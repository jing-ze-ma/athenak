#ifndef EOS_EOS_TABLE_HPP_
#define EOS_EOS_TABLE_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file eos_table.hpp
//! \brief tabulated representation of a general equation of state, on a uniform grid in
//! (log10 rho, log10 T) with cgs units, and the device-side interpolation that every
//! general-EOS kernel call ultimately goes through.
//!
//! WHY A TABLE. The physics (dissociation, Saha ionization) is naturally a function of
//! (rho,T) and costs a nested root find to evaluate directly. Tabulating it collapses
//! that to a handful of flops per call and, more importantly, removes the inner root find
//! entirely: the only iteration left anywhere on the hot path is the outer T(rho,e)
//! inversion, which the code already pays for exactly once per cell per stage.
//!
//! WHAT IS TABULATED, AND WHAT IS NOT. Only the GAS contribution is tabulated, as the
//! specific quantities log10(e_gas/rho) and log10(p_gas/rho). Radiation is added
//! analytically on top, in Eval() below. Keeping radiation out of the table makes the
//! runtime radiation switch exact rather than a second table, and keeps the tabulated
//! surfaces smooth and of modest dynamic range -- an aT^4 term spans sixteen decades
//! across the grid and would dominate the interpolation error everywhere else.
//!
//! INTERPOLATION. Bicubic Hermite, using node values together with node derivatives
//! d/dx, d/dy and the cross derivative d2/dxdy, all evaluated from the analytic model at
//! table build time. The point of tabulating the derivatives rather than differencing the
//! interpolant is CONSISTENCY: chi_rho, chi_T and c_v are returned as the analytic
//! derivatives OF THE INTERPOLANT, so they are exactly the derivatives of the p and e
//! that the solver is actually using, and the identity
//! Gamma_1 = chi_rho + p chi_T^2/(rho T c_v) holds on the interpolated surface and not
//! merely on the underlying physics. Thermodynamic consistency is what the well-balanced
//! scheme and the Riemann solvers depend on; agreement with the true EOS to interpolation
//! accuracy is a separate and weaker requirement.
//!
//! OFF-TABLE STATES are handled by LINEAR extrapolation from the nearest edge, using the
//! edge value and edge derivative. Cubic extrapolation of a Hermite patch diverges
//! violently, so the local coordinates are clamped to the patch and the offset is carried
//! by an explicit first-order term. The result is C1 across the table boundary and
//! preserves both the positivity of p and the monotonicity of e(T) that the root finds
//! below rely on. It is a numerical safety net, not physics: a run that spends time off
//! the table is misconfigured, and BuildEOSTable() says so at startup.

#include <math.h>

#include <string>

#include "athena.hpp"

//----------------------------------------------------------------------------------------
//! \enum EOSTableVar
//! \brief which slice of the table holds what. Each interpolated surface occupies four
//! consecutive slices, in the order (value, d/dx, d/dy, d2/dxdy), where
//! x = log10(rho[cgs]) and y = log10(T[K]). The surfaces are the SPECIFIC quantities
//! log10(e_gas/rho) and log10(p_gas/rho), the mean molecular weight, and log10 of the
//! free-electron fraction n_e/n_tot.
//!
//! ITXE is deliberately NOT read by Eval(): only the resistivity wants it, and charging
//! every EOS call for a fourth Hermite patch would be a pure loss. ElectronFraction()
//! reads it on its own.

enum EOSTableVar {ITE=0, ITP=4, ITMU=8, ITXE=12, ITNVAR=16};

//----------------------------------------------------------------------------------------
//! \struct EOSThermoState
//! \brief everything one table evaluation yields, in cgs. Bundled because callers almost
//! always want several of these and each is nearly free once the patch has been read.

struct EOSThermoState {
  Real e;          // internal energy density, erg/cm^3
  Real p;          // pressure, erg/cm^3
  Real chi_rho;    // dln p/dln rho at constant T
  Real chi_t;      // dln p/dln T at constant rho
  Real cv;         // specific heat at constant volume, erg/g/K
  Real dlne_dlnt;  // dln e/dln T at constant rho; the derivative the T solve needs
  Real mu;         // mean molecular weight, in units of m_u
};

//----------------------------------------------------------------------------------------
//! \struct EOSTable
//! \brief the table itself, plus the interpolation and the inversions built on it. Held
//! by value inside EOS_Data and captured into device lambdas; the Kokkos Views are
//! reference-counted handles, so copying it is cheap.

struct EOSTable {
  bool active = false;       // false leaves every EOS accessor on its analytic branch
  bool radiation = false;    // add the aT^4 terms on top of the tabulated gas

  int nx = 0, ny = 0;        // node counts along x = log10 rho_cgs and y = log10 T_K
  Real xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0;
  Real dx = 1.0, dy = 1.0, dxi = 1.0, dyi = 1.0;

  // Unit conversions, code -> cgs. Duplicated from EOS_Data so that the table is
  // self-contained and can be evaluated without one.
  Real dens_cgs = 1.0, pres_cgs = 1.0, temp_cgs = 1.0;
  Real arad = 7.5657332503e-15;   // radiation constant in cgs
  Real efmax = 0.0;               // global bound on e(rho,pfloor), CODE units

  DvceArray3D<Real> tbl;     // (ITNVAR, ny, nx)
  DvceArray1D<Real> efbnd;   // per-x-cell upper bound on e(rho, pfloor), CODE units

  // ln(10), and the log-space convergence tolerance for the root finds below. The
  // tolerance tracks the working precision: a double-precision run converges to round
  // off in two or three Newton steps from a warm start, a single-precision one has no
  // business chasing digits it does not have.
  static constexpr Real ln10 = 2.302585092994045684;
  static constexpr Real logtol = (sizeof(Real) == 4) ? 1.0e-6 : 1.0e-13;

  //--------------------------------------------------------------------------------------
  //! \fn Real Pow10
  //! \brief 10^a, as an exponential rather than a general power
  KOKKOS_INLINE_FUNCTION
  static Real Pow10(const Real a) { return exp(a*ln10); }

  //--------------------------------------------------------------------------------------
  //! \fn void HermitePatch
  //! \brief evaluate one bicubic Hermite surface and its two first derivatives.
  //!
  //! `iv` is the value slice; the three derivative slices follow it immediately. `ix`
  //! and `iy` index the lower corner of the patch and `u`,`v` are the local coordinates
  //! within it,
  //! already clamped to [0,1] by the caller.
  KOKKOS_INLINE_FUNCTION
  void HermitePatch(const int iv, const int ix, const int iy, const Real u, const Real v,
                    Real &f, Real &fx, Real &fy) const {
    // Hermite basis and its derivative, in each direction, ordered
    // (value at 0, slope at 0, value at 1, slope at 1)
    const Real u2 = u*u, u3 = u2*u;
    const Real v2 = v*v, v3 = v2*v;
    const Real hu[4] = {2.0*u3 - 3.0*u2 + 1.0, u3 - 2.0*u2 + u,
                        -2.0*u3 + 3.0*u2, u3 - u2};
    const Real hv[4] = {2.0*v3 - 3.0*v2 + 1.0, v3 - 2.0*v2 + v,
                        -2.0*v3 + 3.0*v2, v3 - v2};
    const Real du[4] = {6.0*u2 - 6.0*u, 3.0*u2 - 4.0*u + 1.0,
                        -6.0*u2 + 6.0*u, 3.0*u2 - 2.0*u};
    const Real dv[4] = {6.0*v2 - 6.0*v, 3.0*v2 - 4.0*v + 1.0,
                        -6.0*v2 + 6.0*v, 3.0*v2 - 2.0*v};

    // Gather the sixteen corner data into the Hermite coefficient matrix. The derivative
    // entries are scaled by the cell size because the basis above is written in the local
    // coordinates u,v rather than in x,y.
    Real cf[4][4];
    for (int a=0; a<2; ++a) {
      for (int b=0; b<2; ++b) {
        const int i0 = ix + a;
        const int j0 = iy + b;
        cf[2*a][2*b]     = tbl(iv,   j0, i0);
        cf[2*a+1][2*b]   = tbl(iv+1, j0, i0)*dx;
        cf[2*a][2*b+1]   = tbl(iv+2, j0, i0)*dy;
        cf[2*a+1][2*b+1] = tbl(iv+3, j0, i0)*dx*dy;
      }
    }

    f = 0.0; fx = 0.0; fy = 0.0;
    for (int a=0; a<4; ++a) {
      Real sf = 0.0, sd = 0.0;
      for (int b=0; b<4; ++b) {
        sf += cf[a][b]*hv[b];
        sd += cf[a][b]*dv[b];
      }
      f  += hu[a]*sf;
      fx += du[a]*sf;
      fy += hu[a]*sd;
    }
    fx *= dxi;
    fy *= dyi;
    return;
  }

  //--------------------------------------------------------------------------------------
  //! \fn void Interpolate
  //! \brief bicubic interpolation of one surface at an arbitrary (x,y), with linear
  //! extrapolation outside the table
  KOKKOS_INLINE_FUNCTION
  void Interpolate(const int iv, const Real x, const Real y,
                   Real &f, Real &fx, Real &fy) const {
    Real gx = (x - xmin)*dxi;
    Real gy = (y - ymin)*dyi;
    int ix = static_cast<int>(floor(gx));
    int iy = static_cast<int>(floor(gy));
    ix = (ix < 0) ? 0 : ((ix > nx-2) ? nx-2 : ix);
    iy = (iy < 0) ? 0 : ((iy > ny-2) ? ny-2 : iy);
    Real u = gx - static_cast<Real>(ix);
    Real v = gy - static_cast<Real>(iy);
    Real uc = (u < 0.0) ? 0.0 : ((u > 1.0) ? 1.0 : u);
    Real vc = (v < 0.0) ? 0.0 : ((v > 1.0) ? 1.0 : v);
    HermitePatch(iv, ix, iy, uc, vc, f, fx, fy);
    // linear continuation outside the table, at the boundary slope
    if (u != uc) f += (u - uc)*dx*fx;
    if (v != vc) f += (v - vc)*dy*fy;
    return;
  }

  //--------------------------------------------------------------------------------------
  //! \fn void Eval
  //! \brief the full thermodynamic state at a density and temperature given in CGS.
  //!
  //! The gas part comes from the table; the radiation part, if enabled, is added
  //! analytically here. Every derivative returned is the derivative of the value returned
  //! alongside it, so a caller can build Gamma_1 from the standard identity and get an
  //! answer consistent with the p it is about to hand to a Riemann solver.
  KOKKOS_INLINE_FUNCTION
  void Eval(const Real rho, const Real t, EOSThermoState &s) const {
    const Real x = log10(rho);
    const Real y = log10(t);

    Real ev, evx, evy, pv, pvx, pvy, muv, mux, muy;
    Interpolate(ITE, x, y, ev, evx, evy);
    Interpolate(ITP, x, y, pv, pvx, pvy);
    Interpolate(ITMU, x, y, muv, mux, muy);

    // The tabulated surfaces are the SPECIFIC quantities, so the density factor and the
    // derivative it contributes are restored here:
    //   ln p_gas = ln rho + ln10 * P  =>  dln p/dln rho = 1 + dP/dx, dln p/dln T = dP/dy
    // and likewise for e, whose logarithmic temperature derivative is what sets c_v.
    const Real egas = rho*Pow10(ev);
    const Real pgas = rho*Pow10(pv);
    const Real chir_g = 1.0 + pvx;
    const Real chit_g = pvy;
    const Real cv_g = (egas/rho)*evy/t;

    if (radiation) {
      const Real erad = arad*t*t*t*t;
      const Real prad = erad/3.0;
      s.e = egas + erad;
      s.p = pgas + prad;
      // radiation is independent of density at fixed T, so it dilutes chi_rho and pulls
      // chi_T towards its own value of 4
      s.chi_rho = pgas*chir_g/s.p;
      s.chi_t = (pgas*chit_g + 4.0*prad)/s.p;
      s.cv = cv_g + 4.0*erad/(rho*t);
      s.dlne_dlnt = (egas*evy + 4.0*erad)/s.e;
    } else {
      s.e = egas;
      s.p = pgas;
      s.chi_rho = chir_g;
      s.chi_t = chit_g;
      s.cv = cv_g;
      s.dlne_dlnt = evy;
    }
    s.mu = muv;
    return;
  }

  //--------------------------------------------------------------------------------------
  //! \fn Real SolveLog
  //! \brief the safeguarded root find that all three inversions below are built from.
  //!
  //! Solves log10(q(z)) = ltarget for z, where z is log10 of the unknown and q is one of
  //! the tabulated quantities, selected by `mode`:
  //!   0 -- q = e at fixed rho, z = log10 T, derivative dln e/dln T
  //!   1 -- q = p at fixed rho, z = log10 T, derivative chi_T
  //!   2 -- q = p at fixed T,   z = log10 rho, derivative chi_rho
  //! All three targets are monotonically increasing in their unknown, which is what makes
  //! the bracket update valid.
  //!
  //! WHY NOT PLAIN SAFEGUARDED NEWTON. Requiring only that the Newton step stay inside
  //! the bracket is NOT enough. Across an ionization transition dln e/dln T rises sharply
  //! and then falls again, and a Newton step taken from either side overshoots to the
  //! other; both steps land inside the bracket, so a bracket-only guard permits an
  //! endless ping-pong that never converges and never shrinks the bracket. The test below
  //! is the classical one: fall back to bisection whenever the Newton step either leaves
  //! the bracket or fails to at least halve the interval relative to the previous step,
  //! which bounds the iteration count by the bisection rate in the worst case while
  //! keeping quadratic convergence in the common one.
  KOKKOS_INLINE_FUNCTION
  Real SolveLog(const int mode, const Real fixed, const Real ltarget,
                const Real zguess, const Real zlo_in, const Real zhi_in) const {
    Real zlo = zlo_in;
    Real zhi = zhi_in;
    Real z = (zguess > zlo && zguess < zhi) ? zguess : 0.5*(zlo + zhi);

    Real dzold = fabs(zhi - zlo);
    Real dz = dzold;

    Real g, dg;
    EvalResidual(mode, fixed, z, ltarget, g, dg);
    for (int it=0; it<80; ++it) {
      if (g > 0.0) {
        zhi = z;
      } else {
        zlo = z;
      }
      // bisect if Newton would leave the bracket, or is converging too slowly
      bool bisect = (((z - zhi)*dg - g)*((z - zlo)*dg - g) > 0.0) ||
                    (fabs(2.0*g) > fabs(dzold*dg)) || !(dg > 0.0);
      dzold = dz;
      if (bisect) {
        dz = 0.5*(zhi - zlo);
        z = zlo + dz;
      } else {
        dz = g/dg;
        z -= dz;
      }
      if (fabs(dz) < logtol) break;
      EvalResidual(mode, fixed, z, ltarget, g, dg);
    }
    return z;
  }

  //! \fn void EvalResidual
  //! \brief the residual log10(q) - ltarget and its derivative with respect to z, for the
  //! three inversion modes described on SolveLog()
  KOKKOS_INLINE_FUNCTION
  void EvalResidual(const int mode, const Real fixed, const Real z, const Real ltarget,
                    Real &g, Real &dg) const {
    EOSThermoState s;
    if (mode == 2) {
      Eval(Pow10(z), fixed, s);
      g = log10(s.p) - ltarget;
      dg = s.chi_rho;
    } else {
      Eval(fixed, Pow10(z), s);
      // d log10(q)/dz is the same ratio as dln q/dln T, the base of the logarithm
      // cancelling between numerator and denominator
      if (mode == 0) {
        g = log10(s.e) - ltarget;
        dg = s.dlne_dlnt;
      } else {
        g = log10(s.p) - ltarget;
        dg = s.chi_t;
      }
    }
    return;
  }

  //--------------------------------------------------------------------------------------
  //! \fn Real SolveTemperature
  //! \brief invert e(rho,T) for T. Density and energy density in CGS, temperature
  //! returned in Kelvin.
  //!
  //! This is THE expensive operation of the general EOS and the only iteration left on
  //! the hot path. It iterates on y = log10 T, which is the right variable: dln e/dln T
  //! is 3/2 for a monatomic gas and larger inside an ionization zone, but never small, so
  //! the iteration is well conditioned everywhere and a warm start from the previous
  //! stage converges in two or three steps.
  //!
  //! Monotonicity of e in T is what makes the bracket valid, and it is a property of the
  //! interpolant, not merely of the physics: the model is monotonic at every node and the
  //! linear continuation outside the table preserves it.
  KOKKOS_INLINE_FUNCTION
  Real SolveTemperature(const Real rho, const Real etarget, const Real tguess) const {
    // Bracket generously beyond the table: the continuation there is linear in the logs,
    // so a state off the table still has a well defined temperature.
    Real zg = (tguess > 0.0) ? log10(tguess) : -1.0e30;
    return Pow10(SolveLog(0, rho, log10(etarget), zg, ymin - 3.0, ymax + 3.0));
  }

  //--------------------------------------------------------------------------------------
  //! \fn Real SolveTemperatureFromP
  //! \brief invert p(rho,T) for T, in CGS. Driven by chi_T instead of dln e/dln T. Used
  //! by the pressure floor and by problem generators.
  KOKKOS_INLINE_FUNCTION
  Real SolveTemperatureFromP(const Real rho, const Real ptarget,
                             const Real tguess) const {
    Real zg = (tguess > 0.0) ? log10(tguess) : -1.0e30;
    return Pow10(SolveLog(1, rho, log10(ptarget), zg, ymin - 3.0, ymax + 3.0));
  }

  //--------------------------------------------------------------------------------------
  //! \fn Real ElectronFractionCgs
  //! \brief free electrons per particle, n_e/n_tot, at a density and temperature in CGS.
  //!
  //! Kept out of Eval() on purpose: this is a fourth Hermite patch that only the
  //! resistivity needs, and Eval() runs on every EOS call.
  //!
  //! The tabulated surface is log10 of the fraction, because it spans some twenty decades
  //! across a planetary atmosphere and interpolating it linearly would be meaningless.
  KOKKOS_INLINE_FUNCTION
  Real ElectronFractionCgs(const Real rho, const Real t) const {
    Real f, fx, fy;
    Interpolate(ITXE, log10(rho), log10(t), f, fx, fy);
    return Pow10(f);
  }

  //--------------------------------------------------------------------------------------
  //! \fn Real SolveDensity
  //! \brief invert p(rho,T) for rho at fixed T, in CGS.
  //!
  //! Setup-time only. Note that this inversion is ILL POSED when radiation dominates the
  //! pressure, because chi_rho then tends to zero and p stops depending on density at
  //! all; the iteration below still returns a value inside its bracket, but a problem
  //! generator building a background in that regime is asking for something the (p,T)
  //! pair does not determine.
  KOKKOS_INLINE_FUNCTION
  Real SolveDensity(const Real ptarget, const Real t, const Real dguess) const {
    Real zg = (dguess > 0.0) ? log10(dguess) : -1.0e30;
    return Pow10(SolveLog(2, t, log10(ptarget), zg, xmin - 4.0, xmax + 4.0));
  }

  //--------------------------------------------------------------------------------------
  //! \fn Real EnergyFloorBound
  //! \brief a cheap upper bound on the internal energy density at the pressure floor, in
  //! CODE units, for a density given in CODE units.
  //!
  //! Precomputed at table build time: for each x cell the maximum of e(rho,pfloor) over
  //! that cell is stored, so the lookup is one index and one load. A bound rather than
  //! the value, because its only job is to prove that the common case is above the floor
  //! without paying for the inversion. Densities off the table fall back to the global
  //! maximum, which is conservative in the only direction that matters.
  KOKKOS_INLINE_FUNCTION
  Real EnergyFloorBound(const Real d) const {
    Real gx = (log10(d*dens_cgs) - xmin)*dxi;
    if (gx < 0.0 || gx > static_cast<Real>(nx-1)) return efmax;
    int ix = static_cast<int>(floor(gx));
    ix = (ix > nx-2) ? nx-2 : ix;
    return efbnd(ix);
  }
};

//----------------------------------------------------------------------------------------
// host-side table construction; implemented in eos_table.cpp

class ParameterInput;
void BuildEOSTable(EOSTable &tbl, ParameterInput *pin, const std::string &block,
                   const Real dens_cgs, const Real pres_cgs, const Real temp_cgs,
                   const Real pfloor);

#endif // EOS_EOS_TABLE_HPP_
