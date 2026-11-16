#ifndef EOS_EOS_TABLE_HPP_
#define EOS_EOS_TABLE_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file eos_table.hpp
//! \brief tabulated general equation of state on a uniform (log10 rho, log10 T) grid in
//! cgs, and the device interpolation every general-EOS call goes through.
//!
//! Only the gas part is tabulated, as log10(e_gas/rho) and log10(p_gas/rho); radiation is
//! added analytically in Eval().  Interpolation is bicubic Hermite with node derivatives
//! from the analytic model, so chi_rho, chi_T and c_v are the derivatives of the
//! interpolated p and e and Gamma_1 = chi_rho + p chi_T^2/(rho T c_v) holds on the
//! interpolant.  Off the table the surfaces continue linearly from the edge (C1, keeps p
//! > 0 and e(T) monotone); BuildEOSTable reports runs that need it.

#include <math.h>

#include <string>

#include "athena.hpp"
#include "utils/rad_taper.hpp"

//----------------------------------------------------------------------------------------
//! \enum EOSTableVar
//! \brief table slices: each surface occupies four consecutive slices (value, d/dx, d/dy,
//! d2/dxdy), x = log10(rho[cgs]), y = log10(T[K]).  Surfaces: log10(e_gas/rho),
//! log10(p_gas/rho), mu and log10(n_e/n_tot); ITXE is read only by ElectronFractionCgs.

enum EOSTableVar {ITE=0, ITP=4, ITMU=8, ITXE=12, ITNVAR=16};

//----------------------------------------------------------------------------------------
//! \struct EOSThermoState
//! \brief thermodynamic state of one table evaluation, in cgs

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
//! \brief the table with its interpolation and inversions; held by value in EOS_Data and
//! captured into device lambdas (the Views are reference-counted handles)

struct EOSTable {
  bool active = false;       // false leaves every EOS accessor on its analytic branch
  bool radiation = false;    // add the aT^4 terms on top of the tabulated gas

  int nx = 0, ny = 0;        // node counts along x = log10 rho_cgs and y = log10 T_K
  Real xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0;
  Real dx = 1.0, dy = 1.0, dxi = 1.0, dyi = 1.0;

  // code -> cgs unit factors (as in EOS_Data, so the table is self-contained)
  Real dens_cgs = 1.0, pres_cgs = 1.0, temp_cgs = 1.0;
  Real arad = 7.5657332503e-15;   // radiation constant in cgs
  // radiation taper (utils/rad_taper.hpp): the aT^4 terms are multiplied by a smoothstep
  // w(log10 rho), 1 for rho >= 10^rad_lrho_hi and 0 for rho <= 10^rad_lrho_lo
  bool rad_taper = false;
  Real rad_lrho_lo = 0.0, rad_lrho_hi = 0.0;   // log10 of rho_lo, rho_hi, in cgs
  // temperature gate of the taper (<block>/eos_rad_t_lo, eos_rad_t_hi): w = max(w_rho,
  // w_T), so hot cells keep their radiation (rad_taper::WeightGated)
  bool rad_tgate = false;
  Real rad_lt_lo = 0.0, rad_lt_hi = 0.0;       // log10 of T_lo, T_hi, in Kelvin
  // composition the table was built with, for consistency checks against opacity fits
  bool metal_ion = false;         // were the metal electron donors included?
  Real metal_mh = 0.0;            // [M/H] in dex used for those donors
  Real efmax = 0.0;               // global bound on e(rho,pfloor), CODE units

  // table (ny, nx, ITNVAR), variable index fastest so that one Hermite node is contiguous
  DvceArray3D<Real> tbl;
  DvceArray1D<Real> efbnd;   // per-x-cell upper bound on e(rho, pfloor), CODE units
  // HostSpace mirrors of tbl and efbnd for host callers (problem generators)
  DvceArray3D<Real>::HostMirror tbl_h;
  DvceArray1D<Real>::HostMirror efbnd_h;

  // ln(10) and the log-space tolerance of the root finds (tracks the working precision)
  static constexpr Real ln10 = 2.302585092994045684;
  static constexpr Real logtol = (sizeof(Real) == 4) ? 1.0e-6 : 1.0e-13;

  //--------------------------------------------------------------------------------------
  //! \fn Real Pow10
  //! \brief 10^a, as an exponential rather than a general power
  KOKKOS_INLINE_FUNCTION
  static Real Pow10(const Real a) { return exp(a*ln10); }

  //--------------------------------------------------------------------------------------
  //! \fn void GatherPatch
  //! \brief gathers the sixteen corner data of one Hermite patch from table view t (tbl
  //! on the device, tbl_h on the host)
  template <class TV>
  KOKKOS_INLINE_FUNCTION
  void GatherPatch(const TV &t, const int iv, const int ix, const int iy,
                   Real (&cf)[4][4]) const {
    for (int a=0; a<2; ++a) {
      for (int b=0; b<2; ++b) {
        const int i0 = ix + a;
        const int j0 = iy + b;
        cf[2*a][2*b]     = t(j0, i0, iv  );
        cf[2*a+1][2*b]   = t(j0, i0, iv+1)*dx;
        cf[2*a][2*b+1]   = t(j0, i0, iv+2)*dy;
        cf[2*a+1][2*b+1] = t(j0, i0, iv+3)*dx*dy;
      }
    }
  }

  //--------------------------------------------------------------------------------------
  //! \fn void HermitePatch
  //! \brief one bicubic Hermite surface and its first derivatives.  iv is the value slice
  //! (derivative slices follow), (ix, iy) the patch corner, (u, v) in [0, 1] the local
  //! coordinates; WANT_FX / WANT_FY skip derivatives the caller does not read.
  template <bool WANT_FX, bool WANT_FY>
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
    Real du[4] = {0.0, 0.0, 0.0, 0.0};
    Real dv[4] = {0.0, 0.0, 0.0, 0.0};
    if (WANT_FX) {
      du[0] = 6.0*u2 - 6.0*u;   du[1] = 3.0*u2 - 4.0*u + 1.0;
      du[2] = -6.0*u2 + 6.0*u;  du[3] = 3.0*u2 - 2.0*u;
    }
    if (WANT_FY) {
      dv[0] = 6.0*v2 - 6.0*v;   dv[1] = 3.0*v2 - 4.0*v + 1.0;
      dv[2] = -6.0*v2 + 6.0*v;  dv[3] = 3.0*v2 - 2.0*v;
    }

    // corner data -> Hermite coefficients; derivatives scaled by the cell size (the basis
    // is in the local u, v)
    Real cf[4][4];
    KOKKOS_IF_ON_DEVICE((GatherPatch(tbl, iv, ix, iy, cf);))
    KOKKOS_IF_ON_HOST((GatherPatch(tbl_h, iv, ix, iy, cf);))

    f = 0.0; fx = 0.0; fy = 0.0;
    for (int a=0; a<4; ++a) {
      Real sf = 0.0, sd = 0.0;
      for (int b=0; b<4; ++b) {
        sf += cf[a][b]*hv[b];
        if (WANT_FY) sd += cf[a][b]*dv[b];
      }
      f  += hu[a]*sf;
      if (WANT_FX) fx += du[a]*sf;
      if (WANT_FY) fy += hu[a]*sd;
    }
    fx *= dxi;
    fy *= dyi;
    return;
  }

  //--------------------------------------------------------------------------------------
  //! \fn void Interpolate
  //! \brief bicubic interpolation of one surface at an arbitrary (x,y), with linear
  //! extrapolation outside the table
  template <bool WANT_FX = true, bool WANT_FY = true>
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
    // off the table the linear continuation needs both derivatives
    if (u != uc || v != vc) {
      HermitePatch<true, true>(iv, ix, iy, uc, vc, f, fx, fy);
      // linear continuation outside the table, at the boundary slope
      if (u != uc) f += (u - uc)*dx*fx;
      if (v != vc) f += (v - vc)*dy*fy;
    } else {
      HermitePatch<WANT_FX, WANT_FY>(iv, ix, iy, uc, vc, f, fx, fy);
    }
    return;
  }

  //--------------------------------------------------------------------------------------
  //! \fn void EvalImpl
  //! \brief thermodynamic state at (rho, T) in cgs: gas from the table, radiation added
  //! analytically; every derivative is that of the value returned beside it.  WANT_MU
  //! selects the mu surface (not needed on the hot path).
  template <bool WANT_MU>
  KOKKOS_INLINE_FUNCTION
  void EvalImpl(const Real rho, const Real t, EOSThermoState &s) const {
    EvalFromLogs<WANT_MU>(log10(rho), log10(t), rho, t, s);
  }

  //! \fn void EvalFromLogs
  //! \brief body of EvalImpl taking both log10 and linear (rho, T)
  template <bool WANT_MU>
  KOKKOS_INLINE_FUNCTION
  void EvalFromLogs(const Real x, const Real y, const Real rho, const Real t,
                    EOSThermoState &s) const {
    Real ev, evx, evy, pv, pvx, pvy, muv = 0.0, mux, muy;
    Interpolate<false, true>(ITE, x, y, ev, evx, evy);   // evx is never used
    Interpolate<true, true>(ITP, x, y, pv, pvx, pvy);
    if (WANT_MU) Interpolate<false, false>(ITMU, x, y, muv, mux, muy);

    // restore the density factor of the specific surfaces:
    //   ln p_gas = ln rho + ln10 P  =>  dln p/dln rho = 1 + dP/dx, dln p/dln T = dP/dy
    const Real egas = rho*Pow10(ev);
    const Real pgas = rho*Pow10(pv);
    const Real chir_g = 1.0 + pvx;
    const Real chit_g = pvy;
    const Real cv_g = (egas/rho)*evy/t;

    if (radiation) {
      const Real erad0 = arad*t*t*t*t;
      Real w = 1.0, dwdx = 0.0, dwdy = 0.0;
      if (rad_taper) {
        rad_taper::WeightGated(x, rad_lrho_lo, rad_lrho_hi, y, rad_lt_lo, rad_lt_hi,
                               rad_tgate, w, dwdx, dwdy);
      }
      const Real erad = rad_taper ? w*erad0 : erad0;
      const Real prad = erad/3.0;
      s.e = egas + erad;
      s.p = pgas + prad;
      // radiation is independent of density at fixed T, so it dilutes chi_rho and pulls
      // chi_T towards its own value of 4
      s.chi_rho = pgas*chir_g/s.p;
      // with the taper, w depends on density: d(w prad0)/dln rho = prad0 dw/dx/ln10
      if (rad_taper) s.chi_rho += (erad0/3.0)*dwdx*M_LOG10E/s.p;
      s.chi_t = (pgas*chit_g + 4.0*prad)/s.p;
      s.cv = cv_g + 4.0*erad/(rho*t);
      s.dlne_dlnt = (egas*evy + 4.0*erad)/s.e;
      // the temperature gate adds dw/dlnT to the temperature derivatives:
      // d(w erad0)/dlnT = 4 w erad0 + erad0 dw/dy/ln10 (same arithmetic as EvalPOnly and
      // EvalEOnly)
      if (dwdy != 0.0) {
        const Real dgt = erad0*dwdy*M_LOG10E;
        s.chi_t += (erad0/3.0)*dwdy*M_LOG10E/s.p;
        s.cv += dgt/(rho*t);
        s.dlne_dlnt += dgt/s.e;
      }
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

  //! \fn void Eval
  //! \brief the full state, mu included
  KOKKOS_INLINE_FUNCTION
  void Eval(const Real rho, const Real t, EOSThermoState &s) const {
    EvalImpl<true>(rho, t, s);
  }

  //! \fn void EvalNoMu
  //! \brief the full state without mu (s.mu = 0)
  KOKKOS_INLINE_FUNCTION
  void EvalNoMu(const Real rho, const Real t, EOSThermoState &s) const {
    EvalImpl<false>(rho, t, s);
  }

  //--------------------------------------------------------------------------------------
  //! \fn Real SolveLog
  //! \brief safeguarded Newton solve of log10(q(z)) = ltarget, z = log10 of the unknown:
  //!   MODE 0: q = e at fixed rho, z = log10 T
  //!   MODE 1: q = p at fixed rho, z = log10 T
  //!   MODE 2: q = p at fixed T,   z = log10 rho
  //! All targets increase monotonically in z.  The step falls back to bisection when
  //! Newton leaves the bracket or fails to halve the previous step, which stops the
  //! ping-pong of a bracket-only guard across ionization zones.
  template <int MODE>
  KOKKOS_INLINE_FUNCTION
  Real SolveLog(const Real xfixed, const Real ltarget,
                const Real zguess, const Real zlo_in, const Real zhi_in) const {
    Real zlo = zlo_in;
    Real zhi = zhi_in;
    Real z = (zguess > zlo && zguess < zhi) ? zguess : 0.5*(zlo + zhi);

    Real dzold = fabs(zhi - zlo);
    Real dz = dzold;

    Real g, dg;
    EvalResidual<MODE>(xfixed, z, ltarget, g, dg);
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
      EvalResidual<MODE>(xfixed, z, ltarget, g, dg);
    }
    return z;
  }

  //--------------------------------------------------------------------------------------
  //! \fn void EvalEOnly
  //! \brief e and dln e/dln T from the energy patch alone (Eval's arithmetic, bitwise
  //! equal)
  KOKKOS_INLINE_FUNCTION
  void EvalEOnly(const Real rho, const Real t, Real &e, Real &dlne_dlnt) const {
    Real ev, evx, evy;
    const Real xr = log10(rho);
    const Real yr = log10(t);
    Interpolate<false, true>(ITE, xr, yr, ev, evx, evy);
    const Real egas = rho*Pow10(ev);
    if (radiation) {
      Real erad = arad*t*t*t*t;
      Real dgt = 0.0;
      if (rad_taper) {
        Real w, dwdx, dwdy;
        rad_taper::WeightGated(xr, rad_lrho_lo, rad_lrho_hi, yr, rad_lt_lo, rad_lt_hi,
                               rad_tgate, w, dwdx, dwdy);
        dgt = erad*dwdy*M_LOG10E;   // erad0 dw/dlnT, before the weight is applied
        erad *= w;
      }
      e = egas + erad;
      dlne_dlnt = (egas*evy + 4.0*erad)/e;
      if (dgt != 0.0) dlne_dlnt += dgt/e;
    } else {
      e = egas;
      dlne_dlnt = evy;
    }
    return;
  }

  //--------------------------------------------------------------------------------------
  //! \fn void EvalPOnly
  //! \brief p, chi_rho and chi_T from the pressure patch alone (Eval's arithmetic,
  //! bitwise equal)
  KOKKOS_INLINE_FUNCTION
  void EvalPOnly(const Real rho, const Real t,
                 Real &p, Real &chi_rho, Real &chi_t) const {
    Real pv, pvx, pvy;
    const Real xr = log10(rho);
    const Real yr = log10(t);
    Interpolate(ITP, xr, yr, pv, pvx, pvy);
    const Real pgas = rho*Pow10(pv);
    const Real chir_g = 1.0 + pvx;
    const Real chit_g = pvy;
    if (radiation) {
      const Real erad0 = arad*t*t*t*t;
      Real w = 1.0, dwdx = 0.0, dwdy = 0.0;
      if (rad_taper) {
        rad_taper::WeightGated(xr, rad_lrho_lo, rad_lrho_hi, yr, rad_lt_lo, rad_lt_hi,
                               rad_tgate, w, dwdx, dwdy);
      }
      const Real erad = rad_taper ? w*erad0 : erad0;
      const Real prad = erad/3.0;
      p = pgas + prad;
      chi_rho = pgas*chir_g/p;
      if (rad_taper) chi_rho += (erad0/3.0)*dwdx*M_LOG10E/p;
      chi_t = (pgas*chit_g + 4.0*prad)/p;
      if (dwdy != 0.0) chi_t += ((erad0/3.0)*dwdy*M_LOG10E)/p;
    } else {
      p = pgas;
      chi_rho = chir_g;
      chi_t = chit_g;
    }
    return;
  }

  //! \fn void EvalResidual
  //! \brief residual g = log10(q) - ltarget and dg/dz of SolveLog, kept in log space
  //! (xfixed = log10 of the fixed quantity): log10 e = log10 rho + E(x, y); only the
  //! radiation term needs exponentials
  template <int MODE>
  KOKKOS_INLINE_FUNCTION
  void EvalResidual(const Real xfixed, const Real z, const Real ltarget,
                    Real &g, Real &dg) const {
    // d log10(q)/dz is the same ratio as dln q/dln T, the base of the logarithm
    // cancelling between numerator and denominator
    const Real x = (MODE == 2) ? z : xfixed;      // log10 rho
    const Real y = (MODE == 2) ? xfixed : z;      // log10 T
    Real qv, qvx = 0.0, qvy = 0.0;
    // mode 2 differentiates in density, the other two in temperature
    Interpolate<(MODE == 2), (MODE != 2)>((MODE == 0) ? ITE : ITP, x, y, qv, qvx, qvy);
    // log10 of the volumetric quantity, without ever leaving log space
    g = x + qv - ltarget;
    if (MODE == 0) {
      dg = qvy;
    } else {
      dg = (MODE == 2) ? (1.0 + qvx) : qvy;
    }
    if (radiation) {
      // radiation is additive in the linear quantity (arithmetic as in Eval)
      const Real t = Pow10(y);
      Real erad = arad*t*t*t*t;
      Real w = 1.0, dwdx = 0.0, dwdy = 0.0;
      if (rad_taper) {
        rad_taper::WeightGated(x, rad_lrho_lo, rad_lrho_hi, y, rad_lt_lo, rad_lt_hi,
                               rad_tgate, w, dwdx, dwdy);
        // dwdx is with respect to log10 rho, so only the density inversion sees it
        dwdx *= erad*M_LOG10E;
        // ...and dwdy with respect to log10 T, which only the two T inversions see
        dwdy *= erad*M_LOG10E;
        erad *= w;
      }
      const Real qrad = (MODE == 0) ? erad : erad/3.0;
      const Real qgas = Pow10(g + ltarget);
      const Real qtot = qgas + qrad;
      const Real f = (MODE == 0) ? 4.0 : ((MODE == 2) ? 0.0 : 4.0);
      if (MODE == 2) {
        dg = qgas*(1.0 + qvx)/qtot;
        if (rad_taper) dg += (dwdx/3.0)/qtot;
      } else {
        dg = (qgas*qvy + f*qrad)/qtot;
        // gate term dw/dlnT >= 0: the target stays monotone in T
        if (dwdy != 0.0) dg += ((MODE == 0) ? dwdy : dwdy/3.0)/qtot;
      }
      g = log10(qtot) - ltarget;
    }
    return;
  }

  //--------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------
  //! \fn Real ClampLogT
  //! \brief holds a solved log10 T in the tabulated range [ymin, ymax] and flags the
  //! clamp. The root-find bracket extends 3 decades beyond the table, and an out-of-range
  //! or non-finite target pins at its edge (a saturated state with a huge sound speed);
  //! the nearest tabulated row is returned instead.  States inside the table are
  //! unchanged.
  KOKKOS_INLINE_FUNCTION
  Real ClampLogT(const Real z, bool &clamped) const {
    if (!(z >= ymin)) {   // below the table, or NaN
      clamped = true;
      return ymin;
    }
    if (z > ymax) {
      clamped = true;
      return ymax;
    }
    return z;
  }

  //! \fn Real ClampedSolveLogT
  //! \brief mode-0/1 root find, clamped; a non-finite target is caught before the solve
  //! (a NaN residual would walk the bisection to the top of the bracket)
  template <int MODE>
  KOKKOS_INLINE_FUNCTION
  Real ClampedSolveLogT(const Real lfixed, const Real ltarget, const Real zguess,
                        bool &clamped) const {
    if (!Kokkos::isfinite(ltarget)) {
      clamped = true;
      return (ltarget > 0.0) ? ymax : ymin;   // +inf -> top, -inf/NaN -> bottom
    }
    return ClampLogT(SolveLog<MODE>(lfixed, ltarget, zguess, ymin - 3.0, ymax + 3.0),
                     clamped);
  }

  //! \fn Real SolveLogTemperature
  //! \brief inverts e(rho, T) for T in log space (log10 rho in, log10 T out), clamped
  KOKKOS_INLINE_FUNCTION
  Real SolveLogTemperature(const Real lrho, const Real letarget, const Real zguess,
                           bool &clamped) const {
    return ClampedSolveLogT<0>(lrho, letarget, zguess, clamped);
  }

  KOKKOS_INLINE_FUNCTION
  Real SolveLogTemperature(const Real lrho, const Real letarget,
                           const Real zguess) const {
    bool clamped = false;
    return ClampedSolveLogT<0>(lrho, letarget, zguess, clamped);
  }

  //! \fn Real SolveTemperature
  //! \brief inverts e(rho, T) for T in cgs (Kelvin), clamped; iterates on log10 T, a
  //! well-conditioned variable, from the warm start tguess
  KOKKOS_INLINE_FUNCTION
  Real SolveTemperature(const Real rho, const Real etarget, const Real tguess,
                        bool &clamped) const {
    Real zg = (tguess > 0.0) ? log10(tguess) : -1.0e30;
    return Pow10(ClampedSolveLogT<0>(log10(rho), log10(etarget), zg, clamped));
  }

  KOKKOS_INLINE_FUNCTION
  Real SolveTemperature(const Real rho, const Real etarget, const Real tguess) const {
    bool clamped = false;
    return SolveTemperature(rho, etarget, tguess, clamped);
  }

  //--------------------------------------------------------------------------------------
  //! \fn Real SolveTemperatureFromP
  //! \brief inverts p(rho, T) for T in cgs (pressure floor, problem generators); clamped
  //! like the energy inversion, so the floor never puts a cell below the table
  KOKKOS_INLINE_FUNCTION
  Real SolveTemperatureFromP(const Real rho, const Real ptarget, const Real tguess,
                             bool &clamped) const {
    Real zg = (tguess > 0.0) ? log10(tguess) : -1.0e30;
    return Pow10(ClampedSolveLogT<1>(log10(rho), log10(ptarget), zg, clamped));
  }

  KOKKOS_INLINE_FUNCTION
  Real SolveTemperatureFromP(const Real rho, const Real ptarget,
                             const Real tguess) const {
    bool clamped = false;
    return SolveTemperatureFromP(rho, ptarget, tguess, clamped);
  }

  //--------------------------------------------------------------------------------------
  //! \fn Real ElectronFractionCgs
  //! \brief free electrons per particle n_e/n_tot at (rho, T) in cgs (log10 surface)
  KOKKOS_INLINE_FUNCTION
  Real ElectronFractionCgs(const Real rho, const Real t) const {
    Real f, fx, fy;
    Interpolate(ITXE, log10(rho), log10(t), f, fx, fy);
    return Pow10(f);
  }

  //--------------------------------------------------------------------------------------
  //! \fn Real SolveDensity
  //! \brief inverts p(rho, T) for rho at fixed T in cgs (setup time; ill posed where
  //! radiation dominates the pressure)
  KOKKOS_INLINE_FUNCTION
  Real SolveDensity(const Real ptarget, const Real t, const Real dguess) const {
    Real zg = (dguess > 0.0) ? log10(dguess) : -1.0e30;
    return Pow10(SolveLog<2>(log10(t), log10(ptarget), zg, xmin - 4.0, xmax + 4.0));
  }

  //--------------------------------------------------------------------------------------
  //! \fn Real EnergyFloorBound
  //! \brief upper bound on e(rho, pfloor) in code units: the per-x-cell maximum stored at
  //! build time (global maximum off the table)
  KOKKOS_INLINE_FUNCTION
  Real EnergyFloorBound(const Real d) const {
    Real gx = (log10(d*dens_cgs) - xmin)*dxi;
    if (gx < 0.0 || gx > static_cast<Real>(nx-1)) return efmax;
    int ix = static_cast<int>(floor(gx));
    ix = (ix > nx-2) ? nx-2 : ix;
    Real b = efmax;
    KOKKOS_IF_ON_DEVICE((b = efbnd(ix);))
    KOKKOS_IF_ON_HOST((b = efbnd_h(ix);))
    return b;
  }
};

//----------------------------------------------------------------------------------------
// host-side table construction; implemented in eos_table.cpp

class ParameterInput;
void BuildEOSTable(EOSTable &tbl, ParameterInput *pin, const std::string &block,
                   const Real dens_cgs, const Real pres_cgs, const Real temp_cgs,
                   const Real pfloor);

#endif // EOS_EOS_TABLE_HPP_
