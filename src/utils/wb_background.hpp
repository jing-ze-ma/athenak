#ifndef UTILS_WB_BACKGROUND_HPP_
#define UTILS_WB_BACKGROUND_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file wb_background.hpp
//! \brief local hydrostatic background states for the deviation-based well-balanced
//! scheme, for a GENERAL equation of state.
//!
//! The well-balanced scheme of Kappeli & Mishra (2014, 2016) subtracts a local
//! hydrostatic background from the state before reconstruction, reconstructs only the
//! deviation, and adds the background back at the interfaces.  Its ideal-gas
//! implementation (Hydro::getWBerho / MHD::getWBerho) writes that background in closed
//! form, because for p = (gamma-1) e each of the three thermodynamic constraints below
//! makes some power or logarithm of e and d exactly linear in the potential.  None of
//! those closed forms survives a general EOS.
//!
//! What DOES survive is the underlying statement.  Hydrostatic equilibrium is
//!
//!     dp/dPhi = -d,
//!
//! and closing it requires one thermodynamic constraint.  The three options are:
//!
//!   isodensity   d fixed        =>  p is linear in Phi.  Closed form, exact.
//!   isothermal   T fixed        =>  dln d/dPhi = -d/(p chi_rho),  e = e(d,T).
//!   isentropic   s fixed        =>  dd/dPhi = -d/cs^2,  de/dPhi = -(e+p)/cs^2.
//!
//! The isentropic pair is just h + Phi = const (Bernoulli) written differentially, using
//! dh = dp/d and dp = cs^2 dd along an isentrope.  Every coefficient is a function of
//! (d,e) that the EOS already provides, so NO entropy function and NO root find is
//! needed -- the background is obtained by direct integration.
//!
//! Each expression reduces to the ideal-gas closed form when the EOS is a gamma law:
//! chi_rho = 1 and p = d T give dln d/dPhi = -1/T (the exp() branch of getWBerho), and
//! cs^2 = gamma p/d gives dln d/dPhi = -d/(gamma p), whose exact integral is the
//! d^(gamma-1) branch.
//!
//! ACCURACY.  isodensity and isothermal are integrated exactly: the first is a closed
//! form, and the second has a coefficient that is constant along the background by
//! construction (T is the constrained variable), so a single exponential step is exact
//! for any EOS whose chi_rho and T-relation are evaluated at the anchor.  Only the
//! isentropic pair has a genuinely varying coefficient, and it is integrated with RK4.
//!
//! WELL-BALANCEDNESS DOES NOT DEPEND ON THAT ACCURACY.  The gravitational source term is
//! discretised as the background's own pressure difference across the cell -- not as
//! -d g dx -- so a state that equals the discrete background produces a momentum flux
//! difference that cancels the source identically, for ANY background, however it was
//! obtained.  Integration accuracy controls only how much of the true stratification is
//! removed before reconstruction, i.e. the size of the deviation, not whether the
//! balance is exact.  The one hard requirement is that the reconstruction and the source
//! term call this routine with identical arguments.
//!
//! DIFFERENCE FROM THE IDEAL-GAS BRANCH.  getWBerho() reaches the neighbouring cell
//! centres in two steps, using the neighbour's own state for the coefficient on the
//! outer half cell.  The routine here integrates every target point from the anchor cell
//! alone, which keeps the background a single thermodynamically consistent profile.  The
//! two constructions agree exactly whenever the flow actually satisfies the assumed
//! constraint -- which is precisely the regime in which the scheme is well balanced --
//! and differ at O(dPhi * gradient of the constraint variable) otherwise.

#include <math.h>

#include "athena.hpp"
#include "eos/eos.hpp"

//----------------------------------------------------------------------------------------
//! \struct WBState
//! \brief background density, internal energy density and pressure at one point.

struct WBState {
  Real d, e, p;
};

//----------------------------------------------------------------------------------------
//! \enum WBVar
//! \brief which channel of the background a call site wants.  This cannot be inferred
//! from the array index the reconstruction is looping over: the pressure channel lives in
//! wder, whose DerivedIndex values collide numerically with IDN/IEN.

enum WBVar {wb_dens = 0, wb_eint = 1, wb_pres = 2};

//! \fn void WBReadCache
//! \brief the five background states of one channel from the per-cell cache
//! (Hydro/MHD::wbq0, layout (m, 5*var + {im1,imh,i,iph,ip1}, k, j, i)), built once per
//! stage by BuildWBCache() so the density, energy and pressure channels and the source
//! term all read the same walk instead of each repeating it.
KOKKOS_INLINE_FUNCTION
void WBReadCache(const DvceArray5D<Real> &c, const int var, const int m, const int k,
                 const int j, const int i, Real &q0_im1, Real &q0_imh, Real &q0_i,
                 Real &q0_iph, Real &q0_ip1) {
  const int b = 5*var;
  q0_im1 = c(m,b,k,j,i); q0_imh = c(m,b+1,k,j,i); q0_i = c(m,b+2,k,j,i);
  q0_iph = c(m,b+3,k,j,i); q0_ip1 = c(m,b+4,k,j,i);
}

//! \fn Real WBT
//! \brief the cached temperature of a cell (Hydro/MHD::wtemp, solved once by ConsToPrim),
//! or -1 ("unknown, solve for it") when the array is not allocated -- an ideal gas.
KOKKOS_INLINE_FUNCTION
Real WBT(const DvceArray4D<Real> &w, const int m, const int k, const int j, const int i) {
  return (w.extent_int(0) > 0) ? w(m,k,j,i) : -1.0;
}

//----------------------------------------------------------------------------------------
//! \fn int WBOptionNumber
//! \brief maps the WBOption enum onto the integer branch used below, resolving `adaptive`
//! against the local state.  `adaptive` picks whichever of the isothermal and isentropic
//! backgrounds the flow is closer to, by comparing the relative variation of T with the
//! variation of the entropy across the stencil.  Only entropy DIFFERENCES are needed,
//!
//!     ds = c_v dln T - (p chi_T/(d T)) dln d,
//!
//! so no absolute entropy function is required.  For an ideal gas this reduces to the
//! ratio of dT/T to ds/s used by the ideal-gas branch.

//! `tguess` is a temperature to start the inversions from -- the anchor cell's, which
//! WBBackgroundStencil has already solved for. Every state in the stencil is within a
//! half cell of the anchor, so it is a good guess; a bad one costs iterations, never
//! accuracy, and an ideal gas ignores it entirely.

KOKKOS_INLINE_FUNCTION
int WBOptionNumber(const EOS_Data &eos, const WBOption wb_option,
                   const Real d_m, const Real e_m, const Real d_p, const Real e_p,
                   const Real d_i, const Real e_i, const Real tguess = -1.0,
                   const Real t_mc = -1.0, const Real t_pc = -1.0) {
  switch (wb_option) {
    case WBOption::isodensity:
      return 0;
    case WBOption::isothermal:
      return 1;
    case WBOption::isentropic:
      return 2;
    case WBOption::polytropic:
      return 3;
    case WBOption::isentropic_dt:
      return 4;
    case WBOption::adaptive:
    case WBOption::adaptive_fast: {
      // the three temperatures are ConsToPrim's cached ones when the caller has them
      Real t_m = (t_mc > 0.0) ? t_mc : eos.Temperature(d_m, e_m, tguess);
      Real t_p = (t_pc > 0.0) ? t_pc : eos.Temperature(d_p, e_p, tguess);
      Real t_i = (tguess > 0.0) ? tguess : eos.Temperature(d_i, e_i);
      Real dlnt = log(t_p/t_m);
      Real dlnd = log(d_p/d_m);
      // entropy difference across the stencil, in units of c_v, evaluated at cell i. All
      // three coefficients are taken at the temperature already solved for above: asking
      // for them without it would repeat the same root find three more times.
      Real cv = eos.SpecificHeatCv(d_i, e_i, t_i);
      Real p_i = eos.Pressure(d_i, e_i, t_i);
      Real dsdivcv = dlnt - (p_i*eos.ChiT(d_i,e_i,t_i)/(d_i*t_i*cv))*dlnd;
      // adaptive_fast: the isentrope is walked in (d,T) with no root find
      const int isen = (wb_option == WBOption::adaptive_fast) ? 4 : 2;
      return (fabs(dlnt) > fabs(dsdivcv)) ? isen : 1;
    }
    default:
      return 1;
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real WBEnergyFromEnthalpy
//! \brief solves (e + p(d,e))/d = h for e at fixed d, by Newton iteration.
//!
//! This is what turns h + Phi = const from an identity into something usable.  Along a
//! hydrostatic isentrope dh = dp/d and dp/dPhi = -d, so h + Phi is EXACTLY conserved, with
//! no EOS entering at all.  That single relation does not determine the state -- it is one
//! equation in (d,e), and the missing half is which adiabat -- but once d has been obtained
//! some other way, it pins e down exactly, and the background then satisfies hydrostatic
//! equilibrium to machine precision rather than to the accuracy of the density step.
//!
//! f(e) = e + p(d,e) - d h is monotonically increasing in e, and f'(e) = 1 + (dp/de)_d with
//! (dp/de)_d = p chi_T/(d T c_v).  For an ideal gas f is linear and the first iteration is
//! exact, after which the residual is zero and the loop leaves e untouched -- so a gamma
//! law reproduces the closed-form answer bitwise.

//! `t_out` receives the temperature of the last iterate, which the caller reuses instead
//! of solving the root find again on the converged e. The two differ only by the residual
//! the loop stopped at -- 1e-15 relative at most, and exactly zero when it stopped on
//! f == 0, which is what an ideal gas does on its second pass.

KOKKOS_INLINE_FUNCTION
Real WBEnergyFromEnthalpy(const EOS_Data &eos, const Real d, const Real h,
                          const Real e_guess, Real &t_out, const Real tguess = -1.0,
                          Real *p_out = nullptr) {
  Real e = e_guess;
  // Each iterate is a small step from the last, so the previous iterate's temperature is
  // very nearly the answer for the next one. Carrying it forward turns every inversion
  // after the first into a couple of Newton steps instead of a fresh bracket, and the
  // first one starts from the caller's guess.
  Real tg = tguess;
  for (int it=0; it<10; ++it) {
    // One temperature solve per iteration, shared by the residual and the derivative.
    // This loop is the hottest EOS consumer in the dynamic well-balanced path -- up to
    // ten iterations per half cell per direction -- so evaluating p, chi_T and c_v
    // without a temperature would cost four root finds per iteration instead of one.
    Real t = eos.Temperature(d, e, tg);
    tg = t;
    t_out = t;
    // one table evaluation serves p, chi_T and c_v
    Real et, p, chir, chit, cv;
    eos.ThermoAt(d, t, et, p, chir, chit, cv);
    if (p_out != nullptr) {*p_out = p;}
    Real f = e + p - d*h;
    if (f == 0.0) {break;}
    Real dpde = p*chit/(d*t*cv);
    Real de = -f/(1.0 + dpde);
    // keep the iterate positive; a general EOS is undefined at e <= 0
    e = (e + de > 0.0) ? (e + de) : 0.5*e;
    if (fabs(de) <= 1.0e-15*fabs(e)) {break;}
  }
  return e;
}

//----------------------------------------------------------------------------------------
//! \fn void WBAdvance
//! \brief advances the background state by dphi along one segment of the stencil walk,
//! with the thermodynamic coefficient frozen at (d_c,e_c).
//!
//! The ideal-gas getWBerho() reaches a neighbouring cell centre in TWO segments: the inner
//! half cell uses the anchor cell's own coefficient, the outer half cell the neighbour's.
//! That is not an accident -- when the flow does not exactly satisfy the assumed
//! constraint, letting the neighbour speak for its own half cell tracks the real profile
//! more closely and leaves a smaller deviation to reconstruct.  Measured on an isothermal
//! atmosphere reconstructed with the (deliberately wrong) isodensity and isentropic
//! backgrounds, integrating everything from the anchor instead inflates the residual
//! hydrostatic imbalance by about a factor of three.  So the general-EOS background walks
//! the same two segments, and this is the single-segment step.
//!
//! `t` returns the temperature of the state this leaves behind, which every branch knows
//! for free: the isodensity branch gets it out of the pressure inversion, the isothermal
//! branch IS the temperature it started from, and the isentropic branch has it from the
//! last Newton iterate.  It exists so WBBackgroundStencil() can evaluate the pressure
//! channel without solving for T a second time on a state that was just constructed.
//! Handing back the PRESSURE instead would be the obvious move and is the wrong one: the
//! ideal-gas getWBq0() derives its pressure channel as (gamma-1) times the same energy it
//! reconstructs, so the general path must do the same or it stops reproducing the ideal
//! path bit for bit -- and the whole verification of this file rests on that.  A
//! temperature is safe precisely because Pressure(d,e,t) IGNORES t under a gamma law.

//! `dlntdphi` is the local logarithmic temperature gradient with respect to the
//! potential, d ln T / d Phi, measured across the stencil by WBBackgroundStencil(); it is
//! only read by the POLYTROPIC branch (wb_opt == 3) and defaults to zero, which makes
//! that branch the isothermal one.

KOKKOS_INLINE_FUNCTION
void WBAdvance(const EOS_Data &eos, const int wb_opt, const Real d_c, const Real e_c,
               const Real dphi, Real &d, Real &e, Real &t, const Real tguess = -1.0,
               const Real dlntdphi = 0.0, const Real tstart = -1.0,
               const Real t_cc = -1.0, Real *p_end = nullptr) {
  // `p_end` receives the pressure of the state left behind, from the same table
  // evaluation that produced its energy, so the caller does not evaluate it again.
  Real dum;
  // `tstart` is the temperature of the state the segment starts from and `t_cc` that of
  // the coefficient cell, when the caller KNOWS them (the anchor's cached T, the
  // interface T the previous segment handed back, the neighbour's cached T).  Either
  // non-positive means "solve for it".  With both known, the isothermal, polytropic and
  // isodensity branches do no root find at all.
  if (wb_opt == 4) {
    // ISENTROPE IN (d,T), no root find.  Along an adiabat dlnT/dln d = Gamma_3 - 1 =
    // p chi_T/(d T c_v), and hydrostatic balance with dp = p(chi_rho dln d + chi_T dlnT)
    // gives  dln d/dPhi = -d/(p Gamma_1),  Gamma_1 = chi_rho + chi_T (Gamma_3 - 1).
    // Integrated with the midpoint rule like the polytropic branch; every coefficient
    // is one table evaluation at a known (d,T).  The Newton branch (wb_opt == 2) instead
    // projects the energy onto the exact invariant h + Phi = const, at the price of a
    // temperature root find per iteration; the two agree to the walk's O(dphi^2), and
    // neither affects the balance itself, which is exact for the background either way.
    const Real t0 = (tstart > 0.0) ? tstart : eos.Temperature(d, e, tguess);
    Real e0, p0, chir, chit, cv;
    eos.ThermoAt(d, t0, e0, p0, chir, chit, cv);
    Real g3m1 = p0*chit/(d*t0*cv);
    Real k1 = -d/(p0*(chir + chit*g3m1));
    const Real dm = d*exp(0.5*k1*dphi);
    const Real tm = t0*exp(0.5*g3m1*k1*dphi);
    Real em, pm, chirm, chitm, cvm;
    eos.ThermoAt(dm, tm, em, pm, chirm, chitm, cvm);
    const Real g3m1m = pm*chitm/(dm*tm*cvm);
    const Real k2 = -dm/(pm*(chirm + chitm*g3m1m));
    d *= exp(k2*dphi);
    t = t0*exp(g3m1m*k2*dphi);
    Real pe, chire, chite;
    eos.ThermoAt(d, t, e, pe, chire, chite, dum);
    if (p_end != nullptr) {*p_end = pe;}
  } else if (wb_opt == 3) {
    // LOCAL POLYTROPE: T follows the gradient the stencil actually has,
    //     ln T(Phi) = ln T_0 + a (Phi - Phi_0),   a = dlntdphi,
    // and the density follows hydrostatic balance with p = p(d,T):
    //     dp = p chi_rho dln d + p chi_T dln T  and  dp/dPhi = -d
    //     =>  dln d/dPhi = -(d + a p chi_T)/(p chi_rho).
    // For a gamma law this is the polytrope d ~ T^n with n = -1/(a T) - 1, which contains
    // the isothermal (a = 0) and isentropic (a = -(gamma-1)/gamma) backgrounds as
    // members -- so a radiative profile, which is neither, is reproduced to one order
    // higher than either of them.  The walk stays in (d,T): the coefficients p, chi_rho
    // and chi_T at a known temperature are direct table reads, no root find, which is
    // what makes this branch cheaper than the isentropic one despite the midpoint step.
    // The segment is integrated with the midpoint rule (coefficient re-evaluated at the
    // half-way state), so the O(dphi^2) truncation of a single frozen-coefficient
    // exponential does not eat the order the closure just bought.
    const Real t0 = (tstart > 0.0) ? tstart : eos.Temperature(d, e, tguess);
    Real e0, p0, chir0, chit0, em, pm, chirm, chitm, pe, chire, chite;
    eos.ThermoAt(d, t0, e0, p0, chir0, chit0, dum);
    const Real k1 = -(d + dlntdphi*p0*chit0)/(p0*chir0);
    const Real dm = d*exp(0.5*k1*dphi);
    const Real tm = t0*exp(0.5*dlntdphi*dphi);
    eos.ThermoAt(dm, tm, em, pm, chirm, chitm, dum);
    const Real k2 = -(dm + dlntdphi*pm*chitm)/(pm*chirm);
    d *= exp(k2*dphi);
    t = t0*exp(dlntdphi*dphi);
    eos.ThermoAt(d, t, e, pe, chire, chite, dum);
    if (p_end != nullptr) {*p_end = pe;}
  } else if (wb_opt == 0) {
    // ISODENSITY: d is fixed, and dp/dPhi = -d integrates exactly over the segment.
    Real p = eos.Pressure(d, e, (tstart > 0.0) ? tstart : eos.Temperature(d, e, tguess))
             - d_c*dphi;
    e = eos.EnergyFromPressure(d, p, t);
    if (p_end != nullptr) {*p_end = eos.Pressure(d, e, t);}
  } else if (wb_opt == 1) {
    // ISOTHERMAL: dln d/dPhi = -d/(p chi_rho), frozen at the coefficient cell, is exact
    // for an ideal gas (the coefficient is 1/T) and second-order otherwise.  T is the
    // constrained variable, so it is read back off the state the segment starts from,
    // where the constraint already holds -- no separate reference has to be threaded
    // through.  For an ideal gas EnergyFromTemperature makes e scale with d, reproducing
    // the exp() branch of getWBerho exactly.
    Real t_ref = (tstart > 0.0) ? tstart : eos.Temperature(d, e, tguess);
    // p and chi_rho are both wanted at the coefficient cell, so its temperature is solved
    // for once and handed to both
    Real t_c = (t_cc > 0.0) ? t_cc : eos.Temperature(d_c, e_c, t_ref);
    Real e_cc, p_c, chir_c, chit_c;
    eos.ThermoAt(d_c, t_c, e_cc, p_c, chir_c, chit_c, dum);
    d *= exp(-d_c*dphi/(p_c*chir_c));
    Real pe, chire, chite;
    eos.ThermoAt(d, t_ref, e, pe, chire, chite, dum);
    if (p_end != nullptr) {*p_end = pe;}
    t = t_ref;   // the whole point of this branch: T is what is held fixed
  } else {
    // ISENTROPIC.  Two steps, neither of which needs an entropy function.
    //
    // (a) DENSITY, from a local power-law adiabat.  For a gamma law the exact statement is
    // that d^(gamma-1) is linear in Phi with slope -(gamma-1)/gamma * d_c^gamma/p_c, which
    // is getWBerho's isentropic branch.  Gamma_1 is precisely the exponent that plays
    // gamma's role here -- (dln p/dln d) at constant entropy -- so using Gamma_1 of the
    // COEFFICIENT CELL gives the same step for a general EOS.  This is not the discredited
    // "substitute Gamma_1 into the closed form" move: that failed because p = (gamma-1)e
    // is an ideal-gas identity applied globally, whereas Gamma_1 is used here only as a
    // local coefficient over half a cell, exactly the role the adiabat constant plays in
    // the ideal-gas expression.  It reduces to that expression identically when Gamma_1 is
    // constant, and re-anchors to the neighbouring cell in the same way.
    //
    // (b) ENERGY, from the exact invariant h + Phi = const.  The density step carries the
    // error of freezing Gamma_1; projecting the energy onto the invariant keeps that error
    // out of the pressure-gravity balance entirely.
    Real t_h = (tstart > 0.0) ? tstart : eos.Temperature(d, e, tguess);
    Real h_target = eos.Enthalpy(d, e, t_h) - dphi;
    // Gamma_1 and p are both wanted at the coefficient cell; one temperature serves both
    Real t_c = (t_cc > 0.0) ? t_cc : eos.Temperature(d_c, e_c, t_h);
    // Gamma_1 = chi_rho + p chi_T^2/(d T c_v) and p, from one evaluation
    Real e_cc, p_c, chir_c, chit_c, cv_c;
    eos.ThermoAt(d_c, t_c, e_cc, p_c, chir_c, chit_c, cv_c);
    Real g1 = chir_c + p_c*chit_c*chit_c/(d_c*t_c*cv_c);
    Real gm1 = g1 - 1.0;
    Real d_new;
    if (fabs(gm1) > 1.0e-8) {
      Real u = pow(d, gm1) - dphi*gm1/g1*pow(d_c, g1)/p_c;
      // a large enough step can drive the power-law variable through zero; the gamma -> 1
      // exponential below is the smooth limit of the same relation and stays positive
      d_new = (u > 0.0) ? pow(u, 1.0/gm1) : d*exp(-dphi*d_c/(g1*p_c));
    } else {
      d_new = d*exp(-dphi*d_c/(g1*p_c));
    }
    e = WBEnergyFromEnthalpy(eos, d_new, h_target, e*d_new/d, t, t_c, p_end);
    d = d_new;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void WBBackgroundStencil
//! \brief background states at the five points the deviation reconstruction needs:
//! the two neighbouring cell centres, the two interfaces, and the anchor cell itself.
//! Argument order matches the ideal-gas getWBerho() so the two are interchangeable at the
//! call sites.  The anchor state is returned unchanged at i, by construction.

KOKKOS_INLINE_FUNCTION
void WBBackgroundStencil(const EOS_Data &eos, const WBOption wb_option,
                         const Real rho_im1, const Real rho_i, const Real rho_ip1,
                         const Real e_im1, const Real e_i, const Real e_ip1,
                         const Real phi_im1, const Real phi_imh, const Real phi_i,
                         const Real phi_iph, const Real phi_ip1,
                         WBState &q0_im1, WBState &q0_imh, WBState &q0_i,
                         WBState &q0_iph, WBState &q0_ip1,
                         const Real t_im1c = -1.0, const Real t_ic = -1.0,
                         const Real t_ip1c = -1.0) {
  // t_*c: the cached temperatures of the three cells (WBT), or -1 to solve.
  // The four advanced points get their pressure from the temperature WBAdvance() hands
  // back, not from a fresh root find: the states have only just been constructed and
  // their temperature is already known. Only the anchor, which is not advanced, still
  // pays for a solve -- one per stencil instead of five. Under a gamma law
  // Pressure(d,e,t) does not read t at all, so this is bitwise identical there; under the
  // tabulated EOS it is also more self consistent, since t is the temperature that
  // actually produced e rather than one re-inferred from it.
  // ONE temperature inversion per stencil, at the anchor. Every other state here is
  // within a half cell of it, so it seeds every remaining inversion -- the option test,
  // all four segment advances, and the Newton loop inside the isentropic branch. Before
  // this, each of those bracketed from scratch, which on the dhj outer boundary meant
  // dozens of cold-start root finds per ghost cell per stage.
  Real t_i = (t_ic > 0.0) ? t_ic : eos.Temperature(rho_i, e_i);
  Real dum_cv;

  int wb_opt = WBOptionNumber(eos, wb_option, rho_im1, e_im1, rho_ip1, e_ip1, rho_i, e_i,
                              t_i, t_im1c, t_ip1c);
  // the local polytrope's gradient, d ln T / d Phi across the whole stencil.  Two more
  // inversions, seeded by the anchor's temperature.  Zero when the stencil does not
  // span any potential (a tangential sweep), which reduces the branch to isothermal.
  Real dlntdphi = 0.0;
  if (wb_opt == 3) {
    const Real dphis = phi_ip1 - phi_im1;
    if (fabs(dphis) > 0.0) {
      const Real t_m = (t_im1c > 0.0) ? t_im1c : eos.Temperature(rho_im1, e_im1, t_i);
      const Real t_p = (t_ip1c > 0.0) ? t_ip1c : eos.Temperature(rho_ip1, e_ip1, t_i);
      dlntdphi = log(t_p/t_m)/dphis;
    }
  }

  q0_i.d = rho_i;
  q0_i.e = e_i;

  if (wb_opt == 3) {
    // POLYTROPIC FAST PATH: five table evaluations per stencil, one per state.  T(Phi) is
    // prescribed (ln T linear in Phi with the stencil's gradient a), so every state's
    // (d,T) is known once its density is; the density steps use the coefficient
    //     k = -(d + a p chi_T)/(p chi_rho)
    // from the evaluation at each segment's end, which also delivers that state's e and
    // p: anchor (p_i and the anchor coefficient), two interfaces, two neighbours.  A
    // third of the midpoint version's table traffic at the same order (see `seg`).
    Real ei_t, p_i, chir, chit;
    eos.ThermoAt(rho_i, t_i, ei_t, p_i, chir, chit, dum_cv);
    q0_i.p = p_i;
    const Real a = dlntdphi;
    // one segment: predictor with the start coefficient k0, ONE table evaluation at the
    // predicted end point (its e, p and chi's), trapezoidal corrector on ln d with the
    // end coefficient, and e, p moved to the corrected density with the evaluated
    // chi_rho (p ~ d^chi_rho at fixed T) and e ~ d (exact for a gamma law, O(dphi^2)
    // otherwise on a correction that is itself O(dphi^2)).  The end evaluation's chi's
    // seed the next segment.  Second order per segment, five evaluations per stencil.
    auto seg = [&](const Real d0, const Real t0, const Real k0, const Real dph,
                   Real &d1, Real &t1, Real &e1, Real &p1, Real &k1) {
      const Real dpred = d0*exp(k0*dph);
      t1 = t0*exp(a*dph);
      Real ep, pp, cr, ct;
      eos.ThermoAt(dpred, t1, ep, pp, cr, ct, dum_cv);
      const Real kend = -(dpred + a*pp*ct)/(pp*cr);
      d1 = d0*exp(0.5*(k0 + kend)*dph);
      // r - 1 is O(dphi^2), so first order in it is enough (and saves a pow)
      const Real r = d1/dpred;
      p1 = pp*(1.0 + cr*(r - 1.0));
      e1 = ep*r;
      k1 = -(d1 + a*p1*ct)/(p1*cr);
    };
    const Real k0 = -(rho_i + a*p_i*chit)/(p_i*chir);
    Real d1, t1, e1, p1, k1, d2, t2, e2, p2, k2;
    seg(rho_i, t_i, k0, phi_imh - phi_i, d1, t1, e1, p1, k1);
    q0_imh.d = d1; q0_imh.e = e1; q0_imh.p = p1;
    seg(d1, t1, k1, phi_im1 - phi_imh, d2, t2, e2, p2, k2);
    q0_im1.d = d2; q0_im1.e = e2; q0_im1.p = p2;
    seg(rho_i, t_i, k0, phi_iph - phi_i, d1, t1, e1, p1, k1);
    q0_iph.d = d1; q0_iph.e = e1; q0_iph.p = p1;
    seg(d1, t1, k1, phi_ip1 - phi_iph, d2, t2, e2, p2, k2);
    q0_ip1.d = d2; q0_ip1.e = e2; q0_ip1.p = p2;
    return;
  }
  q0_i.p = eos.Pressure(rho_i, e_i, t_i);

  // inner half cells, from the anchor, using the anchor's own coefficient
  Real dm = rho_i, em = e_i, tm = -1.0, pm = 0.0;
  WBAdvance(eos, wb_opt, rho_i, e_i, phi_imh - phi_i, dm, em, tm, t_i, dlntdphi,
            t_i, t_i, &pm);
  q0_imh.d = dm; q0_imh.e = em; q0_imh.p = pm;

  Real dp = rho_i, ep = e_i, tp = -1.0, pp = 0.0;
  WBAdvance(eos, wb_opt, rho_i, e_i, phi_iph - phi_i, dp, ep, tp, t_i, dlntdphi,
            t_i, t_i, &pp);
  q0_iph.d = dp; q0_iph.e = ep; q0_iph.p = pp;

  // outer half cells, continuing from the interfaces, now with the neighbour speaking for
  // its own half cell -- the same two-segment walk the ideal-gas closed forms perform
  // these two continue FROM the interfaces, whose temperatures the calls above just
  // handed back in tm/tp -- a closer guess still than the anchor's
  WBAdvance(eos, wb_opt, rho_im1, e_im1, phi_im1 - phi_imh, dm, em, tm, tm, dlntdphi,
            tm, t_im1c, &pm);
  q0_im1.d = dm; q0_im1.e = em; q0_im1.p = pm;

  WBAdvance(eos, wb_opt, rho_ip1, e_ip1, phi_ip1 - phi_iph, dp, ep, tp, tp, dlntdphi,
            tp, t_ip1c, &pp);
  q0_ip1.d = dp; q0_ip1.e = ep; q0_ip1.p = pp;
  return;
}

#endif // UTILS_WB_BACKGROUND_HPP_
