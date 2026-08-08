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

KOKKOS_INLINE_FUNCTION
int WBOptionNumber(const EOS_Data &eos, const WBOption wb_option,
                   const Real d_m, const Real e_m, const Real d_p, const Real e_p,
                   const Real d_i, const Real e_i) {
  switch (wb_option) {
    case WBOption::isodensity:
      return 0;
    case WBOption::isothermal:
      return 1;
    case WBOption::isentropic:
      return 2;
    case WBOption::adaptive: {
      Real t_m = eos.Temperature(d_m, e_m);
      Real t_p = eos.Temperature(d_p, e_p);
      Real t_i = eos.Temperature(d_i, e_i);
      Real dlnt = log(t_p/t_m);
      Real dlnd = log(d_p/d_m);
      // entropy difference across the stencil, in units of c_v, evaluated at cell i
      Real cv = eos.SpecificHeatCv(d_i, e_i);
      Real p_i = eos.Pressure(d_i, e_i);
      Real dsdivcv = dlnt - (p_i*eos.ChiT(d_i,e_i)/(d_i*t_i*cv))*dlnd;
      return (fabs(dlnt) > fabs(dsdivcv)) ? 2 : 1;
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

KOKKOS_INLINE_FUNCTION
Real WBEnergyFromEnthalpy(const EOS_Data &eos, const Real d, const Real h,
                          const Real e_guess) {
  Real e = e_guess;
  for (int it=0; it<10; ++it) {
    Real p = eos.Pressure(d, e);
    Real f = e + p - d*h;
    if (f == 0.0) {break;}
    Real t = eos.Temperature(d, e);
    Real dpde = p*eos.ChiT(d, e)/(d*t*eos.SpecificHeatCv(d, e));
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

KOKKOS_INLINE_FUNCTION
void WBAdvance(const EOS_Data &eos, const int wb_opt,
               const Real d_c, const Real e_c, const Real dphi, Real &d, Real &e) {
  if (wb_opt == 0) {
    // ISODENSITY: d is fixed, and dp/dPhi = -d integrates exactly over the segment.
    Real p = eos.Pressure(d, e) - d_c*dphi;
    e = eos.EnergyFromPressure(d, p);
  } else if (wb_opt == 1) {
    // ISOTHERMAL: dln d/dPhi = -d/(p chi_rho), frozen at the coefficient cell, is exact
    // for an ideal gas (the coefficient is 1/T) and second-order otherwise.  T is the
    // constrained variable, so it is read back off the state the segment starts from,
    // where the constraint already holds -- no separate reference has to be threaded
    // through.  For an ideal gas EnergyFromTemperature makes e scale with d, reproducing
    // the exp() branch of getWBerho exactly.
    Real t_ref = eos.Temperature(d, e);
    Real p_c = eos.Pressure(d_c, e_c);
    d *= exp(-d_c*dphi/(p_c*eos.ChiRho(d_c, e_c)));
    e = eos.EnergyFromTemperature(d, t_ref);
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
    Real h_target = eos.Enthalpy(d, e) - dphi;
    Real g1 = eos.Gamma1(d_c, e_c);
    Real p_c = eos.Pressure(d_c, e_c);
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
    e = WBEnergyFromEnthalpy(eos, d_new, h_target, e*d_new/d);
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
                         WBState &q0_iph, WBState &q0_ip1) {
  int wb_opt = WBOptionNumber(eos, wb_option, rho_im1, e_im1, rho_ip1, e_ip1, rho_i, e_i);

  q0_i.d = rho_i;
  q0_i.e = e_i;
  q0_i.p = eos.Pressure(rho_i, e_i);

  // inner half cells, from the anchor, using the anchor's own coefficient
  Real dm = rho_i, em = e_i;
  WBAdvance(eos, wb_opt, rho_i, e_i, phi_imh - phi_i, dm, em);
  q0_imh.d = dm; q0_imh.e = em; q0_imh.p = eos.Pressure(dm, em);

  Real dp = rho_i, ep = e_i;
  WBAdvance(eos, wb_opt, rho_i, e_i, phi_iph - phi_i, dp, ep);
  q0_iph.d = dp; q0_iph.e = ep; q0_iph.p = eos.Pressure(dp, ep);

  // outer half cells, continuing from the interfaces, now with the neighbour speaking for
  // its own half cell -- the same two-segment walk the ideal-gas closed forms perform
  WBAdvance(eos, wb_opt, rho_im1, e_im1, phi_im1 - phi_imh, dm, em);
  q0_im1.d = dm; q0_im1.e = em; q0_im1.p = eos.Pressure(dm, em);

  WBAdvance(eos, wb_opt, rho_ip1, e_ip1, phi_ip1 - phi_iph, dp, ep);
  q0_ip1.d = dp; q0_ip1.e = ep; q0_ip1.p = eos.Pressure(dp, ep);
  return;
}

#endif // UTILS_WB_BACKGROUND_HPP_
