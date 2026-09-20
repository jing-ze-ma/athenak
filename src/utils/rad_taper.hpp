#ifndef UTILS_RAD_TAPER_HPP_
#define UTILS_RAD_TAPER_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_taper.hpp
//! \brief the density weight w(rho) that switches the LTE radiation terms of the general
//! EOS off in an optically thin layer, shared by the EOS and by the two-stream solver.
//!
//! WHY. With <eos>/eos_radiation the table EOS adds aT^4 and aT^4/3 to e and p at EVERY
//! density. Deep down that is right: the radiation field is a Planckian in equilibrium
//! with the gas, it is advected with it, and -grad(aT^4/3) is the radiative force. At the
//! top of a stellar-surface box it is wrong twice over. The real radiation field there is
//! the two-stream's, which is anisotropic, streams at c, and is not a property of the
//! local gas at all; and its energy is not in the cell's internal energy, so counting
//! aT^4 there inflates e, p and the sound speed without bound. MEASURED on the relaxed
//! B-star column (bench/bstar_fecz): at tau = 1e-2 the tabulated Prad/Pgas is 8.6 and the
//! radiation energy is 17x the gas thermal energy, in the most violent layer of the box.
//!
//! WHAT. w(rho) is a smoothstep in log10(rho): 1 at and below rho_hi... strictly, 1 for
//! rho >= rho_hi (dense, LTE) and 0 for rho <= rho_lo (thin, the two-stream owns the
//! radiation), with the C1 cubic 3s^2 - 2s^3 in between. The EOS multiplies its radiation
//! terms by w; the momentum source in the two-stream supplies the force the pressure
//! gradient no longer carries, weighted by (1 - w). The two are complementary by
//! construction, so the total force is -grad Prad in the deep limit and kappa rho F/c in
//! the thin limit, with no gap at either end of the ramp.
//!
//! The DERIVATIVE is returned alongside, because both users need it: the EOS to keep
//! chi_rho the true derivative of the p it returns, and the source term to cancel the
//! P_rad grad w piece of that gradient, which is an artefact of the taper and not a
//! force. dw/dx is with respect to x = log10(rho[cgs]), the EOS table's own abscissa.

#include "athena.hpp"

namespace rad_taper {

//----------------------------------------------------------------------------------------
//! \fn void Weight
//! \brief the smoothstep w and its derivative dw/dx at x = log10(rho[cgs]).
//! `xlo` = log10(rho_lo) and `xhi` = log10(rho_hi) with xlo < xhi; outside [xlo,xhi] the
//! weight is constant and the derivative is exactly zero, so a cell that is entirely deep
//! or entirely thin costs nothing and picks up no spurious gradient term.

KOKKOS_INLINE_FUNCTION
void Weight(const Real x, const Real xlo, const Real xhi, Real &w, Real &dwdx) {
  if (x <= xlo) {
    w = 0.0;
    dwdx = 0.0;
  } else if (x >= xhi) {
    w = 1.0;
    dwdx = 0.0;
  } else {
    const Real idx = 1.0/(xhi - xlo);
    const Real s = (x - xlo)*idx;
    w = s*s*(3.0 - 2.0*s);
    dwdx = 6.0*s*(1.0 - s)*idx;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void WeightGated
//! \brief the taper weight with a TEMPERATURE GATE: w = max(w_rho(x), w_T(y)).
//!
//! WHY (he4-presn, tests_r7/tests_r8).  w(rho) alone says "thin" of any cell whose own
//! density is low, and that is not the statement the taper is meant to make.  The taper
//! marks where the two-stream takes the radiation over from the gas, which is an OPTICAL
//! DEPTH statement.  On a star with a FeCZ density inversion the two part company: the
//! He-star convection zone sits only 1.2x above `eos_rad_rho_hi`, so a cell that a
//! seam-row flux error thins by half a decade falls INTO the window at tau ~ 150 and is
//! handed the optically-thin surface force, kappa F/(c g) = 1.23, in the middle of the
//! convection zone.  It then evacuates further; that loop killed every 3-D arm of
//! tests_r6 at 1.04-1.08 turnover and the 1-D column at 1.19 (11 of 11 collapses on a
//! panel-seam row), and `rt_rad_force = false` was the only switch that removed it.
//!
//! The gate closes that loop with the one monotone stratification variable the EOS has in
//! hand at EVERY evaluation point: the temperature.  `y` = log10(T[K]); w_T = 1 at and
//! above `yhi` = log10 T(tau_hi), 0 at and below `ylo` = log10 T(tau_lo), the same cubic
//! in between.  A cell inside the star is hot, so it keeps w = 1: the radiation stays in
//! the EOS and the force is -grad P_rad through the hydro -- WHATEVER its density does.
//!
//! WHY TEMPERATURE AND NOT RADIUS OR tau.  Both would be better statements of "optically
//! thick", and neither is reachable: the EOS is evaluated per cell with no cell context
//! at all (see EOSTable -- it takes rho and T and nothing else), and the radiation terms
//! are read by the well-balanced background, conduction, the reconstruction, the Riemann
//! solvers, RaiseVel and the problem generators as well as by ConsToPrim, ~150 call sites
//! over 28 files.  Threading a coordinate through all of them is the only way to make a
//! radius or tau gate CONSISTENT, and an inconsistent gate is worse than none.  On a
//! star T(r) is monotone where rho(r) is not -- which is the defect itself -- and the two
//! thresholds are read off the SAME two optical depths a radius gate would have used.
//!
//! Exactly one of the two derivatives is non-zero, because max() picks one branch: the
//! caller that needs grad w (the two-stream force) therefore gets d w/d log10 rho in
//! `dwdx` and d w/d log10 T in `dwdy` and can form both gradient terms unconditionally.
//! `gate = false` returns Weight()'s own answer with dwdy = 0, bit for bit.

KOKKOS_INLINE_FUNCTION
void WeightGated(const Real x, const Real xlo, const Real xhi,
                 const Real y, const Real ylo, const Real yhi, const bool gate,
                 Real &w, Real &dwdx, Real &dwdy) {
  Weight(x, xlo, xhi, w, dwdx);
  dwdy = 0.0;
  if (gate) {
    Real wt, dwt;
    Weight(y, ylo, yhi, wt, dwt);
    if (wt > w) {
      w = wt;
      dwdx = 0.0;
      dwdy = dwt;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Real WeightOnly
//! \brief the weight alone, for callers with no use for the derivative.

KOKKOS_INLINE_FUNCTION
Real WeightOnly(const Real x, const Real xlo, const Real xhi) {
  Real w, dwdx;
  Weight(x, xlo, xhi, w, dwdx);
  return w;
}

}  // namespace rad_taper

#endif  // UTILS_RAD_TAPER_HPP_
