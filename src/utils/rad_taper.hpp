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
