#ifndef RAD_M1_RAD_M1_OPACITY_HPP_
#define RAD_M1_RAD_M1_OPACITY_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_opacity.hpp
//! \brief the analytic opacity laws of <rad_m1>, as device-inline functions.
//!
//! Everything here returns an opacity PER UNIT MASS in code units; the caller multiplies
//! by rho and stores rho*kappa (an inverse length) in RadiationM1::opac.  Milestone 1b
//! has no table opacity: `kappa_F` from the existing Rosseland table (design sect. 4) is
//! milestone 2 work.

#include <math.h>
#include "athena.hpp"
#include "rad_m1/rad_m1.hpp"

namespace radm1 {

//----------------------------------------------------------------------------------------
//! \fn M1UserOpacity
//! \brief THE USER HOOK of <rad_m1>/opacity = user.  It is device-callable and is
//! deliberately inline rather than a function pointer, because a host function pointer
//! cannot be dereferenced inside a device lambda (see the CLAUDE.md kernel rules).  A
//! problem that wants its own opacity law edits this one function and rebuilds.
//!
//! The default body returns the constant opacities, so `opacity = user` is a working
//! (if trivial) configuration out of the box and the branch is always compiled.
//!
//! \param[in]  d, t        density and temperature of the cell, code units
//! \param[in]  kp,ke,kf,ks the <rad_m1> constants, as defaults
//! \param[out] op,oe,of,os the four opacities per unit mass

KOKKOS_INLINE_FUNCTION
void M1UserOpacity(const Real d, const Real t, const Real kp, const Real ke,
                   const Real kf, const Real ks,
                   Real &op, Real &oe, Real &of, Real &os) {
  // TODO(@user): replace with the problem's own kappa(rho,T)
  (void) d;
  (void) t;
  op = kp;
  oe = ke;
  of = kf;
  os = ks;
}

//----------------------------------------------------------------------------------------
//! \fn M1Opacities
//! \brief the four opacities per unit mass at (d, t) under the law `otype`:
//!   const     kappa = kappa0
//!   powerlaw  kappa = kappa0 (d/rho_ref)^a (t/t_ref)^b,  ONE pair of exponents for all
//!             four opacities (design sect. 4 wants no more than that in stage 1)
//!   user      M1UserOpacity above

KOKKOS_INLINE_FUNCTION
void M1Opacities(const int otype, const Real d, const Real t,
                 const Real kp, const Real ke, const Real kf, const Real ks,
                 const Real rho_ref, const Real t_ref, const Real aa, const Real bb,
                 Real &op, Real &oe, Real &of, Real &os) {
  if (otype == M1_OPAC_POWERLAW) {
    Real fd = (aa == 0.0) ? 1.0 : pow(fmax(d, 0.0)/rho_ref, aa);
    Real ft = (bb == 0.0) ? 1.0 : pow(fmax(t, 0.0)/t_ref, bb);
    Real sc = fd*ft;
    op = kp*sc;
    oe = ke*sc;
    of = kf*sc;
    os = ks*sc;
  } else if (otype == M1_OPAC_USER) {
    M1UserOpacity(d, t, kp, ke, kf, ks, op, oe, of, os);
  } else {
    op = kp;
    oe = ke;
    of = kf;
    os = ks;
  }
}

} // namespace radm1
#endif // RAD_M1_RAD_M1_OPACITY_HPP_
