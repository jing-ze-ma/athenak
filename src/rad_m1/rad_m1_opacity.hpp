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
#include "diffusion/conduction.hpp"   // RosselandTable: ONE tabulated-opacity lookup

namespace radm1 {

//----------------------------------------------------------------------------------------
//! \fn M1TableOpacities
//! \brief <rad_m1>/opacity = table: the stellar Rosseland and Planck means, looked up in
//! KELVIN and g/cm^3 on the shared (log10 T, log10 rho) grid of M1OpacTab.
//!
//! THE CONSISTENT SET (stage-2 plan sect. 4).  The tabulated Rosseland total ALREADY
//! contains electron scattering, so it must not be added a second time: the transport
//! opacity is kappa_F = kappa_R with kappa_s = 0.  The Planck table is absorption only
//! and serves both the emission and the absorption mean, kappa_P = kappa_E = kappa_Pl.
//! On the He column the two differ by 15x to 527x, which is why the design refuses
//! planck_from_rosseland for any quoted result.
//!
//! \param[in]  tab   the two tables and their unit conversions
//! \param[in]  d, t  density and temperature of the cell, CODE units
//! \param[out] op,oe,of,os  the four opacities per unit mass, CODE units

KOKKOS_INLINE_FUNCTION
void M1TableOpacities(const M1OpacTab &tab, const Real d, const Real t,
                      Real &op, Real &oe, Real &of, Real &os) {
  const Real tk = t*tab.tunit;
  const Real dc = d*tab.dunit;
  // m1-fast4: both tables share the grid, so the logarithms, the two bisections and the
  // weights are made ONCE; the arithmetic is RosselandTable's, term for term (bitwise)
  const Real x = log10(tk), y = log10(dc);
  const int nT = tab.nT, nD = tab.nD;
  int i = 0, j = 0;
  Real fx = 0.0, fy = 0.0;
  if (!(x > tab.lT(0))) {
    i = 0; fx = 0.0;
  } else if (x >= tab.lT(nT-1)) {
    i = nT-2; fx = 1.0;
  } else {
    int hi = nT-2;
    while (i < hi) {
      const int mid = (i + hi + 1)/2;
      if (tab.lT(mid) <= x) {i = mid;} else {hi = mid - 1;}
    }
    fx = (x - tab.lT(i))/(tab.lT(i+1) - tab.lT(i));
  }
  if (!(y > tab.lD(0))) {
    j = 0; fy = 0.0;
  } else if (y >= tab.lD(nD-1)) {
    j = nD-2; fy = 1.0;
  } else {
    int hj = nD-2;
    while (j < hj) {
      const int mid = (j + hj + 1)/2;
      if (tab.lD(mid) <= y) {j = mid;} else {hj = mid - 1;}
    }
    fy = (y - tab.lD(j))/(tab.lD(j+1) - tab.lD(j));
  }
  const Real lr = (1.0-fx)*((1.0-fy)*tab.kr(i,j) + fy*tab.kr(i,j+1))
                +      fx *((1.0-fy)*tab.kr(i+1,j) + fy*tab.kr(i+1,j+1));
  const Real lp = (1.0-fx)*((1.0-fy)*tab.kp(i,j) + fy*tab.kp(i,j+1))
                +      fx *((1.0-fy)*tab.kp(i+1,j) + fy*tab.kp(i+1,j+1));
  const Real kr = pow(10.0, lr);
  const Real kp = pow(10.0, lp);
  op = kp*tab.kunit;
  oe = op;
  of = kr*tab.kunit;
  os = 0.0;
}

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
