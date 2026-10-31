#ifndef UTILS_CB_SENSOR_HPP_
#define UTILS_CB_SENSOR_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cb_sensor.hpp
//! \brief Odd-even (2-cell, "checkerboard") sensor for the low-Mach Riemann solvers
//! LHLLC (hydro) and LHLLD (MHD); keys <hydro>/lhllc_cb_sensor, <mhd>/lhlld_cb_sensor.
//!
//! Both solvers multiply the velocity-jump term rho c (u_L - u_R) of the contact
//! pressure by phi = chi (2 - chi), chi ~ M.  That term is the only part of the flux that
//! damps a 2-cell COMPRESSIVE mode of the normal velocity (the role Rhie-Chow
//! interpolation plays in collocated incompressible codes): the mass flux of a pure
//! checkerboard u_n = (-1)^i du is zero on every face, so the pressure never sees it.
//! With phi = O(M) the damping rate of that mode drops from ~c/dx to ~|u|/dx and the
//! mode survives.  The low-Mach reduction is meant for RESOLVED, nearly incompressible
//! flow, which has no 2dx content; a 2dx compressive mode is an unresolved acoustic wave
//! and should get the full acoustic dissipation.
//!
//! The sensor takes the face-normal mass flux q = rho u_n at the four cell centres
//! c-2..c+1 along the face normal (the face sits between c-1 and c).  With the
//! consecutive differences d1 = q(c-1)-q(c-2), d2 = q(c)-q(c-1), d3 = q(c+1)-q(c),
//!   J = (2 d2 - d1 - d3)/4
//! is the projection of the face difference d2 on the 2dx mode: J = d2 exactly for a pure
//! checkerboard (also on top of any linear q), and J = -dx^3 q'''/4 + ... for a smooth
//! q.  s = clamp(J/d2, 0, 1) is the fraction of the face jump that is 2dx (raised to
//! the power p = the key's value; p = 1 is the plain fraction, p > 1 spares faces that
//! are only partly 2dx, a pure checkerboard keeps s = 1 for any p), and the solver uses
//!   phi_eff = phi + (1 - phi) s,
//! so the face dissipation always lies between the low-Mach solver (s = 0) and plain
//! HLLC/HLLD (s = 1) and keeps their sign.  q is the MASS flux, not u_n, so a resolved
//! stratified flow obeying the anelastic constraint div(rho u) ~ 0 reads as smooth even
//! where div u = u_r/H_rho is not small.
//!
//! The floor is returned as chi_cb with phi_eff = chi_cb (2 - chi_cb), i.e.
//! chi_cb = 1 - (1 - chi) sqrt(1 - s), and is 0 when s = 0, so that
//! chi = max(chi, chi_cb) leaves the solver bitwise unchanged on faces the sensor does
//! not flag.

#include <cmath>

#include "athena.hpp"

KOKKOS_INLINE_FUNCTION
Real CBSensorChiFloor(const DvceArray5D<Real> &w0, const int m, const int k, const int j,
                      const int i, const int ivx, const Real chi, const Real gain) {
  const int cdi = (ivx == IVX) ? 1 : 0;
  const int cdj = (ivx == IVY) ? 1 : 0;
  const int cdk = (ivx == IVZ) ? 1 : 0;
  // stencil clamped to the array: FOFC-widened faces reach one cell beyond the ghosts
  const int nn = (cdi ? static_cast<int>(w0.extent(4)) :
                 (cdj ? static_cast<int>(w0.extent(3)) :
                        static_cast<int>(w0.extent(2)))) - 1;
  const int c0 = cdi*i + cdj*j + cdk*k;
  Real qs[4];
  for (int o = 0; o < 4; ++o) {
    int c = c0 - 2 + o;
    c = (c < 0) ? 0 : ((c > nn) ? nn : c);
    const int kk = cdk ? c : k, jj = cdj ? c : j, ii = cdi ? c : i;
    qs[o] = w0(m,IDN,kk,jj,ii)*w0(m,ivx,kk,jj,ii);
  }
  const Real d2 = qs[2] - qs[1];
  const Real jodd = 0.25*(2.0*d2 - (qs[1] - qs[0]) - (qs[3] - qs[2]));
  Real s = (jodd*d2 > 0.0) ? fmin(1.0, jodd/d2) : 0.0;
  if (gain != 1.0 && s > 0.0) s = pow(s, gain);
  return (s > 0.0) ? 1.0 - (1.0 - chi)*sqrt(1.0 - s) : 0.0;
}

#endif // UTILS_CB_SENSOR_HPP_
