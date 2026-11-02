#ifndef RAD_M1_RAD_M1_VGDTHETA_HPP_
#define RAD_M1_RAD_M1_VGDTHETA_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_vgdtheta.hpp
//! \brief accsp-1102: the device helpers of vet_gd on a STRETCHED theta and on the FULL
//! sphere (RadiationM1::vgd_sth / vgd_pol, rad_m1_vetgd_sph2.inc).
//!
//! The band rows of a block carry the node table tn(m, jb) (the face midpoint of band row
//! jb in the EXTENDED chart, strictly increasing).  An upwind point's fractional band row
//! is linear in theta between the two nodes that bracket it (the uniform index formula
//! j + (thu - th)/dx2 on a uniform grid, up to round-off); outside the table it is
//! extrapolated with the end interval (then clamped by the caller as before).

#include <cmath>

#include "athena.hpp"

// fractional band row of theta thu; jb = a starting guess (the cell's own row: upwind
// points lie a few rows away, so the walk is short)
template <typename T>
KOKKOS_INLINE_FUNCTION
Real VgdThFj(const T &tn, const int m, const int jb, const Real thu, const int c2w) {
  int j0 = (jb < 0) ? 0 : ((jb > c2w - 2) ? (c2w - 2) : jb);
  while (j0 < c2w - 2 && tn(m,j0+1) <= thu) {++j0;}
  while (j0 > 0 && tn(m,j0) > thu) {--j0;}
  return j0 + (thu - tn(m,j0))/(tn(m,j0+1) - tn(m,j0));
}

// the full sphere: an upwind point BEHIND the polar axis (|dphi| > pi/2 from the cell)
// is taken in the extended chart of the cell's hemisphere: (theta, phi) -> (-theta,
// phi + pi) north, (2 pi - theta, phi + pi) south; the band rows beyond the pole hold
// the pole partner's reversed rows at exactly those positions
KOKKOS_INLINE_FUNCTION
void VgdPoleChart(const bool pol, const Real th, Real &thu, Real &dph) {
  if (!pol || !(fabs(dph) > 0.5*M_PI)) {return;}
  thu = (th < 0.5*M_PI) ? -thu : (2.0*M_PI - thu);
  dph += (dph > 0.0) ? -M_PI : M_PI;
}

// the band row a pole-slot read takes in its on-rank partner: rows reversed about the
// pole (north: band row jj < w -> interior row 2 w - 1 - jj; south: jj >= w + nx2 ->
// 2 (w + nx2) - 1 - jj)
KOKKOS_INLINE_FUNCTION
int VgdPoleRow(const int jj, const int w, const int nx2) {
  return (jj < w) ? (2*w - 1 - jj) : (2*(w + nx2) - 1 - jj);
}

#endif  // RAD_M1_RAD_M1_VGDTHETA_HPP_
