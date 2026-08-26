#ifndef COORDINATES_GRID_STRETCH_HPP_
#define COORDINATES_GRID_STRETCH_HPP_
//========================================================================================
// Athena++K astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file grid_stretch.hpp
//! \brief coordinate stretching maps for the spherical-polar grid.
//!
//! Kept in their own header, free of any dependency beyond athena.hpp, because both
//! Coordinates and problem generators need them: Coordinates stretches x1v/x1f/dx1, while
//! a problem generator that rebuilds radii with CellCenterX/LeftEdgeX to lay down an
//! initial condition must apply the identical map or the two end up on different grids.

#include <cmath>

#include "athena.hpp"

// Number of coefficients in the polynomial radial stretch (see StretchRPoly):
//     u(xi) = xi + sum_{k=1}^{NSTRETCH_R_POLY} c_k xi^k (1-xi).
#define NSTRETCH_R_POLY 4


KOKKOS_INLINE_FUNCTION
void StretchR(const Real a, const Real r0, const Real r1, Real &r) {
  Real xi = (r-r0)/(r1-r0);
  Real denom = 1.0 - exp(-a);
  r = r0 + (r1 - r0)*(1.0 - exp(-a*xi))/denom;
//    r = r0*pow(r1/r0,xi);
}
//--------------------------------------------------------------------------------------
//! \fn StretchRPoly
//! \brief polynomial radial grid stretch, an alternative to the exponential StretchR.
//!
//! Maps the uniform coordinate onto
//!     r = r0 + (r1 - r0) * u(xi),   xi = (r - r0)/(r1 - r0),
//!     u(xi) = xi + sum_{k=1}^{NSTRETCH_R_POLY} c_k xi^k (1 - xi),
//! so that u(0) = 0 and u(1) = 1 for ANY coefficients: the domain end points are fixed
//! and only the interior distribution moves. c = 0 recovers the uniform grid.
//!
//! WHY A POLYNOMIAL AND NOT THE EXPONENTIAL. StretchR is monotonic in cell width, so it
//! can only coarsen (or only refine) outwards. A stratified atmosphere does not want
//! that: the pressure scale height H is large at depth, collapses across a dissociation
//! front, then grows steadily aloft, so the cell width that resolves H equally
//! everywhere is NON-MONOTONIC. This family can represent that.
//!
//! HOW TO CHOOSE THE COEFFICIENTS. Constant cells-per-scale-height means dr ~ H(r),
//! which is exactly uniform spacing in ln p. Take a relaxed snapshot, form the
//! cumulative scale-height coordinate N_H(r) = int dr'/H(r'), and least-squares fit
//! u_target(xi) = (r(xi N_H^tot) - r0)/(r1 - r0) in the basis xi^k (1-xi).
//!
//! The mapping must be strictly increasing; Mesh's constructor samples du/dxi over
//! [0,1] and fatals if it is not, since a fold-over gives negative cell widths.
KOKKOS_INLINE_FUNCTION
void StretchRPoly(const Real *c, const Real r0, const Real r1, Real &r) {
  Real xi = (r-r0)/(r1-r0);
  Real u = xi;
  Real xik = xi;                      // xi^k, built up as k increases
  for (int k=1; k<=NSTRETCH_R_POLY; ++k) {
    u += c[k-1]*xik*(1.0-xi);
    xik *= xi;
  }
  r = r0 + (r1-r0)*u;
}
KOKKOS_INLINE_FUNCTION
void StretchTheta(const Real a, Real &t) {
  Real xi = t/M_PI;
  t = M_PI/2.0*(1.0+sinh(a*(2.0*xi-1.0))/sinh(a));
}
  


#endif // COORDINATES_GRID_STRETCH_HPP_
