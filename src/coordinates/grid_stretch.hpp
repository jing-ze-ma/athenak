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
// Inputs give mesh/f_stretch_r_c1..c8; missing ones are zero, which adds exact zeros
// (a 4-coefficient input gives the same grid bit for bit as when this was 4).
// Storage size of Mesh::fStretchRPoly: NSTRETCH_R_PCOEF polynomial coefficients, then
// NSTRETCH_R_BUMP local "bumps" of (amplitude, centre, width) each, then one "plateau"
// of (amplitude, xa, xb, width) (see StretchRPoly).
#define NSTRETCH_R_PCOEF 8
#define NSTRETCH_R_BUMP 2
#define NSTRETCH_R_PLAT 4
#define NSTRETCH_R_POLY (NSTRETCH_R_PCOEF + 3*NSTRETCH_R_BUMP + NSTRETCH_R_PLAT)

//! ln cosh(x) without overflow: |x| + ln(1 + exp(-2|x|)) - ln 2
KOKKOS_INLINE_FUNCTION
Real LogCoshStable(const Real x) {
  const Real ax = fabs(x);
  return ax + log1p(exp(-2.0*ax)) - 0.69314718055994530942;
}

//! the plateau term of StretchRPoly (mesh/f_stretch_r_p_amp, _xa, _xb, _w):
//!   du/dxi += a [g(xi) - G(1)],  g = (tanh((xi - xa)/w) - tanh((xi - xb)/w))/2,
//!   u      += a [G(xi) - xi G(1)], G(xi) = int_0^xi g (closed form in ln cosh),
//! so u(0) = 0 and u(1) = 1 stay fixed.  g is ~1 on [xa, xb] and ~0 outside, with tanh
//! edges of half-width w: du/dxi is FLAT inside and outside (the ratio of the two cell
//! widths is 1 + a/(1 - a G(1))), unlike a sech^2 bump or a polynomial.  a < 0 = finer
//! cells on [xa, xb].
KOKKOS_INLINE_FUNCTION
Real StretchRPlateauG(const Real xa, const Real xb, const Real w, const Real xi) {
  return 0.5*w*(LogCoshStable((xi - xa)/w) - LogCoshStable(-xa/w)
                - LogCoshStable((xi - xb)/w) + LogCoshStable(-xb/w));
}


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
  for (int k=1; k<=NSTRETCH_R_PCOEF; ++k) {
    u += c[k-1]*xik*(1.0-xi);
    xik *= xi;
  }
  // optional local bumps (mesh/f_stretch_r_b<n>_amp, _x, _w; skipped when amp = 0, so a
  // grid without them is bit for bit the plain polynomial): with a < 0 the cells near
  // xi = xb are narrower, du/dxi changes by a sech^2((xi - xb)/w) (minus the small
  // constant a w [tanh((1-xb)/w) + tanh(xb/w)] that keeps u(0) = 0 and u(1) = 1)
  for (int b=0; b<NSTRETCH_R_BUMP; ++b) {
    const Real a = c[NSTRETCH_R_PCOEF + 3*b];
    if (a == 0.0) continue;
    const Real xb = c[NSTRETCH_R_PCOEF + 3*b + 1], w = c[NSTRETCH_R_PCOEF + 3*b + 2];
    u += a*w*(tanh((xi - xb)/w) - (1.0 - xi)*tanh(-xb/w) - xi*tanh((1.0 - xb)/w));
  }
  // optional plateau (mesh/f_stretch_r_p_amp, _xa, _xb, _w; skipped when amp = 0, so a
  // grid without it is bit for bit the one above): see StretchRPlateauG
  {
    const Real *cp = c + NSTRETCH_R_PCOEF + 3*NSTRETCH_R_BUMP;
    if (cp[0] != 0.0) {
      u += cp[0]*(StretchRPlateauG(cp[1], cp[2], cp[3], xi)
                  - xi*StretchRPlateauG(cp[1], cp[2], cp[3], 1.0));
    }
  }
  r = r0 + (r1-r0)*u;
}
//----------------------------------------------------------------------------------------
//! \fn ApplyRStretch()
//! \brief Apply whichever radial stretch is active to a position built by
//! CellCenterX/LeftEdgeX (which always come out UNSTRETCHED). Same definitions as
//! Coordinates uses, so a problem generator and the coordinate arrays cannot drift apart.
KOKKOS_INLINE_FUNCTION
void ApplyRStretch(const bool str_r, const Real fstr_r, const bool str_rp, const Real *c,
                   const Real r0, const Real r1, Real &r) {
  if (str_r) {
    StretchR(fstr_r, r0, r1, r);
  } else if (str_rp) {
    StretchRPoly(c, r0, r1, r);
  }
}

//----------------------------------------------------------------------------------------
//! \fn RadialCentroid()
//! \brief The volume centroid of the shell [r_l, r_r] -- what Coordinates::x1v stores on
//! BOTH spherical grids, and hence where a cell value lives.  NOT the midpoint.
KOKKOS_INLINE_FUNCTION
Real RadialCentroid(const Real r_l, const Real r_r) {
  const Real q = r_l/r_r;
  return 0.25*(q*q + 1.0)/((1.0/3.0)*(q*q + q + 1.0))*(r_r + r_l);
}

KOKKOS_INLINE_FUNCTION
void StretchTheta(const Real a, Real &t) {
  Real xi = t/M_PI;
  t = M_PI/2.0*(1.0+sinh(a*(2.0*xi-1.0))/sinh(a));
}
  


#endif // COORDINATES_GRID_STRETCH_HPP_
