#ifndef DIFFUSION_RADIAL_FACE_WEIGHTS_HPP_
#define DIFFUSION_RADIAL_FACE_WEIGHTS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file radial_face_weights.hpp
//! \brief share-1104 G9: RadialFaceWeights, moved verbatim from resistivity_gnomonic.cpp
//! so that the spherical-polar current (current_density.hpp, <mhd>/sp_resist_x1_centred)
//! shares the cubed sphere's radial-face derivative (<mhd>/cs_resist_x1_centred).

#include "athena.hpp"

//----------------------------------------------------------------------------------------
//! \brief Weights w[0..3] on the cell centres i-2, i-1, i, i+1 such that
//!     sum_s w[s] g(i-2+s) = (r_c(i) - r_c(i-1)) * dg/dr AT THE x1 FACE r_l(i).
//!
//! WHY. The x2/x3 Stokes loops of pass 1 run through the cell CENTRES r_c(i-1), r_c(i),
//! so their radial difference g(i) - g(i-1) (g = dxface*B, i.e. arc*r*B_t) is the
//! derivative at the MIDPOINT of the two centres, while the EMF it produces is used at
//! the edge on the face r_l(i).  On a uniform grid the two coincide to O(dr^2/r); on a
//! STRETCHED grid they differ by ~(dr_i - dr_(i-1))/4, a truncation error that is only
//! as small as the cell-to-cell width ratio is close to 1.  On the WASP-121b grid (8-term
//! polynomial, width ratio up to 1.09 below 1e-6 bar and 1.31 at the top) it dominated
//! the resistive diffusion error: a 1-D replica of this operator (csresist_0928) gave a
//! 1.2e-2 field error below 1e-6 bar against 2.3e-4 with these weights (and the
//! loop normalized at r_l, the rsc factor in pass 1).
//!
//! The derivative at r_l is the MEAN of the two quadratics through (i-1, i, i+1) and
//! (i-2, i-1, i) -- symmetric, second order at r_l -- or one of them alone on a physical
//! radial boundary face (side = +1: inner boundary, (i-1, i, i+1); side = -1: outer,
//! (i-2, i-1, i)), so no stencil reaches further into the ghosts than today's.  The
//! spectrum of the resulting 1-D diffusion operator stays real and negative with eta
//! varying by 1e4 from cell to cell, and its spectral radius (hence the stable dt) moves
//! by < 0.1 %.

KOKKOS_INLINE_FUNCTION
void RadialFaceWeights(const Real r0, const Real r1, const Real r2, const Real r3,
                       const Real rl, const int side, Real w[4]) {
  // derivative of the quadratic through (xa, xb, xc) at rl: weights on the three nodes
  auto dq = [](Real xa, Real xb, Real xc, Real x, Real &wa, Real &wb, Real &wc) {
    wa = ((x - xb) + (x - xc))/((xa - xb)*(xa - xc));
    wb = ((x - xa) + (x - xc))/((xb - xa)*(xb - xc));
    wc = ((x - xa) + (x - xb))/((xc - xa)*(xc - xb));
  };
  Real ra1 = 0.0, ra2 = 0.0, ra3 = 0.0;   // right stencil (i-1, i, i+1)
  Real la0 = 0.0, la1 = 0.0, la2 = 0.0;   // left  stencil (i-2, i-1, i)
  dq(r1, r2, r3, rl, ra1, ra2, ra3);
  dq(r0, r1, r2, rl, la0, la1, la2);
  const Real fr = (side > 0) ? 1.0 : ((side < 0) ? 0.0 : 0.5);
  const Real fl = 1.0 - fr;
  const Real dr = r2 - r1;
  w[0] = dr*(fl*la0);
  w[1] = dr*(fl*la1 + fr*ra1);
  w[2] = dr*(fl*la2 + fr*ra2);
  w[3] = dr*(fr*ra3);
}


#endif // DIFFUSION_RADIAL_FACE_WEIGHTS_HPP_
