//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file resistivity_gnomonic.cpp
//! \brief The resistive EMF on the CUBED SPHERE.
//!
//! WHY THIS NEEDS ITS OWN FILE. `current_density.hpp` branches on `use_spherical_polar`
//! only, so on the cubed sphere it fell through to the Cartesian formula and divided by
//! `size.dx1/dx2/dx3` -- which on this grid are index spacings, not lengths (x2 and x3
//! are equiangular coordinates on [-1,1]). That is wrong by an O(1) factor everywhere.
//! Its curvilinear branch is not usable either, for a reason that is not a matter of
//! filling in different metric coefficients:
//!
//!   THE TANGENT BASIS IS NOT ORTHOGONAL. `b0.x2f` and `b0.x3f` hold B.nhat, the flux
//!   density through the face, and on the cubed sphere nhat_xi is perpendicular to
//!   e_eta, NOT parallel to e_xi. Stokes' theorem needs the component ALONG each edge,
//!   B.e_xi, and the two differ by the metric:
//!       B.e_xi = (B.nhat_xi + c B.nhat_eta)/s,   c = cos_cell, s = sin_cell,
//!   so the correction is O(c) -- O(1) near a panel corner, not O(h). This is the same
//!   trap that made the cell-centred EMF inconsistent in mhd_corner_e.cpp.
//!
//! The curl is therefore taken in two passes. Pass one applies Stokes' theorem on the
//! loop through the four cell CENTRES around each edge (lengths from Coordinates::dxface,
//! enclosed area from Coordinates::areaedge), which yields J in the FACE-NORMAL frame
//! {rhat, nhat_xi, nhat_eta}. Pass two rotates it into the frame the edge EMFs are stored
//! in -- E.that along the edge, which is what mhd_ct.cpp's dxedge*E/area consumes and
//! what GnomonicEquiangleEmfX1 produces for the ideal EMF -- and multiplies by eta. The
//! two passes cannot be fused: the rotation mixes J at x2-edges with J at x3-edges, which
//! are different points.

#include <iostream>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "resistivity.hpp"
#include "mhd/mhd.hpp"

//----------------------------------------------------------------------------------------
//! \brief B.e_xi at the centre of an x2 face, from the face-normal components. The
//! partner component lives on the four x3 faces that share this face's corners.

KOKKOS_INLINE_FUNCTION
Real BcovXi(const DvceFaceFld4D<Real> &b, const DvceArray3D<Real> &cosf,
            const DvceArray3D<Real> &sinf, const int m, const int k, const int j,
            const int i) {
  const Real b3a = 0.25*(b.x3f(m,k,j,i) + b.x3f(m,k+1,j,i)
                       + b.x3f(m,k,j-1,i) + b.x3f(m,k+1,j-1,i));
  return (b.x2f(m,k,j,i) + cosf(m,k,j)*b3a)/sinf(m,k,j);
}

//----------------------------------------------------------------------------------------
//! \brief B.e_eta at the centre of an x3 face.

KOKKOS_INLINE_FUNCTION
Real BcovEta(const DvceFaceFld4D<Real> &b, const DvceArray3D<Real> &cosf,
             const DvceArray3D<Real> &sinf, const int m, const int k, const int j,
             const int i) {
  const Real b2a = 0.25*(b.x2f(m,k,j,i) + b.x2f(m,k,j+1,i)
                       + b.x2f(m,k-1,j,i) + b.x2f(m,k-1,j+1,i));
  return (b.x3f(m,k,j,i) + cosf(m,k,j)*b2a)/sinf(m,k,j);
}

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
//! 1.2e-2 field error below 1e-6 bar against 6e-4 with these weights.
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

//----------------------------------------------------------------------------------------
//! \fn Resistivity::AddEMFGnomonicResist
//! \brief Adds E_resistive = eta J to the edge-centred electric fields, cubed sphere.

void Resistivity::AddEMFGnomonicResist(const DvceFaceFld4D<Real> &b0,
                                       DvceEdgeFld4D<Real> &efld) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  if (!(pmy_pack->pmesh->three_d)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Resistivity on the cubed sphere requires a 3D mesh" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto &cxi = pmy_pack->pcoord->cos_face_xi;
  auto &sxi = pmy_pack->pcoord->sin_face_xi;
  auto &cet = pmy_pack->pcoord->cos_face_eta;
  auto &set = pmy_pack->pcoord->sin_face_eta;
  auto &dxf1 = pmy_pack->pcoord->dxface.x1f;
  auto &dxf2 = pmy_pack->pcoord->dxface.x2f;
  auto &dxf3 = pmy_pack->pcoord->dxface.x3f;
  auto &ae1 = pmy_pack->pcoord->areaedge.x1e;
  auto &ae2 = pmy_pack->pcoord->areaedge.x2e;
  auto &ae3 = pmy_pack->pcoord->areaedge.x3e;
  auto b = b0;
  auto jn1 = jnorm.x1e;
  auto jn2 = jnorm.x2e;
  auto jn3 = jnorm.x3e;

  // ---- pass 1: J in the face-normal frame, by Stokes' theorem on the dual loop --------
  // Each loop is bounded by the four cell centres around the edge, so every side passes
  // through a face centre and carries that face's field. Pass 2 needs J one layer beyond
  // its own edge range in the two transverse directions, hence the -1 / +1 below.
  par_for("res_cs_j1", DevExeSpace(), 0,nmb1, ks,ke+1, js,je+1, is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    jn1(m,k,j,i) = (dxf3(m,k,j,i)*BcovEta(b,cet,set,m,k,j,i)
                  - dxf3(m,k,j-1,i)*BcovEta(b,cet,set,m,k,j-1,i)
                  - dxf2(m,k,j,i)*BcovXi(b,cxi,sxi,m,k,j,i)
                  + dxf2(m,k-1,j,i)*BcovXi(b,cxi,sxi,m,k-1,j,i))/ae1(m,k,j,i);
  });

  // STRETCHED radial grid: the radial difference is taken AT the x1 face (see
  // RadialFaceWeights).  Off, or on an unstretched grid, this is the two-point loop,
  // bit for bit.
  const bool x1c = cs_resist_x1_centred && (pmy_pack->pmesh->use_grid_stretch_r ||
                                            pmy_pack->pmesh->use_grid_stretch_r_poly);
  if (!x1c) {
    par_for("res_cs_j2", DevExeSpace(), 0,nmb1, ks,ke+1, js-1,je+1, is,ie+1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      jn2(m,k,j,i) = (-dxf3(m,k,j,i)*BcovEta(b,cet,set,m,k,j,i)
                      + dxf3(m,k,j,i-1)*BcovEta(b,cet,set,m,k,j,i-1)
                      + dxf1(m,k,j,i)*b.x1f(m,k,j,i)
                      - dxf1(m,k-1,j,i)*b.x1f(m,k-1,j,i))/ae2(m,k,j,i);
    });

    par_for("res_cs_j3", DevExeSpace(), 0,nmb1, ks-1,ke+1, js,je+1, is,ie+1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      jn3(m,k,j,i) = (dxf2(m,k,j,i)*BcovXi(b,cxi,sxi,m,k,j,i)
                    - dxf2(m,k,j,i-1)*BcovXi(b,cxi,sxi,m,k,j,i-1)
                    - dxf1(m,k,j,i)*b.x1f(m,k,j,i)
                    + dxf1(m,k,j-1,i)*b.x1f(m,k,j-1,i))/ae3(m,k,j,i);
    });
  } else {
    auto &rc = pmy_pack->pcoord->x1v;    // cell centroids, ghosts included
    auto &rf = pmy_pack->pcoord->xx1f;   // x1 faces
    auto &mbbcs = pmy_pack->pmb->mb_bcs;
    if (indcs.ng < 2) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<mhd>/cs_resist_x1_centred needs nghost >= 2"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    par_for("res_cs_j2c", DevExeSpace(), 0,nmb1, ks,ke+1, js-1,je+1, is,ie+1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      // +1: physical inner radial boundary on this face; -1: outer; 0: interior
      int sd = 0;
      if (i == is && mbbcs.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::block) {
        sd = 1;
      }
      if (i == ie+1 && mbbcs.d_view(m,BoundaryFace::outer_x1) != BoundaryFlag::block) {
        sd = -1;
      }
      Real w[4];
      RadialFaceWeights(rc(m,i-2), rc(m,i-1), rc(m,i), rc(m,i+1), rf(m,i), sd, w);
      Real dg = 0.0;
      for (int s=0; s<4; ++s) {
        if (w[s] != 0.0) dg += w[s]*dxf3(m,k,j,i-2+s)*BcovEta(b,cet,set,m,k,j,i-2+s);
      }
      jn2(m,k,j,i) = (-dg
                      + dxf1(m,k,j,i)*b.x1f(m,k,j,i)
                      - dxf1(m,k-1,j,i)*b.x1f(m,k-1,j,i))/ae2(m,k,j,i);
    });

    par_for("res_cs_j3c", DevExeSpace(), 0,nmb1, ks-1,ke+1, js,je+1, is,ie+1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      // +1: physical inner radial boundary on this face; -1: outer; 0: interior
      int sd = 0;
      if (i == is && mbbcs.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::block) {
        sd = 1;
      }
      if (i == ie+1 && mbbcs.d_view(m,BoundaryFace::outer_x1) != BoundaryFlag::block) {
        sd = -1;
      }
      Real w[4];
      RadialFaceWeights(rc(m,i-2), rc(m,i-1), rc(m,i), rc(m,i+1), rf(m,i), sd, w);
      Real dg = 0.0;
      for (int s=0; s<4; ++s) {
        if (w[s] != 0.0) dg += w[s]*dxf2(m,k,j,i-2+s)*BcovXi(b,cxi,sxi,m,k,j,i-2+s);
      }
      jn3(m,k,j,i) = (dg
                    - dxf1(m,k,j,i)*b.x1f(m,k,j,i)
                    + dxf1(m,k,j-1,i)*b.x1f(m,k,j-1,i))/ae3(m,k,j,i);
    });
  }

  // ---- pass 2: rotate into the edge (covariant) frame and multiply by eta -------------
  // J.e_xi = (J.nhat_xi + c J.nhat_eta)/s and likewise for eta, the same transform
  // GnomonicEquiangleRaiseVelMHD applies to the field. The radial edge needs no rotation:
  // rhat is orthogonal to both tangent directions.
  auto e1 = efld.x1e;
  auto e2 = efld.x2e;
  auto e3 = efld.x3e;
  auto eta_b_ = eta_b;

  par_for("res_cs_e1", DevExeSpace(), 0,nmb1, ks,ke+1, js,je+1, is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real et = 0.25*(eta_b_(m,k,j,i) + eta_b_(m,k,j-1,i)
                        + eta_b_(m,k-1,j,i) + eta_b_(m,k-1,j-1,i));
    e1(m,k,j,i) += et*jn1(m,k,j,i);
  });

  par_for("res_cs_e2", DevExeSpace(), 0,nmb1, ks,ke+1, js,je, is,ie+1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real et = 0.25*(eta_b_(m,k,j,i) + eta_b_(m,k,j,i-1)
                        + eta_b_(m,k-1,j,i) + eta_b_(m,k-1,j,i-1));
    const Real j3a = 0.25*(jn3(m,k,j,i) + jn3(m,k,j+1,i)
                         + jn3(m,k-1,j,i) + jn3(m,k-1,j+1,i));
    e2(m,k,j,i) += et*(jn2(m,k,j,i) + cet(m,k,j)*j3a)/set(m,k,j);
  });

  par_for("res_cs_e3", DevExeSpace(), 0,nmb1, ks,ke, js,je+1, is,ie+1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real et = 0.25*(eta_b_(m,k,j,i) + eta_b_(m,k,j-1,i)
                        + eta_b_(m,k,j,i-1) + eta_b_(m,k,j-1,i-1));
    const Real j2a = 0.25*(jn2(m,k,j,i) + jn2(m,k+1,j,i)
                         + jn2(m,k,j-1,i) + jn2(m,k+1,j-1,i));
    e3(m,k,j,i) += et*(jn3(m,k,j,i) + cxi(m,k,j)*j2a)/sxi(m,k,j);
  });

  return;
}
