//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file resistivity_gnomonic.cpp
//! \brief The resistive EMF on the cubed sphere. The face fields b0.x2f / b0.x3f hold
//! B.nhat, and the tangent basis is not orthogonal: Stokes' theorem needs the component
//! along each edge, B.e_xi = (B.nhat_xi + c B.nhat_eta)/s (c = cos_cell, s = sin_cell).
//! The curl is taken in two passes. Pass one applies Stokes' theorem on the loop through
//! the four cell centres around each edge (lengths from Coordinates::dxface, area from
//! Coordinates::areaedge), which gives J in the face-normal frame {rhat, nhat_xi,
//! nhat_eta}. Pass two rotates it into the frame of the edge EMFs (E.that along the
//! edge, as GnomonicEquiangleEmfX1 produces for the ideal EMF) and multiplies by eta.
//! The rotation mixes J at x2 edges with J at x3 edges, so the passes are separate.

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "resistivity.hpp"
#include "mhd/mhd.hpp"
#include "diffusion/cs_vertex_res.hpp"
#include "diffusion/radial_face_weights.hpp"

//----------------------------------------------------------------------------------------
//! \fn Real BcovXi
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
//! \fn Real BcovEta
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
//! \fn void Resistivity::CsVertexJFix
//! \brief <mesh>/cs_vertex_3cell: pass 1's dual loops at a cube vertex
//! run through the made-up corner block (ng x ng diagonal ghosts: no physical cells).
//! j1 = true: the vertex RADIAL current from the 3-cell loop through this block's corner
//! cell and the two seam ghost cells.  Each panel forms the part of the loop it owns,
//! its two seam-face segments (Stokes, B.e from the same BcovXi / BcovEta as pass 1)
//! over its quarter of the edge area, and the three panels' values are averaged
//! afterwards with the EMF (the per-panel half-loop H/(A/4) averages to the full 3-cell
//! loop sum H/(3A/4): the three corner cells are congruent).  j1 = false: the x3 edge
//! (kn, jc) of the ghost row and the x2 edge (kc, jn) of the ghost column at the vertex
//! (corner cell edges, read by the transverse J averages of the Ohmic metric term and of
//! ambipolar diffusion) are extrapolated linearly along their own ghost row / column
//! from the two real edges next to them (constant if the block has < 3 cells).

void Resistivity::CsVertexJFix(const DvceFaceFld4D<Real> &b0, const int il,
                               const int iu, const bool j1) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int js = indcs.js, je = indcs.je, ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto &msz = pmy_pack->pmb->mb_size;
  auto &cxi = pmy_pack->pcoord->cos_face_xi;
  auto &sxi = pmy_pack->pcoord->sin_face_xi;
  auto &cet = pmy_pack->pcoord->cos_face_eta;
  auto &set = pmy_pack->pcoord->sin_face_eta;
  auto &dxf2 = pmy_pack->pcoord->dxface.x2f;
  auto &dxf3 = pmy_pack->pcoord->dxface.x3f;
  auto &ae1 = pmy_pack->pcoord->areaedge.x1e;
  auto b = b0;
  auto jn1 = jnorm.x1e, jn2 = jnorm.x2e, jn3 = jnorm.x3e;
  par_for("res_cs_vtxj", DevExeSpace(), 0, nmb1, 0, 3, il, iu,
  KOKKOS_LAMBDA(const int m, const int cn, const int i) {
    const bool jlo = (cn & 1) == 0, klo = (cn & 2) == 0;
    const int jc = jlo ? js : je+1, kc = klo ? ks : ke+1;
    CsVtx3 v;
    if (!CsVertexAt(msz.d_view(m), kc, jc, js, je, ks, ke, v)) return;
    const int dj = jlo ? 1 : -1, dk = klo ? 1 : -1;
    if (j1) {
      // pass 1's dual-loop J.r (res_cs_j) at a radial edge (k, j) of the seam GHOST
      // region: corner-free for the two edges used here
      auto j1at = [&](const int k, const int j) {
        return (dxf3(m,k,j,i)*BcovEta(b,cet,set,m,k,j,i)
              - dxf3(m,k,j-1,i)*BcovEta(b,cet,set,m,k,j-1,i)
              - dxf2(m,k,j,i)*BcovXi(b,cxi,sxi,m,k,j,i)
              + dxf2(m,k-1,j,i)*BcovXi(b,cxi,sxi,m,k-1,j,i))/ae1(m,k,j,i);
      };
      // the vertex is the midpoint of the point-reflected pair (one step into the x3-seam
      // ghost row, one into the x2-seam ghost column, diagonally)
      const Real p1 = j1at(kc-dk, jc+dj), p2 = j1at(kc+dk, jc-dj);
      // the pair's curvature term (half the second difference along the
      // pair's direction, centred on the own interior edge (kc+dk, jc+dj), from the two
      // own seam edges next to the vertex), smoothness-gated and bounded (cs_fc3)
      const Real fa = jn1(m,kc,jc+2*dj,i), fb = jn1(m,kc+dk,jc+dj,i);
      const Real fc = jn1(m,kc+2*dk,jc,i);
      const Real lo = fmin(fmin(fmin(p1, p2), fmin(fa, fb)), fc);
      const Real hi = fmax(fmax(fmax(p1, p2), fmax(fa, fb)), fc);
      jn1(m,kc,jc,i) = cs_fc3::CurvCorrect(0.5*(p1 + p2), fa, fb, fc, 0.5, lo, hi);
      return;
    }
    // the corner cell's own transverse edges (x3 edge at (kn, jc), x2 edge at (kc, jn))
    // loop through the corner block.  They only enter 4-point averages (the Ohmic metric
    // term, ambipolar / mean_j transverse J at the vertex-cell along-seam edges); the
    // parallelogram value makes each of those averages the midpoint of the two real
    // entries on the diagonal through the average's point (linear-exact, convex)
    jn3(m,v.kn,jc,i) = jn3(m,v.kv,jc,i) + jn3(m,v.kn,jc+dj,i) - jn3(m,v.kv,jc+dj,i);
    jn2(m,kc,v.jn,i) = jn2(m,kc,v.jv,i) + jn2(m,kc+dk,v.jn,i) - jn2(m,kc+dk,v.jv,i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void Resistivity::AddEMFGnomonicResist
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

  // <mesh>/cs_vertex_3cell: jn2 one layer further in k and jn3 in j (the vertex
  // point-reflection pairs of the ambipolar / mean_j vertex J, CsEdgeOps::X1)
  const int v3w = pmy_pack->pmesh->cs_vertex_3cell ? 1 : 0;
  // ---- pass 1: J in the face-normal frame, by Stokes' theorem on the dual loop --------
  // Each loop is bounded by the four cell centres around the edge, so every side passes
  // through a face centre and carries that face's field. Pass 2 needs J one layer beyond
  // its own edge range in the two transverse directions, hence the -1 / +1 below.
  auto j1_ = KOKKOS_LAMBDA(int m, int k, int j, int i) {
    jn1(m,k,j,i) = (dxf3(m,k,j,i)*BcovEta(b,cet,set,m,k,j,i)
                  - dxf3(m,k,j-1,i)*BcovEta(b,cet,set,m,k,j-1,i)
                  - dxf2(m,k,j,i)*BcovXi(b,cxi,sxi,m,k,j,i)
                  + dxf2(m,k-1,j,i)*BcovXi(b,cxi,sxi,m,k-1,j,i))/ae1(m,k,j,i);
  };

  // stretched radial grid: the radial difference is taken at the x1 face (see
  // RadialFaceWeights); otherwise the two-point loop
  const bool x1c = cs_resist_x1_centred && (pmy_pack->pmesh->use_grid_stretch_r ||
                                            pmy_pack->pmesh->use_grid_stretch_r_poly);
  if (!x1c) {
    auto j2_ = KOKKOS_LAMBDA(int m, int k, int j, int i) {
      jn2(m,k,j,i) = (-dxf3(m,k,j,i)*BcovEta(b,cet,set,m,k,j,i)
                      + dxf3(m,k,j,i-1)*BcovEta(b,cet,set,m,k,j,i-1)
                      + dxf1(m,k,j,i)*b.x1f(m,k,j,i)
                      - dxf1(m,k-1,j,i)*b.x1f(m,k-1,j,i))/ae2(m,k,j,i);
    };

    auto j3_ = KOKKOS_LAMBDA(int m, int k, int j, int i) {
      jn3(m,k,j,i) = (dxf2(m,k,j,i)*BcovXi(b,cxi,sxi,m,k,j,i)
                    - dxf2(m,k,j,i-1)*BcovXi(b,cxi,sxi,m,k,j,i-1)
                    - dxf1(m,k,j,i)*b.x1f(m,k,j,i)
                    + dxf1(m,k,j-1,i)*b.x1f(m,k,j-1,i))/ae3(m,k,j,i);
    };
    // pass 1 as one launch (independent loops)
    par_for_fused3("res_cs_j", DevExeSpace(), nmb1,
                   LoopBox3{ks, ke+1, js, je+1, is, ie}, j1_,
                   LoopBox3{ks-v3w, ke+1+v3w, js-1, je+1, is, ie+1}, j2_,
                   LoopBox3{ks-1, ke+1, js-v3w, je+1+v3w, is, ie+1}, j3_);
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
    auto j2c_ = KOKKOS_LAMBDA(int m, int k, int j, int i) {
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
      // ae2 = rbar*(r_c(i)-r_c(i-1))*arc is the loop area at the loop's MEAN radius;
      // the curl AT the face divides by r_l instead (curl = (1/r)[...] at r = r_l)
      const Real rsc = 0.5*(rc(m,i) + rc(m,i-1))/rf(m,i);
      jn2(m,k,j,i) = rsc*(-dg
                      + dxf1(m,k,j,i)*b.x1f(m,k,j,i)
                      - dxf1(m,k-1,j,i)*b.x1f(m,k-1,j,i))/ae2(m,k,j,i);
    };

    auto j3c_ = KOKKOS_LAMBDA(int m, int k, int j, int i) {
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
      const Real rsc = 0.5*(rc(m,i) + rc(m,i-1))/rf(m,i);   // as in res_cs_j2c
      jn3(m,k,j,i) = rsc*(dg
                    - dxf1(m,k,j,i)*b.x1f(m,k,j,i)
                    + dxf1(m,k,j-1,i)*b.x1f(m,k,j-1,i))/ae3(m,k,j,i);
    };
    par_for_fused3("res_cs_jc", DevExeSpace(), nmb1,
                   LoopBox3{ks, ke+1, js, je+1, is, ie}, j1_,
                   LoopBox3{ks-v3w, ke+1+v3w, js-1, je+1, is, ie+1}, j2c_,
                   LoopBox3{ks-1, ke+1, js-v3w, je+1+v3w, is, ie+1}, j3c_);
  }

  // <mesh>/cs_vertex_3cell: no J entry may read the cube-vertex corner
  if (pmy_pack->pmesh->cs_vertex_3cell) {
    CsVertexJFix(b0, is, ie+1, false);
    CsVertexJFix(b0, is, ie, true);
  }
  const bool v3_ = pmy_pack->pmesh->cs_vertex_3cell;

  // ---- pass 2: rotate into the edge (covariant) frame and multiply by eta -------------
  // J.e_xi = (J.nhat_xi + c J.nhat_eta)/s and likewise for eta, the same transform
  // GnomonicEquiangleRaiseVelMHD applies to the field. The radial edge needs no rotation:
  // rhat is orthogonal to both tangent directions.
  auto e1 = efld.x1e;
  auto e2 = efld.x2e;
  auto e3 = efld.x3e;
  auto eta_b_ = eta_b;
  // cs_seam_mhd: the seam edges take the two-chart mean of the resistive EMF
  // (MHD::CsSeamResist, after the edge EMFs below)
  const bool rsc_ = (pmy_pack->pmhd != nullptr && pmy_pack->pmhd->cs_seam_mhd);
  auto &msz_ = pmy_pack->pmb->mb_size;

  auto e1_ = KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real et = 0.25*(eta_b_(m,k,j,i) + eta_b_(m,k,j-1,i)
                  + eta_b_(m,k-1,j,i) + eta_b_(m,k-1,j-1,i));
    CsVtx3 vx;   // a cube vertex has three cells, not four
    if (v3_ && CsVertexAt(msz_.d_view(m), k, j, js, je, ks, ke, vx)) {
      et = CsVertex3(eta_b_, m, vx, i);
    }
    e1(m,k,j,i) += et*jn1(m,k,j,i);
  };

  auto e2_ = KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real et = 0.25*(eta_b_(m,k,j,i) + eta_b_(m,k,j,i-1)
                        + eta_b_(m,k-1,j,i) + eta_b_(m,k-1,j,i-1));
    const Real j3a = 0.25*(jn3(m,k,j,i) + jn3(m,k,j+1,i)
                         + jn3(m,k-1,j,i) + jn3(m,k-1,j+1,i));
    e2(m,k,j,i) += et*(jn2(m,k,j,i) + cet(m,k,j)*j3a)/set(m,k,j);
  };

  auto e3_ = KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real et = 0.25*(eta_b_(m,k,j,i) + eta_b_(m,k,j-1,i)
                        + eta_b_(m,k,j,i-1) + eta_b_(m,k,j-1,i-1));
    const Real j2a = 0.25*(jn2(m,k,j,i) + jn2(m,k+1,j,i)
                         + jn2(m,k,j-1,i) + jn2(m,k+1,j-1,i));
    e3(m,k,j,i) += et*(jn3(m,k,j,i) + cxi(m,k,j)*j2a)/sxi(m,k,j);
  };
  // pass 2 as one launch (each edge loop reads only J)
  par_for_fused3("res_cs_e", DevExeSpace(), nmb1,
                 LoopBox3{ks, ke+1, js, je+1, is, ie}, e1_,
                 LoopBox3{ks, ke+1, js, je, is, ie+1}, e2_,
                 LoopBox3{ks, ke, js, je+1, is, ie+1}, e3_);
  // with ambipolar diffusion the two-chart mean is taken once, after the ambipolar EMF
  // has been added too (AddEMFAmbipolarGnomonic); the own-grid Ohmic part is added here
  if (rsc_ && ad_type == 0) {
    pmy_pack->pmhd->CsSeamResist(b0, eta_b, jnorm, efld);
  }

  return;
}
