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

//----------------------------------------------------------------------------------------
//! \fn Resistivity::CsVertexJFix
//! \brief cornerfix-1103 <mesh>/cs_vertex_3cell.  Pass 1's dual loops at a cube vertex
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
      // cornerfix3-1103: the pair's curvature term (half the second difference along the
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

  // cornerfix2-1103 <mesh>/cs_vertex_3cell: jn2 one layer further in k and jn3 in j (the
  // vertex point-reflection pairs of the ambipolar / mean_j vertex J, CsEdgeOps::X1);
  // extra entries only, every value computed before is unchanged
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

  // STRETCHED radial grid: the radial difference is taken AT the x1 face (see
  // RadialFaceWeights).  Off, or on an unstretched grid, this is the two-point loop,
  // bit for bit.
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
    // r1mspeed-1102: pass 1 as ONE launch (independent loops; bitwise)
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

  // cornerfix-1103 <mesh>/cs_vertex_3cell: no J entry may read the cube-vertex corner
  if (pmy_pack->pmesh->cs_vertex_3cell) {
    CsVertexJFix(b0, is, ie+1, false);
    CsVertexJFix(b0, is, ie, true);
  }
  const bool v3_ = pmy_pack->pmesh->cs_vertex_3cell;

  // fcorder3-1103 DIAGNOSTIC (env CS_RESJ_STAT=<every>, off unless set; no data change):
  // max |J.n| of the pass-1 current on the panel-edge (seam) edges vs >= 2 cells inside,
  // radial edges and tangential (x2/x3) edges, printed every <every>-th call by rank 0
  {static const char *rjs_ = std::getenv("CS_RESJ_STAT");
   static int rjcall_ = 0;
   if (rjs_ != nullptr && (rjcall_++ % std::max(1, std::atoi(rjs_))) == 0) {
    const int nx2_ = je - js + 1, nx3_ = ke - ks + 1, ni_ = ie - is + 1;
    const int ntot_ = (nmb1 + 1)*(nx3_ + 1)*(nx2_ + 1)*ni_;
    auto &msd = pmy_pack->pmb->mb_size;
    Real mx[4];
    for (int kd=0; kd<4; ++kd) {
      Real r_ = 0.0;
      Kokkos::parallel_reduce("res_cs_jstat", Kokkos::RangePolicy<>(DevExeSpace(), 0, ntot_),
      KOKKOS_LAMBDA(const int t, Real &rr) {
        const int i = is + t%ni_;
        const int j = js + (t/ni_)%(nx2_ + 1);
        const int k = ks + (t/(ni_*(nx2_ + 1)))%(nx3_ + 1);
        const int m = t/(ni_*(nx2_ + 1)*(nx3_ + 1));
        const auto z = msd.d_view(m);
        const Real h2 = (z.x2max - z.x2min)/nx2_, h3 = (z.x3max - z.x3min)/nx3_;
        // distance (faces) of face j / k from the panel edge
        const Real xf = z.x2min + (j - js)*h2, ef = z.x3min + (k - ks)*h3;
        const int dj = static_cast<int>(floor((1.0 - fabs(xf))/h2 + 0.5));
        const int dk = static_cast<int>(floor((1.0 - fabs(ef))/h3 + 0.5));
        const bool rad = (kd < 2);
        const bool seam = rad ? (dj == 0 || dk == 0) : (dj == 0 || dk == 0);
        const bool intr = (dj >= 2 && dk >= 2);
        if ((kd % 2 == 0) ? !seam : !intr) return;
        Real v = 0.0;
        if (rad) {
          v = fabs(jn1(m,k,j,i));
        } else {
          if (dj == 0 || intr) v = fmax(v, (k <= ke) ? fabs(jn3(m,k,j,i)) : 0.0);
          if (dk == 0 || intr) v = fmax(v, (j <= je) ? fabs(jn2(m,k,j,i)) : 0.0);
        }
        rr = fmax(rr, v);
      }, Kokkos::Max<Real>(r_));
      mx[kd] = r_;
    }
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, mx, 4, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif
    if (global_variable::my_rank == 0) {
      std::printf("### CS_RESJ_STAT call %d |J| max: radial seam %.4e interior %.4e |"
                  " tangential seam %.4e interior %.4e\n", rjcall_ - 1, mx[0], mx[1],
                  mx[2], mx[3]);
    }
  }}

  // ---- pass 2: rotate into the edge (covariant) frame and multiply by eta -------------
  // J.e_xi = (J.nhat_xi + c J.nhat_eta)/s and likewise for eta, the same transform
  // GnomonicEquiangleRaiseVelMHD applies to the field. The radial edge needs no rotation:
  // rhat is orthogonal to both tangent directions.
  auto e1 = efld.x1e;
  auto e2 = efld.x2e;
  auto e3 = efld.x3e;
  auto eta_b_ = eta_b;
  // ucthshare-1102 <mhd>/cs_seam_mhd = shared_g | shared_o: the seam edges take the SHARED
  // resistive EMF (MHD::CsSeamResist, after pass 2).  cs_seam_res = own: pass 2 skips
  // them; mean (default): pass 2 computes them and CsSeamResist takes the two-chart mean
  const bool rsc_ = (pmy_pack->pmhd != nullptr && pmy_pack->pmhd->cs_seam_mhd >= 3);
  const bool rsh_ = rsc_ && pmy_pack->pmhd->cs_seam_res_own;
  auto &msz_ = pmy_pack->pmb->mb_size;

  auto e1_ = KOKKOS_LAMBDA(int m, int k, int j, int i) {
    if (rsh_ && ((k == ks && msz_.d_view(m).x3min <= -1.0 + 1.0e-12) ||
                 (k == ke+1 && msz_.d_view(m).x3max >= 1.0 - 1.0e-12))) return;
    if (rsh_ && ((j == js && msz_.d_view(m).x2min <= -1.0 + 1.0e-12) ||
                 (j == je+1 && msz_.d_view(m).x2max >= 1.0 - 1.0e-12))) return;
    Real et = 0.25*(eta_b_(m,k,j,i) + eta_b_(m,k,j-1,i)
                  + eta_b_(m,k-1,j,i) + eta_b_(m,k-1,j-1,i));
    CsVtx3 vx;   // cornerfix-1103: a cube vertex has three cells, not four
    if (v3_ && CsVertexAt(msz_.d_view(m), k, j, js, je, ks, ke, vx)) {
      et = CsVertex3(eta_b_, m, vx, i);
    }
    e1(m,k,j,i) += et*jn1(m,k,j,i);
  };

  auto e2_ = KOKKOS_LAMBDA(int m, int k, int j, int i) {
    if (rsh_ && ((k == ks && msz_.d_view(m).x3min <= -1.0 + 1.0e-12) ||
                 (k == ke+1 && msz_.d_view(m).x3max >= 1.0 - 1.0e-12))) return;
    const Real et = 0.25*(eta_b_(m,k,j,i) + eta_b_(m,k,j,i-1)
                        + eta_b_(m,k-1,j,i) + eta_b_(m,k-1,j,i-1));
    const Real j3a = 0.25*(jn3(m,k,j,i) + jn3(m,k,j+1,i)
                         + jn3(m,k-1,j,i) + jn3(m,k-1,j+1,i));
    e2(m,k,j,i) += et*(jn2(m,k,j,i) + cet(m,k,j)*j3a)/set(m,k,j);
  };

  auto e3_ = KOKKOS_LAMBDA(int m, int k, int j, int i) {
    if (rsh_ && ((j == js && msz_.d_view(m).x2min <= -1.0 + 1.0e-12) ||
                 (j == je+1 && msz_.d_view(m).x2max >= 1.0 - 1.0e-12))) return;
    const Real et = 0.25*(eta_b_(m,k,j,i) + eta_b_(m,k,j-1,i)
                        + eta_b_(m,k,j,i-1) + eta_b_(m,k,j-1,i-1));
    const Real j2a = 0.25*(jn2(m,k,j,i) + jn2(m,k+1,j,i)
                         + jn2(m,k,j-1,i) + jn2(m,k+1,j-1,i));
    e3(m,k,j,i) += et*(jn3(m,k,j,i) + cxi(m,k,j)*j2a)/sxi(m,k,j);
  };
  // r1mspeed-1102: pass 2 as ONE launch (each edge loop reads only J; bitwise)
  par_for_fused3("res_cs_e", DevExeSpace(), nmb1, LoopBox3{ks, ke+1, js, je+1, is, ie}, e1_,
                 LoopBox3{ks, ke+1, js, je, is, ie+1}, e2_,
                 LoopBox3{ks, ke, js, je+1, is, ie+1}, e3_);
  // with ambipolar diffusion the two-chart mean is taken once, after the ambipolar EMF
  // has been added too (AddEMFAmbipolarGnomonic); the own-grid Ohmic part is added here
  if (rsc_ && (rsh_ || ad_type == 0)) {
    pmy_pack->pmhd->CsSeamResist(b0, eta_b, jnorm, efld);
  }

  return;
}
