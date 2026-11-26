//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_sph_facesum.cpp
//! \brief spherical polar: the r faces of the closure_consistency face sums made in the
//! sweep, as the Cartesian slab's (rays/sc_sph_sweep_body.inc, SLAB): a ray of direction
//! d swept into cell c carries through its r face (between c and its upwind shell) the
//! intensity where its segment crosses that face, an entry the face's inflow, a turned
//! ray (no r face of its own) the cell below, and a ray leaving through the inner or top
//! face its intensity there.  Every r face then holds each direction once.  The lateral
//! faces stay RayFaceSums' upwind-cell rule, which on the curved shells is the slab's
//! foot-shift flux of a shift below one cell (the upwind cell of the face's own shell,
//! the receiving ray's upwind shell) without its S2 / S3 halves: those use the
//! receiving cell's own shift, which differs from its neighbour's and would break the
//! mirror symmetry of the sums.  A uniform field gives every face I (n . A): the closed
//! vector areas make each cell's balance zero for every ray.  The sums are added per
//! cell in ray order (no atomics) as Cartesian (energy, momentum) sums into vgd_cfs;
//! RayFaceSums projects them on the face bases.
//!
//!  called from: rays/sc_sph_sweep_shell.inc, rad_vet_closure_faces.cpp

#include <math.h>

#include <cstdint>

#include "athena.hpp"
// FMA contraction only within one expression in this TU
#if defined(__clang__)
#pragma clang fp contract(on)
#endif
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_moments.hpp"
#include "rad_vet/rays/sc_segment.hpp"
#include "rad_vet/rays/sc_halo.hpp"
#include "rad_vet/rays/sc_sph_facesum.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::ScSpFsOn
//! \brief the spherical-polar sweep makes the closure_consistency face sums (VET_SP_FS)

bool RadiationVet::ScSpFsOn() const {
#if VET_SP_FS
  return cln_consist && vgd_on && vgd_slab != 1 &&
         (pmy_pack->pmesh->use_spherical_polar || vgd_csg);
#else
  return false;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScSpFaceZero
//! \brief before a sweep with the face sums: the sums zeroed, the rays' scalar buffer
//! sized (a fused sweep: both sets')

void RadiationVet::ScSpFaceZero(const DevExeSpace &ex) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack, n = vgd_n;
  if (vgd_nb1 > 1) {
    ScFatal("closure_consistency: the face sums the spherical-polar sweep makes need "
            "one MeshBlock along r (a build with -D VET_SP_FS=0 takes RayFaceSums)");
  }
  for (int x = 0; x < (vgd_xff_now ? 2 : 1); ++x) {
    auto &rg = x ? vgd_rgs_xb : vgd_rgs;
    auto &cf = x ? vgd_cfs_xb : vgd_cfs;
    if (rg.extent_int(0) != 2*nmb || rg.extent_int(1) != VGD_SPFS_NS*n) {
      Kokkos::realloc(rg, 2*nmb, VGD_SPFS_NS*n, indcs.nx3, indcs.nx2);
    }
    if (cf.extent_int(0) != nmb) {
      Kokkos::realloc(cf, nmb, 4, indcs.nx3, indcs.nx2, indcs.nx1 + 1);
    }
    Kokkos::deep_copy(ex, cf, 0.0);
  }
  // the r faces' closed vector areas over r_f^2 (FaceVector; both meshes scale as r^2)
  if (vgd_cfa.extent_int(0) != nmb) {
    Kokkos::realloc(vgd_cfa, nmb, 3, indcs.nx3, indcs.nx2);
    auto fa_ = vgd_cfa;
    auto x3f = pmy_pack->pcoord->xx3f;
    const int ks = indcs.ks, js = indcs.js, is = indcs.is;
    const SolveArgs *px = DevArgs(MakeArgs(), 3);
    par_for("vet_vgd_spfs_area", ex, 0, nmb - 1, 0, indcs.nx3 - 1, 0, indcs.nx2 - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      const SolveArgs &x = *px;
      Real a[3];
      FaceVector(x, x3f, 0, m, ks + k, js + j, is, a);
      const Real r2 = SQR(x.x1f(m,is));
      for (int c = 0; c < 3; ++c) {fa_(m,c,k,j) = a[c]/r2;}
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScSpFacePass
//! \brief at the start of a pass: no ray of either shell slot counts until swept (the
//! pass's rays are the same on every shell; a fused sweep: both sets)

void RadiationVet::ScSpFacePass(const DevExeSpace &ex) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack, n = vgd_n;
  for (int x = 0; x < (vgd_xff_now ? 2 : 1); ++x) {
    auto rg_ = x ? vgd_rgs_xb : vgd_rgs;
    par_for("vet_vgd_spfs_pass", ex, 0, 2*nmb - 1, 0, n - 1, 0, indcs.nx3 - 1, 0,
            indcs.nx2 - 1,
    KOKKOS_LAMBDA(const int b, const int d, const int k, const int j) {
      rg_(b,SPFS_XV*n + d,k,j) = -1.0;
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScSpFaceSum
//! \brief the r-face sums of the rays swept into shell i (pass inw), after its interior
//! and edge items: one thread per cell, the rays in order, the four sums (energy,
//! Cartesian momentum) through the rays' r face (inward: the cell's high face) and
//! through the exits (the inner or top face)

void RadiationVet::ScSpFaceSum(const int i, const int l, const bool inw,
                               const DevExeSpace &ex, const bool alt) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is;
  const int nmb = pmy_pack->nmb_thispack, n = vgd_n, ib = l & 1;
  const int n1 = indcs.nx1, lcut = vgd_scut;
  auto rg_ = alt ? vgd_rgs_xb : vgd_rgs;
  auto dir_ = alt ? vgd_tabd.dir : vgd_dir;
  auto cf_ = alt ? vgd_cfs_xb : vgd_cfs;
  auto fa_ = vgd_cfa;
  auto x1f = pmy_pack->pcoord->xx1f;
  auto x1f_ = vgd_x1f;
  const int ir = inw ? (i + 1) : i, ix = inw ? i : (i + 1);
  // (exits: the first shell's inward rays, the top shell's outward rays)
  const bool xs = inw ? (l == lcut) : (l == n1 - 1);
  par_for("vet_vgd_spfs_sum", ex, 0, nmb - 1, 0, indcs.nx3 - 1, 0, indcs.nx2 - 1,
  KOKKOS_LAMBDA(const int m, const int kc, const int jc) {
    const int b = ib*nmb + m;
    const bool xb = xs && (x1f_(m) & (inw ? 1 : 2));
    const Real rr = SQR(x1f(m,ir)), rx = SQR(x1f(m,ix));
    const Real a0 = fa_(m,0,kc,jc), a1 = fa_(m,1,kc,jc), a2 = fa_(m,2,kc,jc);
    Real sr[4] = {0.0, 0.0, 0.0, 0.0}, sx[4] = {0.0, 0.0, 0.0, 0.0};
    // w (n . A) v and its momentum to the sums s
    auto add = [&](const int r, const Real ra, const Real v, Real s[4]) {
      const Real e0 = dir_(r,0), e1 = dir_(r,1), e2 = dir_(r,2);
      const Real q = dir_(r,3)*((e0*(ra*a0) + e1*(ra*a1) + e2*(ra*a2))*v);
      s[0] += q;
      s[1] += q*e0;
      s[2] += q*e1;
      s[3] += q*e2;
    };
    for (int r = 0; r < n; ++r) {
      const Real xv = rg_(b,SPFS_XV*n + r,kc,jc);
      if (xv < 0.0) {continue;}
      add(r, rr, xv, sr);
      const Real xe = xb ? rg_(b,SPFS_XE*n + r,kc,jc) : -1.0;
      if (xe >= 0.0) {add(r, rx, xe, sx);}
    }
    for (int c = 0; c < 4; ++c) {
      cf_(m,c,kc,jc,ir-is) += sr[c];
      if (xb) {cf_(m,c,kc,jc,ix-is) += sx[c];}
    }
  });
}

}  // namespace radvet
