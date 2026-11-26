//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_sph_sweep.cpp
//! \brief <rad_vet>/vet_gd: the Eddington tensor on the spherical-polar grid from a
//! short-characteristics formal solution along globally fixed (Cartesian) directions.
//!
//! The grey equation n . grad I = chi (S - I) is swept shell by shell for a fixed set of
//! straight-ray directions (HEALPix or product set; weights integrate 1, n, nn exactly).
//! Two passes, one launch per shell, one thread per (block, cell, direction):
//!   inward (n . r_hat < 0), top-down: upwind point on shell i+1 (vacuum at the top);
//!   outward (n . r_hat >= 0), bottom-up: shell i-1, or beyond the tangent point on the
//!     same shell (inward value), or the inner face (diffusion intensity).
//! Upwind I, ln chi, ln S are bilinear in (theta, phi) index space; S is linear in tau.
//! Periodic theta/phi walls: ghost intensities are re-indexed by the local-frame map
//! n -> M n (nearest set direction). Block boundaries: lagged inflow (block-Jacobi).
//! Output: the intensities vgd_i (ClosureFromSc, RayFaceSums).
//!
//!  called from: rad_vet_closure_ring.cpp, rays/sc_build.cpp
//!  calls into: rad_vet_params_closure.cpp, rays/sc_bandring.cpp, rays/sc_cart_setup.cpp,
//!    rays/sc_cs_sweep.cpp, rays/sc_setup.cpp, rays/sc_sph_globalband.cpp,
//!    rays/sc_sph_halo.cpp, rays/sc_sph_halo_exact.cpp, rays/sc_sph_halo_fused.cpp,
//!    rays/sc_sph_walls.cpp, rays/sc_tables.cpp

#include <math.h>

#include <algorithm>
#include <array>
#include <type_traits>
#include <vector>

#include "athena.hpp"
// FMA contraction only within one expression in this TU
#if defined(__clang__)
#pragma clang fp contract(on)
#endif
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "coordinates/coordinates.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rays/sc_segment.hpp"
#include "rad_vet/rays/sc_ray.hpp"
#include "rad_vet/rays/sc_halo.hpp"
#include "rad_vet/rays/sc_cs.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScSweep
//! \brief one sweep (inward then outward) of every direction

void RadiationVet::ScSweep() {
  if (vgd_csg) {   // the cubed-sphere sweep (rays/sc_cs_sweep.cpp)
    ScCsSweep();
    return;
  }
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  // band index space of vgd_i / vgd_cs / the tables: active cells at [wb, wb + nx)
  const int wb = vgd_w;
  const int js = wb, je = wb + indcs.nx2 - 1;
  const int ks = wb, ke = wb + indcs.nx3 - 1;
  const int n1 = indcs.nx1;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n = vgd_n;
  const int lcut = vgd_scut;
  auto cs_ = vgd_cs;
  auto vi_ = vgd_i;
  auto dir_ = vgd_dir;
  auto cnt_ = vray_cnt;
  auto mr_ = vgd_mr;
  auto x1f_ = vgd_x1f;
  // outside the block: on-rank neighbours are read from their interior (re-indexed at
  // a wall), remote ones from the band filled per shell by ScHalo (vgd_raw: owner or
  // band, re-indexed at the read)
  auto hl_ = vgd_hloc;
  const bool raw_ = vgd_raw;
  auto own_ = vgd_own;
  auto wl_ = vgd_wall;
  auto mp_ = vgd_map;
  const int nx2b = indcs.nx2, nx3b = indcs.nx3;
  const bool wint = vgd_wint;
  auto m3_ = vgd_m3;
  auto w3_ = vgd_w3;
  // the wall kernel (bodym) for mirror, vacuum and extrap theta walls
  const bool thm = vgd_thmir || (vgd_thvac > 0) ||
                   (vgd_slab == 1 && (vet_x2open || vet_x3open));
  const bool bfc = vet_bc_face;   // vet_bc_diff = face
  const bool imir = (ibc_x1min == VET_IBC_REFLECT) && lcut == 0;
  constexpr bool spl = true;   // (the segment rule: split, emission and jump terms)
  constexpr int spe = 3;
  constexpr int vim = 3;   // (corner at a material edge, log elsewhere)
  // inner emission at the physical inner face (icut = 0): black sphere of flux F_in
  // (vet_bc_inner, VetBcEmit: emit, or auto at a Marshak face)
  const bool emi = VetBcEmit() && (lcut == 0) && ibc_x1min != VET_IBC_REFLECT;
  const Real iemit = VetBathLo();
  // top bath (vet_bc_bath, Marshak outer x1, implicit_ebath_x1max != 0): carried by the
  // top shell's inward rays, attenuated from the top face (as vet_sc)
  const bool tbh = vgd_tbath && (ibc_x1max == VET_IBC_MARSHAK) && (iebath_x1max != 0.0);
  const Real itop = tbh ? iebath_x1max : 0.0;
  const Real mt0 = pm->mesh_size.x2min, mtl = pm->mesh_size.x2max - mt0;
  auto cx1v = (vgd_slab == 1) ? vgd_sx1v : pmy_pack->pcoord->x1v;
  auto cx1f = (vgd_slab == 1) ? vgd_sx1f : pmy_pack->pcoord->xx1f;
  auto &mbsize = pmy_pack->pmb->mb_size;
  const ScLG lg_ = vgd_lg;
  const int jlo = 0, jhi = je + wb - 1;
  const int klo = 0, khi = ke + wb - 1;
  const Real twopi = 2.0*M_PI;
  // stretched theta / full sphere (sph2.inc): pol_ 2 = vet_gd_band_global (cell-window
  // clamp), 1 = block band
  const bool sth_ = vgd_sth;
  const int pol_ = vgd_pol ? (vgd_gband ? 2 : 1) : 0;
  const bool gphb_ = vgd_gband; auto phcb_ = vgd_phc; const int gpwb_ = vgd_w;
  auto thn_ = vgd_thn;
  auto hp_ = vgd_hpol;
  auto wsd_ = vgd_wsd;
  const int c2w_ = indcs.nx2 + 2*wb;
  const int xr_ = vgd_xr;   // reads beyond the 2x2 stencil (rad_vet.hpp)
  // Cartesian slab (ScSlabSelect): x1 plays r, constant foot offsets per direction;
  // x1 periodic: the pass's entry (ScSweepX1)
  const bool slb_ = (vgd_slab == 1);
  const int slbc_ = slb_ ? ((indcs.nx3 > 1) ? 2 : 1) : 0;   // (ScCand)
  const int x1p_ = slb_ ? vgd_x1p : 0;
  const bool x1si_ = vet_x1src;
  // the fused sweep counts as two sweeps, a list-only preparation as none
  if (!vgd_xff_prep) {vgd_hsweep += vgd_xff_now ? 2 : 1;}
  // vet_gd_shell_list: per pass the (m, k, j, d) whose n . r_hat (face midpoints) lies
  // on the pass's branch, widened by 1e-12 (the kernel's test still decides); remade
  // when the set changes
  if (vgd_shl_alpha != vgd_alpha) {
    const int nkj = (ke - ks + 1)*(je - js + 1), nj = je - js + 1;
    const int ntot = (nmb1 + 1)*nkj*n;
    for (int p = 0; p < 2; ++p) {
      const bool inw = (p == 0);
      if (static_cast<int>(vgd_shl[p].extent(0)) < ntot) {
        Kokkos::realloc(vgd_shl[p], ntot);
      }
      auto L_ = vgd_shl[p];
      auto mo_ = vgd_x1mo;
      int cnt = 0;
      const bool gph_ = vgd_gband; auto phc_ = vgd_phc; const int gpw_ = vgd_w;
      Kokkos::parallel_scan("vet_vgd_shl", Kokkos::RangePolicy<>(vgd_cur, 0, ntot),
      KOKKOS_LAMBDA(const int c, int &acc, const bool fin) {
        int t = c;
        const int d = t % n;
        t /= n;
        const int m = mo_(t/nkj);
        t -= (t/nkj)*nkj;
        const int k = ks + t/nj;
        const int j = js + (t % nj);
        const Real nx = dir_(d,0), ny = dir_(d,1), nz = dir_(d,2);
        const Real th = sth_ ? thn_(m,j)
                             : lg_.Th(mbsize, m, j, js);
        const Real ph = (gph_ ? phc_(m, gpw_ + k - ks) : lg_.Ph(mbsize, m, k, ks));
        const Real st = sin(th), ct = cos(th), sp = sin(ph), cp = cos(ph);
        const Real mr = slb_ ? nx : (nx*st*cp + ny*st*sp + nz*ct);
        const bool keep = inw ? (mr < 1.0e-12) : (mr >= -1.0e-12);
        if (keep) {
          if (fin) {L_(acc) = (m*nkj + t)*n + d;}
          acc += 1;
        }
      }, cnt);
      vgd_shn[p] = cnt;
      // several MeshBlocks along x1: the code range of each level
      vgd_shlv[p].assign(vgd_nb1 + 1, cnt);
      vgd_shlv[p][0] = 0;
      if (vgd_nb1 > 1) {
        auto lh = Kokkos::create_mirror_view_and_copy(HostMemSpace(), L_);
        for (int q = cnt - 1; q >= 0; --q) {vgd_shlv[p][vgd_x1lv[lh(q)/(n*nkj)]] = q;}
        for (int lv = vgd_nb1 - 1; lv > 0; --lv) {
          vgd_shlv[p][lv] = std::min(vgd_shlv[p][lv], vgd_shlv[p][lv + 1]);
        }
      }
    }
    vgd_shl_alpha = vgd_alpha;
  }
  // per-cell sin / cos of the face-midpoint (theta, phi) for the shell kernel (made
  // once; bitwise the kernel's own evaluation)
  if (static_cast<int>(vgd_trig.extent(0)) != nmb1 + 1) {
    Kokkos::realloc(vgd_trig, nmb1 + 1, ke - ks + 1, je - js + 1, 4);
    auto tw_ = vgd_trig;
    const bool gph_ = vgd_gband; auto phc_ = vgd_phc; const int gpw_ = vgd_w;
    par_for("vet_vgd_trig", vgd_cur, 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      const Real th = sth_ ? thn_(m,j)
                           : lg_.Th(mbsize, m, j, js);
      const Real ph = (gph_ ? phc_(m, gpw_ + k - ks) : lg_.Ph(mbsize, m, k, ks));
      // Cartesian slab: r_hat = x1_hat (theta = pi/2, phi = 0 exactly)
      tw_(m,k-ks,j-js,0) = slb_ ? 1.0 : sin(th);
      tw_(m,k-ks,j-js,1) = slb_ ? 0.0 : cos(th);
      tw_(m,k-ks,j-js,2) = slb_ ? 0.0 : sin(ph);
      tw_(m,k-ks,j-js,3) = slb_ ? 1.0 : cos(ph);
    });
  }
  auto tg_ = vgd_trig;
  // the exact lists of every (pass, shell) whose record does not match (new set, cut or
  // depth), in a few batched collectives (ScHxBatch)
  if (vgd_ipc_on && vgd_ipc_state == 0 && vgd_hx_on && vgd_hmpi && (vgd_hcomp > 0)) {
    const size_t cs = std::max<size_t>(vgd_crb.extent(0), 1);
    ScIpcInit(vgd_prk, cs, 2*cs);   // vet_gd_halo_ipc (collective, once)
  }
  if (vgd_hx_on && vgd_hmpi && (vgd_hcomp > 0)) {
    const int n1c = indcs.nx1 + 2*indcs.ng;
    if (static_cast<int>(vgd_hxe.size()) != 2*n1c) {vgd_hxe.resize(2*n1c);}
    std::vector<std::array<int, 4>> jobs;
    const int nshb = n1 - lcut;
    for (int pass = 0; pass < 2; ++pass) {
      const bool inw = (pass == 0);
      for (int q = 0; q < nshb; ++q) {
        const int l = inw ? (n1 - 1 - q) : (lcut + q);
        const int i = is + l;
        const int wsx = inw ? vgd_wsi[i] : vgd_wso[i];
        if (wsx < 1 || wsx > vgd_w) {continue;}
        if (!ScHxKey(pass*n1c + i, wsx)) {jobs.push_back({pass, i, l, wsx});}
      }
    }
    ScHxBatch(jobs);
  }
  // a list-only preparation of this set (shell lists, exact halo lists)
  if (vgd_xff_prep) {return;}
  const bool xffs = vgd_xff_now;
  if (ibc_x1max == VET_IBC_REFLECT) {ScTopMirror();}
  // vet_gd_band_ring: shells are lent a wide slab while their band is needed (the
  // other set's array as the second one, ScRingFlush)
  ScRingSweepBegin(xffs);
  const bool rga = vgd_ring_act;
  for (int pass = 0; pass < 2; ++pass) {
    const bool inw = (pass == 0);
    if (imir && !inw) {ScInnerMirror();}
    auto isc_ = vgd_isc;
#define VET_VGD_BODY_NAME body
#define VET_VGD_BODY_MIR 0
#include "rad_vet/rays/sc_sph_sweep_body.inc"
#undef VET_VGD_BODY_NAME
#undef VET_VGD_BODY_MIR
#define VET_VGD_BODY_NAME bodym
#define VET_VGD_BODY_MIR 1
#include "rad_vet/rays/sc_sph_sweep_body.inc"
#undef VET_VGD_BODY_NAME
#undef VET_VGD_BODY_MIR
#if defined(KOKKOS_ENABLE_CUDA)
    // CUDA gx-false kernels (rays/sc_sph_sweep_body.inc): bodyg / bodymg compile in
    // vpa = gx = false (nvcc has no generic device lambda for std::false_type)
#define VET_VGD_BODY_GX0
#define VET_VGD_BODY_NAME bodyg
#define VET_VGD_BODY_MIR 0
#include "rad_vet/rays/sc_sph_sweep_body.inc"
#undef VET_VGD_BODY_NAME
#undef VET_VGD_BODY_MIR
#define VET_VGD_BODY_NAME bodymg
#define VET_VGD_BODY_MIR 1
#include "rad_vet/rays/sc_sph_sweep_body.inc"
#undef VET_VGD_BODY_NAME
#undef VET_VGD_BODY_MIR
#undef VET_VGD_BODY_GX0
#endif
    // shell bodies of the fused two-set kernel (run_xff), one per set, with its tables
    // and intensities bound to the body's names. Both come from one lambda expression
    // (one closure type), so the kernel holds one copy of the body and stack frame
#if defined(KOKKOS_ENABLE_CUDA)
#define VET_VGD_BODY_GX0
#endif
    using ScDirV = decltype(vgd_dir);
    using ScWallV = decltype(vgd_wall);
    using ScMapV = decltype(vgd_map);
    using ScMrV = decltype(vgd_mr);
    using ScM3V = decltype(vgd_m3);
    using ScW3V = decltype(vgd_w3);
    using ScIV = decltype(vgd_i);
    auto mkx0 = [&](const ScDirV dir_, const ScWallV wl_, const ScMapV mp_,
                    const ScMrV mr_, const ScM3V m3_, const ScW3V w3_,
                    const ScIV vi_) {
#define VET_VGD_BODY_NAME bodyx
#define VET_VGD_BODY_MIR 0
#include "rad_vet/rays/sc_sph_sweep_body.inc"
#undef VET_VGD_BODY_NAME
#undef VET_VGD_BODY_MIR
      return bodyx;
    };
    auto mkx1 = [&](const ScDirV dir_, const ScWallV wl_, const ScMapV mp_,
                    const ScMrV mr_, const ScM3V m3_, const ScW3V w3_,
                    const ScIV vi_) {
#define VET_VGD_BODY_NAME bodyx
#define VET_VGD_BODY_MIR 1
#include "rad_vet/rays/sc_sph_sweep_body.inc"
#undef VET_VGD_BODY_NAME
#undef VET_VGD_BODY_MIR
      return bodyx;
    };
#if defined(KOKKOS_ENABLE_CUDA)
#undef VET_VGD_BODY_GX0
#endif
    // vet_gd_body_global: this pass's closure copied once to device memory; launches
    // carry only its address (bitwise)
    using ScBodyT = decltype(body);
    using ScBodymT = decltype(bodym);
    const ScBodyT *bdp = nullptr;
    const ScBodymT *bmp = nullptr;
    if (thm) {
      bmp = ScClosureToDevice(bodym, vgd_bbuf1);
    } else {
      bdp = ScClosureToDevice(body, vgd_bbuf0);
    }
#if defined(KOKKOS_ENABLE_CUDA)
    using ScBodygT = decltype(bodyg);
    using ScBodymgT = decltype(bodymg);
    const ScBodygT *bgp = nullptr;
    const ScBodymgT *bgq = nullptr;
    if (bmp != nullptr) {bgq = ScClosureToDevice(bodymg, vgd_bbuf3);}
    if (bdp != nullptr) {bgp = ScClosureToDevice(bodyg, vgd_bbuf2);}
#endif
    // both sets' gx-false closures in device memory (ScXffOk)
    using ScXb0T = decltype(mkx0(vgd_dir, vgd_wall, vgd_map, vgd_mr, vgd_m3, vgd_w3,
                                 vgd_i));
    using ScXb1T = decltype(mkx1(vgd_dir, vgd_wall, vgd_map, vgd_mr, vgd_m3, vgd_w3,
                                 vgd_i));
    const ScXb0T *xa0 = nullptr, *xq0 = nullptr;
    const ScXb1T *xa1 = nullptr, *xq1 = nullptr;
    if (xffs) {
#if defined(KOKKOS_ENABLE_CUDA)
      // the current set's closure: the buffer of its gx-false body
      auto &bba = vgd_bbuf2;
#else
      auto &bba = vgd_bbuf0;
#endif
      if (thm) {
        xa1 = ScClosureToDevice(mkx1(vgd_dir, vgd_wall, vgd_map, vgd_mr, vgd_m3, vgd_w3,
                                     vgd_i), bba);
        xq1 = ScClosureToDevice(mkx1(vgd_tabd.dir, vgd_tabd.wall, vgd_tabd.map,
                                     vgd_tabd.mr, vgd_tabd.m3, vgd_tabd.w3, vgd_i_xb),
                                vgd_bbufx);
      } else {
        xa0 = ScClosureToDevice(mkx0(vgd_dir, vgd_wall, vgd_map, vgd_mr, vgd_m3, vgd_w3,
                                     vgd_i), bba);
        xq0 = ScClosureToDevice(mkx0(vgd_tabd.dir, vgd_tabd.wall, vgd_tabd.map,
                                     vgd_tabd.mr, vgd_tabd.m3, vgd_tabd.w3, vgd_i_xb),
                                vgd_bbufx);
      }
    }
    const int nkj = (ke - ks + 1)*(je - js + 1), nj = je - js + 1;
    // fused shell kernel of (i, l): this pass's codes of the current set, then the
    // other set's, one launch; each thread runs its code's set (bitwise two launches)
    auto run_xff = [&](const int i, const int l) {
      const int nA = vgd_shn[pass], nB = vgd_shn_xb[pass];
      if (nA + nB <= 0) {return;}
      auto LA_ = vgd_shl[pass];
      auto LB_ = vgd_shl_xb[pass];
      const std::false_type f_;
      // current set's codes padded to a multiple of 1024: every block (power of two
      // <= 1024) runs one set only
      const int nP = ((nA + 1023)/1024)*1024;
      if (thm) {
        const auto *qa = xa1;
        const auto *qb = xq1;
        Kokkos::parallel_for("vet_vgd_shell_xff", ScShellPol(vgd_cur, 0, nP + nB),
        KOKKOS_LAMBDA(const int c) {
          const bool sb = (c >= nP);
          if (!sb && c >= nA) {return;}
          int t = sb ? LB_(c - nP) : LA_(c);
          const int d = t % n;
          t /= n;
          const int m = t/nkj;
          t -= m*nkj;
          const int k = ks + t/nj;
          const int j = js + (t % nj);
          (*(sb ? qb : qa))(m, k, j, d, i, l, f_, f_);   // one body: one call site
        });
      } else {
        const auto *qa = xa0;
        const auto *qb = xq0;
        Kokkos::parallel_for("vet_vgd_shell_xff", ScShellPol(vgd_cur, 0, nP + nB),
        KOKKOS_LAMBDA(const int c) {
          const bool sb = (c >= nP);
          if (!sb && c >= nA) {return;}
          int t = sb ? LB_(c - nP) : LA_(c);
          const int d = t % n;
          t /= n;
          const int m = t/nkj;
          t -= m*nkj;
          const int k = ks + t/nj;
          const int j = js + (t % nj);
          (*(sb ? qb : qa))(m, k, j, d, i, l, f_, f_);   // one body: one call site
        });
      }
    };
    // the shell kernel over codes cl(c), c in [c0, c1) (vet_gd_shell_list)
    auto run_list = [&](const DvceArray1D<int> &cl, const int c0, const int c1,
                        const int i, const int l, const char *name) {
      if (c1 <= c0) {return;}
      auto L_ = cl;
      // vet_gd_body_global: same codes, closure through a pointer; one kernel per
      // variant
      const std::false_type f_;
#if defined(KOKKOS_ENABLE_CUDA)
      // on CUDA the gx-false closures bodyg / bodymg
      const ScBodygT *bp_ = bgp;
      const ScBodymgT *bq_ = bgq;
#else
      const ScBodyT *bp_ = bdp;
      const ScBodymT *bq_ = bmp;
#endif
#if defined(KOKKOS_ENABLE_CUDA)
      // launch-bounds variants of the CUDA gx-false kernel (VET_VGD_CUDA_LB; 0 = the
      // ScShellPol launch): register allocation and block size only (same bits)
#define VET_VGD_GF_CU(...)                                                    \
      if (bq_ != nullptr) {                                               \
        Kokkos::parallel_for(name, __VA_ARGS__, KOKKOS_LAMBDA(const int c) { \
          int t = L_(c);                                                  \
          const int d = t % n;                                            \
          t /= n;                                                         \
          const int m = t/nkj;                                            \
          t -= m*nkj;                                                     \
          const int k = ks + t/nj;                                        \
          const int j = js + (t % nj);                                    \
          (*bq_)(m, k, j, d, i, l, f_, f_);                               \
        });                                                               \
      } else {                                                            \
        Kokkos::parallel_for(name, __VA_ARGS__, KOKKOS_LAMBDA(const int c) { \
          int t = L_(c);                                                  \
          const int d = t % n;                                            \
          t /= n;                                                         \
          const int m = t/nkj;                                            \
          t -= m*nkj;                                                     \
          const int k = ks + t/nj;                                        \
          const int j = js + (t % nj);                                    \
          (*bp_)(m, k, j, d, i, l, f_, f_);                               \
        });                                                               \
      }
      switch (VET_VGD_CUDA_LB) {
        case 1: VET_VGD_GF_CU(ScLbPol<128, 4>(vgd_cur, c0, c1)); break;
        case 2: VET_VGD_GF_CU(ScLbPol<128, 5>(vgd_cur, c0, c1)); break;
        case 3: VET_VGD_GF_CU(ScLbPol<128, 6>(vgd_cur, c0, c1)); break;
        case 4: VET_VGD_GF_CU(ScLbPol<256, 2>(vgd_cur, c0, c1)); break;
        case 5: VET_VGD_GF_CU(ScLbPol<256, 3>(vgd_cur, c0, c1)); break;
        default: VET_VGD_GF_CU(ScShellPol(vgd_cur, c0, c1)); break;
      }
#undef VET_VGD_GF_CU
#else
      if (bq_ != nullptr) {
        Kokkos::parallel_for(name, ScShellPol(vgd_cur, c0, c1),
        KOKKOS_LAMBDA(const int c) {
          int t = L_(c);
          const int d = t % n;
          t /= n;
          const int m = t/nkj;
          t -= m*nkj;
          const int k = ks + t/nj;
          const int j = js + (t % nj);
          (*bq_)(m, k, j, d, i, l, f_, f_);
        });
      } else {
        Kokkos::parallel_for(name, ScShellPol(vgd_cur, c0, c1),
        KOKKOS_LAMBDA(const int c) {
          int t = L_(c);
          const int d = t % n;
          t /= n;
          const int m = t/nkj;
          t -= m*nkj;
          const int k = ks + t/nj;
          const int j = js + (t % nj);
          (*bp_)(m, k, j, d, i, l, f_, f_);
        });
      }
#endif
    };
    // interior / edge split of a shell (vgd_mid): edge = false runs the codes cl(c),
    // c in [c0, c1), deferring those whose stencil reaches the band to the list vgd_edl;
    // edge = true runs that list (the codes are independent: same bits as one launch)
    auto run_split = [&](const DvceArray1D<int> &cl, const int c0, const int c1,
                         const int i, const int l, const bool edge,
                         const DevExeSpace &ex) {
      if (c1 <= c0) {return;}
      if (vgd_edl.extent_int(0) < c1 - c0) {Kokkos::realloc(vgd_edl, c1 - c0);}
      if (!vgd_edn.is_allocated()) {
        vgd_edn = Kokkos::View<int, DevMemSpace>("vet_vgd_edn");
      }
      auto L_ = cl;
      auto E_ = vgd_edl;
      auto En_ = vgd_edn;
      const std::false_type f_;
#if defined(KOKKOS_ENABLE_CUDA)
      const ScBodygT *bp_ = bgp;
      const ScBodymgT *bq_ = bgq;
#else
      const ScBodyT *bp_ = bdp;
      const ScBodymT *bq_ = bmp;
#endif
      auto item = KOKKOS_LAMBDA(const int t0, int &m, int &k, int &j, int &d) {
        int t = t0;
        d = t % n;
        t /= n;
        m = t/nkj;
        t -= m*nkj;
        k = ks + t/nj;
        j = js + (t % nj);
      };
      if (edge) {
        Kokkos::parallel_for("vet_vgd_shell_edge", ScShellPol(ex, 0, c1 - c0),
        KOKKOS_LAMBDA(const int c) {
          if (c >= En_()) {return;}
          int m, k, j, d;
          item(E_(c), m, k, j, d);
          if (bq_ != nullptr) {
            (*bq_)(m, k, j, d, i, l, f_, f_);
          } else {
            (*bp_)(m, k, j, d, i, l, f_, f_);
          }
        });
        return;
      }
      Kokkos::deep_copy(ex, vgd_edn, 0);
      Kokkos::parallel_for("vet_vgd_shell_int", ScShellPol(ex, c0, c1),
      KOKKOS_LAMBDA(const int c) {
        int m, k, j, d;
        item(L_(c), m, k, j, d);
        const bool df = (bq_ != nullptr) ? (*bq_)(m, k, j, d, i, l, true, f_)
                                         : (*bp_)(m, k, j, d, i, l, true, f_);
        if (df) {E_(Kokkos::atomic_fetch_add(&En_(), 1)) = L_(c);}
      });
    };
    const int nsh = n1 - lcut;
    // vet_gd_halo_exact: read marks of shell i's band in this pass (when its record
    // does not match; same decision on every rank): the ghosts read by inward shell i-1
    // and outward shell i's turned rays (inward exchange) or outward shell i+1 (outward)
    const int c2b = nx2b + 2*wb, c3b = nx3b + 2*wb;
    auto hx_mark = [&](const int i, const int l, const int wsx) -> bool {
      const int n1c = indcs.nx1 + 2*indcs.ng;
      if (static_cast<int>(vgd_hxe.size()) != 2*n1c) {vgd_hxe.resize(2*n1c);}
      if (ScHxKey(pass*n1c + i, wsx)) {return false;}
      const size_t nall = static_cast<size_t>(nmb1 + 1)*n*c3b*c2b;
      // stamped marks (no clearing): a ghost is marked when it holds the build's stamp
      if (vgd_hx_mk.extent(0) < nall) {
        Kokkos::realloc(vgd_hx_mk, nall);
        vgd_hx_stamp = 0;
      }
      if (!vgd_hx_oc.is_allocated()) {
        vgd_hx_oc = Kokkos::View<int, DevMemSpace>("vet_vgd_hx_oc");
      }
      vgd_hx_stamp += 1;
      const int stamp = vgd_hx_stamp;
      Kokkos::deep_copy(vgd_cur, vgd_hx_oc, 0);
      auto oc_ = vgd_hx_oc;
      auto mk_ = vgd_hx_mk;
      auto mark = [&](const int pr, const int ir, const int lr) {
        const bool inwr = (pr == 0);
        auto L_ = vgd_shl[pr];
        Kokkos::parallel_for("vet_vgd_hx_mark", Kokkos::RangePolicy<>(vgd_cur, 0,
                                                                      vgd_shn[pr]),
        KOKKOS_LAMBDA(const int c) {
          int t = L_(c);
          const int d = t % n;
          t /= n;
          const int m = t/nkj;
          t -= m*nkj;
          const int k = ks + t/nj;
          const int j = js + (t % nj);
          ScCand(m, k, j, d, ir, lr, inwr, n1, lcut, js, ks, jlo, jhi, klo, khi, dir_,
                 mbsize, lg_, cx1v, cx1f, sth_, pol_, thn_, c2w_, imir, xr_, slbc_,
          [&](const int iu, const int ka, const int kb, const int ja, const int jb,
              const int dr) {
            if (iu != i) {return;}
            for (int kk = ka; kk <= kb; ++kk) {
              for (int jj = ja; jj <= jb; ++jj) {
                const int ek = (kk < wb) ? -1 : ((kk >= wb + nx3b) ? 1 : 0);
                const int ej = (jj < wb) ? -1 : ((jj >= wb + nx2b) ? 1 : 0);
                if (ek == 0 && ej == 0) {continue;}
                if (ej != 0 && wl_(m,kk,jj) >= 2) {
                 // the main sweep reads the mirror row (flag 2) or the edge row
                 // (extrap, 4) through the phi-edge slot; vacuum (3) reads nothing.
                 // Flags 3, 4: the twin's periodic reads below too
                 const int wf = wl_(m,kk,jj);
                 const int o2 = 3*(ek + 1) + 1;
                 if (wf != 3 && ek != 0 && hl_(8*m + ((o2 < 4) ? o2 : (o2 - 1))) < 0) {
                    const bool hi = (ej > 0);
                    const int g = hi ? (jj - wb - nx2b) : (wb - 1 - jj);
                    const int jm = (wf == 4) ? (hi ? (wb + nx2b - 1) : wb)
                                             : ScMirRow(g, hi, wb, nx2b);
                    for (int t = 0; t < ((wint && wf == 2) ? 3 : 1); ++t) {
                      if (t > 0 && w3_(m,kk,jj,dr,t) == 0.0) {continue;}
                      const int dq = (wf == 4) ? dr : (wint ? m3_(m,kk,jj,dr,t)
                                                           : mp_(m,kk,jj,dr));
                      mk_(((m*n + dq)*c3b + kk)*c2b + jm) = stamp;
                    }
                    if (kk < wb - wsx || kk >= wb + nx3b + wsx) {
                      Kokkos::atomic_add(&oc_(), 1);
                    }
                 }
                 if (wf == 2) {continue;}
                }
                const int oo = 3*(ek + 1) + (ej + 1);
                if (hl_(8*m + ((oo < 4) ? oo : (oo - 1))) >= 0) {continue;}
                mk_(((m*n + dr)*c3b + kk)*c2b + jj) = stamp;
                if (jj < wb - wsx || jj >= wb + nx2b + wsx || kk < wb - wsx ||
                    kk >= wb + nx3b + wsx) {
                 Kokkos::atomic_add(&oc_(), 1);
                }
              }
            }
          });
        });
      };
      if (inw) {
        if (l - 1 >= lcut) {mark(0, i - 1, l - 1);}
        mark(1, i, l);
      } else if (l + 1 <= n1 - 1) {
        mark(1, i + 1, l + 1);
      }
      return true;
    };
    const bool hxm = vgd_hx_on && vgd_hmpi && (vgd_hcomp > 0);
    // several MeshBlocks along x1: levels in the pass's order, each from the boundary
    // shell handed on (wall maps couple all directions: no wavefront)
    const int nst = vgd_nb1;
    // remote partners over MPI: the interior / edge split (vgd_mid), so the next shell's
    // interior runs while this shell's messages travel
    // IPC (CUDA): the interior runs on the second stream vgd_es2, concurrent with the
    // exchange kernel on the sweep stream (events order the two)
#if defined(KOKKOS_ENABLE_CUDA)
    const bool pipc = vgd_hmpi && vgd_ipc_state == 1 && !xffs && nst == 1;
#else
    const bool pipc = false;
#endif
    const bool ppl = vgd_hmpi && (vgd_ipc_state != 1 || pipc) && !xffs && nst == 1;
    if (pipc) {ScEs2Make();}
    ScCsEv ego, eint;   // interior may start (sweep stream) / interior done (vgd_es2)
    for (int qq = 0; qq < nst*nsh; ++qq) {
      const int q = qq % nsh, lv = inw ? (nst - 1 - qq/nsh) : qq/nsh;
      if (nst > 1 && q == 0 && qq > 0) {ScX1Hand(!inw);}
      const int l = inw ? (n1 - 1 - q) : (lcut + q);
      const int i = is + l;
      // the shell exchanged after this one (its mask is prepared while this shell's
      // messages travel)
      if (q + 1 < n1 - lcut) {
        const int in = is + (inw ? (n1 - 2 - q) : (lcut + q + 1));
        vgd_hnext[0] = in;
        vgd_hnext[1] = inw ? 1 : 0;
        vgd_hnext[2] = inw ? vgd_wsi[in] : vgd_wso[in];
      } else if (inw) {
        const int in = is + lcut;
        vgd_hnext[0] = in;
        vgd_hnext[1] = 0;
        vgd_hnext[2] = vgd_wso[in];
      } else {
        vgd_hnext[0] = -1;
      }
      if (rga) {
        // this shell's wide slab (queued during the previous exchange, else now):
        // inward from scratch, outward with its interior and inward band
        ScRingFlush();
        if (vgd_rg_state[i] < 0) {
          ScRingOp o;
          o.kind = 3;
          o.i = i;
          o.copyin = !inw;
          o.restore = !inw;
          vgd_rg_pend.push_back(o);
        }
        ScRingFlush();
      }
      if (xffs) {
        // both cross-fade sets, one launch
        run_xff(i, l);
      } else if (ppl && q > 0) {
        // the interior ran behind the previous exchange: the deferred edge items
        if (pipc) {eint.StreamWait(vgd_cur);}
        run_split(vgd_shl[pass], vgd_shlv[pass][lv], vgd_shlv[pass][lv + 1], i, l, true,
                  vgd_cur);
      } else {
        // shell list: only this pass's (and level's) (m, k, j, d)
        run_list(vgd_shl[pass], vgd_shlv[pass][lv], vgd_shlv[pass][lv + 1], i, l,
                 "vet_vgd_shell");
      }
      if (rga) {
        // the previous shell has been read (stream order): save its inward band
        // (inward pass), its interior back to the compact slab, lend the next shell a
        // slab; in this shell's exchange hook (ScHlBegin) or before the next kernel
        if (q > 0) {
          const int ip = inw ? (i + 1) : (i - 1);
          ScRingOp o;
          if (inw) {
            o.kind = 1;
            o.i = ip;
            vgd_rg_pend.push_back(o);
          }
          o.kind = 2;
          o.i = ip;
          vgd_rg_pend.push_back(o);
        }
        if (q + 1 < nsh) {
          ScRingOp o;
          o.kind = 3;
          o.i = inw ? (i - 1) : (i + 1);
          o.copyin = !inw;
          o.restore = !inw;
          vgd_rg_pend.push_back(o);
        }
      }
      const int wsx = inw ? vgd_wsi[i] : vgd_wso[i];
      if (ppl && q + 1 < nsh) {
        // the next shell's interior, behind this shell's pack (after the ring ops that
        // lend it its slab)
        const int ln = inw ? (l - 1) : (l + 1);
        vgd_mid = [&, ln, lv]() {
          run_split(vgd_shl[pass], vgd_shlv[pass][lv], vgd_shlv[pass][lv + 1], is + ln,
                    ln, false, vgd_cur);
        };
        if (pipc) {
          // the ring ops (they lend the next shell its slab) first, then the interior
          // on vgd_es2 while this shell's exchange runs on the sweep stream
          ScRingFlush();
          ego.Record(vgd_cur);
          ego.StreamWait(vgd_es2);
          run_split(vgd_shl[pass], vgd_shlv[pass][lv], vgd_shlv[pass][lv + 1], is + ln,
                    ln, false, vgd_es2);
          eint.Record(vgd_es2);
          vgd_mid = nullptr;
        }
      }
      // the shell is complete on every block: its lateral band (exact, not lagged)
      Kokkos::Timer th;
      vgd_hinw = inw;
      if (hxm && hx_mark(i, l, wsx)) {
        // vet_gd_halo_exact = 2: lists from the marks alone, before the exchange
        ScHxBuild2(pass*(indcs.nx1 + 2*indcs.ng) + i, i, wsx);
      }
      if (xffs) {
        // one message per partner with both sets' exact lists (vgd_raw: ScGbHaloX);
        // else two exchanges, the other set's first (each writes only its own band)
        vgd_rg_lastx = false;
        const bool hx2 = vgd_raw ? ScGbHaloX(i, wsx)
                              : ScHlXff(i, wsx);
        if (!hx2 && rga && !vgd_raw) {
          ScFatal("vet_gd_rotate_xfade_fuse: the band ring needs the merged exchange");
        }
        if (rga && inw) {   // record of the saved inward band, as below; both sets
          vgd_rg_sk[i] = vgd_raw ? 3 : (vgd_rg_lastx ? 1 : 2);
          vgd_rg_se[i] = vgd_raw ? ScRawKey(i, true, wsx) : (vgd_rg_lastx ? i : -1);
        }
        if (!hx2) {
          ScXffSwap();
          ScHalo(vgd_i, n, i, i, wsx, true);
          ScXffSwap();
          ScHalo(vgd_i, n, i, i, wsx, true);
        }
      } else {
        vgd_rg_lastx = false;
#if defined(KOKKOS_ENABLE_CUDA)
        // IPC: the ring ops of this hook touch only the shells before and after this
        // one (an acquire never takes the slot just released), so they run before the
        // exchange, which is then the one fused kernel vet_vgd_ipc_x (same bits)
        if (rga && vgd_ipc_state == 1) {ScRingFlush();}
#endif
        ScHalo(vgd_i, n, i, i, wsx, true);
        if (rga && inw) {
          // record of the values to save from this shell's inward band: the exchange's
          // exact list (vgd_raw: what it wrote), else the whole band
          const int ex = i;
          const bool xl = vgd_rg_lastx && (vgd_rg_lastx_e == ex);
          vgd_rg_sk[i] = vgd_raw ? 3 : (xl ? 1 : 2);
          vgd_rg_se[i] = vgd_raw ? ScRawKey(i, true, wsx) : (xl ? ex : -1);
        }
      }
      if (rga) {ScRingFlush();}
      ScMidRun();   // (not run by the exchange: after it)
      // the next kernel is stream-ordered behind the unpack: no fence needed
      vgd_thalo += th.seconds();
    }
    if (rga) {
      // last shell of the pass: saved (inward) and back to its compact slab
      const int il = is + (inw ? lcut : (n1 - 1));
      if (vgd_rg_state[il] >= 0) {
        ScRingOp o;
        if (inw) {
          o.kind = 1;
          o.i = il;
          vgd_rg_pend.push_back(o);
        }
        o.kind = 2;
        o.i = il;
        vgd_rg_pend.push_back(o);
      }
      ScRingFlush();
    }
  }
  ScRingSweepEnd();
}

}  // namespace radvet
