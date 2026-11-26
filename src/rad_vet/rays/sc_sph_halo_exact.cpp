//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_sph_halo_exact.cpp
//! \brief spherical-polar rays: exact halo lists (only the ghosts the sweep reads).
//!
//!  called from: rays/sc_bandring.cpp, rays/sc_sph_halo.cpp, rays/sc_sph_halo_fused.cpp,
//!    rays/sc_sph_sweep.cpp

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include "athena.hpp"
// FMA contraction only within one expression in this TU
#if defined(__clang__)
#pragma clang fp contract(on)
#endif
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "coordinates/coordinates.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rays/sc_segment.hpp"
#include "rad_vet/rays/sc_ray.hpp"
#include "rad_vet/rays/sc_halo.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::ScHxKey
//! \brief vet_gd_halo_exact: record e (exact lists or fallback) was made for the
//! current direction set, cut and depth (same answer on every rank)

bool RadiationVet::ScHxKey(const int e, const int ws) const {
  if (e < 0 || e >= static_cast<int>(vgd_hxe.size())) {return false;}
  const ScHlEntry &X = vgd_hxe[e];
  return X.gen != -1 && X.alpha == vgd_alpha && ScHxCut(e, X) && X.ws == ws;
}

// record X of (pass, shell l) = entry e holds for the current ray cut: made for it, or
// the move crosses neither l nor l + 1 (its reads are shells l - 1 .. l + 1), so a cut
// move remakes only the lists of the shells it crosses
bool RadiationVet::ScHxCut(const int e, const ScHlEntry &X) const {
  if (X.scut == vgd_scut) {return true;}
  if (X.scut < 0 || vgd_scut < 0) {return false;}
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int l = e % (indcs.nx1 + 2*indcs.ng) - indcs.is;
  auto sd = [](const int a, const int c) {return (a > c) - (a < c);};
  return sd(l, X.scut) == sd(l, vgd_scut) && sd(l + 1, X.scut) == sd(l + 1, vgd_scut);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScHxBuild2
//! \brief vet_gd_halo_exact = 2: the exact lists of record e from the read marks
//! vgd_hx_mk alone (no masks, no dense exchange); collective. (1) sources of marked
//! ghosts, (2) their compact positions, (3) per partner the count and ascending wanted
//! offsets, (4) a global vote (no mark outside depth ws, byte budget), (5) send
//! addresses and read ghosts. Every source of a read ghost is sent.

void RadiationVet::ScHxBuild2(const int e, const int i, const int ws) {
#if MPI_PARALLEL_ENABLED
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nmb = nmb1 + 1;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int nv = vgd_n;
  const int mx = std::max(nx2, nx3)*ws*nv;
  const int c2 = nx2 + 2*w, c3 = nx3 + 2*w;
  const int np = static_cast<int>(vgd_prk.size());
  const size_t rtot = static_cast<size_t>(vgd_prdsp[ws].empty() ? 0 :
                      (vgd_prdsp[ws].back() + vgd_prcnt[ws].back()))*nv;
  ScHlEntry &E = vgd_hxe[e];
  auto hl_ = vgd_hloc;
  auto ro_ = vgd_roff;
  auto wl_ = vgd_wall;
  auto mp_ = vgd_map;
  auto m3_ = vgd_m3;
  const bool wint = vgd_wint;
  auto mk_ = vgd_hx_mk;
  const int stamp = vgd_hx_stamp;
  const DevExeSpace ex_ = vgd_cur;
  auto grow = [](DvceArray1D<int> &a, const size_t n) {
    if (a.extent(0) < n) {Kokkos::realloc(Kokkos::WithoutInitializing, a, n);}
  };
  grow(vgd_hx_need, rtot + 1);
  grow(vgd_hx_t0, rtot + 1);
  grow(vgd_hx_t1, rtot + 1);
  auto nd_ = vgd_hx_need;
  auto f0_ = vgd_hx_t0;
  auto f1_ = vgd_hx_t1;
  Kokkos::deep_copy(ex_, Kokkos::subview(nd_, std::make_pair(size_t(0), rtot + 1)), 0);
  // (1) sources of the marked ghosts
  par_for("vet_vgd_hx_need", ex_, 0, nmb1, 0, 7, 0, mx - 1,
  KOKKOS_LAMBDA(const int m, const int o, const int t) {
    if (hl_(8*m + o) >= 0) {return;}
    const int oo = (o < 4) ? o : (o + 1);
    const int dk = oo/3 - 1, dj = oo%3 - 1;
    const int jd = (dj < 0) ? (w - ws) : ((dj == 0) ? w : (w + nx2));
    const int kd = (dk < 0) ? (w - ws) : ((dk == 0) ? w : (w + nx3));
    const int jn = (dj == 0) ? nx2 : ws, kn = (dk == 0) ? nx3 : ws;
    if (t >= nv*kn*jn) {return;}
    const int v = t/(kn*jn);
    const int r1 = t - v*kn*jn;
    const int kk = r1/jn;
    const int jj = r1 - kk*jn;
    if (mk_(((m*nv + v)*c3 + kd + kk)*c2 + jd + jj) != stamp) {return;}
    const size_t rb0 = static_cast<size_t>(ro_(ws,8*m + o))*nv;
    if (wint && wl_(m,kd+kk,jd+jj) != 0) {
      for (int u = 0; u < 3; ++u) {
        nd_(rb0 + (m3_(m,kd+kk,jd+jj,v,u)*kn + kk)*jn + jj) = 1;
      }
    } else {
      const int vs = (wl_(m,kd+kk,jd+jj) != 0) ? mp_(m,kd+kk,jd+jj,v) : v;
      nd_(rb0 + (vs*kn + kk)*jn + jj) = 1;
    }
  });
  int nout = 0;           // marked ghosts outside the region
  {
    Kokkos::View<int, HostMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> h(&nout);
    Kokkos::deep_copy(ex_, h, vgd_hx_oc);
  }
  // read ghosts (non-wall f0, wall f1) over the receive layout
  par_for("vet_vgd_hx_dflag", ex_, 0, nmb1, 0, 7, 0, mx - 1,
  KOKKOS_LAMBDA(const int m, const int o, const int t) {
    if (hl_(8*m + o) >= 0) {return;}
    const int oo = (o < 4) ? o : (o + 1);
    const int dk = oo/3 - 1, dj = oo%3 - 1;
    const int jd = (dj < 0) ? (w - ws) : ((dj == 0) ? w : (w + nx2));
    const int kd = (dk < 0) ? (w - ws) : ((dk == 0) ? w : (w + nx3));
    const int jn = (dj == 0) ? nx2 : ws, kn = (dk == 0) ? nx3 : ws;
    if (t >= nv*kn*jn) {return;}
    const int v = t/(kn*jn);
    const int r1 = t - v*kn*jn;
    const int kk = r1/jn;
    const int jj = r1 - kk*jn;
    const size_t q = static_cast<size_t>(ro_(ws,8*m + o))*nv + t;
    const bool rd = mk_(((m*nv + v)*c3 + kd + kk)*c2 + jd + jj) == stamp;
    const bool wa = wint && wl_(m,kd+kk,jd+jj) != 0;
    f0_(q) = (rd && !wa) ? 1 : 0;
    f1_(q) = (rd && wa) ? 1 : 0;
  });
  // (2) compact receive positions of the needed sources and of the read ghosts
  ScScan2(ex_, nd_, rtot, nd_, 0, "vet_vgd_hx_scan2", vgd_scan_sums);
  ScScan2(ex_, f0_, rtot, f1_, rtot, "vet_vgd_hx_scan2", vgd_scan_sums);
  // per partner: needed counts (scanned need at the piece bounds), host
  std::vector<int> rpx(np + 1), spx(np + 1);
  int nrx = 0, nwx = 0;
  {
    if (static_cast<int>(vgd_hx_hb.extent(0)) < 4*np + 6) {
      vgd_hx_hb = Kokkos::View<int*, Kokkos::SharedHostPinnedSpace>("vet_vgd_hx_hb",
                                                                    4*np + 6);
    }
    auto hb_ = vgd_hx_hb;
    const int o2 = 2*np + 4;
    for (int p = 0; p < np; ++p) {hb_(o2 + p) = vgd_prdsp[ws][p]*nv;}
    hb_(o2 + np) = static_cast<int>(rtot);
    const int np_ = np;
    const size_t rt = rtot;
    Kokkos::parallel_for("vet_vgd_hx_bnd", Kokkos::RangePolicy<>(ex_, 0, np + 3),
    KOKKOS_LAMBDA(const int q) {
      if (q <= np_) {
        hb_(q) = nd_(hb_(o2 + q));
      } else if (q == np_ + 1) {
        hb_(q) = f0_(rt);
      } else {
        hb_(q) = f1_(rt);
      }
    });
    ex_.fence();
    for (int p = 0; p <= np; ++p) {rpx[p] = hb_(p);}
    nrx = hb_(np + 1);
    nwx = hb_(np + 2);
  }
  const int nneed = rpx[np];
  // wanted offsets, ascending, relative to the partner's piece start
  grow(vgd_hx_nc, nneed + 1);
  auto rq_ = vgd_hx_nc;
  {
    std::vector<int> pst(np + 1);
    for (int p = 0; p < np; ++p) {pst[p] = vgd_prdsp[ws][p]*nv;}
    pst[np] = static_cast<int>(rtot);
    if (static_cast<int>(vgd_hx_pin.extent(0)) < np + 1) {
      vgd_hx_pin = Kokkos::View<int*, Kokkos::SharedHostPinnedSpace>("vet_vgd_hx_pin",
                                                                     np + 1);
    }
    auto pst_d = vgd_hx_pin;
    for (int p = 0; p <= np; ++p) {pst_d(p) = pst[p];}
    const int np_ = np;
    Kokkos::parallel_for("vet_vgd_hx_req", Kokkos::RangePolicy<>(ex_, 0, rtot),
    KOKKOS_LAMBDA(const size_t q) {
      if (nd_(q + 1) == nd_(q)) {return;}
      int lo = 0, hi = np_ - 1;
      while (lo < hi) {             // partner piece holding q
        const int md = (lo + hi + 1)/2;
        if (static_cast<size_t>(pst_d(md)) <= q) {lo = md;} else {hi = md - 1;}
      }
      rq_(nd_(q)) = static_cast<int>(q) - pst_d(lo);
    });
    ex_.fence();
  }
  // (3) counts, then the offsets
  {
    std::vector<int> rc(np), sc(np);
    for (int p = 0; p < np; ++p) {rc[p] = rpx[p+1] - rpx[p];}
    std::vector<MPI_Request> req;
    for (int p = 0; p < np; ++p) {
      req.emplace_back();
      MPI_Irecv(&sc[p], 1, MPI_INT, vgd_prk[p], 7005, vgd_comm, &req.back());
      req.emplace_back();
      MPI_Isend(&rc[p], 1, MPI_INT, vgd_prk[p], 7005, vgd_comm, &req.back());
    }
    MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
    spx[0] = 0;
    for (int p = 0; p < np; ++p) {spx[p+1] = spx[p] + sc[p];}
  }
  const int nsx = spx[np];
  grow(vgd_hx_ns, nsx + 1);
  auto so2_ = vgd_hx_ns;
  {
    std::vector<MPI_Request> req;
    for (int p = 0; p < np; ++p) {
      req.emplace_back();
      MPI_Irecv(so2_.data() + spx[p], spx[p+1] - spx[p], MPI_INT, vgd_prk[p], 7006,
                vgd_comm, &req.back());
      req.emplace_back();
      MPI_Isend(rq_.data() + rpx[p], rpx[p+1] - rpx[p], MPI_INT, vgd_prk[p], 7006,
                vgd_comm, &req.back());
    }
    MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
  }
  // (4) the vote
  const size_t nsz[5] = {static_cast<size_t>(std::max(nsx, 1)),
                         static_cast<size_t>(std::max(nrx, 1)),
                         static_cast<size_t>(std::max(nrx, 1)),
                         static_cast<size_t>(std::max(nwx, 1)),
                         static_cast<size_t>(std::max(3*nwx, 1))};
  DvceArray1D<int> *arr[5] = {&E.sl, &E.dl, &E.sr, &E.dw, &E.sw};
  size_t have = 0, after = 0;
  for (int q = 0; q < 5; ++q) {
    have += arr[q]->extent(0);
    after += std::max(arr[q]->extent(0), nsz[q]);
  }
  int ok = (nout == 0) ? 1 : 0;
  if (after > have) {
    if (vgd_hl_cap == 0) {vgd_hl_cap = static_cast<size_t>(vgd_hl_mb) << 20;}
    if (vgd_hl_bytes + sizeof(int)*(after - have) > vgd_hl_cap) {ok = 0;}
  }
  int gok = 0;
  MPI_Allreduce(&ok, &gok, 1, MPI_INT, MPI_MIN, vgd_comm);
  vgd_hx_nout += nout;
  E.alpha = vgd_alpha;
  E.scut = vgd_scut;
  E.ws = ws;
  E.exact = false;
  E.gen = -2;
  if (gok == 0) {vgd_hx_nfall += 1.0; return;}
  if (after > have) {
    for (int q = 0; q < 5; ++q) {
      if (arr[q]->extent(0) < nsz[q]) {Kokkos::realloc(*arr[q], nsz[q]);}
    }
    vgd_hl_bytes += sizeof(int)*(after - have);
  }
  auto sl_ = E.sl;
  auto dl_ = E.dl;
  auto sr_ = E.sr;
  auto dw_ = E.dw;
  auto sw_ = E.sw;
  // (5) send addresses: offset within my piece for partner p -> slot -> (v, k, j)
  {
    // my send pieces of depth ws sorted by start: (start, 8 m + o)
    auto so_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_soff);
    std::vector<std::array<int, 2>> pc;
    for (int t = 0; t < 8*nmb; ++t) {
      if (so_h(ws,t) >= 0) {pc.push_back({so_h(ws,t)*nv, t});}
    }
    std::sort(pc.begin(), pc.end());
    const int npc = static_cast<int>(pc.size());
    if (static_cast<int>(vgd_hx_pin2.extent(0)) < 2*npc + 2*np + 2) {
      vgd_hx_pin2 = Kokkos::View<int*, Kokkos::SharedHostPinnedSpace>("vet_vgd_hx_pin2",
                                                                      2*npc + 2*np + 2);
    }
    auto ps_d = vgd_hx_pin2;
    auto ps_h = ps_d;
    for (int q = 0; q < npc; ++q) {ps_h(q) = pc[q][0]; ps_h(npc + q) = pc[q][1];}
    for (int p = 0; p <= np; ++p) {
      ps_h(2*npc + p) = spx[p];
      ps_h(2*npc + np + 1 + p) = (p < np) ? vgd_pdsp[ws][p]*nv : 0;
    }
    const int np_ = np;
    Kokkos::parallel_for("vet_vgd_hx_slist", Kokkos::RangePolicy<>(ex_, 0, nsx),
    KOKKOS_LAMBDA(const int c) {
      int lo = 0, hi = np_ - 1;     // partner of send entry c
      while (lo < hi) {
        const int md = (lo + hi + 1)/2;
        if (ps_d(2*npc + md) <= c) {lo = md;} else {hi = md - 1;}
      }
      const int q = ps_d(2*npc + np_ + 1 + lo) + so2_(c);
      int a = 0, b = npc - 1;       // piece holding q
      while (a < b) {
        const int md = (a + b + 1)/2;
        if (ps_d(md) <= q) {a = md;} else {b = md - 1;}
      }
      const int sid = ps_d(npc + a);
      const int m = sid/8, o = sid % 8;
      const int oo = (o < 4) ? o : (o + 1);
      const int dk = oo/3 - 1, dj = oo%3 - 1;
      const int jn = (dj == 0) ? nx2 : ws;
      const int kn = (dk == 0) ? nx3 : ws;
      const int t = q - ps_d(a);
      const int v = t/(kn*jn);
      const int r1 = t - v*kn*jn;
      const int kk = r1/jn;
      const int jj = r1 - kk*jn;
      const int js2 = (dj > 0) ? (w + nx2 - ws) : w;
      const int ks2 = (dk > 0) ? (w + nx3 - ws) : w;
      sl_(c) = ((m*nv + v)*c3 + ks2 + kk)*c2 + js2 + jj;
    });
  }
  // read ghosts with their compact source positions (every source is sent)
  par_for("vet_vgd_hx_rlist", ex_, 0, nmb1, 0, 7, 0, mx - 1,
  KOKKOS_LAMBDA(const int m, const int o, const int t) {
    if (hl_(8*m + o) >= 0) {return;}
    const int oo = (o < 4) ? o : (o + 1);
    const int dk = oo/3 - 1, dj = oo%3 - 1;
    const int jd = (dj < 0) ? (w - ws) : ((dj == 0) ? w : (w + nx2));
    const int kd = (dk < 0) ? (w - ws) : ((dk == 0) ? w : (w + nx3));
    const int jn = (dj == 0) ? nx2 : ws, kn = (dk == 0) ? nx3 : ws;
    if (t >= nv*kn*jn) {return;}
    const int v = t/(kn*jn);
    const int r1 = t - v*kn*jn;
    const int kk = r1/jn;
    const int jj = r1 - kk*jn;
    const size_t rb0 = static_cast<size_t>(ro_(ws,8*m + o))*nv;
    const size_t q = rb0 + t;
    const int ad = ((m*nv + v)*c3 + kd + kk)*c2 + jd + jj;
    if (f1_(q + 1) > f1_(q)) {
      const int c = f1_(q);
      dw_(c) = ad;
      for (int u = 0; u < 3; ++u) {
        sw_(3*c + u) = nd_(rb0 + (m3_(m,kd+kk,jd+jj,v,u)*kn + kk)*jn + jj);
      }
    } else if (f0_(q + 1) > f0_(q)) {
      const int vs = (wl_(m,kd+kk,jd+jj) != 0) ? mp_(m,kd+kk,jd+jj,v) : v;
      dl_(f0_(q)) = ad;
      sr_(f0_(q)) = nd_(rb0 + (vs*kn + kk)*jn + jj);
    }
  });
  ex_.fence();
  E.sp = spx;
  E.rp = rpx;
  E.dsto.clear();
  E.ns = nsx;
  E.nr = nrx;
  E.nw = nwx;
  E.exact = true;
  E.gen = ++vgd_hl_gen;
  vgd_hx_nmade += 1.0;
  vgd_hx_vexact += rpx[np];
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScHxBatch
//! \brief the exact lists (vet_gd_halo_exact = 2) of every (pass, shell) job, built up
//! front in ScSweep in groups of G shells: per group (1) device-only marks, flags,
//! scans, counts and offsets, one fence; (2) one count message per partner, the
//! offsets, one vote (MPI_Allreduce of G flags); (3) send / receive lists, one fence.
//! Same expressions as ScHxBuild2 (bitwise). G <= 64 from a 256 MB pool (same on
//! every rank); the budget is checked on the group's sum (a shell over it: dense path).

void RadiationVet::ScHxBatch(const std::vector<std::array<int, 4>> &jobs) {
#if MPI_PARALLEL_ENABLED
  if (jobs.empty()) {return;}
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nmb = nmb1 + 1;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int nv = vgd_n, n = vgd_n;
  const int c2 = nx2 + 2*w, c3 = nx3 + 2*w;
  const int np = static_cast<int>(vgd_prk.size());
  const int n1 = indcs.nx1, is = indcs.is;
  const int n1c = indcs.nx1 + 2*indcs.ng;
  const int lcut = vgd_scut;
  const bool imir_ = (ibc_x1min == VET_IBC_REFLECT) && lcut == 0;
  const int wb = w;
  const int js = wb, je = wb + nx2 - 1, ks = wb, ke = wb + nx3 - 1;
  const int jlo = 0, jhi = je + wb - 1, klo = 0, khi = ke + wb - 1;
  const bool sth_ = vgd_sth, pol_ = vgd_pol;
  auto thn_ = vgd_thn;
  const int c2w_ = indcs.nx2 + 2*wb;
  const int xr_ = vgd_xr;
  const int nkj = (ke - ks + 1)*(je - js + 1), nj = je - js + 1;
  const int nx2b = nx2, nx3b = nx3;
  auto hl_ = vgd_hloc;
  auto ro_ = vgd_roff;
  auto wl_ = vgd_wall;
  auto mp_ = vgd_map;
  auto m3_ = vgd_m3;
  auto dir_ = vgd_dir;
  auto &mbsize = pmy_pack->pmb->mb_size;
  const ScLG lg_ = vgd_lg;
  auto cx1v = pmy_pack->pcoord->x1v;
  auto cx1f = pmy_pack->pcoord->xx1f;
  const bool wint = vgd_wint;
  auto w3_ = vgd_w3;           // mirror theta walls: read marks
  const DevExeSpace ex_ = vgd_cur;
  if (static_cast<int>(vgd_hxe.size()) != 2*n1c) {vgd_hxe.resize(2*n1c);}
  Kokkos::Timer tb;
  auto rtot_of = [&](const int ws) -> size_t {
    return static_cast<size_t>(vgd_prdsp[ws].empty() ? 0 :
           (vgd_prdsp[ws].back() + vgd_prcnt[ws].back()))*nv;
  };
  auto stot_of = [&](const int ws) -> size_t {
    return static_cast<size_t>(vgd_pdsp[ws].empty() ? 0 :
           (vgd_pdsp[ws].back() + vgd_pscnt[ws].back()))*nv;
  };
  // ---- setup (once): slot capacities, G, per-depth partner / piece tables ----
  if (vgd_hxb_g == 0) {
    size_t rcap = 1, scap = 1;
    for (int ws = 1; ws <= w; ++ws) {
      rcap = std::max(rcap, rtot_of(ws) + 1);
      scap = std::max(scap, stot_of(ws) + 1);
    }
    // pool 256 MB (exact lists built in groups, only at a new direction set)
    const double mb = 256.0;
    const int gmax = 64;
    const double per = (4.0*rcap + scap)*sizeof(int);
    int g = static_cast<int>(std::floor(mb*1048576.0/per));
    g = std::max(1, std::min(g, gmax));
    MPI_Allreduce(MPI_IN_PLACE, &g, 1, MPI_INT, MPI_MIN, vgd_comm);
    vgd_hxb_g = g;
    vgd_hxb_rcap = rcap;
    vgd_hxb_scap = scap;
    const size_t G = static_cast<size_t>(g);
    Kokkos::realloc(Kokkos::WithoutInitializing, vgd_hxb_nd, G*rcap);
    Kokkos::realloc(Kokkos::WithoutInitializing, vgd_hxb_f0, G*rcap);
    Kokkos::realloc(Kokkos::WithoutInitializing, vgd_hxb_f1, G*rcap);
    Kokkos::realloc(Kokkos::WithoutInitializing, vgd_hxb_rq, G*rcap);
    Kokkos::realloc(Kokkos::WithoutInitializing, vgd_hxb_so, G*scap);
    vgd_hxb_oc = DvceArray1D<int>("vet_vgd_hxb_oc", G);
    vgd_hxb_hb = Kokkos::View<int*, Kokkos::SharedHostPinnedSpace>("vet_vgd_hxb_hb",
                                                                   G*(np + 4));
    vgd_hxb_sp = Kokkos::View<int*, Kokkos::SharedHostPinnedSpace>("vet_vgd_hxb_sp",
                                                                   G*(np + 1));
    // per depth ws: the partners' receive piece starts (np + 1, last = rtot), my send
    // pieces sorted by start (start, 8 m + o), the partners' send starts
    auto so_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_soff);
    std::vector<std::array<int, 2>> pc0;
    for (int t = 0; t < 8*nmb; ++t) {
      if (so_h(1,t) >= 0) {pc0.push_back({0, t});}
    }
    const int npc = static_cast<int>(pc0.size());
    vgd_hxb_npc = npc;
    vgd_hxb_pst = Kokkos::View<int**, LayoutWrapper, DevMemSpace>("vet_vgd_hxb_pst",
                                                                    w + 1, np + 1);
    vgd_hxb_pcs = Kokkos::View<int**, LayoutWrapper, DevMemSpace>("vet_vgd_hxb_pcs",
                                                                    w + 1,
                                                                    std::max(npc, 1));
    vgd_hxb_pci = Kokkos::View<int**, LayoutWrapper, DevMemSpace>("vet_vgd_hxb_pci",
                                                                    w + 1,
                                                                    std::max(npc, 1));
    vgd_hxb_pds = Kokkos::View<int**, LayoutWrapper, DevMemSpace>("vet_vgd_hxb_pds",
                                                                    w + 1, np + 1);
    auto pst_h = Kokkos::create_mirror_view(vgd_hxb_pst);
    auto pcs_h = Kokkos::create_mirror_view(vgd_hxb_pcs);
    auto pci_h = Kokkos::create_mirror_view(vgd_hxb_pci);
    auto pds_h = Kokkos::create_mirror_view(vgd_hxb_pds);
    for (int ws = 1; ws <= w; ++ws) {
      for (int p = 0; p < np; ++p) {
        pst_h(ws,p) = vgd_prdsp[ws][p]*nv;
        pds_h(ws,p) = vgd_pdsp[ws][p]*nv;
      }
      pst_h(ws,np) = static_cast<int>(rtot_of(ws));
      pds_h(ws,np) = 0;
      std::vector<std::array<int, 2>> pc;
      for (int t = 0; t < 8*nmb; ++t) {
        if (so_h(ws,t) >= 0) {pc.push_back({so_h(ws,t)*nv, t});}
      }
      std::sort(pc.begin(), pc.end());
      if (static_cast<int>(pc.size()) != npc) {
        ScFatal("vet_gd_halo_exact batch: the send pieces differ between depths");
      }
      for (int q = 0; q < npc; ++q) {pcs_h(ws,q) = pc[q][0]; pci_h(ws,q) = pc[q][1];}
    }
    Kokkos::deep_copy(vgd_hxb_pst, pst_h);
    Kokkos::deep_copy(vgd_hxb_pcs, pcs_h);
    Kokkos::deep_copy(vgd_hxb_pci, pci_h);
    Kokkos::deep_copy(vgd_hxb_pds, pds_h);
    const size_t nall = static_cast<size_t>(nmb)*n*c3*c2;
    if (vgd_hx_mk.extent(0) < nall) {
      Kokkos::realloc(vgd_hx_mk, nall);
      vgd_hx_stamp = 0;
    }
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_gd_halo_exact lists made in groups of " << g
                << " shells (pool " << (G*per/1048576.0) << " MB per rank)" << std::endl;
    }
  }
  const int G = vgd_hxb_g;
  const size_t rcap = vgd_hxb_rcap, scap = vgd_hxb_scap;
  const int npc = vgd_hxb_npc;
  auto mk_ = vgd_hx_mk;
  auto ocv_ = vgd_hxb_oc;
  auto pst_ = vgd_hxb_pst;
  auto pcs_ = vgd_hxb_pcs;
  auto pci_ = vgd_hxb_pci;
  auto pds_ = vgd_hxb_pds;
  auto hb_ = vgd_hxb_hb;
  auto spp_ = vgd_hxb_sp;
  const int nb = np + 4;
  const int njob = static_cast<int>(jobs.size());
  for (int g0 = 0; g0 < njob; g0 += G) {
    const int gn = std::min(G, njob - g0);
    Kokkos::deep_copy(ex_, ocv_, 0);
    // ---- (1) device: marks, sources, flags, scans, bounds, offsets ----
    for (int s = 0; s < gn; ++s) {
      const int pass = jobs[g0+s][0], i = jobs[g0+s][1], l = jobs[g0+s][2];
      const int ws = jobs[g0+s][3];
      const bool inw = (pass == 0);
      const size_t rtot = rtot_of(ws);
      const int mx = std::max(nx2, nx3)*ws*nv;
      DvceArray1D<int> nd_(vgd_hxb_nd.data() + s*rcap, rtot + 1);
      DvceArray1D<int> f0_(vgd_hxb_f0.data() + s*rcap, rtot + 1);
      DvceArray1D<int> f1_(vgd_hxb_f1.data() + s*rcap, rtot + 1);
      DvceArray1D<int> rq_(vgd_hxb_rq.data() + s*rcap, rtot + 1);
      vgd_hx_stamp += 1;
      const int stamp = vgd_hx_stamp;
      const int wsx = ws;
      // marks (ScSweep hx_mark): ghosts of shell i's band read by the stencils of
      // the shells reading this exchange
      auto mark = [&](const int pr, const int ir, const int lr) {
        const bool inwr = (pr == 0);
        auto L_ = vgd_shl[pr];
        const int so = s;
        Kokkos::parallel_for("vet_vgd_hxb_mark", Kokkos::RangePolicy<>(ex_, 0,
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
                 mbsize, lg_, cx1v, cx1f, sth_, pol_, thn_, c2w_, imir_, xr_,
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
                      mk_(((m*n + dq)*c3 + kk)*c2 + jm) = stamp;
                    }
                    if (kk < wb - wsx || kk >= wb + nx3b + wsx) {
                      Kokkos::atomic_add(&ocv_(so), 1);
                    }
                 }
                 if (wf == 2) {continue;}
                }
                const int oo = 3*(ek + 1) + (ej + 1);
                if (hl_(8*m + ((oo < 4) ? oo : (oo - 1))) >= 0) {continue;}
                mk_(((m*n + dr)*c3 + kk)*c2 + jj) = stamp;
                if (jj < wb - wsx || jj >= wb + nx2b + wsx || kk < wb - wsx ||
                    kk >= wb + nx3b + wsx) {
                 Kokkos::atomic_add(&ocv_(so), 1);
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
      (void) is;
      Kokkos::deep_copy(ex_, nd_, 0);
      // sources of the marked ghosts and read-ghost flags (ScHxBuild2 (1))
      par_for("vet_vgd_hxb_need", ex_, 0, nmb1, 0, 7, 0, mx - 1,
      KOKKOS_LAMBDA(const int m, const int o, const int t) {
        if (hl_(8*m + o) >= 0) {return;}
        const int oo = (o < 4) ? o : (o + 1);
        const int dk = oo/3 - 1, dj = oo%3 - 1;
        const int jd = (dj < 0) ? (w - ws) : ((dj == 0) ? w : (w + nx2));
        const int kd = (dk < 0) ? (w - ws) : ((dk == 0) ? w : (w + nx3));
        const int jn = (dj == 0) ? nx2 : ws, kn = (dk == 0) ? nx3 : ws;
        if (t >= nv*kn*jn) {return;}
        const int v = t/(kn*jn);
        const int r1 = t - v*kn*jn;
        const int kk = r1/jn;
        const int jj = r1 - kk*jn;
        const size_t rb0 = static_cast<size_t>(ro_(ws,8*m + o))*nv;
        const size_t q = rb0 + t;
        const bool rd = mk_(((m*nv + v)*c3 + kd + kk)*c2 + jd + jj) == stamp;
        const bool wa = wint && wl_(m,kd+kk,jd+jj) != 0;
        f0_(q) = (rd && !wa) ? 1 : 0;
        f1_(q) = (rd && wa) ? 1 : 0;
        if (!rd) {return;}
        if (wa) {
          for (int u = 0; u < 3; ++u) {
            nd_(rb0 + (m3_(m,kd+kk,jd+jj,v,u)*kn + kk)*jn + jj) = 1;
          }
        } else {
          const int vs = (wl_(m,kd+kk,jd+jj) != 0) ? mp_(m,kd+kk,jd+jj,v) : v;
          nd_(rb0 + (vs*kn + kk)*jn + jj) = 1;
        }
      });
      ScScanOne(ex_, nd_, rtot, vgd_scan_sums, vgd_hxb_tmp);
      ScScanOne(ex_, f0_, rtot, vgd_scan_sums, vgd_hxb_tmp);
      ScScanOne(ex_, f1_, rtot, vgd_scan_sums, vgd_hxb_tmp);
      // needed counts at the partner bounds, read-ghost counts, out count
      {
        const int np_ = np, so = s, nb_ = nb;
        const size_t rt = rtot;
        Kokkos::parallel_for("vet_vgd_hxb_bnd", Kokkos::RangePolicy<>(ex_, 0, np + 4),
        KOKKOS_LAMBDA(const int q) {
          int v = 0;
          if (q <= np_) {
            v = nd_(pst_(ws,q));
          } else if (q == np_ + 1) {
            v = f0_(rt);
          } else if (q == np_ + 2) {
            v = f1_(rt);
          } else {
            v = ocv_(so);
          }
          hb_(so*nb_ + q) = v;
        });
      }
      // wanted offsets, ascending, relative to the partner's piece start
      {
        const int np_ = np;
        Kokkos::parallel_for("vet_vgd_hxb_req", Kokkos::RangePolicy<>(ex_, 0, rtot),
        KOKKOS_LAMBDA(const size_t q) {
          if (nd_(q + 1) == nd_(q)) {return;}
          int lo = 0, hi = np_ - 1;
          while (lo < hi) {
            const int md = (lo + hi + 1)/2;
            if (static_cast<size_t>(pst_(ws,md)) <= q) {lo = md;} else {hi = md - 1;}
          }
          rq_(nd_(q)) = static_cast<int>(q) - pst_(ws,lo);
        });
      }
    }
    ex_.fence();
    // ---- (2) host + MPI: counts, offsets, vote ----
    std::vector<std::vector<int>> rpx(gn, std::vector<int>(np + 1));
    std::vector<std::vector<int>> spx(gn, std::vector<int>(np + 1));
    std::vector<int> nrx(gn), nwx(gn), nout(gn);
    // per partner one message of 2 gn ints: the counts I need from it and where its
    // values start in my receive layout (rpx: the IPC sender writes into my buffer)
    std::vector<int> rc(static_cast<size_t>(2*np)*gn), sc(static_cast<size_t>(2*np)*gn);
    for (int s = 0; s < gn; ++s) {
      for (int p = 0; p <= np; ++p) {rpx[s][p] = hb_(s*nb + p);}
      nrx[s] = hb_(s*nb + np + 1);
      nwx[s] = hb_(s*nb + np + 2);
      nout[s] = hb_(s*nb + np + 3);
      for (int p = 0; p < np; ++p) {
        rc[2*p*gn + s] = rpx[s][p+1] - rpx[s][p];
        rc[2*p*gn + gn + s] = rpx[s][p];
      }
    }
    {
      std::vector<MPI_Request> req;
      req.reserve(2*np);
      for (int p = 0; p < np; ++p) {
        req.emplace_back();
        MPI_Irecv(sc.data() + 2*p*gn, 2*gn, MPI_INT, vgd_prk[p], 7008, vgd_comm,
                  &req.back());
        req.emplace_back();
        MPI_Isend(rc.data() + 2*p*gn, 2*gn, MPI_INT, vgd_prk[p], 7008, vgd_comm,
                  &req.back());
      }
      MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
    }
    std::vector<std::vector<int>> dsto(gn, std::vector<int>(np));
    for (int s = 0; s < gn; ++s) {
      spx[s][0] = 0;
      for (int p = 0; p < np; ++p) {
        spx[s][p+1] = spx[s][p] + sc[2*p*gn + s];
        dsto[s][p] = sc[2*p*gn + gn + s];
      }
      if (static_cast<size_t>(spx[s][np]) + 1 > scap) {
        ScFatal("vet_gd_halo_exact batch: more requested values than the send layout");
      }
    }
    {
      std::vector<MPI_Request> req;
      for (int s = 0; s < gn; ++s) {
        int *so2 = vgd_hxb_so.data() + s*scap;
        int *rq = vgd_hxb_rq.data() + s*rcap;
        for (int p = 0; p < np; ++p) {
          const int nrc = spx[s][p+1] - spx[s][p], nsd = rpx[s][p+1] - rpx[s][p];
          if (nrc > 0) {
            req.emplace_back();
            MPI_Irecv(so2 + spx[s][p], nrc, MPI_INT, vgd_prk[p], 7009, vgd_comm,
                      &req.back());
          }
          if (nsd > 0) {
            req.emplace_back();
            MPI_Isend(rq + rpx[s][p], nsd, MPI_INT, vgd_prk[p], 7009, vgd_comm,
                      &req.back());
          }
        }
      }
      MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
    }
    // vote: no marked ghost outside the region, and the lists fit the budget (counted
    // over the group as if every earlier shell of it is made)
    std::vector<int> ok(gn), gok(gn);
    std::vector<std::array<size_t, 5>> tgt(gn);
    size_t grow = 0;
    if (vgd_hl_cap == 0) {vgd_hl_cap = static_cast<size_t>(vgd_hl_mb) << 20;}
    for (int s = 0; s < gn; ++s) {
      const int e = jobs[g0+s][0]*n1c + jobs[g0+s][1];
      ScHlEntry &E = vgd_hxe[e];
      const size_t nsz[5] = {static_cast<size_t>(std::max(spx[s][np], 1)),
                             static_cast<size_t>(std::max(nrx[s], 1)),
                             static_cast<size_t>(std::max(nrx[s], 1)),
                             static_cast<size_t>(std::max(nwx[s], 1)),
                             static_cast<size_t>(std::max(3*nwx[s], 1))};
      DvceArray1D<int> *arr[5] = {&E.sl, &E.dl, &E.sr, &E.dw, &E.sw};
      size_t have = 0, after = 0;
      for (int q = 0; q < 5; ++q) {
        // growth with 1/4 headroom (sizes change with the angle)
        tgt[s][q] = (arr[q]->extent(0) >= nsz[q]) ? arr[q]->extent(0) :
                    (nsz[q] + nsz[q]/4 + 64);
        have += arr[q]->extent(0);
        after += tgt[s][q];
      }
      ok[s] = (nout[s] == 0) ? 1 : 0;
      if (after > have) {
        if (vgd_hl_bytes + grow + sizeof(int)*(after - have) > vgd_hl_cap) {
          ok[s] = 0;
        } else {
          grow += sizeof(int)*(after - have);
        }
      }
    }
    MPI_Allreduce(ok.data(), gok.data(), gn, MPI_INT, MPI_MIN, vgd_comm);
    // ---- (3) device: the lists ----
    for (int s = 0; s < gn; ++s) {
      const int pass = jobs[g0+s][0], i = jobs[g0+s][1], ws = jobs[g0+s][3];
      const int e = pass*n1c + i;
      ScHlEntry &E = vgd_hxe[e];
      vgd_hx_nout += nout[s];
      E.alpha = vgd_alpha;
      E.scut = vgd_scut;
      E.ws = ws;
      E.exact = false;
      E.gen = -2;
      if (gok[s] == 0) {vgd_hx_nfall += 1.0; continue;}
      DvceArray1D<int> *arr[5] = {&E.sl, &E.dl, &E.sr, &E.dw, &E.sw};
      for (int q = 0; q < 5; ++q) {
        if (arr[q]->extent(0) < tgt[s][q]) {
          vgd_hl_bytes += sizeof(int)*(tgt[s][q] - arr[q]->extent(0));
          Kokkos::realloc(Kokkos::WithoutInitializing, *arr[q], tgt[s][q]);
          vgd_hxb_nrealloc += 1.0;
        }
      }
      const size_t rtot = rtot_of(ws);
      const int mx = std::max(nx2, nx3)*ws*nv;
      DvceArray1D<int> nd_(vgd_hxb_nd.data() + s*rcap, rtot + 1);
      DvceArray1D<int> f0_(vgd_hxb_f0.data() + s*rcap, rtot + 1);
      DvceArray1D<int> f1_(vgd_hxb_f1.data() + s*rcap, rtot + 1);
      DvceArray1D<int> so2_(vgd_hxb_so.data() + s*scap, scap);
      auto sl_ = E.sl;
      auto dl_ = E.dl;
      auto sr_ = E.sr;
      auto dw_ = E.dw;
      auto sw_ = E.sw;
      const int nsx = spx[s][np];
      for (int p = 0; p <= np; ++p) {spp_(s*(np + 1) + p) = spx[s][p];}
      // send addresses (ScHxBuild2 (5))
      {
        const int np_ = np, so = s;
        Kokkos::parallel_for("vet_vgd_hxb_slist", Kokkos::RangePolicy<>(ex_, 0, nsx),
        KOKKOS_LAMBDA(const int c) {
          int lo = 0, hi = np_ - 1;     // partner of send entry c
          while (lo < hi) {
            const int md = (lo + hi + 1)/2;
            if (spp_(so*(np_ + 1) + md) <= c) {lo = md;} else {hi = md - 1;}
          }
          const int q = pds_(ws,lo) + so2_(c);
          int a = 0, b = npc - 1;       // piece holding q
          while (a < b) {
            const int md = (a + b + 1)/2;
            if (pcs_(ws,md) <= q) {a = md;} else {b = md - 1;}
          }
          const int sid = pci_(ws,a);
          const int m = sid/8, o = sid % 8;
          const int oo = (o < 4) ? o : (o + 1);
          const int dk = oo/3 - 1, dj = oo%3 - 1;
          const int jn = (dj == 0) ? nx2 : ws;
          const int kn = (dk == 0) ? nx3 : ws;
          const int t = q - pcs_(ws,a);
          const int v = t/(kn*jn);
          const int r1 = t - v*kn*jn;
          const int kk = r1/jn;
          const int jj = r1 - kk*jn;
          const int js2 = (dj > 0) ? (w + nx2 - ws) : w;
          const int ks2 = (dk > 0) ? (w + nx3 - ws) : w;
          sl_(c) = ((m*nv + v)*c3 + ks2 + kk)*c2 + js2 + jj;
        });
      }
      // read ghosts with their compact source positions
      par_for("vet_vgd_hxb_rlist", ex_, 0, nmb1, 0, 7, 0, mx - 1,
      KOKKOS_LAMBDA(const int m, const int o, const int t) {
        if (hl_(8*m + o) >= 0) {return;}
        const int oo = (o < 4) ? o : (o + 1);
        const int dk = oo/3 - 1, dj = oo%3 - 1;
        const int jd = (dj < 0) ? (w - ws) : ((dj == 0) ? w : (w + nx2));
        const int kd = (dk < 0) ? (w - ws) : ((dk == 0) ? w : (w + nx3));
        const int jn = (dj == 0) ? nx2 : ws, kn = (dk == 0) ? nx3 : ws;
        if (t >= nv*kn*jn) {return;}
        const int v = t/(kn*jn);
        const int r1 = t - v*kn*jn;
        const int kk = r1/jn;
        const int jj = r1 - kk*jn;
        const size_t rb0 = static_cast<size_t>(ro_(ws,8*m + o))*nv;
        const size_t q = rb0 + t;
        const int ad = ((m*nv + v)*c3 + kd + kk)*c2 + jd + jj;
        if (f1_(q + 1) > f1_(q)) {
          const int c = f1_(q);
          dw_(c) = ad;
          for (int u = 0; u < 3; ++u) {
            sw_(3*c + u) = nd_(rb0 + (m3_(m,kd+kk,jd+jj,v,u)*kn + kk)*jn + jj);
          }
        } else if (f0_(q + 1) > f0_(q)) {
          const int vs = (wl_(m,kd+kk,jd+jj) != 0) ? mp_(m,kd+kk,jd+jj,v) : v;
          dl_(f0_(q)) = ad;
          sr_(f0_(q)) = nd_(rb0 + (vs*kn + kk)*jn + jj);
        }
      });
      E.sp = spx[s];
      E.rp = rpx[s];
      E.dsto = dsto[s];
      E.ns = nsx;
      E.nr = nrx[s];
      E.nw = nwx[s];
      E.exact = true;
      E.gen = ++vgd_hl_gen;
      vgd_hx_nmade += 1.0;
      vgd_hx_vexact += rpx[s][np];
    }
    ex_.fence();
    vgd_hxb_ngroup += 1.0;
  }
  vgd_hxb_t += tb.seconds();
#endif
}

}  // namespace radvet
