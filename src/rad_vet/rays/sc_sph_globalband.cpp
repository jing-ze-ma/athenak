//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_sph_globalband.cpp
//! \brief vet_gd on a stretched theta and on the full sphere (use_polar_boundary);
//! the halo of all bands in one exchange (ScGbHaloC, ScGbHaloX).
//! Stretched theta: band rows at the real face midpoints (node table vgd_thn), the
//! fractional row linear in theta between nodes (ScThFj); vgd_sth false if uniform.
//! Full sphere: rows beyond a pole continue the extended chart (-theta north,
//! 2 pi - theta south, phi + pi); the pole slot takes the partner's reversed rows
//! (ScPoleRow); directions are global Cartesian (no map across the pole); polar rings
//! read beyond the band depth set by vet_gd_band_colat are clamped and counted.
//!
//!  called from: rays/sc_bandring.cpp, rays/sc_setup.cpp, rays/sc_sph_halo.cpp,
//!    rays/sc_sph_sweep.cpp
//!  calls into: rays/sc_setup.cpp, rays/sc_sph_walls.cpp

#include <math.h>

#include <algorithm>
#include <iostream>
#include <map>
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
#include "coordinates/coordinates.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rays/sc_segment.hpp"
#include "rad_vet/rays/sc_halo.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScGlobalBandInit
//! \brief vet_gd_band_global: band lists keyed on the global cell, as ScCsInit.
//! Band (kk, jj) of block (lx2, lx3) is J = lx2 nx2 + jj - w, K = lx3 nx3 + kk - w;
//! beyond a pole J -> -1 - J (north) or 2 nt2 - 1 - J (south), K -> K + nx3t/2.
//! Each position at depth d >= 1 takes its owner's interior value; lists are sorted by
//! depth so the depth-ws halo is a prefix. vgd_hloc = -1: all reads from the band.
//! vgd_raw (wedge): J, K periodic per x1 level; vgd_own: on-rank owners.

void RadiationVet::ScGlobalBandInit() {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int g0 = pmy_pack->gids;
  const int me = global_variable::my_rank;
  const int ngb = pm->nmb_total;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int c2w = nx2 + 2*w, c3w = nx3 + 2*w;
  const int nt2 = pm->mesh_indcs.nx2, nt3 = pm->mesh_indcs.nx3;
  const int nb2 = vgd_nb2, nb3 = vgd_nb3;
  if (pm->multilevel) {ScFatal("vet_gd_band_global: SMR/AMR not supported");}
  if (nb2*nx2 != nt2 || nb3*nx3 != nt3 || (vgd_pol && (nt3 % 2) != 0)) {
    ScFatal("vet_gd_band_global: the MeshBlocks must tile the mesh, nx3 even");
  }
  std::vector<int> gof(vgd_nb1*nb2*nb3, -1);
  for (int g = 0; g < ngb; ++g) {
    const LogicalLocation &l = pm->lloc_eachmb[g];
    gof[(l.lx1*nb3 + l.lx3)*nb2 + l.lx2] = g;
  }
  // owner (gid, band k, band j) of band position (kk, jj) of block g; false: none
  auto owner = [&](const int g, const int kk, const int jj, int &go, int &ko,
                   int &jo) -> bool {
    const LogicalLocation &l = pm->lloc_eachmb[g];
    int J = l.lx2*nx2 + jj - w, K = l.lx3*nx3 + kk - w;
    if (!vgd_pol) {
      J = ((J % nt2) + nt2) % nt2;
    } else if (J < 0) {
      J = -1 - J;
      K += nt3/2;
    } else if (J >= nt2) {
      J = 2*nt2 - 1 - J;
      K += nt3/2;
    }
    if (J < 0 || J >= nt2) {return false;}
    K = ((K % nt3) + nt3) % nt3;
    const int b2 = J/nx2, b3 = K/nx3;
    go = gof[(l.lx1*nb3 + b3)*nb2 + b2];
    ko = K - b3*nx3 + w;
    jo = J - b2*nx2 + w;
    return true;
  };
  auto depth = [&](const int kk, const int jj) -> int {
    const int dk = std::max(std::max(w - kk, kk - (w + nx3 - 1)), 0);
    const int dj = std::max(std::max(w - jj, jj - (w + nx2 - 1)), 0);
    return std::max(dk, dj);
  };
  Kokkos::realloc(vgd_own, nmb, c3w, c2w, 3);
  auto ow_h = Kokkos::create_mirror_view(vgd_own);
  for (int m = 0; m < nmb; ++m) {
    for (int kk = 0; kk < c3w; ++kk) {
      for (int jj = 0; jj < c2w; ++jj) {
        int go = -1, ko = kk, jo = jj;
        if (!owner(g0 + m, kk, jj, go, ko, jo) || pm->rank_eachmb[go] != me) {
          go = g0 - 1;
          ko = kk;
          jo = jj;
        }
        ow_h(m,kk,jj,0) = go - g0;
        ow_h(m,kk,jj,1) = ko;
        ow_h(m,kk,jj,2) = jo;
      }
    }
  }
  Kokkos::deep_copy(vgd_own, ow_h);
  // buckets by depth (1..w): on-rank (6 ints), per partner rank send / recv (3 ints)
  std::vector<std::vector<int>> lb(w + 1);
  std::map<int, std::vector<std::vector<int>>> sb, rb;
  // source cell's global (J, K) per send / receive entry (compact keep test, both sides)
  std::map<int, std::vector<std::vector<int>>> sbj, rbj;
  for (int gr = 0; gr < ngb; ++gr) {
    const int rr = pm->rank_eachmb[gr];
    for (int kk = 0; kk < c3w; ++kk) {
      for (int jj = 0; jj < c2w; ++jj) {
        const int d = depth(kk, jj);
        if (d == 0) {continue;}
        int go, ko, jo;
        if (!owner(gr, kk, jj, go, ko, jo)) {
          ScFatal("vet_gd_band_global: a band position beyond the sphere (band deeper "
                  "than the theta rows)");
        }
        const int ro = pm->rank_eachmb[go];
        if (rr != me && ro != me) {continue;}
        const int jg = pm->lloc_eachmb[go].lx2*nx2 + jo - w;
        const int kg = pm->lloc_eachmb[go].lx3*nx3 + ko - w;
        // tw 1: image across a theta wall (keeps its whole branch); tw 2 (vgd_raw):
        // image across a phi wall only (keeps what its wall map reads, ScGbHaloC)
        const int jr = pm->lloc_eachmb[gr].lx2*nx2 + jj - w;
        const int kr = pm->lloc_eachmb[gr].lx3*nx3 + kk - w;
        const bool pw = (kr < 0 || kr >= nt3);
        const int tw = (vgd_pol || !(pw || jr < 0 || jr >= nt2)) ? 0
                       : ((pw && jr >= 0 && jr < nt2 && vgd_raw && !vgd_thmir &&
                           vgd_thvac == 0) ? 2 : 1);
        if (rr == me && ro == me) {
          lb[d].insert(lb[d].end(), {gr - g0, kk, jj, go - g0, ko, jo, d});
        } else if (rr == me) {
          auto &v = rb[ro];
          if (v.empty()) {v.resize(w + 1);}
          v[d].insert(v[d].end(), {gr - g0, kk, jj});
          auto &vj = rbj[ro];
          if (vj.empty()) {vj.resize(w + 1);}
          vj[d].insert(vj[d].end(), {jg, kg, tw});
        } else {
          auto &v = sb[rr];
          if (v.empty()) {v.resize(w + 1);}
          v[d].insert(v[d].end(), {go - pm->gids_eachrank[me], ko, jo});
          auto &vj = sbj[rr];
          if (vj.empty()) {vj.resize(w + 1);}
          vj[d].insert(vj[d].end(), {jg, kg, tw});
        }
      }
    }
  }
  // on-rank list in depth order, its prefix counts
  std::vector<int> loc;
  vgd_gdnl.assign(w + 1, 0);
  for (int d = 1; d <= w; ++d) {
    for (size_t q = 0; q < lb[d].size(); q += 7) {
      loc.insert(loc.end(), lb[d].begin() + q, lb[d].begin() + q + 6);
    }
    vgd_gdnl[d] = static_cast<int>(loc.size()/6);
  }
  vgd_csnloc = static_cast<int>(loc.size()/6);
  Kokkos::realloc(vgd_csloc, std::max(vgd_csnloc, 1), 6);
  {
    auto l_h = Kokkos::create_mirror_view(vgd_csloc);
    for (int e = 0; e < vgd_csnloc; ++e) {
      for (int c = 0; c < 6; ++c) {l_h(e,c) = loc[6*e + c];}
    }
    Kokkos::deep_copy(vgd_csloc, l_h);
  }
  // partners (sorted), their segments in depth order
  vgd_csprk.clear();
  for (const auto &x : sb) {vgd_csprk.push_back(x.first);}
  for (const auto &x : rb) {
    if (sb.find(x.first) == sb.end()) {vgd_csprk.push_back(x.first);}
  }
  std::sort(vgd_csprk.begin(), vgd_csprk.end());
  const int np = static_cast<int>(vgd_csprk.size());
  vgd_csps.assign(np, 0);
  vgd_cspn.assign(np, 0);
  vgd_cspr.assign(np, 0);
  vgd_cspm.assign(np, 0);
  vgd_gdsn.assign(w + 1, std::vector<int>(np, 0));
  vgd_gdrn.assign(w + 1, std::vector<int>(np, 0));
  std::vector<int> sl, rl, sd, rd, sj, rj;
  for (int p = 0; p < np; ++p) {
    const int rk = vgd_csprk[p];
    vgd_csps[p] = static_cast<int>(sl.size()/3);
    vgd_cspr[p] = static_cast<int>(rl.size()/3);
    for (int d = 1; d <= w; ++d) {
      if (sb.count(rk)) {
        const auto &v = sb[rk][d];
        sl.insert(sl.end(), v.begin(), v.end());
        sd.insert(sd.end(), v.size()/3, d);
        const auto &vj = sbj[rk][d];
        sj.insert(sj.end(), vj.begin(), vj.end());
      }
      if (rb.count(rk)) {
        const auto &v = rb[rk][d];
        rl.insert(rl.end(), v.begin(), v.end());
        rd.insert(rd.end(), v.size()/3, d);
        const auto &vj = rbj[rk][d];
        rj.insert(rj.end(), vj.begin(), vj.end());
      }
      vgd_gdsn[d][p] = static_cast<int>(sl.size()/3) - vgd_csps[p];
      vgd_gdrn[d][p] = static_cast<int>(rl.size()/3) - vgd_cspr[p];
    }
    vgd_cspn[p] = static_cast<int>(sl.size()/3) - vgd_csps[p];
    vgd_cspm[p] = static_cast<int>(rl.size()/3) - vgd_cspr[p];
  }
  vgd_csnsnd = static_cast<int>(sl.size()/3);
  vgd_csnrcv = static_cast<int>(rl.size()/3);
  Kokkos::realloc(vgd_cssl, std::max(vgd_csnsnd, 1), 3);
  Kokkos::realloc(vgd_csrl, std::max(vgd_csnrcv, 1), 3);
  Kokkos::realloc(vgd_gdps, std::max(vgd_csnsnd, 1));
  Kokkos::realloc(vgd_gdpr, std::max(vgd_csnrcv, 1));
  {
    auto s_h = Kokkos::create_mirror_view(vgd_cssl);
    auto r_h = Kokkos::create_mirror_view(vgd_csrl);
    auto ds_h = Kokkos::create_mirror_view(vgd_gdps);
    auto dr_h = Kokkos::create_mirror_view(vgd_gdpr);
    for (int e = 0; e < vgd_csnsnd; ++e) {
      for (int c = 0; c < 3; ++c) {s_h(e,c) = sl[3*e + c];}
      ds_h(e) = sd[e];
    }
    for (int e = 0; e < vgd_csnrcv; ++e) {
      for (int c = 0; c < 3; ++c) {r_h(e,c) = rl[3*e + c];}
      dr_h(e) = rd[e];
    }
    Kokkos::deep_copy(vgd_cssl, s_h);
    Kokkos::deep_copy(vgd_csrl, r_h);
    Kokkos::deep_copy(vgd_gdps, ds_h);
    Kokkos::deep_copy(vgd_gdpr, dr_h);
  }
  // message buffers: one shell of every direction, or ln chi and ln S of some shells
  const int nent = std::max(std::max(vgd_csnsnd, vgd_csnrcv), 1);
  const int nvb = std::max(vgd_n, 2*8);
  // device buffers (CUDA-aware MPI); the cs pinned host buffers get one element
  Kokkos::realloc(vgd_sbuf, static_cast<size_t>(nent)*nvb);
  Kokkos::realloc(vgd_rbuf, static_cast<size_t>(nent)*nvb);
  Kokkos::realloc(vgd_cssb, 1);
  Kokkos::realloc(vgd_csrb, 1);
  vgd_hmpi = (np > 0);
  // band column phi from the global index, P0 + (K + 0.5) D (layout independent; the
  // branch test n . r_hat >= 0 of a tangent cell depends on the last bit)
  {
    const double p0 = pm->mesh_size.x3min, dp = (pm->mesh_size.x3max - p0)/nt3;
    vgd_phc_h.assign(static_cast<size_t>(nmb)*c3w, 0.0);
    Kokkos::realloc(vgd_phc, nmb, c3w);
    auto ph_h = Kokkos::create_mirror_view(vgd_phc);
    for (int m = 0; m < nmb; ++m) {
      const int lx3 = pm->lloc_eachmb[g0 + m].lx3;
      for (int kb = 0; kb < c3w; ++kb) {
        const int K = lx3*nx3 + kb - w;
        const double v = p0 + (K + 0.5)*dp;
        vgd_phc_h[static_cast<size_t>(m)*c3w + kb] = v;
        ph_h(m,kb) = v;
      }
    }
    Kokkos::deep_copy(vgd_phc, ph_h);
  }
  // compact global band (vet_gd_halo_compact > 0, MPI band); kept lists per direction set
  vgd_gbc = (vgd_hcomp > 0) && vgd_hmpi;
  if (vgd_gbc) {
    const int n = vgd_n;
    Kokkos::realloc(vgd_gbjs, std::max(vgd_csnsnd, 1), 3);
    Kokkos::realloc(vgd_gbjr, std::max(vgd_csnrcv, 1), 3);
    {
      auto a_h = Kokkos::create_mirror_view(vgd_gbjs);
      for (int e = 0; e < 3*vgd_csnsnd; ++e) {a_h(e/3,e%3) = sj[e];}
      Kokkos::deep_copy(vgd_gbjs, a_h);
      auto b_h = Kokkos::create_mirror_view(vgd_gbjr);
      for (int e = 0; e < 3*vgd_csnrcv; ++e) {b_h(e/3,e%3) = rj[e];}
      Kokkos::deep_copy(vgd_gbjr, b_h);
    }
    // global rows: node theta, cell width, polar zone (colatitude < vet_gd_band_colat,
    // + w + 1 rows: every row a clamped polar read can reach)
    Kokkos::realloc(vgd_gbth, nt2);
    Kokkos::realloc(vgd_gbdt, nt2);
    Kokkos::realloc(vgd_gbcap, nt2);
    auto th_h = Kokkos::create_mirror_view(vgd_gbth);
    auto dt_h = Kokkos::create_mirror_view(vgd_gbdt);
    auto cp_h = Kokkos::create_mirror_view(vgd_gbcap);
    const double tb = vgd_bcolat*M_PI/180.0;
    int nz = 0;
    for (int J = 0; J < nt2; ++J) {
      const double f0 = ScThetaFace(pm, J), f1 = ScThetaFace(pm, J + 1);
      th_h(J) = 0.5*(f0 + f1);
      dt_h(J) = f1 - f0;
      if (0.5*(f0 + f1) < tb) {nz = J + 1;}
    }
    for (int J = 0; J < nt2; ++J) {
      cp_h(J) = (vgd_pol && (J < nz + w + 1 || J >= nt2 - nz - w - 1)) ? 1 : 0;
    }
    Kokkos::deep_copy(vgd_gbth, th_h);
    Kokkos::deep_copy(vgd_gbdt, dt_h);
    Kokkos::deep_copy(vgd_gbcap, cp_h);
    vgd_gbdmin = static_cast<Real>(ScThetaDmin());
    auto x1v_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmy_pack->pcoord->x1v);
    auto x1f_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmy_pack->pcoord->xx1f);
    vgd_r1v.assign(x1v_h.extent(1), 0.0);
    vgd_r1f.assign(x1f_h.extent(1), 0.0);
    for (size_t q = 0; q < vgd_r1v.size(); ++q) {vgd_r1v[q] = x1v_h(0,q);}
    for (size_t q = 0; q < vgd_r1f.size(); ++q) {vgd_r1f[q] = x1f_h(0,q);}
    Kokkos::realloc(vgd_gbfs, static_cast<size_t>(vgd_csnsnd)*n + 1);
    Kokkos::realloc(vgd_gbfr, static_cast<size_t>(vgd_csnrcv)*n + 1);
    vgd_gbpb = Kokkos::View<int*, Kokkos::SharedHostPinnedSpace>("vet_vgd_gb_pb",
                                                                 8*np + 2);
    if (vgd_scan_sums.extent(0) == 0) {Kokkos::realloc(vgd_scan_sums, 1);}
    if (me == 0) {
      std::cout << "<rad_vet> vet_gd_band_global: COMPACT exchange (vet_gd_halo_compact "
                   "= "
                << vgd_hcomp << "): per shell and pass only the kept (cell, "
                << "direction) values; polar zone rows " << nz + w + 1 << " per pole "
                << "keep their whole branch" << std::endl;
    }
  }
  // vet_gd_halo_ipc (vgd_raw): kept-list exchanges GPU-driven over IPC, one receive
  // region per partner (largest dense receive segment); collective
  if (vgd_raw && vgd_ipc_on && vgd_hcomp > 0 && vgd_ipc_state == 0) {
    // clean_rot_fuse: room for both cross-fade sets (ScGbHaloX)
    size_t cs = 1;
    for (int p = 0; p < np; ++p) {cs = std::max<size_t>(cs, 1ULL*vgd_cspm[p]*vgd_n);}
    if (cln_rfuse) {cs *= 2;}
    ScIpcInit(vgd_csprk, cs, std::max(np, 1)*cs);
  }
  vgd_gbcf.clear();
  vgd_gbcb.clear();
  vgd_gbcf_xb.clear();
  vgd_gbcb_xb.clear();
  // every slot reads the band
  Kokkos::deep_copy(vgd_hloc, -1);
  for (int t = 0; t < 8*nmb; ++t) {vgd_hrank[t] = me;}
  if (me == 0) {
    std::cout << "<rad_vet> vet_gd_band_global: band " << w
              << " cells keyed on the global cell (exact copies from any block; rank 0: "
              << "on-rank " << vgd_csnloc << ", sent " << vgd_csnsnd << ", received "
              << vgd_csnrcv << " positions per shell at full depth, partners " << np
              << "); stencil clamped to the cell's own window (layout independent)"
              << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::ScGbHaloC
//! \brief compact per-shell exchange of the global band (vet_gd_halo_compact > 0):
//! only (entry, direction) values passing ScHaloKeepTP are sent. Both sides compute the
//! same flags from the source cell's (J, K), depth and shell radii, so messages carry no
//! index. Polar-zone rows keep their whole branch; vgd_raw phi-wall image flags come
//! from the receiver. Unsent positions get 0 (never read). Returns false.

template <class V>
bool RadiationVet::ScGbHaloC(V &a, const int i, const int ws) {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  auto a_ = a;
  const int nv = vgd_n;
  // lean (sweep on the default stream): pack / unpack stream-ordered, host waits for pack
  const bool lean = vgd_cur.impl_instance_id() == DevExeSpace().impl_instance_id();
  if (!lean) vgd_cur.fence();   // the shell's intensities are complete
  const int wq = std::min(ws, static_cast<int>(vgd_gdnl.size()) - 1);
  {
    const int nl = vgd_raw ? 0 : vgd_gdnl[wq];   // vgd_raw: on-rank owners read direct
    if (nl > 0) {
      auto L_ = vgd_csloc;
      par_for("vet_vgdgb_loc", DevExeSpace(), 0, nl - 1, 0, nv - 1,
      KOKKOS_LAMBDA(const int e, const int v) {
        a_(L_(e,0),v,L_(e,1),L_(e,2),i) = a_(L_(e,3),v,L_(e,4),L_(e,5),i);
      });
    }
  }
#if MPI_PARALLEL_ENABLED
  if (!vgd_hmpi) {return false;}
  const int np = static_cast<int>(vgd_csprk.size());
  const int ns = vgd_csnsnd, nr = vgd_csnrcv;
  const bool inw = vgd_hinw;
  const int mode = vgd_hcomp;
  const Real p0 = pm->mesh_size.x3min;
  const Real dph = (pm->mesh_size.x3max - p0)/pm->mesh_indcs.nx3;
  const Real dmin = vgd_gbdmin;
  const int lcut = indcs.is + vgd_scut;
  const Real rs = vgd_r1v[i];
  const bool hasm = (i >= lcut);
  const Real rm = (i > lcut) ? vgd_r1v[i-1] : vgd_r1f[i];
  const bool hasp = (i < indcs.ie);
  const Real rp = hasp ? vgd_r1v[i+1] : rs;
  auto dir_ = vgd_dir;
  auto th_ = vgd_gbth;
  auto dt_ = vgd_gbdt;
  auto cp_ = vgd_gbcap;
  auto JS_ = vgd_gbjs;
  auto JR_ = vgd_gbjr;
  auto DS_ = vgd_gdps;
  auto DR_ = vgd_gdpr;
  auto S_ = vgd_cssl;
  auto R_ = vgd_csrl;
  auto fs_ = vgd_gbfs;
  auto fr_ = vgd_gbfr;
  const int nmx = std::max(ns, nr);
  // vgd_raw: kept lists and bounds per (shell, pass, depth) of the current direction
  // set, made at the first exchange; pack / unpack run over kept values only
  const int key = (2*i + (inw ? 1 : 0))*(vgd_w + 1) + wq;
  auto &cf = vgd_gbcf[vgd_raw ? key : -1];
  auto &cb = vgd_gbcb[vgd_raw ? key : -1];
  // the lists depend on the cut vgd_scut (last entry); a cut change remakes the lists
  // of the shells it moves across
  auto side = [&](const int c) {
    return (i > indcs.is + c) ? 1 : ((i < indcs.is + c) ? -1 : 0);
  };
  if (vgd_raw && !cb.empty() && side(cb.back()) != side(vgd_scut)) {cb.clear();}
  auto pv = Kokkos::subview(vgd_gbpb, std::make_pair(4*np, 8*np));
  if (vgd_raw && !cb.empty()) {
    for (int q = 0; q < 4*np; ++q) {pv(q) = cb[q];}
  } else {
    // (1) flags of both sides (send: my interior cells; receive: my band positions)
    par_for("vet_vgdgb_flag", DevExeSpace(), 0, nmx - 1, 0, nv - 1,
    KOKKOS_LAMBDA(const int e, const int v) {
      const Real nx = dir_(v,0), ny = dir_(v,1), nz = dir_(v,2);
      for (int sd = 0; sd < 2; ++sd) {
        if (e >= ((sd == 0) ? ns : nr)) {continue;}
        const int dep = (sd == 0) ? DS_(e) : DR_(e);
        const int tw = (sd == 0) ? JS_(e,2) : JR_(e,2);
        int f = 0;
        if (dep <= wq && tw != 2) {
          const int J = (sd == 0) ? JS_(e,0) : JR_(e,0);
          const int K = (sd == 0) ? JS_(e,1) : JR_(e,1);
          const Real gam = 2.0*fmax(dt_(J), dph);
          const bool all = (cp_(J) != 0) || (tw != 0);
          f = ScHaloKeepTP(nx, ny, nz, th_(J), p0 + (K + 0.5)*dph, inw, mode, all, dep,
                           rs, rm, rp, hasm, hasp, gam, dmin) ? 1 : 0;
        }
        if (sd == 0) {fs_(static_cast<size_t>(e)*nv + v) = f;}
        else {fr_(static_cast<size_t>(e)*nv + v) = f;}
      }
    });
    // (1b) vgd_raw phi-wall image (tw 2): flag the directions its wall map reads for
    // each q kept at the image's unwrapped position; senders take the receivers' flags
    if (vgd_raw) {
      auto phc_ = vgd_phc;
      auto mp_ = vgd_map;
      auto m3_ = vgd_m3;
      const bool wint = vgd_wint;
      par_for("vet_vgdgb_wflag", DevExeSpace(), 0, std::max(nr, 1) - 1, 0, nv - 1,
      KOKKOS_LAMBDA(const int e, const int q) {
        if (e >= nr || DR_(e) > wq || JR_(e,2) != 2) {return;}
        const int m = R_(e,0), kk = R_(e,1), jj = R_(e,2), J = JR_(e,0);
        const Real gam = 2.0*fmax(dt_(J), dph);
        if (!ScHaloKeepTP(dir_(q,0), dir_(q,1), dir_(q,2), th_(J), phc_(m,kk), inw, mode,
                          cp_(J) != 0, DR_(e), rs, rm, rp, hasm, hasp, gam, dmin)) {
          return;
        }
        for (int t = 0; t < (wint ? 3 : 1); ++t) {
          const int v = wint ? m3_(m,kk,jj,q,t) : mp_(m,kk,jj,q);
          Kokkos::atomic_max(&fr_(static_cast<size_t>(e)*nv + v), 1);
        }
      });
      Kokkos::fence();
      std::vector<MPI_Request> fq;
      for (int p = 0; p < np; ++p) {
        const int rk = vgd_csprk[p];
        const int sn = vgd_gdsn[wq][p]*nv, rn = vgd_gdrn[wq][p]*nv;
        if (sn > 0) {
          fq.emplace_back();
          MPI_Irecv(fs_.data() + static_cast<size_t>(vgd_csps[p])*nv, sn, MPI_INT, rk,
                    7113, MPI_COMM_WORLD, &fq.back());
        }
        if (rn > 0) {
          fq.emplace_back();
          MPI_Isend(fr_.data() + static_cast<size_t>(vgd_cspr[p])*nv, rn, MPI_INT, rk,
                    7113, MPI_COMM_WORLD, &fq.back());
        }
      }
      MPI_Waitall(static_cast<int>(fq.size()), fq.data(), MPI_STATUSES_IGNORE);
    }
    // (2) exclusive scans in place (flag -> position; totals at the end)
    ScScan2(DevExeSpace(), vgd_gbfs, static_cast<size_t>(ns)*nv, vgd_gbfr,
            static_cast<size_t>(nr)*nv, "vet_vgdgb_scan2", vgd_scan_sums);
    // (3) positions at the partner segment bounds (depth-wq prefix) -> pinned host
    {
      std::vector<int> bd(4*np);
      for (int p = 0; p < np; ++p) {
        bd[4*p] = vgd_csps[p];
        bd[4*p + 1] = vgd_csps[p] + vgd_gdsn[wq][p];
        bd[4*p + 2] = vgd_cspr[p];
        bd[4*p + 3] = vgd_cspr[p] + vgd_gdrn[wq][p];
      }
      for (int q = 0; q < 4*np; ++q) {vgd_gbpb(q) = bd[q];}
    }
    auto pb_ = vgd_gbpb;
    const int np_ = np;
    Kokkos::fence();
    Kokkos::parallel_for("vet_vgdgb_bnd", Kokkos::RangePolicy<>(DevExeSpace(), 0, 4*np),
    KOKKOS_LAMBDA(const int q) {
      const size_t at = static_cast<size_t>(pb_(q))*nv;
      pb_(4*np_ + q) = ((q & 3) < 2) ? fs_(at) : fr_(at);
    });
    if (vgd_raw) {
      Kokkos::fence();
      vgd_gbrec.erase(key);
      cb.assign(pv.data(), pv.data() + 4*np);
      for (int sd = 0; sd < 2; ++sd) {
        auto f_ = (sd == 0) ? fs_ : fr_;
        const int nq = ((sd == 0) ? ns : nr)*nv;
        int nk = 0;
        Kokkos::deep_copy(nk, Kokkos::subview(f_, nq));
        cf[sd] = DvceArray1D<int>("vet_vgdgb_keep", std::max(nk, 1));
        auto l_ = cf[sd];
        Kokkos::parallel_for("vet_vgdgb_list", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                                                                     nq),
        KOKKOS_LAMBDA(const int q) {if (f_(q + 1) > f_(q)) {l_(f_(q)) = q;}});
        cb.push_back(nk);
      }
      cb.push_back(vgd_scut);
    }
  }
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  // vet_gd_halo_ipc: pack into partner p's receive region of me, set its flag, wait for
  // mine, unpack: stream-ordered, no host wait, no MPI; bitwise the MPI result.
  // Parity buffers: a partner's flag of exchange s follows its unpack of s - 1.
  if (vgd_raw && !cb.empty() && vgd_ipc_state == 1 && lean) {
    ++vgd_ipc_seq;
    const long long seq = vgd_ipc_seq, par = seq & 1;
    const int np_ = np;
    auto pslot = vgd_ipc_pslot.d_view;
    ScIpcIdx b0{}, b1{};
    // bodies (> 512 B) as kernel arguments, not constant memory (avoids a host wait)
    auto lw = [](const int n) {return Kokkos::Experimental::require(Kokkos::RangePolicy<>(
        DevExeSpace(), 0, n), Kokkos::Experimental::WorkItemProperty::HintLightWeight);};
    for (int p = 0; p < np; ++p) {b0.v[p] = cb[4*p]; b1.v[p] = cb[4*p + 1];}
    auto prb = vgd_ipc_prb.d_view;
    auto phalf = vgd_ipc_phalf.d_view;
    auto pcs = vgd_ipc_pcs.d_view;
    auto psig = vgd_ipc_psig.d_view;
    auto ls_ = cf[0];
    Kokkos::parallel_for("vet_vgdgb_ipc_pack", lw(std::max(cb[4*np], 1)),
    KOKKOS_LAMBDA(const int o) {
      int p = 0;
      while (p < np_ && !(o >= b0.v[p] && o < b1.v[p])) {++p;}
      if (p == np_) {return;}
      const int e = ls_(o)/nv;
      reinterpret_cast<Real *>(prb(p))[par*phalf(p) + pslot(p)*pcs(p) + o - b0.v[p]] =
          a_(S_(e,0),ls_(o) - e*nv,S_(e,1),S_(e,2),i);
    });
    Kokkos::parallel_for("vet_vgdgb_ipc_set", Kokkos::RangePolicy<>(DevExeSpace(), 0, 1),
    KOKKOS_LAMBDA(const int) {ScIpcSetN(psig.data(), pslot.data(), np_, seq);});
    for (int p = 0; p < np; ++p) {
      vgd_nbyte += static_cast<Real>(b1.v[p] - b0.v[p])*sizeof(Real);
      vgd_csvsent += static_cast<double>(b1.v[p] - b0.v[p]);
      vgd_csvdense += static_cast<double>(vgd_gdsn[wq][p])*nv;
    }
    vgd_ipc_n += 1.0;
    for (int p = 0; p < np; ++p) {b0.v[p] = cb[4*p + 2]; b1.v[p] = cb[4*p + 3];}
    long long *sig = vgd_ipc_sig;
    int *err = vgd_ipc_err;
    Kokkos::parallel_for("vet_vgdgb_ipc_wait", Kokkos::RangePolicy<>(DevExeSpace(), 0, 1),
    KOKKOS_LAMBDA(const int) {ScIpcWaitN(sig, np_, seq, err);});
    vgd_nexch += 1.0;
    const Real *rb = vgd_ipc_rbuf + par*vgd_ipc_half;
    const size_t cs = vgd_ipc_cs;
    auto lr_ = cf[1];
    Kokkos::parallel_for("vet_vgdgb_ipc_unpack", lw(std::max(cb[4*np + 1], 1)),
    KOKKOS_LAMBDA(const int o) {
      int p = 0;
      while (p < np_ && !(o >= b0.v[p] && o < b1.v[p])) {++p;}
      if (p == np_) {return;}
      const int e = lr_(o)/nv;
      a_(R_(e,0),lr_(o) - e*nv,R_(e,1),R_(e,2),i) = rb[p*cs + o - b0.v[p]];
    });
    return false;
  }
#endif
  // (4) compact pack (my kept values in list order)
  auto sb_ = vgd_sbuf;
  auto rb_ = vgd_rbuf;
  auto ls_ = cf[0], lr_ = cf[1];
  if (vgd_raw && cb[4*np] > 0) {
    Kokkos::parallel_for("vet_vgdgb_kpack",
                         Kokkos::RangePolicy<>(DevExeSpace(), 0, cb[4*np]),
    KOKKOS_LAMBDA(const int o) {
      const int e = ls_(o)/nv;
      sb_(o) = a_(S_(e,0),ls_(o) - e*nv,S_(e,1),S_(e,2),i);
    });
  } else if (ns > 0 && !vgd_raw) {
    par_for("vet_vgdgb_pack", DevExeSpace(), 0, ns - 1, 0, nv - 1,
    KOKKOS_LAMBDA(const int e, const int v) {
      const size_t q = static_cast<size_t>(e)*nv + v;
      if (fs_(q + 1) == fs_(q)) {return;}
      sb_(fs_(q)) = a_(S_(e,0),v,S_(e,1),S_(e,2),i);
    });
  }
  // lean: wait for the pack alone (an event on the sweep's stream, not a fence)
  if (lean) {
    ScPackEvtRecord(DevExeSpace());
    ScPackEvtWait(DevExeSpace());
  } else {
    Kokkos::fence();
  }
  Kokkos::Timer tq;
  std::vector<MPI_Request> req;
  double nsent = 0.0;
  for (int p = 0; p < np; ++p) {
    const int rk = vgd_csprk[p];
    const int s0 = pv(4*p), s1 = pv(4*p + 1), r0 = pv(4*p + 2), r1 = pv(4*p + 3);
    if (r1 > r0) {
      req.emplace_back();
      MPI_Irecv(rb_.data() + r0, r1 - r0, MPI_ATHENA_REAL, rk, 7111, MPI_COMM_WORLD,
                &req.back());
    }
    if (s1 > s0) {
      req.emplace_back();
      MPI_Isend(sb_.data() + s0, s1 - s0, MPI_ATHENA_REAL, rk, 7111, MPI_COMM_WORLD,
                &req.back());
      vgd_nbyte += static_cast<Real>(s1 - s0)*sizeof(Real);
      nsent += static_cast<double>(s1 - s0);
    }
    vgd_csvdense += static_cast<double>(vgd_gdsn[wq][p])*nv;
  }
  MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
  vgd_tmpi += tq.seconds();
  vgd_nexch += 1.0;
  vgd_csvsent += nsent;
  // (5) expand into the band (unsent positions of the prefix: 0)
  if (vgd_raw && cb[4*np + 1] > 0) {
    Kokkos::parallel_for("vet_vgdgb_kunpack",
                         Kokkos::RangePolicy<>(DevExeSpace(), 0, cb[4*np + 1]),
    KOKKOS_LAMBDA(const int o) {
      const int e = lr_(o)/nv;
      a_(R_(e,0),lr_(o) - e*nv,R_(e,1),R_(e,2),i) = rb_(o);
    });
  } else if (nr > 0 && !vgd_raw) {
    par_for("vet_vgdgb_unpack", DevExeSpace(), 0, nr - 1, 0, nv - 1,
    KOKKOS_LAMBDA(const int e, const int v) {
      if (DR_(e) > wq) {return;}
      const size_t q = static_cast<size_t>(e)*nv + v;
      a_(R_(e,0),v,R_(e,1),R_(e,2),i) = (fr_(q + 1) > fr_(q)) ? rb_(fr_(q)) : 0.0;
    });
  }
  if (!lean) Kokkos::fence();
#else
  (void) i;
#endif
  return false;
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::ScGbHaloX
//! \brief clean_rot_fuse (vgd_raw): compact exchange of shell i for both cross-fade
//! sets in one message per partner (current set's piece, then the other's); bitwise two
//! ScGbHaloC calls. false (nothing done) unless both sets hold the kept lists.

bool RadiationVet::ScGbHaloX(const int i, const int ws) {
#if MPI_PARALLEL_ENABLED
  const int np = static_cast<int>(vgd_csprk.size());
  if (!vgd_raw || !vgd_gbc || !vgd_hmpi || np < 1 || np > VGD_IPC_NS) {
    return false;
  }
  const int wq = std::min(ws, static_cast<int>(vgd_gdnl.size()) - 1);
  const int key = (2*i + (vgd_hinw ? 1 : 0))*(vgd_w + 1) + wq;
  auto ia = vgd_gbcb.find(key), ib = vgd_gbcb_xb.find(key);
  if (ia == vgd_gbcb.end() || ib == vgd_gbcb_xb.end()) {return false;}
  const std::vector<int> &A = ia->second, &B = ib->second;
  // each set's pieces contiguous in partner order: piece p at [c(4p), c(4p+1))
  auto ok = [&](const std::vector<int> &c) {
    if (static_cast<int>(c.size()) != 4*np + 3 || c.back() != vgd_scut || c[0] != 0 ||
        c[2] != 0) {
      return false;
    }
    for (int p = 0; p + 1 < np; ++p) {
      if (c[4*p + 1] != c[4*p + 4] || c[4*p + 3] != c[4*p + 6]) {return false;}
    }
    return c[4*np - 3] == c[4*np] && c[4*np - 1] == c[4*np + 1];
  };
  if (!ok(A) || !ok(B)) {return false;}
  ScIpcIdx sa{}, sb{}, sm{}, ra{}, rb{}, rm{};
  for (int p = 0; p <= np; ++p) {
    sa.v[p] = (p < np) ? A[4*p] : A[4*np];
    sb.v[p] = (p < np) ? B[4*p] : B[4*np];
    ra.v[p] = (p < np) ? A[4*p + 2] : A[4*np + 1];
    rb.v[p] = (p < np) ? B[4*p + 2] : B[4*np + 1];
    sm.v[p] = sa.v[p] + sb.v[p];
    rm.v[p] = ra.v[p] + rb.v[p];
  }
  const int nsA = sa.v[np], nS = sm.v[np], nrA = ra.v[np], nR = rm.v[np];
  if (static_cast<int>(vgd_xff_sb.extent(0)) < std::max(nS, 1)) {
    Kokkos::realloc(vgd_xff_sb, std::max(nS, 1) + nS/8);
  }
  if (static_cast<int>(vgd_xff_rb.extent(0)) < std::max(nR, 1)) {
    Kokkos::realloc(vgd_xff_rb, std::max(nR, 1) + nR/8);
  }
  auto sbuf_ = vgd_xff_sb;
  auto rbuf_ = vgd_xff_rb;
  auto aA_ = vgd_i;
  auto aB_ = vgd_i_xb;
  auto S_ = vgd_cssl;
  auto R_ = vgd_csrl;
  auto lsA_ = vgd_gbcf[key][0], lrA_ = vgd_gbcf[key][1];
  auto lsB_ = vgd_gbcf_xb[key][0], lrB_ = vgd_gbcf_xb[key][1];
  const int nv = vgd_n, npk = np;
  const bool lean = vgd_cur.impl_instance_id() == DevExeSpace().impl_instance_id();
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  if (vgd_ipc_state == 1) {
    // vet_gd_halo_ipc, as ScGbHaloC with the merged pieces
    for (int p = 0; p < np; ++p) {
      if (static_cast<size_t>(sm.v[p+1] - sm.v[p]) > vgd_ipc_cs ||
          static_cast<size_t>(rm.v[p+1] - rm.v[p]) > vgd_ipc_cs) {
        return false;   // no room: the two exchanges
      }
    }
    if (!lean) {return false;}
    ++vgd_ipc_seq;
    const long long seq = vgd_ipc_seq, par = seq & 1;
    auto pslot = vgd_ipc_pslot.d_view;
    auto prb = vgd_ipc_prb.d_view;
    auto phalf = vgd_ipc_phalf.d_view;
    auto pcs = vgd_ipc_pcs.d_view;
    auto psig = vgd_ipc_psig.d_view;
    if (nS > 0) {
      Kokkos::parallel_for("vet_vgdgb_ipc_xpack", ScLight(DevExeSpace(), 0, nS),
      KOKKOS_LAMBDA(const int c) {
        const bool b = (c >= nsA);
        const int cc = b ? (c - nsA) : c;
        const ScIpcIdx &sx = b ? sb : sa;
        int p = 0;
        while (p + 1 < npk && sx.v[p+1] <= cc) {++p;}
        const int q = b ? lsB_(cc) : lsA_(cc);
        const int e = q/nv;
        const int o = (b ? (sa.v[p+1] - sa.v[p]) : 0) + (cc - sx.v[p]);
        reinterpret_cast<Real *>(prb(p))[par*phalf(p) + pslot(p)*pcs(p) + o] =
            b ? aB_(S_(e,0),q - e*nv,S_(e,1),S_(e,2),i)
              : aA_(S_(e,0),q - e*nv,S_(e,1),S_(e,2),i);
      });
    }
    Kokkos::parallel_for("vet_vgdgb_ipc_set", Kokkos::RangePolicy<>(DevExeSpace(), 0, 1),
    KOKKOS_LAMBDA(const int) {ScIpcSetN(psig.data(), pslot.data(), npk, seq);});
    long long *sig = vgd_ipc_sig;
    int *err = vgd_ipc_err;
    Kokkos::parallel_for("vet_vgdgb_ipc_wait", Kokkos::RangePolicy<>(DevExeSpace(), 0, 1),
    KOKKOS_LAMBDA(const int) {ScIpcWaitN(sig, npk, seq, err);});
    for (int p = 0; p < np; ++p) {
      vgd_nbyte += static_cast<Real>(sm.v[p+1] - sm.v[p])*sizeof(Real);
      vgd_csvsent += static_cast<double>(sm.v[p+1] - sm.v[p]);
      vgd_csvdense += 2.0*static_cast<double>(vgd_gdsn[wq][p])*nv;
    }
    vgd_ipc_n += 1.0;
    vgd_nexch += 2.0;
    vgd_xff_nmsg += 1.0;
    const Real *rbp = vgd_ipc_rbuf + par*vgd_ipc_half;
    const size_t cs = vgd_ipc_cs;
    if (nR > 0) {
      Kokkos::parallel_for("vet_vgdgb_ipc_xunpack", ScLight(DevExeSpace(), 0, nR),
      KOKKOS_LAMBDA(const int c) {
        const bool b = (c >= nrA);
        const int cc = b ? (c - nrA) : c;
        const ScIpcIdx &rx = b ? rb : ra;
        int p = 0;
        while (p + 1 < npk && rx.v[p+1] <= cc) {++p;}
        const int q = b ? lrB_(cc) : lrA_(cc);
        const int e = q/nv;
        const Real x = rbp[p*cs + (b ? (ra.v[p+1] - ra.v[p]) : 0) + (cc - rx.v[p])];
        if (b) {
          aB_(R_(e,0),q - e*nv,R_(e,1),R_(e,2),i) = x;
        } else {
          aA_(R_(e,0),q - e*nv,R_(e,1),R_(e,2),i) = x;
        }
      });
    }
    return true;
  }
#endif
  if (!lean) vgd_cur.fence();   // the shell's intensities are complete
  if (nS > 0) {
    Kokkos::parallel_for("vet_vgdgb_xpack", ScLight(DevExeSpace(), 0, nS),
    KOKKOS_LAMBDA(const int c) {
      const bool b = (c >= nsA);
      const int cc = b ? (c - nsA) : c;
      const ScIpcIdx &sx = b ? sb : sa;
      int p = 0;
      while (p + 1 < npk && sx.v[p+1] <= cc) {++p;}
      const int q = b ? lsB_(cc) : lsA_(cc);
      const int e = q/nv;
      const int pos = sm.v[p] + (b ? (sa.v[p+1] - sa.v[p]) : 0) + (cc - sx.v[p]);
      sbuf_(pos) = b ? aB_(S_(e,0),q - e*nv,S_(e,1),S_(e,2),i)
                     : aA_(S_(e,0),q - e*nv,S_(e,1),S_(e,2),i);
    });
  }
  if (lean) {
    ScPackEvtRecord(DevExeSpace());
    ScPackEvtWait(DevExeSpace());
  } else {
    Kokkos::fence();
  }
  Kokkos::Timer tq;
  std::vector<MPI_Request> req;
  for (int p = 0; p < np; ++p) {
    const int rk = vgd_csprk[p];
    req.emplace_back();
    MPI_Irecv(rbuf_.data() + rm.v[p], rm.v[p+1] - rm.v[p], MPI_ATHENA_REAL, rk, 7112,
              MPI_COMM_WORLD, &req.back());
    req.emplace_back();
    MPI_Isend(sbuf_.data() + sm.v[p], sm.v[p+1] - sm.v[p], MPI_ATHENA_REAL, rk, 7112,
              MPI_COMM_WORLD, &req.back());
  }
  std::vector<MPI_Status> stat(req.size());
  MPI_Waitall(static_cast<int>(req.size()), req.data(), stat.data());
  vgd_tmpi += tq.seconds();
  vgd_nexch += 2.0;
  vgd_xff_nmsg += 1.0;
  for (int p = 0; p < np; ++p) {
    int got = 0;
    MPI_Get_count(&stat[2*p], MPI_ATHENA_REAL, &got);
    if (got != rm.v[p+1] - rm.v[p]) {
      ScFatal("clean_rot_fuse: sender and receiver kept lists disagree");
    }
    vgd_nbyte += static_cast<Real>(sm.v[p+1] - sm.v[p])*sizeof(Real);
    vgd_csvsent += static_cast<double>(sm.v[p+1] - sm.v[p]);
    vgd_csvdense += 2.0*static_cast<double>(vgd_gdsn[wq][p])*nv;
  }
  if (nR > 0) {
    Kokkos::parallel_for("vet_vgdgb_xunpack", ScLight(DevExeSpace(), 0, nR),
    KOKKOS_LAMBDA(const int c) {
      const bool b = (c >= nrA);
      const int cc = b ? (c - nrA) : c;
      const ScIpcIdx &rx = b ? rb : ra;
      int p = 0;
      while (p + 1 < npk && rx.v[p+1] <= cc) {++p;}
      const int q = b ? lrB_(cc) : lrA_(cc);
      const int e = q/nv;
      const Real x = rbuf_(rm.v[p] + (b ? (ra.v[p+1] - ra.v[p]) : 0) + (cc - rx.v[p]));
      if (b) {
        aB_(R_(e,0),q - e*nv,R_(e,1),R_(e,2),i) = x;
      } else {
        aA_(R_(e,0),q - e*nv,R_(e,1),R_(e,2),i) = x;
      }
    });
  }
  if (!lean) Kokkos::fence();
  return true;
#else
  (void) i;
  (void) ws;
  return false;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationVet::ScRawKey
//! \brief vgd_raw + vet_gd_band_ring: key and record of the band positions the raw
//! exchange of (shell, pass, depth) writes, encoded ((m n + v) c3 + k) c2 + j; the ring
//! saves / restores them (ScRingFlush)

int RadiationVet::ScRawKey(const int i, const bool inw, const int ws) const {
  const int wq = std::min(ws, static_cast<int>(vgd_gdnl.size()) - 1);
  return (2*i + (inw ? 1 : 0))*(vgd_w + 1) + wq;
}

DvceArray1D<int> RadiationVet::ScRawRec(const int key) {
  auto it = vgd_gbrec.find(key);
  if (it != vgd_gbrec.end()) {return it->second;}
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nv = vgd_n, wq = key % (vgd_w + 1), np = static_cast<int>(vgd_csprk.size());
  const int c3 = indcs.nx3 + 2*vgd_w, c2 = indcs.nx2 + 2*vgd_w;
  if (static_cast<double>(pmy_pack->nmb_thispack)*nv*c3*c2 >= INT32_MAX) {
    ScFatal("vet_gd_band_ring: raw band index >= 2^31");
  }
  const bool kc = vgd_gbc && vgd_gbcb.count(key) && vgd_gbcb[key].size() > 4*np + 1;
  auto R_ = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_csrl);
  auto D_ = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_gdpr);
  std::vector<int> q;   // kept (entry, direction), else all within depth wq
  if (kc) {
    q.resize(vgd_gbcb[key][4*np + 1]);
    auto lh = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_gbcf[key][1]);
    for (size_t o = 0; o < q.size(); ++o) {q[o] = lh(o);}
  } else {
    for (int e = 0; e < vgd_csnrcv; ++e) {
      for (int v = 0; v < nv && D_(e) <= wq; ++v) {q.push_back(e*nv + v);}
    }
  }
  DvceArray1D<int> l("vet_vgd_rawrec", q.size());
  auto l_h = Kokkos::create_mirror_view(l);
  for (size_t o = 0; o < q.size(); ++o) {
    const int e = q[o]/nv, v = q[o] - e*nv;
    l_h(o) = ((R_(e,0)*nv + v)*c3 + R_(e,1))*c2 + R_(e,2);
  }
  Kokkos::deep_copy(l, l_h);
  vgd_gbrec[key] = l;
  return l;
}

template bool RadiationVet::ScGbHaloC<DvceArray5D<Real>>(DvceArray5D<Real> &a,
                                                         const int i, const int ws);

template bool RadiationVet::ScGbHaloC<RadiationVet::ScRag>(RadiationVet::ScRag &a,
                                                           const int i, const int ws);

}  // namespace radvet
