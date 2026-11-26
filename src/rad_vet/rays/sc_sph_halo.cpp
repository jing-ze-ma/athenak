//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_sph_halo.cpp
//! \brief spherical-polar rays: halo of a shell, by masks and lists.
//!
//!  called from: rays/sc_build.cpp, rays/sc_sph_sweep.cpp, rays/sc_sph_walls.cpp
//!  calls into: rays/sc_bandring.cpp, rays/sc_cs_halo.cpp, rays/sc_sph_globalband.cpp,
//!    rays/sc_sph_halo_exact.cpp, rays/sc_sph_walls.cpp

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include "athena.hpp"
// FMA contraction only within one expression in this TU
#if defined(__clang__)
#pragma clang fp contract(on)
#endif
#include "mesh/mesh.hpp"
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
//! \fn void RadiationVet::VetGdHalo
//! \brief fill the lateral band (vgd_w cells) of a(m, v, k, j, i), v in [0, nv), i in
//! [i0, i1], from the 8 lateral neighbours' interiors: local ones copied directly (read
//! interiors, write bands: no race), remote ones by one message per (block, slot).
//! vgd_raw: global-cell lists, unmapped; intensities (mapd) only from remote owners

template <class V>
void RadiationVet::VetGdHalo(V &a, const int nv, const int i0, const int i1, const int ws,
                             const bool mapd) {
  // vet_gd_band_global: every band position within depth ws from its owner's interior
  // (global-cell lists, sph2.inc); no direction re-indexing at the full sphere's walls
  if (vgd_gband || vgd_raw) {
    // compact per-shell exchange of the global band
    if (mapd && vgd_gbc && i1 == i0 && nv == vgd_n) {
      VetGdGbHaloC(a, i0, ws);
      return;
    }
    VetGdCsHalo(a, nv, i0, i1, !(vgd_raw && mapd), ws);
    return;
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int ni = i1 - i0 + 1;
  const int mx = std::max(nx2, nx3)*ws*nv*ni;
  auto hl_ = vgd_hloc;
  auto hp_ = vgd_hpol;
  auto a_ = a;
  auto wl_ = vgd_wall;
  auto mp_ = vgd_map;
  const bool wint = vgd_wint;
  auto m3_ = vgd_m3;
  auto w3_ = vgd_w3;
  // region of slot o along one direction d = -1, 0, +1 with depth ws: (dest start, len,
  // source start) = (w - ws, ws, w + nx - ws), (w, nx, w), (w + nx, ws, w); wall
  // ghosts (mapd) take direction mp(v) of the source (new(d) = old(map(d)))
  const bool mpi = vgd_hmpi;
  if (mapd && !mpi) {return;}
  auto sb_ = vgd_sbuf;
  auto so_ = vgd_soff;
  auto ro_ = vgd_roff;
  const int nvi = nv*ni;
  const bool hcomp = mapd && (vgd_hcomp > 0) && (ni == 1);
  // true when the compact halo's list path wrote the band itself
  bool done = false;
  if (hcomp) {done = VetGdHaloCompact(a, nv, i0, ws);}
  if (!hcomp) {
  par_for("vet_vgd_halo", vgd_cur, 0, nmb1, 0, 7, 0, mx - 1,
  KOKKOS_LAMBDA(const int m, const int o, const int t) {
    const int oo = (o < 4) ? o : (o + 1);
    const int dk = oo/3 - 1, dj = oo%3 - 1;
    const int jn = (dj == 0) ? nx2 : ws;
    const int kn = (dk == 0) ? nx3 : ws;
    const int cnt = nv*kn*jn*ni;
    if (t >= cnt) {return;}
    const int v = t/(kn*jn*ni);
    const int r1 = t - v*kn*jn*ni;
    const int kk = r1/(jn*ni);
    const int r2 = r1 - kk*jn*ni;
    const int jj = r2/ni;
    const int i = i0 + (r2 - jj*ni);
    const int nl = hl_(8*m + o);
    // a slot beyond a pole (full sphere): the partner's rows reversed
    const int pl = hp_(8*m + o);
    if (nl >= 0 && mapd) {return;}    // the sweep reads on-rank neighbours directly
    if (nl >= 0) {
      const int jd = (dj < 0) ? (w - ws) : ((dj == 0) ? w : (w + nx2));
      const int kd = (dk < 0) ? (w - ws) : ((dk == 0) ? w : (w + nx3));
      const int js_ = (dj < 0) ? (w + nx2 - ws) : w;
      // nx3 = 1: every phi band column is a copy of the one column
      const int ks_ = (nx3 == 1) ? w : (((dk < 0) ? (w + nx3 - ws) : w) + kk);
      const int vs = (mapd && wl_(m,kd+kk,jd+jj) != 0) ? mp_(m,kd+kk,jd+jj,v) : v;
      const int jsrc = (pl != 0) ? VgdPoleRow(jd + jj, w, nx2) : (js_ + jj);
      a_(m,v,kd+kk,jd+jj,i) = a_(nl,vs,ks_,jsrc,i);
    } else if (mpi) {
      // my interior that the remote neighbour needs in its opposite slot (across a
      // pole: my rows reversed)
      const int js2 = (dj > 0) ? (w + nx2 - ws) : w;
      const int ks2 = (nx3 == 1) ? w : (((dk > 0) ? (w + nx3 - ws) : w) + kk);
      const int jsrc = (pl == 0) ? (js2 + jj)
                       : ((pl == 1) ? (w + ws - 1 - jj) : (w + nx2 - 1 - jj));
      sb_(static_cast<size_t>(so_(ws,8*m + o))*nvi + t) = a_(m,v,ks2,jsrc,i);
    }
  });
  }  // !hcomp
#if MPI_PARALLEL_ENABLED
  if (mpi && !hcomp) {
    vgd_cur.fence();
    Kokkos::Timer tq;
    std::vector<MPI_Request> req;
    auto rb_ = vgd_rbuf;
    // one message per partner rank and direction
    for (size_t p = 0; p < vgd_prk.size(); ++p) {
      const int rk = vgd_prk[p];
      req.emplace_back();
      MPI_Irecv(rb_.data() + static_cast<size_t>(vgd_prdsp[ws][p])*nvi,
                vgd_prcnt[ws][p]*nvi, MPI_ATHENA_REAL, rk, 7001, vgd_comm,
                &req.back());
      req.emplace_back();
      MPI_Isend(sb_.data() + static_cast<size_t>(vgd_pdsp[ws][p])*nvi,
                vgd_pscnt[ws][p]*nvi, MPI_ATHENA_REAL, rk, 7001, vgd_comm,
                &req.back());
    }
    vgd_tpost += tq.seconds();
    MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
    vgd_tmpi += tq.seconds();
    vgd_nexch += 1.0;
    for (size_t p = 0; p < vgd_prk.size(); ++p) {
      vgd_nbyte += static_cast<Real>(vgd_pscnt[ws][p])*nvi*sizeof(Real);
    }
  }
  if (mpi && !done) {
    auto rb_ = vgd_rbuf;
    par_for("vet_vgd_unpack", vgd_cur, 0, nmb1, 0, 7, 0, mx - 1,
    KOKKOS_LAMBDA(const int m, const int o, const int t) {
      if (hl_(8*m + o) >= 0) {return;}
      const int oo = (o < 4) ? o : (o + 1);
      const int dk = oo/3 - 1, dj = oo%3 - 1;
      const int jd = (dj < 0) ? (w - ws) : ((dj == 0) ? w : (w + nx2));
      const int kd = (dk < 0) ? (w - ws) : ((dk == 0) ? w : (w + nx3));
      const int jn = (dj == 0) ? nx2 : ws, kn = (dk == 0) ? nx3 : ws;
      const int cnt = nv*kn*jn*ni;
      if (t >= cnt) {return;}
      const int v = t/(kn*jn*ni);
      const int r1 = t - v*kn*jn*ni;
      const int kk = r1/(jn*ni);
      const int r2 = r1 - kk*jn*ni;
      const int jj = r2/ni;
      const int ii = r2 - jj*ni;
      const size_t rb0 = static_cast<size_t>(ro_(ws,8*m + o))*nvi;
      if (mapd && wint && wl_(m,kd+kk,jd+jj) != 0) {
        Real val = 0.0;
        for (int q = 0; q < 3; ++q) {
          const int vq = m3_(m,kd+kk,jd+jj,v,q);
          const size_t at = rb0 + ((vq*kn + kk)*jn + jj)*ni + ii;
          val += w3_(m,kd+kk,jd+jj,v,q)*rb_(at);
        }
        a_(m,v,kd+kk,jd+jj,i0+ii) = val;
        return;
      }
      const int vs = (mapd && wl_(m,kd+kk,jd+jj) != 0) ? mp_(m,kd+kk,jd+jj,v) : v;
      const size_t at = rb0 + ((vs*kn + kk)*jn + jj)*ni + ii;
      a_(m,v,kd+kk,jd+jj,i0+ii) = rb_(at);
    });
  }
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdHaloCs
//! \brief VetGdHalo of ln chi, ln S over shells [i0, i1] at full band depth, in chunks
//! that fit vgd_maxcnt per slot (bitwise the one-call exchange)

void RadiationVet::VetGdHaloCs(DvceArray5D<Real> &a, const int i0, const int i1) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int per = std::max(indcs.nx2, indcs.nx3)*vgd_w*2;   // one shell, one slot
  const int nch = std::max(1, (per > 0) ? vgd_maxcnt/per : 1);
  for (int c0 = i0; c0 <= i1; c0 += nch) {
    VetGdHalo(a, 2, c0, std::min(c0 + nch - 1, i1), vgd_w, false);
  }
  if (vgd_thmir) {VetGdThMirCs(a, i0, i1);}
  if (vgd_thvac > 0) {VetGdThVacCs(a, i0, i1);}
}

namespace {

KOKKOS_INLINE_FUNCTION
bool VgdHaloKeep(const Real nx, const Real ny, const Real nz, const int gj, const int gk,
                 const Real t0, const Real dth, const Real p0, const Real dph,
                 const bool inw, const int mode, const bool wall, const int dep,
                 const Real rs, const Real rm, const Real rp, const bool hasm,
                 const bool hasp, const Real gam, const Real dmin) {
  const Real th = t0 + (gj + 0.5)*dth, ph = p0 + (gk + 0.5)*dph;
  return VgdHaloKeepTP(nx, ny, nz, th, ph, inw, mode, wall, dep, rs, rm, rp, hasm, hasp,
                       gam, dmin);
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdHcPrep
//! \brief vet_gd_halo_compact = 1, 2: the keep masks of (shell, pass, depth) on the
//! dense layouts of both sides (VgdHaloKeep), their scans and partner bounds
void RadiationVet::VetGdHcPrep(const int slot, const int i0, const bool inw0,
                               const int ws) {
#if MPI_PARALLEL_ENABLED
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int nv = vgd_n;
  const int mx = std::max(nx2, nx3)*ws*nv;
  const int i = i0;
  const bool inw = inw0;
  const int mode = vgd_hcomp;
  const int nb2 = vgd_nb2, nb3 = vgd_nb3;
  // mirror / extrap theta walls: rows within ws of a wall keep their branch (gj global)
  const bool thm = vgd_thmir || (vgd_thvac == 2);
  const int nx2m = pm->mesh_indcs.nx2;
  const Real t0 = pm->mesh_size.x2min, p0 = pm->mesh_size.x3min;
  const Real dth = (pm->mesh_size.x2max - t0)/pm->mesh_indcs.nx2;
  const Real dph = (pm->mesh_size.x3max - p0)/pm->mesh_indcs.nx3;
  const Real smin = fmin(sin(t0), sin(pm->mesh_size.x2max));
  const Real dmin = fmin(dth, smin*dph);
  const Real gam = 2.0*fmax(dth, dph);
  const int lcut = indcs.is + vgd_scut;
  const Real rs = vgd_r1v[i];
  const bool hasm = (i >= lcut);
  const Real rm = (i > lcut) ? vgd_r1v[i-1] : vgd_r1f[i];
  const bool hasp = (i < indcs.ie);
  const Real rp = hasp ? vgd_r1v[i+1] : rs;
  auto hl_ = vgd_hloc;
  auto lx_ = vgd_lxy;
  auto dir_ = vgd_dir;
  auto so_ = vgd_soff;
  auto ro_ = vgd_roff;
  auto fs_ = vgd_hfs[slot];
  auto fr_ = vgd_hfr[slot];
  const int nvi = nv;
  const size_t stot = static_cast<size_t>(vgd_pdsp[ws].empty() ? 0 :
                      (vgd_pdsp[ws].back() + vgd_pscnt[ws].back()))*nvi;
  const size_t rtot = static_cast<size_t>(vgd_prdsp[ws].empty() ? 0 :
                      (vgd_prdsp[ws].back() + vgd_prcnt[ws].back()))*nvi;
  // (1) flags of both dense layouts (sender: my interior; receiver: my band)
  par_for("vet_vgd_hc_flag", vgd_cur, 0, nmb1, 0, 7, 0, mx - 1,
  KOKKOS_LAMBDA(const int m, const int o, const int t) {
    if (hl_(8*m + o) >= 0) {return;}
    const int oo = (o < 4) ? o : (o + 1);
    const int dk = oo/3 - 1, dj = oo%3 - 1;
    const int jn = (dj == 0) ? nx2 : ws;
    const int kn = (dk == 0) ? nx3 : ws;
    if (t >= nv*kn*jn) {return;}
    const int v = t/(kn*jn);
    const int r1 = t - v*kn*jn;
    const int kk = r1/jn;
    const int jj = r1 - kk*jn;
    const int lx2 = lx_(m,0), lx3 = lx_(m,1);
    const Real nx = dir_(v,0), ny = dir_(v,1), nz = dir_(v,2);
    {
      // send (my slot o, receiver slot direction -dj, -dk)
      const int js2 = (dj > 0) ? (w + nx2 - ws) : w;
      const int ks2 = (dk > 0) ? (w + nx3 - ws) : w;
      const int gj = lx2*nx2 + (js2 + jj - w), gk = lx3*nx3 + (ks2 + kk - w);
      const int rdj = -dj, rdk = -dk;
      const int dpj = (rdj < 0) ? (ws - jj) : ((rdj > 0) ? (jj + 1) : 0);
      const int dpk = (rdk < 0) ? (ws - kk) : ((rdk > 0) ? (kk + 1) : 0);
      const int lr2 = (lx2 + dj + nb2) % nb2, lr3 = (lx3 + dk + nb3) % nb3;
      // mirror theta walls: rows within ws of a wall are read by mirrored rays too;
      // keep their whole branch (decided from gj alone)
      const bool wall = (lr2 + rdj < 0) || (lr2 + rdj >= nb2) || (lr3 + rdk < 0)
                        || (lr3 + rdk >= nb3) || (thm && (gj < ws || gj >= nx2m - ws));
      fs_(static_cast<size_t>(so_(ws,8*m + o))*nvi + t) =
        VgdHaloKeep(nx, ny, nz, gj, gk, t0, dth, p0, dph, inw, mode, wall,
                    (dpj > dpk) ? dpj : dpk, rs, rm, rp, hasm, hasp, gam, dmin) ? 1 : 0;
    }
    {
      // receive (my slot o): the source cell in the neighbour, global indices
      const int js_ = (dj < 0) ? (w + nx2 - ws) : w;
      const int ks_ = (dk < 0) ? (w + nx3 - ws) : w;
      const int ln2 = (lx2 + dj + nb2) % nb2, ln3 = (lx3 + dk + nb3) % nb3;
      const int gj = ln2*nx2 + (js_ + jj - w), gk = ln3*nx3 + (ks_ + kk - w);
      const int dpj = (dj < 0) ? (ws - jj) : ((dj > 0) ? (jj + 1) : 0);
      const int dpk = (dk < 0) ? (ws - kk) : ((dk > 0) ? (kk + 1) : 0);
      const bool wall = (lx2 + dj < 0) || (lx2 + dj >= nb2) || (lx3 + dk < 0)
                        || (lx3 + dk >= nb3) || (thm && (gj < ws || gj >= nx2m - ws));
      fr_(static_cast<size_t>(ro_(ws,8*m + o))*nvi + t) =
        VgdHaloKeep(nx, ny, nz, gj, gk, t0, dth, p0, dph, inw, mode, wall,
                    (dpj > dpk) ? dpj : dpk, rs, rm, rp, hasm, hasp, gam, dmin) ? 1 : 0;
    }
  });
  // (2) exclusive scans in place (flag -> position, total at n); flag = pos(q+1) -
  // pos(q); both in one launch (two 32-bit counts in one 64-bit sum)
  VgdScan2(vgd_cur, vgd_hfs[slot], stot, vgd_hfr[slot], rtot, "vet_vgd_hc_scan2",
           vgd_scan_sums);
  // (3) scan values at the partner boundaries, straight into pinned host memory (read
  // after the exchange's fence)
  const int np = static_cast<int>(vgd_prk.size());
  {
    auto pbd_ = vgd_pbd;
    auto pbv_ = vgd_hpb[slot];
    Kokkos::parallel_for("vet_vgd_hc_bnd",
                         Kokkos::RangePolicy<>(vgd_cur, 0, 2*(np + 1)),
    KOKKOS_LAMBDA(const int q) {
      const size_t at = static_cast<size_t>(pbd_(ws,q))*nvi;
      pbv_(q) = (q <= np) ? fs_(at) : fr_(at);
    });
  }
  vgd_htag[slot][0] = i0;
  vgd_htag[slot][1] = inw0 ? 1 : 0;
  vgd_htag[slot][2] = ws;
  vgd_htag[slot][3] = vgd_scut;
  vgd_htag[slot][4] = vgd_hsweep;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::VetGdHlBuild
//! \brief vet_gd_halo_list_mb: entry e (pass inw, shell i, depth ws) from the compact
//! masks fs / fr and partner bounds pbh:
//!   sl     band address (m, v, k, j) of the c-th sent value,
//!   dl, sr band address of every receiver ghost with a sent source, and its position,
//!   dw, sw the same for vet_gd_wall_interp wall ghosts (3 sources, -1 = unsent).
//! Address = ((m nv + v) c3 + k) c2 + j at depth vgd_w. False when over the budget.

bool RadiationVet::VetGdHlBuild(const int e, const int i, const bool inw, const int ws,
                                const DvceArray1D<int> &fs, const DvceArray1D<int> &fr,
                                const std::vector<int> &pbh) {
#if MPI_PARALLEL_ENABLED
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int nv = vgd_n;
  const int mx = std::max(nx2, nx3)*ws*nv;
  const int c2 = nx2 + 2*w, c3 = nx3 + 2*w;
  const int np = static_cast<int>(vgd_prk.size());
  VgdHlEntry &E = vgd_hle[e];
  E.gen = -1;
  const size_t rtot = static_cast<size_t>(vgd_prdsp[ws].empty() ? 0 :
                      (vgd_prdsp[ws].back() + vgd_prcnt[ws].back()))*nv;
  if (vgd_hl_fn.extent(0) < rtot + 1) {
    Kokkos::realloc(vgd_hl_fn, rtot + 1);
    Kokkos::realloc(vgd_hl_fw, rtot + 1);
  }
  auto hl_ = vgd_hloc;
  auto so_ = vgd_soff;
  auto ro_ = vgd_roff;
  auto wl_ = vgd_wall;
  auto mp_ = vgd_map;
  auto m3_ = vgd_m3;
  const bool wint = vgd_wint;
  auto fs_ = fs;
  auto fr_ = fr;
  auto fn_ = vgd_hl_fn;
  auto fw_ = vgd_hl_fw;
  // (1) receiver ghosts that the dense unpack fills from at least one sent value
  par_for("vet_vgd_hl_dflag", vgd_cur, 0, nmb1, 0, 7, 0, mx - 1,
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
    int fn = 0, fw = 0;
    if (wint && wl_(m,kd+kk,jd+jj) != 0) {
      for (int u = 0; u < 3; ++u) {
        const size_t at = rb0 + (m3_(m,kd+kk,jd+jj,v,u)*kn + kk)*jn + jj;
        if (fr_(at + 1) > fr_(at)) {fw = 1;}
      }
    } else {
      const int vs = (wl_(m,kd+kk,jd+jj) != 0) ? mp_(m,kd+kk,jd+jj,v) : v;
      const size_t at = rb0 + (vs*kn + kk)*jn + jj;
      fn = (fr_(at + 1) > fr_(at)) ? 1 : 0;
    }
    fn_(q) = fn;
    fw_(q) = fw;
  });
  const DevExeSpace ex_ = vgd_cur;
  VgdScan2(ex_, vgd_hl_fn, rtot, vgd_hl_fw, rtot, "vet_vgd_hl_scan2", vgd_scan_sums);
  int nrw[2] = {0, 0};
  {
    Kokkos::View<int*, HostMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> h(nrw, 1);
    Kokkos::deep_copy(ex_, h, Kokkos::subview(vgd_hl_fn, std::make_pair(rtot, rtot + 1)));
    Kokkos::View<int*, HostMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> h2(nrw + 1,
                                                                                 1);
    Kokkos::deep_copy(ex_, h2,
                      Kokkos::subview(vgd_hl_fw, std::make_pair(rtot, rtot + 1)));
    ex_.fence();
  }
  const int ns = pbh[np] - pbh[0];
  const int nr = nrw[0], nw = nrw[1];
  // every array at least its new size (capacity kept); budget on the sum
  const size_t nsz[5] = {static_cast<size_t>(std::max(ns, 1)),
                         static_cast<size_t>(std::max(nr, 1)),
                         static_cast<size_t>(std::max(nr, 1)),
                         static_cast<size_t>(std::max(nw, 1)),
                         static_cast<size_t>(std::max(3*nw, 1))};
  DvceArray1D<int> *arr[5] = {&E.sl, &E.dl, &E.sr, &E.dw, &E.sw};
  size_t have = 0, after = 0;
  for (int q = 0; q < 5; ++q) {
    have += arr[q]->extent(0);
    after += std::max(arr[q]->extent(0), nsz[q]);
  }
  if (after > have) {
    // budget: vet_gd_halo_list_mb, capped on a GPU at half the device memory free at
    // the first build
    if (vgd_hl_cap == 0) {
      vgd_hl_cap = static_cast<size_t>(vgd_hl_mb) << 20;
      size_t fr = 0, tt = 0;
#if defined(KOKKOS_ENABLE_CUDA)
      if (cudaMemGetInfo(&fr, &tt) == cudaSuccess) {
        vgd_hl_cap = std::min(vgd_hl_cap, fr/2);
      }
#elif defined(KOKKOS_ENABLE_HIP)
      if (hipMemGetInfo(&fr, &tt) == hipSuccess) {
        vgd_hl_cap = std::min(vgd_hl_cap, fr/2);
      }
#endif
      (void) fr;
      (void) tt;
    }
    const size_t cap = vgd_hl_cap;
    if (vgd_hl_bytes + sizeof(int)*(after - have) > cap) {
      // over the budget: no new attempt for this mask set (dense path)
      E.alpha = vgd_alpha;
      E.scut = vgd_scut;
      E.ws = ws;
      E.nofit = true;
      return false;
    }
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
  // (2) the send list (compact pack addresses)
  par_for("vet_vgd_hl_slist", vgd_cur, 0, nmb1, 0, 7, 0, mx - 1,
  KOKKOS_LAMBDA(const int m, const int o, const int t) {
    if (hl_(8*m + o) >= 0) {return;}
    const int oo = (o < 4) ? o : (o + 1);
    const int dk = oo/3 - 1, dj = oo%3 - 1;
    const int jn = (dj == 0) ? nx2 : ws;
    const int kn = (dk == 0) ? nx3 : ws;
    if (t >= nv*kn*jn) {return;}
    const size_t q = static_cast<size_t>(so_(ws,8*m + o))*nv + t;
    if (fs_(q + 1) == fs_(q)) {return;}
    const int v = t/(kn*jn);
    const int r1 = t - v*kn*jn;
    const int kk = r1/jn;
    const int jj = r1 - kk*jn;
    const int js2 = (dj > 0) ? (w + nx2 - ws) : w;
    const int ks2 = (dk > 0) ? (w + nx3 - ws) : w;
    sl_(fs_(q)) = ((m*nv + v)*c3 + ks2 + kk)*c2 + js2 + jj;
  });
  // (3) the receive lists
  par_for("vet_vgd_hl_rlist", vgd_cur, 0, nmb1, 0, 7, 0, mx - 1,
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
    if (fw_(q + 1) > fw_(q)) {
      const int c = fw_(q);
      dw_(c) = ad;
      for (int u = 0; u < 3; ++u) {
        const size_t at = rb0 + (m3_(m,kd+kk,jd+jj,v,u)*kn + kk)*jn + jj;
        sw_(3*c + u) = (fr_(at + 1) > fr_(at)) ? fr_(at) : -1;
      }
    } else if (fn_(q + 1) > fn_(q)) {
      const int vs = (wl_(m,kd+kk,jd+jj) != 0) ? mp_(m,kd+kk,jd+jj,v) : v;
      const size_t at = rb0 + (vs*kn + kk)*jn + jj;
      dl_(fn_(q)) = ad;
      sr_(fn_(q)) = fr_(at);
    }
  });
  E.sp.assign(pbh.begin(), pbh.begin() + np + 1);
  E.rp.assign(pbh.begin() + np + 1, pbh.begin() + 2*np + 2);
  E.ns = ns;
  E.nr = nr;
  E.nw = nw;
  E.alpha = vgd_alpha;
  E.scut = vgd_scut;
  E.ws = ws;
  E.gen = ++vgd_hl_gen;
  E.nofit = false;
  vgd_hl_nmade += 1.0;
  return true;
#else
  return false;
#endif
}

//----------------------------------------------------------------------------------------
// an event behind the per-shell halo pack on the sweep stream, so the host
// waits for the pack alone while the ring copies queued after it run (VetGdHlBegin)
#if defined(KOKKOS_ENABLE_HIP)
hipEvent_t g_vgd_pkev = nullptr;
void VgdPackEvtRecord(const DevExeSpace &ex) {
  if (g_vgd_pkev == nullptr) {
    (void)hipEventCreateWithFlags(&g_vgd_pkev, hipEventDisableTiming);
  }
  (void)hipEventRecord(g_vgd_pkev, ex.hip_stream());
}
void VgdPackEvtWait(const DevExeSpace &) {(void)hipEventSynchronize(g_vgd_pkev);}
#elif defined(KOKKOS_ENABLE_CUDA)
cudaEvent_t g_vgd_pkev = nullptr;
void VgdPackEvtRecord(const DevExeSpace &ex) {
  if (g_vgd_pkev == nullptr) {
    (void)cudaEventCreateWithFlags(&g_vgd_pkev, cudaEventDisableTiming);
  }
  (void)cudaEventRecord(g_vgd_pkev, ex.cuda_stream());
}
void VgdPackEvtWait(const DevExeSpace &) {(void)cudaEventSynchronize(g_vgd_pkev);}
#else
void VgdPackEvtRecord(const DevExeSpace &) {}
void VgdPackEvtWait(const DevExeSpace &ex) {ex.fence();}
#endif

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::VetGdHlBegin
//! \brief first half of the list path of VetGdHaloCompact (gather pack, fence, posted
//! messages, zeros on Y's entries); false when the list path does not apply.
//! VetGdHlEnd completes it; the sweep overlap runs local items in between.

template <class V>
bool RadiationVet::VetGdHlBegin(V &a, const int nv, const int i0, const int ws) {
#if MPI_PARALLEL_ENABLED
  if (vgd_hlp.on) {VgdFatal("VetGdHlBegin: an exchange is still pending");}
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int i = i0;
  const bool inw = vgd_hinw;
  auto a_ = a;
  auto csb_ = vgd_csb;
  auto crb_ = vgd_crb;
  const int np = static_cast<int>(vgd_prk.size());
  // vet_gd_halo_list_mb: list entries X (this pass) and Y (the other pass) of shell i,
  // and the last writers of the band per pass
  const int n1c = indcs.nx1 + 2*indcs.ng;
  const int c2 = nx2 + 2*w, c3 = nx3 + 2*w;
  bool lst = (nv == vgd_n) && (i >= 0) && (i < n1c) &&
             (static_cast<int64_t>(nmb1 + 1)*nv*c3*c2 <
              static_cast<int64_t>(std::numeric_limits<int>::max()));
  const int px = inw ? 0 : 1;
  const int ex = px*n1c + i, ey = (1 - px)*n1c + i;
  std::array<std::vector<int64_t>, 2> *gwa = nullptr;
  if (lst) {
    if (static_cast<int>(vgd_hle.size()) != 2*n1c) {vgd_hle.resize(2*n1c);}
    auto gwget = [&](const void *p) {
      auto &g = vgd_hl_gw[p];
      if (static_cast<int>(g[0].size()) != n1c) {
        g[0].assign(n1c, 0);
        g[1].assign(n1c, 0);
      }
      return &g;
    };
    if constexpr (std::is_same<V, VgdIView>::value) {
      gwa = gwget(static_cast<const void *>(a.d.data()));
    } else {
      lst = false;
    }
  }
  auto xvalid = [&]() {
    const VgdHlEntry &X = vgd_hle[ex];
    return X.gen > 0 && X.alpha == vgd_alpha && X.scut == vgd_scut && X.ws == ws;
  };
  auto gwok = [&](std::array<std::vector<int64_t>, 2> *g) {
    const int64_t gx = (*g)[px][i], gy = (*g)[1 - px][i];
    return (gx == 0 || gx == vgd_hle[ex].gen) && (gy == 0 || (vgd_hle[ey].gen > 0 &&
                                                               gy == vgd_hle[ey].gen));
  };
  // vet_gd_halo_exact: only an exact record of this key (no last-writer condition, no
  // zeros: every ghost the sweep reads is on the list)
  const bool xe = vgd_hx_on;
  if (xe) {
    if (!(lst && VetGdHxKey(ex, ws) && vgd_hxe[ex].exact)) {return false;}
  } else if (!(lst && xvalid() && gwok(gwa))) {
    return false;
  }
  if constexpr (std::is_same<V, VgdIView>::value) {
    const VgdHlEntry &X = xe ? vgd_hxe[ex] : vgd_hle[ex];
    auto sl_ = X.sl;
    const int c2_ = c2, c3_ = c3, nv_ = nv;
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    // vet_gd_halo_ipc: node-local and GPU-driven: wait for the partners' acks of this
    // parity, gather my values into their receive buffers, set their arrival flags,
    // wait for mine; all on the sweep stream (bitwise the MPI path)
    if (xe && vgd_ipc_state == 1 && static_cast<int>(X.dsto.size()) == np) {
      const long long seq = ++vgd_ipc_seq;
      const long long par = seq & 1;
      long long *sig = vgd_ipc_sig;
      int *err = vgd_ipc_err;
      auto psig = vgd_ipc_psig.d_view;
      auto pslot = vgd_ipc_pslot.d_view;
      auto prb = vgd_ipc_prb.d_view;
      auto phalf = vgd_ipc_phalf.d_view;
      const int np_ = np;
      // lean protocol, no ack kernel: partner p sets my arrival flag of exchange s in
      // its k2(s), after its scatter of s-1; so once my k2(s) has seen it, p has read
      // its parity (s-1)&1 buffer and my pack of s+1 may write there
      VgdIpcIdx spi, doi;
      for (int p = 0; p <= np; ++p) {spi.v[p] = X.sp[p];}
      for (int p = 0; p < np; ++p) {doi.v[p] = X.dsto[p];}
      // fused exchange (CUDA lean, no band-ring copy in between): teams pack as
      // vet_vgd_ipc_pack, the last one sets the partners' flags, every team waits for
      // mine and scatters as vet_vgd_hl_scatter (bitwise). At most 216 teams of 128: all
      // resident, so early spinning teams cannot starve the last one's pack
#if defined(KOKKOS_ENABLE_CUDA)
      if (vgd_hlp.now && vgd_rg_pend.empty() && X.ns + X.nr + X.nw > 0) {
        if (!vgd_hlp.cnt.is_allocated()) {
          vgd_hlp.cnt = Kokkos::View<int, DevMemSpace>("vet_vgd_ipc_cnt");
        }
        auto cnt_ = vgd_hlp.cnt;
        long long *sig = vgd_ipc_sig;
        int *err = vgd_ipc_err;
        const Real *rb_ = vgd_ipc_rbuf + par*vgd_ipc_half;
        const int ns_ = X.ns, nr = X.nr, nw = X.nw;
        auto dl_ = X.dl;
        auto sr_ = X.sr;
        auto dw_ = X.dw;
        auto sw_ = X.sw;
        auto w3_ = vgd_w3;
        constexpr int kT = 128;
        const int nt = std::max(ns_, nr + nw);
        const int nl = std::min((nt + kT - 1)/kT, 216);
        // light-weight launch: body as kernel arguments, no host wait
        Kokkos::parallel_for("vet_vgd_ipc_x", Kokkos::Experimental::require(
            Kokkos::TeamPolicy<>(vgd_cur, nl, kT),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight),
        KOKKOS_LAMBDA(const TeamMember_t &tm) {
          const int c0 = tm.league_rank()*kT + tm.team_rank(), dc = tm.league_size()*kT;
          for (int c = c0; c < ns_; c += dc) {
            int p = 0;
            while (p + 1 < np_ && spi.v[p+1] <= c) {++p;}
            int ad = sl_(c);
            const int jb = ad % c2_;
            ad /= c2_;
            const int kb = ad % c3_;
            ad /= c3_;
            const int v = ad % nv_;
            const int m = ad/nv_;
            Real *dst = reinterpret_cast<Real *>(prb(p)) + par*phalf(p) + doi.v[p] +
                        (c - spi.v[p]);
            *dst = a_(m,v,kb,jb,i);
          }
#if defined(__CUDA_ARCH__)
          __threadfence_system();
#endif
          tm.team_barrier();
          if (tm.team_rank() == 0) {
            if (Kokkos::atomic_fetch_add(&cnt_(), 1) == tm.league_size() - 1) {
              cnt_() = 0;
              VgdIpcSetN(psig.data(), pslot.data(), np_, seq);
            }
            VgdIpcWaitN(sig, np_, seq, err);
          }
          tm.team_barrier();
          for (int c = c0; c < nr + nw; c += dc) {
            int ad = (c < nr) ? dl_(c) : dw_(c - nr);
            const int jb = ad % c2_;
            ad /= c2_;
            const int kb = ad % c3_;
            ad /= c3_;
            const int v = ad % nv_;
            const int m = ad/nv_;
            if (c < nr) {
              const int s = sr_(c);
              a_(m,v,kb,jb,i) = (s >= 0) ? rb_[s] : 0.0;
              continue;
            }
            const int cw = c - nr;
            Real val = 0.0;
            for (int u = 0; u < 3; ++u) {
              const int s = sw_(3*cw + u);
              val += w3_(m,kb,jb,v,u)*((s >= 0) ? rb_[s] : 0.0);
            }
            a_(m,v,kb,jb,i) = val;
          }
        });
        vgd_hlp.done = true;
      } else
#endif
      if (X.ns > 0) {
        Kokkos::parallel_for("vet_vgd_ipc_pack", VgdLight(vgd_cur, 0, X.ns),
        KOKKOS_LAMBDA(const int c) {
          int p = 0;
          while (p + 1 < np_ && spi.v[p+1] <= c) {++p;}
          int ad = sl_(c);
          const int jb = ad % c2_;
          ad /= c2_;
          const int kb = ad % c3_;
          ad /= c3_;
          const int v = ad % nv_;
          const int m = ad/nv_;
          Real *dst = reinterpret_cast<Real *>(prb(p)) + par*phalf(p) + doi.v[p] +
                      (c - spi.v[p]);
          *dst = a_(m,v,kb,jb,i);
        });
      }
      // band-ring copies for this hook: after the pack that reads the slabs, before
      // the scatter that writes the lent slab
      VetGdRingFlush();
      if (!vgd_hlp.done) {
      Kokkos::parallel_for("vet_vgd_ipc_k2", Kokkos::RangePolicy<>(vgd_cur, 0, 1),
      KOKKOS_LAMBDA(const int) {
        VgdIpcSetN(psig.data(), pslot.data(), np_, seq);
        VgdIpcWaitN(sig, np_, seq, err);
      });
      }
      vgd_hlp.on = true;
      vgd_hlp.ipc = true;
      vgd_hlp.seq = seq;
      vgd_hlp.xe = xe;
      vgd_hlp.ex = ex;
      vgd_hlp.px = px;
      vgd_hlp.i = i;
      vgd_hlp.gwa = gwa;
      vgd_hlp.a = a;
      vgd_ipc_n += 1.0;
      return true;
    }
#endif
    if (X.ns > 0) {
      Kokkos::parallel_for("vet_vgd_hl_pack", VgdLight(vgd_cur, 0, X.ns),
      KOKKOS_LAMBDA(const int c) {
        int ad = sl_(c);
        const int jb = ad % c2_;
        ad /= c2_;
        const int kb = ad % c3_;
        ad /= c3_;
        const int v = ad % nv_;
        const int m = ad/nv_;
        csb_(c) = a_(m,v,kb,jb,i);
      });
    }
    const std::vector<int> &sp = X.sp, &rp_ = X.rp;
    // ring copies (other shells' slabs; stream order pack, ring, zero, scatter) are
    // queued on an event behind the pack, so they run while the host posts the
    // messages (the send buffer is complete at the event: bitwise)
    VgdPackEvtRecord(vgd_cur);
    VetGdRingFlush();
    VgdPackEvtWait(vgd_cur);
    vgd_hlp.tq.reset();
    auto &tq = vgd_hlp.tq;
    auto &req = vgd_hlp.req;
    req.clear();
    for (int p = 0; p < np; ++p) {
      const int rk = vgd_prk[p];
      req.emplace_back();
      MPI_Irecv(crb_.data() + rp_[p], rp_[p+1] - rp_[p], MPI_ATHENA_REAL, rk, 7002,
                vgd_comm, &req.back());
      req.emplace_back();
      MPI_Isend(csb_.data() + sp[p], sp[p+1] - sp[p], MPI_ATHENA_REAL, rk, 7002,
                vgd_comm, &req.back());
    }
    vgd_tpost += tq.seconds();
    // vet_gd_band_ring: the sweep's slab copies (other shells only) while the messages
    // travel
    VetGdRingFlush();
    // zeros on Y's receive entries inside X's region (the dense unpack writes 0 there)
    const VgdHlEntry &Y = vgd_hle[ey];
    const bool zy = !xe && ((*gwa)[1 - px][i] != 0);
    if (zy && (Y.nr > 0 || Y.nw > 0)) {
      const int ny = Y.nr, nyw = Y.nw;
      auto ydl_ = Y.dl;
      auto ydw_ = Y.dw;
      Kokkos::parallel_for("vet_vgd_hl_zero", VgdLight(vgd_cur, 0, ny + nyw),
      KOKKOS_LAMBDA(const int c) {
        int ad = (c < ny) ? ydl_(c) : ydw_(c - ny);
        const int jb = ad % c2_;
        ad /= c2_;
        const int kb = ad % c3_;
        ad /= c3_;
        const int v = ad % nv_;
        const int m = ad/nv_;
        if (jb < w - ws || jb >= w + nx2 + ws || kb < w - ws || kb >= w + nx3 + ws) {
          return;
        }
        a_(m,v,kb,jb,i) = 0.0;
      });
    }

    vgd_hlp.on = true;
    vgd_hlp.ipc = false;
    vgd_hlp.xe = xe;
    vgd_hlp.ex = ex;
    vgd_hlp.px = px;
    vgd_hlp.i = i;
    vgd_hlp.gwa = gwa;
    vgd_hlp.a = a;
    return true;
  }
#endif
  return false;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdHlEnd
//! \brief second half of a list-path exchange started by VetGdHlBegin: MPI wait, size
//! check, scatter, last-writer ids

void RadiationVet::VetGdHlEnd() {
#if MPI_PARALLEL_ENABLED
  if (!vgd_hlp.on) {return;}
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int c2_ = nx2 + 2*w, c3_ = nx3 + 2*w, nv_ = vgd_n;
  const int i = vgd_hlp.i, px = vgd_hlp.px;
  auto a_ = vgd_hlp.a;
  DvceArray1D<Real> crb_ = vgd_crb;
  auto &req = vgd_hlp.req;
  auto &tq = vgd_hlp.tq;
  auto *gwa = vgd_hlp.gwa;
  const int np = static_cast<int>(vgd_prk.size());
  const VgdHlEntry &X = vgd_hlp.xe ? vgd_hxe[vgd_hlp.ex] : vgd_hle[vgd_hlp.ex];
  const std::vector<int> &sp = X.sp, &rp_ = X.rp;
  if (vgd_hlp.ipc) {
    // vet_gd_halo_ipc: the values are in my IPC buffer of this parity (k2 waited for
    // them); the scatter reads it in place of crb
    const long long par = vgd_hlp.seq & 1;
    crb_ = DvceArray1D<Real>(vgd_ipc_rbuf + par*vgd_ipc_half, vgd_ipc_cs);
    vgd_nexch += 1.0;
    for (int p = 0; p < np; ++p) {
      vgd_nbyte += static_cast<Real>(sp[p+1] - sp[p])*sizeof(Real);
    }
  } else {
  std::vector<MPI_Status> stat(req.size());
  MPI_Waitall(static_cast<int>(req.size()), req.data(), stat.data());
  vgd_tmpi += tq.seconds();
  vgd_nexch += 1.0;
  for (int p = 0; p < np; ++p) {
    int got = 0;
    MPI_Get_count(&stat[2*p], MPI_ATHENA_REAL, &got);
    if (got != rp_[p+1] - rp_[p]) {
      VgdFatal("vet_gd_halo_list: sender and receiver lists disagree");
    }
    vgd_nbyte += static_cast<Real>(sp[p+1] - sp[p])*sizeof(Real);
  }
  }   // !ipc
  // the scatter (dense unpack of the expanded buffer, on the written entries)
  const int nr = X.nr, nw = X.nw;
  auto dl_ = X.dl;
  auto sr_ = X.sr;
  auto dw_ = X.dw;
  auto sw_ = X.sw;
  auto w3_ = vgd_w3;
  if (nr + nw > 0 && !vgd_hlp.done) {
    Kokkos::parallel_for("vet_vgd_hl_scatter",
                         VgdLight(vgd_cur, 0, nr + nw),
    KOKKOS_LAMBDA(const int c) {
      int ad = (c < nr) ? dl_(c) : dw_(c - nr);
      const int jb = ad % c2_;
      ad /= c2_;
      const int kb = ad % c3_;
      ad /= c3_;
      const int v = ad % nv_;
      const int m = ad/nv_;
      if (c < nr) {
        // exact lists: sr = -1 is a read ghost the masks do not send: 0, as the dense
        // path
        const int s = sr_(c);
        a_(m,v,kb,jb,i) = (s >= 0) ? crb_(s) : 0.0;
        return;
      }
      const int cw = c - nr;
      Real val = 0.0;
      for (int u = 0; u < 3; ++u) {
        const int s = sw_(3*cw + u);
        val += w3_(m,kb,jb,v,u)*((s >= 0) ? crb_(s) : 0.0);
      }
      a_(m,v,kb,jb,i) = val;
    });
  }
  vgd_rg_lastx = vgd_hlp.xe;          // the band of this exchange came by list
  vgd_rg_lastx_e = vgd_hlp.ex;
  if (!vgd_hlp.xe) {
    (*gwa)[px][i] = X.gen;
  }
  vgd_hl_nuse += 1.0;
  vgd_hlp.on = false;
  vgd_hlp.ipc = false;
  vgd_hlp.done = false;
  vgd_hlp.a = VgdRag();
#endif
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::VetGdHaloCompact
//! \brief exchange of one shell: mask from VetGdHcPrep (prepared during the previous
//! shell's wait), compact pack, MPI, expansion; queue the next mask. false: the caller
//! unpacks the dense buffer. With a valid list entry X (vet_gd_halo_list_mb): gather
//! pack of X, MPI, zeros on the other pass's list Y in X's region, scatter of X (same
//! band values as the dense path; true). Defined below VetGdHlBegin / VetGdHlEnd.

template <class V>
bool RadiationVet::VetGdHaloCompact(V &a, const int nv, const int i0, const int ws) {
#if MPI_PARALLEL_ENABLED
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int mx = std::max(nx2, nx3)*ws*nv;
  const int i = i0;
  const bool inw = vgd_hinw;
  auto hl_ = vgd_hloc;
  auto a_ = a;
  auto so_ = vgd_soff;
  auto csb_ = vgd_csb;
  auto crb_ = vgd_crb;
  auto rb_ = vgd_rbuf;
  const int nvi = nv;
  const size_t rtot = static_cast<size_t>(vgd_prdsp[ws].empty() ? 0 :
                      (vgd_prdsp[ws].back() + vgd_prcnt[ws].back()))*nvi;
  const int np = static_cast<int>(vgd_prk.size());
  // vet_gd_halo_list_mb: list entries X (this pass) and Y (the other pass) of shell i,
  // and the last writers of the band per pass
  const int n1c = indcs.nx1 + 2*indcs.ng;
  const int c2 = nx2 + 2*w, c3 = nx3 + 2*w;
  bool lst = (nv == vgd_n) && (i >= 0) && (i < n1c) &&
             (static_cast<int64_t>(nmb1 + 1)*nv*c3*c2 <
              static_cast<int64_t>(std::numeric_limits<int>::max()));
  const int px = inw ? 0 : 1;
  const int ex = px*n1c + i;
  std::array<std::vector<int64_t>, 2> *gwa = nullptr;
  if (lst) {
    if (static_cast<int>(vgd_hle.size()) != 2*n1c) {vgd_hle.resize(2*n1c);}
    auto gwget = [&](const void *p) {
      auto &g = vgd_hl_gw[p];
      if (static_cast<int>(g[0].size()) != n1c) {
        g[0].assign(n1c, 0);
        g[1].assign(n1c, 0);
      }
      return &g;
    };
    if constexpr (std::is_same<V, VgdIView>::value) {
      gwa = gwget(static_cast<const void *>(a.d.data()));
    } else {
      lst = false;
    }
  }
  auto xvalid = [&]() {
    const VgdHlEntry &X = vgd_hle[ex];
    return X.gen > 0 && X.alpha == vgd_alpha && X.scut == vgd_scut && X.ws == ws;
  };
  // ---- list path (begin + end; Begin decides: the mask lists X, or with
  // vet_gd_halo_exact X's exact record) ----
  vgd_hlp.now = true;
  if (lst && VetGdHlBegin(a, nv, i0, ws)) {
    vgd_hlp.now = false;
    VetGdHlEnd();
    return true;
  }
  vgd_hlp.now = false;
  // ---- the dense path ----
  // the slot holding this shell's mask (prepared by the previous exchange), else now
  auto tagged = [&](const int sl) {
    return vgd_htag[sl][0] == i && vgd_htag[sl][1] == (inw ? 1 : 0) &&
           vgd_htag[sl][2] == ws && vgd_htag[sl][3] == vgd_scut &&
           vgd_htag[sl][4] == vgd_hsweep;
  };
  int slot = -1;
  if (tagged(0)) {slot = 0;}
  if (tagged(1)) {slot = 1;}
  if (slot < 0) {
    slot = 1 - vgd_hlast;
    VetGdHcPrep(slot, i, inw, ws);
  }
  vgd_hlast = slot;
  auto fs_ = vgd_hfs[slot];
  auto fr_ = vgd_hfr[slot];
  // (3) compact pack
  par_for("vet_vgd_hc_pack", vgd_cur, 0, nmb1, 0, 7, 0, mx - 1,
  KOKKOS_LAMBDA(const int m, const int o, const int t) {
    if (hl_(8*m + o) >= 0) {return;}
    const int oo = (o < 4) ? o : (o + 1);
    const int dk = oo/3 - 1, dj = oo%3 - 1;
    const int jn = (dj == 0) ? nx2 : ws;
    const int kn = (dk == 0) ? nx3 : ws;
    if (t >= nv*kn*jn) {return;}
    const size_t q = static_cast<size_t>(so_(ws,8*m + o))*nvi + t;
    if (fs_(q + 1) == fs_(q)) {return;}
    const int v = t/(kn*jn);
    const int r1 = t - v*kn*jn;
    const int kk = r1/jn;
    const int jj = r1 - kk*jn;
    const int js2 = (dj > 0) ? (w + nx2 - ws) : w;
    const int ks2 = (dk > 0) ? (w + nx3 - ws) : w;
    csb_(fs_(q)) = a_(m,v,ks2+kk,js2+jj,i);
  });
  vgd_cur.fence();
  std::vector<int> sp(np + 1), rp_(np + 1);
  {
    auto pbh = vgd_hpb[slot];
    for (int p = 0; p <= np; ++p) {
      sp[p] = pbh(p);
      rp_[p] = pbh(np + 1 + p);
    }
  }
  // vet_gd_halo_list_mb: X's lists from this exchange's masks (when X is not valid for
  // them), before the pipelined prep below reuses the masks' slot.
  // vet_gd_halo_exact: the exact lists when the record does not match (collective:
  // no rank has a usable record for this key)
  if (vgd_hx_on && lst) {
    if (static_cast<int>(vgd_hxe.size()) != 2*n1c) {vgd_hxe.resize(2*n1c);}
    vgd_hl_ndense += 1.0;
  }
  if (lst && !vgd_hx_on) {
    const VgdHlEntry &X0 = vgd_hle[ex];
    const bool nofit = X0.nofit && X0.alpha == vgd_alpha && X0.scut == vgd_scut &&
                       X0.ws == ws;
    if (!xvalid() && !nofit) {
      std::vector<int> pbv(sp);
      pbv.insert(pbv.end(), rp_.begin(), rp_.end());
      VetGdHlBuild(ex, i, inw, ws, fs_, fr_, pbv);
    }
    // the dense unpack writes X's whole region: the band follows X's lists when X is
    // valid, else it is unknown to the list path
    const int64_t g = xvalid() ? vgd_hle[ex].gen : -1;
    (*gwa)[px][i] = g;
    vgd_hl_ndense += 1.0;
  }
  Kokkos::Timer tq;
  std::vector<MPI_Request> req;
  for (int p = 0; p < np; ++p) {
    const int rk = vgd_prk[p];
    req.emplace_back();
    MPI_Irecv(crb_.data() + rp_[p], rp_[p+1] - rp_[p], MPI_ATHENA_REAL, rk, 7002,
              vgd_comm, &req.back());
    req.emplace_back();
    MPI_Isend(csb_.data() + sp[p], sp[p+1] - sp[p], MPI_ATHENA_REAL, rk, 7002,
              vgd_comm, &req.back());
  }
  vgd_tpost += tq.seconds();
  // the next shell's mask on the device while the messages travel
  if (vgd_hnext[0] >= 0) {
    VetGdHcPrep(1 - slot, vgd_hnext[0], vgd_hnext[1] != 0, vgd_hnext[2]);
  }
  std::vector<MPI_Status> stat(req.size());
  MPI_Waitall(static_cast<int>(req.size()), req.data(), stat.data());
  vgd_tmpi += tq.seconds();
  vgd_nexch += 1.0;
  for (int p = 0; p < np; ++p) {
    int got = 0;
    MPI_Get_count(&stat[2*p], MPI_ATHENA_REAL, &got);
    if (got != rp_[p+1] - rp_[p]) {
      VgdFatal("vet_gd_halo_compact: sender and receiver masks disagree");
    }
    vgd_nbyte += static_cast<Real>(sp[p+1] - sp[p])*sizeof(Real);
  }
  // (5) expansion to the dense receive layout
  Kokkos::parallel_for("vet_vgd_hc_expand", Kokkos::RangePolicy<>(vgd_cur, 0, rtot),
  KOKKOS_LAMBDA(const size_t q) {
    rb_(q) = (fr_(q + 1) > fr_(q)) ? crb_(fr_(q)) : 0.0;
  });
  return false;
#else
  return false;
#endif
}

// band exchanges called by VetGdSweep (its own TU)
template void RadiationVet::VetGdHalo<RadiationVet::VgdRag>(VgdRag &a, const int nv,
    const int i0, const int i1, const int ws, const bool mapd);

template bool RadiationVet::VetGdHlBegin<RadiationVet::VgdRag>(VgdRag &a, const int nv,
    const int i0, const int ws);

}  // namespace radvet
