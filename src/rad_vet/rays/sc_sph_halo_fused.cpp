//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_sph_halo_fused.cpp
//! \brief spherical-polar rays: the merged halo of two direction sets swept together.
//!
//!  called from: rad_vet_closure_ring.cpp, rays/sc_sph_sweep.cpp
//!  calls into: rays/sc_bandring.cpp, rays/sc_sph_halo_exact.cpp

#include <algorithm>
#include <cstdint>
#include <limits>
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
#include "rad_vet/rays/sc_halo.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::ScXffOk
//! \brief does the fused two-set sweep cover this configuration (same answer on every
//! rank)? Otherwise the two sweeps run

bool RadiationVet::ScXffOk() {
  if (!vgd_xff || !vgd_xf) {return false;}
  if (vray_icut < vgd_rag_cut) {return false;}   // (a cut below the shared slabs)
  if (vgd_csg) {
    // the cubed-sphere sweep (ScCsSweep): pass lists, no early send
    return !vgd_gband && (vgd_i_xb.dummy == vgd_i.dummy);
  }
  if (vgd_gband && !vgd_raw) {
    return false;
  }
  // IPC halo: one set per exchange (raw band: ScGbHaloX; CUDA lean: ScHlXff's
  // IPC form, no early send)
#if !defined(KOKKOS_ENABLE_CUDA)
  if (vgd_ipc_state == 1 && !vgd_raw) {return false;}
#endif
  if (vgd_i_xb.dummy != vgd_i.dummy) {return false;}
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::ScHlXff
//! \brief exchange of shell i (pass vgd_hinw, depth ws) for both fused sets in one
//! message per partner (piece p = [A_p, B_p]); bitwise the two exact-list exchanges.
//! false (nothing done) unless both sets hold exact records of this exchange.

bool RadiationVet::ScHlXff(const int i, const int ws) {
#if MPI_PARALLEL_ENABLED
  if (!vgd_hmpi || vgd_hlp.on) {return false;}
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int n1c = indcs.nx1 + 2*indcs.ng;
  const int c2_ = nx2 + 2*w, c3_ = nx3 + 2*w, nv_ = vgd_n;
  if (i < 0 || i >= n1c ||
      static_cast<int64_t>(nmb1 + 1)*nv_*c3_*c2_ >=
      static_cast<int64_t>(std::numeric_limits<int>::max())) {
    return false;
  }
  const int px = vgd_hinw ? 0 : 1;
  const int ex = px*n1c + i;
  if (ex >= static_cast<int>(vgd_hxe.size()) ||
      ex >= static_cast<int>(vgd_hxe_alt.size())) {
    return false;
  }
  const ScHlEntry &XA = vgd_hxe[ex];
  const ScHlEntry &XB = vgd_hxe_alt[ex];
  auto okx = [&](const ScHlEntry &X, const Real al) {
    return X.gen != -1 && X.alpha == al && ScHxCut(ex, X) && X.ws == ws && X.exact;
  };
  if (!okx(XA, vgd_alpha) || !okx(XB, vgd_xf_alt)) {return false;}
  const int np = static_cast<int>(vgd_prk.size());
  if (np > VGD_IPC_NS || static_cast<int>(XA.sp.size()) != np + 1 ||
      static_cast<int>(XB.sp.size()) != np + 1 ||
      static_cast<int>(XA.rp.size()) != np + 1 ||
      static_cast<int>(XB.rp.size()) != np + 1) {
    return false;
  }
  // merged layouts: piece p at so[p] (send) / ro[p] (receive), A first
  std::vector<int> so(np + 1), ro(np + 1);
  ScIpcIdx sa{}, sb{}, ra{}, rb{}, sm{}, rm{};
  for (int p = 0; p <= np; ++p) {
    so[p] = XA.sp[p] + XB.sp[p];
    ro[p] = XA.rp[p] + XB.rp[p];
    sa.v[p] = XA.sp[p];
    sb.v[p] = XB.sp[p];
    ra.v[p] = XA.rp[p];
    rb.v[p] = XB.rp[p];
    sm.v[p] = so[p];
    rm.v[p] = ro[p];
  }
  if (static_cast<int>(vgd_xff_sb.extent(0)) < std::max(so[np], 1)) {
    Kokkos::realloc(vgd_xff_sb, std::max(so[np], 1) + so[np]/8);
  }
  if (static_cast<int>(vgd_xff_rb.extent(0)) < std::max(ro[np], 1)) {
    Kokkos::realloc(vgd_xff_rb, std::max(ro[np], 1) + ro[np]/8);
  }
  auto sbuf_ = vgd_xff_sb;
  auto rbuf_ = vgd_xff_rb;
  auto aA_ = vgd_i;
  auto aB_ = vgd_i_xb;
  const int nsA = XA.ns, nsB = XB.ns;
  const int npk = np;
#if defined(KOKKOS_ENABLE_CUDA)
  // IPC form (CUDA lean, no early send): one kernel packs both sets' kept values of
  // partner p into its receive buffer, the last team sets the partners' flags, every
  // team waits for mine and scatters both sets (bitwise vet_vgd_xff_scatter)
  if (vgd_ipc_state == 1 &&
      vgd_rg_pend.empty() && static_cast<int>(XA.dsto.size()) == np &&
      static_cast<int>(XB.dsto.size()) == np) {
    if (!vgd_hlp.cnt.is_allocated()) {
      vgd_hlp.cnt = Kokkos::View<int, DevMemSpace>("vet_vgd_ipc_cnt");
    }
    const long long seq = ++vgd_ipc_seq, par = seq & 1;
    auto cnt_ = vgd_hlp.cnt;
    long long *sig = vgd_ipc_sig;
    int *err = vgd_ipc_err;
    auto psig = vgd_ipc_psig.d_view;
    auto pslot = vgd_ipc_pslot.d_view;
    auto prb = vgd_ipc_prb.d_view;
    auto phalf = vgd_ipc_phalf.d_view;
    auto pcs = vgd_ipc_pcs.d_view;
    const Real *rb_ = vgd_ipc_rbuf + par*vgd_ipc_half;
    const size_t cs_ = vgd_ipc_cs;
    ScIpcIdx da{}, db{};
    for (int p = 0; p < np; ++p) {da.v[p] = XA.dsto[p]; db.v[p] = XB.dsto[p];}
    auto aA_ = vgd_i;
    auto aB_ = vgd_i_xb;
    auto slA_ = XA.sl;
    auto slB_ = XB.sl;
    const int nsA = XA.ns, nsB = XB.ns;
    const int nrA = XA.nr, nA = XA.nr + XA.nw, nrB = XB.nr, nB = XB.nr + XB.nw;
    auto dlA_ = XA.dl;
    auto srA_ = XA.sr;
    auto dwA_ = XA.dw;
    auto swA_ = XA.sw;
    auto dlB_ = XB.dl;
    auto srB_ = XB.sr;
    auto dwB_ = XB.dw;
    auto swB_ = XB.sw;
    auto w3A_ = vgd_w3;
    auto w3B_ = vgd_tabd.w3;
    constexpr int kT = 128;
    const int nt = std::max(nsA + nsB, nA + nB);
    const int nl = std::min((nt + kT - 1)/kT, 216);
    Kokkos::parallel_for("vet_vgd_ipc_xff", Kokkos::Experimental::require(
        Kokkos::TeamPolicy<>(vgd_cur, std::max(nl, 1), kT),
        Kokkos::Experimental::WorkItemProperty::HintLightWeight),
    KOKKOS_LAMBDA(const TeamMember_t &tm) {
      const int e0 = tm.league_rank()*kT + tm.team_rank(), de = tm.league_size()*kT;
      for (int c = e0; c < nsA + nsB; c += de) {
        const bool b = (c >= nsA);
        const int cc = b ? (c - nsA) : c;
        const ScIpcIdx &sx = b ? sb : sa;
        int p = 0;
        while (p + 1 < npk && sx.v[p+1] <= cc) {++p;}
        int ad = b ? slB_(cc) : slA_(cc);
        const int jb = ad % c2_;
        ad /= c2_;
        const int kb = ad % c3_;
        ad /= c3_;
        const int v = ad % nv_;
        const int m = ad/nv_;
        Real *dst = reinterpret_cast<Real *>(prb(p)) + par*phalf(p) + (b ? pcs(p) : 0) +
                    (b ? db.v[p] : da.v[p]) + (cc - sx.v[p]);
        *dst = b ? aB_(m,v,kb,jb,i) : aA_(m,v,kb,jb,i);
      }
#if defined(__CUDA_ARCH__)
      __threadfence_system();
#endif
      tm.team_barrier();
      if (tm.team_rank() == 0) {
        if (Kokkos::atomic_fetch_add(&cnt_(), 1) == tm.league_size() - 1) {
          cnt_() = 0;
          sig[VGD_IPC_SEQ] = seq;   // the exchange number (graphs read it)
          ScIpcSetN(psig.data(), pslot.data(), npk, seq);
        }
        ScIpcWaitN(sig, npk, seq, err);
      }
      tm.team_barrier();
      for (int c0 = e0; c0 < nA + nB; c0 += de) {
        const bool b = (c0 >= nA);
        const int c = b ? (c0 - nA) : c0;
        const int nr = b ? nrB : nrA;
        const Real *rx = rb_ + (b ? cs_ : 0);
        int ad = (c < nr) ? (b ? dlB_(c) : dlA_(c)) : (b ? dwB_(c - nr) : dwA_(c - nr));
        const int jb = ad % c2_;
        ad /= c2_;
        const int kb = ad % c3_;
        ad /= c3_;
        const int v = ad % nv_;
        const int m = ad/nv_;
        Real val = 0.0;
        if (c < nr) {
          const int s = b ? srB_(c) : srA_(c);
          val = (s >= 0) ? rx[s] : 0.0;
        } else {
          const int cw = c - nr;
          for (int u = 0; u < 3; ++u) {
            const int s = b ? swB_(3*cw + u) : swA_(3*cw + u);
            const Real wq = b ? w3B_(m,kb,jb,v,u) : w3A_(m,kb,jb,v,u);
            val += wq*((s >= 0) ? rx[s] : 0.0);
          }
        }
        if (b) {
          aB_(m,v,kb,jb,i) = val;
        } else {
          aA_(m,v,kb,jb,i) = val;
        }
      }
    });
    vgd_ipc_n += 1.0;
    vgd_nexch += 2.0;
    vgd_hl_nuse += 2.0;
    vgd_xff_nmsg += 1.0;
    vgd_rg_lastx = true;
    vgd_rg_lastx_e = ex;
    return true;
  }
#endif
  if (nsA + nsB > 0) {
    auto slA_ = XA.sl;
    auto slB_ = XB.sl;
    Kokkos::parallel_for("vet_vgd_xff_pack", ScLight(vgd_cur, 0, nsA + nsB),
    KOKKOS_LAMBDA(const int c) {
      const bool b = (c >= nsA);
      const int cc = b ? (c - nsA) : c;
      const ScIpcIdx &sx = b ? sb : sa;
      int p = 0;
      while (p + 1 < npk && sx.v[p+1] <= cc) {++p;}
      int ad = b ? slB_(cc) : slA_(cc);
      const int jb = ad % c2_;
      ad /= c2_;
      const int kb = ad % c3_;
      ad /= c3_;
      const int v = ad % nv_;
      const int m = ad/nv_;
      const int pos = sm.v[p] + (b ? (sa.v[p+1] - sa.v[p]) : 0) + (cc - sx.v[p]);
      sbuf_(pos) = b ? aB_(m,v,kb,jb,i) : aA_(m,v,kb,jb,i);
    });
  }
  // band-ring copies of other shells go between pack and scatter (as ScHlBegin)
  ScPackEvtRecord(vgd_cur);
  ScRingFlush();
  ScPackEvtWait(vgd_cur);
  Kokkos::Timer tq;
  std::vector<MPI_Request> req;
  for (int p = 0; p < np; ++p) {
    const int rk = vgd_prk[p];
    req.emplace_back();
    MPI_Irecv(rbuf_.data() + ro[p], ro[p+1] - ro[p], MPI_ATHENA_REAL, rk, 7005, vgd_comm,
              &req.back());
    req.emplace_back();
    MPI_Isend(sbuf_.data() + so[p], so[p+1] - so[p], MPI_ATHENA_REAL, rk, 7005, vgd_comm,
              &req.back());
  }
  vgd_tpost += tq.seconds();
  std::vector<MPI_Status> stat(req.size());
  MPI_Waitall(static_cast<int>(req.size()), req.data(), stat.data());
  vgd_tmpi += tq.seconds();
  vgd_nexch += 2.0;
  for (int p = 0; p < np; ++p) {
    int got = 0;
    MPI_Get_count(&stat[2*p], MPI_ATHENA_REAL, &got);
    if (got != ro[p+1] - ro[p]) {
      ScFatal("vet_gd_rotate_xfade_fuse: sender and receiver lists disagree");
    }
    vgd_nbyte += static_cast<Real>(so[p+1] - so[p])*sizeof(Real);
  }
  // scatter both sets (A, then B), as ScHlEnd per set
  const int nrA = XA.nr, nwA = XA.nw, nrB = XB.nr, nwB = XB.nw;
  const int nA = nrA + nwA, nB = nrB + nwB;
  if (nA + nB > 0) {
    auto dlA_ = XA.dl;
    auto srA_ = XA.sr;
    auto dwA_ = XA.dw;
    auto swA_ = XA.sw;
    auto dlB_ = XB.dl;
    auto srB_ = XB.sr;
    auto dwB_ = XB.dw;
    auto swB_ = XB.sw;
    auto w3A_ = vgd_w3;
    auto w3B_ = vgd_tabd.w3;
    Kokkos::parallel_for("vet_vgd_xff_scatter", ScLight(vgd_cur, 0, nA + nB),
    KOKKOS_LAMBDA(const int c0) {
      const bool b = (c0 >= nA);
      const int c = b ? (c0 - nA) : c0;
      const int nr = b ? nrB : nrA;
      const ScIpcIdx &rx = b ? rb : ra;
      // merged receive position of compact source s of this set
      auto at = [&](const int s) -> int {
        int p = 0;
        while (p + 1 < npk && rx.v[p+1] <= s) {++p;}
        return rm.v[p] + (b ? (ra.v[p+1] - ra.v[p]) : 0) + (s - rx.v[p]);
      };
      int ad = (c < nr) ? (b ? dlB_(c) : dlA_(c)) : (b ? dwB_(c - nr) : dwA_(c - nr));
      const int jb = ad % c2_;
      ad /= c2_;
      const int kb = ad % c3_;
      ad /= c3_;
      const int v = ad % nv_;
      const int m = ad/nv_;
      if (c < nr) {
        const int s = b ? srB_(c) : srA_(c);
        const Real x = (s >= 0) ? rbuf_(at(s)) : 0.0;
        if (b) {
          aB_(m,v,kb,jb,i) = x;
        } else {
          aA_(m,v,kb,jb,i) = x;
        }
        return;
      }
      const int cw = c - nr;
      Real val = 0.0;
      for (int u = 0; u < 3; ++u) {
        const int s = b ? swB_(3*cw + u) : swA_(3*cw + u);
        const Real wq = b ? w3B_(m,kb,jb,v,u) : w3A_(m,kb,jb,v,u);
        val += wq*((s >= 0) ? rbuf_(at(s)) : 0.0);
      }
      if (b) {
        aB_(m,v,kb,jb,i) = val;
      } else {
        aA_(m,v,kb,jb,i) = val;
      }
    });
  }
  vgd_hl_nuse += 2.0;
  vgd_xff_nmsg += 1.0;
  vgd_rg_lastx = true;   // the band of this exchange came from both sets' lists
  vgd_rg_lastx_e = ex;
  return true;
#else
  (void) i;
  (void) ws;
  return false;
#endif
}

}  // namespace radvet
