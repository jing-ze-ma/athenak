//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_closure_ring.cpp
//! \brief ring: the closure of rotating direction sets, the ring of offsets of
//! the older sets and the mixture.
//!
//!  called from: rad_vet_closure.cpp, rad_vet_restart.cpp, rad_vet_tasks.cpp,
//!    driver/dt_retry.cpp
//!  calls into: rad_vet_bvals.cpp, rad_vet_closure_faces.cpp, rad_vet_closure_rays.cpp,
//!    rad_vet_init.cpp, rays/sc_build.cpp, rays/sc_cart_setup.cpp,
//!    rays/sc_cart_sweep.cpp, rays/sc_directions.cpp, rays/sc_sph_halo_fused.cpp,
//!    rays/sc_sph_sweep.cpp, rays/sc_sph_walls.cpp, rays/sc_tables.cpp

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_args.hpp"
#include "rad_vet/rad_vet_moments.hpp"

namespace radvet {

extern bool g_vra_adv;   // vet_rotate_advance (rays/sc_directions.cpp)

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RingUnion
//! \brief ring, cubed-sphere sc rays: the direction union (cln_fvd) of sets
//! b, b+1 (slots blk mod 2) with cross-fade weights 1 - s, s

void RadiationVet::RingUnion() {
  const int n = vgd_n;
  const int64_t b = cln_rbld;
  if (cln_fvd.extent_int(0) != 2*n) cln_fvd = DvceArray2D<Real>("cln_fvd", 2*n, 4);
  auto d_h = Kokkos::create_mirror_view(cln_fvd);
  for (int64_t blk = b; blk <= b + 1; ++blk) {
    const int q = static_cast<int>(blk & 1);
    const std::vector<double> d = ScRotate(RotAngleBlk(blk));
    const Real w = (blk > b) ? cln_rs : 1.0 - cln_rs;
    for (int e = 0; e < n; ++e) {
      for (int c = 0; c < 3; ++c) d_h(q*n+e,c) = d[4*e+c];
      d_h(q*n+e,3) = w*d[4*e+3];
    }
  }
  Kokkos::deep_copy(cln_fvd, d_h);
  cln_fvn = 2*n;
}

//----------------------------------------------------------------------------------------
//! \fn bool RingMix
//! \brief ring: raw moments q of the ring mix: current set (slot 0), next (qn,
//! weight s), older sets r (weights 1, the last 1 - s), J-weighted with r_r =
//! J_{b-r}/J_b.  Invalid offsets (r_r < 0) skipped; false if the current set has no rays

struct RingL {
  DvceArray5D<Real> rl, cs;
  Kokkos::Array<int, CLN_N> fm;
  int nf, ol, ks, js, is;
  // (slot o's normalised moment n at interior cell (k,j,i); slot ol is read from cln_s)
  KOKKOS_INLINE_FUNCTION
  Real operator()(const int m, const int o, const int n, const int k, const int j,
                  const int i) const {
    if (o == ol) {
      const Real jr = cs(m,CLN_J,k+ks,j+js,i+is), ij = (jr > 0.0) ? 1.0/jr : 0.0;
      return (CellMoment(n) ? 1.0 : ij)*cs(m,n,k+ks,j+js,i+is);
    }
    return rl(m,((ol < 0) ? o*nf : 0)+fm[n],k,j,i);
  }
};

KOKKOS_INLINE_FUNCTION
bool RingMix(const RingL &rl, const DvceArray5D<Real> &rr,
             const Kokkos::Array<int, CLN_N> &fr, const int nfr, const int qn,
             const int nr, const Real s, const int m, const int k, const int j,
             const int i, Real *q) {
  const Real jb = rl(m,0,CLN_J,k,j,i);
  if (!(jb > 0.0)) return false;
  Real a[RadiationVet::VET_RA_NF];
  const Real wn = s*rl(m,qn,CLN_J,k,j,i)/jb;
  Real na = 1.0 + s, nw = 1.0 + wn;
  for (int r = 0; r < nr; ++r) {
    const Real rj = rr(m,r*nfr+fr[int{CLN_J}],k,j,i), w = (r == nr - 1) ? 1.0 - s : 1.0;
    a[r] = (rj >= 0.0) ? w*rj : 0.0;
    na += (rj >= 0.0) ? w : 0.0;
    nw += a[r];
  }
  q[CLN_J] = jb*nw/na;
  for (int n = 0; n < CLN_N; ++n) {
    if (n == CLN_J) continue;
    if (rl.fm[n] < 0) {q[n] = 0.0; continue;}
    const Real gb = rl(m,0,n,k,j,i);
    if (fr[n] < 0) {   // (not on the ring: the swept sets', as the fallback below)
      q[n] = (1.0 - s)*(gb*jb) + s*(rl(m,qn,n,k,j,i)*rl(m,qn,CLN_J,k,j,i));
      continue;
    }
    Real v = gb + wn*rl(m,qn,n,k,j,i);
    for (int r = 0; r < nr; ++r) v += a[r]*(gb + rr(m,r*nfr+fr[n],k,j,i));
    q[n] = v/nw*q[CLN_J];
  }
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool RingEma
//! \brief ring with vet_rotate_avg_ema (production's form, RotAvgEma*,
//! RotXfBlend): the raw moments q of (1 - s) A + s T in the normalised fields and J. A =
//! the current set with the EMA offsets (e = rr's first nfr, its J entry the ratio r =
//! J_avg/J_b), T = the next set with the EMA the coming boundary makes of this state's
//! Delta (RotAvgEmaShift): aw = max(1/N, a_ch, 1/(n + 1)), r' = (1 - aw) r/R + aw, e' =
//! (1 - aw)(r/R)/r' (e - D), R = J_{b+1}/J_b; a_ch is this build's (pd's entry nfr), n
//! rr's entry nfr + 2 (0: no history, A and T the sets themselves).  The face sums (not
//! cell moments) as (1 - s) A + s T in raw moments, linear in the faces' own sums: no
//! cell's J ratio weights its low face, so that a face and its mirror face mix alike

KOKKOS_INLINE_FUNCTION
bool RingEma(const RingL &rl, const DvceArray5D<Real> &rr,
             const DvceArray5D<Real> &pd, const Kokkos::Array<int, CLN_N> &fr,
             const int nfr, const int qn, const Real s, const Real a0, const int m,
             const int k, const int j, const int i, Real *q) {
  const Real jb = rl(m,0,CLN_J,k,j,i), jn = rl(m,qn,CLN_J,k,j,i);
  if (!(jb > 0.0)) return false;
  const Real nset = rr(m,nfr+2,k,j,i), r = rr(m,fr[int{CLN_J}],k,j,i);
  const bool hist = nset > 0.5 && r > 0.0, sh = hist && jn > 0.0;
  Real jt = jn, fac = 0.0;
  if (sh) {
    const Real aw = (a0 < 0.0) ? fmax(0.0, fmin(pd(m,nfr,k,j,i), 1.0)) :
                     fmax(fmax(a0, fmin(pd(m,nfr,k,j,i), 1.0)), 1.0/(nset + 1.0));
    const Real rq = r*jb/jn, r1 = (1.0 - aw)*rq + aw;
    fac = (1.0 - aw)*rq/r1;
    jt = jn*r1;
  }
  q[CLN_J] = (1.0 - s)*jb*(hist ? r : 1.0) + s*jt;
  for (int n = 0; n < CLN_N; ++n) {
    if (n == CLN_J) continue;
    if (rl.fm[n] < 0) {q[n] = 0.0; continue;}
    const Real gb = rl(m,0,n,k,j,i), gn = rl(m,qn,n,k,j,i);
    if (fr[n] < 0) {   // (not on the ring: the swept sets', as RingMix)
      q[n] = (1.0 - s)*(gb*jb) + s*(gn*jn);
      continue;
    }
    const Real e = rr(m,fr[n],k,j,i);
    const Real ga = gb + (hist ? e : 0.0), gt = sh ? gn + fac*(e - (gn - gb)) : gn;
    q[n] = CellMoment(n) ? ((1.0 - s)*ga + s*gt)*q[CLN_J]
                              : (1.0 - s)*ga*jb*(hist ? r : 1.0) + s*gt*jt;
  }
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationVet::RingFields
//! \brief ring: the fields the ring carries, fm[n] = slot or -1; returns their
//! number.  Left out: fields cln_s lacks, the stress parts of the face sums and (cubed
//! sphere) CLN_G (mixed as 0).  file: the restart layout; files before version 7 keep the
//! stress parts (and CLN_G on the cubed sphere)

int RadiationVet::RingFields(int *fm, const bool file) const {
  const bool cub = pmy_pack->pmesh->use_cubed_sphere;
  const int ns = NumClosureSlots(cub, cln_rfh);   // (the slots of cln_s)
  int nf = 0;
  for (int n = 0; n < CLN_N; ++n) {
    const bool fs = (n >= CLN_F2 && n < CLN_G && (n - CLN_F2)%4 != 0) ||   // (stress
                    (n >= CLN_FH && n < CLN_F1 && (n - CLN_FH)%4 != 0);    // parts)
    const bool un = fs || (cub && (n == CLN_G || n == CLN_G + 1));
    const bool f6 = file && t2_rstv < 7;
    fm[n] = (n < ns && !(f6 ? fs && !cub : un)) ? nf++ : -1;
  }
  return nf;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RingShift
//! \brief ring: block b's first build shifts the ring by the pending Delta:
//! d_1 = -Delta, d_j -= Delta, r_j' = r_{j-1}/R, r_1' = 1/R

void RadiationVet::RingShift(const int64_t b) {
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int nr = vra_k - 1;
  Kokkos::Array<int, CLN_N> fm;
  const int nf = RingFields(fm.data());
  const Real tk = static_cast<Real>(b)*vgd_rot_dt;
  if (std::abs(cln_rpt - tk) > 1.0e-9*vgd_rot_dt) cln_roff_sum += 1.0;
  vra_gen += 1;
  auto rr = cln_rr;
  auto pd = cln_rpd;
  if (vra_ema) {
    // (EMA update with the boundary's Delta and weight a_ch; no history or no Delta:
    // the current set alone)
    const Real a0 = vra_ema_wb ? -1.0 : 1.0/vra_ema_n;   // (< 0 = a_ch alone)
    par_for("cln_ring_eshift", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1, 0,
            id.nx3 - 1, 0, id.nx2 - 1, 0, id.nx1 - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int fj = fm[int{CLN_J}];
      const Real R = pd(m,fj,k,j,i), nset = rr(m,nf+2,k,j,i), r = rr(m,fj,k,j,i);
      rr(m,nf,k,j,i) = pd(m,nf+1,k,j,i);
      rr(m,nf+1,k,j,i) = pd(m,nf+2,k,j,i);
      if (!(R > 0.0) || !(nset > 0.5) || !(r > 0.0)) {
        for (int f = 0; f < nf; ++f) rr(m,f,k,j,i) = (f == fj) ? 1.0 : 0.0;
        rr(m,nf+2,k,j,i) = (R > 0.0) ? 1.0 : 0.0;
        return;
      }
      const Real aw = (a0 < 0.0) ? fmax(0.0, fmin(pd(m,nf,k,j,i), 1.0)) :
                      fmax(fmax(a0, fmin(pd(m,nf,k,j,i), 1.0)), 1.0/(nset + 1.0));
      const Real r1 = (1.0 - aw)*r/R + aw, fac = (1.0 - aw)*(r/R)/r1;
      for (int f = 0; f < nf; ++f) {
        rr(m,f,k,j,i) = (f == fj) ? r1 : fac*(rr(m,f,k,j,i) - pd(m,f,k,j,i));
      }
      rr(m,nf+2,k,j,i) = (aw >= 1.0) ? 1.0 : fmin(nset + 1.0, 1.0e9);
    });
  } else {
  par_for("cln_ring_shift", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1, 0, id.nx3 - 1,
          0, id.nx2 - 1, 0, id.nx1 - 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int fj = fm[int{CLN_J}];
    const Real rp = pd(m,fj,k,j,i);
    for (int f = 0; f < nf; ++f) {
      for (int r = nr - 1; r > 0; --r) {
        const Real o = rr(m,(r-1)*nf+f,k,j,i);
        rr(m,r*nf+f,k,j,i) = (f != fj) ? o - pd(m,f,k,j,i) :
                             ((rp > 0.0 && o >= 0.0) ? o/rp : -1.0);
      }
      rr(m,f,k,j,i) = (f != fj) ? -pd(m,f,k,j,i) : ((rp > 0.0) ? 1.0/rp : -1.0);
    }
  });
  }
  cln_rok = true;
  cln_rb = b;
}

int64_t RadiationVet::RingBlk(const Real dt) const {
  const Real tau = (pmy_pack->pmesh->time + dt)/vgd_rot_dt;
  const int64_t b = static_cast<int64_t>(std::floor(tau));
  return (tau - static_cast<Real>(b) < 1.0e-9) ? b - 1 : b;   // (RingSets)
}

// (live with the cross-fade: every build sweeps the sets b and b+1)
bool RadiationVet::RingTwo() const {
  return cln_ring == 1 && vgd_xf && !(vra_on && !g_vra_adv);
}

// (cubed-sphere sc rays: true if the face half ranges cln_fh are rebuilt before every
// solve; else they are state, kept by the dt_retry checkpoint and restart files)
bool RadiationVet::FhBuilt() const {
  return cln_rfh || (cln_cssc && cln_extrap == 0 && !cln_closit && !cln_edd &&
                     !cln_sphere && (!cln_ring || RingTwo()));
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RingPreShift
//! \brief before the dt_retry checkpoint: with ring = live, cross-fade and
//! vet_gd_rotate_align every retry's build lies in one block, so a step starting on a
//! block boundary shifts the ring here and the checkpoint skips offsets and Delta
//! (cln_rfix).  RingBlk: the block of a build at time + dt

void RadiationVet::RingPreShift() {
  const Real dt = pmy_pack->pmesh->dt;
  cln_rfix = RingTwo() && vgd_rot_dt > 0.0 && dt < vgd_rot_dt &&
             RingBlk(dt) == RingBlk(1.0e-6*dt);   // (any retry's dt)
  if (cln_rfix && cln_rok && cln_rpb == cln_rb && RingBlk(dt) == cln_rb + 1) {
    RingShift(cln_rb + 1);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RingMoments
//! \brief ring: closure of the rotating union of K = vet_rotate_avg sets with
//! cross-fade, sweeping only the current set b and the next b+1:
//!   M = (1 - s)/K M_{b-K+1} + 1/K sum_{j=1}^{K-2} M_{b-j} + 1/K M_b + s/K M_{b+1}
//! older sets from the current g = M/J plus offsets d_j = g_{b-j} - g_b,
//! r_j = J_{b-j}/J_b (RingMix).  The build closing block b (s = 1) gives Delta =
//! M_{b+1} - M_b; block b+1's first build shifts the ring (RingShift).  Offsets are
//! piecewise constant, so the closure stays continuous in t.  A fresh start fills the
//! ring with one full sweep.  Unrealisable mixtures fall back to (1 - s) M_b + s M_{b+1}
//! (counted); linear in the moments, so face sums stay single-valued and conservative

void RadiationVet::RingMoments() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack, nk = RingK(), nr = nk - 1, ns = cln_rns;
  const int64_t b = cln_rbld;
  const Real s = cln_rs;
  const int c3 = ke - ks + 1, c2 = je - js + 1, c1 = ie - is + 1;
  // (swept moments on interior cells: current set in cln_rl slot 0, next in slot 1; a
  // fill's older sets straight into cln_rr.  fr: fields of offsets and Delta)
  Kokkos::Array<int, CLN_N> fm, fr;
  const int nf = RingFields(fm.data()), nfr = RingFields(fr.data());
  // (sweep order: cubed-sphere sc rays the older sets by block, then b, b+1 by slot
  // blk mod 2, as FaceHalf folds them; elsewhere slot order (keeps the restart
  // bitwise).  The last swept set stays in cln_s, ol its slot (-1: rl holds both))
  std::vector<int> ord(ns);
  for (int q = 0; q < ns; ++q) ord[q] = q;
  auto okey = [&](int u) {
    const int64_t bu = cln_fvblk[u];
    return (bu >= b) ? std::make_pair(int64_t{1}, bu & 1)
                     : std::make_pair(int64_t{0}, bu);
  };
  if (cln_cssc) {
    std::sort(ord.begin(), ord.end(), [&](int u, int v) {return okey(u) < okey(v);});
  }
  const int64_t pl = cln_fvblk[ord[ns - 1]] - (b - nk + 1);
  const int ol = (cln_ring != 1) ? -1 : ((pl == nr) ? 0 : ((pl == nk) ? 1 : -1));
  // (over the extrap history's free set, cln_mc while it holds none, or cln_csy off the
  // arena, where they fit: rewritten before every read, dead across a build)
  const int nrl = ((ol < 0) ? 2 : 1)*nf;
  const std::size_t nrs = static_cast<std::size_t>(nmb)*nrl*c3*c2*c1;
  if (ExtrapSx() && !cln_have_c && cln_mc.span() >= nrs) {
    cln_rl = DvceArray5D<Real>(cln_mc.data(), nmb, nrl, c3, c2, c1);
  } else if (cln_csy.span() >= nrs && !MemArenaIn(cln_csy.data())) {
    cln_rl = DvceArray5D<Real>(cln_csy.data(), nmb, nrl, c3, c2, c1);
  } else if (cln_rl.extent_int(0) != nmb || cln_rl.extent_int(1) != nrl) {
    Kokkos::realloc(cln_rl, nmb, nrl, c3, c2, c1);
  }
  if (cln_rr.extent_int(0) != nmb) {   // (a restart has read them)
    vra_gen += 1;   // (cln_rr's writers: the dt_retry checkpoint skips a clean copy)
    Kokkos::realloc(cln_rr, nmb, RingNr(nfr), c3, c2, c1);
    Kokkos::realloc(cln_rpd, nmb, RingNd(nfr), c3, c2, c1);
    cln_rok = false;
  }
  // vet_rotate_avg_ema: rr = EMA offsets, aux J and D11 (last boundary), set count n;
  // pd also a_ch and the next set's J and D11
  const bool ema = vra_ema, estep = vra_ema_chstep;
  const Real a0 = ema ? (vra_ema_wb ? -1.0 : 1.0/vra_ema_n) : 1.0, eps = kEmaEps;
  std::vector<int> slot(nk + 2, -1);   // slot of block b - nk + 1 + p
  for (int q = 0; q < ns; ++q) {
    const int64_t p = cln_fvblk[q] - (b - nk + 1);
    if (p >= 0 && p <= nk) slot[p] = q;
  }
  bool full = true;
  for (int p = 0; p < nr; ++p) full = full && slot[p] >= 0;
  const bool shift = cln_rok && b == cln_rb + 1 && cln_rpb == cln_rb;
  const bool fill = (!cln_rok || b != cln_rb) && !shift && full;
  // (the pending Delta is read again only after this build has rewritten it)
  cln_rpfree = cln_ring == 1 && slot[nk] >= 0 && !shift;
  // moments of each swept set (ray_min tested on the mixture); cubed-sphere sc rays: the
  // union's directions first
  if (cln_cssc) RingUnion();
  auto dfull = cln_fvd;
  bool fh0 = true;
  const Real ray_min = cln_raymin;
  cln_raymin = 0.0;
  // (J ratios only where both sets' J >= closure_ray_min E, so a source
  // hole's J ~ 1e-200 gives no huge ratio)
  const Real jmin = ray_min;
  auto iwr = iw;
  auto cs = cln_s;
  const RingL rl{cln_rl, cln_s, fm, nf, ol, ks, js, is};
  auto rr = cln_rr;
  // ring_fuse: sweep the live sets b, b+1 in one fused ScSweep (vgd_xff_now;
  // own tables and intensities each); the loop below then takes each set's moments
  const bool fuse = vgd_on && cln_rfuse && !cln_scmem && ns == 2 && !fill && ScXffOk();
  // (the cut below the sets' shared slabs, ScRagThick: the same cut on every rank)
  const bool low = vgd_rag_cut > vray_icut;
  if (vgd_on && cln_rfuse && !fuse && ns == 2 && !fill) {
    vgd_xff_nfb += 1.0;
    if (low && !vgd_rag_said && global_variable::my_rank == 0) {
      std::cout << "<rad_vet> ring_fuse: the ray cut (shell " << vray_icut << ") is "
                << "below the shells whose slabs the two sets share (" << vgd_rag_cut
                << "): two sweeps per build while it does; a restart shares fewer shells"
                << std::endl;
    }
    vgd_rag_said = vgd_rag_said || low;
  }
  if (vgd_on && cln_rfuse && !fuse && !low && ns == 2 && !fill && !cln_rfuse_said &&
      global_variable::my_rank == 0) {
    std::cout << "<rad_vet> ring_fuse: not taken (needs the sc rays, no "
              << "reflecting top and a fused-sweep configuration, ScXffOk): two "
              << "sweeps per build" << std::endl;
    cln_rfuse_said = true;
  }
  if (fuse) {
    const int qa = (cln_fvblk[0] < cln_fvblk[1]) ? 0 : 1, qb = 1 - qa;
    vgd_scut = vray_icut;
    ScXfTables(cln_fvkey[qb]);
    vgd_xff_prep = true;
    ScSweep();
    ScXfTables(cln_fvkey[qa]);
    ScSweep();
    vgd_xff_prep = false;
    vgd_xff_now = true;
    ScSweep();
    vgd_xff_now = false;
    vgd_xff_n += 1.0;   // (fused build count)
  }
  for (int qo = 0; qo < ns; ++qo) {
    const int q = ord[qo];
    const int64_t p = cln_fvblk[q] - (b - nk + 1);
    const bool older = (p >= 0 && p < nr);
    const bool keep = (older && fill) || p == nr || p == nk;
    if (older && fill) vra_gen += 1;
    if (!vgd_on && vgd_slab != 1) {   // (vet_sc: the set of block cln_fvblk[q] alone)
      VetRotSet(VetRotBeta(cln_fvblk[q]));
      ScCartSweep();
      ClosureFromRays();
    } else {
      cln_fvd = cln_fvdm[q];
      // the sc rays: set q swept alone (its tables, its own weights)
      cln_rq = q;
      cln_rqm = static_cast<int>(cln_fvblk[q] & 1);
      ScXfTables(cln_fvkey[q]);
      vgd_scut = vray_icut;
      if (!fuse) ScSweep();
      if (cln_scmem) ScTopCopy();
      ClosureFromSc();
      if (cln_consist || cln_cssc) RayFaceSums();
      if (cln_cssc && !cln_rfh && (p == nr || p == nk)) {
        FaceHalf(dfull, cln_rqm*vgd_n, fh0, 0);
        fh0 = false;
      }
    }
    const int sl = (p == nk) ? 1 : 0;   // (the slot of a live set)
    if (!keep || (!older && sl == ol)) continue;
    const int o = static_cast<int>(older ? (nr - 1 - p)*nfr : ((ol < 0) ? sl*nf : 0));
    auto ra = older ? rr : cln_rl;
    const auto fa = older ? fr : fm;
    par_for("cln_ring_raw", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real jr = cs(m,CLN_J,k,j,i), ij = (jr > 0.0) ? 1.0/jr : 0.0;
      for (int n = 0; n < CLN_N; ++n) {
        if (fa[n] < 0) continue;
        ra(m,o+fa[n],k-ks,j-js,i-is) = (CellMoment(n) ? 1.0 : ij)*cs(m,n,k,j,i);
      }
    });
  }
  if (vgd_on || vgd_slab == 1) ScPost();
  cln_rpfree = false;
  cln_fvd = dfull;
  cln_raymin = ray_min;
  // (jump: the next set only at the build that closes the block; its weight is s = 0)
  const int qb = 0, qn = (slot[nk] >= 0) ? 1 : 0;
  auto pd = cln_rpd;
  if (!cln_rok || b != cln_rb) {
    // (f0: fill outside the steady start; a restart's first
    // fill from an older-format file fills as the steady start does)
    const bool f0 = !cln_insteady && !(cln_rstold && !cln_rfdone && full);
    cln_rfdone = cln_rfdone || (cln_rstold && full);
    if (shift) {
      RingShift(b);
    } else if (full && ema) {
      // EMA init: J-weighted equal mix of the current and previous set, aux := current.
      // f0: the current set alone, count 1; a restart with converged EMA history
      // (cln_rhist): count N, no warm-up
      const Real n0 = cln_rhist ? vra_ema_n : 1.0;
      cln_rhist = false;
      par_for("cln_ring_einit", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k0, const int j0, const int i0) {
        const int k = k0 - ks, j = j0 - js, i = i0 - is;
        const Real jb = rl(m,qb,CLN_J,k,j,i), jp = rr(m,fr[int{CLN_J}],k,j,i);
        const Real jm = jmin*iwr(m,VET_IW_EN,k0,j0,i0);
        const bool v = !f0 && jb > 0.0 && jb >= jm && jp > 0.0 && jp >= jm;
        const Real ws = v ? jp/jb : 1.0;
        for (int n = 0; n < CLN_N; ++n) {
          const int f = fr[n];
          if (f < 0) continue;
          rr(m,f,k,j,i) = (n == CLN_J) ? (v ? 0.5*(1.0 + ws) : 1.0) :
                          (v ? ws*(rr(m,f,k,j,i) - rl(m,qb,n,k,j,i))/(1.0 + ws) : 0.0);
        }
        rr(m,nfr,k,j,i) = jb;
        rr(m,nfr+1,k,j,i) = rl(m,qb,CLN_D11,k,j,i);
        rr(m,nfr+2,k,j,i) = v ? fmax(2.0, n0) : ((f0 && jb > 0.0 && jb >= jm) ? n0 : 0.0);
      });
    } else if (full) {
      // a fill: in SteadyStart the older sets' offsets from their sweeps; elsewhere
      // zero offsets (J ratio 1), the older sets entering one per block via the shifts
      // (sweeping them at a step's build would add an O(dt) error)
      par_for("cln_ring_fill", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k0, const int j0, const int i0) {
        const int k = k0 - ks, j = j0 - js, i = i0 - is;
        const Real jb = rl(m,qb,CLN_J,k,j,i);
        const Real jm = jmin*iwr(m,VET_IW_EN,k0,j0,i0);
        for (int r = 0; r < nr; ++r) {
          for (int n = 0; n < CLN_N; ++n) {
            const int f = fr[n];
            if (f < 0) continue;
            const Real o = rr(m,r*nfr+f,k,j,i);
            rr(m,r*nfr+f,k,j,i) = f0 ? ((n != CLN_J) ? 0.0 : 1.0) : (n != CLN_J) ?
                                  o - rl(m,qb,n,k,j,i) :
                                  ((jb > 0.0 && jb >= jm && o >= jm) ? o/jb : -1.0);
          }
        }
      });
    } else {
      // a boundary missed: this build takes the swept sets, the next one the union
      cln_roff_sum += 1.0;
      cln_rok = false;
    }
    if (shift || full) {
      cln_rok = true;
      cln_rb = b;
    }
  }
  // the pending Delta of this build
  if (slot[nk] >= 0) {
    par_for("cln_ring_delta", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k0, const int j0, const int i0) {
      const int k = k0 - ks, j = j0 - js, i = i0 - is;
      const Real jb = rl(m,qb,CLN_J,k,j,i), jn = rl(m,qn,CLN_J,k,j,i);
      const Real jm = jmin*iwr(m,VET_IW_EN,k0,j0,i0);
      for (int n = 0; n < CLN_N; ++n) {
        if (fr[n] < 0) continue;
        pd(m,fr[n],k,j,i) = (n != CLN_J) ? rl(m,qn,n,k,j,i) - rl(m,qb,n,k,j,i) :
                                ((jb > 0.0 && jb >= jm && jn >= jm) ? jn/jb : -1.0);
      }
      if (!ema) return;
      // (EMA change weight a_ch from c = max |ln J_b/J_aux|, |D11_b - D11_aux|; the next
      // set's J and D11 for the aux)
      Real ach = -1.0;
      if (eps > 0.0) {
        const Real jo = rr(m,nfr,k,j,i);
        Real c = (jb > 0.0 && jo > 0.0) ? fabs(log(jb/jo)) : 1.0e30;
        c = fmax(c, fabs(rl(m,qb,CLN_D11,k,j,i) - rr(m,nfr+1,k,j,i)));
        ach = estep ? fmin(1.0, fmax(0.0, c/eps - 1.0)) : fmin(1.0, c/eps);
      }
      pd(m,nfr,k,j,i) = ach;
      pd(m,nfr+1,k,j,i) = jn;
      pd(m,nfr+2,k,j,i) = rl(m,qn,CLN_D11,k,j,i);
    });
    cln_rpb = b;
    cln_rpt = (static_cast<Real>(b) + (cln_rlast ? 1.0 : s))*vgd_rot_dt;
  }
  // the mixture; jump: at the block-closing build also the next block's window (cln_rnew)
  const bool ring = cln_rok;
  const bool nxt = ring && cln_rlast && slot[nk] >= 0;
  if (nxt && cln_rnew.extent_int(0) != nmb) {   // (the layout of cln_mc it becomes)
    Kokkos::realloc(cln_rnew, nmb, cln_s.extent_int(1), c3, c2, c1);
  }
  auto rw = cln_rnew;
  auto iw_ = iw;
  auto valid = cln_v;
  SolveArgs x = MakeArgs();
  int nguard = 0;
  const int ni = ie - is + 1, nj = je - js + 1, nkk = ke - ks + 1;
  Kokkos::parallel_reduce("cln_ring_mix",
                          Kokkos::RangePolicy<DevExeSpace>(0, nmb*nkk*nj*ni),
  KOKKOS_LAMBDA(const int t, int &ng) {
    const int i = is + t%ni;
    const int j = js + (t/ni)%nj;
    const int k = ks + (t/(ni*nj))%nkk;
    const int m = t/(ni*nj*nkk);
    // (the swept sets' raw moments M = g J)
    const int kr = k - ks, jr = j - js, ir = i - is;
    // (slot o; fields the ring does not carry: 0, never read)
    auto raw = [&](const int o, const int n) {
      const Real jo = rl(m,o,CLN_J,kr,jr,ir);
      return (n == CLN_J) ? jo : ((fm[n] < 0) ? 0.0 : rl(m,o,n,kr,jr,ir)*jo);
    };
    Real q[CLN_N];
    const bool mix = ring && (ema ? RingEma(rl, rr, pd, fr, nfr, qn, s, a0, m, kr,
                                            jr, ir, q)
                                  : RingMix(rl, rr, fr, nfr, qn, nr, s, m, kr, jr,
                                            ir, q));
    if (!mix) {
      for (int n = 0; n < CLN_N; ++n) q[n] = (1.0 - s)*raw(0, n) + s*raw(qn, n);
    }
    if (mix && !(q[CLN_J] > 0.0 && Realizable(x, q, m, k, j, i))) {
      if ((1.0 - s)*raw(0, CLN_J) + s*raw(qn, CLN_J) > 0.0) ++ng;
      for (int n = 0; n < CLN_N; ++n) {
        if (!x.cub || CellMoment(n)) q[n] = (1.0 - s)*raw(0, n) + s*raw(qn, n);
      }
    }
    if (q[CLN_J] < ray_min*iw_(m,VET_IW_EN,k,j,i)) q[CLN_J] = 0.0;
    StoreMoments(cs, valid, q, m, k, j, i);
    if (!nxt) return;
    // (the next block's window: the union at s = 1)
    if (!(RingMix(rl, rr, fr, nfr, qn, nr, 1.0, m, kr, jr, ir, q) &&
          q[CLN_J] > 0.0 && Realizable(x, q, m, k, j, i))) {
      if (rl(m,qn,CLN_J,kr,jr,ir) > 0.0) ++ng;
      for (int n = 0; n < CLN_N; ++n) {
        if (!x.cub || CellMoment(n)) q[n] = raw(qn, n);
      }
    }
    if (q[CLN_J] < ray_min*iw_(m,VET_IW_EN,k,j,i)) {
      for (int n = 0; n < CLN_N; ++n) {
        if (CellMoment(n)) q[n] = 0.0;
      }
    }
    for (int n = 0; n < cs.extent_int(1); ++n) rw(m,n,kr,jr,ir) = q[n];
  }, Kokkos::Sum<int>(nguard));
  cln_rjump = nxt;
  cln_rguard_sum += GlobalReduce(nguard, false);
  cln_rl = DvceArray5D<Real>();   // (build scratch: freed, every build writes it first)
  if (cln_cssc) cln_fvi = DvceArray5D<Real>();   // (and so cln_fvi)
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RingSets
//! \brief the swept directions (cln_fvd, cln_fvn) for the closure at time tc: the union
//! of the rotating sets, the set weights folded into the quadrature weights (moments are
//! weighted means: single-valued, conservative, realisable).  vet_rotate_avg = K: sets
//! b-K+1 .. b, 1/K each; vet_gd_rotate_xfade adds b+1, weights (1 - s)/K oldest and s/K
//! next, s the fraction of block b passed.  Slot blk mod (K + xfade); a function of t

void RadiationVet::RingSets(const Real tc, const int cyc) {
  const int n = vgd_n;
  const Mesh *pm = pmy_pack->pmesh;
  const int nk = vra_on ? RingK() : 1;
  const bool xf = vgd_xf && !(vra_on && !g_vra_adv);
  int64_t b = RotBlock(cyc);
  Real s = 0.0;
  if (vgd_rot_dt > 0.0) {
    // (without the cross-fade the block of the step's midpoint, as RotBlock)
    const Real dt = (pm->dt < vgd_rot_dt) ? pm->dt : 0.0;
    const Real tau = (tc - (xf ? 0.0 : 0.5*dt))/vgd_rot_dt;
    b = static_cast<int64_t>(std::floor(tau));
    if (xf) s = tau - static_cast<Real>(b);
    // ring = live: a build on a block boundary closes the block (s = 1)
    if (cln_ring == 1 && s < 1.0e-9) {
      b -= 1;
      s = 1.0;
    }
  }
  // ring with its history: the current set and the next one, (1 - s) and s;
  // jump: the next one (weight 0) only at the build that closes the block
  const bool two = cln_ring && cln_rok;
  cln_rlast = (cln_ring == 2) &&
              tc >= (static_cast<Real>(b) + 1.0 - 1.0e-9)*vgd_rot_dt;
  const int nx = (xf || cln_rlast) ? 1 : 0;
  const int ns = (two ? 1 : nk) + nx;
  // (the sets' tables by their own count: RingUnion resets cln_fvn to two sets)
  if (cln_fvn != ns*n || static_cast<int>(cln_fvkey.size()) != ns) {
    cln_fvn = ns*n;
    cln_fvd = DvceArray2D<Real>("cln_fvd", ns*n, 4);
    cln_fvkey.assign(ns, std::numeric_limits<Real>::quiet_NaN());
    cln_fvblk.assign(ns, 0);
    if (cln_ring) {
      cln_fvdm.resize(ns);
      for (int q = 0; q < ns; ++q) cln_fvdm[q] = DvceArray2D<Real>("cln_fvdm", ns*n, 4);
    }
  }
  auto d_h = Kokkos::create_mirror_view(cln_fvd);
  for (int64_t blk = (two ? b : b - nk + 1); blk <= b + nx; ++blk) {
    const int q = static_cast<int>(((blk % ns) + ns) % ns);
    // (Cartesian slab: the turn of the vet_sc set about x1)
    const Real key = (vgd_slab == 1) ? ScSlabKey(blk) : RotAngleBlk(blk);
    Real w = ((xf && blk == b - nk + 1) ? 1.0 - s : ((blk > b) ? s : 1.0))/nk;
    if (two) w = (blk > b) ? s : 1.0 - s;
    const std::vector<double> d = ScRotate(key);
    for (int e = 0; e < n; ++e) {
      for (int c = 0; c < 3; ++c) d_h(q*n+e,c) = d[4*e+c];
      d_h(q*n+e,3) = w*d[4*e+3];
    }
    cln_fvkey[q] = key;
    cln_fvblk[q] = blk;
  }
  Kokkos::deep_copy(cln_fvd, d_h);
  if (cln_ring) {
    // the directions of each slot alone, with the set's own weights
    for (int q = 0; q < ns; ++q) {
      auto m_h = Kokkos::create_mirror_view(cln_fvdm[q]);
      const std::vector<double> d = ScRotate(cln_fvkey[q]);
      for (int e = 0; e < ns*n; ++e) {
        for (int c = 0; c < 3; ++c) m_h(e,c) = d_h(e,c);
        m_h(e,3) = (e/n == q) ? d[4*(e%n)+3] : 0.0;
      }
      Kokkos::deep_copy(cln_fvdm[q], m_h);
    }
    cln_rns = ns;
    cln_rs = s;
    cln_rbld = b;
  }
}

}  // namespace radvet
