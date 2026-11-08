//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_rotavg.cpp
//! \brief rotavg-1012, <rad_m1>/vet_rotate_avg = K >= 2: the closure averaged over the
//! last K rotated direction sets (closure-agnostic; vet_gd registers its fields in
//! VetGdRotAvgSetup, rad_m1_vetgd.cpp).
//!
//! WHY.  With the set rotated every rotate_dt and c dt/dx >> 1, every set relaxes to its
//! own steady state within a step or two, so the closure in use is the CURRENT set's,
//! ray effects included ("last set wins"): rotation does not average the ray effects.
//!
//! WHAT.  g_s = the normalised hand-over fields of set s (ratios to J: the Eddington
//! tensor parts, the half-range ratios, C1/J; and J itself).  The moment step gets
//!   used = g_b(now) + (1/K) sum_{j=1}^{K-1} d_j        (additive fields)
//!   used = g_b(now) (1 + sum_j r_j)/K                  (J: ratio fields)
//! with the OFFSETS d_j = g_{b-j} - g_b (r_j = J_{b-j}/J_b) of the K-1 previous sets from
//! the current one.  Every offset is a difference of two sweeps of the SAME state: at the
//! first build of a new block the old set is swept once more (before the tables switch),
//! Delta = g_new - g_old, and the ring is shifted: d_1' = -Delta, d_j' = d_{j-1} - Delta
//! (r_1' = 1/R, r_j' = r_{j-1}/R with R = J_new/J_old).  On a fresh start or a restart
//! the K-1 previous sets are swept once at the first build.  In a steady state used =
//! (1/K) sum_s g_s exactly: the normalised moments of the equal mixture of the K sets'
//! normalised intensities, i.e. a realizable closure with the tensor and the half-range
//! ratios from the SAME mixture.  In a changing state the current set is live (no lag of
//! the closure); only the ray-pattern differences are carried, aged <= K rotate_dt.
//! Clamps (kind 1): the half-range ratios to their ranges, C1 >= 0, the outer q to its
//! range; the tensor's D_rr is clamped by the fold (VetLatBuild).  Field arrays: two
//! (m, n, k, j, i) arrays and one column array (m, k, j) held at the top cell i = ie.
//! Kinds: 0 additive, 1 additive clamped to [lo, hi], 2 ratio (J itself), 3 a raw moment
//! normalised by the J field (vet_sc's K_ab, H_a, the half-range sums, C1 in E units):
//! averaged as X/J, handed back as (X/J)_avg J_avg.
//!
//! VALIDITY.  A cell is valid in a set when i >= the closure's first shell (ilo) and the
//! flag field (if any) is >= 0 (and the J field > 0 when kind-3 fields exist); an offset is kept only where both sweeps were valid, and
//! the mean runs over the valid offsets (+ the current set).
//!
//! COST.  One extra sweep per rotation block (+ K-1 sweeps at the first build); memory
//! (K-1)(nf+1) reals per active cell.


#include <algorithm>
#include <cstdio>
#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "rad_m1/rad_m1.hpp"

namespace radm1 {

namespace {
// the field arrays of one call (device copies)
struct RaArr {
  DvceArray5D<Real> a0, a1;
  DvceArray3D<Real> c;
  KOKKOS_INLINE_FUNCTION
  Real &At(const int arr, const int comp, const int m, const int k, const int j,
           const int i) const {
    return (arr == 2) ? c(m,k,j) : ((arr == 0) ? a0(m,comp,k,j,i) : a1(m,comp,k,j,i));
  }
};

using Tab = RadiationM1::RotAvgTab;

KOKKOS_INLINE_FUNCTION
Real RaJ(const RaArr &a, const Tab &t, const int m, const int k, const int j,
         const int i) {
  return (t.jf >= 0) ? a.At(t.arr[t.jf], t.comp[t.jf], m, k, j, i) : 1.0;
}

KOKKOS_INLINE_FUNCTION
bool RaValid(const RaArr &a, const Tab &t, const int ilo, const int m, const int k,
             const int j, const int i) {
  if (i < ilo) return false;
  if (t.farr >= 0 && a.At(t.farr, t.fcomp, m, k, j, i) < 0.0) return false;
  if (t.jf >= 0 && !(RaJ(a, t, m, k, j, i) > 0.0)) return false;
  return true;
}

// field f of the cell in its averaged form (kind 3: / J), or false when not held here
KOKKOS_INLINE_FUNCTION
bool RaGet(const RaArr &a, const Tab &t, const int f, const int ie, const int m,
           const int k, const int j, const int i, const Real jr, Real &g) {
  if (t.arr[f] == 2 && i != ie) return false;
  g = a.At(t.arr[f], t.comp[f], m, k, j, i);
  if (t.kind[f] == 3) g = (jr > 0.0) ? g/jr : 0.0;
  return true;
}
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::RotAvgSnap
//! \brief the registered fields as they are (raw, one set; kind 3 as X/J) -> ring slot
//! `slot`, with the validity of each cell

void RadiationM1::RotAvgSnap(const int slot, const int ilo) {
  vra_gen += 1;   // gpumem-1022: vra_d may change (dt_retry_host skips clean saves)
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nf = vra_tab.nf;
  const int nalloc = vra_k;   // the ring (K-1 slots) + the Delta slot K-1
  if (vra_d.extent_int(0) < nmb1 + 1 || vra_d.extent_int(1) != nalloc*(nf + 1)) {
    Kokkos::realloc(vra_d, nmb1 + 1, nalloc*(nf + 1), ke - ks + 1, je - js + 1,
                    ie - is + 1);
    Kokkos::deep_copy(vra_d, 0.0);
  }
  auto d_ = vra_d;
  const RaArr a{vra_a0, vra_a1, vra_c};
  const Tab t = vra_tab;
  const int o = slot*(nf + 1);
  par_for("m1_ra_snap", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const bool v = RaValid(a, t, ilo, m, k, j, i);
    const Real jr = RaJ(a, t, m, k, j, i);
    for (int f = 0; f < t.nf; ++f) {
      Real g;
      if (RaGet(a, t, f, ie, m, k, j, i, jr, g)) d_(m,o+f,k-ks,j-js,i-is) = g;
    }
    d_(m,o+t.nf,k-ks,j-js,i-is) = v ? 1.0 : 0.0;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::RotAvgUpdate
//! \brief init (the first build of a run): with the CURRENT set's raw fields in place, the
//! raw snapshots S_j of the K-1 previous sets -> offsets S_j - g (ratio fields S_j/g)

void RadiationM1::RotAvgUpdate(const int ilo) {
  if (vra_pend != 1) return;
  if (vra_ema) {
    RotAvgEmaInit(ilo);
    return;
  }
  vra_gen += 1;   // gpumem-1022: vra_d may change (dt_retry_host skips clean saves)
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nslot = vra_k - 1;
  auto d_ = vra_d;
  const RaArr a{vra_a0, vra_a1, vra_c};
  const Tab t = vra_tab;
  const int nf1 = t.nf + 1;
  const bool init = (vra_pend == 1);
  par_for("m1_ra_upd", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int kk = k - ks, jj = j - js, ii = i - is;
    const bool vg = RaValid(a, t, ilo, m, k, j, i);
    const Real jr = RaJ(a, t, m, k, j, i);
    if (init) {
      for (int s = 0; s < nslot; ++s) {
        const int o = s*nf1;
        const bool v = vg && (d_(m,o+t.nf,kk,jj,ii) > 0.5);
        for (int f = 0; f < t.nf; ++f) {
          Real g;
          if (!RaGet(a, t, f, ie, m, k, j, i, jr, g)) continue;
          const Real sv = d_(m,o+f,kk,jj,ii);
          d_(m,o+f,kk,jj,ii) = (t.kind[f] == 2) ? ((v && g > 0.0 && sv > 0.0) ? sv/g : 1.0)
                                                : (v ? (sv - g) : 0.0);
        }
        d_(m,o+t.nf,kk,jj,ii) = v ? 1.0 : 0.0;
      }
    }
  });
  vra_pend = 0;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::RotAvgDelta
//! \brief the LAST build of a block: slot holds the current set's raw fields g_b, the
//! fields now hold the NEXT set's raw S at the same state -> slot := Delta = S - g_b
//! (ratio fields S/g_b), valid where both are.  The state of the last build is the
//! predicted P ~ U(t_k) + O(dt^2) at the block boundary t_k on every dt (the first build
//! of the next block would sample U(t_k + dt): an O(dt) offset on a dt ladder)

void RadiationM1::RotAvgDelta(const int slot, const int ilo) {
  if (vra_ema) {RotAvgEmaChange(ilo);}
  vra_gen += 1;   // gpumem-1022: vra_d may change (dt_retry_host skips clean saves)
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto d_ = vra_d;
  const RaArr a{vra_a0, vra_a1, vra_c};
  const Tab t = vra_tab;
  const int o = slot*(t.nf + 1);
  par_for("m1_ra_delta", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int kk = k - ks, jj = j - js, ii = i - is;
    const bool v = RaValid(a, t, ilo, m, k, j, i) && (d_(m,o+t.nf,kk,jj,ii) > 0.5);
    const Real jr = RaJ(a, t, m, k, j, i);
    for (int f = 0; f < t.nf; ++f) {
      Real sn;
      if (!RaGet(a, t, f, ie, m, k, j, i, jr, sn)) continue;
      const Real g = d_(m,o+f,kk,jj,ii);
      d_(m,o+f,kk,jj,ii) = (t.kind[f] == 2) ? ((v && g > 0.0 && sn > 0.0) ? sn/g : 1.0)
                                            : (v ? (sn - g) : 0.0);
    }
    d_(m,o+t.nf,kk,jj,ii) = v ? 1.0 : 0.0;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::RotAvgShift
//! \brief the first build of a new block: the ring moves on by one set with the Delta of
//! the last build (slot ds): the dropped oldest slot becomes d_1 = -Delta (1/R), every
//! other d_j -= Delta (/= R)

void RadiationM1::RotAvgShift(const int ds) {
  if (vra_ema) {
    RotAvgEmaShift(ds);
    return;
  }
  vra_gen += 1;   // gpumem-1022: vra_d may change (dt_retry_host skips clean saves)
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nslot = vra_k - 1;
  vra_head = (vra_head + nslot - 1) % nslot;
  const int h = vra_head;
  auto d_ = vra_d;
  const Tab t = vra_tab;
  const int nf1 = t.nf + 1;
  const int od = ds*nf1;
  par_for("m1_ra_shift", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int kk = k - ks, jj = j - js, ii = i - is;
    const bool vd = d_(m,od+t.nf,kk,jj,ii) > 0.5;
    for (int s = 0; s < nslot; ++s) {
      const int o = s*nf1;
      const bool v = vd && ((s == h) || d_(m,o+t.nf,kk,jj,ii) > 0.5);
      for (int f = 0; f < t.nf; ++f) {
        const bool mul = (t.kind[f] == 2);
        const Real x = d_(m,od+f,kk,jj,ii);
        Real y;
        if (!v) {
          y = mul ? 1.0 : 0.0;
        } else if (s == h) {
          y = mul ? 1.0/x : -x;
        } else {
          y = mul ? d_(m,o+f,kk,jj,ii)/x : d_(m,o+f,kk,jj,ii) - x;
        }
        d_(m,o+f,kk,jj,ii) = y;
      }
      d_(m,o+t.nf,kk,jj,ii) = v ? 1.0 : 0.0;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::RotAvgApply
//! \brief the registered fields := the mean over the current set and the valid offsets
//! (ratio fields: the plain mean; the others J-weighted with the ratio field jw, i.e.
//! the normalised moments of the union of the K sets; kind 3: (X/J)_avg times the mean J)

void RadiationM1::RotAvgApply(const int ilo) {
  if (vra_ema) {
    RotAvgEmaApply(ilo);
    return;
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nslot = vra_k - 1;
  auto d_ = vra_d;
  const RaArr a{vra_a0, vra_a1, vra_c};
  const Tab t = vra_tab;
  const int nf1 = t.nf + 1;
  par_for("m1_ra_apply", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (!RaValid(a, t, ilo, m, k, j, i)) return;
    const int kk = k - ks, jj = j - js, ii = i - is;
    Real nv = 1.0;
    for (int s = 0; s < nslot; ++s) {nv += d_(m,s*nf1+t.nf,kk,jj,ii);}
    if (nv < 1.5) return;
    const Real jr = RaJ(a, t, m, k, j, i);
    Real ju = jr;   // the averaged J (kind 3 fields are handed back with it)
    if (t.jf >= 0) {
      Real sm = 1.0;
      for (int s = 0; s < nslot; ++s) {
        const int o = s*nf1;
        if (d_(m,o+t.nf,kk,jj,ii) > 0.5) sm += d_(m,o+t.jf,kk,jj,ii);
      }
      ju = jr*sm/nv;
    }
    // the weight of each older set: its J over the current set's (field jw, a ratio
    // field), so the normalised fields are those of the UNION of the sets' intensities,
    // sum_s M_s / sum_s J_s (equal weights 1 without a J field or vet_rotate_avg_weight =
    // equal)
    Real ws[M1_RA_NF], wsum = 1.0;
    for (int s = 0; s < nslot; ++s) {
      const int o = s*nf1;
      ws[s] = (d_(m,o+t.nf,kk,jj,ii) > 0.5) ? ((t.jw >= 0) ? d_(m,o+t.jw,kk,jj,ii) : 1.0)
                                             : 0.0;
      wsum += ws[s];
    }
    for (int f = 0; f < t.nf; ++f) {
      Real g;
      if (!RaGet(a, t, f, ie, m, k, j, i, jr, g)) continue;
      const bool mul = (t.kind[f] == 2);
      Real sm = mul ? 1.0 : 0.0;
      for (int s = 0; s < nslot; ++s) {
        const int o = s*nf1;
        if (d_(m,o+t.nf,kk,jj,ii) > 0.5) sm += mul ? d_(m,o+f,kk,jj,ii)
                                                    : ws[s]*d_(m,o+f,kk,jj,ii);
      }
      Real v;
      if (mul) {
        v = g*sm/nv;
      } else {
        v = g + sm/wsum;
        if (t.kind[f] == 1 || (t.kind[f] == 3 && t.lo[f] < t.hi[f])) {
          v = fmin(fmax(v, t.lo[f]), t.hi[f]);
        }
        if (t.kind[f] == 3) v *= ju;
      }
      a.At(t.arr[f], t.comp[f], m, k, j, i) = v;
    }
  });
}

//----------------------------------------------------------------------------------------
//! raycorr-1108 vet_rotate_avg_ema = N: the exponential moving average over the rotated
//! sets in the offset form of the K-window (one ring slot, vra_k = 2).  With weights
//! w_j = a (1 - a)^j over the sets b - j (J-weighted as the K-window: the normalised
//! moments of the weighted union of the sets' intensities), slot 0 holds
//!   e_f = avg_f - g_f(now)  (additive fields),  r_f = avg_f/g_f(now)  (ratio fields: J)
//! and the boundary b-1 -> b with the same-state Delta (D_f = g_new - g_old, R_f =
//! g_new/g_old) updates them as
//!   r' = (1 - a) r/R + a,   e' = (1 - a) (r_jw/R_jw) (e - D)/r'_jw
//! (the factor r_jw/R_jw / r'_jw is 1 with equal weights).  In a steady state used = the
//! EMA of the sets; the current set stays live, only pattern differences are carried
//! (they decay as (1 - a)^blocks).  Init (the first build of a run): the equal mix of the
//! current and the one previous set (as K = 2).

namespace {
// the EMA weight of block boundary for a cell at radius rr: a = max(1/N, 1 - exp(-dtb/tau))
// with tau = tau0 (physical) and/or lc r/c (local light crossing); 1/N without either
KOKKOS_INLINE_FUNCTION
Real RaEmaA(const Real a0, const Real dtb, const Real tau0, const Real lc, const Real rr,
            const Real cl) {
  Real a = a0;
  if (tau0 > 0.0) a = fmax(a, 1.0 - exp(-dtb/tau0));
  if (lc > 0.0) a = fmax(a, 1.0 - exp(-dtb*cl/(lc*fmax(rr, 1.0e-30))));
  return fmin(a, 1.0);
}
}  // namespace

void RadiationM1::RotAvgEmaInit(const int ilo) {
  vra_gen += 1;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto d_ = vra_d;
  const RaArr a{vra_a0, vra_a1, vra_c};
  const Tab t = vra_tab;
  par_for("m1_ra_ema_init", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int kk = k - ks, jj = j - js, ii = i - is;
    const bool v = RaValid(a, t, ilo, m, k, j, i) && (d_(m,t.nf,kk,jj,ii) > 0.5);
    const Real jr = RaJ(a, t, m, k, j, i);
    // the previous set's weight relative to the current one (its J ratio, or 1)
    Real ws = 1.0;
    if (v && t.jw >= 0) {
      Real g;
      if (RaGet(a, t, t.jw, ie, m, k, j, i, jr, g)) {
        const Real sv = d_(m,t.jw,kk,jj,ii);
        ws = (g > 0.0 && sv > 0.0) ? sv/g : 1.0;
      }
    }
    for (int f = 0; f < t.nf; ++f) {
      Real g;
      if (!RaGet(a, t, f, ie, m, k, j, i, jr, g)) continue;
      const Real sv = d_(m,f,kk,jj,ii);
      if (t.kind[f] == 2) {
        d_(m,f,kk,jj,ii) = (v && g > 0.0 && sv > 0.0) ? 0.5*(1.0 + sv/g) : 1.0;
      } else {
        d_(m,f,kk,jj,ii) = v ? ws*(sv - g)/(1.0 + ws) : 0.0;
      }
    }
    d_(m,t.nf,kk,jj,ii) = v ? 2.0 : 0.0;   // the flag counts the sets in the EMA (warm-up)
    // the aux slot 1: the current set's raw fields (the change test of the next
    // boundary compares the same set at two times); its flag entry holds the change
    // weight of the coming boundary (-1: none measured)
    const int o1 = t.nf + 1;
    for (int f = 0; f < t.nf; ++f) {
      Real g;
      if (RaGet(a, t, f, ie, m, k, j, i, jr, g)) d_(m,o1+f,kk,jj,ii) = g;
    }
    d_(m,o1+t.nf,kk,jj,ii) = -1.0;
  });
  vra_pend = 0;
}

//! raycorr-1108: the change test on the LAST build of a block, before RotAvgDelta turns
//! the Delta slot (vra_k - 1, holding the current set b's raw fields g_b(now)) into the
//! Delta: the same set at the previous boundary is in the aux slot 1 (it was the "next"
//! set there).  c = max(|ln J_b(now)/J_b(then)|, |D_b(now) - D_b(then)|) (D = field 0, the
//! D_rr part) gives the boundary's weight a_ch = min(1, c/eps) (vet_rotate_avg_ema_change
//! = eps, 0 = off), stored in the aux flag entry; then the aux slot := the NEXT set's raw
//! fields (now in place) for the following boundary.  A changed geometry (a source
//! switched on, a shadow moved) thereby drops the history that belongs to the old one.
void RadiationM1::RotAvgEmaChange(const int ilo) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto d_ = vra_d;
  const RaArr a{vra_a0, vra_a1, vra_c};
  const Tab t = vra_tab;
  const int o1 = t.nf + 1, od = (vra_k - 1)*(t.nf + 1);
  const Real eps = vra_ema_eps;
  const bool usej = vra_ema_chj, step = vra_ema_chstep;
  par_for("m1_ra_ema_chg", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int kk = k - ks, jj = j - js, ii = i - is;
    Real ach = -1.0;
    if (eps > 0.0) {
      Real c = 0.0;
      if (usej && t.jf >= 0 && t.arr[t.jf] != 2) {
        const Real jn = d_(m,od+t.jf,kk,jj,ii), jo = d_(m,o1+t.jf,kk,jj,ii);
        c = (jn > 0.0 && jo > 0.0) ? fabs(log(jn/jo)) : 1.0e30;
      }
      if (t.nf > 0 && t.arr[0] != 2 && t.kind[0] != 2) {
        c = fmax(c, fabs(d_(m,od,kk,jj,ii) - d_(m,o1,kk,jj,ii)));
      }
      // linear: a_ch = c/eps (any drift shortens the memory); step: 0 below eps, 1 at
      // 2 eps (only a change larger than eps per block drops the history)
      ach = step ? fmin(1.0, fmax(0.0, c/eps - 1.0)) : fmin(1.0, c/eps);
    }
    const Real jr = RaJ(a, t, m, k, j, i);
    for (int f = 0; f < t.nf; ++f) {
      Real g;
      if (RaGet(a, t, f, ie, m, k, j, i, jr, g)) d_(m,o1+f,kk,jj,ii) = g;
    }
    d_(m,o1+t.nf,kk,jj,ii) = ach;
  });
}

void RadiationM1::RotAvgEmaShift(const int ds) {
  vra_gen += 1;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto d_ = vra_d;
  const Tab t = vra_tab;
  const int od = ds*(t.nf + 1);
  const Real a0 = 1.0/vra_ema_n;
  const Real tau0 = vra_ema_tau, lc = vra_ema_lc, cl = c_light;
  // the block length: the physical rotation clock, else the cycles of a block x dt
  Mesh *pm = pmy_pack->pmesh;
  const Real dtb = (vgd_rot_dt > 0.0) ? vgd_rot_dt
                                      : static_cast<Real>(std::max(vgd_rot, 1))*pm->dt;
  auto &size = pmy_pack->pmb->mb_size;
  const int nx1 = indcs.nx1;
  // raycorr-1108 vet_rotate_avg_ema_diag = M > 0: every M-th boundary, the change
  // weights of the cells with a history (mean, fraction > 1/N, fraction >= 0.5), printed
  vra_ema_nsh += 1;
  if (vra_ema_diag > 0 && vra_ema_nsh % vra_ema_diag == 0) {
    const int nkji = (ke-ks+1)*(je-js+1)*(ie-is+1);
    const int nji = (je-js+1)*(ie-is+1), ni = ie-is+1;
    Real sn = 0.0, sa = 0.0, s1 = 0.0, s2 = 0.0;
    Kokkos::parallel_reduce("m1_ra_ema_diag",
    Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1+1)*nkji),
    KOKKOS_LAMBDA(const int idx, Real &rn, Real &ra, Real &r1, Real &r2) {
      const int m = idx/nkji;
      const int kk = (idx - m*nkji)/nji;
      const int jj = (idx - m*nkji - kk*nji)/ni;
      const int ii = idx - m*nkji - kk*nji - jj*ni;
      if (!(d_(m,od+t.nf,kk,jj,ii) > 0.5) || !(d_(m,t.nf,kk,jj,ii) > 0.5)) return;
      const Real ach = fmax(d_(m,(t.nf+1)+t.nf,kk,jj,ii), 0.0);
      rn += 1.0;
      ra += ach;
      if (ach > a0) r1 += 1.0;
      if (ach >= 0.5) r2 += 1.0;
    }, Kokkos::Sum<Real>(sn), Kokkos::Sum<Real>(sa), Kokkos::Sum<Real>(s1),
       Kokkos::Sum<Real>(s2));
    Real v[4] = {sn, sa, s1, s2};
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, v, 4, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (global_variable::my_rank == 0) {
      const Real n = fmax(v[0], 1.0);
      std::printf("<rad_m1> ema_diag cycle %d: cells %.0f, mean a_ch %.4e, "
                  "f(a_ch > 1/N) %.4f, f(a_ch >= 0.5) %.4f\n", pm->ncycle, v[0],
                  v[1]/n, v[2]/n, v[3]/n);
    }
  }
  par_for("m1_ra_ema_shift", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int kk = k - ks, jj = j - js, ii = i - is;
    const bool vd = d_(m,od+t.nf,kk,jj,ii) > 0.5;
    const bool vh = d_(m,t.nf,kk,jj,ii) > 0.5;
    if (!vd || !vh) {   // no Delta or no history: restart from the current set alone
      for (int f = 0; f < t.nf; ++f) {d_(m,f,kk,jj,ii) = (t.kind[f] == 2) ? 1.0 : 0.0;}
      d_(m,t.nf,kk,jj,ii) = vd ? 1.0 : 0.0;   // count 1: the current set alone
      return;
    }
    Real rr = 1.0;
    if (lc > 0.0) {   // the cell radius (x1 = r on the sp wedge and the cubed sphere)
      const Real x1min = size.d_view(m).x1min, x1max = size.d_view(m).x1max;
      rr = x1min + (static_cast<Real>(ii) + 0.5)*(x1max - x1min)/nx1;
    }
    const Real ach = d_(m,(t.nf+1)+t.nf,kk,jj,ii);   // the change weight (aux slot)
    // warm-up: with n sets seen so far the weight is at least 1/(n+1) (the plain mean of
    // all sets until n reaches N; the flag entry counts them)
    const Real nset = d_(m,t.nf,kk,jj,ii);
    const Real aw = fmax(fmax(RaEmaA(a0, dtb, tau0, lc, rr, cl), fmin(ach, 1.0)),
                         1.0/(nset + 1.0));
    Real fac = 1.0 - aw;
    if (t.jw >= 0) {
      const Real r0 = d_(m,t.jw,kk,jj,ii), R = d_(m,od+t.jw,kk,jj,ii);
      const Real r1 = (1.0 - aw)*r0/R + aw;
      fac = (1.0 - aw)*(r0/R)/r1;
    }
    for (int f = 0; f < t.nf; ++f) {
      const Real x = d_(m,od+f,kk,jj,ii);
      const Real y = d_(m,f,kk,jj,ii);
      d_(m,f,kk,jj,ii) = (t.kind[f] == 2) ? ((1.0 - aw)*y/x + aw) : fac*(y - x);
    }
    d_(m,t.nf,kk,jj,ii) = (aw >= 1.0) ? 1.0 : fmin(nset + 1.0, 1.0e9);
  });
}

void RadiationM1::RotAvgEmaApply(const int ilo) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto d_ = vra_d;
  const RaArr a{vra_a0, vra_a1, vra_c};
  const Tab t = vra_tab;
  par_for("m1_ra_ema_apply", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (!RaValid(a, t, ilo, m, k, j, i)) return;
    const int kk = k - ks, jj = j - js, ii = i - is;
    if (!(d_(m,t.nf,kk,jj,ii) > 0.5)) return;
    const Real jr = RaJ(a, t, m, k, j, i);
    const Real ju = (t.jf >= 0) ? jr*d_(m,t.jf,kk,jj,ii) : jr;
    for (int f = 0; f < t.nf; ++f) {
      Real g;
      if (!RaGet(a, t, f, ie, m, k, j, i, jr, g)) continue;
      Real v;
      if (t.kind[f] == 2) {
        v = g*d_(m,f,kk,jj,ii);
      } else {
        v = g + d_(m,f,kk,jj,ii);
        if (t.kind[f] == 1 || (t.kind[f] == 3 && t.lo[f] < t.hi[f])) {
          v = fmin(fmax(v, t.lo[f]), t.hi[f]);
        }
        if (t.kind[f] == 3) v *= ju;
      }
      a.At(t.arr[f], t.comp[f], m, k, j, i) = v;
    }
  });
}

//----------------------------------------------------------------------------------------
//! rotavg-rst-1012: the averaging state in the restart file (radm1::kM1RaRstMagic), so a
//! restart with vet_rotate_avg = K is a bitwise continuation.  The header (int64 x
//! M1_RA_RST_NH): have = 1, nslab, K, nf, closure (1 vet_gd, 2 vet_sc), head, pend,
//! nslot (ring slots written: K-1, or K with the Delta slot when pend = 3), blk, the
//! interior cells per MeshBlock.  The slots' interior values (nslot (nf+1) nx3 nx2 nx1
//! per MeshBlock) are packed densely into nslab slabs of the record's (ghosted) slab
//! size, the LAST of each MeshBlock record, so the history moves with its MeshBlock (any
//! rank count / distribution).  Written only with a complete history (vra_state = 1).
//! rotavg-rst2-1013: only the radial range i >= i0 is written, i0 = the first interior
//! i (over all MeshBlocks and ranks) where a written slot is valid: the closure's shells
//! (vet_gd: i >= is + vlat_icut).  Cells below are invalid in every slot and never read;
//! they load as zero (flag 0).  Header entry 10 = i0 (a 10-entry header: i0 = 0, the
//! full range, as written by rotavg-rst-1012).

// the registered field count of the history: the table's, or (a restart written before
// the first build of the run, e.g. a 0-step restart) the one loaded from the file
int RadiationM1::RotAvgNf() const {
  return (vra_tab.nf > 0) ? vra_tab.nf : vra_rst_nf;
}

// the first interior i where a written slot is valid, over all MeshBlocks and ranks
// (nx1 - 1 at most, so at least one radial cell is written)
int RadiationM1::RotAvgRstI0() const {
  const int nf = RotAvgNf();
  const int nslot = (vra_pend == 3) ? vra_k : vra_k - 1;
  auto d_ = vra_d;
  const int nmb = vra_d.extent_int(0), nk = vra_d.extent_int(2);
  const int nj = vra_d.extent_int(3), ni = vra_d.extent_int(4);
  int i0 = ni - 1;
  if (nmb > 0 && nf > 0) {
    const int nf1 = nf + 1;
    Kokkos::parallel_reduce("m1_ra_rst_i0",
        Kokkos::RangePolicy<DevExeSpace>(0, nmb*nk*nj),
    KOKKOS_LAMBDA(const int q, int &lmin) {
      const int m = q/(nk*nj);
      const int k = (q/nj) % nk;
      const int j = q % nj;
      for (int i = 0; i < lmin; ++i) {
        bool v = false;
        for (int s = 0; s < nslot; ++s) {v = v || (d_(m,s*nf1+nf,k,j,i) > 0.5);}
        if (v) {lmin = i; break;}
      }
    }, Kokkos::Min<int>(i0));
  }
  i0 = std::min(i0, ni - 1);
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &i0, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif
  return i0;
}

int RadiationM1::RotAvgRstNslab(const int ncell) const {
  const int nf = RotAvgNf();
  if (!vra_on || vra_state != 1 || nf <= 0 || vra_d.extent_int(0) == 0 || ncell <= 0) {
    return 0;
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int64_t nslot = (vra_pend == 3) ? vra_k : vra_k - 1;
  const int i0 = RotAvgRstI0();
  const int64_t nv = nslot*(nf + 1)*static_cast<int64_t>(indcs.nx1 - i0)*indcs.nx2*
                     indcs.nx3;
  return static_cast<int>((nv + ncell - 1)/ncell);
}

void RadiationM1::RotAvgRstHdr(int64_t *h) const {
  h[10] = 0;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  h[0] = 1;
  h[2] = vra_k;
  h[3] = RotAvgNf();
  h[4] = (vgd_on ? 1 : 2) + (vra_ema ? 16 : 0);   // raycorr-1108: EMA history
  h[5] = vra_head;
  h[6] = vra_pend;
  h[7] = (vra_pend == 3) ? vra_k : vra_k - 1;
  h[8] = vra_blk;
  h[9] = static_cast<int64_t>(indcs.nx1)*indcs.nx2*indcs.nx3;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  h[1] = RotAvgRstNslab((indcs.nx1 + 2*indcs.ng)*n2*n3);
  if (h[1] > 0) {h[10] = RotAvgRstI0();}
}

namespace {
// the flat copy between vra_d's first nslot (nf+1) components (interior, radial range
// i >= i0) and the dense slabs of the record (dir 0: vra_d -> slabs, zero-padded; 1:
// slabs -> vra_d)
void RaRstCopy(const DvceArray5D<Real> &d, const DvceArray5D<Real> &sl, const int nmb,
               const int ncomp, const int dir, const int i0) {
  const int ni = d.extent_int(4) - i0, nj = d.extent_int(3), nk = d.extent_int(2);
  const int n1 = sl.extent_int(4), n2 = sl.extent_int(3), n3 = sl.extent_int(2);
  const int64_t ncint = static_cast<int64_t>(ni)*nj*nk;
  const int64_t ncout = static_cast<int64_t>(n1)*n2*n3;
  const int64_t nv = ncomp*ncint;
  const int64_t ntot = static_cast<int64_t>(sl.extent_int(1))*ncout;
  // fin51-1014: one MeshBlock per launch, and launches of at most 2^30 values (the
  // index range is int; a block's record can exceed 2^31 values); one launch = as before
  const int64_t nch = int64_t(1) << 30;
  for (int m = 0; m < nmb; ++m) {
  for (int64_t q0 = 0; q0 < ntot; q0 += nch) {
  const int64_t nq = std::min(nch, ntot - q0);
  par_for("m1_ra_rst", DevExeSpace(), 0, static_cast<int>(nq - 1),
  KOKKOS_LAMBDA(const int qi) {
    const int64_t q = q0 + qi;
    const int s = static_cast<int>(q/ncout);
    const int64_t r2 = q - s*ncout;
    const int k = static_cast<int>(r2/(static_cast<int64_t>(n2)*n1));
    const int j = static_cast<int>((r2/n1) % n2);
    const int i = static_cast<int>(r2 % n1);
    if (q >= nv) {
      if (dir == 0) {sl(m,s,k,j,i) = 0.0;}
      return;
    }
    const int n = static_cast<int>(q/ncint);
    const int64_t r = q - n*ncint;
    const int kk = static_cast<int>(r/(static_cast<int64_t>(nj)*ni));
    const int jj = static_cast<int>((r/ni) % nj);
    const int ii = static_cast<int>(r % ni) + i0;
    if (dir == 0) {
      sl(m,s,k,j,i) = d(m,n,kk,jj,ii);
    } else {
      d(m,n,kk,jj,ii) = sl(m,s,k,j,i);
    }
  });
  }
  }
}
}  // namespace

void RadiationM1::RotAvgRstPack(DvceArray5D<Real> &dst, const int nmb) {
  RaRstCopy(vra_d, dst, nmb, ((vra_pend == 3) ? vra_k : vra_k - 1)*(RotAvgNf() + 1), 0,
            RotAvgRstI0());
}

//! mem-1014: the same record straight into the host array of the restart output
//! (nmb, nslab, n3, n2, n1), through a device buffer of at most nsc slabs of one
//! MeshBlock instead of the whole (nmb, nslab) record: every element is the value
//! RaRstCopy(dir 0) puts there (the same index map), so the file is byte-identical
void RadiationM1::RotAvgRstPackHost(HostArray5D<Real> &out, const int nmb, const int nsc) {
  const DvceArray5D<Real> d = vra_d;
  const int ncomp = ((vra_pend == 3) ? vra_k : vra_k - 1)*(RotAvgNf() + 1);
  const int i0 = RotAvgRstI0();
  const int ni = d.extent_int(4) - i0, nj = d.extent_int(3), nk = d.extent_int(2);
  const int n1 = out.extent_int(4), n2 = out.extent_int(3), n3 = out.extent_int(2);
  const int nsl = out.extent_int(1);
  const int64_t ncint = static_cast<int64_t>(ni)*nj*nk;
  const int64_t ncout = static_cast<int64_t>(n1)*n2*n3;
  const int64_t nv = ncomp*ncint;
  const int cs = std::max(1, std::min(nsc, nsl));
  DvceArray5D<Real> sl("rst-m1r-chunk", 1, cs, n3, n2, n1);
  for (int m = 0; m < nmb; ++m) {
    for (int s0 = 0; s0 < nsl; s0 += cs) {
      const int ns = std::min(cs, nsl - s0);
      const int64_t q0 = static_cast<int64_t>(s0)*ncout;
      par_for("m1_ra_rst", DevExeSpace(), 0, static_cast<int>(ns*ncout - 1),
      KOKKOS_LAMBDA(const int qi) {
        const int64_t q = q0 + qi;
        const int s = static_cast<int>(q/ncout);
        const int64_t r2 = q - s*ncout;
        const int k = static_cast<int>(r2/(static_cast<int64_t>(n2)*n1));
        const int j = static_cast<int>((r2/n1) % n2);
        const int i = static_cast<int>(r2 % n1);
        if (q >= nv) {
          sl(0,s-s0,k,j,i) = 0.0;
          return;
        }
        const int n = static_cast<int>(q/ncint);
        const int64_t r = q - n*ncint;
        const int kk = static_cast<int>(r/(static_cast<int64_t>(nj)*ni));
        const int jj = static_cast<int>((r/ni) % nj);
        const int ii = static_cast<int>(r % ni) + i0;
        sl(0,s-s0,k,j,i) = d(m,n,kk,jj,ii);
      });
      // (4-D subviews, both contiguous with equal strides: a device -> host copy)
      Kokkos::deep_copy(Kokkos::subview(out, m, std::make_pair(s0, s0+ns), Kokkos::ALL,
                                        Kokkos::ALL, Kokkos::ALL),
                        Kokkos::subview(sl, 0, std::make_pair(0, ns), Kokkos::ALL,
                                        Kokkos::ALL, Kokkos::ALL));
    }
  }
}

//! the history of the file into vra_d + the scalars, when it fits this run (the same K
//! and closure, the same MeshBlock size); else nothing (the first build re-measures)
bool RadiationM1::RotAvgRstLoad(const int64_t *h, const DvceArray5D<Real> &src,
                                const int nmb) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int64_t ncint = static_cast<int64_t>(indcs.nx1)*indcs.nx2*indcs.nx3;
  const int64_t close = (vgd_on ? 1 : 2) + (vra_ema ? 16 : 0);
  const bool ok = vra_on && h[0] == 1 && h[2] == vra_k && h[4] == close && h[3] > 0 &&
                  h[3] < M1_RA_NF && h[9] == ncint && h[5] >= 0 && h[5] < vra_k - 1 &&
                  (h[6] == 0 || h[6] == 3) && h[7] == ((h[6] == 3) ? vra_k : vra_k - 1) &&
                  h[10] >= 0 && h[10] < indcs.nx1;
  if (!ok) {
    if (global_variable::my_rank == 0) {
      std::cout << "### WARNING: the <rad_m1> vet_rotate_avg history of this restart file "
                << "does not fit this run (K " << h[2] << ", closure " << h[4] << "); it "
                << "is re-measured at the first build and this restart is not bitwise."
                << std::endl;
    }
    return false;
  }
  const int nf = static_cast<int>(h[3]);
  vra_gen += 1;   // gpumem-1022: vra_d may change (dt_retry_host skips clean saves)
  Kokkos::realloc(vra_d, nmb, vra_k*(nf + 1), indcs.nx3, indcs.nx2, indcs.nx1);
  Kokkos::deep_copy(vra_d, 0.0);
  RaRstCopy(vra_d, src, nmb, static_cast<int>(h[7])*(nf + 1), 1, static_cast<int>(h[10]));
  vra_head = static_cast<int>(h[5]);
  vra_pend = static_cast<int>(h[6]);
  vra_blk = h[8];
  vra_state = 1;
  vra_rst_nf = nf;
  return true;
}

//! the first averaging call after a load (the table is complete there): the registered
//! fields must be those of the file, else the history is dropped (re-measured)
bool RadiationM1::RotAvgRstCheck() {
  if (vra_rst_nf < 0) {return true;}
  const bool ok = (vra_tab.nf == vra_rst_nf);
  vra_rst_nf = -1;
  if (!ok) {
    if (global_variable::my_rank == 0) {
      std::cout << "### WARNING: <rad_m1> vet_rotate_avg: the restart file's history has "
                << "other fields than this run's closure; re-measured at the next build "
                << "(not bitwise)." << std::endl;
    }
    vra_state = 0;
    vra_pend = 0;
  }
  return ok;
}

}  // namespace radm1
