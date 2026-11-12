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
// the EMA weight of a block boundary for a cell at radius rr:
// a = max(1/N, 1 - exp(-dtb/tau))
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
  const bool dema = vra_ema_dema;
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
    d_(m,t.nf,kk,jj,ii) = v ? 2.0 : 0.0;   // the flag counts the sets (warm-up)
    // the aux slot 1: the current set's raw fields (the change test of the next
    // boundary compares the same set at two times); its flag entry holds the change
    // weight of the coming boundary (-1: none measured)
    const int o1 = t.nf + 1;
    for (int f = 0; f < t.nf; ++f) {
      Real g;
      if (RaGet(a, t, f, ie, m, k, j, i, jr, g)) d_(m,o1+f,kk,jj,ii) = g;
    }
    d_(m,o1+t.nf,kk,jj,ii) = -1.0;
    // raycheap-1109 dema: the second EMA starts equal to the first (2 e1 - e2 = e1)
    if (dema) {
      const int o2 = 2*(t.nf + 1);
      for (int f = 0; f <= t.nf; ++f) {d_(m,o2+f,kk,jj,ii) = d_(m,f,kk,jj,ii);}
    }
  });
  vra_pend = 0;
}

//! raycorr-1108: the change test on the LAST build of a block, before RotAvgDelta turns
//! the Delta slot (vra_k - 1, holding the current set b's raw fields g_b(now)) into the
//! Delta: the same set at the previous boundary is in the aux slot 1 (it was the "next"
//! set there).  c = max(|ln J_b(now)/J_b(then)|, |D_b(now) - D_b(then)|) (J = the
//! kind-2 J field, vet_rotate_avg_ema_change_j; D = field 0, the
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
  // the J field of the test: jf, else the ratio field (kind 2, the set's J; vet_gd
  // registers it as the weight jw only, so jf alone left the J test dead)
  int jt = t.jf;
  for (int f = 0; f < t.nf && jt < 0; ++f) {
    if (t.kind[f] == 2 && t.arr[f] != 2) jt = f;
  }
  par_for("m1_ra_ema_chg", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int kk = k - ks, jj = j - js, ii = i - is;
    Real ach = -1.0;
    if (eps > 0.0) {
      Real c = 0.0;
      if (usej && jt >= 0 && t.arr[jt] != 2) {
        const Real jn = d_(m,od+jt,kk,jj,ii), jo = d_(m,o1+jt,kk,jj,ii);
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
  const bool dema = vra_ema_dema;
  const int o2 = 2*(t.nf + 1);
  par_for("m1_ra_ema_shift", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int kk = k - ks, jj = j - js, ii = i - is;
    const bool vd = d_(m,od+t.nf,kk,jj,ii) > 0.5;
    const bool vh = d_(m,t.nf,kk,jj,ii) > 0.5;
    if (!vd || !vh) {   // no Delta or no history: restart from the current set alone
      for (int f = 0; f < t.nf; ++f) {d_(m,f,kk,jj,ii) = (t.kind[f] == 2) ? 1.0 : 0.0;}
      d_(m,t.nf,kk,jj,ii) = vd ? 1.0 : 0.0;   // count 1: the current set alone
      if (dema) {
        for (int f = 0; f <= t.nf; ++f) {d_(m,o2+f,kk,jj,ii) = d_(m,f,kk,jj,ii);}
      }
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
    // raycheap-1109 dema: the second EMA, of the new first one, in the same frame, with
    // the same (J-) weights: r2' = (1 - aw) r2/R + aw r1', e2' = fac2 (e2 - D) + w1 e1'
    Real fac2 = 1.0 - aw, w1 = aw;
    if (dema && t.jw >= 0) {
      const Real r2 = d_(m,o2+t.jw,kk,jj,ii), R = d_(m,od+t.jw,kk,jj,ii);
      const Real r0 = d_(m,t.jw,kk,jj,ii);
      const Real r1n = (1.0 - aw)*r0/R + aw;
      const Real r2n = (1.0 - aw)*r2/R + aw*r1n;
      fac2 = (1.0 - aw)*(r2/R)/r2n;
      w1 = aw*r1n/r2n;
    }
    for (int f = 0; f < t.nf; ++f) {
      const Real x = d_(m,od+f,kk,jj,ii);
      const Real y = d_(m,f,kk,jj,ii);
      const Real y1 = (t.kind[f] == 2) ? ((1.0 - aw)*y/x + aw) : fac*(y - x);
      d_(m,f,kk,jj,ii) = y1;
      if (dema) {
        const Real y2 = d_(m,o2+f,kk,jj,ii);
        d_(m,o2+f,kk,jj,ii) = (t.kind[f] == 2) ? ((1.0 - aw)*y2/x + aw*y1)
                                                : fac2*(y2 - x) + w1*y1;
      }
    }
    d_(m,t.nf,kk,jj,ii) = (aw >= 1.0) ? 1.0 : fmin(nset + 1.0, 1.0e9);
    if (dema) {d_(m,o2+t.nf,kk,jj,ii) = d_(m,t.nf,kk,jj,ii);}
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
  const bool dema = vra_ema_dema;
  const int o2 = 2*(t.nf + 1);
  par_for("m1_ra_ema_apply", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (!RaValid(a, t, ilo, m, k, j, i)) return;
    const int kk = k - ks, jj = j - js, ii = i - is;
    if (!(d_(m,t.nf,kk,jj,ii) > 0.5)) return;
    const Real jr = RaJ(a, t, m, k, j, i);
    // raycheap-1109 dema: the lag-corrected offset 2 e1 - e2 (ratios 2 r1 - r2 in
    // [r1/4, 4 r1]); without dema the EMA itself
    auto off = [&](const int f) -> Real {
      const Real e1 = d_(m,f,kk,jj,ii);
      if (!dema) return e1;
      const Real e2 = d_(m,o2+f,kk,jj,ii);
      if (t.kind[f] == 2) return fmin(fmax(2.0*e1 - e2, 0.25*e1), 4.0*e1);
      return 2.0*e1 - e2;
    };
    const Real ju = (t.jf >= 0) ? jr*off(t.jf) : jr;
    for (int f = 0; f < t.nf; ++f) {
      Real g;
      if (!RaGet(a, t, f, ie, m, k, j, i, jr, g)) continue;
      Real v;
      if (t.kind[f] == 2) {
        v = g*off(f);
      } else {
        v = g + off(f);
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

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::RotXfSnap
//! \brief vgdxf-1112 (vet_gd_rotate_xfade): the registered fields as they are (no kind-3
//! transform: vet_gd registers kinds 0-2) -> vgd_xf_d slot (0: the current set's raw
//! fields G, before the EMA; 1: the closure A this build would hand over), + validity

void RadiationM1::RotXfSnap(const int slot, const int ilo) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const Tab t = vra_tab;
  const int nf1 = t.nf + 1;
  // (vgd_xf_m > 0: slot 4 = the Delta before the last; lean: slot 5 = the previous
  // build's G, RotXfPreLean)
  const int nsl = (vgd_xf_m > 0) ? (vgd_xf_lean ? 7 : 5) : 4;
  if (vgd_xf_d.extent_int(0) != nmb1 + 1 || vgd_xf_d.extent_int(1) != nsl*nf1 ||
      vgd_xf_d.extent_int(2) != ke - ks + 1 || vgd_xf_d.extent_int(3) != je - js + 1 ||
      vgd_xf_d.extent_int(4) != ie - is + 1) {
    Kokkos::realloc(vgd_xf_d, nmb1 + 1, nsl*nf1, ke - ks + 1, je - js + 1, ie - is + 1);
  }
  auto x_ = vgd_xf_d;
  const RaArr a{vra_a0, vra_a1, vra_c};
  const int o = slot*nf1;
  // (slot 5 / the G times move on only when the time does: a repeated sweep at the
  // same time, e.g. a retry or a second sweep of the build, replaces the latest G)
  const Real tnow = pmy_pack->pmesh->time;
  const bool adv = (slot == 0) && !(tnow == vgd_xf_gt[1] &&
                                    vgd_xf_gb[1] == VetGdBlock(pmy_pack->pmesh->ncycle));
  const bool k5 = adv && (nsl == 7);
  const int o5 = 5*nf1;
  par_for("m1_xf_snap", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int kk = k - ks, jj = j - js, ii = i - is;
    for (int f = 0; f <= t.nf; ++f) {
      if (k5 && (f == t.nf || !(t.arr[f] == 2 && i != ie))) {
        x_(m,o5+f,kk,jj,ii) = x_(m,o+f,kk,jj,ii);
      }
    }
    for (int f = 0; f < t.nf; ++f) {
      if (t.arr[f] == 2 && i != ie) continue;
      x_(m,o+f,kk,jj,ii) = a.At(t.arr[f], t.comp[f], m, k, j, i);
    }
    x_(m,o+t.nf,kk,jj,ii) = RaValid(a, t, ilo, m, k, j, i) ? 1.0 : 0.0;
  });
  if (adv) {
    vgd_xf_gb[0] = vgd_xf_gb[1];
    vgd_xf_gt[0] = vgd_xf_gt[1];
    vgd_xf_gb[1] = VetGdBlock(pmy_pack->pmesh->ncycle);
    vgd_xf_gt[1] = pmy_pack->pmesh->time;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::RotXfBlend
//! \brief vgdxf-1112: the registered fields hold the NEXT set's raw sweep S at this
//! state; vgd_xf_d holds G (slot 0) and A (slot 1).  T = the closure of the first build
//! of the next block at this state: S, or with the EMA history S + fac (e - Delta) (ratio
//! fields S ((1 - aw) r/x + aw)), Delta = S - G (x = S/G), exactly RotAvgEmaShift's
//! update; the change weight is the aux slot's (last = the last build of the block: the
//! RotAvgEmaChange of this step) or measured here as RotAvgEmaChange does (G against the
//! aux slot = this set at the previous boundary).  Fields := (1 - s) A + s T where A and
//! S are valid, else A

void RadiationM1::RotXfBlend(const Real s, const bool last, const int ilo,
                              const Real rx, const bool store) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const Tab t = vra_tab;
  const int nf1 = t.nf + 1;
  auto x_ = vgd_xf_d;
  const RaArr a{vra_a0, vra_a1, vra_c};
  const bool ema = vra_ema && vra_on && vra_state == 1 && vra_d.extent_int(0) > nmb1 &&
                   vra_d.extent_int(1) >= 2*nf1;
  // the K-window (vet_rotate_avg = K): the ring after RotAvgShift of the coming boundary
  const int nslot = vra_k - 1;
  const bool kwin = !vra_ema && vra_on && vra_state == 1 && nslot >= 1 &&
                    vra_d.extent_int(0) > nmb1 && vra_d.extent_int(1) >= nslot*nf1;
  const int hn = (nslot >= 1) ? (vra_head + nslot - 1) % nslot : 0;
  auto d_ = vra_d;
  const Real a0 = (vra_ema_n > 0.0) ? 1.0/vra_ema_n : 1.0;
  const Real tau0 = vra_ema_tau, lc = vra_ema_lc, cl = c_light;
  const Real dtb = vgd_rot_dt;
  auto &size = pmy_pack->pmb->mb_size;
  const int nx1 = indcs.nx1;
  const Real eps = vra_ema_eps;
  const bool usej = vra_ema_chj, step = vra_ema_chstep;
  int jt = t.jf;
  for (int f = 0; f < t.nf && jt < 0; ++f) {
    if (t.kind[f] == 2 && t.arr[f] != 2) jt = f;
  }
  const Real sb = fmin(fmax(s, 0.0), 1.0);
  // the boundary Delta at the boundary time t_k (rx = (t_k - t)/(t - t_prev) >= 0; < 0:
  // no previous build of this block): the last build's Delta (state t_k - dt, an O(dt)
  // offset carried into the whole next block) extrapolated with the previous build's
  const int od = (vra_k - 1)*nf1;
  const bool dext = last && rx >= 0.0 && (ema || kwin) && vra_d.extent_int(1) >= od + nf1;
  // vgd_xf_m > 0: slot 4 := the Delta of the previous live sweep (RotXfEst)
  const bool keep4 = store && vgd_xf_d.extent_int(1) >= 5*nf1;
  const int o4 = 4*nf1;
  par_for("m1_xf_blend", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int kk = k - ks, jj = j - js, ii = i - is;
    const bool vs = RaValid(a, t, ilo, m, k, j, i);
    const bool vg = x_(m,t.nf,kk,jj,ii) > 0.5;
    const bool va = x_(m,nf1+t.nf,kk,jj,ii) > 0.5;
    {
      const int o2 = 2*nf1;
      const bool vn = vs && vg;
      const bool vp = x_(m,o2+t.nf,kk,jj,ii) > 0.5;
      const bool ext = dext && vn && vp && (d_(m,od+t.nf,kk,jj,ii) > 0.5);
      for (int f = 0; f < t.nf; ++f) {
        if (t.arr[f] == 2 && i != ie) continue;
        const Real sv = a.At(t.arr[f], t.comp[f], m, k, j, i);
        const Real gv = x_(m,f,kk,jj,ii);
        const bool mul = (t.kind[f] == 2);
        const Real dn = mul ? ((gv > 0.0 && sv > 0.0) ? sv/gv : 1.0) : (sv - gv);
        if (ext) {
          const Real de = dn + rx*(dn - x_(m,o2+f,kk,jj,ii));
          d_(m,od+f,kk,jj,ii) = (mul && !(de > 0.0)) ? dn : de;
        }
        if (store) {
          if (keep4) {x_(m,o4+f,kk,jj,ii) = x_(m,o2+f,kk,jj,ii);}
          x_(m,o2+f,kk,jj,ii) = dn;
          x_(m,o2+nf1+f,kk,jj,ii) = gv;   // slot 3: G of this build (RotXfPre)
        }
      }
      if (store) {
        if (keep4) {x_(m,o4+t.nf,kk,jj,ii) = x_(m,o2+t.nf,kk,jj,ii);}
        x_(m,o2+t.nf,kk,jj,ii) = vn ? 1.0 : 0.0;
        x_(m,o2+nf1+t.nf,kk,jj,ii) = vg ? 1.0 : 0.0;
      }
    }
    if (!vs || !va) {   // no blend: A as it was
      for (int f = 0; f < t.nf; ++f) {
        if (t.arr[f] == 2 && i != ie) continue;
        a.At(t.arr[f], t.comp[f], m, k, j, i) = x_(m,nf1+f,kk,jj,ii);
      }
      return;
    }
    // the EMA shift of the coming boundary (RotAvgEmaShift), with Delta from this state
    bool sh = ema && vg && (d_(m,t.nf,kk,jj,ii) > 0.5);
    Real aw = 1.0, fac = 0.0;
    if (sh) {
      Real ach = -1.0;
      if (last) {
        ach = d_(m,nf1+t.nf,kk,jj,ii);
      } else if (eps > 0.0) {
        Real c = 0.0;
        if (usej && jt >= 0 && t.arr[jt] != 2) {
          const Real jn = x_(m,jt,kk,jj,ii), jo = d_(m,nf1+jt,kk,jj,ii);
          c = (jn > 0.0 && jo > 0.0) ? fabs(log(jn/jo)) : 1.0e30;
        }
        if (t.nf > 0 && t.arr[0] != 2 && t.kind[0] != 2) {
          c = fmax(c, fabs(x_(m,0,kk,jj,ii) - d_(m,nf1,kk,jj,ii)));
        }
        ach = step ? fmin(1.0, fmax(0.0, c/eps - 1.0)) : fmin(1.0, c/eps);
      }
      Real rr = 1.0;
      if (lc > 0.0) {
        const Real x1min = size.d_view(m).x1min, x1max = size.d_view(m).x1max;
        rr = x1min + (static_cast<Real>(ii) + 0.5)*(x1max - x1min)/nx1;
      }
      const Real nset = d_(m,t.nf,kk,jj,ii);
      aw = fmax(fmax(RaEmaA(a0, dtb, tau0, lc, rr, cl), fmin(ach, 1.0)),
                1.0/(nset + 1.0));
      fac = 1.0 - aw;
      if (t.jw >= 0) {
        const Real sj = a.At(t.arr[t.jw], t.comp[t.jw], m, k, j, i);
        const Real gj = x_(m,t.jw,kk,jj,ii);
        const Real R = (gj > 0.0 && sj > 0.0) ? sj/gj : 1.0;
        const Real r0 = d_(m,t.jw,kk,jj,ii);
        const Real r1 = (1.0 - aw)*r0/R + aw;
        fac = (1.0 - aw)*(r0/R)/r1;
      }
    }
    // K-window: the shifted offsets d'_s (slot hn: -Delta, the others d_s - Delta; ratio
    // fields 1/x, d_s/x), valid where Delta and the slot are; then RotAvgApply's mean
    const bool kw = kwin && vg;
    Real kws[M1_RA_NF], kwsum = 1.0, knv = 1.0;
    if (kw) {
      Real xj = 1.0;
      if (t.jw >= 0) {
        const Real sj = a.At(t.arr[t.jw], t.comp[t.jw], m, k, j, i);
        const Real gj = x_(m,t.jw,kk,jj,ii);
        xj = (gj > 0.0 && sj > 0.0) ? sj/gj : 1.0;
      }
      for (int q = 0; q < nslot; ++q) {
        const bool vq = (q == hn) || (d_(m,q*nf1+t.nf,kk,jj,ii) > 0.5);
        knv += vq ? 1.0 : 0.0;
        kws[q] = vq ? ((t.jw >= 0) ? ((q == hn) ? 1.0/xj : d_(m,q*nf1+t.jw,kk,jj,ii)/xj)
                                   : 1.0) : 0.0;
        kwsum += kws[q];
      }
    }
    for (int f = 0; f < t.nf; ++f) {
      if (t.arr[f] == 2 && i != ie) continue;
      Real &fv = a.At(t.arr[f], t.comp[f], m, k, j, i);
      const Real sv = fv;
      Real tv = sv;
      if (kw && knv > 1.5) {
        const Real gv = x_(m,f,kk,jj,ii);
        const bool mul = (t.kind[f] == 2);
        const Real x = mul ? ((gv > 0.0 && sv > 0.0) ? sv/gv : 1.0) : (sv - gv);
        Real sm = mul ? 1.0 : 0.0;
        for (int q = 0; q < nslot; ++q) {
          if (!((q == hn) || d_(m,q*nf1+t.nf,kk,jj,ii) > 0.5)) continue;
          const Real dq = (q == hn) ? (mul ? 1.0/x : -x)
                                    : (mul ? d_(m,q*nf1+f,kk,jj,ii)/x
                                           : d_(m,q*nf1+f,kk,jj,ii) - x);
          sm += mul ? dq : kws[q]*dq;
        }
        if (mul) {
          tv = sv*sm/knv;
        } else {
          tv = sv + sm/kwsum;
          if (t.kind[f] == 1) tv = fmin(fmax(tv, t.lo[f]), t.hi[f]);
        }
      } else if (sh) {
        const Real gv = x_(m,f,kk,jj,ii);
        const Real e = d_(m,f,kk,jj,ii);
        if (t.kind[f] == 2) {
          const Real xr = (gv > 0.0 && sv > 0.0) ? sv/gv : 1.0;
          tv = sv*((1.0 - aw)*e/xr + aw);
        } else {
          tv = sv + fac*(e - (sv - gv));
          if (t.kind[f] == 1) tv = fmin(fmax(tv, t.lo[f]), t.hi[f]);
        }
      }
      fv = (1.0 - sb)*x_(m,nf1+f,kk,jj,ii) + sb*tv;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::RotXfPre
//! \brief vgdxf-1112, the LAST build of a block, before RotAvgDelta: the boundary
//! quantities at the boundary time t_k instead of this build's t_k - dt (an O(dt) offset
//! the EMA / K-window would carry into the whole next block): the current set's raw
//! fields (the Delta slot and G) and the next set's sweep (the fields) extrapolated
//! linearly with the previous build of this block, rx = (t_k - t)/(t - t_prev).  The
//! EMA change test, the aux slot and the Delta then all refer to t_k; the fields are
//! put back to the live sweep afterwards (VetGdBuild) and G (slot 0) is not changed

void RadiationM1::RotXfPre(const Real rx, const int ilo) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const Tab t = vra_tab;
  const int nf1 = t.nf + 1;
  if (!(rx >= 0.0) || vgd_xf_d.extent_int(1) < 4*nf1 ||
      vra_d.extent_int(1) < vra_k*nf1) {return;}
  auto x_ = vgd_xf_d;
  auto d_ = vra_d;
  const RaArr a{vra_a0, vra_a1, vra_c};
  const int od = (vra_k - 1)*nf1, o2 = 2*nf1, o3 = 3*nf1;
  par_for("m1_xf_pre", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int kk = k - ks, jj = j - js, ii = i - is;
    if (!RaValid(a, t, ilo, m, k, j, i) || !(x_(m,t.nf,kk,jj,ii) > 0.5) ||
        !(x_(m,o2+t.nf,kk,jj,ii) > 0.5) || !(x_(m,o3+t.nf,kk,jj,ii) > 0.5) ||
        !(d_(m,od+t.nf,kk,jj,ii) > 0.5)) {return;}
    for (int f = 0; f < t.nf; ++f) {
      if (t.arr[f] == 2 && i != ie) continue;
      const bool mul = (t.kind[f] == 2);
      Real &fv = a.At(t.arr[f], t.comp[f], m, k, j, i);
      const Real sv = fv, gv = x_(m,f,kk,jj,ii);
      const Real gp = x_(m,o3+f,kk,jj,ii), dp = x_(m,o2+f,kk,jj,ii);
      const Real sp = mul ? gp*dp : gp + dp;
      Real gk = gv + rx*(gv - gp), sk = sv + rx*(sv - sp);
      if (mul) {
        if (!(gk > 0.0)) gk = gv;
        if (!(sk > 0.0)) sk = sv;
      } else if (t.kind[f] == 1) {
        gk = fmin(fmax(gk, t.lo[f]), t.hi[f]);
        sk = fmin(fmax(sk, t.lo[f]), t.hi[f]);
      }
      d_(m,od+f,kk,jj,ii) = gk;
      fv = sk;   // put back after RotAvgDelta (VetGdBuild): the blend uses the live S
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::RotXfPreLean
//! \brief vgdxf-1112 vet_gd_rotate_xfade_lean: RotXfPre without a live next-set sweep at
//! the previous build: G_k = G + rg (G - G_prev build) (slot 5, rg = (t_k - t)/(t -
//! t_prev)); S_k = G_k (+) (D + rs (P - L)), D = the live sweep's Delta at this state,
//! (P - L) the estimate line's rise (RotXfSample; rs = (t_k - t)/its length; < 0: D held)

void RadiationM1::RotXfPreLean(const Real rg, const Real rs, const int ilo) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const Tab t = vra_tab;
  const int nf1 = t.nf + 1;
  if (!(rg >= 0.0) || vgd_xf_d.extent_int(1) < 7*nf1 ||
      vra_d.extent_int(1) < vra_k*nf1) {return;}
  auto x_ = vgd_xf_d;
  auto d_ = vra_d;
  const RaArr a{vra_a0, vra_a1, vra_c};
  const int od = (vra_k - 1)*nf1, o4 = 4*nf1, o5 = 5*nf1, o6 = 6*nf1;
  const bool sl = (rs >= 0.0);
  par_for("m1_xf_prel", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int kk = k - ks, jj = j - js, ii = i - is;
    if (!RaValid(a, t, ilo, m, k, j, i) || !(x_(m,t.nf,kk,jj,ii) > 0.5) ||
        !(x_(m,o5+t.nf,kk,jj,ii) > 0.5) || !(d_(m,od+t.nf,kk,jj,ii) > 0.5)) {return;}
    const bool ds = sl && (x_(m,o4+t.nf,kk,jj,ii) > 0.5);
    for (int f = 0; f < t.nf; ++f) {
      if (t.arr[f] == 2 && i != ie) continue;
      const bool mul = (t.kind[f] == 2);
      Real &fv = a.At(t.arr[f], t.comp[f], m, k, j, i);
      const Real sv = fv, gv = x_(m,f,kk,jj,ii), gp = x_(m,o5+f,kk,jj,ii);
      Real gk = gv + rg*(gv - gp);
      if (mul && !(gk > 0.0)) gk = gv;
      const Real dn = mul ? ((gv > 0.0 && sv > 0.0) ? sv/gv : 1.0) : (sv - gv);
      Real dk = dn;
      if (ds) {
        dk = dn + rs*(x_(m,o6+f,kk,jj,ii) - x_(m,o4+f,kk,jj,ii));
        if (mul && !(dk > 0.0)) dk = dn;
      }
      Real sk = mul ? gk*dk : gk + dk;
      if (t.kind[f] == 1) {
        gk = fmin(fmax(gk, t.lo[f]), t.hi[f]);
        sk = fmin(fmax(sk, t.lo[f]), t.hi[f]);
      }
      d_(m,od+f,kk,jj,ii) = gk;
      fv = sk;   // put back after RotAvgDelta (VetGdBuild)
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::RotXfKeep
//! \brief vgdxf-1112 lean + vet_gd_rotate_xfade_pen: the live next-set sweep of the
//! second-to-last build only feeds RotXfPre: slot 2 := its Delta, slot 3 := G (+ flags);
//! (the line slots 4, 6 are untouched)

void RadiationM1::RotXfKeep(const int ilo) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const Tab t = vra_tab;
  const int nf1 = t.nf + 1;
  auto x_ = vgd_xf_d;
  const RaArr a{vra_a0, vra_a1, vra_c};
  const int o2 = 2*nf1, o3 = 3*nf1;
  par_for("m1_xf_keep", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int kk = k - ks, jj = j - js, ii = i - is;
    const bool vg = x_(m,t.nf,kk,jj,ii) > 0.5;
    const bool vn = RaValid(a, t, ilo, m, k, j, i) && vg;
    for (int f = 0; f < t.nf; ++f) {
      if (t.arr[f] == 2 && i != ie) continue;
      const Real sv = a.At(t.arr[f], t.comp[f], m, k, j, i), gv = x_(m,f,kk,jj,ii);
      x_(m,o2+f,kk,jj,ii) = (t.kind[f] == 2) ? ((gv > 0.0 && sv > 0.0) ? sv/gv : 1.0)
                                             : (sv - gv);
      x_(m,o3+f,kk,jj,ii) = gv;
    }
    x_(m,o2+t.nf,kk,jj,ii) = vn ? 1.0 : 0.0;
    x_(m,o3+t.nf,kk,jj,ii) = vg ? 1.0 : 0.0;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::RotXfSample
//! \brief vgdxf-1112 vet_gd_rotate_xfade_lean, a build with a live next-set sweep (the
//! registered fields): D = its Delta (slot 2).  The estimate stays continuous in time:
//! the new line starts at the old line's value L = L' + q (P' - L') (first sample of a
//! block: L = D) and aims at P = 2 D - D_prev (first: D) at the end of this part (slots
//! 4, 6; ratio fields P <= 0: D).  The fields := G (+) L, the estimate at this time

void RadiationM1::RotXfSample(const Real q, const int ilo) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const Tab t = vra_tab;
  const int nf1 = t.nf + 1;
  if (vgd_xf_d.extent_int(1) < 7*nf1) {return;}
  auto x_ = vgd_xf_d;
  const RaArr a{vra_a0, vra_a1, vra_c};
  const int o2 = 2*nf1, o4 = 4*nf1, o6 = 6*nf1;
  const bool first = !(q >= 0.0);
  par_for("m1_xf_samp", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int kk = k - ks, jj = j - js, ii = i - is;
    const bool vg = x_(m,t.nf,kk,jj,ii) > 0.5;
    const bool vn = RaValid(a, t, ilo, m, k, j, i) && vg;
    const bool vo = !first && (x_(m,o2+t.nf,kk,jj,ii) > 0.5) &&
                    (x_(m,o4+t.nf,kk,jj,ii) > 0.5);
    for (int f = 0; f < t.nf; ++f) {
      if (t.arr[f] == 2 && i != ie) continue;
      const bool mul = (t.kind[f] == 2);
      Real &fv = a.At(t.arr[f], t.comp[f], m, k, j, i);
      const Real sv = fv, gv = x_(m,f,kk,jj,ii);
      const Real dn = mul ? ((gv > 0.0 && sv > 0.0) ? sv/gv : 1.0) : (sv - gv);
      Real lv = dn, pv = dn;
      if (vo) {
        const Real lo = x_(m,o4+f,kk,jj,ii);
        lv = lo + q*(x_(m,o6+f,kk,jj,ii) - lo);
        pv = 2.0*dn - x_(m,o2+f,kk,jj,ii);
        if (mul) {
          if (!(lv > 0.0)) lv = dn;
          if (!(pv > 0.0)) pv = dn;
        }
      }
      x_(m,o2+f,kk,jj,ii) = dn;
      x_(m,o4+f,kk,jj,ii) = lv;
      x_(m,o6+f,kk,jj,ii) = pv;
      Real v = mul ? gv*lv : gv + lv;
      if (t.kind[f] == 1) v = fmin(fmax(v, t.lo[f]), t.hi[f]);
      fv = v;
    }
    x_(m,o2+t.nf,kk,jj,ii) = vn ? 1.0 : 0.0;
    x_(m,o4+t.nf,kk,jj,ii) = (vn && (first || vo)) ? 1.0 : 0.0;
    if (!vn && t.farr >= 0) {a.At(t.farr, t.fcomp, m, k, j, i) = -1.0;}
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::RotXfEst
//! \brief vgdxf-1112 vet_gd_rotate_xfade_m = M > 0 (the amortized cross-fade): a build
//! without a next-set sweep: the fields := the next set's estimate G (+) Delta^ (ratio
//! fields G Delta^), Delta^ = D2 + r (D2 - D4): the last two live sweeps' same-state
//! differences (slots 2, 4) extrapolated linearly in time, r = (t - t_2)/(t_2 - t_4)
//! (r < 0 or slot 4 invalid: D2 held; a ratio extrapolation <= 0: D2); where D2 is
//! missing the flag field is set invalid (RotXfBlend then hands over A)

void RadiationM1::RotXfEst(const int ilo, const Real r) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const Tab t = vra_tab;
  const int nf1 = t.nf + 1;
  auto x_ = vgd_xf_d;
  const RaArr a{vra_a0, vra_a1, vra_c};
  const int o2 = 2*nf1, o4 = 4*nf1, o6 = 6*nf1;
  // lean: Delta^ = L + r (P - L) (slots 4, 6; r = the line's fraction, RotXfSample)
  const bool ln = vgd_xf_lean && vgd_xf_d.extent_int(1) >= 7*nf1;
  const bool lin = !ln && (r >= 0.0) && vgd_xf_d.extent_int(1) >= 5*nf1;
  par_for("m1_xf_est", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int kk = k - ks, jj = j - js, ii = i - is;
    const bool v = (i >= ilo) && (x_(m,t.nf,kk,jj,ii) > 0.5) &&
                   (x_(m,(ln ? o4 : o2)+t.nf,kk,jj,ii) > 0.5);
    const bool l2 = lin && (x_(m,o4+t.nf,kk,jj,ii) > 0.5);
    for (int f = 0; f < t.nf; ++f) {
      if (t.arr[f] == 2 && i != ie) continue;
      const Real gv = x_(m,f,kk,jj,ii), d2 = x_(m,o2+f,kk,jj,ii);
      Real dv = d2;
      if (ln) {
        const Real lv = x_(m,o4+f,kk,jj,ii);
        dv = lv + r*(x_(m,o6+f,kk,jj,ii) - lv);
        if (t.kind[f] == 2 && !(dv > 0.0)) dv = lv;
      } else if (l2) {
        dv = d2 + r*(d2 - x_(m,o4+f,kk,jj,ii));
        if (t.kind[f] == 2 && !(dv > 0.0)) dv = d2;
      }
      Real v = (t.kind[f] == 2) ? gv*dv : gv + dv;
      if (t.kind[f] == 1) v = fmin(fmax(v, t.lo[f]), t.hi[f]);
      a.At(t.arr[f], t.comp[f], m, k, j, i) = v;
    }
    if (!v && t.farr >= 0) {a.At(t.farr, t.fcomp, m, k, j, i) = -1.0;}
  });
}

}  // namespace radm1
