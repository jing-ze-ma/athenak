//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_cart_rotavg.cpp
//! \brief vet_rotate_avg = K >= 2: the closure averaged over the last K rotated direction
//! sets (K-window or EMA), and its restart record.
//!
//! With c dt/dx >> 1 each rotated set relaxes to its own steady state, so plain rotation
//! does not average the ray effects.  g_s = normalised hand-over fields of set s (ratios
//! to J, and J).  The moment step gets
//!   used = g_b(now) + (1/K) sum_{j=1}^{K-1} d_j      (additive fields)
//!   used = g_b(now) (1 + sum_j r_j)/K                (J: ratio fields)
//! with offsets d_j = g_{b-j} - g_b (r_j = J_{b-j}/J_b).  Each offset is a difference of
//! two sweeps of the SAME state: at a block change the old set is swept once more,
//! Delta = g_new - g_old, and the ring shifts d_1' = -Delta, d_j' = d_{j-1} - Delta
//! (ratios: r_1' = 1/R, r_j' = r_{j-1}/R).  Steady state: used = (1/K) sum_s g_s exactly
//! (a realizable mixture); the current set stays live, only pattern differences age.
//! Kinds: 0 additive, 1 additive clamped to [lo, hi], 2 ratio (J), 3 raw moment averaged
//! as X/J and handed back times the mean J.  A cell is valid in a set when i >= ilo, the
//! flag field is >= 0 and J > 0; offsets count only where both sweeps were valid.
//! Memory: (K-1)(nf+1) reals per active cell.
//!
//!  called from: rays/sc_cart_sweep.cpp, outputs/restart.cpp, pgen/pgen.cpp

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <utility>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

namespace {

// field arrays of one call (device copies)
struct RaArr {
  DvceArray5D<Real> a0, a1;
  DvceArray3D<Real> c;
  KOKKOS_INLINE_FUNCTION
  Real &At(const int arr, const int comp, const int m, const int k, const int j,
           const int i) const {
    return (arr == 2) ? c(m,k,j) : ((arr == 0) ? a0(m,comp,k,j,i) : a1(m,comp,k,j,i));
  }
};

using Tab = RadiationVet::RotAvgTab;

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
//! \fn void RadiationVet::RotAvgSnap
//! \brief registered fields as they are (kind 3 as X/J) -> ring slot `slot`, with the
//! validity of each cell

void RadiationVet::RotAvgSnap(const int slot, const int ilo) {
  vra_gen += 1;   // vra_d may change (dt_retry_host skips clean saves)
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nf = vra_tab.nf;
  const int nalloc = vra_k;   // ring (K-1 slots) + Delta slot K-1
  if (vra_d.extent_int(0) < nmb1 + 1 || vra_d.extent_int(1) != nalloc*(nf + 1)) {
    Kokkos::realloc(vra_d, nmb1 + 1, nalloc*(nf + 1), ke - ks + 1, je - js + 1,
                    ie - is + 1);
    Kokkos::deep_copy(vra_d, 0.0);
  }
  auto d_ = vra_d;
  const RaArr a{vra_a0, vra_a1, vra_c};
  const Tab t = vra_tab;
  const int o = slot*(nf + 1);
  par_for("vet_ra_snap", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
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
//! \fn void RadiationVet::RotAvgUpdate
//! \brief first build of a run: with the current set's raw fields in place, the
//! snapshots S_j of the K-1 previous sets -> offsets S_j - g (ratio fields S_j/g)

void RadiationVet::RotAvgUpdate(const int ilo) {
  if (vra_pend != 1) return;
  if (vra_ema) {
    RotAvgEmaInit(ilo);
    return;
  }
  vra_gen += 1;   // vra_d may change (dt_retry_host skips clean saves)
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
  par_for("vet_ra_upd", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
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
          d_(m,o+f,kk,jj,ii) = (t.kind[f] == 2)
                               ? ((v && g > 0.0 && sv > 0.0) ? sv/g : 1.0)
                               : (v ? (sv - g) : 0.0);
        }
        d_(m,o+t.nf,kk,jj,ii) = v ? 1.0 : 0.0;
      }
    }
  });
  vra_pend = 0;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RotAvgDelta
//! \brief last build of a block: slot holds g_b, the fields the NEXT set's S at the same
//! state -> slot := Delta = S - g_b (ratio: S/g_b), valid where both are.  The last build
//! samples the block boundary to O(dt^2) on every dt.

void RadiationVet::RotAvgDelta(const int slot, const int ilo) {
  if (vra_ema) {RotAvgEmaChange(ilo);}
  vra_gen += 1;   // vra_d may change (dt_retry_host skips clean saves)
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto d_ = vra_d;
  const RaArr a{vra_a0, vra_a1, vra_c};
  const Tab t = vra_tab;
  const int o = slot*(t.nf + 1);
  par_for("vet_ra_delta", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
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
//! \fn void RadiationVet::RotAvgShift
//! \brief first build of a new block: the ring moves on by one set with the Delta of
//! slot ds: the dropped oldest slot becomes d_1 = -Delta (1/R), every other d_j -= Delta
//! (/= R)

void RadiationVet::RotAvgShift(const int ds) {
  if (vra_ema) {
    RotAvgEmaShift(ds);
    return;
  }
  vra_gen += 1;   // vra_d may change (dt_retry_host skips clean saves)
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
  par_for("vet_ra_shift", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
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
//! \fn void RadiationVet::RotAvgApply
//! \brief registered fields := mean over the current set and the valid offsets (ratio
//! fields plain, others J-weighted by jw; kind 3: (X/J)_avg times the mean J)

void RadiationVet::RotAvgApply(const int ilo) {
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
  par_for("vet_ra_apply", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (!RaValid(a, t, ilo, m, k, j, i)) return;
    const int kk = k - ks, jj = j - js, ii = i - is;
    Real nv = 1.0;
    for (int s = 0; s < nslot; ++s) {nv += d_(m,s*nf1+t.nf,kk,jj,ii);}
    if (nv < 1.5) return;
    const Real jr = RaJ(a, t, m, k, j, i);
    Real ju = jr;   // averaged J (kind 3 fields are handed back with it)
    if (t.jf >= 0) {
      Real sm = 1.0;
      for (int s = 0; s < nslot; ++s) {
        const int o = s*nf1;
        if (d_(m,o+t.nf,kk,jj,ii) > 0.5) sm += d_(m,o+t.jf,kk,jj,ii);
      }
      ju = jr*sm/nv;
    }
    // weight of each older set: its J over the current one's (jw), so the fields are the
    // normalised moments of the union of the sets (1 without jw or with weight = equal)
    Real ws[VET_RA_NF], wsum = 1.0;
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
//! vet_rotate_avg_ema = N: exponential moving average in the offset form (one ring slot,
//! vra_k = 2), weights a (1 - a)^j over the sets b - j (J-weighted).  Slot 0 holds
//!   e_f = avg_f - g_f(now) (additive),  r_f = avg_f/g_f(now) (ratio: J)
//! and the boundary Delta (D_f, R_f) updates them as
//!   r' = (1 - a) r/R + a,   e' = (1 - a) (r_jw/R_jw) (e - D)/r'_jw.
//! Init: the equal mix of the current and the previous set (as K = 2).

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RotAvgEmaInit
//! \brief EMA init on the first build (equal mix of two sets) and the aux slot

void RadiationVet::RotAvgEmaInit(const int ilo) {
  vra_gen += 1;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto d_ = vra_d;
  const RaArr a{vra_a0, vra_a1, vra_c};
  const Tab t = vra_tab;
  par_for("vet_ra_ema_init", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int kk = k - ks, jj = j - js, ii = i - is;
    const bool v = RaValid(a, t, ilo, m, k, j, i) && (d_(m,t.nf,kk,jj,ii) > 0.5);
    const Real jr = RaJ(a, t, m, k, j, i);
    // previous set's weight relative to the current one (its J ratio, or 1)
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
    d_(m,t.nf,kk,jj,ii) = v ? 2.0 : 0.0;   // flag counts the sets (warm-up)
    // aux slot 1: the current set's raw fields for the next boundary's change test; its
    // flag entry holds that boundary's change weight (-1: none measured)
    const int o1 = t.nf + 1;
    for (int f = 0; f < t.nf; ++f) {
      Real g;
      if (RaGet(a, t, f, ie, m, k, j, i, jr, g)) d_(m,o1+f,kk,jj,ii) = g;
    }
    d_(m,o1+t.nf,kk,jj,ii) = -1.0;
  });
  vra_pend = 0;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RotAvgEmaChange
//! \brief change test on the last build of a block, before RotAvgDelta: compares set b
//! now with the same set at the previous boundary (aux slot 1).
//! c = max(|ln J_b(now)/J_b(then)|, |D_b(now) - D_b(then)|) (vet_rotate_avg_ema_change_j,
//! field 0 = D_rr) gives a_ch = min(1, c/eps) (vet_rotate_avg_ema_change = eps, 0 off),
//! stored in the aux flag; then aux := the next set's raw fields.  A changed geometry
//! thus drops the old history.
void RadiationVet::RotAvgEmaChange(const int ilo) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto d_ = vra_d;
  const RaArr a{vra_a0, vra_a1, vra_c};
  const Tab t = vra_tab;
  const int o1 = t.nf + 1, od = (vra_k - 1)*(t.nf + 1);
  const Real eps = kEmaEps;
  const bool step = vra_ema_chstep;
  // J field of the test: jf, else the kind-2 ratio field (vet_gd registers J as jw only)
  int jt = t.jf;
  for (int f = 0; f < t.nf && jt < 0; ++f) {
    if (t.kind[f] == 2 && t.arr[f] != 2) jt = f;
  }
  par_for("vet_ra_ema_chg", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int kk = k - ks, jj = j - js, ii = i - is;
    Real c = 0.0;
    if (jt >= 0 && t.arr[jt] != 2) {
      const Real jn = d_(m,od+jt,kk,jj,ii), jo = d_(m,o1+jt,kk,jj,ii);
      c = (jn > 0.0 && jo > 0.0) ? fabs(log(jn/jo)) : 1.0e30;
    }
    if (t.nf > 0 && t.arr[0] != 2 && t.kind[0] != 2) {
      c = fmax(c, fabs(d_(m,od,kk,jj,ii) - d_(m,o1,kk,jj,ii)));
    }
    // linear: a_ch = c/eps; step: 0 below eps, 1 at 2 eps
    const Real ach = step ? fmin(1.0, fmax(0.0, c/eps - 1.0)) : fmin(1.0, c/eps);
    const Real jr = RaJ(a, t, m, k, j, i);
    for (int f = 0; f < t.nf; ++f) {
      Real g;
      if (RaGet(a, t, f, ie, m, k, j, i, jr, g)) d_(m,o1+f,kk,jj,ii) = g;
    }
    d_(m,o1+t.nf,kk,jj,ii) = ach;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RotAvgEmaShift
//! \brief EMA update at a block boundary (first build of the new block)

void RadiationVet::RotAvgEmaShift(const int ds) {
  if (vra_ema_wb) {
    std::cout << "### FATAL ERROR <rad_vet>/vet_rotate_avg_ema_wb needs the clean ring "
              << "(ring = live)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  vra_gen += 1;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto d_ = vra_d;
  const Tab t = vra_tab;
  const int od = ds*(t.nf + 1);
  const Real a0 = 1.0/vra_ema_n;
  par_for("vet_ra_ema_shift", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int kk = k - ks, jj = j - js, ii = i - is;
    const bool vd = d_(m,od+t.nf,kk,jj,ii) > 0.5;
    const bool vh = d_(m,t.nf,kk,jj,ii) > 0.5;
    if (!vd || !vh) {   // no Delta or no history: restart from the current set alone
      for (int f = 0; f < t.nf; ++f) {d_(m,f,kk,jj,ii) = (t.kind[f] == 2) ? 1.0 : 0.0;}
      d_(m,t.nf,kk,jj,ii) = vd ? 1.0 : 0.0;   // count 1: current set alone
      return;
    }
    const Real ach = d_(m,(t.nf+1)+t.nf,kk,jj,ii);   // change weight (aux slot)
    // warm-up: with n sets seen the weight is >= 1/(n+1) (plain mean until n = N)
    const Real nset = d_(m,t.nf,kk,jj,ii);
    const Real aw = fmax(fmax(fmin(a0, 1.0), fmin(ach, 1.0)), 1.0/(nset + 1.0));
    Real fac = 1.0 - aw;
    if (t.jw >= 0) {
      const Real r0 = d_(m,t.jw,kk,jj,ii), R = d_(m,od+t.jw,kk,jj,ii);
      const Real r1 = (1.0 - aw)*r0/R + aw;
      fac = (1.0 - aw)*(r0/R)/r1;
    }
    for (int f = 0; f < t.nf; ++f) {
      const Real x = d_(m,od+f,kk,jj,ii);
      const Real y = d_(m,f,kk,jj,ii);
      const Real y1 = (t.kind[f] == 2) ? ((1.0 - aw)*y/x + aw) : fac*(y - x);
      d_(m,f,kk,jj,ii) = y1;
    }
    d_(m,t.nf,kk,jj,ii) = (aw >= 1.0) ? 1.0 : fmin(nset + 1.0, 1.0e9);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RotAvgEmaApply
//! \brief registered fields := current set + EMA offset

void RadiationVet::RotAvgEmaApply(const int ilo) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto d_ = vra_d;
  const RaArr a{vra_a0, vra_a1, vra_c};
  const Tab t = vra_tab;
  par_for("vet_ra_ema_apply", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (!RaValid(a, t, ilo, m, k, j, i)) return;
    const int kk = k - ks, jj = j - js, ii = i - is;
    if (!(d_(m,t.nf,kk,jj,ii) > 0.5)) return;
    const Real jr = RaJ(a, t, m, k, j, i);
    auto off = [&](const int f) -> Real {return d_(m,f,kk,jj,ii);};
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
//! Restart record of the averaging state (kVetRaRstMagic), for a bitwise continuation.
//! Header (int64 x VET_RA_RST_NH): have = 1, nslab, K, nf, closure (1 vet_gd, 2 vet_sc,
//! +16 EMA), head, pend, nslot (K-1, or K with the Delta slot when pend = 3), blk,
//! interior cells per MeshBlock, i0.  The slots' interior values for i >= i0 (first i
//! valid in any slot, over all ranks) are packed densely into nslab slabs at the end of
//! each MeshBlock record, so the history moves with its block.  Written only with a
//! complete history (vra_state = 1); cells below i0 load as zero.  A 10-entry header
//! reads i0 = 0.

// registered field count: the table's, or the file's (restart before the first build)
int RadiationVet::RotAvgNf() const {
  return (vra_tab.nf > 0) ? vra_tab.nf : vra_rst_nf;
}

// first interior i where a written slot is valid, over all MeshBlocks and ranks
// (at most nx1 - 1)
int RadiationVet::RotAvgRstI0() const {
  const int nf = RotAvgNf();
  const int nslot = (vra_pend == 3) ? vra_k : vra_k - 1;
  auto d_ = vra_d;
  const int nmb = vra_d.extent_int(0), nk = vra_d.extent_int(2);
  const int nj = vra_d.extent_int(3), ni = vra_d.extent_int(4);
  int i0 = ni - 1;
  if (nmb > 0 && nf > 0) {
    const int nf1 = nf + 1;
    Kokkos::parallel_reduce("vet_ra_rst_i0",
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

// slabs of the restart record (0: nothing to write)
int RadiationVet::RotAvgRstNslab(const int ncell) const {
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

// restart header h[0..10]
void RadiationVet::RotAvgRstHdr(int64_t *h) const {
  h[10] = 0;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  h[0] = 1;
  h[2] = vra_k;
  h[3] = RotAvgNf();
  h[4] = (vgd_on ? 1 : 2) + (vra_ema ? 16 : 0);   // EMA history
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

// flat copy between vra_d's first nslot (nf+1) components (interior, i >= i0) and the
// dense slabs (dir 0: vra_d -> slabs, zero-padded; 1: slabs -> vra_d)
void RaRstCopy(const DvceArray5D<Real> &d, const DvceArray5D<Real> &sl, const int nmb,
               const int ncomp, const int dir, const int i0) {
  const int ni = d.extent_int(4) - i0, nj = d.extent_int(3), nk = d.extent_int(2);
  const int n1 = sl.extent_int(4), n2 = sl.extent_int(3), n3 = sl.extent_int(2);
  const int64_t ncint = static_cast<int64_t>(ni)*nj*nk;
  const int64_t ncout = static_cast<int64_t>(n1)*n2*n3;
  const int64_t nv = ncomp*ncint;
  const int64_t ntot = static_cast<int64_t>(sl.extent_int(1))*ncout;
  // one MeshBlock per launch, at most 2^30 values each (int index range)
  const int64_t nch = int64_t(1) << 30;
  for (int m = 0; m < nmb; ++m) {
  for (int64_t q0 = 0; q0 < ntot; q0 += nch) {
  const int64_t nq = std::min(nch, ntot - q0);
  par_for("vet_ra_rst", DevExeSpace(), 0, static_cast<int>(nq - 1),
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

//! the same record straight into the host output array (nmb, nslab, n3, n2, n1) through
//! a device buffer of at most nsc slabs of one MeshBlock; byte-identical to RaRstCopy
void RadiationVet::RotAvgRstPackHost(HostArray5D<Real> &out, const int nmb,
                                     const int nsc) {
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
  DvceArray5D<Real> sl("rst-vetr-chunk", 1, cs, n3, n2, n1);
  for (int m = 0; m < nmb; ++m) {
    for (int s0 = 0; s0 < nsl; s0 += cs) {
      const int ns = std::min(cs, nsl - s0);
      const int64_t q0 = static_cast<int64_t>(s0)*ncout;
      par_for("vet_ra_rst", DevExeSpace(), 0, static_cast<int>(ns*ncout - 1),
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
      // 4-D subviews, contiguous with equal strides
      Kokkos::deep_copy(Kokkos::subview(out, m, std::make_pair(s0, s0+ns), Kokkos::ALL,
                                        Kokkos::ALL, Kokkos::ALL),
                        Kokkos::subview(sl, 0, std::make_pair(0, ns), Kokkos::ALL,
                                        Kokkos::ALL, Kokkos::ALL));
    }
  }
}

//! load the file's history into vra_d when it fits this run (same K, closure, MeshBlock
//! size); else nothing (the first build re-measures)
bool RadiationVet::RotAvgRstLoad(const int64_t *h, const DvceArray5D<Real> &src,
                                 const int nmb) {
  // not averaged through vra_d (no rotation, or the clean ring): skip; cln_rhist records
  // an EMA history under the clean ring (its closure is converged)
  if (!ScRotOn() || cln_ring) {
    cln_rhist = cln_ring && vra_ema && h[0] == 1;
    return false;
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int64_t ncint = static_cast<int64_t>(indcs.nx1)*indcs.nx2*indcs.nx3;
  const int64_t close = (vgd_on ? 1 : 2) + (vra_ema ? 16 : 0);
  const bool ok = vra_on && h[0] == 1 && h[2] == vra_k && h[4] == close && h[3] > 0 &&
                  h[3] < VET_RA_NF && h[9] == ncint && h[5] >= 0 && h[5] < vra_k - 1 &&
                  (h[6] == 0 || h[6] == 3) && h[7] == ((h[6] == 3) ? vra_k : vra_k - 1) &&
                  h[10] >= 0 && h[10] < indcs.nx1;
  if (!ok) {
    if (global_variable::my_rank == 0) {
      std::cout << "### WARNING: the <rad_vet> vet_rotate_avg history of this restart "
                << "file does not fit this run (K " << h[2] << ", closure " << h[4]
                << "); it is re-measured at the first build and this restart is not "
                << "bitwise."
                << std::endl;
    }
    return false;
  }
  const int nf = static_cast<int>(h[3]);
  vra_gen += 1;   // vra_d may change (dt_retry_host skips clean saves)
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

//! first averaging call after a load: the registered fields must match the file's,
//! else the history is dropped (re-measured)
bool RadiationVet::RotAvgRstCheck() {
  if (vra_rst_nf < 0) {return true;}
  const bool ok = (vra_tab.nf == vra_rst_nf);
  vra_rst_nf = -1;
  if (!ok) {
    if (global_variable::my_rank == 0) {
      std::cout << "### WARNING: <rad_vet> vet_rotate_avg: the restart file's history "
                   "has "
                << "other fields than this run's closure; re-measured at the next build "
                << "(not bitwise)." << std::endl;
    }
    vra_state = 0;
    vra_pend = 0;
  }
  return ok;
}

}  // namespace radvet
