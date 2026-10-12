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
//! Clamps (kind 1): the half-range ratios to their ranges, C1 >= 0; the tensor's D_rr is
//! clamped by the fold (VetLatBuild).
//!
//! VALIDITY.  A cell is valid in a set when i >= the closure's first shell (ilo) and the
//! flag field (if any) is >= 0; an offset is kept only where both sweeps were valid, and
//! the mean runs over the valid offsets (+ the current set).
//!
//! COST.  One extra sweep per rotation block (+ K-1 sweeps at the first build); memory
//! (K-1)(nf+1) reals per active cell.

#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "rad_m1/rad_m1.hpp"

namespace radm1 {

namespace {
KOKKOS_INLINE_FUNCTION
bool RaValid(const DvceArray5D<Real> &a0, const DvceArray5D<Real> &a1, const int farr,
             const int fcomp, const int ilo, const int m, const int k, const int j,
             const int i) {
  if (i < ilo) return false;
  if (farr < 0) return true;
  const Real f = (farr == 0) ? a0(m,fcomp,k,j,i) : a1(m,fcomp,k,j,i);
  return f >= 0.0;
}
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::RotAvgSnap
//! \brief the registered fields as they are (raw, one set) -> ring slot `slot`, with the
//! validity of each cell

void RadiationM1::RotAvgSnap(const int slot, const int ilo) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nf = vra_tab.nf;
  const int nslot = vra_k - 1;
  if (vra_d.extent_int(0) < nmb1 + 1 || vra_d.extent_int(1) != nslot*(nf + 1)) {
    Kokkos::realloc(vra_d, nmb1 + 1, nslot*(nf + 1), ke - ks + 1, je - js + 1,
                    ie - is + 1);
    Kokkos::deep_copy(vra_d, 0.0);
  }
  auto d_ = vra_d;
  auto a0_ = vra_a0;
  auto a1_ = vra_a1;
  const RotAvgTab t = vra_tab;
  const int o = slot*(nf + 1);
  par_for("m1_ra_snap", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const bool v = RaValid(a0_, a1_, t.farr, t.fcomp, ilo, m, k, j, i);
    for (int f = 0; f < t.nf; ++f) {
      d_(m,o+f,k-ks,j-js,i-is) = (t.arr[f] == 0) ? a0_(m,t.comp[f],k,j,i)
                                                 : a1_(m,t.comp[f],k,j,i);
    }
    d_(m,o+t.nf,k-ks,j-js,i-is) = v ? 1.0 : 0.0;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::RotAvgUpdate
//! \brief with the CURRENT set's raw fields in place: turn the pending raw snapshots into
//! offsets (init: every slot S_j -> S_j - g; rotation: the head slot S_old -> S_old - g =
//! -Delta, every other slot += -Delta; ratio fields multiplicatively)

void RadiationM1::RotAvgUpdate(const int ilo) {
  if (vra_pend == 0) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nslot = vra_k - 1;
  auto d_ = vra_d;
  auto a0_ = vra_a0;
  auto a1_ = vra_a1;
  const RotAvgTab t = vra_tab;
  const int nf1 = t.nf + 1;
  const int h = vra_head;
  const bool init = (vra_pend == 1);
  par_for("m1_ra_upd", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int kk = k - ks, jj = j - js, ii = i - is;
    const bool vg = RaValid(a0_, a1_, t.farr, t.fcomp, ilo, m, k, j, i);
    if (init) {
      for (int s = 0; s < nslot; ++s) {
        const int o = s*nf1;
        const bool v = vg && (d_(m,o+t.nf,kk,jj,ii) > 0.5);
        for (int f = 0; f < t.nf; ++f) {
          const Real g = (t.arr[f] == 0) ? a0_(m,t.comp[f],k,j,i) : a1_(m,t.comp[f],k,j,i);
          const Real sv = d_(m,o+f,kk,jj,ii);
          Real x;
          if (t.kind[f] == 2) {
            x = (v && g > 0.0 && sv > 0.0) ? sv/g : 1.0;
          } else {
            x = v ? (sv - g) : 0.0;
          }
          d_(m,o+f,kk,jj,ii) = x;
        }
        d_(m,o+t.nf,kk,jj,ii) = v ? 1.0 : 0.0;
      }
      return;
    }
    // rotation: the head slot holds the old set's raw fields of this state
    const int oh = h*nf1;
    const bool vh = vg && (d_(m,oh+t.nf,kk,jj,ii) > 0.5);
    for (int f = 0; f < t.nf; ++f) {
      const Real g = (t.arr[f] == 0) ? a0_(m,t.comp[f],k,j,i) : a1_(m,t.comp[f],k,j,i);
      const Real sv = d_(m,oh+f,kk,jj,ii);
      const bool mul = (t.kind[f] == 2);
      Real x;
      if (mul) {
        x = (vh && g > 0.0 && sv > 0.0) ? sv/g : 1.0;
      } else {
        x = vh ? (sv - g) : 0.0;
      }
      d_(m,oh+f,kk,jj,ii) = x;
      for (int s = 0; s < nslot; ++s) {
        if (s == h) continue;
        const int o = s*nf1;
        if (!(vh && d_(m,o+t.nf,kk,jj,ii) > 0.5)) {
          d_(m,o+f,kk,jj,ii) = mul ? 1.0 : 0.0;
        } else if (mul) {
          d_(m,o+f,kk,jj,ii) *= x;
        } else {
          d_(m,o+f,kk,jj,ii) += x;
        }
      }
    }
    for (int s = 0; s < nslot; ++s) {
      if (s == h) continue;
      const int o = s*nf1;
      d_(m,o+t.nf,kk,jj,ii) = (vh && d_(m,o+t.nf,kk,jj,ii) > 0.5) ? 1.0 : 0.0;
    }
    d_(m,oh+t.nf,kk,jj,ii) = vh ? 1.0 : 0.0;
  });
  vra_pend = 0;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::RotAvgApply
//! \brief the registered fields := the mean over the current set and the valid offsets

void RadiationM1::RotAvgApply(const int ilo) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nslot = vra_k - 1;
  auto d_ = vra_d;
  auto a0_ = vra_a0;
  auto a1_ = vra_a1;
  const RotAvgTab t = vra_tab;
  const int nf1 = t.nf + 1;
  par_for("m1_ra_apply", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (!RaValid(a0_, a1_, t.farr, t.fcomp, ilo, m, k, j, i)) return;
    const int kk = k - ks, jj = j - js, ii = i - is;
    Real nv = 1.0;
    for (int s = 0; s < nslot; ++s) {nv += d_(m,s*nf1+t.nf,kk,jj,ii);}
    if (nv < 1.5) return;
    for (int f = 0; f < t.nf; ++f) {
      Real sm = (t.kind[f] == 2) ? 1.0 : 0.0;
      for (int s = 0; s < nslot; ++s) {
        const int o = s*nf1;
        if (d_(m,o+t.nf,kk,jj,ii) > 0.5) sm += d_(m,o+f,kk,jj,ii);
      }
      Real &g = (t.arr[f] == 0) ? a0_(m,t.comp[f],k,j,i) : a1_(m,t.comp[f],k,j,i);
      if (t.kind[f] == 2) {
        g *= sm/nv;
      } else {
        Real v = g + sm/nv;
        if (t.kind[f] == 1) v = fmin(fmax(v, t.lo[f]), t.hi[f]);
        g = v;
      }
    }
  });
}

}  // namespace radm1
