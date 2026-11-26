//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_bandring.cpp
//! \brief curved-mesh rays: wide band slabs lent to the shells that need them.
//!
//!  called from: rays/sc_cs_sweep.cpp, rays/sc_sph_halo.cpp, rays/sc_sph_halo_fused.cpp,
//!    rays/sc_sph_sweep.cpp
//!  calls into: rays/sc_sph_globalband.cpp, rays/sc_sph_halo_exact.cpp

#include <cstdint>

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

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScRingSweepBegin / ScRingSweepEnd / ScRingFlush
//! \brief vet_gd_band_ring (rad_vet.hpp): per-sweep bookkeeping and the kernel of a
//! group of slab operations (save a shell's inward band; release = interior back to
//! the compact slab; acquire = interior in and saved band back, then both shells'
//! tables). Pure copies: the sweep reads the values of the dense-band layout.

// band_ring copy kernel (one thread per value; the last thread sets the tables)
struct ScRingK {
  Real *da, *db, *st;
  int64_t *aoff, *boff;
  int *awi, *bwi;
  const int64_t *oc, *ocb, *so;
  const int *sb, *svdl, *svdw, *rsdl, *rsdw;
  int64_t offs0, offs1, ocs0, ocs1, sos0, sos1;
  int64_t nsv, nrl, ncp, nrs, svso, rsso, svn, rsn;
  int sbs0, sbs1, n, nx2, nx3, w, c2b, c3b, nm, svk, rsk, svnr, rsnr;
  int svs, svw, rli, rls, rlw, aqi, aqs, aqw;
  bool two, hrl, haq;
  // xb: db is the other cross-fade set's array; its saved / restored values follow its
  // own records (lists 2, counts svn2 / rsn2, after the first set's)
  const int *svdl2, *svdw2, *rsdl2, *rsdw2;
  int64_t svn2, rsn2;
  int svnr2, rsnr2;
  bool xb;
  // value of the wide slab of slot s (depth wl) of block m at band index (kb, jb);
  // -1 outside the shell's rectangle
  KOKKOS_INLINE_FUNCTION
  int64_t Wadr(const int s, const int wl, const int m, const int v, const int kb,
               const int jb) const {
    const int bk0 = sb[m*sbs0]*wl, bj0 = sb[m*sbs0 + 2*sbs1]*wl;
    const int c3 = nx3 + bk0 + sb[m*sbs0 + sbs1]*wl;
    const int c2 = nx2 + bj0 + sb[m*sbs0 + 3*sbs1]*wl;
    const int kk = kb - (w - bk0), jj = jb - (w - bj0);
    if (kk < 0 || kk >= c3 || jj < 0 || jj >= c2) {return -1;}
    return so[s*sos0 + m*sos1] + v +
           static_cast<int64_t>(n)*(kk + static_cast<int64_t>(c3)*jj);
  }
  // the t-th saved value of a record: its band index (m, v, kb, jb); false: none
  KOKKOS_INLINE_FUNCTION
  bool Rec(const int k, const int wl, const int *dl, const int *dw, const int nr,
           const int64_t t, int &m, int &v, int &kb, int &jb) const {
    if (k == 1) {
      int a = (t < nr) ? dl[t] : dw[t - nr];
      jb = a % c2b;
      a /= c2b;
      kb = a % c3b;
      a /= c3b;
      v = a % n;
      m = a/n;
      return true;
    }
    // dense: (m, v, kb, jb) over the band index rectangle of depth wl, interior out
    const int64_t c3w = nx3 + 2*wl, c2w = nx2 + 2*wl;
    int64_t r = t;
    v = static_cast<int>(r % n);
    r /= n;
    const int kw = static_cast<int>(r % c3w);
    r /= c3w;
    const int jw = static_cast<int>(r % c2w);
    m = static_cast<int>(r/c2w);
    kb = kw + (w - wl);
    jb = jw + (w - wl);
    return !(kb >= w && kb < w + nx3 && jb >= w && jb < w + nx2);
  }
  KOKKOS_INLINE_FUNCTION
  void operator()(const int64_t tt) const {
    int64_t t = tt;
    if (xb && t < nsv) {
      const bool b2 = (t >= svn);
      const int64_t tl = b2 ? (t - svn) : t;
      int m, v, kb, jb;
      if (!Rec(svk, svw, b2 ? svdl2 : svdl, b2 ? svdw2 : svdw, b2 ? svnr2 : svnr, tl,
               m, v, kb, jb)) {return;}
      const int64_t q = Wadr(svs, svw, m, v, kb, jb);
      if (q < 0) {return;}
      st[svso + t] = b2 ? db[q] : da[q];
      return;
    }
    if (t < nsv) {
      int m, v, kb, jb;
      if (!Rec(svk, svw, svdl, svdw, svnr, t, m, v, kb, jb)) {return;}
      const int64_t q = Wadr(svs, svw, m, v, kb, jb);
      if (q < 0) {return;}
      st[svso + t] = da[q];
      if (two) {st[svso + svn + t] = db[q];}
      return;
    }
    t -= nsv;
    if (t < nrl + ncp) {
      const bool out = (t < nrl);
      if (!out) {t -= nrl;}
      const int i = out ? rli : aqi;
      const int s = out ? rls : aqs;
      const int wl = out ? rlw : aqw;
      // 32-bit index arithmetic (nrl, ncp < 2^31, checked on the host)
      int r = static_cast<int>(t);
      const int v = r % n;
      r /= n;
      const int kk = r % nx3;
      r /= nx3;
      const int jj = r % nx2;
      const int m = r/nx2;
      const int64_t qi = v + static_cast<int64_t>(n)*(kk + static_cast<int64_t>(nx3)*jj);
      const int64_t qc = oc[m*ocs0 + i*ocs1] + qi, qb = ocb[m*ocs0 + i*ocs1] + qi;
      const int64_t qw = Wadr(s, wl, m, v, kk + w, jj + w);
      if (out) {
        da[qc] = da[qw];
        if (two) {db[qb] = db[qw];}
      } else {
        da[qw] = da[qc];
        if (two) {db[qw] = db[qb];}
      }
      return;
    }
    t -= nrl + ncp;
    if (xb && t < nrs) {
      const bool b2 = (t >= rsn);
      const int64_t tl = b2 ? (t - rsn) : t;
      int m, v, kb, jb;
      if (!Rec(rsk, aqw, b2 ? rsdl2 : rsdl, b2 ? rsdw2 : rsdw, b2 ? rsnr2 : rsnr, tl,
               m, v, kb, jb)) {return;}
      const int64_t q = Wadr(aqs, aqw, m, v, kb, jb);
      if (q < 0) {return;}
      if (b2) {
        db[q] = st[rsso + t];
      } else {
        da[q] = st[rsso + t];
      }
      return;
    }
    if (t < nrs) {
      int m, v, kb, jb;
      if (!Rec(rsk, aqw, rsdl, rsdw, rsnr, t, m, v, kb, jb)) {return;}
      const int64_t q = Wadr(aqs, aqw, m, v, kb, jb);
      if (q < 0) {return;}
      da[q] = st[rsso + t];
      if (two) {db[q] = st[rsso + rsn + t];}
      return;
    }
    // the tables (no other thread of this kernel reads them)
    for (int m = 0; m < nm; ++m) {
      if (hrl) {
        aoff[m*offs0 + rli*offs1] = oc[m*ocs0 + rli*ocs1];
        if (two) {boff[m*offs0 + rli*offs1] = ocb[m*ocs0 + rli*ocs1];}
      }
      if (haq) {
        aoff[m*offs0 + aqi*offs1] = so[aqs*sos0 + m*sos1];
        if (two) {boff[m*offs0 + aqi*offs1] = so[aqs*sos0 + m*sos1];}
      }
    }
    if (hrl) {
      awi[rli] = 0;
      if (two) {bwi[rli] = 0;}
    }
    if (haq) {
      awi[aqi] = aqw;
      if (two) {bwi[aqi] = aqw;}
    }
  }
};

static_assert(sizeof(ScRingK) <= 512, "ScRingK: keep below the by-value launch limit");

void RadiationVet::ScRingSweepBegin(const bool two) {
  vgd_ring_act = vgd_ring;
  if (!vgd_ring_act) {return;}
  for (size_t i = 0; i < vgd_rg_state.size(); ++i) {
    if (vgd_rg_state[i] >= 0) {ScFatal("vet_gd_band_ring: a shell still holds a slot");}
  }
  for (int s = 0; s < vgd_rg_k; ++s) {vgd_rg_sown[s] = -1;}
  vgd_rg_pend.clear();
  vgd_rg_two = two;
  vgd_rg_xb = vgd_xff_now;   // the second array is the other set's
  vgd_rg_stop = 0;
  vgd_rg_nk += 1.0;
  // store of the saved inward band values: the exact lists of this sweep's inward
  // exchanges (made by the batch; grown below if needed)
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int n1c = indcs.nx1 + 2*indcs.ng;
  int64_t need = 0;
  if (static_cast<int>(vgd_hxe.size()) == 2*n1c) {
    for (int i = indcs.is; i <= indcs.ie; ++i) {
      const ScHlEntry &X = vgd_hxe[i];
      if (X.gen != -1 && X.exact && X.alpha == vgd_alpha && ScHxCut(i, X)) {
        need += static_cast<int64_t>(X.nr) + X.nw;
      }
    }
  }
  need *= (two ? 2 : 1);
  if (static_cast<int64_t>(vgd_rg_store.extent(0)) < need) {
    Kokkos::fence();
    vgd_rg_store = decltype(vgd_rg_store)();
    Kokkos::realloc(vgd_rg_store, need + need/16 + 1024);
  }
}

void RadiationVet::ScRingSweepEnd() {
  if (!vgd_ring_act) {return;}
  ScRingFlush();
  for (size_t i = 0; i < vgd_rg_state.size(); ++i) {
    if (vgd_rg_state[i] >= 0) {ScFatal("vet_gd_band_ring: a shell kept its slot");}
  }
  vgd_ring_act = false;
}

int RadiationVet::ScRingFreeSlot(const int avoid) const {
  for (int s = 0; s < vgd_rg_k; ++s) {
    if (vgd_rg_sown[s] < 0 && s != avoid) {return s;}
  }
  return -1;
}

void RadiationVet::ScRingFlush() {
  if (!vgd_ring_act || vgd_rg_pend.empty()) {return;}
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int n = vgd_n, nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int c2b = nx2 + 2*w, c3b = nx3 + 2*w;
  const bool two = vgd_rg_two;
  const bool xb = two && vgd_rg_xb;
  // one group: at most one save, one release, one acquire
  size_t p = 0;
  while (p < vgd_rg_pend.size()) {
    ScRingOp sv, rl, aq;
    while (p < vgd_rg_pend.size()) {
      const ScRingOp &o = vgd_rg_pend[p];
      if ((o.kind == 1 && sv.kind) || (o.kind == 2 && rl.kind) ||
          (o.kind == 3 && aq.kind)) {
        break;
      }
      if (o.kind == 1) {
        sv = o;
        sv.slot = vgd_rg_state[o.i];
        if (sv.slot < 0) {ScFatal("vet_gd_band_ring: save of a compact shell");}
      } else if (o.kind == 2) {
        rl = o;
        rl.slot = vgd_rg_state[o.i];
        if (rl.slot < 0) {ScFatal("vet_gd_band_ring: release of a compact shell");}
        if (sv.kind && sv.i != o.i) {ScFatal("vet_gd_band_ring: save / release differ");}
      } else {
        aq = o;
        if (vgd_rg_state[o.i] >= 0) {ScFatal("vet_gd_band_ring: shell already lent");}
        aq.slot = ScRingFreeSlot(rl.kind ? rl.slot : -1);
        if (aq.slot < 0) {ScFatal("vet_gd_band_ring: no free slot");}
      }
      ++p;
    }
    // save record of the shell saved here, and the store
    int64_t sv_so = 0, sv_n = 0, sv_n2 = 0;
    int sv_k = 0;
    DvceArray1D<int> sv_dl, sv_dw, sv_dl2, sv_dw2;
    int sv_nr = 0, sv_nr2 = 0;
    if (sv.kind) {
      const int i = sv.i;
      const int wl = vgd_rg_wfull[i];
      if (vgd_rg_sk[i] != 1 && vgd_rg_sk[i] != 3) {vgd_rg_sk[i] = 2;}   // whole band
      sv_k = vgd_rg_sk[i];
      if (xb && sv_k == 3) {   // the other set's own kept list
        ScGbSwap();
        sv_dl2 = ScRawRec(vgd_rg_se[i]);
        ScGbSwap();
        sv_nr2 = sv_dl2.extent_int(0);
      }
      if (sv_k == 3) {   // vgd_raw: the positions the raw exchange wrote
        sv_dl = ScRawRec(vgd_rg_se[i]);
        sv_nr = sv_dl.extent_int(0);
        sv_n = sv_nr;
        sv_k = 1;
      } else if (sv_k == 1) {
        const ScHlEntry &X = vgd_hxe[vgd_rg_se[i]];
        sv_dl = X.dl;
        sv_dw = X.dw;
        sv_nr = X.nr;
        sv_n = static_cast<int64_t>(X.nr) + X.nw;
      } else {
        sv_n = static_cast<int64_t>(nmb)*n*(nx3 + 2*wl)*(nx2 + 2*wl);
        vgd_rg_ndense += 1.0;
      }
      sv_n2 = (xb && vgd_rg_sk[i] == 3) ? sv_nr2 : sv_n;
      if (xb && sv_k == 1 && vgd_rg_sk[i] != 3) {
        const ScHlEntry &X2 = vgd_hxe_alt[vgd_rg_se[i]];
        sv_dl2 = X2.dl;
        sv_dw2 = X2.dw;
        sv_nr2 = X2.nr;
        sv_n2 = static_cast<int64_t>(X2.nr) + X2.nw;
      }
      if (static_cast<int>(vgd_rg_sn2.size()) != static_cast<int>(vgd_rg_sn.size())) {
        vgd_rg_sn2.assign(vgd_rg_sn.size(), 0);
      }
      vgd_rg_sn2[i] = sv_n2;
      const int64_t need = vgd_rg_stop + (xb ? (sv_n + sv_n2) : sv_n*(two ? 2 : 1));
      if (static_cast<int64_t>(vgd_rg_store.extent(0)) < need) {
        Kokkos::fence();
        Kokkos::resize(vgd_rg_store, need + need/8 + 1024);
      }
      sv_so = vgd_rg_stop;
      vgd_rg_so[i] = sv_so;
      vgd_rg_sn[i] = sv_n;
      vgd_rg_stop = need;
      // keep the saved values' lists with the record (the restore reads them)
    }
    // restore record of the acquired shell (outward pass)
    int64_t rs_so = 0, rs_n = 0, rs_n2 = 0;
    int rs_k = 0, rs_nr = 0, rs_nr2 = 0;
    DvceArray1D<int> rs_dl, rs_dw, rs_dl2, rs_dw2;
    if (aq.kind && aq.restore) {
      const int i = aq.i;
      rs_k = vgd_rg_sk[i];
      rs_so = vgd_rg_so[i];
      rs_n = vgd_rg_sn[i];
      if (xb && rs_k == 3) {   // the other set's raw record
        ScGbSwap();
        rs_dl2 = ScRawRec(vgd_rg_se[i]);
        ScGbSwap();
        rs_nr2 = rs_dl2.extent_int(0);
      }
      if (rs_k == 3) {
        rs_dl = ScRawRec(vgd_rg_se[i]);
        rs_nr = rs_dl.extent_int(0);
        rs_k = 1;
        if (rs_nr != rs_n) {ScFatal("vet_gd_band_ring: the raw record changed");}
      } else if (rs_k == 1) {
        const ScHlEntry &X = vgd_hxe[vgd_rg_se[i]];
        rs_dl = X.dl;
        rs_dw = X.dw;
        rs_nr = X.nr;
        if (static_cast<int64_t>(X.nr) + X.nw != rs_n) {
          ScFatal("vet_gd_band_ring: the exact list changed between save and restore");
        }
      } else if (rs_k != 2) {
        ScFatal("vet_gd_band_ring: restore without a save");
      }
      if (xb) {
        rs_n2 = vgd_rg_sn2[i];
        if (vgd_rg_sk[i] == 3 && rs_nr2 != rs_n2) {
          ScFatal("vet_gd_band_ring: the other set's raw record changed");
        }
        if (rs_k == 1 && vgd_rg_sk[i] != 3) {
          const ScHlEntry &X2 = vgd_hxe_alt[vgd_rg_se[i]];
          rs_dl2 = X2.dl;
          rs_dw2 = X2.dw;
          rs_nr2 = X2.nr;
          if (static_cast<int64_t>(X2.nr) + X2.nw != rs_n2) {
            ScFatal("vet_gd_band_ring: the other set's exact list changed between save "
                    "and restore");
          }
        }
      }
    }
    auto a_ = vgd_i;
    auto b_ = xb ? vgd_i_xb : vgd_i;
    auto sb_ = vgd_i.sb;
    auto oc_ = a_.offc, ocb_ = b_.offc;
    auto so_ = vgd_rg_dsoff;
    auto st_ = vgd_rg_store;
    const int64_t nsv = sv.kind ? (xb ? (sv_n + sv_n2) : sv_n) : 0;
    const int64_t nin = static_cast<int64_t>(nmb)*n*nx3*nx2;
    const int64_t nrl = rl.kind ? nin : 0;
    const int64_t ncp = (aq.kind && aq.copyin) ? nin : 0;
    const int64_t nrs = (aq.kind && aq.restore) ? (xb ? (rs_n + rs_n2) : rs_n) : 0;
    const int64_t ntot = nsv + nrl + ncp + nrs + 1;
    if (nin >= static_cast<int64_t>(INT32_MAX)) {
      ScFatal("vet_gd_band_ring: shell interior >= 2^31 values (32-bit copy index)");
    }
    // the copy kernel: a small POD functor of raw pointers (a functor of Views exceeds
    // the 512 B by-value launch limit and would go through constant memory)
    ScRingK k;
    k.da = a_.d.data();
    k.db = b_.d.data();
    k.st = st_.data();
    k.aoff = a_.off.data();
    k.boff = b_.off.data();
    k.awi = a_.wi.data();
    k.bwi = b_.wi.data();
    k.oc = oc_.data();
    k.ocb = ocb_.data();
    k.so = so_.data();
    k.sb = sb_.data();
    k.offs0 = a_.off.stride(0);
    k.offs1 = a_.off.stride(1);
    k.ocs0 = oc_.stride(0);
    k.ocs1 = oc_.stride(1);
    k.sos0 = so_.stride(0);
    k.sos1 = so_.stride(1);
    k.sbs0 = static_cast<int>(sb_.stride(0));
    k.sbs1 = static_cast<int>(sb_.stride(1));
    if (b_.off.stride(0) != a_.off.stride(0) || b_.off.stride(1) != a_.off.stride(1)) {
      ScFatal("vet_gd_band_ring: main / twin offset tables differ in layout");
    }
    k.svdl = sv_dl.data();
    k.svdw = sv_dw.data();
    k.rsdl = rs_dl.data();
    k.rsdw = rs_dw.data();
    k.nsv = nsv;
    k.nrl = nrl;
    k.ncp = ncp;
    k.nrs = nrs;
    k.svso = sv_so;
    k.rsso = rs_so;
    k.svn = sv_n;
    k.rsn = rs_n;
    k.n = n;
    k.nx2 = nx2;
    k.nx3 = nx3;
    k.w = w;
    k.c2b = c2b;
    k.c3b = c3b;
    k.nm = nmb;
    k.svk = sv_k;
    k.rsk = rs_k;
    k.svnr = sv_nr;
    k.rsnr = rs_nr;
    k.svs = sv.slot;
    k.svw = sv.kind ? vgd_rg_wfull[sv.i] : 0;
    k.rli = rl.i;
    k.rls = rl.slot;
    k.rlw = rl.kind ? vgd_rg_wfull[rl.i] : 0;
    k.aqi = aq.i;
    k.aqs = aq.slot;
    k.aqw = aq.kind ? vgd_rg_wfull[aq.i] : 0;
    k.two = two;
    k.hrl = rl.kind;
    k.haq = aq.kind;
    k.xb = xb;
    k.svdl2 = sv_dl2.data();
    k.svdw2 = sv_dw2.data();
    k.rsdl2 = rs_dl2.data();
    k.rsdw2 = rs_dw2.data();
    k.svn2 = sv_n2;
    k.rsn2 = rs_n2;
    k.svnr2 = sv_nr2;
    k.rsnr2 = rs_nr2;
    Kokkos::parallel_for("vet_vgd_ring", Kokkos::RangePolicy<>(vgd_cur, 0, ntot), k);
    if (rl.kind) {
      vgd_rg_sown[rl.slot] = -1;
      vgd_rg_state[rl.i] = -1;
    }
    if (aq.kind) {
      vgd_rg_sown[aq.slot] = aq.i;
      vgd_rg_state[aq.i] = aq.slot;
    }
  }
  vgd_rg_pend.clear();
}

}  // namespace radvet
