#ifndef UTILS_TWO_STREAM_CK_RST_STATE_HPP_
#define UTILS_TWO_STREAM_CK_RST_STATE_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file two_stream_ck_rst_state.hpp
//! \brief the heavy half of utils/two_stream_ck_rst.hpp: collect the implicit ck
//! solver's cross-call state for a restart file, put it back on the first call(s) of a
//! restarted run, and REBUILD the ck_impl_xstep stored operator from the state its
//! storing call saw.  Included by two_stream_rt.hpp right behind two_stream_column_ck.hpp
//! (it needs the solver's pointers, not the sweep).  Hooks in two_stream_rt.hpp:
//!   CkRstAfterAlloc  after CkImplAlloc (lazy allocation of the solver arrays)
//!   CkXsSnap         in the ck_impl_xstep store branch of pass 0
//!   CkRstBeginCall   at the top of every implicit call (picket_fence_two_stream_RT)
//!   CkRstCadRestore  at the top of CkCadStep
//! Every hook is a no-op without a restart block, so a run that never reads one is
//! unchanged; the snapshot copies of CkXsSnap only exist with ck_impl_xstep > 0 and do
//! not touch the solution.

#include <cstdint>
#include <cstring>
#include <iostream>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "utils/two_stream_ck_rst.hpp"
#include "utils/two_stream_column_ck.hpp"
#include "utils/deep_copy_across.hpp"

namespace two_stream_rt {

using deep_copy_across::DeepCopyAcross;

inline void picket_fence_two_stream_RT_pass(Mesh *pm, Real bdt);
inline void CkRstScalars();

// ---- the ck_impl_xstep snapshot: the state pass 0 of the last STORING call read, and
// that call's bdt.  Refreshed on every store (a device-to-device copy of u0, wtemp,
// bcc0 and ck_thk, i.e. ~7 Reals per cell every xstep cycles).
inline DvceArray5D<Real> *ck_xsU_ptr = nullptr;
inline DvceArray4D<Real> *ck_xsW_ptr = nullptr;
inline DvceArray4D<Real> *ck_xsK_ptr = nullptr;
inline DvceArray5D<Real> *ck_xsB_ptr = nullptr;
inline Real ck_xs_bdt = 0.0;
inline bool ck_xs_have = false;
// restart bookkeeping
inline bool ck_rst_scalars_done = false;
inline bool ck_rst_rebuild_active = false;
inline bool ck_rst_warned = false;
inline bool ck_rst_first_done = false;      // the first implicit call has begun

//----------------------------------------------------------------------------------------
//! \fn void CkXsSnap
//! \brief called by pass 0 of a call that STORES the xstep operator, before rt_pre_opac
//! and the apply read or rewrite ck_thk: everything that pass reads from the gas, and
//! the thin/thick verdict it starts from.

inline void CkXsSnap(Mesh *pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  const bool hyd = (pmbp->phydro != nullptr);
  DvceArray5D<Real> u0 = hyd ? pmbp->phydro->u0 : pmbp->pmhd->u0;
  DvceArray4D<Real> wt = hyd ? pmbp->phydro->wtemp : pmbp->pmhd->wtemp;
  if (ck_xsU_ptr == nullptr || ck_xsU_ptr->extent(0) != u0.extent(0)) {
    if (ck_xsU_ptr != nullptr) {
      delete ck_xsU_ptr;
      delete ck_xsW_ptr;
      delete ck_xsK_ptr;
      if (ck_xsB_ptr != nullptr) delete ck_xsB_ptr;
      ck_xsB_ptr = nullptr;
    }
    ck_xsU_ptr = new DvceArray5D<Real>("ck_xsU", u0.extent(0), u0.extent(1),
                                       u0.extent(2), u0.extent(3), u0.extent(4));
    ck_xsW_ptr = new DvceArray4D<Real>("ck_xsW", wt.extent(0), wt.extent(1),
                                       wt.extent(2), wt.extent(3));
    auto &thk = *ck_thk_ptr;
    ck_xsK_ptr = new DvceArray4D<Real>("ck_xsK", thk.extent(0), thk.extent(1),
                                       thk.extent(2), thk.extent(3));
    if (!hyd) {
      auto &bc = pmbp->pmhd->bcc0;
      ck_xsB_ptr = new DvceArray5D<Real>("ck_xsB", bc.extent(0), bc.extent(1),
                                         bc.extent(2), bc.extent(3), bc.extent(4));
    }
  }
  Kokkos::deep_copy(DevExeSpace(), *ck_xsU_ptr, u0);
  if (wt.size() > 0) Kokkos::deep_copy(DevExeSpace(), *ck_xsW_ptr, wt);
  Kokkos::deep_copy(DevExeSpace(), *ck_xsK_ptr, *ck_thk_ptr);
  if (!hyd) Kokkos::deep_copy(DevExeSpace(), *ck_xsB_ptr, pmbp->pmhd->bcc0);
  ck_xs_bdt = bdt;
  ck_xs_have = true;
}

//----------------------------------------------------------------------------------------
//! \fn int CkRstCollect
//! \brief restart.cpp's LoadOutputData: the slabs (nslab, nmb, n3, n2, n1) and header.
//! Must return the same slab list on every rank (rank 0 writes the header).

inline int CkRstCollect(Mesh *pm, HostArray5D<Real> &out, CkRstHdr &h, int nmb,
                        int n3, int n2, int n1) {
  std::memset(&h, 0, sizeof(h));
  // a restarted run that has not made an implicit call yet (e.g. only linearised
  // cadence steps since the restart) still holds its state in the staging: pass it on
  if (ck_rst_have) CkRstScalars();
  const bool thk_stg = (ck_thk_ptr == nullptr) && CkRstStaged(kCkSlabThk);
  if (!ck_implicit || (ck_thk_ptr == nullptr && !thk_stg)) return 0;
  MeshBlockPack *pmbp = pm->pmb_pack;
  const bool hyd = (pmbp->phydro != nullptr);
  DvceArray5D<Real> u0 = hyd ? pmbp->phydro->u0 : pmbp->pmhd->u0;
  // does the next call re-apply a stored operator on some rank?  Then every rank writes
  // the snapshot slabs (a rank without one writes zeros, and its XSCYC says -1)
  const bool xs_stg = !ck_xs_have && CkRstStaged(kCkSlabXsK);
  int need = (ck_impl_xstep > 0 && (ck_xs_have || xs_stg) && ck_impl_xs_cyc >= 0 &&
              pm->ncycle - ck_impl_xs_cyc < ck_impl_xstep) ? 1 : 0;
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &need, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
  std::vector<int> ids;
  ids.push_back(kCkSlabThk);
  const bool cad_stg = (ck_cad_ptr == nullptr) && CkRstStaged(kCkSlabCad0);
  if (ck_cad_ptr != nullptr || cad_stg) {
    for (int s=0; s<4; ++s) ids.push_back(kCkSlabCad0 + s);
  }
  if (ck_impl_cvkeep && (ck_cv_ptr != nullptr || CkRstStaged(kCkSlabCv))) {
    ids.push_back(kCkSlabCv);
  }
  if (ck_impl_xstep > 0) {
    ids.push_back(kCkSlabXsCyc);
    ids.push_back(kCkSlabXsBdt);
  }
  if (need) {
    ids.push_back(kCkSlabXsW);
    ids.push_back(kCkSlabXsK);
    const int nv = static_cast<int>(u0.extent(1));
    for (int n=0; n<nv; ++n) ids.push_back(kCkSlabXsU0 + n);
    if (!hyd) {
      for (int n=0; n<3; ++n) ids.push_back(kCkSlabXsB0 + n);
    }
  }
  const int ns = static_cast<int>(ids.size());
  if (ns > kCkRstMaxSlab) return 0;
  Kokkos::realloc(out, ns, nmb, n3, n2, n1);
  auto nm = std::make_pair(0, nmb);
  using Kokkos::ALL;
  for (int s=0; s<ns; ++s) {
    const int id = ids[s];
    auto dst = Kokkos::subview(out, s, ALL, ALL, ALL, ALL);
    // still staged (not consumed since the restart): copy the staged bytes
    const bool stg = CkRstStaged(id) &&
        ((id == kCkSlabThk && thk_stg) || (id == kCkSlabCv && ck_cv_ptr == nullptr) ||
         (id >= kCkSlabCad0 && id < kCkSlabCad0 + 4 && cad_stg) ||
         (id >= kCkSlabXsW && xs_stg));
    if (stg) {
      const std::size_t nd = dst.size();
      if (ck_rst_stage[id].size() == nd) {
        std::memcpy(dst.data(), ck_rst_stage[id].data(), nd*sizeof(Real));
      } else {
        Kokkos::deep_copy(dst, 0.0);
      }
    } else if (id == kCkSlabThk) {
      DeepCopyAcross(dst, Kokkos::subview(*ck_thk_ptr, nm, ALL, ALL, ALL));
    } else if (id >= kCkSlabCad0 && id < kCkSlabCad0 + 4) {
      DeepCopyAcross(dst, Kokkos::subview(*ck_cad_ptr, nm, id - kCkSlabCad0,
                                          ALL, ALL, ALL));
    } else if (id == kCkSlabCv) {
      DeepCopyAcross(dst, Kokkos::subview(*ck_cv_ptr, nm, ALL, ALL, ALL));
    } else if (id == kCkSlabXsCyc) {
      Kokkos::deep_copy(dst, static_cast<Real>(ck_impl_xs_cyc));
    } else if (id == kCkSlabXsBdt) {
      Kokkos::deep_copy(dst, ck_xs_bdt);
    } else if (!ck_xs_have) {
      Kokkos::deep_copy(dst, 0.0);
    } else if (id == kCkSlabXsW) {
      if (ck_xsW_ptr->size() > 0) {
        DeepCopyAcross(dst, Kokkos::subview(*ck_xsW_ptr, nm, ALL, ALL, ALL));
      } else {
        Kokkos::deep_copy(dst, 0.0);
      }
    } else if (id == kCkSlabXsK) {
      DeepCopyAcross(dst, Kokkos::subview(*ck_xsK_ptr, nm, ALL, ALL, ALL));
    } else if (id >= kCkSlabXsU0 && id < kCkSlabXsU0 + 16) {
      DeepCopyAcross(dst, Kokkos::subview(*ck_xsU_ptr, nm, id - kCkSlabXsU0,
                                          ALL, ALL, ALL));
    } else if (id >= kCkSlabXsB0 && id < kCkSlabXsB0 + 3) {
      DeepCopyAcross(dst, Kokkos::subview(*ck_xsB_ptr, nm, id - kCkSlabXsB0,
                                          ALL, ALL, ALL));
    }
  }
  h.version = 1;
  h.nslab = ns;
  for (int s=0; s<ns; ++s) h.id[s] = ids[s];
  h.xs_cyc = ck_impl_xs_cyc;
  h.jac_built = ck_impl_jac_built ? 1 : 0;
  h.nstore = ck_impl_nstore;
  h.nreuse = ck_impl_nreuse;
  h.nsweep = ck_impl_nsweep;
  h.nfull = ck_cad_nfull;
  h.nlin = ck_cad_nlin;
  h.nguard = ck_cad_nguard;
  h.fref = ck_cad_fref;
  h.xs_bdt = static_cast<double>(ck_xs_bdt);
  return ns;
}
inline const bool ck_rst_registered_ = (ck_rst_collect_fn = &CkRstCollect, true);

//----------------------------------------------------------------------------------------
//! \fn bool CkRstPut
//! \brief copy staged slab id into dst (a device View or subview of the slab's shape)
//! and free it.  False (and the slab dropped) when it does not fit.

template <typename V>
inline bool CkRstPut(const int id, const V &dst) {
  if (!CkRstStaged(id)) return false;
  using UnmanagedHost4D = Kokkos::View<Real****, Kokkos::LayoutRight,
                                       Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
  const bool fit = (dst.extent(0) >= static_cast<std::size_t>(ck_rst_nmb)) &&
                   dst.extent(1) == static_cast<std::size_t>(ck_rst_n3) &&
                   dst.extent(2) == static_cast<std::size_t>(ck_rst_n2) &&
                   dst.extent(3) == static_cast<std::size_t>(ck_rst_n1);
  if (!fit) {
    if (global_variable::my_rank == 0) {
      std::cout << "### WARNING: ck restart slab " << id << " does not fit this "
                << "decomposition; dropped (this restart is not bitwise)." << std::endl;
    }
    CkRstStageFree(id);
    return false;
  }
  UnmanagedHost4D h(ck_rst_stage[id].data(), ck_rst_nmb, ck_rst_n3, ck_rst_n2,
                    ck_rst_n1);
  DeepCopyAcross(Kokkos::subview(dst, std::make_pair(0, ck_rst_nmb), Kokkos::ALL,
                                 Kokkos::ALL, Kokkos::ALL), h);
  CkRstStageFree(id);
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void CkRstScalars
//! \brief once, on the first call of a restarted run that carried a block.

inline void CkRstScalars() {
  if (ck_rst_scalars_done) return;
  ck_rst_scalars_done = true;
  if (!ck_rst_have) {
    if (ck_rst_restarted && ck_implicit && !ck_rst_warned &&
        global_variable::my_rank == 0) {
      ck_rst_warned = true;
      std::cout << "### WARNING: restart file has no implicit-ck state block (written "
                << "before it was added, or with ck_implicit off): ck_thk starts all "
                << "thin, no xstep operator is stored and the cadence starts with a full "
                << "call; this restart is not bitwise." << std::endl;
    }
    return;
  }
  const CkRstHdr &h = ck_rst_hdr;
  ck_impl_xs_cyc = static_cast<int>(h.xs_cyc);
  if (CkRstStaged(kCkSlabXsCyc)) {
    ck_impl_xs_cyc = static_cast<int>(ck_rst_stage[kCkSlabXsCyc][0]);
    CkRstStageFree(kCkSlabXsCyc);
  }
  ck_xs_bdt = static_cast<Real>(h.xs_bdt);
  if (CkRstStaged(kCkSlabXsBdt)) {
    ck_xs_bdt = ck_rst_stage[kCkSlabXsBdt][0];
    CkRstStageFree(kCkSlabXsBdt);
  }
  ck_impl_jac_built = (h.jac_built != 0);
  ck_impl_nstore = h.nstore;
  ck_impl_nreuse = h.nreuse;
  ck_impl_nsweep = h.nsweep;
  ck_cad_nfull = h.nfull;
  ck_cad_nlin = h.nlin;
  ck_cad_nguard = h.nguard;
  ck_cad_fref = h.fref;
}

//----------------------------------------------------------------------------------------
//! \fn void CkRstAfterAlloc
//! \brief right after CkImplAlloc: ck_thk (the snapshot's during a rebuild) and ck_cv.

inline void CkRstAfterAlloc() {
  // a restarted run whose first implicit call is a cadence GUARD call: that call's
  // ck_done = mask copy ran before ck_done existed, so it is made here, or the masked
  // columns (which already took the linearised step) would be solved as well
  if (ck_cad_partial && ck_done_ptr != nullptr && ck_cad_mask_ptr != nullptr) {
    Kokkos::deep_copy(*ck_done_ptr, *ck_cad_mask_ptr);
  }
  if (!ck_rst_have) return;
  if (ck_rst_rebuild_active) {
    CkRstPut(kCkSlabXsK, *ck_thk_ptr);
  } else {
    CkRstPut(kCkSlabThk, *ck_thk_ptr);
  }
  if (ck_cv_ptr != nullptr) CkRstPut(kCkSlabCv, *ck_cv_ptr);
}

//----------------------------------------------------------------------------------------
//! \fn void CkRstRebuild
//! \brief re-run pass 0 of the storing call on the state it saw (the XS slabs), which
//! stores exactly the operator (and cut, frozen opacity, beam, ...) that call stored;
//! then put the restarted state back.  Nothing of the gas is kept from the rebuild.

inline void CkRstRebuild(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  const bool hyd = (pmbp->phydro != nullptr);
  DvceArray5D<Real> u0 = hyd ? pmbp->phydro->u0 : pmbp->pmhd->u0;
  DvceArray4D<Real> wt = hyd ? pmbp->phydro->wtemp : pmbp->pmhd->wtemp;
  DvceArray5D<Real> bc;
  if (!hyd) bc = pmbp->pmhd->bcc0;
  using Kokkos::ALL;
  // keep the restarted state
  DvceArray5D<Real> u0k("ck_rst_u0k", u0.extent(0), u0.extent(1), u0.extent(2),
                        u0.extent(3), u0.extent(4));
  Kokkos::deep_copy(u0k, u0);
  DvceArray4D<Real> wtk("ck_rst_wtk", wt.extent(0), wt.extent(1), wt.extent(2),
                        wt.extent(3));
  if (wt.size() > 0) Kokkos::deep_copy(wtk, wt);
  DvceArray5D<Real> bck;
  if (!hyd) {
    bck = DvceArray5D<Real>("ck_rst_bck", bc.extent(0), bc.extent(1), bc.extent(2),
                            bc.extent(3), bc.extent(4));
    Kokkos::deep_copy(bck, bc);
  }
  // the storing call's state
  const int nv = static_cast<int>(u0.extent(1));
  for (int n=0; n<nv && n<16; ++n) {
    CkRstPut(kCkSlabXsU0 + n, Kokkos::subview(u0, ALL, n, ALL, ALL, ALL));
  }
  if (wt.size() > 0) {
    CkRstPut(kCkSlabXsW, wt);
  } else {
    CkRstStageFree(kCkSlabXsW);
  }
  if (!hyd) {
    for (int n=0; n<3; ++n) {
      CkRstPut(kCkSlabXsB0 + n, Kokkos::subview(bc, ALL, n, ALL, ALL, ALL));
    }
  }
  if (ck_thk_ptr != nullptr) CkRstPut(kCkSlabXsK, *ck_thk_ptr);   // else: AfterAlloc
  const int xs_keep = ck_impl_xs_cyc;
  const std::int64_t nst = ck_impl_nstore, nre = ck_impl_nreuse;
  const bool jb = ck_impl_jac_built;
  ck_impl_xs_cyc = -1;                  // pass 0 stores
  ck_impl_pass = 0;
  ck_impl_prev_res = -1.0;
  ck_impl_jac_again = false;
  if (ck_done_ptr != nullptr) Kokkos::deep_copy(*ck_done_ptr, 0.0);
  ck_rst_rebuild_active = true;
  picket_fence_two_stream_RT_pass(pm, ck_xs_bdt);
  ck_rst_rebuild_active = false;
  ck_impl_pass = -1;
  ck_impl_xs_cyc = xs_keep;
  ck_impl_nstore = nst;
  ck_impl_nreuse = nre;
  ck_impl_jac_built = jb;
  CkRstStageFree(kCkSlabXsK);
  // the restarted state back, and the thin/thick verdict the last call ended on
  Kokkos::deep_copy(u0, u0k);
  if (wt.size() > 0) Kokkos::deep_copy(wt, wtk);
  if (!hyd) Kokkos::deep_copy(bc, bck);
  CkRstPut(kCkSlabThk, *ck_thk_ptr);
  if (global_variable::my_rank == 0) {
    std::cout << "ck restart: rebuilt the ck_impl_xstep operator stored on cycle "
              << xs_keep << " (bdt " << ck_xs_bdt << ")" << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CkRstBeginCall
//! \brief the top of every implicit call (picket_fence_two_stream_RT): on the first one
//! of a restarted run, restore the scalars and, if this call is going to re-apply the
//! stored operator, rebuild it.

inline void CkRstBeginCall(Mesh *pm) {
  CkRstScalars();
  if (!ck_rst_have || ck_rst_first_done) return;
  ck_rst_first_done = true;
  if (ck_thk_ptr != nullptr) {        // already allocated: deliver now
    CkRstPut(kCkSlabThk, *ck_thk_ptr);
    if (ck_cv_ptr != nullptr) CkRstPut(kCkSlabCv, *ck_cv_ptr);
  }
  const bool snap = CkRstStaged(kCkSlabXsK);
  const bool reuse = (ck_impl_xstep > 0) && (ck_impl_xs_cyc >= 0) &&
                     (pm->ncycle - ck_impl_xs_cyc < ck_impl_xstep) && !ck_cad_partial;
  if (reuse && snap) {
    CkRstRebuild(pm);
  } else if (reuse) {
    ck_impl_xs_cyc = -1;              // no operator to re-apply: store afresh
    if (global_variable::my_rank == 0) {
      std::cout << "### WARNING: restart file has no ck_impl_xstep snapshot; the first "
                << "call stores a new operator and this restart is not bitwise."
                << std::endl;
    }
  }
  if (!reuse || !snap) {
    for (int id=kCkSlabXsW; id<kCkRstMaxSlab; ++id) CkRstStageFree(id);
  }
  // once the first call has run, nothing staged is wanted any more
  if (!CkRstStaged(kCkSlabThk) && !CkRstStaged(kCkSlabCv) &&
      !CkRstStaged(kCkSlabCad0)) {
    ck_rst_have = false;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CkRstCadRestore
//! \brief the top of CkCadStep: the cadence's Q0, D, T0, rho0, so that the first step
//! of a restarted run follows the schedule instead of a forced full call.

inline void CkRstCadRestore(Mesh *pm) {
  CkRstScalars();
  if (!ck_rst_have || ck_cad_ptr != nullptr || !CkRstStaged(kCkSlabCad0)) return;
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int n1 = indcs.nx1 + 2*indcs.ng;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
  ck_cad_ptr = new DvceArray5D<Real>("ck_cad", pmbp->nmb_thispack, 6, n3, n2, n1);
  ck_cad_mask_ptr = new DvceArray3D<Real>("ck_cad_mask", pmbp->nmb_thispack, n3, n2);
  Kokkos::deep_copy(*ck_cad_mask_ptr, 2.0);
  bool ok = true;
  for (int s=0; s<4; ++s) {
    ok = CkRstPut(kCkSlabCad0 + s, Kokkos::subview(*ck_cad_ptr, Kokkos::ALL, s,
                  Kokkos::ALL, Kokkos::ALL, Kokkos::ALL)) && ok;
  }
  if (!ok) {                          // does not fit: back to the forced full call
    delete ck_cad_ptr;
    delete ck_cad_mask_ptr;
    ck_cad_ptr = nullptr;
    ck_cad_mask_ptr = nullptr;
  }
}

}  // namespace two_stream_rt

#endif  // UTILS_TWO_STREAM_CK_RST_STATE_HPP_
