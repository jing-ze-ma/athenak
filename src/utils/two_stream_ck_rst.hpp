#ifndef UTILS_TWO_STREAM_CK_RST_HPP_
#define UTILS_TWO_STREAM_CK_RST_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file two_stream_ck_rst.hpp
//! \brief THE IMPLICIT CORRELATED-K SOLVER'S CROSS-CALL STATE AS RESTART STATE (the
//! light half: file format, header, staging; restart.cpp and pgen.cpp include this one).
//! The heavy half -- what is collected, how it is put back, and the rebuild of the
//! ck_impl_xstep operator -- is utils/two_stream_ck_rst_state.hpp, included by
//! two_stream_rt.hpp.  See tests_ck_implicit/README_restart.md.
//!
//! What persists from one RT call to the next under problem/ck_implicit (T4), the ck-fast
//! levers (c2) and the cadence (ck_impl_every), and how each piece travels:
//!   * ck_thk, the thin/thick verdict of the last apply, read by pass 0 of the next call
//!     (rt_pre_opac)                                                   -> slab THK
//!   * the cadence's Q0, D, T0, rho0 (ck_cad slots 0-3)                 -> slabs CAD0-3
//!   * ck_cv, the secant cv the call ended on (ck_impl_cvkeep only)     -> slab CV
//!   * the ck_impl_xstep STORED OPERATOR (kappa rho, coefficient triple, tm factorisation,
//!     cut, frozen beam, ...): O(100) Reals per cell, far too large to write.  It is a
//!     deterministic function of the state the storing call's pass 0 saw, so THAT state
//!     is written (u0, wtemp, bcc0 for MHD, and ck_thk at the start of that call) with
//!     the call's bdt, and the first call after the restart re-runs that pass 0 on it
//!     before it starts (CkRstRebuild).  Only written when the next call can re-apply
//!     the operator (xs_cyc >= 0 and ncycle - xs_cyc < xstep).       -> slabs XS*
//!   * ck_impl_xs_cyc per rank (a rank-local quantity under ck_impl_xstep_thr and the
//!     cadence guard)                                                   -> slab XSCYC
//!   * host scalars: ck_impl_jac_built, the xstep store bdt, and the run counters that
//!     the verbose lines print                                          -> header
//!   * problem/flux_hst_rst: the face-flux cache of the flux history (diagnostic only;
//!     restored into rt_Fb/rt_Qb chain block 0 as block sums)          -> slabs FH*
//!   * problem/ck_sph_face = 5: the VEF factors (O(100) Reals per cell) are, like the
//!     xstep operator, re-formed from the state the last REFRESH read (u0, wtemp, bcc0,
//!     ck_thk), its bdt and its cycle ck_vef_last (CkRstVefRebuild)    -> slabs VEF*
//! Everything else the solver keeps (ck_rprev, ck_ep/ck_tp, ck_jac, ck_done, ck_dep,
//! the per-call scalars) is reset or rewritten before it is read in every call.
//! NOT covered: ck_impl_warm (ck_dep + ck_ei carried; a failed lever, refused with the
//! cadence), ck_impl_reuse_jac = 2, ck_impl_glob, ck_impl_aa.
//!
//! THE FILE FORMAT, as the other marked blocks (kPgenRstMagic, kEintRstMagic): a marked
//! header behind every other marked header (magic, IOWrapperSizeT length, CkRstHdr),
//! and nslab per-MeshBlock slabs behind every other tail slab of each MeshBlock record.
//! No block is written when there is nothing to carry (ck_implicit off), so such files
//! are byte-identical to the ones written before this existed.  A file without the block
//! (older, or ck_implicit off when it was written) restarts as before: ck_thk all thin,
//! no stored operator, and under the cadence a full call first.

#include <cstdint>
#include <cstring>
#include <vector>

#include "athena.hpp"

class Mesh;

namespace two_stream_rt {

constexpr char kCkRstMagic[8] = {'C', 'K', 'R', 'S', 'T', '0', '0', '1'};
constexpr int kCkRstMaxSlab = 48;

// slab ids
enum CkRstSlab : std::int32_t {
  kCkSlabThk = 1,
  kCkSlabCad0 = 2,    // Q0, D, T0, rho0 = 2 .. 5
  kCkSlabCv = 6,
  kCkSlabXsCyc = 7,   // ck_impl_xs_cyc of the rank, as a Real, in every cell
  kCkSlabXsW = 8,     // wtemp at the store
  kCkSlabXsK = 9,     // ck_thk at the start of the storing call
  kCkSlabXsBdt = 10,  // the storing call's bdt, per rank, in every cell
  // problem/flux_hst_rst (rthst-1101): the radiative face-flux cache the flux history
  // columns read (rt_Fb, rt_Qb, rt_icut, rt_cf slot 3, rt_hist_fstar/albedo), as the
  // last RT call left it.  Without it the first rows after a restart, written before
  // the cadence's next full call, had no radiative part (see RtFhCollect).
  kCkSlabFhFb = 11,   // rt_Fb summed over the chain blocks, (m,k,j,i), faces i
  kCkSlabFhQb = 12,   // rt_Qb summed over the chain blocks
  kCkSlabFhCol = 13,  // per column, i = 0..3: rt_icut, rt_cf(3), fstar, albedo
  // ckhesdirk-1102 problem/ck_time_scheme = hesdirk2: the FSAL slope K1 and its
  // validity (a Real, 1 or 0, in every cell)
  kCkSlabT2K1 = 14,
  kCkSlabT2Ok = 15,
  kCkSlabXsU0 = 16,   // u0 variable n at the store = 16 + n (n < 10)
  // ck2ndt-1102 problem/ck_cad_cont (utils/ck_2nd.hpp (8b)): the continuous model's
  // ck_cad2 slots G (0), t0 (1), S_n (6), C (7) = 26 .. 29 and its ramp length L = 30
  // (ids 26..31 are XsU0 + 10.. that no gas variable count reaches)
  kCkSlabCd2 = 26,
  kCkSlabCd2L = 30,
  kCkSlabXsB0 = 32,   // bcc0 component n at the store = 32 + n (MHD)
  // problem/ck_sph_face = 5: the state the last refresh of the VEF factors read, from
  // which a restart re-forms them (CkRstVefRebuild)
  kCkSlabVefCyc = 35,  // ck_vef_last, as a Real, in every cell
  kCkSlabVefBdt = 36,  // the refreshing call's bdt
  kCkSlabVefW = 37,    // wtemp
  kCkSlabVefK = 38,    // ck_thk at the start of the refreshing call
  kCkSlabVefB0 = 39,   // bcc0 component n = 39 + n (MHD)
  kCkSlabVefU0 = 42    // u0 variable n = 42 + n (n < 6)
};

struct CkRstHdr {
  std::int32_t version;
  std::int32_t nslab;
  std::int32_t id[kCkRstMaxSlab];
  std::int64_t xs_cyc;       // rank 0's; the XSCYC slab has every rank's
  std::int64_t jac_built;
  std::int64_t nstore, nreuse, nsweep, nfull, nlin, nguard;
  double fref;
  double xs_bdt;             // the bdt of the storing call (the XS slabs)
};

// what restart.cpp calls to collect the slabs (nullptr: the solver is not linked in).
// Fills out (nslab, nmb, n3, n2, n1) and the header; returns nslab (0 = no block).
using CkRstCollectFn = int (*)(Mesh *pm, HostArray5D<Real> &out, CkRstHdr &hdr,
                               int nmb, int n3, int n2, int n1);
inline CkRstCollectFn ck_rst_collect_fn = nullptr;

// problem/flux_hst_rst (deep_hot_jupiter_rt.cpp, default flux_hst && ck_implicit):
// carry the flux history's face-flux cache in this block (slabs kCkSlabFh*)
inline bool rt_fhst_rst = false;

// the restart staging, filled by pgen.cpp, consumed by the first RT call(s).  Plain
// std::vectors (not Views) so that an unconsumed staging does not outlive Kokkos.
inline bool ck_rst_restarted = false;     // a restart was read (pgen.cpp)
inline bool ck_rst_have = false;          // the file carried a block
inline CkRstHdr ck_rst_hdr;
inline std::vector<std::vector<Real>> ck_rst_stage;   // indexed by slab id
inline int ck_rst_nmb = 0, ck_rst_n3 = 0, ck_rst_n2 = 0, ck_rst_n1 = 0;

inline void CkRstStageFree(const int id) {
  if (id >= 0 && id < static_cast<int>(ck_rst_stage.size())) {
    std::vector<Real>().swap(ck_rst_stage[id]);
  }
}

inline bool CkRstStaged(const int id) {
  return id >= 0 && id < static_cast<int>(ck_rst_stage.size()) &&
         !ck_rst_stage[id].empty();
}

}  // namespace two_stream_rt

#endif  // UTILS_TWO_STREAM_CK_RST_HPP_
