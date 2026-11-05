//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hydro_splitbval.cpp
//! \brief ckbval-1102 <hydro>/split_bval_fast: the check that lets the ghost refresh after
//! the operator-split source (Hydro::RTOpSplitBvals) exchange u0(IEN) alone.

#include <type_traits>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "driver/driver.hpp"
#include "pgen/pgen.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace hydro {
// share-1104 P6: the test and the snapshot copy are shared with MHD
#define SBV_CLASS Hydro
#define SBV_NVAR nhydro
#include "hydro/splitbval_impl.hpp"
#undef SBV_CLASS
#undef SBV_NVAR

//----------------------------------------------------------------------------------------
//! \fn void Hydro::SplitBvalSnap
//! \brief u0 at the input of the first ConToPrim after a ghost refresh (called from
//! ConToPrim when sbf_want), i.e. right after the exchange and the physical BCs.  Kept
//! only for the conversions a split refresh can follow: the last RK stage, or the split
//! refresh itself (no driver); any other refresh leaves no valid snapshot.

void Hydro::SplitBvalSnap(Driver *pdrive, int stage) {
  sbf_want = false;
  if (pdrive == nullptr || (!c2p_split && stage == pdrive->nexp_stages)) {
    SbfSnapCopy(sbf_s, u0, sbf_valid);
  }
}

//----------------------------------------------------------------------------------------
//! \fn bool Hydro::SplitC2POnce
//! \brief <hydro>/split_c2p_once: true when the ConToPrimSplit tasks are redundant.
//! Each ConToPrimSplit converts a state that has already been converted with nothing
//! writing u0 in between:
//!   - user_split_once: "after_timeintegrator" = RTStrangSplit (source, then
//!     RTOpSplitBvals = exchange, BCs, ConToPrim over the whole mesh) -> ConToPrimSplit;
//!     the next "before_timeintegrator" = (RTStrangSplit returns at stage 0) ->
//!     ConToPrimSplit;
//!   - no split source at all (e.g. <problem>/ck_time_scheme = hesdirk2, ck inside the
//!     RK stages): both lists hold only the ConToPrimSplit after the last stage's own
//!     ConToPrim.
//! Between those conversions and the next RK stage nothing writes u0 or w0 (NewTimeStep
//! and the outputs only read them).  The repeats are not bitwise no-ops -- each one redoes
//! u0(IEN) = (E - rho Phi) + rho Phi under etotgrav and re-solves T from the previous T
//! -- so this is a ROUND-OFF change, off unless named.  The Strang split (a split source
//! without user_split_once) keeps both tasks.

bool Hydro::SplitC2POnce() const {
  if (!scp_once) return false;
  const Mesh *pm = pmy_pack->pmesh;
  if (pm->pgen == nullptr || pm->pgen->user_split_func == nullptr) return true;
  return pm->pgen->user_split_once;
}

}  // namespace hydro
