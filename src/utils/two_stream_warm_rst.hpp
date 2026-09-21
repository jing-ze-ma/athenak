#ifndef UTILS_TWO_STREAM_WARM_RST_HPP_
#define UTILS_TWO_STREAM_WARM_RST_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file two_stream_warm_rst.hpp
//! \brief THE MODE-3 NEWTON WARM START AS RESTART STATE.
//!
//! problem/rt_impl_warm starts the mode-3 column Newton at the PREVIOUS call's converged
//! Planck function instead of at the entry state's, held per cell in rt_c3bp (level n)
//! and, for rt_impl_warm = 2, rt_c3bp2 (level n-1).  That history is real state: the
//! iterate it sets decides which of the round-off-equivalent answers inside rt_impl_tol
//! Newton stops on, so a restart that cold starts it is NOT a bitwise continuation --
//! measured on the G1 box configuration, the first post-restart history row differs by
//! 2e-15 in the mass and 2-5e-11 in the top flux.
//!
//! This header holds the four namespace-scope variables the history needs (moved out of
//! two_stream_rt.hpp, which is far too heavy to include in the I/O code) plus the STAGING
//! the restart reader fills.  The ordering problem it solves: the reader runs inside
//! ProblemGenerator's restart constructor, long before the problem generator's first RT
//! call allocates rt_c3bp, so the bytes cannot be deep_copy'd to their destination when
//! they are read.  They are parked in rt_warm_stage_ptr / rt_warm_stage2_ptr, and the
//! lazy allocation in two_stream_rt.hpp calls RtWarmStageConsume() immediately after it
//! creates the Views: the staging is copied in and freed, and the first solve of the
//! restarted run takes exactly the code path the straight run's solve took.
//!
//! THE FILE FORMAT.  Two pieces, because the per-cell slabs and the bookkeeping have
//! different homes:
//!   * a MARKED header block (kRtWarmRstMagic), written right behind the optional pgen
//!     state block and read by the same eight-byte peek: int32 nlev, int32 pad, Real
//!     bdt_prev.  Marked rather than inferred, because the per-MeshBlock tail alone is
//!     ambiguous -- one extra slab could be a warm level or the mhd wtemp, two could be
//!     two warm levels or the hydro derived cache -- and because bdt_prev, the warm = 2
//!     extrapolation ratio's denominator, has to travel with it.
//!   * nlev per-MeshBlock slabs appended BEHIND the wtemp/wder tail, laid out exactly
//!     like them, so a reader that does not want them simply does not read them and a
//!     file written without them differs only in its tail.
//! A file with no header block has no tail either: nlev = 0, one warning line, cold
//! start, i.e. today's behaviour.

#include <cstdint>
#include <cstring>
#include <iostream>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"

namespace two_stream_rt {

//----------------------------------------------------------------------------------------
// the warm start itself.  Defined here and not in two_stream_rt.hpp so that restart.cpp
// and pgen.cpp can see them without pulling in the whole solver.

// problem/rt_impl_warm: 0 = off (bitwise the old code), 1 = the previous converged b per
// cell, 2 = linear extrapolation in time from the last two.
inline int rt_impl_warm = 0;
inline DvceArray4D<Real> *rt_c3bp_ptr = nullptr;   // the warm-start history, level n
inline DvceArray4D<Real> *rt_c3bp2_ptr = nullptr;  // the warm-start history, level n-1
inline Real rt_c3_bdt_prev = 0.0;                  // the previous call's bdt

//----------------------------------------------------------------------------------------
// the restart staging: host copies of the history levels a restart file carried, and the
// number of them (0, 1 or 2), plus the extents they were read with.  PLAIN std::vectors
// and not Kokkos Views on purpose: staging that is never consumed (a file with a history
// read by a run that does not use one) must not still be holding a View when
// Kokkos::finalize runs.

inline std::vector<Real> rt_warm_stage;
inline std::vector<Real> rt_warm_stage2;
inline int rt_warm_stage_nlev = 0;
inline int rt_warm_stage_nmb = 0;
inline int rt_warm_stage_n3 = 0;
inline int rt_warm_stage_n2 = 0;
inline int rt_warm_stage_n1 = 0;

// the eight-byte marker of the header block, exactly as long as the IOWrapperSizeT it
// stands in front of so that the reader's peek can fall back (see kPgenRstMagic)
constexpr char kRtWarmRstMagic[8] = {'R', 'T', 'W', 'A', 'R', 'M', '0', '1'};

//----------------------------------------------------------------------------------------
//! \fn int RtWarmLevels()
//! \brief how many history levels this run keeps, and therefore writes: 0 when the warm
//! start is off or has not allocated yet (a restart file written before the first RT call
//! of the run carries no history, and must not claim to).

inline int RtWarmLevels() {
  if (rt_impl_warm <= 0 || rt_c3bp_ptr == nullptr) return 0;
  return (rt_impl_warm > 1) ? 2 : 1;
}

//----------------------------------------------------------------------------------------
//! \fn void RtWarmStageFree()
//! \brief drop the staging without using it (nothing to consume it, or it does not fit).

inline void RtWarmStageFree() {
  std::vector<Real>().swap(rt_warm_stage);
  std::vector<Real>().swap(rt_warm_stage2);
  rt_warm_stage_nlev = 0;
}

//----------------------------------------------------------------------------------------
//! \fn void RtWarmStageConsume()
//! \brief copy the staged history into the freshly allocated rt_c3bp / rt_c3bp2 and free
//! it.  Called from the lazy allocation in two_stream_rt.hpp, once, right after the Views
//! exist.  A mismatch in the extents (a restart onto a different mesh decomposition,
//! which this history cannot follow) drops the staging, leaving the cold start in place.

inline void RtWarmStageConsume() {
  if (rt_warm_stage_nlev <= 0 || rt_c3bp_ptr == nullptr) {
    RtWarmStageFree();
    return;
  }
  using UnmanagedHost4D = Kokkos::View<Real****, Kokkos::LayoutRight,
                                       Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
  auto &bp = *rt_c3bp_ptr;
  const std::size_t nm = static_cast<std::size_t>(rt_warm_stage_nmb);
  if (bp.extent(0) != nm ||
      bp.extent(1) != static_cast<std::size_t>(rt_warm_stage_n3) ||
      bp.extent(2) != static_cast<std::size_t>(rt_warm_stage_n2) ||
      bp.extent(3) != static_cast<std::size_t>(rt_warm_stage_n1)) {
    if (global_variable::my_rank == 0) {
      std::cout << "### WARNING: the mode-3 warm-start history in the restart file does "
                << "not fit this decomposition; the first column solve is cold started "
                << "and this restart is not bitwise." << std::endl;
    }
    RtWarmStageFree();
    return;
  }
  UnmanagedHost4D h1(rt_warm_stage.data(), rt_warm_stage_nmb, rt_warm_stage_n3,
                     rt_warm_stage_n2, rt_warm_stage_n1);
  Kokkos::deep_copy(bp, h1);
  if (rt_warm_stage_nlev > 1 && rt_c3bp2_ptr != nullptr &&
      rt_c3bp2_ptr->extent(3) == static_cast<std::size_t>(rt_warm_stage_n1)) {
    UnmanagedHost4D h2(rt_warm_stage2.data(), rt_warm_stage_nmb, rt_warm_stage_n3,
                       rt_warm_stage_n2, rt_warm_stage_n1);
    Kokkos::deep_copy(*rt_c3bp2_ptr, h2);
  }
  RtWarmStageFree();
}

}  // namespace two_stream_rt

#endif  // UTILS_TWO_STREAM_WARM_RST_HPP_
