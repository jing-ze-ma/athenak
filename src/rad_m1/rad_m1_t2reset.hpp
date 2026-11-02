#ifndef RADIATION_M1_RAD_M1_T2RESET_HPP_
#define RADIATION_M1_RAD_M1_T2RESET_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_t2reset.hpp
//! \brief csacc-1031 hesdirk2 RESET-SECTOR exemption (interface between a problem generator
//! and the rad_m1 stage admissibility check).  A pgen that OVERWRITES the gas (and maybe
//! the radiation) of some cells after every stage -- the accretor's stream nozzle, donor
//! mask and back sink -- allocates g_t2_reset (m, k, j, i; the hold_m extents) and sets
//! g_t2_reset_on; every stage it writes 1 in the cells it resets (0 elsewhere) BEFORE the
//! implicit solve.  The admissibility minimum of a hesdirk2 stage (solved E, solved T,
//! written-back gas eint) and the time2_posfix marking then skip those cells: their
//! values are replaced by the pgen anyway (gas at the floor or the nozzle state; the E of
//! a dark sink cell is floored by the write-back as in a backward-Euler step), so a
//! non-positive stage value there is no reason to redo the whole step with backward Euler.
//! Off (g_t2_reset_on false, nothing allocated): bitwise.  Skipped negatives are counted
//! (t2_reset_nskip, printed at the end).  Works on spherical polar and the cubed sphere
//! alike (a per-cell mask; the pgen's own geometry decides the cells).

#include "athena.hpp"

namespace radm1 {
extern bool g_t2_reset_on;
extern DvceArray4D<Real> g_t2_reset;
extern double g_t2_reset_nskip;   // admissibility values <= 0 skipped (all stages, rank)
}  // namespace radm1

#endif  // RADIATION_M1_RAD_M1_T2RESET_HPP_
