//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_update.cpp
//! \brief Performs explicit update of MHD conserved variables (u0) for each stage of the
//! SSP RK integrators (e.g. RK1, RK2, RK3) implemented in AthenaK, using weighted average
//! and partial time update of flux divergence. Source terms are added in the
//! MHDSrcTerms() function.

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "mhd.hpp"
#include "hydro/rk_update_shared.hpp"
#include "dyn_grmhd/dyn_grmhd.hpp"

namespace mhd {
//----------------------------------------------------------------------------------------
//! \fn  void MHD::Update
//  \brief Explicit RK update including flux divergence terms (share-1104 P6: the stage
//  update itself is fluid::RKUpdateCore, shared with Hydro::RKUpdate)

TaskStatus MHD::RKUpdate(Driver *pdriver, int stage) {
  if (boris_on) BorisSaveRef(pdriver, stage);   // boris-1102 (mhd_boris.cpp)
  fluid::RKUpdateCore(pmy_pack, pdriver, stage, nmhd + nscalars, uflx, u0, u1);
  return TaskStatus::complete;
}
} // namespace mhd
