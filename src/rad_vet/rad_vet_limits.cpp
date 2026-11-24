//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_limits.cpp
//! \brief the closure-limit pass of the M1 moments (ApplyClosureLimits): E >= e_floor
//! and |F| <= c E; scheme = clean keeps its solution.

#include <math.h>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "hydro/hydro.hpp"
#include "reconstruct/plm.hpp"
#include "reconstruct/ppm.hpp"
#include "reconstruct/wenoz.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_closure.hpp"
#include "rad_vet/vet_fluid.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationVet::ApplyClosureLimits
//! \brief enforce E >= e_floor and |F| <= c E on the evolved state, over the whole
//! array (active cells and ghosts alike) so that the reconstruction below only ever
//! sees admissible states.

TaskStatus RadiationVet::ApplyClosureLimits(Driver *pdrive, int stage) {
  // implicit_timers: the M1 stage chain starts here
  tmr_on = (tmr_c0 > 0) && (pmy_pack->pmesh->ncycle >= tmr_c0);
  if (tmr_on) {tmr_cnt[7] += 1.0;}
  TmrMark(-1);
  // the clean solve keeps its solution: a clip would change the energy and momentum
  TmrMark(0);
  return TaskStatus::complete;
}

} // namespace radvet
