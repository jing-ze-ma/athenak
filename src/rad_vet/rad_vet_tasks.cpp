//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_tasks.cpp
//! \brief assembly of the rad_vet task lists and the thin task wrappers.
//!
//! The module runs in its OWN sub-cycled task lists, "vet_before_stagen" / "vet_stagen" /
//! "vet_after_stagen", which the Driver executes after "after_timeintegrator" (design
//! sect. 5, modelled on the resistive RKG loop).  Nothing is added to the lists of the
//! main time integrator, so a run without a <rad_vet> block is bitwise unchanged.

#include <map>
#include <memory>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"

namespace radvet {
//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::AssembleRadVetTasks
//! \brief Adds rad_vet tasks to the m1_* task lists.  Called by
//! MeshBlockPack::AddPhysics.

void RadiationVet::AssembleRadVetTasks(std::map<std::string,
                                       std::shared_ptr<TaskList>> tl) {
  TaskID none(0);

  // "vet_before_stagen": work that must complete over all MeshBlocks before the stage
  id.irecv = tl["vet_before_stagen"]->AddTask(&RadiationVet::InitRecv, this, none);

  // "vet_stagen": the clean stage solve (scheme = clean, transport = implicit), then the
  // hydro inversion and the ghost-zone exchange of the moments
  id.closure = tl["vet_stagen"]->AddTask(&RadiationVet::ApplyClosureLimits,this,none);
  id.opac    = tl["vet_stagen"]->AddTask(&RadiationVet::Opacity, this, id.closure);
  id.update  = tl["vet_stagen"]->AddTask(&RadiationVet::CleanSolve, this, id.opac);
  id.c2p     = tl["vet_stagen"]->AddTask(&RadiationVet::HydroConToPrim,this,id.update);
  id.restu   = tl["vet_stagen"]->AddTask(&RadiationVet::RestrictU, this, id.c2p);
  id.sendu   = tl["vet_stagen"]->AddTask(&RadiationVet::SendU, this, id.restu);
  id.recvu   = tl["vet_stagen"]->AddTask(&RadiationVet::RecvU, this, id.sendu);
  id.bcs     = tl["vet_stagen"]->AddTask(&RadiationVet::ApplyPhysicalBCs,this,id.recvu);
  id.prol    = tl["vet_stagen"]->AddTask(&RadiationVet::Prolongate, this, id.bcs);
  // "vet_after_stagen"
  id.csend = tl["vet_after_stagen"]->AddTask(&RadiationVet::ClearSend, this, none);
  id.crecv = tl["vet_after_stagen"]->AddTask(&RadiationVet::ClearRecv, this, id.csend);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationVet::InitRecv
//! \brief post non-blocking receives (with MPI) and reset boundary receive flags

TaskStatus RadiationVet::InitRecv(Driver *pdrive, int stage) {
  TaskStatus tstat = pbval_u->InitRecv(M1_NVAR);
  if (tstat != TaskStatus::complete) return tstat;

  // do not post receives for fluxes when stage < 0 (i.e. ICs)
  if (stage >= 0) {
    if (pmy_pack->pmesh->multilevel) {
      tstat = pbval_u->InitFluxRecv(M1_NVAR);
      if (tstat != TaskStatus::complete) return tstat;
    }
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationVet::RestrictU

TaskStatus RadiationVet::RestrictU(Driver *pdrive, int stage) {
  if (pmy_pack->pmesh->multilevel) {
    pmy_pack->pmesh->pmr->RestrictCC(u0, coarse_u0);
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationVet::SendU

TaskStatus RadiationVet::SendU(Driver *pdrive, int stage) {
  return pbval_u->PackAndSendCC(u0, coarse_u0);
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationVet::RecvU

TaskStatus RadiationVet::RecvU(Driver *pdrive, int stage) {
  return pbval_u->RecvAndUnpackCC(u0, coarse_u0);
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationVet::ApplyPhysicalBCs

TaskStatus RadiationVet::ApplyPhysicalBCs(Driver *pdrive, int stage) {
  // do not apply BCs if domain is strictly periodic
  if (pmy_pack->pmesh->strictly_periodic) return TaskStatus::complete;

  pbval_u->RadVetBCs((pmy_pack), (pbval_u->u_in), u0, c_light, e_floor);

  // user BCs (the beam test patches the inflow face here)
  if (pmy_pack->pmesh->pgen->user_bcs) {
    (pmy_pack->pmesh->pgen->user_bcs_func)(pmy_pack->pmesh);
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationVet::Prolongate

TaskStatus RadiationVet::Prolongate(Driver *pdrive, int stage) {
  if (pmy_pack->pmesh->multilevel) {
    pbval_u->FillCoarseInBndryCC(u0, coarse_u0);
    pbval_u->ProlongateCC(u0, coarse_u0);
  }
  TmrMark(10);
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationVet::ClearSend

TaskStatus RadiationVet::ClearSend(Driver *pdrive, int stage) {
  TaskStatus tstat = pbval_u->ClearSend();
  if (tstat != TaskStatus::complete) return tstat;
  if (stage >= 0) {
    if (pmy_pack->pmesh->multilevel) {
      tstat = pbval_u->ClearFluxSend();
      if (tstat != TaskStatus::complete) return tstat;
    }
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationVet::ClearRecv

TaskStatus RadiationVet::ClearRecv(Driver *pdrive, int stage) {
  TaskStatus tstat = pbval_u->ClearRecv();
  if (tstat != TaskStatus::complete) return tstat;
  if (stage >= 0) {
    if (pmy_pack->pmesh->multilevel) {
      tstat = pbval_u->ClearFluxRecv();
      if (tstat != TaskStatus::complete) return tstat;
    }
  }
  return TaskStatus::complete;
}

} // namespace radvet
