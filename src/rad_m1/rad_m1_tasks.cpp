//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_tasks.cpp
//! \brief assembly of the rad_m1 task lists and the thin task wrappers.
//!
//! The module runs in its OWN sub-cycled task lists, "m1_before_stagen" / "m1_stagen" /
//! "m1_after_stagen", which the Driver executes after "after_timeintegrator" (design
//! sect. 5, modelled on the resistive RKG loop).  Nothing is added to the lists of the
//! main time integrator, so a run without a <rad_m1> block is bitwise unchanged.

#include <map>
#include <memory>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

namespace radm1 {
//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::AssembleRadM1Tasks
//! \brief Adds rad_m1 tasks to the m1_* task lists.  Called by MeshBlockPack::AddPhysics.

void RadiationM1::AssembleRadM1Tasks(std::map<std::string,
                                     std::shared_ptr<TaskList>> tl) {
  TaskID none(0);

  // "m1_before_stagen": work that must complete over all MeshBlocks before the stage
  id.irecv = tl["m1_before_stagen"]->AddTask(&RadiationM1::InitRecv, this, none);

  // ---- MILESTONE 3a hook: the implicit column solve replaces the whole explicit stage
  // chain (closure limits -> opacity -> fluxes -> update -> coupling).  Everything after
  // the solve -- the hydro inversion and the ghost-zone exchange of the moments -- is
  // the same as in the explicit path.
  if (transport >= M1_TRANSPORT_IMPLICIT_X1) {
    id.closure = tl["m1_stagen"]->AddTask(&RadiationM1::ApplyClosureLimits,this,none);
    id.opac    = tl["m1_stagen"]->AddTask(&RadiationM1::Opacity, this, id.closure);
    id.update  = tl["m1_stagen"]->AddTask(&RadiationM1::ImplicitSolve, this, id.opac);
    id.c2p     = tl["m1_stagen"]->AddTask(&RadiationM1::HydroConToPrim,this,id.update);
    id.restu   = tl["m1_stagen"]->AddTask(&RadiationM1::RestrictU, this, id.c2p);
    id.sendu   = tl["m1_stagen"]->AddTask(&RadiationM1::SendU, this, id.restu);
    id.recvu   = tl["m1_stagen"]->AddTask(&RadiationM1::RecvU, this, id.sendu);
    id.bcs     = tl["m1_stagen"]->AddTask(&RadiationM1::ApplyPhysicalBCs,this,id.recvu);
    id.prol    = tl["m1_stagen"]->AddTask(&RadiationM1::Prolongate, this, id.bcs);
    id.csend = tl["m1_after_stagen"]->AddTask(&RadiationM1::ClearSend, this, none);
    id.crecv = tl["m1_after_stagen"]->AddTask(&RadiationM1::ClearRecv, this, id.csend);
    return;
  }
  // ---- end of the 3a hook

  // "m1_stagen"
  id.copyu   = tl["m1_stagen"]->AddTask(&RadiationM1::CopyCons, this, none);
  id.closure = tl["m1_stagen"]->AddTask(&RadiationM1::ApplyClosureLimits,this,id.copyu);
  id.opac    = tl["m1_stagen"]->AddTask(&RadiationM1::Opacity, this, id.closure);
  id.flux    = tl["m1_stagen"]->AddTask(&RadiationM1::CalculateFluxes, this, id.opac);
  id.sendf   = tl["m1_stagen"]->AddTask(&RadiationM1::SendFlux, this, id.flux);
  id.recvf   = tl["m1_stagen"]->AddTask(&RadiationM1::RecvFlux, this, id.sendf);
  id.update  = tl["m1_stagen"]->AddTask(&RadiationM1::Update, this, id.recvf);
  id.coupl   = tl["m1_stagen"]->AddTask(&RadiationM1::Coupling, this, id.update);
  id.c2p     = tl["m1_stagen"]->AddTask(&RadiationM1::HydroConToPrim,this,id.coupl);
  id.restu   = tl["m1_stagen"]->AddTask(&RadiationM1::RestrictU, this, id.c2p);
  id.sendu   = tl["m1_stagen"]->AddTask(&RadiationM1::SendU, this, id.restu);
  id.recvu   = tl["m1_stagen"]->AddTask(&RadiationM1::RecvU, this, id.sendu);
  id.bcs     = tl["m1_stagen"]->AddTask(&RadiationM1::ApplyPhysicalBCs,this,id.recvu);
  id.prol    = tl["m1_stagen"]->AddTask(&RadiationM1::Prolongate, this, id.bcs);

  // "m1_after_stagen"
  id.csend = tl["m1_after_stagen"]->AddTask(&RadiationM1::ClearSend, this, none);
  id.crecv = tl["m1_after_stagen"]->AddTask(&RadiationM1::ClearRecv, this, id.csend);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationM1::InitRecv
//! \brief post non-blocking receives (with MPI) and reset boundary receive flags

TaskStatus RadiationM1::InitRecv(Driver *pdrive, int stage) {
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
//! \fn TaskStatus RadiationM1::CopyCons
//! \brief copy u0 --> u1 in the first stage of each substep

TaskStatus RadiationM1::CopyCons(Driver *pdrive, int stage) {
  if (stage == 1) {
    Kokkos::deep_copy(DevExeSpace(), u1, u0);
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationM1::SendFlux
//! \brief pack/send restricted face fluxes at fine/coarse boundaries

TaskStatus RadiationM1::SendFlux(Driver *pdrive, int stage) {
  TaskStatus tstat = TaskStatus::complete;
  if (pmy_pack->pmesh->multilevel) {
    tstat = pbval_u->PackAndSendFluxCC(uflx);
  }
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationM1::RecvFlux
//! \brief recv/unpack restricted face fluxes at fine/coarse boundaries

TaskStatus RadiationM1::RecvFlux(Driver *pdrive, int stage) {
  TaskStatus tstat = TaskStatus::complete;
  if (pmy_pack->pmesh->multilevel) {
    tstat = pbval_u->RecvAndUnpackFluxCC(uflx);
  }
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationM1::RestrictU

TaskStatus RadiationM1::RestrictU(Driver *pdrive, int stage) {
  if (pmy_pack->pmesh->multilevel) {
    pmy_pack->pmesh->pmr->RestrictCC(u0, coarse_u0);
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationM1::SendU

TaskStatus RadiationM1::SendU(Driver *pdrive, int stage) {
  return pbval_u->PackAndSendCC(u0, coarse_u0);
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationM1::RecvU

TaskStatus RadiationM1::RecvU(Driver *pdrive, int stage) {
  return pbval_u->RecvAndUnpackCC(u0, coarse_u0);
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationM1::ApplyPhysicalBCs

TaskStatus RadiationM1::ApplyPhysicalBCs(Driver *pdrive, int stage) {
  // do not apply BCs if domain is strictly periodic
  if (pmy_pack->pmesh->strictly_periodic) return TaskStatus::complete;

  pbval_u->RadM1BCs((pmy_pack), (pbval_u->u_in), u0, c_light, e_floor);

  // user BCs (the beam test patches the inflow face here)
  if (pmy_pack->pmesh->pgen->user_bcs) {
    (pmy_pack->pmesh->pgen->user_bcs_func)(pmy_pack->pmesh);
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationM1::Prolongate

TaskStatus RadiationM1::Prolongate(Driver *pdrive, int stage) {
  if (pmy_pack->pmesh->multilevel) {
    pbval_u->FillCoarseInBndryCC(u0, coarse_u0);
    pbval_u->ProlongateCC(u0, coarse_u0);
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationM1::ClearSend

TaskStatus RadiationM1::ClearSend(Driver *pdrive, int stage) {
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
//! \fn TaskStatus RadiationM1::ClearRecv

TaskStatus RadiationM1::ClearRecv(Driver *pdrive, int stage) {
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

} // namespace radm1
