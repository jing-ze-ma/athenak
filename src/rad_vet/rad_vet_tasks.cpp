//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_tasks.cpp
//! \brief the workflow of <rad_vet>, read this file first: the task lists, and the
//! implicit stage solve (CleanSolve) as a driver of named steps.
//!
//!  called from: mesh/meshblock_pack.cpp (task lists), driver/driver.cpp
//!  calls into: rad_vet_bvals.cpp, rad_vet_closure.cpp, rad_vet_closure_faces.cpp,
//!    rad_vet_closure_ring.cpp, rad_vet_flux_seam.cpp, rad_vet_init.cpp,
//!    rad_vet_ledger.cpp, rad_vet_newton.cpp, rad_vet_oldvector.cpp, rad_vet_rhs.cpp,
//!    rad_vet_solver.cpp, rad_vet_steady.cpp, rad_vet_update.cpp, rays/sc_source.cpp

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "tasklist/task_list.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "driver/driver.hpp"
#include "pgen/pgen.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_fluid.hpp"
#include "rad_vet/rad_vet_args.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationVet::CleanSolve
//! \brief the implicit solve of one stage (E in cells, F on faces, gas), as named steps:
//! closure -> old vector -> passes { rows, E, gas T, F and momentum } -> new state

TaskStatus RadiationVet::CleanSolve(Driver *pdriver, int stage) {
  CleanSolveBegin();                  // arrays, direction sets, first faces  (_init)

  // 1. what the stage starts from
  CleanStageClosure();                // D, h+-, face sums from the rays      (_closure)
  CleanWoBind();                      // work vectors over the rays' memory   (_init)
  cln_xdh = false;                    // no donor hold before this solve's passes
  CleanOldVector();                   // E_old and gas + stage increments     (_oldvector)
  const bool opu = cln_opt && cln_feedback;   // CleanFirstPass: both, with its opacities
  if (!opu) CleanFaceCoefficients();  // chi_f, theta, blend weight           (_solver)
  if (cln_gflux || cln_consist) CleanGhostFill(cln_w, CLW_EO, 1);   // donors' E_old
  if (cln_consist && !opu) CleanStagBracket();                      // gate q (_newton)

  // 2. passes: Newton on the gas coupling, Picard on the lagged slopes and the closure
  SolvePass p = CleanPassSetup();     // what is iterated, tolerances         (_newton)
  CleanFirstPass();                   // first iterate, drag, donors          (_newton)
  for (p.n = 0; p.n < p.nmax; ++p.n) {
    CleanPassClosure(p);              // closure, opacity, donors of E_k      (_newton)
    CleanPassRows(p);                 // right-hand side, operator        (_rhs, _solver)
    if (!CleanPassSolve(p)) break;    // E: BiCGStab, multigrid               (_solver)
    CleanPassUpdate(p);               // gas T; face F, momentum, drag  (_newton, _flux)
    if (CleanPassSettled(p)) break;   // stop test                            (_newton)
  }
  CleanSolveEnd(p);                   // counts: realizability, E < 0         (_newton)
  CleanVerdictSolve(p);               // unconverged: dt_retry, FATAL         (_newton)

  // 3. the new state: face fluxes, E, cell fluxes, gas, stage slopes, L_z books
  CleanUpdateState(p);                // CleanWriteBack                      (_update)
  CleanVerdictState(p);               // gas eint: t2_fail, dt_retry          (_update)
  cln_xd_on = false;                  // the donor q is of this solve's passes only

  // 4. clean_closure = extrap(2): ray moments of the new state, for the next closure
  if (cln_extrap > 0 && !t2_fail) CleanClosureAfterSolve();         //        (_closure)
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::AssembleRadVetTasks
//! \brief adds the rad_vet tasks to the vet_* task lists (called by AddPhysics)

void RadiationVet::AssembleRadVetTasks(std::map<std::string,
                                       std::shared_ptr<TaskList>> tl) {
  TaskID none(0);

  // "vet_before_stagen"
  id.irecv = tl["vet_before_stagen"]->AddTask(&RadiationVet::InitRecv, this, none);

  // "vet_stagen": clean stage solve, hydro inversion, ghost exchange of the moments
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
  TaskStatus tstat = pbval_u->InitRecv(VET_NVAR);
  if (tstat != TaskStatus::complete) return tstat;

  // do not post receives for fluxes when stage < 0 (i.e. ICs)
  if (stage >= 0) {
    if (pmy_pack->pmesh->multilevel) {
      tstat = pbval_u->InitFluxRecv(VET_NVAR);
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

  // user BCs
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

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationVet::ApplyClosureLimits
//! \brief first task of the stage chain: starts the implicit_timers; applies no clip

TaskStatus RadiationVet::ApplyClosureLimits(Driver *pdrive, int stage) {
  tmr_on = (tmr_c0 > 0) && (pmy_pack->pmesh->ncycle >= tmr_c0);
  if (tmr_on) {tmr_cnt[7] += 1.0;}
  TmrMark(-1);
  // no clip: it would break energy and momentum conservation
  TmrMark(0);
  return TaskStatus::complete;
}

}  // namespace radvet
