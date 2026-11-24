//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_tasks.cpp
//! \brief rad_vet task lists, the thin task wrappers and the time step.

#include <float.h>

#include <algorithm>
#include <limits>
#include <map>
#include <memory>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "coordinates/coordinates.hpp"
#include "driver/driver.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"

namespace radvet {
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

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationVet::NewTimeStep

TaskStatus RadiationVet::NewTimeStep(Driver *pdriver, int stage) {
  // no light-speed CFL unless implicit_cfl > 0 (dt = implicit_cfl dx/c)
  dtnew = static_cast<Real>(std::numeric_limits<float>::max());
  if (!(impl_cfl > 0.0)) return TaskStatus::complete;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  int nx3 = indcs.nx3;

  Real dt1 = std::numeric_limits<float>::max();
  Real dt2 = std::numeric_limits<float>::max();
  Real dt3 = std::numeric_limits<float>::max();

  auto &size = pmy_pack->pmb->mb_size;
  const int nmkji = (pmy_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;

  Kokkos::parallel_reduce("RadVetNudt", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &min_dt1, Real &min_dt2, Real &min_dt3) {
    int m = (idx)/nkji;
    min_dt1 = fmin((size.d_view(m).dx1), min_dt1);
    min_dt2 = fmin((size.d_view(m).dx2), min_dt2);
    min_dt3 = fmin((size.d_view(m).dx3), min_dt3);
  }, Kokkos::Min<Real>(dt1), Kokkos::Min<Real>(dt2), Kokkos::Min<Real>(dt3));

  Real dxmin = dt1;
  if (pmy_pack->pmesh->multi_d) {dxmin = std::min(dxmin, dt2);}
  if (pmy_pack->pmesh->three_d) {dxmin = std::min(dxmin, dt3);}
  if (sph_geom) {
    // spherical: mb_size dx2/dx3 are angles, so take the smallest physical width
    // (dr, r dtheta, r sin(theta) dphi) over the active cells
    auto cx1 = pmy_pack->pcoord->dx1;
    auto cx2 = pmy_pack->pcoord->dx2;
    auto cx3 = pmy_pack->pcoord->dx3;
    const int is = indcs.is, js = indcs.js, ks = indcs.ks;
    const int nji = nx2*nx1;
    const bool md = pmy_pack->pmesh->multi_d, td = pmy_pack->pmesh->three_d;
    Real smin = std::numeric_limits<float>::max();
    Kokkos::parallel_reduce("RadVetNudtSph",
                            Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int &idx, Real &lmin) {
      int m = idx/nkji;
      int r = idx - m*nkji;
      int k = r/nji;
      r -= k*nji;
      int j = r/nx1;
      int i = r - j*nx1;
      k += ks; j += js; i += is;
      Real d = cx1(m,k,j,i);
      if (md) {d = fmin(d, cx2(m,k,j,i));}
      if (td) {d = fmin(d, cx3(m,k,j,i));}
      lmin = fmin(d, lmin);
    }, Kokkos::Min<Real>(smin));
    dxmin = smin;
  }

  dtnew = impl_cfl*dxmin/c_light;
  return TaskStatus::complete;
}

} // namespace radvet
