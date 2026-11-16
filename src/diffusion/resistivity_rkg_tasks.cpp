//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file resistivity_rkg_tasks.cpp
//! \brief functions that control the tasks of the resistive RKG super-time-stepping
//! (Runge-Kutta-Gegenbauer, second order) stored in tasklists in MeshBlockPack

#include <map>
#include <memory>
#include <string>
#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "coordinates/coordinates.hpp"
#include "eos/eos.hpp"
#include "diffusion/viscosity.hpp"
#include "diffusion/resistivity.hpp"
#include "diffusion/conduction.hpp"
#include "srcterms/srcterms.hpp"
#include "bvals/bvals.hpp"
#include "shearing_box/shearing_box.hpp"
#include "shearing_box/orbital_advection.hpp"
#include "mhd/mhd.hpp"
#include "dyn_grmhd/dyn_grmhd.hpp"

//----------------------------------------------------------------------------------------
//! \fn void Resistivity::AssembleResistRKGTasks
//! \brief Adds the tasks of one RKG super-stage into the "*_rkg_*" task lists

void Resistivity::AssembleResistRKGTasks(
    std::map<std::string, std::shared_ptr<TaskList>> tl) {
  using mhd::MHD;
  TaskID none(0);
  MHD *pmhd = pmy_pack->pmhd;
  auto &bti = tl["before_rkg_timeintegrator"];
  auto &bst = tl["before_rkg_stagen"];
  auto &stg = tl["rkg_stagen"];
  auto &ast = tl["after_rkg_stagen"];
  auto &ati = tl["after_rkg_timeintegrator"];

  // assemble "before_timeintegrator" task list
  id.res_totstage = bti->AddTask(&Resistivity::TotStage, this, none);

  // assemble "before_stagen" task list; InitRecv after RKGCoeff, which sets
  // pmhd->rkg_skip_u for this stage
  id.res_coeff = bst->AddTask(&Resistivity::RKGCoeff, this, none);
  id.mhd_irecv = bst->AddTask(&MHD::InitRecv, pmhd, id.res_coeff);

  // assemble "stagen" task list
  id.res_flux      = stg->AddTask(&Resistivity::Fluxes, this, none);
  id.mhd_sendf     = stg->AddTask(&Resistivity::SendFluxRKG, this, id.res_flux);
  id.mhd_recvf     = stg->AddTask(&Resistivity::RecvFluxRKG, this, id.mhd_sendf);
  id.res_copyiniuf = stg->AddTask(&Resistivity::CopyIniConsAndFluxes, this, id.mhd_recvf);
  id.res_rkupdt    = stg->AddTask(&Resistivity::RKUpdate, this, id.res_copyiniuf);
  id.mhd_restu     = stg->AddTask(&MHD::RestrictU, pmhd, id.res_rkupdt);
  id.mhd_sendu     = stg->AddTask(&MHD::SendU, pmhd, id.mhd_restu);
  id.mhd_recvu     = stg->AddTask(&MHD::RecvU, pmhd, id.mhd_sendu);
  id.res_efld      = stg->AddTask(&Resistivity::EField, this, id.mhd_recvu);
  id.mhd_sende     = stg->AddTask(&MHD::SendE, pmhd, id.res_efld);
  id.mhd_recve     = stg->AddTask(&MHD::RecvE, pmhd, id.mhd_sende);
  id.res_copyinie  = stg->AddTask(&Resistivity::CopyIniE, this, id.mhd_recve);
  id.res_ct        = stg->AddTask(&Resistivity::CT, this, id.res_copyinie);
  id.mhd_restb     = stg->AddTask(&MHD::RestrictB, pmhd, id.res_ct);
  id.mhd_sendb     = stg->AddTask(&MHD::SendB, pmhd, id.mhd_restb);
  id.mhd_recvb     = stg->AddTask(&MHD::RecvB, pmhd, id.mhd_sendb);
  id.mhd_bcs       = stg->AddTask(&MHD::ApplyPhysicalBCs, pmhd, id.mhd_recvb);
  id.mhd_prol      = stg->AddTask(&MHD::Prolongate, pmhd, id.mhd_bcs);
  id.res_copyu     = stg->AddTask(&Resistivity::CopyCons, this, id.mhd_prol);
  id.mhd_c2p       = stg->AddTask(&Resistivity::ConToPrimRKG, this, id.res_copyu);
  id.mhd_newdt     = stg->AddTask(&Resistivity::NewTimeStepRKG, this, id.mhd_c2p);

  // assemble "after_stagen" task list
  id.mhd_csend = ast->AddTask(&MHD::ClearSend, pmhd, none);
  // although RecvFlux/U/E/B functions check that all recvs complete, add ClearRecv to
  // task list anyways to catch potential bugs in MPI communication logic
  id.mhd_crecv = ast->AddTask(&MHD::ClearRecv, pmhd, id.mhd_csend);

  // assemble "after_timeintegrator" task list
  id.res_eta        = ati->AddTask(&Resistivity::UpdateResistivity, this, none);
  id.mhd_finalnewdt = ati->AddTask(&MHD::NewTimeStep, pmhd, id.res_eta);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Resistivity::TotStage
//! \brief Number of RKG super-stages s for this super-step (dt x rkg_dtfac over the
//! diffusive dt_diff) and the coefficients of the first stages. A Strang half that
//! starts from the state the RK stages left (rkg_half = 2 or 3) refreshes eta first.

TaskStatus Resistivity::TotStage(Driver *pdrive, int stage) {
  if (rkg_half == 2 || rkg_half == 3) RefreshEta();
  Real dtp = pmy_pack->pmesh->dt_diff;
  Real dth = pmy_pack->pmesh->dt*rkg_dtfac;
  tau = dth/dtp;
  alpha = 0.5;
  Real sreal = -alpha + sqrt(SQR(alpha+1.0)+tau*(3.0+2.0*alpha));
  s = ceil(sreal); // second order
  s = (s < 3) ? 3 : s;
  bjm2 = 1.0/3.0;
  bjm1 = 1.0/3.0;
  w1 = (3.0+2.0*alpha)/((s+2.0*alpha+1.0)*(s-1.0));
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Resistivity::RKGCoeff
//! \brief Recurrence coefficients (b_j, mu, nu, mu~, gamma~) of super-stage `stage`

TaskStatus Resistivity::RKGCoeff(Driver *pdrive, int stage) {
  if (stage > 2) {
    bjm2 = bjm1;
    bjm1 = bj;
  }
  if (stage > 1) {
    bj = (2.0*alpha+1.0)/(3.0+2.0*alpha)*(stage+2.0*alpha+1.0)*(stage-1.0)
         /(stage*(stage+2.0*alpha));
    mu = 2.0*(alpha+stage-1.0)/(stage+2.0*alpha-1.0)*bj/bjm1;
    nu = -(stage-1.0)/(stage+2.0*alpha-1.0)*bj/bjm2;
    mut = mu*w1;
    Real ajm1 = 1.0 - bjm1;
    gat = -mut*ajm1;
  }
  pmy_pack->pmhd->rkg_skip_u = rkg_lean && (stage < s);
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Resistivity::Fluxes
//! \brief Resistive energy fluxes (and the resistive edge EMF) of this super-stage

TaskStatus Resistivity::Fluxes(Driver *pdrive, int stage) {
  AddResistiveFluxes(pmy_pack->pmhd->b0, pmy_pack->pmhd->bcc0, pmy_pack->pmhd->uflx);
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn Resistivity::SendFluxRKG / RecvFluxRKG
//! \brief MHD::SendFlux / RecvFlux for a lean super-stage: the cubed-sphere seam flux
//! reconciliation moves and averages only the energy flux, the one slot RKUpdate reads.
//! The variable filter is set around the call only.

TaskStatus Resistivity::SendFluxRKG(Driver *pdrive, int stage) {
  auto *pm = pmy_pack->pmhd;
  pm->pbval_u->seam_flux_var = rkg_lean ? static_cast<int>(IEN) : -1;
  const TaskStatus st = pm->SendFlux(pdrive, stage);
  pm->pbval_u->seam_flux_var = -1;
  return st;
}

TaskStatus Resistivity::RecvFluxRKG(Driver *pdrive, int stage) {
  auto *pm = pmy_pack->pmhd;
  pm->pbval_u->seam_flux_var = rkg_lean ? static_cast<int>(IEN) : -1;
  const TaskStatus st = pm->RecvFlux(pdrive, stage);
  pm->pbval_u->seam_flux_var = -1;
  return st;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Resistivity::CopyIniConsAndFluxes
//! \brief At super-stage 1, saves the initial state and fluxes (u_ideal, b_ideal,
//! uflx_ideal) and initializes the fourth register (u2, b2); with rkg_lean only the
//! energy slot of the conserved variables and fluxes.

TaskStatus Resistivity::CopyIniConsAndFluxes(Driver *pdrive, int stage) {
  if (stage == 1 && rkg_lean) {
    // only the energy slot of the conserved registers is read (RKUpdate, n = IEN)
    auto &u0 = pmy_pack->pmhd->u0;
    auto &uf = pmy_pack->pmhd->uflx;
    const int n = IEN;
    auto A = Kokkos::ALL;
    Kokkos::deep_copy(DevExeSpace(), Kokkos::subview(u_ideal,A,n,A,A,A),
                      Kokkos::subview(u0,A,n,A,A,A));
    Kokkos::deep_copy(DevExeSpace(), b_ideal.x1f, pmy_pack->pmhd->b0.x1f);
    Kokkos::deep_copy(DevExeSpace(), b_ideal.x2f, pmy_pack->pmhd->b0.x2f);
    Kokkos::deep_copy(DevExeSpace(), b_ideal.x3f, pmy_pack->pmhd->b0.x3f);
    Kokkos::deep_copy(DevExeSpace(), Kokkos::subview(uflx_ideal.x1f,A,n,A,A,A),
                      Kokkos::subview(uf.x1f,A,n,A,A,A));
    Kokkos::deep_copy(DevExeSpace(), Kokkos::subview(uflx_ideal.x2f,A,n,A,A,A),
                      Kokkos::subview(uf.x2f,A,n,A,A,A));
    Kokkos::deep_copy(DevExeSpace(), Kokkos::subview(uflx_ideal.x3f,A,n,A,A,A),
                      Kokkos::subview(uf.x3f,A,n,A,A,A));
    Kokkos::deep_copy(DevExeSpace(), Kokkos::subview(u2,A,n,A,A,A),
                      Kokkos::subview(u0,A,n,A,A,A));
    Kokkos::deep_copy(DevExeSpace(), b2.x1f, pmy_pack->pmhd->b0.x1f);
    Kokkos::deep_copy(DevExeSpace(), b2.x2f, pmy_pack->pmhd->b0.x2f);
    Kokkos::deep_copy(DevExeSpace(), b2.x3f, pmy_pack->pmhd->b0.x3f);
  } else if (stage == 1) {
    Kokkos::deep_copy(DevExeSpace(), u_ideal, pmy_pack->pmhd->u0);
    Kokkos::deep_copy(DevExeSpace(), b_ideal.x1f, pmy_pack->pmhd->b0.x1f);
    Kokkos::deep_copy(DevExeSpace(), b_ideal.x2f, pmy_pack->pmhd->b0.x2f);
    Kokkos::deep_copy(DevExeSpace(), b_ideal.x3f, pmy_pack->pmhd->b0.x3f);
    Kokkos::deep_copy(DevExeSpace(), uflx_ideal.x1f, pmy_pack->pmhd->uflx.x1f);
    Kokkos::deep_copy(DevExeSpace(), uflx_ideal.x2f, pmy_pack->pmhd->uflx.x2f);
    Kokkos::deep_copy(DevExeSpace(), uflx_ideal.x3f, pmy_pack->pmhd->uflx.x3f);
    Kokkos::deep_copy(DevExeSpace(), u2, pmy_pack->pmhd->u0);
    Kokkos::deep_copy(DevExeSpace(), b2.x1f, pmy_pack->pmhd->b0.x1f);
    Kokkos::deep_copy(DevExeSpace(), b2.x2f, pmy_pack->pmhd->b0.x2f);
    Kokkos::deep_copy(DevExeSpace(), b2.x3f, pmy_pack->pmhd->b0.x3f);
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Resistivity::EField
//! \brief Resistive corner EMF of this super-stage

TaskStatus Resistivity::EField(Driver *pdrive, int stage) {
  AddResistiveEMFs(pmy_pack->pmhd->b0, pmy_pack->pmhd->efld);
  // spherical polar: the resistive E_r of the pole x1 edge is formed per phi sector;
  // average it over phi as the explicit path does (<mesh>/use_polar_average_eresist)
  if (pmy_pack->pmesh->use_polar_boundary && pmy_pack->pmesh->use_polar_average_eresist) {
    pmy_pack->pmhd->PolarAzimuthalAverageEr();
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Resistivity::CopyIniE
//! \brief At super-stage 1, saves the initial corner EMF (efld_ideal)

TaskStatus Resistivity::CopyIniE(Driver *pdrive, int stage) {
  if (stage == 1) {
    Kokkos::deep_copy(DevExeSpace(), efld_ideal.x1e, pmy_pack->pmhd->efld.x1e);
    Kokkos::deep_copy(DevExeSpace(), efld_ideal.x2e, pmy_pack->pmhd->efld.x2e);
    Kokkos::deep_copy(DevExeSpace(), efld_ideal.x3e, pmy_pack->pmhd->efld.x3e);
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Resistivity::CopyCons
//! \brief Shifts the registers of the recurrence (u1 -> u2, u0 -> u1). b2 is rotated in
//! place by CT and, with rkg_lean, the energy slot of u2 by RKUpdate; b1 (and with
//! rkg_lean u1) is read by no RKG task and is refreshed at the last super-stage only.

TaskStatus Resistivity::CopyCons(Driver *pdrive, int stage) {
  auto *pm = pmy_pack->pmhd;
  if (rkg_lean) {
    // u2 is read only in its energy slot; u1 is fully refreshed at stage 1 (it held the
    // RK stage state) and afterwards differs from u0 only in the energy slot
    const int n = IEN;
    auto A = Kokkos::ALL;
    if (stage == 1) {
      Kokkos::deep_copy(DevExeSpace(), pm->u1, pm->u0);
    } else if (stage == s) {
      Kokkos::deep_copy(DevExeSpace(), Kokkos::subview(pm->u1,A,n,A,A,A),
                        Kokkos::subview(pm->u0,A,n,A,A,A));
    }
  } else {
    // the full u1 -> u2 shift stays: the ConToPrim of an intermediate super-stage
    // (rkg_c2p_shift) reads u1 and u2 in every cell
    if (stage > 1) {
      Kokkos::deep_copy(DevExeSpace(), u2, pm->u1);
    }
    Kokkos::deep_copy(DevExeSpace(), pm->u1, pm->u0);
  }
  if (stage == s) {
    Kokkos::deep_copy(DevExeSpace(), pm->b1.x1f, pm->b0.x1f);
    Kokkos::deep_copy(DevExeSpace(), pm->b1.x2f, pm->b0.x2f);
    Kokkos::deep_copy(DevExeSpace(), pm->b1.x3f, pm->b0.x3f);
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn Resistivity::ConToPrimRKG / NewTimeStepRKG
//! \brief Without rkg_lean: MHD::ConToPrim and MHD::NewTimeStep (the latter fires at
//! super-stage nexp_stages, since it gates on that number). With rkg_lean: the
//! intermediate super-stages rebuild only the cell-centred field; the last one runs the
//! full inversion and then the timestep, on the final state.

TaskStatus Resistivity::ConToPrimRKG(Driver *pdrive, int stage) {
  if (rkg_lean && stage < s) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    int &ng = indcs.ng;
    int n1m1 = indcs.nx1 + 2*ng - 1;
    int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
    int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
    pmy_pack->pcoord->GnomonicCellCenteredB(pmy_pack->pmhd->b0, pmy_pack->pmhd->bcc0,
                                            0, n1m1, 0, n2m1, 0, n3m1);
    return TaskStatus::complete;
  }
  // the ConToPrim of an intermediate super-stage (floors, fallbacks) changes u0(IEN)
  // only, after CopyCons saved u_j to u1. The change d is added to every register the
  // recurrence reads (u1 = u_j, u2 = u_{j-1}, u_ideal = u_0); since mu + nu + (1-mu-nu)
  // = 1 it survives in the final state with weight 1, as if the floor had been applied
  // at the start of the super-step. rkg_lean converts only at the last stage.
  auto *pmhd = pmy_pack->pmhd;
  pmhd->c2p_wt_rkg = 1.0;
  const TaskStatus st = pmhd->ConToPrim(pdrive, stage);
  pmhd->c2p_wt_rkg = -1.0;
  if (!rkg_lean && stage < s) {
    auto u0_ = pmhd->u0;
    auto u1_ = pmhd->u1;
    auto u2_ = u2;
    auto ui_ = u_ideal;
    const int n = IEN;
    const int nmb1 = pmy_pack->nmb_thispack - 1;
    const int n3 = static_cast<int>(u0_.extent_int(2)) - 1;
    const int n2 = static_cast<int>(u0_.extent_int(3)) - 1;
    const int n1 = static_cast<int>(u0_.extent_int(4)) - 1;
    par_for("rkg_c2p_shift", DevExeSpace(), 0, nmb1, 0, n3, 0, n2, 0, n1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real d = u0_(m,n,k,j,i) - u1_(m,n,k,j,i);
      if (d != 0.0) {
        u1_(m,n,k,j,i) = u0_(m,n,k,j,i);
        u2_(m,n,k,j,i) += d;
        ui_(m,n,k,j,i) += d;
      }
    });
  }
  return st;
}

TaskStatus Resistivity::NewTimeStepRKG(Driver *pdrive, int stage) {
  if (!rkg_lean) return pmy_pack->pmhd->NewTimeStep(pdrive, stage);
  if (stage < s) return TaskStatus::complete;
  return pmy_pack->pmhd->NewTimeStep(pdrive, pdrive->nexp_stages);
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Resistivity::UpdateResistivity
//! \brief After the super-step: combines the x1-face energy fluxes of the Strang halves
//! and refreshes eta from the new state

TaskStatus Resistivity::UpdateResistivity(Driver *pdrive, int stage) {
  pmy_pack->pmhd->rkg_skip_u = false;
  // the x1-face energy flux of the cycle is the mean of the two halves' fluxes, each
  // applied over dt/2
  if (track_x1_bflx && rkg_half != 0 && rkg_bflx.extent_int(1) >= 8) {
    auto bf = rkg_bflx;
    const int h = rkg_half;
    const int nmb1 = rkg_bflx.extent_int(0) - 1;
    const int n3 = rkg_bflx.extent_int(2) - 1, n2 = rkg_bflx.extent_int(3) - 1;
    par_for("rkg_bflx_half", DevExeSpace(), 0, nmb1, 0, n3, 0, n2,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      for (int s2=0; s2<2; ++s2) {
        if (h == 1 || h == 3) {   // 3 = the leading super-step of a fused cycle
          bf(m,6+s2,k,j) = bf(m,4+s2,k,j);
        } else {
          bf(m,4+s2,k,j) = 0.5*(bf(m,6+s2,k,j) + bf(m,4+s2,k,j));
        }
      }
    });
  }
  RefreshEta();
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void Resistivity::RefreshEta
//! \brief eta_b (perna/eos) and ad_coef (eos) from w0 in every cell, ghosts included

void Resistivity::RefreshEta() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int &ng = indcs.ng;
  int n1m1 = indcs.nx1 + 2*ng - 1;
  int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
  int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
  if (OhmicNeedsUpdate()) {
    SetResistivity(pmy_pack->pmhd->w0, pmy_pack->pmhd->peos->eos_data,
                   pmy_pack->pmesh->pgen->hot_jupiter_param.Rgas, eta_b,
                   0, n1m1, 0, n2m1, 0, n3m1);
  }
  if (AmbipolarNeedsUpdate()) {
    SetAmbipolar(pmy_pack->pmhd->w0, pmy_pack->pmhd->peos->eos_data,
                 pmy_pack->pmesh->pgen->hot_jupiter_param.Rgas,
                 0, n1m1, 0, n2m1, 0, n3m1);
  }
}
