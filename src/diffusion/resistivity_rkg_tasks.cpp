//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file resistivity_rkg_tasks.cpp
//! \brief functions that control resistivity rkg tasks stored in tasklists in MeshBlockPack

#include <map>
#include <memory>
#include <string>
#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"
#include "mesh/mesh.hpp"
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

void Resistivity::AssembleResistRKGTasks(std::map<std::string, std::shared_ptr<TaskList>> tl) {
  TaskID none(0);
  mhd::MHD *pmhd = pmy_pack->pmhd;

  // assemble "before_timeintegrator" task list
  id.res_totstage = tl["before_rkg_timeintegrator"]->AddTask(&Resistivity::TotStage, this, none);

  // assemble "before_stagen" task list
  id.res_coeff = tl["before_rkg_stagen"]->AddTask(&Resistivity::RKGCoeff, this, none);
  id.mhd_irecv = tl["before_rkg_stagen"]->AddTask(&mhd::MHD::InitRecv, pmhd, none);

  // assemble "stagen" task list
  id.res_flux      = tl["rkg_stagen"]->AddTask(&Resistivity::Fluxes, this, none);
  id.mhd_sendf     = tl["rkg_stagen"]->AddTask(&mhd::MHD::SendFlux, pmhd, id.res_flux);
  id.mhd_recvf     = tl["rkg_stagen"]->AddTask(&mhd::MHD::RecvFlux, pmhd, id.mhd_sendf);
  id.res_copyiniuf = tl["rkg_stagen"]->AddTask(&Resistivity::CopyIniConsAndFluxes, this, id.mhd_recvf);
  id.res_rkupdt    = tl["rkg_stagen"]->AddTask(&Resistivity::RKUpdate, this, id.res_copyiniuf);
  id.mhd_restu     = tl["rkg_stagen"]->AddTask(&mhd::MHD::RestrictU, pmhd, id.res_rkupdt);
  id.mhd_sendu     = tl["rkg_stagen"]->AddTask(&mhd::MHD::SendU, pmhd, id.mhd_restu);
  id.mhd_recvu     = tl["rkg_stagen"]->AddTask(&mhd::MHD::RecvU, pmhd, id.mhd_sendu);
  id.res_efld      = tl["rkg_stagen"]->AddTask(&Resistivity::EField, this, id.mhd_recvu);
  id.mhd_sende     = tl["rkg_stagen"]->AddTask(&mhd::MHD::SendE, pmhd, id.res_efld);
  id.mhd_recve     = tl["rkg_stagen"]->AddTask(&mhd::MHD::RecvE, pmhd, id.mhd_sende);
  id.res_copyinie = tl["rkg_stagen"]->AddTask(&Resistivity::CopyIniE, this, id.mhd_recve);
  id.res_ct        = tl["rkg_stagen"]->AddTask(&Resistivity::CT, this, id.res_copyinie);
  id.mhd_restb     = tl["rkg_stagen"]->AddTask(&mhd::MHD::RestrictB, pmhd, id.res_ct);
  id.mhd_sendb     = tl["rkg_stagen"]->AddTask(&mhd::MHD::SendB, pmhd, id.mhd_restb);
  id.mhd_recvb     = tl["rkg_stagen"]->AddTask(&mhd::MHD::RecvB, pmhd, id.mhd_sendb);
  id.mhd_bcs       = tl["rkg_stagen"]->AddTask(&mhd::MHD::ApplyPhysicalBCs, pmhd, id.mhd_recvb);
  id.mhd_prol      = tl["rkg_stagen"]->AddTask(&mhd::MHD::Prolongate, pmhd, id.mhd_bcs);
  id.res_copyu     = tl["rkg_stagen"]->AddTask(&Resistivity::CopyCons, this, id.mhd_prol);
  id.mhd_c2p       = tl["rkg_stagen"]->AddTask(&mhd::MHD::ConToPrim, pmhd, id.res_copyu);
  id.mhd_newdt     = tl["rkg_stagen"]->AddTask(&mhd::MHD::NewTimeStep, pmhd, id.mhd_c2p);

  // assemble "after_stagen" task list
  id.mhd_csend = tl["after_rkg_stagen"]->AddTask(&mhd::MHD::ClearSend, pmhd, none);
  // although RecvFlux/U/E/B functions check that all recvs complete, add ClearRecv to
  // task list anyways to catch potential bugs in MPI communication logic
  id.mhd_crecv = tl["after_rkg_stagen"]->AddTask(&mhd::MHD::ClearRecv, pmhd, id.mhd_csend);
    
  id.res_eta = tl["after_rkg_timeintegrator"]->AddTask(&Resistivity::UpdateResistivity, this, none);
  id.mhd_finalnewdt     = tl["after_rkg_timeintegrator"]->AddTask(&mhd::MHD::NewTimeStep, pmhd, id.res_eta);

  return;
}

TaskStatus Resistivity::TotStage(Driver *pdrive, int stage) {
  Real dtp = pmy_pack->pmesh->dt_diff;
  Real dth = pmy_pack->pmesh->dt;
  tau = dth/dtp;
  alpha = 0.5;
  Real sreal = -alpha + sqrt(SQR(alpha+1.0)+tau*(3.0+2.0*alpha));
  s = ceil(sreal); // second order
  s = (s < 3) ? 3 : s;
  bjm2 = 1.0/3.0;
  bjm1 = 1.0/3.0;
  w1 = (3.0+2.0*alpha)/((s+2.0*alpha+1.0)*(s-1.0));
//  tau = 1.0/w1;
//    std::cout << "RKG: dt=" << dth
//                << " dt_diff=" << dtp
//                << " taureal=" << (dth/dtp)
//                << " sreal=" << sreal
//                << " s=" << s
//                << std::endl;
  return TaskStatus::complete;
}

TaskStatus Resistivity::RKGCoeff(Driver *pdrive, int stage) {
  if (stage > 2) {
    bjm2 = bjm1;
    bjm1 = bj;
  }
  if (stage > 1) {
    bj = (2.0*alpha+1.0)/(3.0+2.0*alpha)*(stage+2.0*alpha+1.0)*(stage-1.0)/(stage*(stage+2.0*alpha));
    mu = 2.0*(alpha+stage-1.0)/(stage+2.0*alpha-1.0)*bj/bjm1;
    nu = -(stage-1.0)/(stage+2.0*alpha-1.0)*bj/bjm2;
    mut = mu*w1;
    Real ajm1 = 1.0 - bjm1;
    gat = -mut*ajm1;
  }
  return TaskStatus::complete;
}

TaskStatus Resistivity::Fluxes(Driver *pdrive, int stage) {
  AddResistiveFluxes(pmy_pack->pmhd->b0, pmy_pack->pmhd->bcc0, pmy_pack->pmhd->uflx);
  return TaskStatus::complete;
}

TaskStatus Resistivity::CopyIniConsAndFluxes(Driver *pdrive, int stage) {
  if (stage == 1) {
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

TaskStatus Resistivity::EField(Driver *pdrive, int stage) {
  // Use CT to compute corner E
  AddResistiveEMFs(pmy_pack->pmhd->b0, pmy_pack->pmhd->efld);
  return TaskStatus::complete;
}

TaskStatus Resistivity::CopyIniE(Driver *pdrive, int stage) {
  if (stage == 1) {
    Kokkos::deep_copy(DevExeSpace(), efld_ideal.x1e, pmy_pack->pmhd->efld.x1e);
    Kokkos::deep_copy(DevExeSpace(), efld_ideal.x2e, pmy_pack->pmhd->efld.x2e);
    Kokkos::deep_copy(DevExeSpace(), efld_ideal.x3e, pmy_pack->pmhd->efld.x3e);
  }
  return TaskStatus::complete;
}

TaskStatus Resistivity::CopyCons(Driver *pdrive, int stage) {
  if (stage > 1) {
    Kokkos::deep_copy(DevExeSpace(), u2, pmy_pack->pmhd->u1);
    Kokkos::deep_copy(DevExeSpace(), b2.x1f, pmy_pack->pmhd->b1.x1f);
    Kokkos::deep_copy(DevExeSpace(), b2.x2f, pmy_pack->pmhd->b1.x2f);
    Kokkos::deep_copy(DevExeSpace(), b2.x3f, pmy_pack->pmhd->b1.x3f);
  }
  Kokkos::deep_copy(DevExeSpace(), pmy_pack->pmhd->u1, pmy_pack->pmhd->u0);
  Kokkos::deep_copy(DevExeSpace(), pmy_pack->pmhd->b1.x1f, pmy_pack->pmhd->b0.x1f);
  Kokkos::deep_copy(DevExeSpace(), pmy_pack->pmhd->b1.x2f, pmy_pack->pmhd->b0.x2f);
  Kokkos::deep_copy(DevExeSpace(), pmy_pack->pmhd->b1.x3f, pmy_pack->pmhd->b0.x3f);
  return TaskStatus::complete;
}

TaskStatus Resistivity::UpdateResistivity(Driver *pdrive, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int &ng = indcs.ng;
  int n1m1 = indcs.nx1 + 2*ng - 1;
  int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
  int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
  if (iso_resist_type.compare("constant") != 0) {
    SetResistivity(pmy_pack->pmhd->w0, pmy_pack->pmhd->peos->eos_data, pmy_pack->pmesh->pgen->hot_jupiter_param.Rgas, eta_b, 0, n1m1, 0, n2m1, 0, n3m1);
  }
  return TaskStatus::complete;
}
