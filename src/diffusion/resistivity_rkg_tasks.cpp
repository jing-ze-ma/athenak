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

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

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
#include "mesh/nghbr_index.hpp"

void Resistivity::AssembleResistRKGTasks(std::map<std::string, std::shared_ptr<TaskList>> tl) {
  TaskID none(0);
  mhd::MHD *pmhd = pmy_pack->pmhd;

  // assemble "before_timeintegrator" task list
  id.res_totstage = tl["before_rkg_timeintegrator"]->AddTask(&Resistivity::TotStage, this, none);

  // assemble "before_stagen" task list
  id.res_coeff = tl["before_rkg_stagen"]->AddTask(&Resistivity::RKGCoeff, this, none);
  // after RKGCoeff: that is where rkg_lean sets pmhd->rkg_skip_u for this stage
  id.mhd_irecv = tl["before_rkg_stagen"]->AddTask(&mhd::MHD::InitRecv, pmhd, id.res_coeff);

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
  id.mhd_c2p       = tl["rkg_stagen"]->AddTask(&Resistivity::ConToPrimRKG, this, id.res_copyu);
  id.mhd_newdt     = tl["rkg_stagen"]->AddTask(&Resistivity::NewTimeStepRKG, this, id.mhd_c2p);

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
  // rkgdt-1101 <mhd>/rkg_lmax_diag: once, on the state the first super-step starts from
  if (rkg_lmax_diag > 0 && !rkg_lmax_done) {
    rkg_lmax_done = true;
    LmaxDiag();
  }
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
  pmy_pack->pmhd->rkg_skip_u = rkg_lean && (stage < s);
  if (rkg_halo_ien_only) {
    // rkghalo-1101: super-stage 1 exchanges everything (the state left by the RK stages
    // and the operator-split work before the super-step); stages >= 2 move u0(IEN) and
    // the slots the previous super-stage's ConToPrim wrote on any rank
    auto *pm = pmy_pack->pmhd;
    const int nvar = pm->nmhd + pm->nscalars;
    int vm = 0;
    if (stage > 1 && nvar <= 32 && rkg_halo_changed >= 0) {
      const int all = (nvar == 32) ? -1 : ((1 << nvar) - 1);
      const int mom = (1 << IM1) | (1 << IM2) | (1 << IM3);
      vm = rkg_halo_changed | (1 << IEN);
      if (vm & mom) {vm |= mom;}
      if (vm == all) {vm = 0;}
    }
    if (stage > 1) {
      rkg_halo_nstage++;
      if (vm == 0) {rkg_halo_nfull++;}
      for (int v=0; v<nvar && v<32; ++v) {
        if (vm == 0 || ((vm >> v) & 1)) {rkg_halo_nslot[v]++;}
      }
    }
    pm->rkg_u_vmask = vm;
  }
  return TaskStatus::complete;
}

TaskStatus Resistivity::Fluxes(Driver *pdrive, int stage) {
  AddResistiveFluxes(pmy_pack->pmhd->b0, pmy_pack->pmhd->bcc0, pmy_pack->pmhd->uflx);
  return TaskStatus::complete;
}

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
  if (rkg_lean) {
    // u2 is read only in its energy slot; u1 is fully refreshed at stage 1 (it held the
    // RK stage state) and afterwards differs from u0 only in the energy slot
    auto &pm = pmy_pack->pmhd;
    const int n = IEN;
    auto A = Kokkos::ALL;
    if (stage > 1) {
      Kokkos::deep_copy(DevExeSpace(), Kokkos::subview(u2,A,n,A,A,A),
                        Kokkos::subview(pm->u1,A,n,A,A,A));
      Kokkos::deep_copy(DevExeSpace(), b2.x1f, pm->b1.x1f);
      Kokkos::deep_copy(DevExeSpace(), b2.x2f, pm->b1.x2f);
      Kokkos::deep_copy(DevExeSpace(), b2.x3f, pm->b1.x3f);
      Kokkos::deep_copy(DevExeSpace(), Kokkos::subview(pm->u1,A,n,A,A,A),
                        Kokkos::subview(pm->u0,A,n,A,A,A));
    } else {
      Kokkos::deep_copy(DevExeSpace(), pm->u1, pm->u0);
    }
    Kokkos::deep_copy(DevExeSpace(), pm->b1.x1f, pm->b0.x1f);
    Kokkos::deep_copy(DevExeSpace(), pm->b1.x2f, pm->b0.x2f);
    Kokkos::deep_copy(DevExeSpace(), pm->b1.x3f, pm->b0.x3f);
    return TaskStatus::complete;
  }
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

//----------------------------------------------------------------------------------------
//! \fn Resistivity::ConToPrimRKG / NewTimeStepRKG
//! \brief Without rkg_lean: MHD::ConToPrim and MHD::NewTimeStep, as before (the latter
//! then fires at super-stage nexp_stages, i.e. 2, since it gates on that number).
//! With rkg_lean: the intermediate super-stages rebuild only the cell-centred field; the
//! last one runs the full inversion and then the timestep, on the final state.

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
  // stswall-1029: the ConToPrim of an intermediate super-stage (floors, fallbacks)
  // changes u0(IEN) only, after CopyCons saved u_j to u1.  The recurrence
  // u_{j+1} = mu u_j + nu u_{j-1} + (1-mu-nu) u_ideal - ... then carried that change with
  // the weight X_s of a perturbation of the u0 slot alone (14..22 for 34 super-stages:
  // 5.8e-6 Lsw of floor energy became 8.9e-5 Lsw on the WASP-121b s10 state,
  // stswall_1029/RESULTS.md).  The change d is now added to every register the
  // recurrence reads (u1 = u_j, u2 = u_{j-1}, u_ideal = u_0): mu + nu + (1-mu-nu) = 1,
  // so it survives in the final state with weight exactly 1, as if the floor had been
  // applied at the start of the super-step, and it is booked with weight 1.  Where no
  // floor fires d = 0 and nothing changes.  rkg_lean converts only at the last stage.
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
    if (rkg_halo_ien_only) {
      // rkghalo-1101: did this ConToPrim write any slot but IEN anywhere (ghosts too)?
      // u1 = u0 before it (CopyCons); RKUpdate between two super-stages writes only
      // IEN, so if nothing else changed the other slots' ghosts from the last full
      // exchange are still exactly what a full exchange would deliver.
      // ghosts across a PHYSICAL face (no neighbour) are never written by the exchange
      // (sends read active cells; the cube-vertex widening is tangential), so what the
      // ConToPrim does there is the same with or without the key: skip them
      if (rkg_halo_phys.extent_int(0) != nmb1 + 1) {
        Kokkos::realloc(rkg_halo_phys, nmb1 + 1);
        auto ph = Kokkos::create_mirror_view(rkg_halo_phys);
        auto &ngh = pmy_pack->pmb->nghbr;
        auto *pmsh = pmy_pack->pmesh;
        for (int m=0; m<=nmb1; ++m) {
          int b = 0;
          for (int f=0; f<6; ++f) {
            const int ax = f/2, sg = (f%2 == 0) ? -1 : 1;
            if ((ax == 1 && !pmsh->multi_d) || (ax == 2 && !pmsh->three_d)) continue;
            const int nb = NeighborIndex((ax == 0) ? sg : 0, (ax == 1) ? sg : 0,
                                         (ax == 2) ? sg : 0, 0, 0);
            if (ngh.h_view(m,nb).gid < 0) {b |= (1 << f);}
          }
          ph(m) = b;
        }
        Kokkos::deep_copy(rkg_halo_phys, ph);
      }
      auto phys_ = rkg_halo_phys;
      auto &mi_ = pmy_pack->pmesh->mb_indcs;
      const int is_ = mi_.is, ie_ = mi_.ie;
      const int js_ = mi_.js, je_ = mi_.je, ks_ = mi_.ks, ke_ = mi_.ke;
      const bool md_ = pmy_pack->pmesh->multi_d, td_ = pmy_pack->pmesh->three_d;
      const int nv = u0_.extent_int(1);
      const int64_t nc1 = n1 + 1, nc2 = n2 + 1, nc3 = n3 + 1;
      const int64_t ncell = nc1*nc2*nc3;
      const int64_t ntot = static_cast<int64_t>(nmb1 + 1)*nv*ncell;
      int ndiff = 0;   // bit v: slot v changed somewhere
      Kokkos::parallel_reduce("rkg_halo_check",
        Kokkos::RangePolicy<DevExeSpace, Kokkos::IndexType<int64_t>>(DevExeSpace(), 0,
                                                                     ntot),
      KOKKOS_LAMBDA(const int64_t idx, int &nd) {
        const int64_t mv = idx/ncell;
        int64_t r = idx - mv*ncell;
        const int m = static_cast<int>(mv/nv);
        const int v = static_cast<int>(mv - static_cast<int64_t>(m)*nv);
        const int k = static_cast<int>(r/(nc2*nc1));
        r -= static_cast<int64_t>(k)*nc2*nc1;
        const int j = static_cast<int>(r/nc1);
        const int i = static_cast<int>(r - static_cast<int64_t>(j)*nc1);
        const int b = phys_(m);
        const bool skip = ((b & 1) && i < is_) || ((b & 2) && i > ie_) ||
                          (md_ && (((b & 4) && j < js_) || ((b & 8) && j > je_))) ||
                          (td_ && (((b & 16) && k < ks_) || ((b & 32) && k > ke_)));
        if (v != n && v < 32 && !skip && u0_(m,v,k,j,i) != u1_(m,v,k,j,i)) {
          nd |= (1 << v);
        }
      }, Kokkos::BOr<int>(ndiff));
#if MPI_PARALLEL_ENABLED
      MPI_Allreduce(MPI_IN_PLACE, &ndiff, 1, MPI_INT, MPI_BOR, MPI_COMM_WORLD);
#endif
      rkg_halo_changed = ndiff;
    }
  }
  return st;
}

TaskStatus Resistivity::NewTimeStepRKG(Driver *pdrive, int stage) {
  if (!rkg_lean) return pmy_pack->pmhd->NewTimeStep(pdrive, stage);
  if (stage < s) return TaskStatus::complete;
  return pmy_pack->pmhd->NewTimeStep(pdrive, pdrive->nexp_stages);
}

TaskStatus Resistivity::UpdateResistivity(Driver *pdrive, int stage) {
  pmy_pack->pmhd->rkg_skip_u = false;
  pmy_pack->pmhd->rkg_u_vmask = 0;   // rkghalo-1101: the RK stages exchange everything
  rkg_halo_changed = -1;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int &ng = indcs.ng;
  int n1m1 = indcs.nx1 + 2*ng - 1;
  int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
  int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
  if (OhmicNeedsUpdate()) {
    SetResistivity(pmy_pack->pmhd->w0, pmy_pack->pmhd->peos->eos_data, pmy_pack->pmesh->pgen->hot_jupiter_param.Rgas, eta_b, 0, n1m1, 0, n2m1, 0, n3m1);
  }
  if (AmbipolarNeedsUpdate()) {
    SetAmbipolar(pmy_pack->pmhd->w0, pmy_pack->pmhd->peos->eos_data,
                 pmy_pack->pmesh->pgen->hot_jupiter_param.Rgas,
                 0, n1m1, 0, n2m1, 0, n3m1);
  }
  return TaskStatus::complete;
}
