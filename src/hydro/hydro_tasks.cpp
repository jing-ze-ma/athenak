//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hydro_tasks.cpp
//! \brief functions that control Hydro tasks stored in tasklists in MeshBlockPack

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "eos/eos.hpp"
#include "diffusion/viscosity.hpp"
#include "diffusion/conduction.hpp"
#include "utils/runaway_scan.hpp"
#include "srcterms/srcterms.hpp"
#include "bvals/bvals.hpp"
#include "shearing_box/shearing_box.hpp"
#include "shearing_box/orbital_advection.hpp"
#include "hydro/hydro.hpp"
#include "hydro/hydro_am.hpp"
#include "rad_m1/rad_m1.hpp"
#include "utils/c2p_track.hpp"

namespace hydro {

//----------------------------------------------------------------------------------------
// <problem>/nan_report diagnostics.  Both helpers are no-ops unless Hydro::nan_report is
// set, so a default run executes not one extra kernel and is bit-identical.  They exist
// to say WHICH stage-level operator first writes a non-finite value; see the call sites
// in Hydro::Fluxes and Hydro::HydroSrcTerms.
namespace {
DvceArray1D<int> *nanrep_cnt = nullptr;
DvceArray1D<Real> *nanrep_rec = nullptr;
int nanrep_lines = 0;
const int nanrep_maxlines = 400;

void NanRepAlloc() {
  if (nanrep_cnt == nullptr) {
    nanrep_cnt = new DvceArray1D<int>("nanrep_cnt", 1);
    nanrep_rec = new DvceArray1D<Real>("nanrep_rec", 16);
  }
  Kokkos::deep_copy(*nanrep_cnt, 0);
  Kokkos::deep_copy(*nanrep_rec, 0.0);
}

//! \brief scan one face-centered flux component for a non-finite IDN/IM1/IEN and report
//! the first offender with the two cell states that made it.  dir = 1, 2 or 3.
void NanScanFlux(MeshBlockPack *pmbp, const DvceArray5D<Real> &flx, const int dir,
                 const DvceArray5D<Real> &w0, const EOS_Data &eos, const char *tag) {
  if (nanrep_lines >= nanrep_maxlines) return;
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  const int di = (dir == 1), dj = (dir == 2), dk = (dir == 3);
  NanRepAlloc();
  auto ncnt = *nanrep_cnt;
  auto nrec = *nanrep_rec;
  const bool gen = eos.IsGeneral();
  const Real tunit = eos.temp_cgs, gm1 = eos.gamma - 1.0;
  par_for("nanscan_flx", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real fd = flx(m,IDN,k,j,i);
    const Real f1 = flx(m,IM1,k,j,i);
    const Real fe = flx(m,IEN,k,j,i);
    if (isfinite(fd) && isfinite(f1) && isfinite(fe)) return;
    if (Kokkos::atomic_fetch_add(&ncnt(0), 1) != 0) return;
    const int kl = k-dk, jl = j-dj, il = i-di;
    const Real dl = w0(m,IDN,kl,jl,il), el = w0(m,IEN,kl,jl,il);
    const Real dr = w0(m,IDN,k,j,i), er = w0(m,IEN,k,j,i);
    nrec(0) = static_cast<Real>(m);
    nrec(1) = static_cast<Real>(k);
    nrec(2) = static_cast<Real>(j);
    nrec(3) = static_cast<Real>(i);
    nrec(4) = fd;
    nrec(5) = f1;
    nrec(6) = fe;
    nrec(7) = dl;
    nrec(8) = el;
    nrec(9) = gen ? eos.Temperature(dl, el)*tunit : el/dl*gm1*tunit;
    nrec(10) = w0(m,IVX,kl,jl,il);
    nrec(11) = dr;
    nrec(12) = er;
    nrec(13) = gen ? eos.Temperature(dr, er)*tunit : er/dr*gm1*tunit;
    nrec(14) = w0(m,IVX,k,j,i);
  });
  auto hc = Kokkos::create_mirror_view(ncnt);
  Kokkos::deep_copy(hc, ncnt);
  if (hc(0) <= 0) return;
  auto hr = Kokkos::create_mirror_view(nrec);
  Kokkos::deep_copy(hr, nrec);
  ++nanrep_lines;
  const int mb = static_cast<int>(hr(0));
  std::cout << "### nan_report [" << tag << "] rank " << global_variable::my_rank
            << " cycle " << pmbp->pmesh->ncycle << " t = " << pmbp->pmesh->time
            << ": " << hc(0) << " bad face(s); first (m,k,j,i) = (" << mb << ","
            << static_cast<int>(hr(1)) << "," << static_cast<int>(hr(2)) << ","
            << static_cast<int>(hr(3)) << ") gid = " << (pmbp->gids + mb)
            << " flx(IDN) = " << hr(4) << " flx(IM1) = " << hr(5)
            << " flx(IEN) = " << hr(6)
            << " | L: d = " << hr(7) << " e = " << hr(8) << " T = " << hr(9)
            << " v1 = " << hr(10)
            << " | R: d = " << hr(11) << " e = " << hr(12) << " T = " << hr(13)
            << " v1 = " << hr(14) << std::endl;
}

//! \brief scan the conserved variables for a non-finite or non-positive state.
void NanScanCons(MeshBlockPack *pmbp, const DvceArray5D<Real> &u0, const char *tag) {
  if (nanrep_lines >= nanrep_maxlines) return;
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  NanRepAlloc();
  auto ncnt = *nanrep_cnt;
  auto nrec = *nanrep_rec;
  par_for("nanscan_u", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real d = u0(m,IDN,k,j,i), e = u0(m,IEN,k,j,i);
    const Real m1 = u0(m,IM1,k,j,i), m2 = u0(m,IM2,k,j,i), m3 = u0(m,IM3,k,j,i);
    const bool bad = !isfinite(d) || !isfinite(m1) || !isfinite(m2) || !isfinite(m3)
                     || !isfinite(e) || !(e > 0.0) || !(d > 0.0);
    if (!bad) return;
    if (Kokkos::atomic_fetch_add(&ncnt(0), 1) != 0) return;
    nrec(0) = static_cast<Real>(m);
    nrec(1) = static_cast<Real>(k);
    nrec(2) = static_cast<Real>(j);
    nrec(3) = static_cast<Real>(i);
    nrec(4) = d;
    nrec(5) = e;
    nrec(6) = m1;
    nrec(7) = m2;
    nrec(8) = m3;
  });
  auto hc = Kokkos::create_mirror_view(ncnt);
  Kokkos::deep_copy(hc, ncnt);
  if (hc(0) <= 0) return;
  auto hr = Kokkos::create_mirror_view(nrec);
  Kokkos::deep_copy(hr, nrec);
  ++nanrep_lines;
  const int mb = static_cast<int>(hr(0));
  std::cout << "### nan_report [" << tag << "] rank " << global_variable::my_rank
            << " cycle " << pmbp->pmesh->ncycle << " t = " << pmbp->pmesh->time
            << ": " << hc(0) << " bad cell(s); first (m,k,j,i) = (" << mb << ","
            << static_cast<int>(hr(1)) << "," << static_cast<int>(hr(2)) << ","
            << static_cast<int>(hr(3)) << ") gid = " << (pmbp->gids + mb)
            << " u(IDN) = " << hr(4) << " u(IEN) = " << hr(5)
            << " u(IM1) = " << hr(6) << " u(IM2) = " << hr(7)
            << " u(IM3) = " << hr(8) << std::endl;
}
}  // namespace
//----------------------------------------------------------------------------------------
//! \fn  void Hydro::AssembleHydroTasks
//! \brief Adds hydro tasks to appropriate task lists used by time integrators.
//! Called by MeshBlockPack::AddPhysics() function directly after Hydro constructor.
//! Many of the functions in the task list are implemented in this file because they are
//! simple, or they are wrappers that call one or more other functions.
//!
//! "before_stagen" tasks are those that must be cmpleted over all MeshBlocks BEFORE each
//! stage can be run (such as posting MPI receives, setting BoundaryCommStatus flags, etc)
//!
//! "stagen" tasks are those performed DURING each stage
//!
//! "after_stagen" tasks are those that can only be completed AFTER all the "stagen" tasks
//! are completed over ALL MeshBlocks for each stage, such as clearing all MPI calls, etc.
//!
//! In addition there are "before_timeintegrator" and "after_timeintegrator" task lists
//! in the tl map, which are generally used for operator split tasks.

void Hydro::AssembleHydroTasks(std::map<std::string, std::shared_ptr<TaskList>> tl) {
  TaskID none(0);

  // assemble "before_stagen" task list
  // ---- the Strang-split operator-split source (ProblemGenerator::user_split_func)
  // Half the CYCLE dt before the integrator and half after, each followed by a
  // ConToPrim so the hydro stage (and the next cycle) sees consistent primitives.
  // The hook is null unless a problem generator enrols it, and both tasks then return
  // immediately, so this is bitwise inert for every other run.
  id.splitpre  = tl["before_timeintegrator"]->AddTask(&Hydro::RTStrangSplit, this,
                                                      none);
  id.splitpc2p = tl["before_timeintegrator"]->AddTask(&Hydro::ConToPrimSplit, this,
                                                      id.splitpre);
  id.splitpst  = tl["after_timeintegrator"]->AddTask(&Hydro::RTStrangSplit, this, none);
  id.splitsc2p = tl["after_timeintegrator"]->AddTask(&Hydro::ConToPrimSplit, this,
                                                     id.splitpst);
  pmy_pack->split_hook_tasks = true;

  id.irecv = tl["before_stagen"]->AddTask(&Hydro::InitRecv, this, none);

  // assemble "stagen" task list
  id.copyu     = tl["stagen"]->AddTask(&Hydro::CopyCons, this, none);
  id.flux      = tl["stagen"]->AddTask(&Hydro::Fluxes,this,id.copyu);
  id.sendf     = tl["stagen"]->AddTask(&Hydro::SendFlux, this, id.flux);
  id.recvf     = tl["stagen"]->AddTask(&Hydro::RecvFlux, this, id.sendf);
  id.rkupdt    = tl["stagen"]->AddTask(&Hydro::RKUpdate, this, id.recvf);
  id.srctrms   = tl["stagen"]->AddTask(&Hydro::HydroSrcTerms, this, id.rkupdt);
  // the implicit radial radiative diffusion (<hydro>/rad_implicit_x1) sits between the
  // explicit update and the ghost exchange, so what it writes is what is communicated
  id.impcnd    = tl["stagen"]->AddTask(&Hydro::ImplicitConduction, this, id.srctrms);
  // ... and the implicit TRANSVERSE radiative diffusion (<hydro>/rad_implicit_ang) right
  // after it: the two are operator-split from each other and from the hydro
  id.imptrc    = tl["stagen"]->AddTask(&Hydro::ImplicitTransverseConduction, this,
                                       id.impcnd);
  id.sendu_oa  = tl["stagen"]->AddTask(&Hydro::SendU_OA, this, id.imptrc);
  id.recvu_oa  = tl["stagen"]->AddTask(&Hydro::RecvU_OA, this, id.sendu_oa);
  id.restu     = tl["stagen"]->AddTask(&Hydro::RestrictU, this, id.recvu_oa);
  id.sendu     = tl["stagen"]->AddTask(&Hydro::SendU, this, id.restu);
  id.recvu     = tl["stagen"]->AddTask(&Hydro::RecvU, this, id.sendu);
  id.sendu_shr = tl["stagen"]->AddTask(&Hydro::SendU_Shr, this, id.recvu);
  id.recvu_shr = tl["stagen"]->AddTask(&Hydro::RecvU_Shr, this, id.sendu_shr);
  id.bcs       = tl["stagen"]->AddTask(&Hydro::ApplyPhysicalBCs, this, id.recvu_shr);
  id.prol      = tl["stagen"]->AddTask(&Hydro::Prolongate, this, id.bcs);
  id.c2p       = tl["stagen"]->AddTask(&Hydro::ConToPrim, this, id.prol);
  id.newdt     = tl["stagen"]->AddTask(&Hydro::NewTimeStep, this, id.c2p);

  // assemble "after_stagen" task list
  id.csend = tl["after_stagen"]->AddTask(&Hydro::ClearSend, this, none);
  // although RecvFlux/U functions check that all recvs complete, add ClearRecv to
  // task list anyways to catch potential bugs in MPI communication logic
  id.crecv = tl["after_stagen"]->AddTask(&Hydro::ClearRecv, this, id.csend);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Hydro::InitRecv
//! \brief Wrapper task list function to post non-blocking receives (with MPI), and
//! initialize all boundary receive status flags to waiting (with or without MPI).

TaskStatus Hydro::InitRecv(Driver *pdrive, int stage) {
  // post receives for U
  TaskStatus tstat = pbval_u->InitRecv(nhydro+nscalars);
  if (tstat != TaskStatus::complete) return tstat;

  // with SMR/AMR post receives for fluxes of U
  // do not post receives for fluxes when stage < 0 (i.e. ICs)
  if (pmy_pack->pmesh->multilevel && (stage >= 0)) {
    tstat = pbval_u->InitFluxRecv(nhydro+nscalars);
  }
  // see the note in SendFlux
  if (pmy_pack->pmesh->use_cubed_sphere && (stage >= 0)) {
    tstat = pbval_u->InitFluxSeamRecv(nhydro+nscalars);
  }
  if (tstat != TaskStatus::complete) return tstat;

  // with orbital advection post receives for U
  // only execute for (last stage) AND (3D OR 2d_r_phi)
  if (porb_u != nullptr) {
    if ((stage == pdrive->nexp_stages) &&
        (pmy_pack->pmesh->three_d || porb_u->shearing_box_r_phi)) {
      tstat = porb_u->InitRecv();
    }
  }
  if (tstat != TaskStatus::complete) return tstat;

  // with shearing box boundaries calculate x2-distance x1-boundarues have sheared and
  // with MPI post receives for U.
  // only execute if (3D OR 2d_r_phi)
  if (psbox_u != nullptr) {
    if (pmy_pack->pmesh->three_d || psbox_u->shearing_box_r_phi) {
      Real time = pmy_pack->pmesh->time;
      if (stage == pdrive->nexp_stages) {
        time += pmy_pack->pmesh->dt;
      }
      tstat = psbox_u->InitRecv(time);
    }
  }

  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::CopyCons
//! \brief Simple task list function that copies u0 --> u1 in first stage.  Extended to
//!  handle RK register logic at given stage

TaskStatus Hydro::CopyCons(Driver *pdrive, int stage) {
  if (stage == 1) {
    Kokkos::deep_copy(DevExeSpace(), u1, u0);
  } else {
    if (pdrive->integrator == "rk4") {
      // parallel loop to update u1 with u0 at later stages, only for rk4
      auto &indcs = pmy_pack->pmesh->mb_indcs;
      int is = indcs.is, ie = indcs.ie;
      int js = indcs.js, je = indcs.je;
      int ks = indcs.ks, ke = indcs.ke;
      int nmb1 = pmy_pack->nmb_thispack - 1;
      int nvar = nhydro + nscalars;
      auto &u0 = pmy_pack->phydro->u0;
      auto &u1 = pmy_pack->phydro->u1;
      Real &delta = pdrive->delta[stage-1];
      par_for("rk4_copy_cons", DevExeSpace(),0, nmb1, 0, nvar-1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(int m, int n, int k, int j, int i) {
        u1(m,n,k,j,i) += delta*u0(m,n,k,j,i);
      });
    }
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Hydro::Fluxes
//! \brief Wrapper task list function that calls everything necessary to compute fluxes
//! of conserved variables

TaskStatus Hydro::Fluxes(Driver *pdrive, int stage) {
  if (use_wellbalance_static_reconst_perturb) {
    RemoveWbVar(w0wb,w0);
  }
  // select which calculate_flux function to call based on rsolver_method
  if (rsolver_method == Hydro_RSolver::advect) {
    CalculateFluxes<Hydro_RSolver::advect>(pdrive, stage);
  } else if (rsolver_method == Hydro_RSolver::llf) {
    CalculateFluxes<Hydro_RSolver::llf>(pdrive, stage);
  } else if (rsolver_method == Hydro_RSolver::hlle) {
    CalculateFluxes<Hydro_RSolver::hlle>(pdrive, stage);
  } else if (rsolver_method == Hydro_RSolver::hllc) {
    CalculateFluxes<Hydro_RSolver::hllc>(pdrive, stage);
  } else if (rsolver_method == Hydro_RSolver::lhllc) {
    CalculateFluxes<Hydro_RSolver::lhllc>(pdrive, stage);
  } else if (rsolver_method == Hydro_RSolver::hllclm) {
    CalculateFluxes<Hydro_RSolver::hllclm>(pdrive, stage);
  } else if (rsolver_method == Hydro_RSolver::ausmpup) {
    CalculateFluxes<Hydro_RSolver::ausmpup>(pdrive, stage);
  } else if (rsolver_method == Hydro_RSolver::roe) {
    CalculateFluxes<Hydro_RSolver::roe>(pdrive, stage);
  } else if (rsolver_method == Hydro_RSolver::llf_sr) {
    CalculateFluxes<Hydro_RSolver::llf_sr>(pdrive, stage);
  } else if (rsolver_method == Hydro_RSolver::hlle_sr) {
    CalculateFluxes<Hydro_RSolver::hlle_sr>(pdrive, stage);
  } else if (rsolver_method == Hydro_RSolver::hllc_sr) {
    CalculateFluxes<Hydro_RSolver::hllc_sr>(pdrive, stage);
  } else if (rsolver_method == Hydro_RSolver::llf_gr) {
    CalculateFluxes<Hydro_RSolver::llf_gr>(pdrive, stage);
  } else if (rsolver_method == Hydro_RSolver::hlle_gr) {
    CalculateFluxes<Hydro_RSolver::hlle_gr>(pdrive, stage);
  }

  // <problem>/nan_report: the x1 WB reconstruction, captured before the pfloor mask
  if (nan_report && wbrec_nanrep_cnt != nullptr && wbrec_nanrep_lines < 400) {
    auto hwc = Kokkos::create_mirror_view(*wbrec_nanrep_cnt);
    Kokkos::deep_copy(hwc, *wbrec_nanrep_cnt);
    if (hwc(0) > 0) {
      auto h = Kokkos::create_mirror_view(*wbrec_nanrep_rec);
      Kokkos::deep_copy(h, *wbrec_nanrep_rec);
      ++wbrec_nanrep_lines;
      const int mb = static_cast<int>(h(0));
      std::cout << "### nan_report [wb_recon_x1] rank " << global_variable::my_rank
                << " cycle " << pmy_pack->pmesh->ncycle
                << " t = " << pmy_pack->pmesh->time << ": " << hwc(0)
                << " bad interface(s); first (m,k,j,i) = (" << mb << ","
                << static_cast<int>(h(1)) << "," << static_cast<int>(h(2)) << ","
                << static_cast<int>(h(3)) << ") gid = " << (pmy_pack->gids + mb)
                << " stage = " << static_cast<int>(h(4))
                << " wb_cache_rebuilt = " << static_cast<int>(h(5))
                << " | w0 d[i-1,i,i+1] = " << h(6) << "," << h(7) << "," << h(8)
                << " e = " << h(9) << "," << h(10) << "," << h(11)
                << " T = " << h(12) << "," << h(13) << "," << h(14)
                << " wder_p = " << h(15) << "," << h(16) << "," << h(17)
                << " phicc = " << h(18) << "," << h(19) << "," << h(20)
                << " | bg_d[im1,imh,i,iph,ip1] = " << h(21) << "," << h(22) << ","
                << h(23) << "," << h(24) << "," << h(25)
                << " bg_e = " << h(26) << "," << h(27) << "," << h(28) << ","
                << h(29) << "," << h(30)
                << " bg_p = " << h(31) << "," << h(32) << "," << h(33) << ","
                << h(34) << "," << h(35)
                << " | dl(IDPR) = " << h(36) << " dr(IDPR) = " << h(37)
                << " wl(IEN) = " << h(38) << " wr(IEN) = " << h(39)
                << " dl(IDG1) = " << h(40) << " dr(IDG1) = " << h(41)
                << " wl(IDN) = " << h(42) << " wr(IDN) = " << h(43)
                << " phif = " << h(44) << std::endl;
    }
  }
  // <problem>/nan_report: what the HLLC solver itself captured, then the Riemann fluxes
  // alone, before any diffusive flux is added
  if (nan_report && rsolv_nanrep_cnt != nullptr &&
      rsolv_nanrep_lines < 400) {
    auto hrc = Kokkos::create_mirror_view(*rsolv_nanrep_cnt);
    Kokkos::deep_copy(hrc, *rsolv_nanrep_cnt);
    if (hrc(0) > 0) {
      auto hrr = Kokkos::create_mirror_view(*rsolv_nanrep_rec);
      Kokkos::deep_copy(hrr, *rsolv_nanrep_rec);
      ++rsolv_nanrep_lines;
      const int mb = static_cast<int>(hrr(0));
      std::cout << "### nan_report [hllc_solver] rank " << global_variable::my_rank
                << " cycle " << pmy_pack->pmesh->ncycle
                << " t = " << pmy_pack->pmesh->time << ": " << hrc(0)
                << " bad face(s); first (m,k,j,i) = (" << mb << ","
                << static_cast<int>(hrr(1)) << "," << static_cast<int>(hrr(2)) << ","
                << static_cast<int>(hrr(3)) << ") gid = " << (pmy_pack->gids + mb)
                << " | L: d = " << hrr(4) << " vx = " << hrr(5) << " vy = " << hrr(6)
                << " vz = " << hrr(7) << " e = " << hrr(8) << " p = " << hrr(9)
                << " G1 = " << hrr(10)
                << " | R: d = " << hrr(11) << " vx = " << hrr(12) << " vy = " << hrr(13)
                << " vz = " << hrr(14) << " e = " << hrr(15) << " p = " << hrr(16)
                << " G1 = " << hrr(17)
                << " | cs_l = " << hrr(18) << " cs_r = " << hrr(19)
                << " S_L = " << hrr(20) << " S_R = " << hrr(21) << " S_M = " << hrr(22)
                << " p_star = " << hrr(23) << " ml = " << hrr(24) << " mr = " << hrr(25)
                << " E_l = " << hrr(26) << " E_r = " << hrr(27)
                << " fl.e = " << hrr(28) << " fr.e = " << hrr(29)
                << " | wL = " << hrr(30) << " wR = " << hrr(31) << " wC = " << hrr(32)
                << " flx(mx) = " << hrr(33) << " flx(d) = " << hrr(34)
                << " flx(E) = " << hrr(35) << std::endl;
    }
  }
  if (nan_report) {
    NanScanFlux(pmy_pack, uflx.x1f, 1, w0, peos->eos_data, "hydro_flux_x1");
    if (pmy_pack->pmesh->multi_d) {
      NanScanFlux(pmy_pack, uflx.x2f, 2, w0, peos->eos_data, "hydro_flux_x2");
    }
    if (pmy_pack->pmesh->three_d) {
      NanScanFlux(pmy_pack, uflx.x3f, 3, w0, peos->eos_data, "hydro_flux_x3");
    }
  }

  // Call FOFC if necessary.  It must run on the RIEMANN fluxes ALONE, before any
  // additive correction (diffusion, gravity flux, the well-balanced flux removal, and
  // the resistive fluxes in MHD): the published algorithm forms the trial update from
  // the Riemann fluxes only, and the fallback is meant to replace that flux, with the
  // corrections then applied on top of the first-order flux.
  if (pmy_pack->pradm1 != nullptr) {pmy_pack->pradm1->DetTrace("h_flxraw");}
  if (use_fofc) {
    FOFC(pdrive, stage);
    if (pmy_pack->pradm1 != nullptr) {pmy_pack->pradm1->DetTrace("h_flxfofc");}
  } else if (pmy_pack->pcoord->is_general_relativistic) {
    if (pmy_pack->pcoord->coord_data.bh_excise) {
      FOFC(pdrive, stage);
    }
  }

  // Restore the FULL primitives before any operator that reads them as a physical
  // state.  Under wellbalance_static_reconst the background was removed from w0 at the
  // top of this task so the Riemann reconstruction sees the deviation; everything below
  // -- the conduction and viscous fluxes, and every consumer downstream -- needs the
  // real state.  (FOFC above is the one operator that must run on the deviation: it adds
  // the background back itself, see hydro_fofc.cpp.)  A temperature computed from a
  // deviation density and energy is meaningless, and on a stratified star that is the
  // whole heat flux.
  if (use_wellbalance_static_reconst_perturb) AddWbVar(w0wb,w0);

  // Add diffusion fluxes
  if (pcond != nullptr) {
    // the angular cap (<hydro>/rad_cap_ang) needs the step it is capping, and the
    // angular fluxes are formed here, before the RK update ever sees beta_dt
    pcond->stage_beta_dt = (stage >= 1)
        ? (pdrive->beta[stage-1])*(pmy_pack->pmesh->dt) : 0.0;
    pcond->AddHeatFluxes(w0, peos->eos_data, uflx);
    if (nan_report) {
      NanScanFlux(pmy_pack, uflx.x1f, 1, w0, peos->eos_data, "after_conduction_flux");
    }
  }
  if (pvisc != nullptr) {
    pvisc->AddViscousFluxes(w0, peos->eos_data, uflx);
  }
  if (use_etotgrav) {
    AddGravFlux(phi0,uflx);
  }
  if (use_wellbalance_static) {
    RemoveWbFlux(pfacewb,uflx);
  }

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Hydro::SendFlux
//! \brief Wrapper task list function to pack/send restricted values of fluxes of
//! conserved variables at fine/coarse boundaries

TaskStatus Hydro::SendFlux(Driver *pdrive, int stage) {
  TaskStatus tstat = TaskStatus::complete;
  // Only execute BoundaryVaLUES function with SMR/SMR
  if (pmy_pack->pmesh->multilevel) {
    tstat = pbval_u->PackAndSendFluxCC(uflx);
  }
  // On the cubed sphere, reconcile the flux through a PANEL SEAM instead. A seam face is
  // one physical face held by two panels, each computing its own flux there from
  // interpolated ghosts; without this the update does not telescope across a seam and
  // mass and energy drift secularly. Mutually exclusive with the above: refinement and
  // the cubed sphere cannot both be on.
  if (pmy_pack->pmesh->use_cubed_sphere) {
    tstat = pbval_u->PackAndSendFluxSeamCC(uflx);
  }
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Hydro::RecvFlux
//! \brief Wrapper task list function to recv/unpack restricted values of fluxes of
//! conserved variables at fine/coarse boundaries

TaskStatus Hydro::RecvFlux(Driver *pdrive, int stage) {
  TaskStatus tstat = TaskStatus::complete;
  // Only execute BoundaryValues function with SMR/SMR
  if (pmy_pack->pmesh->multilevel) {
    tstat = pbval_u->RecvAndUnpackFluxCC(uflx);
  }
  // see the note in SendFlux
  if (pmy_pack->pmesh->use_cubed_sphere) {
    tstat = pbval_u->RecvAndUnpackFluxSeamCC(uflx);
  }
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Hydro::HydroSrcTerms
//! \brief Wrapper task list function to apply source terms to conservative vars
//! Note source terms must be computed using only primitives (w0), as the conserved
//! variables (u0) have already been partially updated when this fn called.

TaskStatus Hydro::HydroSrcTerms(Driver *pdrive, int stage) {
  if (pmy_pack->pradm1 != nullptr) {pmy_pack->pradm1->DetTrace("h_presrc");}
  Real beta_dt = (pdrive->beta[stage-1])*(pmy_pack->pmesh->dt);

  // <problem>/nan_report: this task runs immediately after RKUpdate, so u0 here is what
  // the flux divergence produced
  if (nan_report) NanScanCons(pmy_pack, u0, "after_RKUpdate");

  // Add physics source terms (must be computed from primitives)
  if (psrc != nullptr) psrc->ApplySrcTerms(w0, peos->eos_data,  beta_dt, u0);
  if (nan_report && psrc != nullptr) NanScanCons(pmy_pack, u0, "after_ApplySrcTerms");

  // Add shearing box source terms for cell-centered hydro variables
  if (psbox_u != nullptr) psbox_u->SourceTermsCC(w0, peos->eos_data, beta_dt, u0);

  // Add coordinate source terms in GR.  Again, must be computed with only primitives.
  if (pmy_pack->pcoord->is_general_relativistic) {
    pmy_pack->pcoord->CoordSrcTerms(w0, peos->eos_data, beta_dt, u0);
  }
    
  // Add coordinate source terms in curvi-linear grid.  Again, must be computed with only primitives.
  if (pmy_pack->pmesh->use_cubed_sphere) {
    pmy_pack->pcoord->SrcTermsGnomonicEquiangle(w0, wder, pwb, uflx, peos->eos_data,
                                                beta_dt, u0);
    if (nan_report) NanScanCons(pmy_pack, u0, "after_GnomonicSrc");
  }
  if (pmy_pack->pmesh->use_spherical_polar) {
    pmy_pack->pcoord->SrcTermsSphericalPolarHydro(w0, pwb, uflx, peos->eos_data,
                                                 beta_dt, u0);
  }

  // cs-open-1018 <hydro>/am_transport: the axial angular-momentum-conservative transport
  // correction (hydro_am.hpp), from this stage's fluxes and w0
  if (am_transport) {AmTransport(beta_dt, pdrive, stage);}
  if (frame_on) {FrameRotation(beta_dt, pdrive, stage);}
  if (am_vector && am_eta > 0.0 && stage == pdrive->nexp_stages) {AmSpinReturn();}

  // Add user source terms
  if (pmy_pack->pmesh->pgen->user_srcs) {
    (pmy_pack->pmesh->pgen->user_srcs_func)(pmy_pack->pmesh, beta_dt);
  }

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void Hydro::AmTransport
//! \brief cs-open-1018 <hydro>/am_transport: add T_cons - T_code (hydro_am.hpp) of every
//! active cell to its zonal momentum.  The source variants it does not replicate are a
//! FATAL (as the deep_hot_jupiter_rt problem/am_transport).

void Hydro::AmTransport(Real bdt, Driver *pdrive, int stage) {
  Mesh *pm = pmy_pack->pmesh;
  Coordinates *pco = pmy_pack->pcoord;
  if (!am_checked) {
    std::string bad = "";
    if (!(pm->use_cubed_sphere || pm->use_spherical_polar)) {
      bad = "a mesh that is neither cubed sphere nor spherical polar";
    }
    if (use_wellbalance_static) bad = "<hydro>/wellbalance_static";
    if (pm->use_cubed_sphere && pco->cs_wellbalanced_src) bad = "cs_wellbalanced_src";
    if (pm->use_spherical_polar && pco->sp_wellbalanced_src) bad = "sp_wellbalanced_src";
    if (pm->use_spherical_polar && pco->sp_cart_polar_momentum) {
      bad = "sp_cart_polar_momentum";
    }
    if (pm->use_spherical_polar && pco->sp_face_avg) bad = "sp_face_avg";
    if (!bad.empty()) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<hydro>/am_transport does not replicate the angular "
                << "source with " << bad << std::endl;
      std::exit(EXIT_FAILURE);
    }
    am_checked = true;
  }
  hydro::HydroAmV v;
  AmFillV(v);
  auto &ic = pm->mb_indcs;
  auto u0_ = u0;
  if (am_vector) {
    // the spin reservoir follows the RK register logic of u0 (u1 = u0 at stage 1)
    if (pdrive->integrator == "rk4") {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<hydro>/am_vector: integrator rk4 not supported"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (stage == 1) {Kokkos::deep_copy(DevExeSpace(), am_spin1, am_spin);}
    const Real g0 = pdrive->gam0[stage-1], g1 = pdrive->gam1[stage-1];
    auto sp_ = am_spin;
    auto sp1_ = am_spin1;
    auto sps_ = am_spins;
    par_for("hydro_am_vector", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1,
            ic.ks, ic.ke, ic.js, ic.je, ic.is, ic.ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real dT[3], rc[3], b1[3], b2[3], cc;
      hydro::HydroAmVecTorque(v, m, k, j, i, dT, rc, b1, b2, cc);
      const Real r2 = rc[0]*rc[0] + rc[1]*rc[1] + rc[2]*rc[2];
      const Real ir = 1.0/sqrt(r2);
      const Real sr = (dT[0]*rc[0] + dT[1]*rc[1] + dT[2]*rc[2])*ir;   // rhat . dT
      // the tangential momentum whose torque is the perpendicular part of dT:
      // dp = (dT x rc)/(|rc|^2 V)
      Real dp[3];
      hydro::HydroAmCross(dT, rc, dp);
      const Real iv = 1.0/v.vol(m,k,j,i);
      for (int q = 0; q < 3; ++q) {dp[q] *= iv/r2;}
      const Real p1 = dp[0]*b1[0] + dp[1]*b1[1] + dp[2]*b1[2];
      const Real p2 = dp[0]*b2[0] + dp[1]*b2[1] + dp[2]*b2[2];
      // covariant on the cell basis (cs: m_a = p.e_a; sp: orthonormal)
      u0_(m,IM2,k,j,i) += p1*bdt;
      u0_(m,IM3,k,j,i) += p2*bdt;
      sps_(m,k,j,i) = sp_(m,k,j,i);   // the stage's own state (frame_exact)
      sp_(m,k,j,i) = g0*sp_(m,k,j,i) + g1*sp1_(m,k,j,i) + sr*bdt;
      (void) cc;
    });
    return;
  }
  par_for("hydro_am_transport", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1,
          ic.ks, ic.ke, ic.js, ic.je, ic.is, ic.ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real dtq, c2, c3;
    if (!hydro::HydroAmTorque(v, m, k, j, i, dtq, c2, c3)) return;
    u0_(m,IM2,k,j,i) += dtq*c2*bdt;
    u0_(m,IM3,k,j,i) += dtq*c3*bdt;
  });
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Hydro::ImplicitConduction
//! \brief Wrapper task that applies the implicit radial radiative diffusion
//! (<hydro>/rad_implicit_x1).  A no-op unless that flag is set; see
//! Conduction::ImplicitRadialUpdate for what it solves and why.

TaskStatus Hydro::ImplicitConduction(Driver *pdrive, int stage) {
  if (pcond == nullptr) return TaskStatus::complete;
  if (!(pcond->rad_implicit_x1)) return TaskStatus::complete;
  if (stage < 1) return TaskStatus::complete;
  Real beta_dt = (pdrive->beta[stage-1])*(pmy_pack->pmesh->dt);
  pcond->ImplicitRadialUpdate(u0, peos->eos_data, beta_dt);
  runaway_scan::Scan(pmy_pack->pmesh, "implicit_conduction");
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Hydro::ImplicitTransverseConduction
//! \brief Wrapper task that applies the implicit transverse radiative diffusion
//! (<hydro>/rad_implicit_ang).  A no-op unless that flag is set; see
//! Conduction::ImplicitTransverseUpdate for what it solves and why.  With
//! <hydro>/rad_sts_all it instead runs the UNIFIED operator over x1, x2 and x3 in one
//! RKL1 loop (Conduction::StsConductionUpdate), and the radial task above does nothing.

TaskStatus Hydro::ImplicitTransverseConduction(Driver *pdrive, int stage) {
  if (pcond == nullptr) return TaskStatus::complete;
  if (!(pcond->rad_implicit_ang)) return TaskStatus::complete;
  if (stage < 1) return TaskStatus::complete;
  // the operator has been moved out of the stage, to run with the radiation column in
  // the same operator-split step (Conduction::rad_tr_split_out)
  if (pcond->rad_tr_split_out) return TaskStatus::complete;
  Real beta_dt = (pdrive->beta[stage-1])*(pmy_pack->pmesh->dt);
  // <hydro>/rad_sts_once: ONE call per cycle, after the LAST stage and over the
  // FULL dt, instead of one per stage over beta_dt.  An RKL1 super-step of tau costs
  // sqrt(tau) substages, so one call over dt is cheaper than two over dt/2 by ~30 %.
  // The splitting between this operator and the hydro is FIRST ORDER either way -- the
  // operator is not part of the RK right-hand side, so running it every stage does not
  // raise the order, it only halves the size of the first-order term.
  if (pcond->rad_sts_once) {
    if (stage != pdrive->nexp_stages) return TaskStatus::complete;
    beta_dt = pmy_pack->pmesh->dt;
  }
  // <hydro>/rad_sts_all replaces BOTH implicit conduction tasks with this one: the radial
  // task above is a no-op (rad_implicit_x1 is refused alongside it) and the RKL1 loop
  // here runs over all three directions at once.  Same slot in the stage either way.
  if (pcond->rad_sts_all) {
    pcond->StsConductionUpdate(u0, peos->eos_data, beta_dt);
  } else {
    pcond->ImplicitTransverseUpdate(u0, peos->eos_data, beta_dt);
  }
  runaway_scan::Scan(pmy_pack->pmesh, "implicit_transverse_conduction");
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Hydro::SendU_OA
//! \brief Wrapper task list function to pack/send data for orbital advection

TaskStatus Hydro::SendU_OA(Driver *pdrive, int stage) {
  TaskStatus tstat = TaskStatus::complete;
  if (porb_u != nullptr) {
    // only execute if (last stage) AND (3D OR 2d_r_phi)
    if ((stage == pdrive->nexp_stages) &&
        (pmy_pack->pmesh->three_d || porb_u->shearing_box_r_phi)) {
      tstat = porb_u->PackAndSendCC(u0);
    }
  }
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Hydro::RecvU_OA
//! \brief Wrapper task list function to recv/unpack data for orbital advection
//! Orbital remap is performed in this step.

TaskStatus Hydro::RecvU_OA(Driver *pdrive, int stage) {
  TaskStatus tstat = TaskStatus::complete;
  if (porb_u != nullptr) {
    // only execute if (last stage) AND (3D OR 2d_r_phi)
    if ((stage == pdrive->nexp_stages) &&
        (pmy_pack->pmesh->three_d || porb_u->shearing_box_r_phi)) {
      tstat = porb_u->RecvAndUnpackCC(u0, recon_method);
    }
  }
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Hydro::RestrictU
//! \brief Wrapper task list function to restrict conserved vars

TaskStatus Hydro::RestrictU(Driver *pdrive, int stage) {
  if (pmy_pack->pradm1 != nullptr) {pmy_pack->pradm1->DetTrace("h_postsrc");}
  // Only execute Mesh function with SMR/SMR
  if (pmy_pack->pmesh->multilevel) {
    pmy_pack->pmesh->pmr->RestrictCC(u0, coarse_u0);
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Hydro::SendU
//! \brief Wrapper task list function to pack/send cell-centered conserved variables

TaskStatus Hydro::SendU(Driver *pdrive, int stage) {
  TaskStatus tstat = pbval_u->PackAndSendCC(u0, coarse_u0);
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Hydro::RecvU
//! \brief Wrapper task list function to receive/unpack cell-centered conserved variables

TaskStatus Hydro::RecvU(Driver *pdrive, int stage) {
  TaskStatus tstat = pbval_u->RecvAndUnpackCC(u0, coarse_u0);
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Hydro::SendU_Shr
//! \brief Wrapper task list function to pack/send data for shearing box boundaries

TaskStatus Hydro::SendU_Shr(Driver *pdrive, int stage) {
  TaskStatus tstat = TaskStatus::complete;
  if (psbox_u != nullptr) {
    // only execute if (3D OR 2d_r_phi)
    if (pmy_pack->pmesh->three_d || psbox_u->shearing_box_r_phi) {
      tstat = psbox_u->PackAndSendCC(u0, recon_method);
    }
  }
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Hydro::RecvU_Shr
//! \brief Wrapper task list function to recv/unpack data for shearing box boundaries
//! Orbital remap is performed in this step.

TaskStatus Hydro::RecvU_Shr(Driver *pdrive, int stage) {
  TaskStatus tstat = TaskStatus::complete;
  if (psbox_u != nullptr) {
    // only execute if (3D OR 2d_r_phi)
    if (pmy_pack->pmesh->three_d || psbox_u->shearing_box_r_phi) {
      tstat = psbox_u->RecvAndUnpackCC(u0);
    }
  }
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Hydro::ApplyPhysicalBCs
//! \brief Wrapper task list function to call funtions that set physical and user BCs,

TaskStatus Hydro::ApplyPhysicalBCs(Driver *pdrive, int stage) {
  // do not apply BCs if domain is strictly periodic
  if (pmy_pack->pmesh->strictly_periodic) return TaskStatus::complete;

  // physical BCs
  // with etotgrav the copy BCs re-reference rho*Phi to the ghost's own potential
  pbval_u->HydroBCs((pmy_pack), (pbval_u->u_in), u0,
                    use_etotgrav ? phicc0 : DvceArray4D<Real>());

  // user BCs
  if (pmy_pack->pmesh->pgen->user_bcs) {
    (pmy_pack->pmesh->pgen->user_bcs_func)(pmy_pack->pmesh);
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Hydro::Prolongate
//! \brief Wrapper task list function to prolongate conserved (or primitive) variables
//! at fine/coarse boundaries with SMR/AMR

TaskStatus Hydro::Prolongate(Driver *pdrive, int stage) {
  if (pmy_pack->pmesh->multilevel) {  // only prolongate with SMR/AMR
    pbval_u->FillCoarseInBndryCC(u0, coarse_u0);
    if (pmy_pack->pmesh->pmr->prolong_prims) {
      pbval_u->ConsToPrimCoarseBndry(coarse_u0, coarse_w0);
      pbval_u->ProlongateCC(w0, coarse_w0);
      pbval_u->PrimToConsFineBndry(w0, u0);
    } else {
      pbval_u->ProlongateCC(u0, coarse_u0);
    }
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void Hydro::RTOpSplitBvals
//! \brief a COMPLETE ghost-zone update of u0 for an operator that ran OUTSIDE the stage's
//! own communication window (the Strang half-steps of RTStrangSplit).
//!
//! WHY THIS EXISTS.  Every in-stage source writes ACTIVE cells only and is followed,
//! inside the same stage, by RestrictU -> SendU -> RecvU -> ApplyPhysicalBCs ->
//! Prolongate -> ConToPrim, so what the next flux evaluation reads in the ghost zones is
//! the state the source produced.  An operator applied in "before_/after_timeintegrator",
//! or at the head of stage 1 before the stage's own exchange, has no such tail: the
//! MeshBlock halos and -- decisively for a radiating box -- the PHYSICAL x1 wall and top
//! ghosts still hold the state from before it ran.  The first flux divergence of the
//! cycle is then evaluated across a discontinuity of one whole operator step at every
//! block face and at both x1 boundaries.  In the optically thin lid the two-stream
//! changes the top cell by order unity in a step, so the top ghost (built by the user BC
//! from the last active cell) is order-unity wrong once per cycle.
//!
//! The sequence is exactly the stage's own tail, run synchronously: the receives are
//! spun on rather than retried by the task machinery, which costs the MPI overlap once
//! per cycle and nothing else.
//!
//! \param recv_posted  true when "before_stagen" has ALREADY posted the receives for U
//!                     (i.e. this is being called from inside "stagen"); the posted
//!                     receives are then consumed here instead of being posted twice.
//! \param repost       true to post a fresh set of receives afterwards, for the SendU /
//!                     RecvU the rest of the stage still has to run.

void Hydro::RTOpSplitBvals(bool recv_posted, bool repost) {
  if (!recv_posted) {
    while (pbval_u->InitRecv(nhydro+nscalars) != TaskStatus::complete) {}
  }
  if (pmy_pack->pmesh->multilevel) {
    pmy_pack->pmesh->pmr->RestrictCC(u0, coarse_u0);
  }
  while (pbval_u->PackAndSendCC(u0, coarse_u0) != TaskStatus::complete) {}
  while (pbval_u->RecvAndUnpackCC(u0, coarse_u0) != TaskStatus::complete) {}
  (void) ApplyPhysicalBCs(nullptr, 0);
  (void) Prolongate(nullptr, 0);
  (void) ConToPrim(nullptr, 0);
  while (pbval_u->ClearSend() != TaskStatus::complete) {}
  while (pbval_u->ClearRecv() != TaskStatus::complete) {}
  if (repost) {
    while (pbval_u->InitRecv(nhydro+nscalars) != TaskStatus::complete) {}
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Hydro::EnableC2PTrack / C2PTrack
//! \brief floor bookkeeping of ConToPrim (see c2p_track in hydro.hpp and
//! utils/c2p_track.hpp).  Enabled by a problem generator; c2p_acc is zeroed by its owner.

void Hydro::EnableC2PTrack() {
  c2p_track::Enable(pmy_pack, c2p_acc, c2p_ref_e, c2p_ref_d);
  c2p_track = true;
}

void Hydro::C2PTrack(Driver *pdrive, int stage, bool before) {
  if (before) {
    c2p_track::Before(pmy_pack, u0, c2p_ref_e, c2p_ref_d);
  } else {
    // the split radiation's own conversion runs from RTOpSplitBvals (no driver)
    const bool split = c2p_split || (pdrive == nullptr);
    c2p_track::After(pmy_pack, u0, c2p_ref_e, c2p_ref_d, c2p_acc,
                     c2p_track::Weight(pdrive, stage, split), split);
  }
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Hydro::ConToPrimSplit
//! \brief ConToPrim as the operator-split task pair runs it (after RTStrangSplit, outside
//! the RK stages): identical, except that the floor bookkeeping weights it by 1.

TaskStatus Hydro::ConToPrimSplit(Driver *pdrive, int stage) {
  c2p_split = true;
  TaskStatus ts = ConToPrim(pdrive, stage);
  c2p_split = false;
  return ts;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Hydro::RTStrangSplit
//! \brief the Strang half-step of ProblemGenerator::user_split_func.  It runs once in
//! "before_timeintegrator" and once in "after_timeintegrator", each with HALF the
//! CYCLE dt -- pmesh->dt is fixed for the whole cycle, NewTimeStep running only after
//! the second call -- so the two halves are equal and the split is second order.  Each
//! is followed by ConToPrim in the same list so the stage (and the next cycle) sees
//! primitives consistent with the conserved state the source just changed.
//!
//! A null hook returns immediately: every run that does not enrol one is bitwise
//! unchanged.

TaskStatus Hydro::RTStrangSplit(Driver *pdrive, int stage) {
  Mesh *pm = pmy_pack->pmesh;
  if (pm->pgen == nullptr || pm->pgen->user_split_func == nullptr) {
    return TaskStatus::complete;
  }
  // user_split_once: ONE full-dt call, after the last RK stage.  The
  // "before_timeintegrator" half (stage = 0, see Driver::Execute) is skipped entirely,
  // so the source is applied once per cycle on the final state instead of being blended
  // by the SSP-RK stage weights.
  if (pm->pgen->user_split_once) {
    if (stage == 0) return TaskStatus::complete;
    (pm->pgen->user_split_func)(pm, pm->dt);
    RTOpSplitBvals(false, false);
    return TaskStatus::complete;
  }
  (pm->pgen->user_split_func)(pm, 0.5*(pm->dt));
  // the half-step wrote ACTIVE cells only, outside any stage: bring the MeshBlock halos
  // and the physical x1 ghosts up to the state it produced before the stage reads them
  RTOpSplitBvals(false, false);
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Hydro::ConToPrim
//! \brief Wrapper task list function to call ConsToPrim over entire mesh (including gz)

TaskStatus Hydro::ConToPrim(Driver *pdrive, int stage) {
  if (c2p_track) C2PTrack(pdrive, stage, true);
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int &ng = indcs.ng;
  int n1m1 = indcs.nx1 + 2*ng - 1;
  int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
  int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
  // A restart that carried the internal energy (kEintRstMagic, pgen/pgen.hpp): the
  // frozen first conversion takes w0(IEN) from the file and hands u0(IEN) back exactly
  // as the file had it.  Re-deriving the internal energy from the file's u0 goes through
  // (E - rho phi) + rho phi - rho phi with etotgrav, which is not the E - rho phi of the
  // straight run's last inversion; the restart then stopped being bitwise.
  const bool eint_fix = c2p_freeze_derived && c2p_eint_rst;
  DvceArray4D<Real> ekeep;
  if (eint_fix) {
    Kokkos::realloc(ekeep, eint_rst.extent(0), eint_rst.extent(1), eint_rst.extent(2),
                    eint_rst.extent(3));
    Kokkos::deep_copy(ekeep, Kokkos::subview(u0, std::make_pair(0,
                      static_cast<int>(eint_rst.extent(0))), static_cast<int>(IEN),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (use_etotgrav) {
    RemoveGravEtot(phicc0, u0, 0, n1m1, 0, n2m1, 0, n3m1);
  }
  // remembered across the branch below because the cubed-sphere half of the inversion,
  // which owns the cache on that grid, has to be frozen by the SAME call
  const bool frozen = c2p_freeze_derived;
  if (c2p_freeze_derived) {
    // the first conversion of a RESTARTED run under a general EOS: the temperature and
    // the derived (p, Gamma_1) came out of the restart file and must survive this call
    // untouched, or the restart is not a bitwise continuation.  One call only.
    static_cast<GeneralHydro*>(peos)->ConsToPrimFrozen(u0, w0, 0, n1m1, 0, n2m1,
                                                       0, n3m1);
    c2p_freeze_derived = false;
  } else {
    peos->ConsToPrim(u0, w0, false, 0, n1m1, 0, n2m1, 0, n3m1);
  }
  // On the cubed sphere the conserved momentum is COVARIANT while the primitive velocity
  // must be CONTRAVARIANT, and the two differ by the non-orthogonal metric of the
  // gnomonic tangent basis. ConsToPrim cannot know that, so redo the velocity and the
  // internal energy with the metric. See Coordinates::GnomonicEquiangleRaiseVel.
  if (pmy_pack->pmesh->use_cubed_sphere) {
    if (frozen) {
      // ...and on this grid THIS routine, not ConsToPrim, is what re-solves wtemp/wder
      // and re-applies the deferred floors, so the frozen conversion has to reach here
      // too or the cache the restart file restored is overwritten a moment later.  See
      // coordinates/gnomonic_raisevel_frozen.cpp.
      pmy_pack->pcoord->GnomonicEquiangleRaiseVelFrozen(u0, w0, 0, n1m1, 0, n2m1,
                                                        0, n3m1);
    } else {
      pmy_pack->pcoord->GnomonicEquiangleRaiseVel(u0, w0, peos->eos_data, wder, wtemp,
                                                  0, n1m1, 0, n2m1, 0, n3m1);
    }
  }
  if (use_etotgrav) {
    AddGravEtot(phicc0, u0, 0, n1m1, 0, n2m1, 0, n3m1);
  }
  if (eint_fix) {
    const int nmbe = static_cast<int>(eint_rst.extent(0));
    Kokkos::deep_copy(Kokkos::subview(w0, std::make_pair(0,nmbe), static_cast<int>(IEN),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), eint_rst);
    Kokkos::deep_copy(Kokkos::subview(u0, std::make_pair(0,nmbe), static_cast<int>(IEN),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ekeep);
    eint_rst = DvceArray4D<Real>();
    c2p_eint_rst = false;
  }
  if (c2p_track) C2PTrack(pdrive, stage, false);
  runaway_scan::Scan(pmy_pack->pmesh, "ConToPrim_floors");
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Hydro::ClearSend
//! \brief Wrapper task list function that checks all MPI sends have completed. Used in
//! TaskList and in Driver::InitBoundaryValuesAndPrimitives()
//! If stage=(last stage):      clears sends of U, Flx_U, U_OA, U_Shr
//! If (last stage)>stage>=(0): clears sends of U, Flx_U,       U_Shr
//! If stage=(-1):              clears sends of U
//! If stage=(-4):              clears sends of                 U_Shr

TaskStatus Hydro::ClearSend(Driver *pdrive, int stage) {
  TaskStatus tstat;
  // check sends of U complete
  if ((stage >= 0) || (stage == -1)) {
    tstat = pbval_u->ClearSend();
    if (tstat != TaskStatus::complete) return tstat;
  }

  // with SMR/AMR check sends of restricted fluxes of U complete
  // do not check flux send for ICs (stage < 0)
  // On the cubed sphere the PANEL-SEAM flux exchange (SendFlux) uses the same flux_req
  // and must be waited on too: unwaited, every stage leaked one MPI request per off-rank
  // seam neighbour and Cray MPICH aborted once its request pool ran out (DeltaAI, dhj
  // at 2 ranks after ~8000 cycles: "Assertion failed ... req != NULL" in MPI_Isend).
  if ((pmy_pack->pmesh->multilevel || pmy_pack->pmesh->use_cubed_sphere) &&
      (stage >= 0)) {
    tstat = pbval_u->ClearFluxSend();
    if (tstat != TaskStatus::complete) return tstat;
  }

  // with orbital advection check sends of U complete
  // only execute when (shearing box defined) AND (last stage) AND (3D OR 2d_r_phi)
  if (porb_u != nullptr) {
    if ((stage == pdrive->nexp_stages) &&
        (pmy_pack->pmesh->three_d || porb_u->shearing_box_r_phi)) {
      tstat = porb_u->ClearSend();
      if (tstat != TaskStatus::complete) return tstat;
    }
  }

  // with shearing box boundaries check sends of U complete
  // only execute when (shearing box defined) AND (stage>=0 or -4) AND (3D OR 2d_r_phi)
  if (psbox_u != nullptr) {
    if (((stage >= 0) || (stage == -4)) &&
        (pmy_pack->pmesh->three_d || psbox_u->shearing_box_r_phi)) {
      tstat = psbox_u->ClearSend();
      if (tstat != TaskStatus::complete) return tstat;
    }
  }

  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Hydro::ClearRecv
//! \brief Wrapper task list function that checks all MPI receives have completed. Used in
//! TaskList and in Driver::InitBoundaryValuesAndPrimitives()
//! If stage=(last stage):      clears recvs of U, Flx_U, U_OA, U_Shr
//! If (last stage)>stage>=(0): clears recvs of U, Flx_U,       U_Shr
//! If stage=(-1):              clears recvs of U
//! If stage=(-4):              clears recvs of                 U_Shr

TaskStatus Hydro::ClearRecv(Driver *pdrive, int stage) {
  TaskStatus tstat;
  // check receives of U complete
  if ((stage >= 0) || (stage == -1)) {
    tstat = pbval_u->ClearRecv();
    if (tstat != TaskStatus::complete) return tstat;
  }

  // with SMR/AMR check receives of restricted fluxes of U complete
  // do not check flux receives when stage < 0 (i.e. ICs)
  // (and on the cubed sphere, the panel-seam flux receives; see ClearSend)
  if ((pmy_pack->pmesh->multilevel || pmy_pack->pmesh->use_cubed_sphere) &&
      (stage >= 0)) {
    tstat = pbval_u->ClearFluxRecv();
    if (tstat != TaskStatus::complete) return tstat;
  }

  // with orbital advection check receives of U complete
  // only execute when (shearing box defined) AND (last stage) AND (3D OR 2d_r_phi)
  if (porb_u != nullptr) {
    if ((stage == pdrive->nexp_stages) &&
        (pmy_pack->pmesh->three_d || porb_u->shearing_box_r_phi)) {
      tstat = porb_u->ClearRecv();
      if (tstat != TaskStatus::complete) return tstat;
    }
  }

  // with shearing box boundaries check receives of U complete
  // only execute when (shearing box defined) AND (stage>=0 or -4) AND (3D OR 2d_r_phi)
  if (psbox_u != nullptr) {
    if (((stage >= 0) || (stage == -4)) &&
        (pmy_pack->pmesh->three_d || psbox_u->shearing_box_r_phi)) {
      tstat = psbox_u->ClearRecv();
      if (tstat != TaskStatus::complete) return tstat;
    }
  }

  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn void Hydro::AmSpinReturn
//! \brief cs-open-1018 <hydro>/am_spin_return (see hydro.hpp).  Face f between cell c
//! and its lateral neighbour n (directions rhat_c, rhat_n, d = |rhat_c - rhat_n|,
//! s = (rhat_c - rhat_n)/d): c takes dL = lam s, n takes -lam s, lam = -eta (S_c + S_n)
//! d / (16 (1 - rhat_c.rhat_n)) (S the spin), so that each face lowers both spins by
//! ~eta (S_c + S_n)/16.  A cell carries dL exactly: rhat.dL to its spin, the rest as the
//! tangential momentum (dL x rc)/(|rc|^2 V).  Faces at a physical or periodic lateral
//! boundary (where the two sides do not share a geometry) take nothing.

void Hydro::AmSpinReturn() {
  Mesh *pm = pmy_pack->pmesh;
  Coordinates *pco = pmy_pack->pcoord;
  auto &ic = pm->mb_indcs;
  const int is = ic.is, ie = ic.ie, js = ic.js, je = ic.je, ks = ic.ks, ke = ic.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  // the spins in the ghosts (mirror cells at a seam)
  auto sp_ = am_spin;
  auto sw_ = spw;
  par_for("hydro_am_spin_cp", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    sw_(m,0,k,j,i) = sp_(m,k,j,i);
  });
  while (pbval_spin->InitRecv(1) == TaskStatus::incomplete) {}
  while (pbval_spin->PackAndSendCC(spw, spw_c) == TaskStatus::incomplete) {}
  while (pbval_spin->RecvAndUnpackCC(spw, spw_c) == TaskStatus::incomplete) {}
  while (pbval_spin->ClearSend() == TaskStatus::incomplete) {}
  while (pbval_spin->ClearRecv() == TaskStatus::incomplete) {}
  const bool cs = pm->use_cubed_sphere;
  const Real eta = am_eta;
  auto u0_ = u0;
  auto vol = pco->volume;
  auto x1v = pco->x1v;
  auto x2v = pco->x2v;
  auto x3v = pco->x3v;
  auto &size = pmy_pack->pmb->mb_size;
  auto &pan = pmy_pack->pmb->mb_panel;
  auto &bcs = pmy_pack->pmb->mb_bcs;
  const int nx2 = ic.nx2, nx3 = ic.nx3;
  par_for("hydro_am_spin_ret", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // the direction of cell (k', j') of this block (extended panel / sp coordinates)
    auto dirn = [&](const int kq, const int jq, Real *q, Real *b1, Real *b2) {
      if (cs) {
        const Real xi = 0.25*M_PI*CellCenterX(jq-js, nx2, size.d_view(m).x2min,
                                              size.d_view(m).x2max);
        const Real et = 0.25*M_PI*CellCenterX(kq-ks, nx3, size.d_view(m).x3min,
                                              size.d_view(m).x3max);
        cubed_sphere::PanelToCart(pan.d_view(m), xi, et, q);
        cubed_sphere::PanelTangents(pan.d_view(m), xi, et, b1, b2);
      } else {
        const Real th = x2v(m,jq), ph = x3v(m,kq);
        q[0] = sin(th)*cos(ph); q[1] = sin(th)*sin(ph); q[2] = cos(th);
        b1[0] = cos(th)*cos(ph); b1[1] = cos(th)*sin(ph); b1[2] = -sin(th);
        b2[0] = -sin(ph); b2[1] = cos(ph); b2[2] = 0.0;
      }
    };
    Real rh[3], b1[3], b2[3];
    dirn(k, j, rh, b1, b2);
    const Real r = x1v(m,i);
    Real dL[3] = {0.0, 0.0, 0.0};
    const Real sc = sw_(m,0,k,j,i);
    for (int f = 0; f < 4; ++f) {
      const int dj = (f == 0) ? -1 : ((f == 1) ? 1 : 0);
      const int dk = (f == 2) ? -1 : ((f == 3) ? 1 : 0);
      const bool edge = (f == 0 && j == js) || (f == 1 && j == je) ||
                        (f == 2 && k == ks) || (f == 3 && k == ke);
      bool seam = false;
      if (edge) {
        const BoundaryFace bf = (f == 0) ? BoundaryFace::inner_x2
                              : (f == 1) ? BoundaryFace::outer_x2
                              : (f == 2) ? BoundaryFace::inner_x3 : BoundaryFace::outer_x3;
        const BoundaryFlag fl = bcs.d_view(m,bf);
        if (fl == BoundaryFlag::panel) {
          seam = true;
        } else if (fl != BoundaryFlag::block) {
          continue;   // physical / periodic: nothing crosses
        }
      }
      Real rn[3], t1[3], t2[3];
      if (seam) {
        // the mirror cell: rhat_c reflected across the seam plane (normal = gradient of
        // the panel coordinate at the face, which holds the origin)
        const int jf = (f == 1) ? (j + 1) : j, kf = (f == 3) ? (k + 1) : k;
        Real a[3], b[3], n0[3];
        cubed_sphere::PanelFrame(pan.d_view(m), a, b, n0);
        Real nf[3];
        if (f < 2) {
          const Real x = tan(0.25*M_PI*LeftEdgeX(jf-js, nx2, size.d_view(m).x2min,
                                                 size.d_view(m).x2max));
          for (int q = 0; q < 3; ++q) {nf[q] = a[q] - x*n0[q];}
        } else {
          const Real y = tan(0.25*M_PI*LeftEdgeX(kf-ks, nx3, size.d_view(m).x3min,
                                                 size.d_view(m).x3max));
          for (int q = 0; q < 3; ++q) {nf[q] = b[q] - y*n0[q];}
        }
        const Real nn = 1.0/sqrt(nf[0]*nf[0] + nf[1]*nf[1] + nf[2]*nf[2]);
        for (int q = 0; q < 3; ++q) {nf[q] *= nn;}
        const Real pr = rh[0]*nf[0] + rh[1]*nf[1] + rh[2]*nf[2];
        for (int q = 0; q < 3; ++q) {rn[q] = rh[q] - 2.0*pr*nf[q];}
      } else {
        dirn(k + dk, j + dj, rn, t1, t2);
      }
      Real sv[3];
      for (int q = 0; q < 3; ++q) {sv[q] = rh[q] - rn[q];}
      const Real d = sqrt(sv[0]*sv[0] + sv[1]*sv[1] + sv[2]*sv[2]);
      const Real omc = 1.0 - (rh[0]*rn[0] + rh[1]*rn[1] + rh[2]*rn[2]);
      if (!(d > 0.0) || !(omc > 0.0)) continue;
      const Real sn = sw_(m,0,k+dk,j+dj,i);
      const Real lam = -eta*(sc + sn)*d/(16.0*omc);
      for (int q = 0; q < 3; ++q) {dL[q] += lam*sv[q]/d;}
    }
    // this cell's share: radial part to the spin, the rest to the tangential momentum
    const Real srad = dL[0]*rh[0] + dL[1]*rh[1] + dL[2]*rh[2];
    Real rc[3] = {r*rh[0], r*rh[1], r*rh[2]};
    Real dp[3];
    hydro::HydroAmCross(dL, rc, dp);
    const Real iv = 1.0/(vol(m,k,j,i)*r*r);
    Real p1 = 0.0, p2 = 0.0;
    for (int q = 0; q < 3; ++q) {p1 += dp[q]*iv*b1[q]; p2 += dp[q]*iv*b2[q];}
    u0_(m,IM2,k,j,i) += p1;
    u0_(m,IM3,k,j,i) += p2;
    sp_(m,k,j,i) += srad;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void Hydro::FrameRotation
//! \brief cs-open-1018 <hydro>/frame_omega: the rotating-frame force (HydroAmFrame) on the
//! cell momentum and its work (the centrifugal part) on the total energy; under
//! <hydro>/am_vector also the L-conserving torque correction (perpendicular part in the
//! tangential momentum, radial part in the spin; the energy is not touched by it).

void Hydro::FrameRotation(Real bdt, Driver *pdrive, int stage) {
  Mesh *pm = pmy_pack->pmesh;
  if (!(pm->use_cubed_sphere || pm->use_spherical_polar)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<hydro>/frame_omega needs the cubed sphere or spherical "
              << "polar" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const int nst = pdrive->nexp_stages;
  const bool fex = frame_exact && am_vector;
  if (frame_exact && (!am_vector || nst > 2)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<hydro>/frame_exact needs <hydro>/am_vector and integrator "
              << "rk1 | rk2" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  hydro::HydroAmV v;
  AmFillV(v);
  auto &ic = pm->mb_indcs;
  auto u0_ = u0;
  auto u1_ = u1;
  const Real o0 = frame_om[0], o1 = frame_om[1], o2 = frame_om[2];
  const bool amc = am_vector;
  auto sp_ = am_spin;
  auto sp1_ = am_spin1;
  auto sps_ = am_spins;
  // the radiation's L (F / (c chat)) of the cell is part of L_abs: the frame terms book
  // its rotation on the gas (the M1 equations carry no frame terms); u0 F covariant on
  // the cs panel basis (raised with cos_cell), orthonormal on sp
  const bool hrad = amc && (pmy_pack->pradm1 != nullptr);
  auto ur_ = hrad ? pmy_pack->pradm1->u0 : DvceArray5D<Real>();
  const Real icc = hrad ? 1.0/(pmy_pack->pradm1->c_light*pmy_pack->pradm1->chat) : 0.0;
  if (fex && frame_ln.extent_int(0) == 0) {
    Kokkos::realloc(frame_ln, pmy_pack->nmb_thispack, 3, ic.nx3 + 2*ic.ng,
                    ic.nx2 + 2*ic.ng, ic.nx1 + 2*ic.ng);
  }
  auto ln_ = frame_ln;
  const bool rec = fex && (stage == 1);       // record L^n
  const bool fin = fex && (stage == nst);     // add the exact-rotation remainder
  const Real dt = pm->dt;
  // the exact rotation R(-Omega dt) (Rodrigues) and A = -Omega x
  const Real om = sqrt(o0*o0 + o1*o1 + o2*o2);
  const Real kx = (om > 0.0) ? o0/om : 0.0, ky = (om > 0.0) ? o1/om : 0.0;
  const Real kz = (om > 0.0) ? o2/om : 0.0;
  const Real cth = cos(-om*dt), sth = sin(-om*dt);
  const bool one = (nst == 1);
  par_for("hydro_frame_rot", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1,
          ic.ks, ic.ke, ic.js, ic.je, ic.is, ic.ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real omv[3] = {o0, o1, o2};
    Real f[3], dT[3], rc[3], rh[3], b1[3], b2[3], cc;
    hydro::HydroAmFrame(v, omv, amc, m, k, j, i, f, dT, rc, rh, b1, b2, cc);
    const Real v1 = v.w0(m,IVX,k,j,i), v2 = v.w0(m,IVY,k,j,i), v3 = v.w0(m,IVZ,k,j,i);
    Real pw = 0.0;
    for (int q = 0; q < 3; ++q) {pw += f[q]*(v1*rh[q] + v2*b1[q] + v3*b2[q]);}
    Real fr = 0.0, f1 = 0.0, f2 = 0.0;
    for (int q = 0; q < 3; ++q) {fr += f[q]*rh[q]; f1 += f[q]*b1[q]; f2 += f[q]*b2[q];}
    u0_(m,IM1,k,j,i) += fr*bdt;
    u0_(m,IM2,k,j,i) += f1*bdt;
    u0_(m,IM3,k,j,i) += f2*bdt;
    u0_(m,IEN,k,j,i) += pw*bdt;
    if (amc) {
      const Real dv = v.vol(m,k,j,i);
      const Real s2 = 1.0 - cc*cc;
      // a covariant cell vector (radial, x2, x3) -> Cartesian
      auto cart = [&](const Real a1, const Real b2c, const Real b3c, Real *o) {
        const Real a2 = (b2c - cc*b3c)/s2, a3 = (b3c - cc*b2c)/s2;
        for (int q = 0; q < 3; ++q) {o[q] = a1*rh[q] + a2*b1[q] + a3*b2[q];}
      };
      // the spin and radiation parts of L_abs (current): -Omega x (S rhat + L_rad)
      Real lx[3];
      // (the stage's state: the spin before this stage's RK combination)
      for (int q = 0; q < 3; ++q) {lx[q] = sps_(m,k,j,i)*rh[q];}
      if (hrad) {
        Real fc[3], lr[3];
        cart(ur_(m,radm1::M1_F1,k,j,i), ur_(m,radm1::M1_F2,k,j,i),
             ur_(m,radm1::M1_F3,k,j,i), fc);
        for (int q = 0; q < 3; ++q) {fc[q] *= dv*icc;}
        hydro::HydroAmCross(rc, fc, lr);
        for (int q = 0; q < 3; ++q) {lx[q] += lr[q];}
      }
      Real ol[3];
      hydro::HydroAmCross(omv, lx, ol);
      for (int q = 0; q < 3; ++q) {dT[q] -= ol[q];}
      Real dl[3];
      for (int q = 0; q < 3; ++q) {dl[q] = dT[q]*bdt;}
      if (rec) {
        // L_abs^n of the cell: gas V rc x (m + rho Omega x rc) from the step-start
        // conserved state u1 (covariant m), spin^n (am_spin1), radiation (u0 before the
        // step's solves)
        Real mc[3], orc[3], pa[3], la[3];
        cart(u1_(m,IM1,k,j,i), u1_(m,IM2,k,j,i), u1_(m,IM3,k,j,i), mc);
        hydro::HydroAmCross(omv, rc, orc);
        const Real rho = u1_(m,IDN,k,j,i);
        for (int q = 0; q < 3; ++q) {pa[q] = (mc[q] + rho*orc[q])*dv;}
        hydro::HydroAmCross(rc, pa, la);
        for (int q = 0; q < 3; ++q) {
          la[q] += sp1_(m,k,j,i)*rh[q] + (lx[q] - sps_(m,k,j,i)*rh[q]);
          ln_(m,q,k,j,i) = la[q];
        }
      }
      if (fin) {
        Real l[3] = {ln_(m,0,k,j,i), ln_(m,1,k,j,i), ln_(m,2,k,j,i)};
        // R l (Rodrigues about k by -|Omega| dt)
        const Real kl = kx*l[0] + ky*l[1] + kz*l[2];
        const Real kxl[3] = {ky*l[2] - kz*l[1], kz*l[0] - kx*l[2], kx*l[1] - ky*l[0]};
        Real al[3], aal[3];
        {Real t[3]; hydro::HydroAmCross(omv, l, t);
         for (int q = 0; q < 3; ++q) {al[q] = -t[q];}}
        {Real t[3]; hydro::HydroAmCross(omv, al, t);
         for (int q = 0; q < 3; ++q) {aal[q] = -t[q];}}
        for (int q = 0; q < 3; ++q) {
          const Real kq = (q == 0) ? kx : ((q == 1) ? ky : kz);
          const Real rl = l[q]*cth + kxl[q]*sth + kq*kl*(1.0 - cth);
          dl[q] += rl - l[q] - dt*al[q] - (one ? 0.0 : 0.5*dt*dt*aal[q]);
        }
      }
      const Real r2 = rc[0]*rc[0] + rc[1]*rc[1] + rc[2]*rc[2];
      Real dp[3];
      hydro::HydroAmCross(dl, rc, dp);
      const Real iv = 1.0/(dv*r2);
      Real p1 = 0.0, p2 = 0.0;
      for (int q = 0; q < 3; ++q) {p1 += dp[q]*iv*b1[q]; p2 += dp[q]*iv*b2[q];}
      u0_(m,IM2,k,j,i) += p1;
      u0_(m,IM3,k,j,i) += p2;
      sp_(m,k,j,i) += dl[0]*rh[0] + dl[1]*rh[1] + dl[2]*rh[2];
    }
  });
}

//! \brief cs-open-1018: the views of the hydro angular-momentum kernels (hydro_am.hpp)
void Hydro::AmFillV(hydro::HydroAmV &v) {
  Mesh *pm = pmy_pack->pmesh;
  Coordinates *pco = pmy_pack->pcoord;
  auto &ic = pm->mb_indcs;
  v.w0 = w0;
  v.eos = peos->eos_data;
  v.gen = v.eos.IsGeneral();
  v.wder = wder;
  v.f1 = uflx.x1f;  v.f2 = uflx.x2f;  v.f3 = uflx.x3f;
  v.a1 = pco->area.x1f;  v.a2 = pco->area.x2f;  v.a3 = pco->area.x3f;
  v.vol = pco->volume;
  v.xovr = pco->x_ov_rD;
  v.yovr = pco->y_ov_rC;
  v.ccell = pco->cos_cell;
  v.scell = pco->sin_cell;
  v.x1v = pco->x1v;  v.x1f = pco->xx1f;
  v.x2v = pco->x2v;  v.x2f = pco->xx2f;
  v.x3v = pco->x3v;  v.x3f = pco->xx3f;
  v.size = pmy_pack->pmb->mb_size;
  v.panel = pmy_pack->pmb->mb_panel;
  v.is = ic.is;  v.js = ic.js;  v.ks = ic.ks;  v.nx2 = ic.nx2;  v.nx3 = ic.nx3;
  v.cs = pm->use_cubed_sphere;
}

} // namespace hydro
