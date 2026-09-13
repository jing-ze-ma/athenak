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
  id.sendu_oa  = tl["stagen"]->AddTask(&Hydro::SendU_OA, this, id.impcnd);
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
  if (use_fofc) {
    FOFC(pdrive, stage);
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

  // Add user source terms
  if (pmy_pack->pmesh->pgen->user_srcs) {
    (pmy_pack->pmesh->pgen->user_srcs_func)(pmy_pack->pmesh, beta_dt);
  }

  return TaskStatus::complete;
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
  pbval_u->HydroBCs((pmy_pack), (pbval_u->u_in), u0);

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
//! \fn TaskList Hydro::ConToPrim
//! \brief Wrapper task list function to call ConsToPrim over entire mesh (including gz)

TaskStatus Hydro::ConToPrim(Driver *pdrive, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int &ng = indcs.ng;
  int n1m1 = indcs.nx1 + 2*ng - 1;
  int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
  int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
  if (use_etotgrav) {
    RemoveGravEtot(phicc0, u0, 0, n1m1, 0, n2m1, 0, n3m1);
  }
  peos->ConsToPrim(u0, w0, false, 0, n1m1, 0, n2m1, 0, n3m1);
  // On the cubed sphere the conserved momentum is COVARIANT while the primitive velocity
  // must be CONTRAVARIANT, and the two differ by the non-orthogonal metric of the
  // gnomonic tangent basis. ConsToPrim cannot know that, so redo the velocity and the
  // internal energy with the metric. See Coordinates::GnomonicEquiangleRaiseVel.
  if (pmy_pack->pmesh->use_cubed_sphere) {
    pmy_pack->pcoord->GnomonicEquiangleRaiseVel(u0, w0, peos->eos_data, wder, wtemp,
                                                0, n1m1, 0, n2m1, 0, n3m1);
  }
  if (use_etotgrav) {
    AddGravEtot(phicc0, u0, 0, n1m1, 0, n2m1, 0, n3m1);
  }
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
  if (pmy_pack->pmesh->multilevel && (stage >= 0)) {
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
  if (pmy_pack->pmesh->multilevel && (stage >= 0)) {
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

} // namespace hydro
