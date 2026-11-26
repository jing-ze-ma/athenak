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
//! \brief the implicit stage solve of the clean scheme, in the staggered layout (E in
//! cells, F on faces)

TaskStatus RadiationVet::CleanSolve(Driver *pdriver, int stage) {
  if (!cln_init) CleanInit();
  // the vet_gd direction set first (a restart's first solve reads it before any build)
  if (vgd_on && !vray_ready) VetRayInit();
  // sc rays' ring on the cubed sphere: the union's directions and weights too
  if (vgd_on && cln_cssc && cln_ring && cln_fvn == 0) {
    CleanFvSets(pmy_pack->pmesh->time, pmy_pack->pmesh->ncycle - 1);
    // sc ring: its two sets' intensities from the restart file (zero from a file of
    // another ray method, until the first build)
    if (cln_cssc) CleanScUnion2();
  }
  if (cln_cssc && !cln_ring) {   // clean_rays = sc: the sweep's set, its cln_fvi
    cln_fvd = vgd_dir;
    cln_fvn = vgd_n;
  }
  // restart from a file with the intensities (version 5): their face half ranges, set by
  // set
  if (cln_cssc && cln_fh.extent_int(0) == 0 && cln_fvi.extent_int(0) > 0) {
    for (int q = 0; q < cln_fvi.extent_int(1)/vgd_n; ++q) {
      CleanFaceHalf(cln_fvd, q*vgd_n, q == 0, q*vgd_n);
    }
    cln_fvi = DvceArray5D<Real>();
  }
  FluidRef fl = FluidRef::Get(pmy_pack);
  cln_feedback = gas_feedback && fl.on;

  // 1. starting face fluxes on the first cycle (restarts carry them)
  if (!cln_faces_set) {
    cln_faces_set = true;
    if (pmy_pack->pmesh->ncycle == 0) CleanInitFaceFluxes();
  }

  // 2. closure (from the gas before the stage increments), old vector, face coefficients.
  //    clean_closure = iterate with no closure in place takes the predicted one
  if (cln_extrap > 0) {
    CleanClosureExtrap();
  } else if (!cln_closit || !cln_cip) {
    CleanClosure();
  }
  if (cln_ring == 2) CleanRingSolveSet();
  CleanWoBind();
  cln_xdh = false;   // clean_stag_xdonor_freeze: no hold before this solve's passes
  CleanOldVector();
  // (clean_stag_opac_update with gas coupling: CleanStageOpacity takes both, with the
  // first iterate's opacities)
  const bool opu = cln_opt && cln_feedback;
  if (!opu) CleanFaceCoefficients();
  // donors' E_old; the bracket's gate q
  if (cln_gflux || cln_consist) CleanGhostFill(cln_w, CLW_EO, 1);
  if (cln_consist && !opu) CleanStagBracket();
  // 3. Newton iteration on the gas coupling, Picard passes on the lagged slopes and (with
  //    clean_closure = iterate, at most clean_maxpic) the closure; one linear solve
  //    without these.  clean_closure = extrap(2): the stage-1 solve takes the closure of
  //    its own first iterate once (pass 1), else the step is first order in time.
  const bool clos1 = (cln_extrap > 0) && (t2_solve == VET_T2S_STAGE1);
  const bool iterate = cln_feedback || cln_closit || clos1;
  const bool t2stage = (t2_solve == VET_T2S_STAGE1 || t2_solve == VET_T2S_STAGE2);
  int nit = 0, nit_sum = 0, nsolve = 0;
  std::array<int, 2> nfirst = {0, 0};
  bool converged = !iterate;
  const int maxnl = iterate ? cln_maxnl : 1;
  // lagged slopes settle only to the linear tolerance.  clean_conv = budget: the pass's
  // tolerance (CleanNewtonUpdate), the first solve's from the floor (restarts bitwise)
  Real tol_nl = fmax(cln_tol_nl, 10.0*cln_tol);
  cln_tolb = tol_nl;
  // budget norm and time error only in gas-coupled Newton passes with K1 = f(U^n) known;
  // otherwise the strict test
  cln_bon = cln_budget && iterate && cln_feedback && (t2_ok || cln_k1ok) &&
            t2_solve != VET_T2S_NONE && !cln_vzero &&
            t2k1.extent_int(0) >= pmy_pack->nmb_thispack;
  cln_ewr = -1.0;
  // a change that stopped decreasing within 100 x both its round-off level (cln_fe,
  // cln_fv) and tol_nl is converged too
  auto settled = [&](const Real d, const Real d0, const Real f) {
    return d < tol_nl || (d > 0.5*d0 && d < 100.0*fmin(f, tol_nl));
  };
  Real dnl = 0.0, dnl0 = -1.0;   // the last pass's change and the one before
  {
    CleanFirstIterate();
    cln_fsk = false;
    cln_csyok = false;   // lateral sums' coefficients once per closure
    if (opu) CleanStageOpacity();
    if (cln_feedback) CleanAdvectiveCoefficients();
    // the drag's face flux from the first iterate's beta and emission (a restart, failed
    // stage or discarded attempt does not carry the last solve's)
    if (cln_feedback) CleanRightHandSide(true);
    if (cln_feedback) CleanStagDrag(true);
    // clean_stag_xdonor: q of the first iterate with its drag share b (CLW_L holds E_k's
    // lateral sums after the first drag)
    CleanStagXDonor(true, !cln_feedback);
    converged = !iterate;
    cln_gpic = false;   // clean_gas_newton_switch: the solve's history
    std::vector<Real> rhist;
    Real de0 = 1.0e300, dv0 = 1.0e300, trel = 0.0;
    for (int inl = 0; inl < maxnl; ++inl) {
      const bool clos = (cln_closit && inl < cln_maxpic) || (clos1 && inl == 1);
      // opacities of the iterate's T, in passes < clean_stag_opac_passes
      const bool opr = opu && inl > 0 && inl < cln_optn;
      if (clos) {
        cln_csyok = false;
        CleanClosureIterate();
        if (!opr) CleanFaceCoefficients();   // (else in CleanStageOpacity)
        if (cln_consist && !opr) CleanStagBracket();
        if (cln_feedback) CleanAdvectiveCoefficients();
      }
      if (opr) {
        CleanStageOpacity();
        CleanStagDrag(false);
        CleanAdvectiveCoefficients();
      }
      // clean_stag_xdonor: q of the iterate; CLW_L holds E_k's lateral sums unless the
      // closure, bracket or opacity refresh rewrote them (and on the cubed sphere)
      if (inl > 0) CleanStagXDonor(false, clos || opr || !cln_feedback,
                                   cln_xdfz > 0 && inl >= cln_xdfz);
      // after the xdonor: the lateral sums of E_k are its
      CleanRightHandSide(false, !(inl > 0 && cln_xd_on));
      // the first iterate's T from the pass's gas row, then its RHS
      if (inl == 0 && cln_feedback) {
        CleanNewtonUpdate();
        CleanRightHandSide();
      }
      CleanPreconditionerSetup();
      // budget: a pass whose iterate gives the next closure is solved in the strict norm
      const bool cnext = (cln_closit && inl + 1 < cln_maxpic) || (clos1 && inl == 0);
      if (cnext) cln_ewr = (cln_c1b && !dtr_redo) ? -3.0 : -2.0;
      nit = CleanBiCGStab();
      if (cnext) cln_ewr = -1.0;
      nit_sum += std::abs(nit);
      ++nsolve;
      cln_nsolve += 1.0;
      cln_nit += std::abs(nit);
      if (nit < 0) dnl = std::numeric_limits<Real>::infinity();
      // implicit_lin_maxit: a linear solve that ends at its cap leaves its iterate to the
      // pass, which is then not the last one
      const bool lp = nit < 0 && inl + 1 < maxnl;
      if (!iterate || (nit < 0 && !lp)) break;
      const bool fbl = lp;
      if (!t2stage && inl + 1 == maxnl) {   // for CleanWorstCell
        auto sv = Kokkos::subview(cln_w.s, Kokkos::ALL, std::make_pair(CLW_EK,
                                  CLW_EK + 1), Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        Kokkos::realloc(cln_wci, sv.extent(0), 1, sv.extent(2), sv.extent(3),
                        sv.extent(4));
        Kokkos::deep_copy(cln_wci, sv);
      }
      const Real de = CleanNewtonUpdate();
      if (cln_bon) tol_nl = cln_tolb;
      Real dv = 0.0;
      if (cln_feedback) dv = CleanVelocityUpdate();
      dnl0 = (inl + 1 == maxnl) ? dnl : -1.0;
      dnl = fmax(de, dv);
      // clean_gas_newton_switch = W: a stalled solve (pass + 1 >= _min and > W, the
      // minimum over the last W passes >= clean_stall_fac x the one before, up at least
      // twice in the window) drops the gas Newton slope for its remaining passes: the gas
      // T takes the bracketed root of its row, a contraction (Newton cycles where c_v
      // jumps, e.g. H ionisation).  Deterministic in the passes: restarts bitwise.
      if (!cln_gpic) {
        rhist.push_back(dnl);
        constexpr int nw = 4;
        if (inl + 1 >= 6) {
          Real mwin = 1.0e300, mpre = 1.0e300;
          int nup = 0;
          for (int q = 0; q <= inl; ++q) {
            if (q > inl - nw) {
              mwin = fmin(mwin, rhist[q]);
              if (rhist[q] > rhist[q-1]) ++nup;
            } else {
              mpre = fmin(mpre, rhist[q]);
            }
          }
          if (mwin >= 0.5*mpre && nup >= 2) {
            cln_gpic = true;
            cln_ngnsw += 1.0;
            if (global_variable::my_rank == 0) {
              std::printf("CLEAN GNSW cycle %d t2s %d pass %d: stalled, the gas Newton "
                          "slope off (%.0f so far)\n", pmy_pack->pmesh->ncycle, t2_solve,
                          inl, cln_ngnsw);
            }
          }
        }
      }
      // strict test: stop at 0.1 x the first pass's change, at most 1e-4
      if (!cln_bon && inl == 0) trel = fmin(0.1*dnl, 1.0e-4);
      tol_nl = fmax(tol_nl, trel);
      // clean_verbose >= 2: each pass's changes against the tolerance it is tested with
      if (cln_verbose > 1 && global_variable::my_rank == 0) {
        std::printf("CLEAN PASS cycle %d t2s %d pass %d it %d de %.3e dv %.3e tol %.3e "
                    "fe %.3e fv %.3e%s%s\n", pmy_pack->pmesh->ncycle, t2_solve, inl,
                    nit, de, dv, tol_nl, cln_fe, cln_fv, clos ? " clos" : "",
                    opr ? " opac" : "");
      }
      if (fbl) {
        de0 = de;
        dv0 = dv;
        continue;
      }
      if (settled(de, de0, cln_fe) && settled(dv, dv0, cln_fv) && !(clos1 && inl == 0)) {
        converged = true;
        break;
      }
      de0 = de;
      dv0 = dv;
    }
    if (iterate) CleanCopyIterate();
    // R1, R2 of the solution, counted
    CleanStagRealizable();
    nfirst = {cln_nr1, cln_nr2};
  }
  if (!converged) {
    nit = -std::abs(nit) - 1;
  }

  // 4. positivity of the solution
  const int nneg = CleanCountNegative();
  if (nit < 0) cln_nfail += 1.0;
  if (nneg > 0) cln_nneg += 1.0;
  // an unconverged backward-Euler solve fails the cycle: with <time>/dt_retry it is
  // retaken at a smaller dt (a failed stage is redone by backward Euler via t2_fail);
  // without, the run stops if the last change exceeds implicit_resid_fatal
  if (nit < 0 && !t2stage) {
    const bool fatal = !dtr_on && impl_res_fatal > 0.0 && !(dnl <= impl_res_fatal);
    if (dtr_on) {
      if (!dtr_fail) {
        dtr_why = 1;
        dtr_resid = -1.0;
      }
      dtr_fail = true;
    }
    if (global_variable::my_rank == 0) {
      std::printf("CLEAN WARNING cycle %d: backward-Euler solve not converged (%s, %d "
                  "so far; last change %.3e)\n", pmy_pack->pmesh->ncycle,
                  dtr_on ? "dt_retry" : (fatal ? "no dt_retry: > implicit_resid_fatal" :
                  "no dt_retry: accepted"), static_cast<int>(cln_nfail), dnl);
    }
    CleanWorstCell(dnl, dnl0);
    if (fatal) CleanFatal("<rad_vet> clean: a backward-Euler solve did not converge and "
                          "its last change exceeds <rad_vet>/implicit_resid_fatal");
  }
  if (global_variable::my_rank == 0) {
    std::printf("CLEAN cycle %d t2s %d it %d (sum %d) pos 0 mark 0 all 0 neg %d "
                "r1 %d r2 %d first %d %d guard %d solves %d\n",
                pmy_pack->pmesh->ncycle, t2_solve, nit, nit_sum, nneg, cln_nr1, cln_nr2,
                nfirst[0], nfirst[1], cln_nguard, nsolve);
  }
  if (cln_mclc.extent_int(0) == 1) {   // last pass's face momentum terms at the cap
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cln_mclc);
    const Real n = CleanReduce(h(0), false);
    if (global_variable::my_rank == 0 && cln_verbose > 0 && n > 0.0) {
      std::printf(
          "CLEAN MCAP cycle %d t2s %d: %.0f face momentum terms at the cap h chi_f "
          "E in the last pass\n", pmy_pack->pmesh->ncycle, t2_solve, n);
    }
  }

  // 5. new state: face fluxes, E, cell fluxes, gas, stage slopes
  // a failed solve or a probe does not change the energy flux's memory
  cln_dok = (nit >= 0) && !t2_inprobe;
  CleanWriteBack();
  // a write-back that leaves gas eint < 0 is rejected like a failed solve (stage redone
  // by backward Euler, a backward-Euler step by dt_retry)
  if (cln_feedback && cln_neint > 0 && !t2stage && dtr_on) {
    if (!dtr_fail) {
      dtr_why = 4;   // reason eint<=0
      dtr_resid = cln_neint;
    }
    dtr_fail = true;
  }
  if (cln_feedback && cln_neint > 0 && !t2stage && !dtr_on && !t2_inprobe &&
      global_variable::my_rank == 0) {
    std::printf("CLEAN WARNING cycle %d t2s %d: backward-Euler step accepted with %d "
                "cells of gas eint below 0 / its floor, E < 0 or not finite (no "
                "dt_retry)\n", pmy_pack->pmesh->ncycle, t2_solve, cln_neint);
  }
  if (t2_solve == VET_T2S_STAGE1 || t2_solve == VET_T2S_STAGE2) {
    // implicit_pos_floor(_cons): an E < 0 raised by the conservative write-back floor
    // does not fail the stage (counted)
    const bool flr = impl_pos_floor;
    t2_fail = (nit < 0) || (nneg > 0 && !flr) || (cln_feedback && cln_neint > 0) ||
              (t2_solve == VET_T2S_STAGE1 && pmy_pack->pmesh->ncycle == t2_dbg_fail);
    if (cln_lzbook && cln_feedback && t2_solve == VET_T2S_STAGE2 && !t2_fail) {
      CleanLzBook(0.0, false);   // the accepted stage 2's L_z, pending until here
    }
  }
  // K1 = f(U^n) for the next predictor and time error: from an accepted backward-Euler
  // or stage-2 solve; a failed one leaves none (BE redo and auto-BE steps keep the last)
  if (t2_solve == VET_T2S_BESTORE || t2_solve == VET_T2S_STAGE2) {
    cln_k1ok = (nit >= 0) && !dtr_fail && !(t2_solve == VET_T2S_STAGE2 && t2_fail);
  }
  // force_reference = clean: the reference is the force of the converged, accepted first
  // coupled solve (each iterate recaptures it; a failed one is redone)
  if (force_ref == VET_FREF_CLEAN && cln_feedback && !cln_vzero && nit >= 0 &&
      nneg == 0 && !(t2stage && t2_fail) && !dtr_fail && !t2_inprobe) cln_fref_set = true;

  cln_xd_on = false;   // clean_stag_xdonor: q is of this solve's passes only

  // 6. clean_closure = extrap(2): ray moments of the new state
  if (cln_extrap > 0 && !t2_fail) CleanClosureAfterSolve();
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
