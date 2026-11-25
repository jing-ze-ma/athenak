//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_time2.cpp
//! \brief time_scheme = hesdirk2: the IMEX stages (Heun hydro, ESDIRK implicit
//! radiation), their state and its restart channels.
//! Implicit part: A = [[0,0,0],[1-g,g,0],[1/2,b2,g]], g = 1 - 1/sqrt 2, b2 = g/(2(1-g)),
//! first slope FSAL.  A stage is the backward-Euler ImplicitSolve with dt -> g dt and the
//! old vector rhs = (stage start) + t2inc, t2inc = (1-g) dt K1 (stage 1) or
//! dt (g/2 K1 + (b2 - g/2) K2) (stage 2); its slope is K = (Y - rhs)/(g dt).  Without a
//! valid K1 (first step, restart without the slope, failed stage) the step is backward
//! Euler and stores K1 = (Y - rhs)/dt.

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "coordinates/coordinates.hpp"
#include "driver/driver.hpp"
#include "hydro/hydro.hpp"
#include "hydro/hydro_am.hpp"
#include "eos/eos.hpp"
#include "rad_vet/vet_fluid.hpp"
#include "utils/deep_copy_across.hpp"
#include "rad_vet/rad_vet.hpp"
#include "utils/recipe.hpp"
#include "rad_vet/rad_vet_parfor.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_clean.hpp"
#include "rad_vet/rad_vet_opacity.hpp"
#include "utils/eint_from_cons.hpp"
#include "utils/imex_sdirk2.hpp"

namespace radvet {

using deep_copy_across::DeepCopyAcross;

namespace {
const Real kT2G = imex_sdirk2::Gamma();   // shared tableau
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::Time2Init
//! \brief reads time_scheme (its default set by ImplicitInit) and the time2_* keys, and
//! allocates the stage state

void RadiationVet::Time2Init(ParameterInput *pin) {
  time_scheme = VET_TIME_BE;
  if (pin->DoesParameterExist("rad_vet", "time_scheme")) {
    std::string s = pin->GetString("rad_vet", "time_scheme");
    if (s.compare("be") == 0) {
      time_scheme = VET_TIME_BE;
    } else if (s.compare("hesdirk2") == 0) {
      time_scheme = VET_TIME_HESDIRK2;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/time_scheme = '" << s
                << "' is not a choice (be | hesdirk2)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // vet_every > 1 is not accepted under vet_sc
  const std::string kev = VetKeyName(pin, "vet_every", {"vet_sc_every"});
  if (vet_sc && pin->DoesParameterExist("rad_vet", kev) &&
      pin->GetInteger("rad_vet", kev) > 1) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<rad_vet>/" << kev << " > 1 (vet_sc) is removed; only 1 "
              << "is accepted" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // be with t2_beonly runs as hesdirk2 with every step on the backward-Euler path
  if (time_scheme == VET_TIME_BE && t2_beonly) {
    time_scheme = VET_TIME_HESDIRK2;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> time_scheme = be: the hesdirk2 backward-Euler path on "
                   "every "
                << "step (be_legacy = true for the old step)" << std::endl;
    }
  } else if (time_scheme == VET_TIME_HESDIRK2 && t2_beonly &&
             global_variable::my_rank == 0) {
    std::cout << "<rad_vet> time2_be_only = true: every step takes the backward-Euler "
                 "path"
              << std::endl;
  }
  if (time_scheme == VET_TIME_BE) {t2_beonly = false; return;}

  std::string integ = pin->GetString("time", "integrator");
  if (integ.compare("rk2") != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<rad_vet>/time_scheme = hesdirk2 needs <time>/integrator "
              << "= rk2 (its explicit part IS the Heun hydro)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // stage coefficients (t2_tab, utils/imex_sdirk2.hpp):
  //   hesdirk2: a21 = 1-g, Heun average, then g/2 K1 + (b2 - g/2) K2
  //   trbdf2 (c = (0, 2g, 1), w = sqrt(2)/4): a21 = g, hydro weights 1/4 (Y2) and
  //     3/4 (U^n), then (w - g/4)(K1 + K2)
  //   sdirk2: no explicit K1 (Time2GasTime)
  {
    const Real g = kT2G;
    t2_tab = 0;
    {
      // default: sdirk2 (mhd: not wired, hesdirk2)
      const bool sd = pmy_pack->pmhd == nullptr &&
          recipe::Fix(pin, "rad_vet", "time2_tableau", "sdirk2", "hesdirk2");
      std::string tb = pin->GetOrAddString("rad_vet", "time2_tableau",
                                           sd ? "sdirk2" : "hesdirk2");
      if (tb.compare("hesdirk2") == 0) {
        t2_tab = 0;
      } else if (tb.compare("trbdf2") == 0) {
        t2_tab = 1;
      } else if (tb.compare("sdirk2") == 0) {
        t2_tab = 2;
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<rad_vet>/time2_tableau = '" << tb
                  << "' is not a choice (hesdirk2 | trbdf2)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    const imex_sdirk2::Tableau tb = imex_sdirk2::Make(t2_tab);
    (void) g;
    t2_a21 = tb.a21;
    t2_w0 = tb.w0;
    t2_w1 = tb.w1;
    t2_c31 = tb.c31;
    t2_c32 = tb.c32;
    if (t2_tab == 2) {
      t2_gtm = pin->GetOrAddBoolean("rad_vet", "time2_gtime_mom", true);
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_vet> time2_tableau = sdirk2: no explicit K1, Heun hydro"
                  << (t2_gtm ? "" : ", time2_gtime_mom = false") << std::endl;
      }
      // time2_gtime_slope: Time2GasTime's slope-limited exchange.  Default true (with
      // time2_gtime_mom)
      const bool nm = pin->DoesParameterExist("rad_vet", "time2_gtime_slope");
      t2_gts = pin->GetOrAddBoolean("rad_vet", "time2_gtime_slope", t2_gtm &&
          recipe::Fix(pin, "rad_vet", "time2_gtime_slope", "true", "false"));
      recipe::OldStored(pin, "rad_vet", "time2_gtime_slope", "false", "true");
      if (t2_gts && !t2_gtm) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<rad_vet>/time2_gtime_slope = true needs "
                  << "time2_gtime_mom = true (a restart file may store false): name "
                  << "time2_gtime_mom=true, or time2_gtime_slope=false" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_vet> time2_gtime_slope = " << (t2_gts ? "true" : "false")
                  << (nm ? "" : " (default)") << std::endl;
      }
      // time2_gtime_pos: the positive form of the flux state (always with
      // time2_gtime_slope).  Default true
      const bool np = pin->DoesParameterExist("rad_vet", "time2_gtime_pos");
      t2_gtp = pin->GetOrAddBoolean("rad_vet", "time2_gtime_pos",
          recipe::Fix(pin, "rad_vet", "time2_gtime_pos", "true", "false"));
      recipe::OldStored(pin, "rad_vet", "time2_gtime_pos", "false", "true");
      if (t2_gts && !t2_gtp) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<rad_vet>/time2_gtime_slope = true needs "
                  << "time2_gtime_pos = true" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_vet> time2_gtime_pos = " << (t2_gtp ? "true" : "false")
                  << (np ? "" : " (default)") << std::endl;
      }
    }
    if (t2_tab == 2 && coupling && pmy_pack->pmhd != nullptr) {  // no mhd path
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/time2_tableau = sdirk2 is not wired for mhd"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (t2_tab == 1) {
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_vet> time2_tableau = trbdf2: TR-BDF2 implicit stages (stage "
                  << "order 2), RK2 hydro with c2 = 2 - sqrt(2)" << std::endl;
      }
    }
  }
  // closure = vet_col sets tau_closure too; its tensor is built at U^n
  if (tau_closure && !vet_col) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<rad_vet>/time_scheme = hesdirk2 is not wired for "
              << "closure = tau" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // time2_enth_be: the backward-Euler solves of a hesdirk2 run take the stage solves'
  // enthalpy face rule too (else the operator changes for one step: an O(dt) error)
  t2_enth_be = pin->GetOrAddBoolean("rad_vet", "time2_enth_be",
      recipe::Fix(pin, "rad_vet", "time2_enth_be", "true", "false"));
  // time2_vet_extrap: stage solves use D* = D^n + (dt/dt_prev)(D^n - D^{n-1})
  t2_vext = false;
  if (pin->DoesParameterExist("rad_vet", "time2_vet_extrap")) {
    t2_vext = pin->GetBoolean("rad_vet", "time2_vet_extrap");
  }
  // time2_vet (aliases time2_vet_sc / time2_vet_col by closure): when the formal-solution
  // tensor of the two stage solves is built.  Both stages sit at t^{n+1}, so D(U^n) is an
  // O(dt) lag (first order in time wherever the tensor changes).
  //   lag     : D(U^n) for both stages (vet_sc extrap: D* of time2_vet_extrap)
  //   predict : one build per step at P = (stage-1 start) + dt K1, t^{n+1} + O(dt^2)
  //   rebuild : predict for stage 1, a second build at Y1 for stage 2
  t2_vsmode = 0;
  const std::string kvs = VetKeyName(pin, "time2_vet", {"time2_vet_sc"});
  const bool vs_named = pin->DoesParameterExist("rad_vet", kvs);
  if (vet_sc && vs_named && pin->GetString("rad_vet", kvs).compare("extrap") == 0) {
    t2_vext = true;
  }
  if (vet_sc && (vs_named || !t2_vext)) {
    std::string vm = pin->GetOrAddString("rad_vet", kvs,
                                         global_variable::restart_run ? "lag"
                                                                      : "predict");
    if (vm.compare("lag") == 0 || vm.compare("extrap") == 0) {
      t2_vsmode = 0;
    } else if (vm.compare("predict") == 0) {
      t2_vsmode = 1;
    } else if (vm.compare("rebuild") == 0) {
      t2_vsmode = 2;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/" << kvs << " = '" << vm
                << "' is not a choice (lag | predict | rebuild | extrap)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (t2_vsmode != 0 && t2_vext) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/" << kvs << " = " << vm << " excludes "
                << "time2_vet_extrap = true" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // vet_col: the same modes, plus extrap (D^n + (dt/dt_prev)(D^n - D^{n-1}); its history
  // is not in the restart file) and extrap2
  t2_vcmode = 0;
  if (vet_col) {
    const bool vdef = sph_geom && !global_variable::restart_run;
    const std::string kvc = VetKeyName(pin, "time2_vet", {"time2_vet_col"});
    std::string vm = (sph_geom || pin->DoesParameterExist("rad_vet", kvc)) ?
        pin->GetOrAddString("rad_vet", kvc, vdef ? "extrap2" : "lag") :
        std::string("lag");
    if (vm.compare("lag") == 0) {
      t2_vcmode = 0;
    } else if (vm.compare("predict") == 0) {
      t2_vcmode = 1;
    } else if (vm.compare("rebuild") == 0) {
      t2_vcmode = 2;
    } else if (vm.compare("extrap") == 0) {
      t2_vcmode = 3;
    } else if (vm.compare("extrap2") == 0) {
      t2_vcmode = 4;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/" << kvc << " = '" << vm
                << "' is not a choice (lag | predict | rebuild | extrap | extrap2)"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (t2_vcmode != 0 && vcol_every != 1) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/" << kvc << " = " << vm
                << " needs vet_every (vet_col_every) = 1" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // time2_vet_source (alias time2_vet_col_source)
    Time2VetSourceInit(pin, t2_vcmode == 1 || t2_vcmode == 2, true);
  }
  // vet_sc reads the same time2_vet_source
  if (vet_sc) {
    Time2VetSourceInit(pin, t2_vsmode == 1 || t2_vsmode == 2, false);
  }
  // time2_pred_cap: the predicted state is P = U^n + theta dt K1 per cell,
  //   theta = min(1, cap/r),  r = max(|dt K1_E|/E^n, |eint(P) - eint_b|/eint_b),
  // eint_b the stage-start eint without the coupling's K1.  It stops the previous step's
  // slope from overshooting where the state relaxes within a step; theta = 1 (bitwise the
  // unlimited build) where r <= cap.  0 = off.
  t2_pcap = 0.0;
  if (t2_vsmode != 0 || t2_vcmode == 1 || t2_vcmode == 2) {
    t2_pcap = pin->GetOrAddReal("rad_vet", "time2_pred_cap",
                                global_variable::restart_run ? 0.0 : 0.5);
    if (!(t2_pcap >= 0.0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/time2_pred_cap must be >= 0" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // time2_lin_tol / time2_lin_tol_fac: Krylov tolerance of the stage solves
  // (time2_lin_tol, else implicit_lin_tol x time2_lin_tol_fac)
  const bool rs = global_variable::restart_run;
  // time2_one_pass_safety: implicit_one_pass safety factor of the stage solves (0 =
  // implicit_one_pass_safety); one-pass acceptance leaves a Picard error ~ tol/safety
  t2_onep_s = pin->GetOrAddReal("rad_vet", "time2_one_pass_safety", rs ? 0.0 : 30.0);
  if (t2_onep_s != 0.0 && !(t2_onep_s >= 1.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<rad_vet>/time2_one_pass_safety must be 0 or >= 1"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  t2_dbg_fail = -1;
  if (pin->DoesParameterExist("rad_vet", "time2_dbg_fail")) {
    t2_dbg_fail = pin->GetInteger("rad_vet", "time2_dbg_fail");
  }
  // time2_reject: be (backward-Euler redo) | retry (needs <time>/dt_retry = true; without
  // it the step falls back to backward Euler)
  t2_reject = 0;
  if (pin->DoesParameterExist("rad_vet", "time2_reject")) {
    const std::string tr = pin->GetString("rad_vet", "time2_reject");
    if (tr == "retry") {
      t2_reject = 1;
    } else if (tr != "be") {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/time2_reject = '" << tr << "' is not a choice "
                << "(be | retry)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // time2_k1_stiff = s: in cells with relaxation rate k dt > s the exchange part of K1
  // moves to E (Time2FormStage)
  t2_k1stiff = pin->GetOrAddReal("rad_vet", "time2_k1_stiff", 0.0);
  // time2_k1_redist: only 0 is accepted (old restart files may store 0.1)
  if (pin->DoesParameterExist("rad_vet", "time2_k1_redist") &&
      pin->GetReal("rad_vet", "time2_k1_redist") != 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<rad_vet>/time2_k1_redist is removed; only 0 is accepted. Restart "
              << "files of the old solver may store 0.1: name time2_k1_redist = 0"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // time2_face_fct: Zalesak limiter of the stage increment (Time2FaceFct); on for a fresh
  // run.  The clean solve stores the energy flux Phi that moved E (CleanWriteBack)
  t2_ffct = pin->GetOrAddBoolean("rad_vet", "time2_face_fct",
      recipe::Fix(pin, "rad_vet", "time2_face_fct", "true", "false"));
  t2_pfst = true;
  t2_ffct_emin = 0.1;   // (rad_vet.hpp)
  // time2_accept_fallback (read only when named)
  t2_accfb = true;
  if (pin->DoesParameterExist("rad_vet", "time2_accept_fallback")) {
    t2_accfb = pin->GetBoolean("rad_vet", "time2_accept_fallback");
  }
  // time2_adm_poslim (read only when named)
  if (pin->DoesParameterExist("rad_vet", "time2_adm_poslim")) {
    t2_adm_pl = pin->GetBoolean("rad_vet", "time2_adm_poslim");
  }
  t2_cfl_post = pin->GetOrAddBoolean("rad_vet", "time2_cfl_post", t2_cfl_post);
  t2_abe = pin->GetOrAddInteger("rad_vet", "time2_auto_be", t2_abe);
  if (t2_abe < 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<rad_vet>/time2_auto_be must be >= 0" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (global_variable::my_rank == 0 && t2_abe > 0) {
    std::cout << "<rad_vet> time2_auto_be=" << t2_abe << ": after a stage fallback the "
              << "next 1, 2, 4, ... " << t2_abe << " steps take the backward-Euler redo "
              << "path at once" << std::endl;
  }
  if (pin->DoesParameterExist("rad_vet", "restart_pack_per_mb")) {
    rst_pack_mb = pin->GetBoolean("rad_vet", "restart_pack_per_mb");
  }

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int nmb = std::max((pmy_pack->nmb_thispack), (pmy_pack->pmesh->nmb_maxperrank));
  int n1 = indcs.nx1 + 2*(indcs.ng);
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  Kokkos::realloc(t2k1, nmb, VET_T2_NK, n3, n2, n1);
  Kokkos::realloc(t2k2, nmb, VET_T2_NK, n3, n2, n1);
  Kokkos::realloc(t2inc, nmb, VET_T2_NK, n3, n2, n1);
  Kokkos::deep_copy(t2k1, 0.0);
  Kokkos::deep_copy(t2k2, 0.0);
  Kokkos::deep_copy(t2inc, 0.0);
  if (t2_pfst) {  // zero: no face flux known (start, restart)
    Kokkos::realloc(t2pf1, nmb, 3, n3, n2, n1);
    Kokkos::realloc(t2pf2, nmb, 3, n3, n2, n1);
    Kokkos::deep_copy(t2pf1, 0.0);
    Kokkos::deep_copy(t2pf2, 0.0);
    t2pf1_ok = false;
    t2pf2_ok = false;
  }
  Kokkos::realloc(t2f1, nmb, n3, n2, n1+1);
  if (trans_on) {
    Kokkos::realloc(t2f2, nmb, n3, n2+1, n1);
    if (trans_x3) {Kokkos::realloc(t2f3, nmb, n3+1, n2, n1);}
  }
  if (vet_sc) {
    Kokkos::realloc(vet_prev, nmb, VET_T2_NVET, n3, n2, n1);
    Kokkos::realloc(vet_now, nmb, VET_T2_NVET, n3, n2, n1);
    Kokkos::realloc(vet_opac, nmb, vscat ? (VET_NOPAC + 1) : VET_NOPAC, n3, n2,
                    n1);
    Kokkos::deep_copy(vet_prev, 0.0);
  }
  if (vet_col) {
    // time2_vet_col = extrap: the previous build's D and q
    if (t2_vcmode == 3) {
      Kokkos::realloc(vcol_prev, nmb, n3, n2, n1);
      Kokkos::realloc(vcol_qprev, nmb, n3, n2);
      t2_vcprev = false;
    }
  }
  if (t2_pcap > 0.0) {
    Kokkos::realloc(t2pth, nmb, n3, n2, n1);
    Kokkos::deep_copy(t2pth, 1.0);
  }
  t2_ok = false;
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_vet> time_scheme=hesdirk2 (g=" << kT2G << "; first step, restarts "
              << "without the slope and failed stages take backward Euler)" << std::endl;
    std::cout << "         time2_enth_be=" << (t2_enth_be ? "true" : "false")
              << " implicit_cs2_vimp_npass=" << cs2_vnp << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::Time2Active
//! \brief true when the next step runs the two stage solves

bool RadiationVet::Time2Active() {
  if (time_scheme != VET_TIME_HESDIRK2) return false;
  if (static_cast<int>(t2k1.extent(0)) < pmy_pack->nmb_thispack) {t2_ok = false;}
  return t2_ok;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::Time2GasTime
//! \brief time2_tableau = sdirk2: one time per stage for gas and radiation, which keeps
//! the IMEX coupling conditions b_E.c_I = b_I.c_E = 1/2 with the Heun hydro.
//!   which = 1 (before the stage-1 solve): the Heun-stage-1 gas U* is kept in t2hs; the
//!     solve takes S = U^n + g (U* - U^n), the gas at t^n + g dt.
//!   which = 2 (first in HydroConToPrim, gas Y = S + d): the gas becomes U* + d.
//!   which = 3 (after the stage tasks): w0 is refreshed from the hydro's stage-2 flux
//!     state U* + f d (f = 1/g non-stiff, -> 1 stiff), then u0 is U* + d again.
//! time2_gtime_mom: momentum's f = 1/g in every cell; IEN takes its kinetic-energy
//! change.  The flux state's eint increment (hs_ IEN, raw) is (f - 1) d's, or with
//! time2_gtime_slope (1-g) dt minmod(K1, K2), K1 the last step's stage-2 slope of eint
//! and K2 = d/(g dt) (f not used for the energy).  At which = 3 it is added to the gas
//! the caps and floors left, in the Patankar form where it is negative (positive).  Its
//! floors act on the primitives only: u0 is restored to U* + d as it was, so they
//! create no energy and nothing is paid.
//! time2_gtime_pos = false: the old flux state U* + f d,
//! formed at which = 2 on the gas before the caps, its floors kept in u0.

//----------------------------------------------------------------------------------------

void RadiationVet::Time2GasTime(Driver *pdrive, int which) {
  hydro::Hydro *ph = pmy_pack->phydro;
  if (t2_tab != 2 || !coupling || ph == nullptr || (which > 1 && t2_gt != which - 1))
    return;
  t2_gt = (which == 3) ? 0 : which;
  auto u0_ = ph->u0;
  auto u1_ = ph->u1;
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int nmb = u0_.extent_int(0), nv = u0_.extent_int(1);
  const int n3 = u0_.extent_int(2), n2 = u0_.extent_int(3), n1 = u0_.extent_int(4);
  const bool all = (which == 1);
  const int il = all ? 0 : id.is, iu = all ? n1-1 : id.ie, jl = all ? 0 : id.js;
  const int ju = all ? n2-1 : id.je, kl = all ? 0 : id.ks, ku = all ? n3-1 : id.ke;
  if (t2hs.extent_int(0) != nmb || t2hs.extent_int(1) != nv) {
    Kokkos::realloc(t2hs, nmb, nv, n3, n2, n1);
  }
  auto hs_ = t2hs;
  // flux-state factor on d: 1/g + O(k dt), -> 1 where the exchange is stiff, k = ch rho
  // kappa_E + 4 cl rho kappa_P a T^3/(rho c_v)
  auto op_ = opac;
  auto eos = ph->peos->eos_data;
  auto phc = ph->phicc0;
  auto ccl = pmy_pack->pcoord->cos_cell;
  const bool etg = ph->use_etotgrav, csk = cs_geom;
  const Real g = kT2G, cl = c_light, ch = chat, ar = arad, dt = pmy_pack->pmesh->dt;
  const bool gm = t2_gtm, gs = t2_gts, ps = t2_gtp;
  auto k1_ = t2k1;
  par_for("vet_t2_gtime", DevExeSpace(), 0, nmb-1, kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real e = 0.0, xs = 0.0;
    if (which == 2) {
      const Real dd = u0_(m,IDN,k,j,i);
      const Real eg = u0_(m,IEN,k,j,i) - CsKeMom(u0_(m,IM1,k,j,i), u0_(m,IM2,k,j,i),
                      u0_(m,IM3,k,j,i), csk ? ccl(m,k,j) : 0.0)/dd -
                      (etg ? dd*phc(m,k,j,i) : 0.0);
      const Real tk = (eg > 0.0) ? eos.Temperature(dd, eg) : 0.0;
      const Real z = dt*(ch*op_(m,VET_OP_E,k,j,i) + ((eg > 0.0) ?
                     4.0*cl*op_(m,VET_OP_P,k,j,i)*ar*SQR(SQR(tk))/eg : 0.0));
      e = (1.0 + g*z)/(g*(1.0 + z)) - 1.0;
      if (gs) {
        const Real cc = csk ? ccl(m,k,j) : 0.0, r = u1_(m,IDN,k,j,i);
        Real p[3], q[3], s[3];   // momenta of U^n, U^n + g dt K1 and S
        for (int n = 0; n < 3; ++n) {
          p[n] = u1_(m,IM1+n,k,j,i);
          q[n] = p[n] + g*dt*k1_(m,VET_T2_M1+n,k,j,i);
          s[n] = p[n] + g*(hs_(m,IM1+n,k,j,i) - p[n]);
        }
        // the eint changes g dt K1 and g dt K2 = d
        const Real a = g*dt*k1_(m,VET_T2_EN,k,j,i) -
                       (CsKeMom(q[0], q[1], q[2], cc) - CsKeMom(p[0], p[1], p[2], cc))/r;
        const Real es = u1_(m,IEN,k,j,i) + g*(hs_(m,IEN,k,j,i) - u1_(m,IEN,k,j,i)) -
                        CsKeMom(s[0], s[1], s[2], cc)/dd -
                        (etg ? dd*phc(m,k,j,i) : 0.0);
        const Real b = eg - es;
        xs = (a*b > 0.0) ? (1.0/g - 1.0)*((fabs(a) < fabs(b)) ? a : b) : 0.0;
        e = 0.0;
      }
    }
    if (which == 3 && ps) {
      const Real c = csk ? ccl(m,k,j) : 0.0, r = u0_(m,IDN,k,j,i);
      const Real k0 = CsKeMom(u0_(m,IM1,k,j,i), u0_(m,IM2,k,j,i), u0_(m,IM3,k,j,i), c)/r;
      const Real k1 = CsKeMom(u0_(m,IM1,k,j,i) + hs_(m,IM1,k,j,i), u0_(m,IM2,k,j,i) +
                              hs_(m,IM2,k,j,i), u0_(m,IM3,k,j,i) + hs_(m,IM3,k,j,i), c)/r;
      const Real eb = u0_(m,IEN,k,j,i) - k0 - (etg ? r*phc(m,k,j,i) : 0.0);
      xs = hs_(m,IEN,k,j,i);
      hs_(m,IEN,k,j,i) = k1 - k0 + ((eb > 0.0) ? ((xs < 0.0) ? xs*eb/(eb - xs) : xs) :
                                    0.0);
    }
    for (int n = 0; n < nv; ++n) {
      const Real un = u1_(m,n,k,j,i), us = hs_(m,n,k,j,i), u = u0_(m,n,k,j,i);
      if (which == 1) {
        hs_(m,n,k,j,i) = u;
        u0_(m,n,k,j,i) = un + g*(u - un);
      } else if (which == 2) {
        const Real d = u - (un + g*(us - un));   // (0 where the solve left it)
        u0_(m,n,k,j,i) = us + d;
        hs_(m,n,k,j,i) = ((gm && n >= IM1 && n <= IM3) ? 1.0/g - 1.0 : e)*d;
      } else {
        u0_(m,n,k,j,i) = u + us;
        if (ps) {hs_(m,n,k,j,i) = u;}   // (the register to restore)
      }
    }
    if (which == 2 && gs) {
      hs_(m,IEN,k,j,i) = xs;
    } else if (which == 2 && (ps || gm)) {
      // ps: (f - 1) d of the total energy less that of the kinetic; else IEN's with the
      // momentum's f
      const Real c = csk ? ccl(m,k,j) : 0.0, r = u0_(m,IDN,k,j,i) + hs_(m,IDN,k,j,i);
      const Real s = gm ? e*g/(1.0 - g) : 1.0;
      auto ke = [&](const Real q) {
        return CsKeMom(u0_(m,IM1,k,j,i) + q*hs_(m,IM1,k,j,i), u0_(m,IM2,k,j,i) +
                       q*hs_(m,IM2,k,j,i), u0_(m,IM3,k,j,i) + q*hs_(m,IM3,k,j,i), c)/r;
      };
      hs_(m,IEN,k,j,i) += ps ? ke(0.0) - ke(s) : ke(1.0) - ke(s);
    }
  });
  if (which != 3) return;
  const bool dpay = ph->dual.pay;
  if (ps) {ph->dual.pay = ph->dual.book = false;}
  FluidRefresh(pdrive, false);   // ghosts + w0 of the flux state
  ph->dual.pay = dpay;
  ph->dual.book = true;
  // u0 back to U* + d, in the ghosts the neighbour's (exchanged; the flux state at the
  // physical boundaries): a ghost is then its neighbour's cell bitwise, so the stage-2
  // dfloor_flux_limit takes the same donor state on both sides of a MeshBlock face (the
  // face flux does not depend on the MeshBlock layout)
  par_for("vet_t2_gtime3g", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (k >= kl && k <= ku && j >= jl && j <= ju && i >= il && i <= iu) return;
    for (int n = 0; n < nv; ++n) {hs_(m,n,k,j,i) = ps ? u0_(m,n,k,j,i) : 0.0;}
  });
  (void) ph->pbval_u->InitRecv(nv);
  (void) ph->pbval_u->PackAndSendCC(hs_, ph->coarse_u0);
  (void) ph->pbval_u->ClearSend();
  (void) ph->pbval_u->ClearRecv();
  (void) ph->pbval_u->RecvAndUnpackCC(hs_, ph->coarse_u0);
  par_for("vet_t2_gtime3", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int n = 0; n < nv; ++n) {
      u0_(m,n,k,j,i) = ps ? hs_(m,n,k,j,i) : u0_(m,n,k,j,i) - hs_(m,n,k,j,i);
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::Time2FormStage
//! \brief after the hydro's stage `stage`: the radiation stage start (stage 1 U^n, stage
//! 2 the Heun average of U^n and Y2) and t2inc, which ImplicitSolve adds to form the old
//! vector of the stage solve

void RadiationVet::Time2FormStage(int stage, Real dt) {
  const Real g = kT2G;
  const Real w0 = t2_w0, w1 = t2_w1;
  // hesdirk2 and sdirk2 take the Heun stage start; trbdf2 the w0/w1 average
  const bool heun = (t2_tab == 0 || t2_tab == 2);
  Real c1, c2;
  auto u0_ = u0;
  auto u1_ = u1;
  const int nf1 = static_cast<int>(f0x1.extent(3));
  if (stage == 1) {
    Kokkos::deep_copy(DevExeSpace(), u1, u0);
    Kokkos::deep_copy(DevExeSpace(), t2f1, f0x1);
    if (trans_on) {
      Kokkos::deep_copy(DevExeSpace(), t2f2, f0x2);
      if (trans_x3) {Kokkos::deep_copy(DevExeSpace(), t2f3, f0x3);}
    }
    c1 = t2_a21*dt;
    c2 = 0.0;
    // time2_k1_stiff: where k dt > s (k as in Time2GasTime) the exchange part of K1 moves
    // from the gas energy to E (ch/cl weighted, total energy kept), so the L-stable
    // stages alone take the stiff exchange; K1_EN keeps the kinetic-energy change
    if (t2_k1stiff > 0.0) {
      FluidRef fls = FluidRef::Get(pmy_pack);
      if (fls.on) {
        auto &id = pmy_pack->pmesh->mb_indcs;
        const int is = id.is, ie = id.ie, js = id.js, je = id.je, ks = id.ks, ke = id.ke;
        const int nmbr = pmy_pack->nmb_thispack - 1;
        auto ug = fls.u0;
        auto eos = fls.eos;
        const bool etg = fls.use_etotgrav;
        auto phc = etg ? fls.phicc0 : DvceArray4D<Real>();
        const bool emh = fls.mhd;
        auto emg = emag0;
        auto op_ = opac;
        auto k1s_ = t2k1;
        const bool csk = cs_geom;
        auto ccl = pmy_pack->pcoord->cos_cell;
        const Real cl = c_light, ch = chat, ar = arad, sth = t2_k1stiff, dtl = dt;
        const Real ca = c1;
        Real ns = 0.0;
        Kokkos::parallel_reduce("vet_t2_k1stiff", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
            DevExeSpace(), {0, ks, js, is}, {nmbr + 1, ke + 1, je + 1, ie + 1}),
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &sum) {
          const Real dd = ug(m,IDN,k,j,i), idd = 1.0/fmax(dd, 1.0e-300);
          const Real m1 = ug(m,IM1,k,j,i), m2 = ug(m,IM2,k,j,i), m3 = ug(m,IM3,k,j,i);
          const Real cc = csk ? ccl(m,k,j) : 0.0;
          const Real ek0 = CsKeMom(m1, m2, m3, cc)*idd;
          Real egr = etg ? dd*phc(m,k,j,i) : 0.0;
          if (emh) egr += emg(m,k,j,i);
          const Real eg = ug(m,IEN,k,j,i) - ek0 - egr;
          if (!(eg > 0.0)) return;
          const Real tk = eos.Temperature(dd, eg);
          if (!(tk > 0.0)) return;
          const Real rkp = op_(m,VET_OP_P,k,j,i), rke = op_(m,VET_OP_E,k,j,i);
          const Real rate = ch*rke + 4.0*cl*rkp*ar*tk*tk*tk*tk/eg;
          if (!(rate*dtl > sth)) return;
          const Real n1 = m1 + ca*k1s_(m,VET_T2_M1,k,j,i);
          const Real n2 = m2 + ca*k1s_(m,VET_T2_M1+1,k,j,i);
          const Real n3 = m3 + ca*k1s_(m,VET_T2_M1+2,k,j,i);
          const Real ek1 = CsKeMom(n1, n2, n3, cc)*idd;
          const Real kw = (ca != 0.0) ? (ek1 - ek0)/ca : 0.0;
          const Real dex = k1s_(m,VET_T2_EN,k,j,i) - kw;
          k1s_(m,VET_T2_EN,k,j,i) = kw;
          k1s_(m,VET_T2_E,k,j,i) += (ch/cl)*dex;
          sum += 1.0;
        }, Kokkos::Sum<Real>(ns));
        t2_k1stiff_n += ns;
      }
    }
  } else {
    int nmb1 = static_cast<int>(u0.extent(0)) - 1;
    int n3 = static_cast<int>(u0.extent(2)), n2 = static_cast<int>(u0.extent(3));
    int n1 = static_cast<int>(u0.extent(4));
    par_for("vet_t2_avg", DevExeSpace(), 0, nmb1, 0, VET_NVAR-1, 0, n3-1, 0, n2-1,
            0, n1-1,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
      u0_(m,n,k,j,i) = heun ? 0.5*(u1_(m,n,k,j,i) + u0_(m,n,k,j,i))
                            : (w1*u1_(m,n,k,j,i) + w0*u0_(m,n,k,j,i));
    });
    auto avgf = [&](DvceArray4D<Real> &a, DvceArray4D<Real> &b) {
      auto a_ = a;
      auto b_ = b;
      par_for("vet_t2_avgf", DevExeSpace(), 0, static_cast<int>(a.extent(0))-1,
              0, static_cast<int>(a.extent(1))-1, 0, static_cast<int>(a.extent(2))-1,
              0, static_cast<int>(a.extent(3))-1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        a_(m,k,j,i) = heun ? 0.5*(b_(m,k,j,i) + a_(m,k,j,i))
                           : (w1*b_(m,k,j,i) + w0*a_(m,k,j,i));
      });
    };
    avgf(f0x1, t2f1);
    if (trans_on) {
      avgf(f0x2, t2f2);
      if (trans_x3) {avgf(f0x3, t2f3);}
    }
    c1 = t2_c31*dt;
    c2 = t2_c32*dt;
  }
  (void) nf1;
  auto k1_ = t2k1;
  auto k2_ = t2k2;
  auto in_ = t2inc;
  int nmb1 = static_cast<int>(t2inc.extent(0)) - 1;
  int n3 = static_cast<int>(t2inc.extent(2)), n2 = static_cast<int>(t2inc.extent(3));
  int n1 = static_cast<int>(t2inc.extent(4));
  const bool two = (stage != 1);
  par_for("vet_t2_inc", DevExeSpace(), 0, nmb1, 0, VET_T2_NK-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
    Real v = c1*k1_(m,n,k,j,i);
    if (two) {v += c2*k2_(m,n,k,j,i);}
    in_(m,n,k,j,i) = v;
  });
  cln_lzc = c1*cln_lzs1 + (two ? c2*cln_lzs2 : 0.0);   // (clean_cs_lz_abs = 2)
  CleanDInc(two, c1, c2, heun ? 0.5 : w1, heun ? 0.5 : w0);
  // am_vector: the spin reservoir's increment from its slopes
  if (am_vec && aminc.extent_int(0) > 0) {
    auto a1_ = amk1;
    auto a2_ = amk2;
    auto ai_ = aminc;
    par_for("vet_t2_aminc", DevExeSpace(), 0, static_cast<int>(aminc.extent(0)) - 1, 0,
            static_cast<int>(aminc.extent(1)) - 1, 0,
            static_cast<int>(aminc.extent(2)) - 1,
            0, static_cast<int>(aminc.extent(3)) - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real v = c1*a1_(m,k,j,i);
      if (two) {v += c2*a2_(m,k,j,i);}
      ai_(m,k,j,i) = v;
    });
  }
  if (t2_ffct) {Time2FaceFct(stage, c1, c2);}
  CleanBndAcc(stage, c1, c2, (stage == 1) ? 0.0 : (heun ? 0.5 : w0), dt);
  t2_solve = (stage == 1) ? VET_T2S_STAGE1 : VET_T2S_STAGE2;
  t2_fail = false;
  nsub = 1;
  dt_sub = g*dt;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::Time2Stash
//! \brief stage-loop start (before Heun stage 1): checkpoints the clean energy-flux
//! memory and ray history

void RadiationVet::Time2Stash() {
  CleanDMem(true);
  CleanRayHistory(true);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanRayHistory
//! \brief save: keep the rays' history (CleanRayKeep); else restore it (a stage failed)

void RadiationVet::CleanRayHistory(const bool save) {
  DvceArray5D<Real> *a[6] = {&cln_fvi, &cln_fh, &cln_top, &cln_rr, &cln_rpd, &cln_rnew};
  auto &k = cln_rk;
  auto cp = [](auto &d, const auto &s) {
    if (s.size() == 0) {
      d = std::remove_reference_t<decltype(d)>();
      return;
    }
    bool same = d.size() == s.size();
    for (unsigned r = 0; r < s.rank(); ++r) same = same && d.extent(r) == s.extent(r);
    if (!same) {
      d = std::remove_reference_t<decltype(d)>(
          Kokkos::view_alloc(std::string("cln_rk"), Kokkos::WithoutInitializing),
          s.layout());
    }
    Kokkos::deep_copy(d, s);
  };
  for (int q = 0; q < 6; ++q) save ? cp(k.v[q], *a[q]) : cp(*a[q], k.v[q]);
  save ? cp(k.fvd, cln_fvd) : cp(cln_fvd, k.fvd);
  save ? k.fvdm.resize(cln_fvdm.size()) : cln_fvdm.resize(k.fvdm.size());
  for (std::size_t q = 0; q < k.fvdm.size(); ++q) {
    save ? cp(k.fvdm[q], cln_fvdm[q]) : cp(cln_fvdm[q], k.fvdm[q]);
  }
  auto sw = [&](auto &x, auto &y) {if (save) {y = x;} else {x = y;}};
  sw(cln_fvblk, k.fvblk); sw(cln_rb, k.rb); sw(cln_rpb, k.rpb); sw(cln_rbld, k.rbld);
  sw(cln_rpt, k.rpt); sw(cln_rs, k.rs); sw(cln_rguard_sum, k.rg); sw(cln_roff_sum, k.ro);
  sw(cln_rns, k.rns); sw(cln_rok, k.rok); sw(cln_rfix, k.rfix); sw(cln_rlast, k.rlast);
  sw(cln_rjump, k.rjump); sw(cln_fvn, k.fvn); sw(cln_fvkey, k.fvkey);
  sw(cln_rcut, k.rcut); sw(cln_rcut_n, k.rcut_n); sw(cln_rfdone, k.rfdone);
  sw(cln_rhist, k.rhist);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::Time2Restore
//! \brief a stage failed: put U^n back (radiation, faces, gas, B, spin) and refresh the
//! hydro; the Driver reruns the Heun stages and takes the step with backward Euler

void RadiationVet::Time2Restore(Driver *pdrive) {
  t2_gt = 0;   // (Time2GasTime)
  Kokkos::deep_copy(DevExeSpace(), u0, u1);
  Kokkos::deep_copy(DevExeSpace(), f0x1, t2f1);
  if (trans_on) {
    Kokkos::deep_copy(DevExeSpace(), f0x2, t2f2);
    if (trans_x3) {Kokkos::deep_copy(DevExeSpace(), f0x3, t2f3);}
  }
  FluidRef fl = FluidRef::Get(pmy_pack);
  if (fl.on) {
    Kokkos::deep_copy(DevExeSpace(), fl.u0, fl.u1);
    // mhd: the field back to U^n (b1 is the stage-1 copy of b0 incl. ghosts)
    if (fl.mhd) {
      mhd::MHD *pm = pmy_pack->pmhd;
      Kokkos::deep_copy(DevExeSpace(), pm->b0.x1f, pm->b1.x1f);
      Kokkos::deep_copy(DevExeSpace(), pm->b0.x2f, pm->b1.x2f);
      Kokkos::deep_copy(DevExeSpace(), pm->b0.x3f, pm->b1.x3f);
    }
    FluidRefresh(pdrive);
    if (fl.mhd) EmagBuild(false);
  }
  // am_vector: the spin reservoir back to U^n (am_spin1)
  if (am_vec) {
    hydro::AmModule *pam = hydro::AmModuleOf(pmy_pack);
    Kokkos::deep_copy(DevExeSpace(), pam->am_spin, pam->am_spin1);
  }
  CleanDMem(false);
  CleanRayHistory(false);
  cln_lzp = 0.0;   // (the L_z ledger: the failed stages' pending part, CleanLzBook)
  t2_fail = false;
  t2_ok = false;
  t2_nfall += 1.0;
  // the backward-Euler redo is reported, never silent
  if (global_variable::my_rank == 0) {
    std::printf("CLEAN WARNING cycle %d t2s %d: stage solve failed, the step is redone "
                "with backward Euler (stage fallback %.0f)\n", pmy_pack->pmesh->ncycle,
                t2_solve, t2_nfall);
    if (cln_flrun > 0.0 && cln_flrid == 4*pmy_pack->pmesh->ncycle + t2_solve) {
      std::printf(
          "CLEAN WARNING cycle %d t2s %d: the implicit_pos_floor creation %.6e of "
                  "the failed attempt is discarded with it\n", pmy_pack->pmesh->ncycle,
                  t2_solve, cln_flrun);
    }
  }
  cln_flrun = 0.0;
  // the predictor increments were stored by the failed solve: the redo starts cold from
  // U^n (extrapolating them again can grow without bound)
  pred_ok = false;
  pred2_ok = false;
  t2_solve = VET_T2S_NONE;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::Time2FctHalo
//! \brief exchange of the limiter R (iw slot c): first ghost layer only, the clean direct
//! lists, without the deferred clean exchanges

void RadiationVet::Time2FctHalo(const int c) {
  std::vector<CleanGSeg<Real>> hp;
  hp.swap(cln_hpend);
  CleanGhostFill(iw, c, 1, false, 0, 1);
  cln_hpend.swap(hp);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::Time2FaceFct
//! \brief time2_face_fct: Zalesak positivity limiter on the face fluxes of the stage
//! increment t2inc (E and face F channels).  G_f = A_f (c1 t2pf1 + c2 t2pf2) is the
//! energy the increment carries across face f; only cells whose old vector E + dE is not
//! > 0 are limited, by R_i = (1 - 1e-12) Q_i/P_i, and a face takes its donor's R.  Rounds
//! repeat until no cell is newly limited (R only decreases; cap 4096).  Bitwise where no
//! face is limited.

void RadiationVet::Time2FaceFct(int stage, Real c1, Real c2) {
  const bool two = (stage != 1);
  auto rhalo = [&](const int c) {Time2FctHalo(c);};
  // no stored face flux for K1 (first step, restart): K1 is the last stage-2 slope, whose
  // solution is U^n, so its face flux is the stage-start t2f* = f0x*
  if (!two && !t2pf1_ok && t2pf1.extent_int(0) >= pmy_pack->nmb_thispack &&
      t2f1.extent_int(0) >= pmy_pack->nmb_thispack) {
    auto pf_ = t2pf1;
    auto g1_ = t2f1, g2_ = t2f2, g3_ = t2f3;
    const bool q2 = trans_on && t2f2.extent_int(0) > 0;
    const bool q3 = trans_on && trans_x3 && t2f3.extent_int(0) > 0;
    const Real crp = chat/c_light;
    par_for("vet_t2_ffct_fill", DevExeSpace(), 0, pf_.extent_int(0) - 1,
            0, pf_.extent_int(2) - 1, 0, pf_.extent_int(3) - 1, 0, pf_.extent_int(4) - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      pf_(m,0,k,j,i) = crp*g1_(m,k,j,i);
      pf_(m,1,k,j,i) = q2 ? crp*g2_(m,k,j,i) : 0.0;
      pf_(m,2,k,j,i) = q3 ? crp*g3_(m,k,j,i) : 0.0;
    });
    t2pf1_ok = true;
    t2_ffct_fill += 1.0;
  }
  if (!t2pf1_ok || (two && !t2pf2_ok) || t2pf1.extent_int(0) < pmy_pack->nmb_thispack) {
    t2_ffct_nop += 1.0;  // no stored face flux for this slope
    return;
  }
  if (iw.extent_int(1) <= Time2FctSlot(0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<rad_vet>/time2_face_fct needs the Krylov work slots "
              << "(implicit_solver = bicgstab)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  t2_ffct_calls += 1.0;
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int is = id.is, ie = id.ie, js = id.js, je = id.je, ks = id.ks, ke = id.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto in_ = t2inc;
  auto u0_ = u0;
  auto p1_ = t2pf1;
  auto p2_ = t2pf2;
  auto iw_ = iw;
  const int ca = Time2FctSlot(0), cb = Time2FctSlot(1);
  const bool tr2 = trans_on, tr3 = trans_on && trans_x3;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto vol = pmy_pack->pcoord->volume;
  auto ar1 = pmy_pack->pcoord->area.x1f;
  auto ar2 = pmy_pack->pcoord->area.x2f;
  auto ar3 = pmy_pack->pcoord->area.x3f;
  const bool curv = (cs_geom || sph_geom) && vol.extent_int(0) > nmb1 &&
                    vol.extent_int(3) > 1 && ar1.extent_int(0) > nmb1;
  const bool hoi = impl_hold && (hold_m.extent_int(0) > nmb1);
  auto hm_ = hold_m;
  const Real cc1 = c1, cc2 = two ? c2 : 0.0;
  // energy across the lower face of cell (k,j,i) in direction d (+ = towards +x_d)
  auto gf = KOKKOS_LAMBDA(const int m, const int d, const int k, const int j,
                          const int i) {
    const Real phi = cc1*p1_(m,d,k,j,i) + cc2*p2_(m,d,k,j,i);
    Real a;
    if (curv) {
      a = (d == 0) ? ar1(m,k,j,i) : ((d == 1) ? ar2(m,k,j,i) : ar3(m,k,j,i));
    } else {
      const auto &s = mbsize.d_view(m);
      a = (d == 0) ? s.dx2*s.dx3 : ((d == 1) ? s.dx1*s.dx3 : s.dx1*s.dx2);
    }
    return a*phi;
  };
  auto cvol = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const auto &s = mbsize.d_view(m);
    return curv ? vol(m,k,j,i) : s.dx1*s.dx2*s.dx3;
  };
  auto un_ = u1;   // (E^n: the limiter keeps th V E^n in a limited cell)
  const Real eth = t2_ffct_emin;
  // w = V (E + dE) unlimited, outflow po, inflow pi, the inflow pl the donors' R let
  // through, the change dw of V E_new by R (slot sa); lim: an adjacent face is limited
  auto cell = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i,
                            Real &w, Real &po, Real &pi, Real &pl, Real &dw,
                            bool &lim, const int sa) {
    const Real rc = fabs(iw_(m,sa,k,j,i));
    po = 0.0;
    pi = 0.0;
    pl = 0.0;
    dw = 0.0;
    lim = false;
    for (int d = 0; d < 3; ++d) {
      if ((d == 1 && !tr2) || (d == 2 && !tr3)) {continue;}
      const int dk = (d == 2) ? 1 : 0, dj = (d == 1) ? 1 : 0, di = (d == 0) ? 1 : 0;
      // lower face: + = inflow, donor the lower neighbour
      Real g = gf(m, d, k, j, i);
      if (g > 0.0) {
        pi += g;
        const Real rn = fabs(iw_(m,sa,k-dk,j-dj,i-di));
        pl += rn*g;
        if (rn < 1.0) {dw -= (1.0 - rn)*g; lim = true;}
      } else if (g < 0.0) {
        po -= g;
        if (rc < 1.0) {dw -= (1.0 - rc)*g; lim = true;}
      }
      // upper face: + = outflow, donor this cell
      g = gf(m, d, k+dk, j+dj, i+di);
      if (g > 0.0) {
        po += g;
        if (rc < 1.0) {dw += (1.0 - rc)*g; lim = true;}
      } else if (g < 0.0) {
        pi -= g;
        const Real rn = fabs(iw_(m,sa,k+dk,j+dj,i+di));
        pl -= rn*g;
        if (rn < 1.0) {dw += (1.0 - rn)*g; lim = true;}
      }
    }
    w = cvol(m,k,j,i)*(u0_(m,VET_E,k,j,i) + in_(m,VET_T2_E,k,j,i));
  };
  // R = 1 everywhere incl. ghosts (physical ghosts stay unlimited); rounds read slot sa
  // and write sb, then swap (both slots keep the same ghosts)
  for (const int c : {ca, cb}) {
    Kokkos::deep_copy(DevExeSpace(), Kokkos::subview(iw, Kokkos::ALL, c, Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL), 1.0);
  }
  int sa = ca, sb = cb;
  Real ncell = 0.0;
  int nround = 0;
  // zalesak_halo_skip: exchange R only after a round that changed a cell within 2 ng of a
  // MeshBlock face; the first exchange always runs
  const int zd = 2*id.ng;
  const bool zx2 = pmy_pack->pmesh->multi_d, zx3 = pmy_pack->pmesh->three_d;
  const bool zsk = zl_hskip;
  bool zfirst = true;
  // zalesak_active: an R that fell is written with its sign bit set; the next round
  // evaluates only cells with such a face neighbour, the others keep |R|
  const bool zac = zl_active;  // zalesak_active
  // zalesak_active, no exchange: the face ghosts of the slot read next are those of the
  // other slot
  auto fghost = [&](const int c, const int cf) {
    const int n1 = ie - is + 1, n2 = je - js + 1, n3 = ke - ks + 1;
    const int f1 = n3*n2, f2 = tr2 ? n3*n1 : 0, f3 = tr3 ? n2*n1 : 0;
    const int nf = 2*(f1 + f2 + f3);
    Kokkos::parallel_for("vet_t2_ffct_fg", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                         (nmb1 + 1)*nf),
    KOKKOS_LAMBDA(const int t) {
      const int m = t/nf;
      int u = t - m*nf, k, j, i;
      if (u < 2*f1) {
        i = (u < f1) ? is - 1 : ie + 1;
        u %= f1;
        k = ks + u/n2;
        j = js + u % n2;
      } else if (u < 2*(f1 + f2)) {
        u -= 2*f1;
        j = (u < f2) ? js - 1 : je + 1;
        u %= f2;
        k = ks + u/n1;
        i = is + u % n1;
      } else {
        u -= 2*(f1 + f2);
        k = (u < f3) ? ks - 1 : ke + 1;
        u %= f3;
        j = js + u/n1;
        i = is + u % n1;
      }
      iw_(m,c,k,j,i) = fabs(iw_(m,cf,k,j,i));
    });
  };
  for (int r = 0; r < 4096; ++r) {   // until none newly limited
    Real nv[2] = {0.0, 0.0};
    // (zalesak_active: every cell in round 0)
    const bool zsel = zac && r > 0;
    const int s0 = sa, s1 = sb;
    auto rnd = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i,
                             Real &a, Real &ab) {
      const Real rc = fabs(iw_(m,s0,k,j,i));
      Real rn = rc;
      bool act = !zsel;   // (zalesak_active: a face neighbour's R fell)
      for (int d = 0; d < 3 && !act; ++d) {
        if ((d == 1 && !tr2) || (d == 2 && !tr3)) {continue;}
        const int dk = (d == 2), dj = (d == 1), di = (d == 0);
        act = Kokkos::signbit(iw_(m,s0,k-dk,j-dj,i-di)) ||
              Kokkos::signbit(iw_(m,s0,k+dk,j+dj,i+di));
      }
      if (act && !(hoi && hm_(m,k,j,i) > 0.0)) {
        Real w, po, pi, pl, dw;
        bool lim;
        cell(m, k, j, i, w, po, pi, pl, dw, lim, s0);
        const Real wf = eth*cvol(m,k,j,i)*un_(m,VET_E,k,j,i);
        if (!(w + dw > wf)) {
          // held without face transport + remaining inflow, above th V E^n
          const Real q = w + po - pi + pl - wf;
          const Real rq = (q > 0.0 && po > 0.0) ? fmin(1.0, (1.0 - 1.0e-12)*q/po) : 0.0;
          if (rq < rc) {
            rn = zac ? -rq : rq;
            a += 1.0;
            if (zsk && (i < is + zd || i > ie - zd ||
                        (zx2 && (j < js + zd || j > je - zd)) ||
                        (zx3 && (k < ks + zd || k > ke - zd)))) {ab += 1.0;}
          }
        }
      }
      iw_(m,s1,k,j,i) = rn;
    };
    par_reduce_flat2("vet_t2_ffct_r", 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &a,
                  Real &ab) {
      rnd(m, k, j, i, a, ab);
    }, Kokkos::Sum<Real>(nv[0]), Kokkos::Sum<Real>(nv[1]));
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, nv, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
    const Real nn = nv[0];
    if (!(nn > 0.0)) {break;}
    ncell += nn;
    nround = r + 1;
    std::swap(sa, sb);   // (the round's R is read next)
    if (!zsk || zfirst || nv[1] > 0.0) {
      rhalo(sa);
      if (!zac) {   // (the other slot keeps the same ghosts)
        Kokkos::deep_copy(DevExeSpace(), Kokkos::subview(iw, Kokkos::ALL, sb,
                          Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), Kokkos::subview(iw,
                          Kokkos::ALL, sa, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
      }
      zl_halos += 1.0;
      zfirst = false;
    } else if (zac) {
      fghost(sa, sb);
    }
  }
  if (zac) {   // (the last round changed nothing: |R| of its read slot, ghosts too)
    auto r0 = Kokkos::subview(iw, Kokkos::ALL, sa, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
    auto r1 = Kokkos::subview(iw, Kokkos::ALL, ca, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
    par_for("vet_t2_ffct_abs", DevExeSpace(), 0, nmb1, 0, r0.extent_int(1) - 1, 0,
            r0.extent_int(2) - 1, 0, r0.extent_int(3) - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      r1(m,k,j,i) = fabs(r0(m,k,j,i));
    });
  } else if (sa != ca) {   // (the last round changed nothing: the R of both slots)
    Kokkos::deep_copy(DevExeSpace(), Kokkos::subview(iw, Kokkos::ALL, ca, Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL), Kokkos::subview(iw, Kokkos::ALL, sa,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (nround >= 4096) {
    ZalesakCapWarn("time2_face_fct", pmy_pack->pmesh->ncycle, t2_solve);
  }
  zl_rounds += static_cast<Real>(nround);
  t2_ffct_cells += ncell;
  t2_ffct_rmax = std::max(t2_ffct_rmax, static_cast<Real>(nround));
  if (nround == 0) {return;}   // nothing limited: t2inc untouched (bitwise)
  // apply: the E channel of every cell next to a limited face, then the face F channels
  Real lv[3] = {0.0, 0.0, 0.0}, qmin = 1.0e300;
  // a flat range (faster than a 4-D MDRange reduction)
  const int n1 = ie - is + 1, nji = (je - js + 1)*n1, nkji = (ke - ks + 1)*nji;
  Kokkos::parallel_reduce("vet_t2_ffct_e", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int idx, Real &sl, Real &sw, Real &se, Real &sq) {
    const int m = idx/nkji, k = ks + (idx%nkji)/nji, j = js + (idx%nji)/n1;
    const int i = is + idx%n1;
    Real w, po, pi, pl, dw;
    bool lim;
    cell(m, k, j, i, w, po, pi, pl, dw, lim, ca);
    const Real rc = iw_(m,ca,k,j,i);
    if (lim) {in_(m,VET_T2_E,k,j,i) += dw/cvol(m,k,j,i);}
    if (rc < 1.0) {sw += (1.0 - rc)*po;}
    if (!(hoi && hm_(m,k,j,i) > 0.0) && !(w + dw > 0.0)) {sl += 1.0;}
    if (lim) {   // (the limited old vector against E^n)
      const Real q = (u0_(m,VET_E,k,j,i) + in_(m,VET_T2_E,k,j,i))/un_(m,VET_E,k,j,i);
      if (q < eth) {se += 1.0;}
      sq = fmin(sq, q);
    }
  }, Kokkos::Sum<Real>(lv[0]), Kokkos::Sum<Real>(lv[1]), Kokkos::Sum<Real>(lv[2]),
     Kokkos::Min<Real>(qmin));
  // clean staggered: face F is the radiation momentum, left alone by the limiter, so
  // momentum and L_z stay exact
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, lv, 3, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &qmin, 1, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
#endif
  t2_ffct_qmin = std::min(t2_ffct_qmin, qmin);
  t2_ffct_nemp += lv[2];
  t2_ffct_left += lv[0];
  t2_ffct_de += lv[1];
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanBndAcc
//! \brief the history's x1 boundary energy cln_lacc, booking the flux that moved E: stage
//! 0 starts a backward-Euler step (its solves add h A Phi); stage 1 / 2 the old vector of
//! a hesdirk2 stage, w x (energy of Y2) + A (c1 t2pf1 + c2 t2pf2) at the donors' limiter
//! R; the stage solve adds g dt A Phi

void RadiationVet::CleanBndAcc(int stage, Real c1, Real c2, Real w, Real dt) {
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  if (!(sph_geom || cs_geom)) return;
  if (cln_lacc.extent_int(0) <= nmb1) {
    Kokkos::realloc(cln_lacc, nmb1 + 1, 4, u0.extent_int(2), u0.extent_int(3));
  }
  if (stage == 0 || stage == 1) {cln_ldt = dt;}
  if (stage == 0) {Kokkos::deep_copy(DevExeSpace(), cln_lacc, 0.0); return;}
  const bool two = (stage != 1), ok1 = t2pf1_ok;
  if (t2pf1.extent_int(0) <= nmb1 || (two && !t2pf2_ok)) {cln_ldt = 0.0; return;}
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int is = id.is, ie = id.ie;
  auto la = cln_lacc;
  auto p1 = t2pf1, p2 = t2pf2;
  auto g1 = t2f1;
  auto iw_ = iw;
  auto ar = pmy_pack->pcoord->area.x1f;   // (the x1 face area of the clean fluxes)
  const bool lim = t2_ffct && ok1;   // (the limiter's R in iw slot Time2FctSlot(0))
  const int sr = Time2FctSlot(0);
  const Real crp = chat/c_light, cc2 = two ? c2 : 0.0;
  par_for("cln_bndacc", DevExeSpace(), 0, nmb1, 0, 1, id.ks, id.ke, id.js, id.je,
  KOKKOS_LAMBDA(const int m, const int s, const int k, const int j) {
    const int i = s ? ie + 1 : is;
    const Real g = ar(m,k,j,i)*(c1*(ok1 ? p1(m,0,k,j,i) : crp*g1(m,k,j,i)) +
                                cc2*p2(m,0,k,j,i));
    const Real rd = lim ? iw_(m,sr,k,j,(g > 0.0) ? i - 1 : i) : 1.0;
    la(m,2+s,k,j) = w*la(m,s,k,j) + rd*g;
  });
}

//! \fn Real RadiationVet::EbndEsrc
//! \brief ebnd_split: (chat/c) sum esrc V over the active cells, all ranks (0 without
//! esrc)

Real RadiationVet::EbndEsrc() {
  if (!esrc_on) return 0.0;
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int is = id.is, js = id.js, ks = id.ks;
  const int ni = id.nx1, nj = id.nx2, nk = id.nx3;
  const int nmb = pmy_pack->nmb_thispack;
  auto es_ = esrc;
  auto &vol = pmy_pack->pcoord->volume;
  Real l = 0.0;
  Kokkos::parallel_reduce("vet_ebnd_src", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          nmb*nk*nj*ni),
  KOKKOS_LAMBDA(const int n, Real &acc) {
    const int m = n/(nk*nj*ni);
    const int k = (n/(nj*ni))%nk + ks, j = (n/ni)%nj + js, i = n%ni + is;
    acc += es_(m,k,j,i)*vol(m,k,j,i);
  }, Kokkos::Sum<Real>(l));
#if MPI_PARALLEL_ENABLED
  Real g;
  MPI_Allreduce(&l, &g, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  l = g;
#endif
  return (chat/c_light)*l;
}

//! \fn void RadiationVet::EbndStart / EbndStage1 / EbndAccept
//! \brief Ebnd: L at the step start (face flux of K1), of the accepted stage-1 solve
//! (K2), and the step's applied energy at acceptance (after the dt_retry check)

void RadiationVet::EbndStart() {
  if (!ebnd_on) return;
  if (ebnd_split) ebnd_s1 = ebnd_s3_ok ? ebnd_s3 : EbndEsrc();
}

void RadiationVet::EbndStage1() {
  if (!ebnd_on) return;
  if (ebnd_split) ebnd_s2 = EbndEsrc();
}

void RadiationVet::EbndAccept(Real dt, bool t2step) {
  if (!ebnd_on) return;
  // the step's applied energy of the slopes l1 (K1), l2 (K2), l3 (K3)
  auto applied = [&](const Real l1, const Real l2, const Real l3) {
    if (t2step) {
      const Real g = kT2G;
      const Real wy = (t2_tab != 1) ? 0.5 : t2_w0;   // weight of Y1 in the stage-2 start
      return dt*((wy*t2_a21 + t2_c31)*l1 + (wy*g + t2_c32)*l2 + g*l3);
    }
    return dt*l3;
  };
  // clean: f0x1 is the face momentum; the energy through the physical x1 faces is
  // CleanBndAcc's
  Real lc[2] = {0.0, 0.0};
  if (!(cln_ldt > 0.0) && global_variable::my_rank == 0) {
    std::printf("CLEAN WARNING cycle %d: Ebnd without the step's boundary energy "
                "(CleanBndAcc), its faces take 0\n", pmy_pack->pmesh->ncycle);
  }
  if (cln_ldt > 0.0) {
    auto la = cln_lacc;
    auto &bcs = pmy_pack->pmb->mb_bcs;
    auto &id = pmy_pack->pmesh->mb_indcs;
    Kokkos::parallel_reduce("vet_ebnd_cln", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
        DevExeSpace(), {0, id.ks, id.js}, {pmy_pack->nmb_thispack, id.ke + 1, id.je + 1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, Real &si, Real &so) {
      for (int s = 0; s < 2; ++s) {
        const BoundaryFlag f = bcs.d_view(m, s ? BoundaryFace::outer_x1 :
                                                 BoundaryFace::inner_x1);
        if (f == BoundaryFlag::block || f == BoundaryFlag::periodic) continue;
        (s ? so : si) += la(m,s,k,j);
      }
    }, Kokkos::Sum<Real>(lc[0]), Kokkos::Sum<Real>(lc[1]));
    for (int s = 0; s < 2; ++s) lc[s] = CleanReduce(lc[s], false);
  }
  if (ebnd_split) {
    ebnd_acc_in += lc[0];
    ebnd_acc_out += lc[1];
    const Real s3 = EbndEsrc();
    ebnd_acc_src += applied(ebnd_s1, ebnd_s2, s3);
    ebnd_s3 = s3;
    ebnd_s3_ok = true;
    ebnd_acc = ebnd_acc_in - ebnd_acc_out;
    return;
  }
  ebnd_acc += lc[0] - lc[1];
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::Time2RotCorrect
//! \brief vet_gd_rotate_slope = correct, on the first step of a new vet_gd rotation block
//! with a valid K1 of the old set: one stage-type solve of dt_c =
//! vet_gd_rotate_slope_frac dt from U^n with the new set and old vector U^n - dt_c K1.
//! Its slope K1 + (I - dt_c A)^{-1} Delta f takes the closure jump where dt_c A is small
//! and keeps K1 where stiff; a failed solve keeps K1.  Returns true when a solve ran.

bool RadiationVet::Time2RotCorrect(Driver *pdrive, Real dt) {
  if (!t2_rotcorr) return false;
  t2_rotcorr = false;
  if (time_scheme != VET_TIME_HESDIRK2 || !Time2Active()) return false;
  FluidRef fl = FluidRef::Get(pmy_pack);
  const Real dtc = vgd_rot_cfrac*dt;
  // U^n into the stage-start copies (Time2FormStage(1) overwrites them again)
  Kokkos::deep_copy(DevExeSpace(), u1, u0);
  Kokkos::deep_copy(DevExeSpace(), t2f1, f0x1);
  if (trans_on) {
    Kokkos::deep_copy(DevExeSpace(), t2f2, f0x2);
    if (trans_x3) {Kokkos::deep_copy(DevExeSpace(), t2f3, f0x3);}
  }
  if (fl.on) {Kokkos::deep_copy(DevExeSpace(), fl.u1, fl.u0);}
  // mhd: U^n = (u1, b1) for the solve; b1 still holds the last step's start
  if (fl.on && fl.mhd) {
    Kokkos::deep_copy(DevExeSpace(), fl.pm->b1.x1f, fl.pm->b0.x1f);
    Kokkos::deep_copy(DevExeSpace(), fl.pm->b1.x2f, fl.pm->b0.x2f);
    Kokkos::deep_copy(DevExeSpace(), fl.pm->b1.x3f, fl.pm->b0.x3f);
  }
  // the old vector U^n - dt_c K1 (every channel)
  {
    auto k1_ = t2k1;
    auto in_ = t2inc;
    const Real c = -dtc;
    int nmb1 = static_cast<int>(t2inc.extent(0)) - 1;
    int n3 = static_cast<int>(t2inc.extent(2)), n2 = static_cast<int>(t2inc.extent(3));
    int n1 = static_cast<int>(t2inc.extent(4));
    par_for("vet_t2_rcinc", DevExeSpace(), 0, nmb1, 0, VET_T2_NK-1, 0, n3-1, 0, n2-1,
            0, n1-1,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
      in_(m,n,k,j,i) = c*k1_(m,n,k,j,i);
    });
    cln_lzc = c*cln_lzs1;   // (clean_cs_lz_abs = 2)
  }
  const bool pred_save = impl_pred;
  impl_pred = false;
  Mesh *pmesh = pmy_pack->pmesh;
  t2_inprobe = true;
  t2_solve = VET_T2S_STAGE1;
  t2_fail = false;
  nsub = 1;
  dt_sub = dtc;
  const Real its0 = impl_itsum;
  pdrive->ExecuteTaskList(pmesh, "vet_before_stagen", 1);
  pdrive->ExecuteTaskList(pmesh, "vet_stagen", 1);
  pdrive->ExecuteTaskList(pmesh, "vet_after_stagen", 1);
  vgd_rot_ncpass += static_cast<double>(impl_itsum - its0);
  t2_inprobe = false;
  impl_pred = pred_save;
  const bool ok = !t2_fail;
  if (ok) {
    Kokkos::deep_copy(DevExeSpace(), t2k1, t2k2);
    cln_lzs1 = cln_lzs2;   // (clean_cs_lz_abs = 2)
    if (t2_pfst && t2pf1.extent_int(0) > 0) {
      Kokkos::deep_copy(DevExeSpace(), t2pf1, t2pf2);
      t2pf1_ok = t2pf2_ok;
    }
    vgd_rot_ncorr += 1.0;
  } else {
    vgd_rot_ncfail += 1.0;
  }
  // U^n back
  Kokkos::deep_copy(DevExeSpace(), u0, u1);
  Kokkos::deep_copy(DevExeSpace(), f0x1, t2f1);
  if (trans_on) {
    Kokkos::deep_copy(DevExeSpace(), f0x2, t2f2);
    if (trans_x3) {Kokkos::deep_copy(DevExeSpace(), f0x3, t2f3);}
  }
  if (fl.on) {
    Kokkos::deep_copy(DevExeSpace(), fl.u0, fl.u1);
    FluidRefresh(pdrive);
    if (fl.mhd) EmagBuild(false);   // b0 untouched; same faces
  }
  t2_solve = VET_T2S_NONE;
  t2_fail = false;
  // the stage-1 solve keeps this build (vet_col predict | rebuild only)
  t2_vcreuse = ok && vgd_rot_creuse && vet_col && !vet_sc &&
               (t2_vcmode == 1 || t2_vcmode == 2);
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::Time2VetSourceInit
//! \brief time2_vet_source = relaxed | state | blend: VET source of the predict / rebuild
//! builds.  relaxed: a local backward-Euler exchange over dt (one step ahead of the
//! state: first order at P); state: the source of P itself; blend: S*(0) + w (S*(dt) -
//! S*(0)), w = x^2/(x0^2 + x^2), x = c dt rho kappa (time2_vet_source_x0).

void RadiationVet::Time2VetSourceInit(ParameterInput *pin, const bool act,
                                      const bool col) {
  const std::string key = col ? VetKeyName(pin, "time2_vet_source",
                                          {"time2_vet_col_source"})
                              : std::string("time2_vet_source");
  if (!act && !pin->DoesParameterExist("rad_vet", key)) {return;}
  // vet_sc: state; vet_col / vet_gd: blend (relaxed on a restart without the key)
  const char *dflt = col ? (global_variable::restart_run ? "relaxed" : "blend") : "state";
  std::string vs = act ?
      pin->GetOrAddString("rad_vet", key, dflt) :
      pin->GetString("rad_vet", key);
  if (vs.compare("relaxed") == 0) {
    t2_vcsrc = false;
  } else if (vs.compare("state") == 0) {
    t2_vcsrc = true;
  } else if (vs.compare("blend") == 0) {
    t2_vcsbl = true;
    const std::string kx = col ? VetKeyName(pin, "time2_vet_source_x0",
                                           {"time2_vet_col_source_x0"})
                               : std::string("time2_vet_source_x0");
    vcol_sx0 = pin->GetOrAddReal("rad_vet", kx, 1.0);
    if (!(vcol_sx0 > 0.0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/" << kx << " must be > 0" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<rad_vet>/" << key << " = '" << vs
              << "' is not a choice (relaxed | state | blend)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

int RadiationVet::Time2FctSlot(const int q) const {
  return q ? VET_IW_S2 : VET_IW_S1;
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationVet::CleanRstNq
//! \brief closure components per record of the clean extrap history and clean_rot_ring
//! (CleanRstSlot order): all CLN_N

int RadiationVet::CleanRstNq() const {
  return CLN_N;
}

// clean_closure = iterate: the closure in place is restart state (none in a file of
// another closure, Time2RstAccept)
int RadiationVet::CleanRstNcl() const {
  if (t2_rcl >= 0) return t2_rcl;  // a file's (Time2RstAccept)
  return cln_closit ? CleanNs(pmy_pack->pmesh->use_cubed_sphere, cln_rfh) : 0;
}

int RadiationVet::CleanRstNfv() const {
  // cubed sphere sc rays: their cln_fvi slots (before the first sweep: the file's);
  // version >= 6: the face half ranges cln_fh instead (none where every solve rebuilds
  // them)
  if (cln_cssc && t2_rstv >= 6) return CleanFhBuilt() ? 0 : 12;
  if (cln_cssc && cln_fvi.extent_int(0) > 0) return cln_fvi.extent_int(1);
  if (cln_cssc) return (cln_fvn > 0) ? cln_fvn : t2_nfv;
  return cln_scmem ? cln_top.extent_int(1) : 0;
}

//----------------------------------------------------------------------------------------
//! \fn RadiationVet::Time2Rst*
//! \brief the restart channels: K1 (VET_T2_NK), vet_prev (VET_T2_NVET, vet_sc), then the
//! clean scheme's state; files of version < 4 also hold ipred2 and extrap2 channels
//! (read past)

int RadiationVet::Time2RstNchWant() {
  // version >= 4: no ipred2 / extrap2 channels, ring offsets dense (CleanRingRstCopy)
  const bool sk = t2_rstv >= 4;
  const int np = (impl_pord == 2) ? 5 : 3;
  const int nt = (time_scheme != VET_TIME_HESDIRK2) ? 0 :
         VET_T2_NK + ((impl_pred && !sk) ? np : 0) + (vet_sc ? VET_T2_NVET : 0)
         + ((vet_col && t2_vcmode == 4 && !sk) ? VET_X2_NCH : 0);   // extrap2
  return nt + ((force_ref == VET_FREF_CLEAN) ? 4 + cln_lzt : 0)   // the reference, set
         + ((force_ref == VET_FREF_WB_MD) ? 10 : 0)           // wb_arad_md: 9, set
         + ((cln_extrap > 0) ? 2*CleanRstNq() + 1 : 0)    // extrap history
         + CleanRstNcl()                                     // iterate: the closure
         + (cln_ring ? (t2_rstv < 4 ? vra_k*CleanRstNq() : CleanRingRstN()) + 1
                                : 0)                          // clean_rot_ring
         + ((cln_fsolve && cln_fkin) ? 6 : 0)       // stag fkin: d and its K1
         + CleanRstNfv();                                    // the rays' intensities
}

bool RadiationVet::Time2RstAccept(const int nfile, const int ver) {
  // version 1 holds d without its slope K1_d: K1_d = 0, backward-Euler first step
  t2_rstd3 = (ver < 2) && cln_fsolve && cln_fkin;
  // version 5 on the cubed sphere (sc ray intensities) is refused
  if (cln_cssc && ver == 5) {
    std::cout << "### FATAL ERROR <rad_vet> restart: a version-5 clean file on the cubed "
              << "sphere (the sc rays' intensities) is no longer read" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  t2_rstv = ver;
  t2_rcl = -1;
  int nb = Time2RstNchWant() - CleanRstNfv() - (t2_rstd3 ? 3 : 0);
  t2_nfv = (cln_cssc || cln_scmem) ? nfile - nb : 0;
  // a file of another clean_closure: the iterate closure is skipped, or missing (the
  // first solve takes the predicted closure)
  const int ns = CleanNs(pmy_pack->pmesh->use_cubed_sphere, cln_rfh);
  auto fits = [&](const int n) {
    return cln_cssc ? (nfile - n == 0 || nfile - n == 12) : (nfile == n);
  };
  if (!cln_scmem && !fits(nb) && fits(cln_closit ? nb - ns : nb + ns)) {
    t2_rcl = cln_closit ? 0 : ns;
    nb += cln_closit ? -ns : ns;
    t2_nfv = cln_cssc ? nfile - nb : 0;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> NOTE: restart from a file of another clean_closure: "
                << (cln_closit ? "no closure in it, the first solve takes the predicted "
                    "one" : "its closure (clean_closure = iterate) is not read")
                << std::endl;
    }
  }
  // (the ring's half ranges, cln_rfh: no cln_fh channels; a file of a run with them is
  // refused without them, and one without them with them: 1 .. 11 channels over)
  const bool fho = cln_rfh ? t2_nfv != 0 : (t2_nfv > 0 && t2_nfv < 12);
  if (cln_cssc && ver >= 6 && fho) {
    std::cout << "### FATAL ERROR <rad_vet> restart: a file " << (cln_rfh ? "not " : "")
              << "written with the ring's half ranges (vet_rotate_avg_ema_wb on the "
              << "cubed sphere; " << t2_nfv << " channels over)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // cubed sphere, version >= 6: the face half ranges are the last 0 or 12 channels; a
  // surplus is a clean_rot_ring file: with the ring off its channels are skipped, else
  // refused
  t2_rskip = 0;
  if (cln_cssc && ver >= 6 && t2_nfv > 12) {
    t2_rskip = cln_ring ? -1 : CleanRingRstN() + 1;
    const int nh = t2_nfv - t2_rskip;
    if (t2_rskip < 0 || (nh != 0 && nh != 12)) {
      std::cout << "### FATAL ERROR <rad_vet> restart: the clean block has " << nfile
                << " channels, not this run's layout (" << nb << " + 0 or 12) nor it "
                << "with the ring's offsets" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    t2_nfv = nh;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> NOTE: restart with clean_rot_ring off: the file's ring "
                << "offsets (" << t2_rskip << " channels) are not read" << std::endl;
    }
  }
  return (nfile == nb + t2_nfv + t2_rskip) && (t2_nfv >= 0);
}
int RadiationVet::Time2RstNch() {
  t2_rstv = 7;   // (the writer: version 7, after a restart from an older file too)
  t2_rcl = -1;
  // the clean state (and K1 with cln_k1ok) is written in backward-Euler steps too;
  // restart.cpp stores t2_ok
  return Time2RstNchWant();
}

// m0 > 0 packs MeshBlocks m0 .. m0+nmb-1 into a(0 .. nmb-1) (per-MB staging)
void RadiationVet::Time2RstPack(DvceArray5D<Real> &a, int nmb, int m0, int w0) {
  auto mb = std::make_pair(0, nmb);
  auto ms = std::make_pair(m0, m0 + nmb);
  auto AL = Kokkos::ALL;
  const bool h2 = (time_scheme == VET_TIME_HESDIRK2);
  // a holds the record's channels w0 .. w1-1 (the writer streams the record in windows);
  // put copies the record channels c0 .. c0+n-1 = src's s0 .. s0+n-1 that fall inside,
  // one sets channel cc to v
  const int w1 = w0 + a.extent_int(1);
  auto put = [&](const int c0, const int n, const auto &src, const int s0) {
    const int lo = std::max(c0, w0), hi = std::min(c0 + n, w1);
    if (lo < hi) {
      Kokkos::deep_copy(Kokkos::subview(a, mb, std::make_pair(lo - w0, hi - w0), AL, AL,
                                        AL),
                        Kokkos::subview(src, ms,
                                        std::make_pair(s0 + lo - c0, s0 + hi - c0),
                                        AL, AL, AL));
    }
  };
  auto one = [&](const int cc, const Real v) {
    if (cc >= w0 && cc < w1) {
      Kokkos::deep_copy(Kokkos::subview(a, mb, cc - w0, AL, AL, AL), v);
    }
  };
  const int na3 = static_cast<int>(a.extent(2)), na2 = static_cast<int>(a.extent(3));
  const int na1 = static_cast<int>(a.extent(4)), nw = a.extent_int(1);
  int c = 0;
  if (h2) {
    put(c, VET_T2_NK, t2k1, 0);
    c += VET_T2_NK;
  }
  // version >= 4 layout: no predictor / extrap2 channels
  if (h2 && vet_sc) {
    put(c, VET_T2_NVET, vet_prev, 0);
    c += VET_T2_NVET;
  }
  // force_reference = clean (wb_arad_md): the reference force (and potentials) and
  // whether it is set
  if (force_ref == VET_FREF_CLEAN || force_ref == VET_FREF_WB_MD) {
    const int nr = (force_ref == VET_FREF_WB_MD) ? 9 : 3 + cln_lzt;
    const bool have = (cln_fref.extent_int(0) >= m0 + nmb);
    if (have) put(c, nr, cln_fref, 0);
    one(c + nr, (have && cln_fref_set) ? 1.0 : 0.0);
    c += nr + 1;
  }
  // clean_closure = extrap*: the ray moments at t^n and at the end of the step, then
  // t^n, that end time and the flags at i = 0..4
  const int nq = CleanRstNq();
  if (cln_extrap > 0) {
    const bool have = (cln_mn.extent_int(0) >= m0 + nmb);
    for (int n = 0; n < nq && have; ++n) {
      const int s = CleanRstSlot(n);
      if (s >= cln_mn.extent_int(1)) {   // (a slot this mesh lacks: 0)
        one(c + n, 0.0);
        one(c + nq + n, 0.0);
        continue;
      }
      put(c + n, 1, cln_mn, s);
      put(c + nq + n, 1, cln_mc, s);
    }
    auto a_ = a;
    const int cq = c + 2*nq - w0;
    const Real v0 = cln_tn, v1 = cln_tm, v2 = (have && cln_have_n) ? 1.0 : 0.0;
    const Real v3 = (have && cln_have_c) ? 1.0 : 0.0, v4 = cln_cand_be ? 1.0 : 0.0;
    if (cq >= 0 && cq < nw) {
      par_for("vet_t2_clnrst", DevExeSpace(), 0, nmb - 1, 0, na3 - 1, 0, na2 - 1, 0,
              na1 - 1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        a_(m,cq,k,j,i) = (i == 0) ? v0 : ((i == 1) ? v1 : ((i == 2) ? v2 :
                         ((i == 3) ? v3 : ((i == 4) ? v4 : 0.0))));
      });
    }
    c += 2*nq + 1;
  }
  // clean_closure = iterate: the closure in place, the last iterate's (CleanRstNcl)
  if (CleanRstNcl() > 0) {
    if (cln_s.extent_int(0) >= m0 + nmb) put(c, CleanRstNcl(), cln_s, 0);
    c += CleanRstNcl();
  }
  // clean_rot_ring: the offsets and the pending Delta of the fields the ring carries
  // (interior cells, dense: CleanRingRstCopy), then the ring's flag, block, the pending
  // block and its time at i = 0..3
  if (cln_ring) {
    const int nd = CleanRingRstN();
    const bool have = cln_rok && (cln_rr.extent_int(0) >= m0 + nmb);
    CleanRingRstCopy(a, c - w0, nmb, m0, true, have);
    auto a_ = a;
    const int cq = c + nd - w0;
    const Real v0 = have ? 1.0 : 0.0, v1 = static_cast<Real>(cln_rb);
    const Real v2 = static_cast<Real>(cln_rpb), v3 = cln_rpt;
    if (cq >= 0 && cq < nw) {
      par_for("vet_t2_ringrst", DevExeSpace(), 0, nmb - 1, 0, na3 - 1, 0, na2 - 1, 0,
              na1 - 1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        a_(m,cq,k,j,i) = (i == 0) ? v0 : ((i == 1) ? v1 : ((i == 2) ? v2 :
                         ((i == 3) ? v3 : 0.0)));
      });
    }
    c += nd + 1;
  }
  // clean_stag_fkin: the energy flux memory d of the three face families (CLF_N + 3),
  // then its slope K1 (CLF_KD); before the first solve, those read from the file
  if (cln_fsolve && cln_fkin) {
    for (int d = 0; d < 3; ++d) {
      auto fc = (d == 0) ? cln_fc1 : ((d == 1) ? cln_fc2 : cln_fc3);
      if (!cln_init && cln_rstd.extent_int(0) >= m0 + nmb) {
        put(c, 6, cln_rstd, 0);
        break;
      }
      if (fc.extent_int(0) < m0 + nmb) continue;
      put(c + d, 1, fc, CLF_N + 3);
      put(c + 3 + d, 1, fc, CLF_KD);
    }
    c += 6;
  }
  // (a window that does not reach the rays' intensities: restart.cpp windows only a
  // record without them)
  if (w0 > 0 || w1 < c + CleanRstNfv()) return;
  // the sc rays' intensities the next build reads, last (unless a has no room for them:
  // restart.cpp copies them itself)
  if (a.extent_int(1) < c + CleanRstNfv()) return;
  if (cln_cssc && CleanRstNfv() > 0) {   // (version 6: cln_fh)
    auto af = Kokkos::subview(a, mb, std::make_pair(c, c + 12), AL, AL, AL);
    if (cln_fh.extent_int(0) >= m0 + nmb) {
      Kokkos::deep_copy(af, Kokkos::subview(cln_fh, ms, AL, AL, AL, AL));
    } else {
      Kokkos::deep_copy(af, 0.0);
    }
  } else if (CleanRstNfv() > 0 && CleanRstFv().extent_int(0) >= m0 + nmb) {
    Kokkos::deep_copy(Kokkos::subview(a, mb, std::make_pair(c, c + CleanRstNfv()), AL, AL,
                                      AL),
                      Kokkos::subview(CleanRstFv(), ms, AL, AL, AL, AL));
  }
}

void RadiationVet::Time2RstSet(int ch, const HostArray4D<Real> &w, int nmb) {
  auto mb = std::make_pair(0, nmb);
  auto AL = Kokkos::ALL;
  const bool h2 = (time_scheme == VET_TIME_HESDIRK2);
  // a file's predictor / extrap2 channels are read past (a NOTE once if nonzero)
  auto skip = [&w]() {
    static bool said = false;
    for (size_t q = 0; q < w.span() && !said && global_variable::my_rank == 0; ++q) {
      if (w.data()[q] == 0.0) continue;
      std::cout << "<rad_vet> NOTE: the restart file's implicit predictor / extrap2 "
                << "channels are skipped" << std::endl;
      said = true;
    }
  };
  // the last t2_nfv channels: the sc rays' intensities (their number is not
  // known before the first sweep, Time2RstAccept)
  const int nfv0 = Time2RstNchWant() - CleanRstNfv() - (t2_rstd3 ? 3 : 0) + t2_rskip;
  if (cln_cssc && t2_rstv >= 6 && ch >= nfv0) {   // (cln_fh)
    if (cln_fh.extent_int(0) != nmb) {
      Kokkos::realloc(cln_fh, nmb, 12, w.extent_int(1), w.extent_int(2), w.extent_int(3));
    }
    DeepCopyAcross(Kokkos::subview(cln_fh, mb, ch - nfv0, AL, AL, AL), w);
    return;
  }
  if ((cln_cssc || cln_scmem) && ch >= nfv0) {
    auto &fv = CleanRstFv();
    if (fv.extent_int(0) != nmb || fv.extent_int(1) != t2_nfv) {
      Kokkos::realloc(fv, nmb, t2_nfv, w.extent_int(1), w.extent_int(2), w.extent_int(3));
    }
    DeepCopyAcross(Kokkos::subview(fv, mb, ch - nfv0, AL, AL, AL), w);
    return;
  }
  if (h2 && ch < VET_T2_NK) {
    DeepCopyAcross(Kokkos::subview(t2k1, mb, ch, AL, AL, AL), w);
    return;
  }
  if (h2) ch -= VET_T2_NK;
  if (h2 && impl_pred && t2_rstv < 4) {
    const int np = (impl_pord == 2) ? 5 : 3;
    if (ch < np) {
      skip();
      return;
    }
    ch -= np;
  }
  if (h2 && vet_sc) {
    if (ch < VET_T2_NVET) {
      DeepCopyAcross(Kokkos::subview(vet_prev, mb, ch, AL, AL, AL), w);
      return;
    }
    ch -= VET_T2_NVET;
  }
  if (force_ref == VET_FREF_CLEAN || force_ref == VET_FREF_WB_MD) {
    const int nr = (force_ref == VET_FREF_WB_MD) ? 9 : 3 + cln_lzt;
    if (ch < nr) {
      if (cln_fref.extent_int(0) != nmb) {
        Kokkos::realloc(cln_fref, nmb, nr, w.extent_int(1), w.extent_int(2),
                        w.extent_int(3));
      }
      DeepCopyAcross(Kokkos::subview(cln_fref, mb, ch, AL, AL, AL), w);
      return;
    }
    if (ch == nr) {
      cln_fref_set = (w(0,0,0,0) > 0.5);
      return;
    }
    ch -= nr + 1;
  }
  const int nq = CleanRstNq();
  if (cln_extrap > 0) {
    if (ch < 2*nq) {
      auto &mx = (ch < nq) ? cln_mn : cln_mc;
      const int ns = CleanNs(pmy_pack->pmesh->use_cubed_sphere, cln_rfh);
      const int s = CleanRstSlot(ch%nq);
      if (mx.extent_int(0) != nmb) {
        Kokkos::realloc(mx, nmb, ns, w.extent_int(1), w.extent_int(2), w.extent_int(3));
      }
      if (s < ns) DeepCopyAcross(Kokkos::subview(mx, mb, s, AL, AL, AL), w);
      return;
    }
    if (ch == 2*nq) {
      cln_tn = w(0,0,0,0);
      cln_tm = w(0,0,0,1);
      cln_have_n = (w(0,0,0,2) > 0.5);
      cln_have_c = (w(0,0,0,3) > 0.5);
      cln_cand_be = (w(0,0,0,4) > 0.5);
      return;
    }
    ch -= 2*nq + 1;
  }
  if (CleanRstNcl() > 0) {   // (clean_closure = iterate: cln_s, with its ghosts)
    if (ch < CleanRstNcl()) {
      if (!cln_closit) return;   // (an iterate run's file, read by another closure)
      if (cln_s.extent_int(0) != nmb) {
        Kokkos::realloc(cln_s, nmb, CleanRstNcl(), w.extent_int(1), w.extent_int(2),
                        w.extent_int(3));
      }
      DeepCopyAcross(Kokkos::subview(cln_s, mb, ch, AL, AL, AL), w);
      cln_cip = true;
      return;
    }
    ch -= CleanRstNcl();
  }
  if (t2_rskip > 0) {  // the file's ring with the ring off: not read
    if (ch < t2_rskip) return;
    ch -= t2_rskip;
  }
  if (cln_ring && t2_rstv >= 4) {   // (version 4, CleanRingRstCopy)
    if (ch < CleanRingRstN()) {
      DvceArray5D<Real> t("vet_t2_ringrd", nmb, 1, w.extent_int(1), w.extent_int(2),
                          w.extent_int(3));
      DeepCopyAcross(Kokkos::subview(t, AL, 0, AL, AL, AL), w);
      CleanRingRstCopy(t, -ch, nmb, 0, false, true);
      return;
    }
    ch -= CleanRingRstN() - vra_k*nq;  // the flag slab
  }
  if (cln_ring) {
    const int nr = (vra_k - 1)*nq;
    if (ch < nr + nq) {
      // (the interior cells of the carried fields, CleanRingMoments)
      auto &rx = (ch < nr) ? cln_rr : cln_rpd;
      vra_gen += 1;
      auto &id = pmy_pack->pmesh->mb_indcs;
      const int is = id.is, js = id.js, ks = id.ks;
      Kokkos::Array<int, CLN_N> fm;
      const int nf = CleanRingFields(fm.data());
      const int f = fm[CleanRstSlot((ch < nr) ? ch%nq : ch - nr)];
      const int nx = (ch < nr) ? (vra_k - 1)*nf : nf;
      if (rx.extent_int(0) != nmb) {
        Kokkos::realloc(rx, nmb, nx, id.nx3, id.nx2, id.nx1);
      }
      if (f < 0) return;
      const int q = (ch < nr) ? (ch/nq)*nf + f : f;
      DvceArray4D<Real> t("vet_t2_ringrd", nmb, w.extent_int(1), w.extent_int(2),
                          w.extent_int(3));
      DeepCopyAcross(t, w);
      auto rx_ = rx;
      par_for("vet_t2_ringset", DevExeSpace(), 0, nmb - 1, 0, id.nx3 - 1, 0, id.nx2 - 1,
              0, id.nx1 - 1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        rx_(m,q,k,j,i) = t(m,k+ks,j+js,i+is);
      });
      return;
    }
    if (ch == nr + nq) {
      cln_rok = (w(0,0,0,0) > 0.5);
      cln_rb = static_cast<int64_t>(std::llround(w(0,0,0,1)));
      cln_rpb = static_cast<int64_t>(std::llround(w(0,0,0,2)));
      cln_rpt = w(0,0,0,3);
      return;
    }
    ch -= nr + nq + 1;
  }
  if (cln_fsolve && cln_fkin) {   // (CleanInit copies it to the face arrays)
    const int nd = t2_rstd3 ? 3 : 6;
    if (ch < nd) {
      if (cln_rstd.extent_int(0) != nmb) {   // (zeroed: K1_d = 0 for a version-1 file)
        Kokkos::realloc(cln_rstd, nmb, 6, w.extent_int(1), w.extent_int(2),
                        w.extent_int(3));
      }
      DeepCopyAcross(Kokkos::subview(cln_rstd, mb, ch, AL, AL, AL), w);
      return;
    }
    ch -= nd;
  }
  // extrap2 (version < 4): read past
  if (h2 && vet_col && t2_vcmode == 4) skip();
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::Time2AutoBeBegin
//! \brief time2_auto_be: called by the Driver at the cycle start (after VetGdRotStep).
//! A skip step is put in the state the stage fallback leaves (Time2Restore): no valid
//! slope (Heun hydro, then the BESTORE solve), cold predictor, the dt_retry arm of a
//! time2_reject redo.  A retried cycle calls it again with the same state.

bool RadiationVet::Time2AutoBeBegin() {
  t2_abe_now = false;
  // t2_beonly: every step (the first too) is a skip step
  const bool bo = t2_beonly && (time_scheme == VET_TIME_HESDIRK2);
  if (!bo && (t2_abe <= 0 || time_scheme != VET_TIME_HESDIRK2 || t2_abe_left <= 0)) {
    return false;
  }
  if (!bo && !Time2Active()) return false;   // no slope anyway (the BE path follows)
  t2_abe_now = true;
  t2_ok = false;
  t2_fail = false;
  t2_solve = VET_T2S_NONE;
  pred_ok = false;
  pred2_ok = false;
  // as the redo of a rejected stage under time2_reject = retry (arms dt_retry_t2_be)
  if (t2_reject == 1 && dtr_on) {dtr_t2hit = true;}
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::Time2AutoBeEnd
//! \brief an accepted cycle of a hesdirk2 run.  kind 0 stage step, 1 stage fallback +
//! BE redo, 2 time2_auto_be step, 3 other BE step, 4 discarded dt_retry attempt (timing
//! only); tw_a / tw_b: fenced wall seconds of the attempt or step, and of the redo

void RadiationVet::Time2AutoBeEnd(int kind, Real tw_a, Real tw_b) {
  if (kind == 0) {t2w_clean += tw_a;}
  if (kind == 1) {t2w_fail += tw_a; t2w_redo += tw_b;}
  if (kind == 2) {t2w_abe += tw_a;}
  if (kind == 3) {t2w_be += tw_a;}
  if (kind == 4) {t2w_retry += tw_a; return;}
  const bool was_skip = t2_abe_now;
  t2_abe_now = false;
  if (t2_abe <= 0) return;
  if (kind == 0) {
    if (t2_abe_len > 0) {t2_abe_ntry += 1.0; t2_abe_nok += 1.0;}
    t2_abe_len = 0;
    t2_abe_left = 0;
  } else if (kind == 1) {
    if (t2_abe_len > 0) {t2_abe_ntry += 1.0;}
    t2_abe_len = (t2_abe_len == 0) ? 1 : std::min(2*t2_abe_len, t2_abe);
    t2_abe_left = t2_abe_len;
  } else if (kind == 2 && was_skip) {
    t2_abe_n += 1.0;
    t2_abe_left -= 1;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::Time2Report
//! \brief end-of-run counters of the hesdirk2 steps (rank 0)

void RadiationVet::Time2Report() {
  if (time_scheme != VET_TIME_HESDIRK2) return;
  if (global_variable::my_rank != 0) return;
  std::cout << "<rad_vet> time_scheme=" << (t2_beonly ? "hesdirk2(be_only)" : "hesdirk2")
            << ": stage steps=" << t2_nstep
            << " backward-Euler steps=" << t2_nbe << " stage fallbacks=" << t2_nfall
            << " vet clips=" << t2_nclip;
  if (t2_k1stiff > 0.0) {
    std::cout << " k1_stiff(s=" << t2_k1stiff << ") cells(rank 0)=" << t2_k1stiff_n;
  }
  if (t2_ffct) {  // global sums
    std::cout << " face_fct calls=" << t2_ffct_calls << " skipped(no face flux)="
              << t2_ffct_nop << " K1 flux from U^n=" << t2_ffct_fill << " limited cells="
              << t2_ffct_cells << " max rounds="
              << t2_ffct_rmax << " withheld E x vol=" << t2_ffct_de << " left (Q<=0)="
              << t2_ffct_left << " emptied (E < " << t2_ffct_emin << " E^n) cells="
              << t2_ffct_nemp << " min E/E^n="
              << t2_ffct_qmin << " (rank 0: rounds=" << zl_rounds << " R exchanges="
              << zl_halos;
    std::cout << ")";
  }
  if (vet_col && t2_vcmode != 0) {
    std::cout << " time2_vet_col=" << ((t2_vcmode == 1) ? "predict" :
                                       ((t2_vcmode == 2) ? "rebuild" :
                                        ((t2_vcmode == 4) ? "extrap2" : "extrap")))
              << " fallbacks=" << t2_vcnfb;
    if (t2_vcsrc && t2_vcmode < 3) {std::cout << " source=state";}
    if (t2_vcsbl && t2_vcmode < 3) {std::cout << " source=blend x0=" << vcol_sx0;}
  }
  if (vet_sc && t2_vsmode != 0) {
    std::cout << " time2_vet_sc=" << ((t2_vsmode == 1) ? "predict" : "rebuild")
              << " fallbacks=" << t2_vcnfb;
    if (t2_vcsrc) {std::cout << " source=state";}
    if (t2_vcsbl) {std::cout << " source=blend x0=" << vcol_sx0;}
    if (!t2_vcsrc && !t2_vcsbl) {
      std::cout << (vcol_srelax ? (vcol_sesc ? " source=relaxed_esc" : " source=relaxed")
                                : " source=gas");
    }
  }
  if (t2_pcap > 0.0) {
    std::cout << " time2_pred_cap=" << t2_pcap << " limited cell-builds=" << t2_pnlim
              << " of " << t2_pncell;
  }
  if (vgd_rot_slope == 3) {
    std::cout << " vet_gd_rotate_slope=correct solves=" << vgd_rot_ncorr << " failed="
              << vgd_rot_ncfail << " Picard passes=" << vgd_rot_ncpass << " (frac "
              << vgd_rot_cfrac << ", reuse "
              << (vgd_rot_creuse ? "on" : "off") << ")";
  }
  if (t2_accfb) {std::cout << " accepted-with-fallback stages=" << t2_accfb_n;}
  if (t2_adm_pl) {std::cout << " time2_adm_poslim admitted stages=" << t2_adm_pl_n;}
  // VetGdAlignDt calls that shortened dt (NewTimeStep calls, incl. the start-up ones)
  if (vgd_rot_align) {std::cout << " rotate_align shortened dt=" << vgd_rot_nalign;}
  if (cln_ring) {
    std::cout << " clean_rot_ring guard cells=" << cln_rguard_sum
              << " off-boundary shifts=" << cln_roff_sum;
  }
  std::cout << std::endl;
  if (t2_abe > 0) {
    std::cout << "<rad_vet> time2_auto_be=" << t2_abe << ": skipped-stage BE steps="
              << t2_abe_n << " stage retries=" << t2_abe_ntry << " (succeeded "
              << t2_abe_nok << ")" << std::endl;
  }
  // the wall-time split (rank 0, fenced; accepted cycles by kind)
  std::cout << "<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps " << t2w_clean
            << " | failed stage attempts " << t2w_fail << " | their BE redo " << t2w_redo
            << " | auto_be steps " << t2w_abe << " | other BE steps " << t2w_be
            << " | discarded dt_retry attempts " << t2w_retry << std::endl;
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::VsrcNoes
//! \brief vet_source_noesrc in force for the current solve: every backward-Euler solve
//! (under hesdirk2 unless vet_source_noesrc_h2be = false), never a stage solve

bool RadiationVet::VsrcNoes() const {
  if (!vsrc_noes) return false;
  if (time_scheme == VET_TIME_BE) return true;
  return vsrc_noes_h2be && (t2_solve != VET_T2S_STAGE1) && (t2_solve != VET_T2S_STAGE2);
}

}  // namespace radvet
