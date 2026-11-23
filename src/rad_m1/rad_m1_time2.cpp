//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_time2.cpp
//! \brief <rad_m1>/time_scheme = hesdirk2: second order in time for the coupled gas +
//! implicit M1 system (docs/dev/rad_m1_time2_design.md sect. 5,
//! tests_m1/runs_3x_hesdirk2).
//!
//! H-ESDIRK2 is an IMEX Runge-Kutta scheme whose explicit part is AthenaK's Heun (rk2)
//! hydro, unchanged, and whose implicit part is the stiffly accurate, L-stable ESDIRK
//!   A = [[0,0,0],[1-g,g,0],[1/2,b2,g]],  g = 1 - 1/sqrt 2,  b2 = g/(2(1-g)),
//! with the first slope first-same-as-last (FSAL).  Each implicit stage is the existing
//! backward-Euler ImplicitSolve with dt -> g dt and a SUPPLIED old vector
//!   rhs = (stage start state) + t2inc,
//! where the stage start state is what the Heun stage produced (it stays the solve's
//! first iterate, and the only state the EOS and the opacities are evaluated at), and
//!   stage 1: t2inc = (1-g) dt K1,   stage 2: t2inc = dt (g/2 K1 + (b2 - g/2) K2).
//! The slope of a solve is K = (Y - rhs)/(g dt); the stage-2 one is the next K1.
//!
//! A step is backward Euler (with K1 = (Y - rhs)/dt stored) when no valid K1 exists:
//! the first step, a restart from a file without the slope, and the step after a failed
//! stage, which is redone from U^n.

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
#include "rad_m1/m1_fluid.hpp"
#include "utils/deep_copy_across.hpp"
#include "rad_m1/rad_m1.hpp"
#include "utils/recipe.hpp"
#include "rad_m1/rad_m1_parfor.hpp"
#include "rad_m1/rad_m1_implicit.hpp"
#include "rad_m1/rad_m1_clean.hpp"
#include "rad_m1/rad_m1_opacity.hpp"
#include "utils/eint_from_cons.hpp"
#include "utils/imex_sdirk2.hpp"

namespace radm1 {

using deep_copy_across::DeepCopyAcross;

namespace {
const Real kT2G = imex_sdirk2::Gamma();   // shared tableau
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::Time2Init
//! \brief read <rad_m1>/time_scheme (be | hesdirk2) and allocate the stage state.  The
//! key is always present here: ImplicitInit resolves its default (hesdirk2 where it is
//! accepted, be elsewhere and on restarts whose file lacks it; m1-defaults2).

void RadiationM1::Time2Init(ParameterInput *pin) {
  time_scheme = M1_TIME_BE;
  if (pin->DoesParameterExist("rad_m1", "time_scheme")) {
    std::string s = pin->GetString("rad_m1", "time_scheme");
    if (s.compare("be") == 0) {
      time_scheme = M1_TIME_BE;
    } else if (s.compare("hesdirk2") == 0) {
      time_scheme = M1_TIME_HESDIRK2;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/time_scheme = '" << s
                << "' is not a choice (be | hesdirk2)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // the flux table is applied at the stage times (Time2FormStage)
  if (ftab_on && time_scheme != M1_TIME_HESDIRK2 && !t2_beonly) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<rad_m1>/implicit_flux_x1min_table needs time_scheme = "
              << "hesdirk2 (the flux is set at the stage times)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // plain be with t2_beonly (ImplicitInit) runs as hesdirk2, every step on
  // the backward-Euler path (Time2AutoBeBegin)
  if (time_scheme == M1_TIME_BE && t2_beonly) {
    time_scheme = M1_TIME_HESDIRK2;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> time_scheme = be: the hesdirk2 backward-Euler path on every "
                << "step (bematch-1030; be_legacy = true for the old step)" << std::endl;
    }
  } else if (time_scheme == M1_TIME_HESDIRK2 && t2_beonly &&
             global_variable::my_rank == 0) {
    std::cout << "<rad_m1> time2_be_only = true: every step takes the backward-Euler path"
              << std::endl;
  }
  if (time_scheme == M1_TIME_BE) {t2_beonly = false; return;}

  std::string integ = pin->GetString("time", "integrator");
  if (integ.compare("rk2") != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<rad_m1>/time_scheme = hesdirk2 needs <time>/integrator "
              << "= rk2 (its explicit part IS the Heun hydro)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // the stage coefficients (see rad_m1.hpp t2_tab).  hesdirk2: a21 = 1-g, Heun average,
  // then g/2 K1 + (b2 - g/2) K2.  trbdf2 (ESDIRK form of TR-BDF2, c = (0, 2g, 1),
  // A = [[0,0,0],[g,g,0],[w,w,g]], w = sqrt(2)/4): a21 = g, the hydro's stage-2 weights
  // 1/4 (Y2) and 3/4 (U^n), then (w - g/4)(K1 + K2).
  {
    const Real g = kT2G;
    t2_tab = 0;
    if (pin->DoesParameterExist("rad_m1", "time2_tableau")) {
      std::string tb = pin->GetString("rad_m1", "time2_tableau");
      if (tb.compare("hesdirk2") == 0) {
        t2_tab = 0;
      } else if (tb.compare("trbdf2") == 0) {
        t2_tab = 1;
      } else if (tb.compare("sdirk2") == 0) {   // seamE-1106
        t2_tab = 2;
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<rad_m1>/time2_tableau = '" << tb
                  << "' is not a choice (hesdirk2 | trbdf2)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    // the coefficients from the shared tableau (utils/imex_sdirk2.hpp,
    // the same expressions), also used by <problem>/ck_time_scheme = hesdirk2
    const imex_sdirk2::Tableau tb = imex_sdirk2::Make(t2_tab);
    (void) g;
    t2_a21 = tb.a21;
    t2_w0 = tb.w0;
    t2_w1 = tb.w1;
    t2_c31 = tb.c31;
    t2_c32 = tb.c32;
    if (t2_tab == 2 && global_variable::my_rank == 0) {
      std::cout << "<rad_m1> time2_tableau = sdirk2: no explicit K1, Heun hydro"
                << std::endl;
    }
    if (t2_tab == 2 && coupling && pmy_pack->pmhd != nullptr) {   // (Time2GasTime)
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/time2_tableau = sdirk2 is not wired for mhd"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (t2_tab == 1) {
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_m1> time2_tableau = trbdf2: TR-BDF2 implicit stages (stage "
                  << "order 2), RK2 hydro with c2 = 2 - sqrt(2)" << std::endl;
      }
    }
  }
  // closure = vet_col (m1-sph2) sets tau_closure too: its tensor is a formal solution
  // like vet_sc's, built at U^n by Time2VetStart and kept for both stages
  if (tau_closure && !vet_col) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<rad_m1>/time_scheme = hesdirk2 is not wired for "
              << "closure = tau" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // time2_enth_vel = old | start | central (default central): how the enthalpy face
  // coefficient a_f of a stage solve is built under implicit_vimp + implicit_enthalpy =
  // plm (runs_3x_hesdirk2: `old` = plm a_f of the old-vector velocity, unstable at
  // P = 100, tau_lambda = 1e5; `start` = plm a_f of the stage-start a plus the rest as a
  // face mean, also unstable there; `central` = a_f the mean of the two cells, E_f plm)
  t2_afmode = 2;
  // time2_enth_be: with time2_enth_vel = central the stage solves build a_f as the mean
  // of the two cells, but the backward-Euler solves of the same run (the start step, the
  // time2_start slope probe, the rotation probe, the step after a failed stage) take
  // implicit_enthalpy's own a_f (plm: the 4-point plm mean).  That changes the operator
  // for one step of length dt, an O(dt) term of the T ladders (rw_p10t100 E / eint 1.78 /
  // 1.84 at l 16-32).  true: those solves take the stage rule too. DEFAULT true on a
  // fresh run; a restart whose file lacks the key keeps false (FIN3, bitwise) -- the same
  // rule as implicit_cs2_vimp_npass (fresh 3, restart 0).
  t2_enth_be = pin->GetOrAddBoolean("rad_m1", "time2_enth_be",
                                    !global_variable::restart_run);
  if (pin->DoesParameterExist("rad_m1", "time2_enth_vel")) {
    std::string ev = pin->GetString("rad_m1", "time2_enth_vel");
    if (ev.compare("old") == 0) {
      t2_afmode = 0;
    } else if (ev.compare("start") == 0) {
      t2_afmode = 1;
    } else if (ev.compare("central") == 0) {
      t2_afmode = 2;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/time2_enth_vel = '" << ev
                << "' is not a choice (old | start | central)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // time2_vet_extrap (default false): the stage solves use D^n of the formal solution at
  // U^n; true = D* = D^n + (dt/dt_prev)(D^n - D^{n-1}) (runs_3x: fails G1 at P=100,
  // tau=10)
  t2_vext = false;
  if (pin->DoesParameterExist("rad_m1", "time2_vet_extrap")) {
    t2_vext = pin->GetBoolean("rad_m1", "time2_vet_extrap");
  }
  // time2_vet_sc (read only when named; default lag): when the vet_sc
  // tensor of the two stage solves is built.  Both stages sit at t^{n+1} (c = 1), so
  // D(U^n) is an O(dt) lag and the step is FIRST order in time wherever the tensor
  // changes (torder_1010: pulse_k12.8 T 1.05 1.03 1.02; eddington or the extrapolated
  // D* give 2.00).  The modes are those of time2_vet_col (Time2VetColAt):
  //   lag     : D(U^n) for both stages (time2_vet_extrap = true: D* extrapolated)
  //   predict : ONE formal solution per step at P = (stage-1 start) + dt K1, i.e. at
  //             t^{n+1} + O(dt^2); same cost as lag
  //   rebuild : predict for stage 1, a second formal solution at Y1 for stage 2
  // DEFAULT predict on a fresh run (as time2_vet_col on the wedge), unless
  // time2_vet_extrap = true or vet_sc_every > 1 (which predict excludes: lag, silently);
  // a restart whose file lacks the key keeps lag (echoed).  lag is first
  // order in time wherever the tensor changes (pulse_k12.8 T 1.05), predict 2.0.
  // ONE key time2_vet for both formal solutions (lag | predict |
  // rebuild | extrap); the old spellings time2_vet_sc (vet_sc) and time2_vet_col (vet_col
  // / vet_gd) are read as aliases of their own closure (M1KeyName), so an input or a
  // restart file that names only the old key is unchanged.  vet_sc: extrap = lag with
  // time2_vet_extrap = true (the D* extrapolation, as before).  vet_every likewise for
  // vet_sc_every / vet_col_every.
  t2_vsmode = 0;
  const std::string kvs = M1KeyName(pin, "time2_vet", {"time2_vet_sc"});
  const std::string kev = M1KeyName(pin, "vet_every", {"vet_sc_every"});
  const bool vs_named = pin->DoesParameterExist("rad_m1", kvs);
  const bool vs_every = pin->DoesParameterExist("rad_m1", kev) &&
                        (pin->GetInteger("rad_m1", kev) > 1);
  if (vet_sc && vs_named && pin->GetString("rad_m1", kvs).compare("extrap") == 0) {
    t2_vext = true;
  }
  if (vet_sc && (vs_named || (!t2_vext && !vs_every))) {
    std::string vm = pin->GetOrAddString("rad_m1", kvs,
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
                << std::endl << "<rad_m1>/" << kvs << " = '" << vm
                << "' is not a choice (lag | predict | rebuild | extrap)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    const bool every = vs_every;
    if (t2_vsmode != 0 && (t2_vext || every)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/" << kvs << " = " << vm << " excludes "
                << "time2_vet_extrap = true and " << kev << " > 1" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // time2_vet_col (closure = vet_col; m1-sp-order2b, tests_m1/runs_5q_sporder2b): when
  // the vet_col tensor (and surface q) of the two stage solves is built.  Both stages
  // sit at t^{n+1} (c = 1), so D(U^n) is an O(dt) lag and the step first order in time.
  //   lag     : D(U^n) for both stages (the m1-sph2 form)
  //   predict : ONE build per step at P = (stage-1 start) + dt K1 (radiation U^n + dt
  //             K1, gas the Heun predictor + dt K1 of the coupling): D(t^{n+1}) + O(dt^2)
  //   rebuild : predict for stage 1, then a second build at the stage-1 solution Y1
  //   extrap  : DIAGNOSTIC, D^n + (dt/dt_prev)(D^n - D^{n-1}); the history is not in
  //             the restart file
  // DEFAULT predict on the spherical-polar wedge; lag elsewhere and on a restart whose
  // file lacks the key (echoed).  K1 is restart state: predict restarts bitwise.
  // time2_vstage (m1-sp-order2b, tests_m1/runs_5q_sporder2b): in a stage solve under
  // implicit_vimp the lagged v-dependent terms take the stage's own velocity instead of
  // the stage-start one (which misses the stage's radiative kick, an O(dt) lag): the
  // comoving correction E0 - E at v_old + dv^k of the iterate, and the derived cell flux
  // F = F0 + a E at the velocity the write-back gives the gas.  Without it F is first
  // order in time (rsw_u_T: 1.00, E/T/gas 1.9-2.0).  Default true on the spherical-polar
  // wedge (a restart whose file lacks the key keeps false, echoed); elsewhere false
  // unless named.
  // default true on every mesh (was the wedge only): the Cartesian vet_sc
  // ladders need it as much (cenPvs: rw_t10a F1 1.58 2.32 1.97 vs 1.07 1.04
  // 1.02 without); a restart whose file lacks the key keeps false (echoed)
  t2_vcmode = 0;
  if (vet_col) {
    // (USER 10-25): the sph fresh-start default is extrap2 (was predict;
    // RSG E / HG production value); a restart whose file lacks the key keeps lag
    const bool vdef = sph_geom && !global_variable::restart_run;
    const std::string kvc = M1KeyName(pin, "time2_vet", {"time2_vet_col"});
    std::string vm = (sph_geom || pin->DoesParameterExist("rad_m1", kvc)) ?
        pin->GetOrAddString("rad_m1", kvc, vdef ? "extrap2" : "lag") :
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
      t2_vcmode = 4;   // (rad_m1.hpp Time2VetColExtrap2)
      x2_fields = pin->GetOrAddInteger("rad_m1", "time2_vet_x2_fields", 9);
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/" << kvc << " = '" << vm
                << "' is not a choice (lag | predict | rebuild | extrap | extrap2)"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (t2_vcmode != 0 && vcol_every != 1) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/" << kvc << " = " << vm
                << " needs vet_every (vet_col_every) = 1" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // time2_vet_source (was time2_vet_col_source): Time2VetSourceInit
    Time2VetSourceInit(pin, t2_vcmode == 1 || t2_vcmode == 2, true);
  }
  // vet_sc reads the same time2_vet_source (its builds take the
  // VET source rule of vet_col / vet_gd, rad_m1_vet.cpp).  A restart whose file lacks the
  // key keeps vet_sc's old source of the predicted state (state); the old spelling
  // time2_vet_col_source is vet_col's own (a vet_sc file may name it, inert there).
  if (vet_sc) {
    Time2VetSourceInit(pin, t2_vsmode == 1 || t2_vsmode == 2, false);
  }
  // time2_pred_cap: the predicted state of time2_vet_sc / time2_vet_col =
  // predict | rebuild is P = U^n + theta dt K1 with, per cell,
  //   theta = min(1, cap/r),  r = max(|dt K1_E|/E^n, |eint(P) - eint_b|/eint_b),
  // eint_b the gas internal energy of the stage start without the coupling's K1.  K1 is
  // the PREVIOUS step's slope: where the state relaxes within a step (a stiff exchange
  // far from equilibrium, a light front at c dt/dx >> 1) the linear extrapolation
  // overshoots (eq2: E^n + dt K1 = 2 E^n, T(P) = 0 from a negative eint), and the
  // formal solution of that unphysical P puts a closure and a half-range/thin-switch
  // flux into the solve that the state does not have.  theta caps the relative change of
  // the predictor at cap (P stays positive for cap < 1): where r <= cap (every cell of a
  // resolved evolution as dt -> 0, r = O(dt)) theta = 1 and the build is bitwise the
  // unlimited one, so the order is kept; where r > cap the build moves toward D(U^n),
  // the first-order lag, which is what the step can trust there.  0 = off.  Default
  // 0.5 on fresh runs; a restart whose file lacks the key keeps 0 (bitwise).
  t2_pcap = 0.0;
  if (t2_vsmode != 0 || t2_vcmode == 1 || t2_vcmode == 2) {
    t2_pcap = pin->GetOrAddReal("rad_m1", "time2_pred_cap",
                                global_variable::restart_run ? 0.0 : 0.5);
    if (!(t2_pcap >= 0.0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/time2_pred_cap must be >= 0" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // DEBUG dbg_t2_pred = N (read only when named): for the first N predict builds print
  // the predicted state (min E, cells with E <= 0 / eint <= 0, max r, theta), and for
  // vet_sc the formal solution's mean and max H/J and the sweep's clip count
  t2_dbg_pred = 0;
  t2_dbg_pred_n = 0;
  if (pin->DoesParameterExist("rad_m1", "dbg_t2_pred")) {
    t2_dbg_pred = pin->GetInteger("rad_m1", "dbg_t2_pred");
  }
  // time2_lin_tol / time2_lin_tol_fac (tests_m1/runs_5f_h2fast): the linear (Krylov)
  // tolerance of the two stage solves: time2_lin_tol when named, else implicit_lin_tol
  // x time2_lin_tol_fac.  Default 1 since audit-m1-0928 (was 10 on fresh runs from
  // tests_m1/runs_5f_h2fast: 3-D He box 21 -> 15 Krylov iterations per stage solve).
  // Factor 10 degrades the stiff radiation-wave time convergence
  // (/viper/ptmp2/jinma/audit_m1_0928/rw/RESULTS_time.txt: error 1.6-4x larger, order
  // 1.37 -> 0.48 at tau 10), for a gain of -0.4..+7.5 % per simulated second (He box,
  // sph_wedge, MI300A).  Opt-in by name.  A restart whose file records the key keeps its
  // value; one that does not gets 1 (as before).
  const bool rs = global_variable::restart_run;
  // time2_one_pass_safety (tests_m1/runs_5f_h2fast): the implicit_one_pass safety factor
  // of the stage solves.  One-pass acceptance leaves a Picard error ~ tol/safety where
  // the two-pass test leaves ~ q tol; the radwave G1 order at implicit_tol = 1e-11 needs
  // the stage solves closer to the latter (safety 3: median order 1.74 at P=100,
  // tau=1e3; 30: >= 1.94 in all 24 cases).  Default 30; 0 = implicit_one_pass_safety,
  // the default of a restart whose file does not carry the key.
  t2_onep_s = pin->GetOrAddReal("rad_m1", "time2_one_pass_safety", rs ? 0.0 : 30.0);
  if (t2_onep_s != 0.0 && !(t2_onep_s >= 1.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<rad_m1>/time2_one_pass_safety must be 0 or >= 1"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // time2_start / time2_start_frac (rad_m1.hpp; read only when named):
  // be (default, bitwise) | slope.  The Lie-split backward-Euler first step leaves an
  // O(dt) start-up error that the second-order steps after it carry along as a free
  // mode (h2order_1011: P = 10 linear wave, tau 100, v 1.1 1.2 1.0 -> see NOTE).
  if (pin->DoesParameterExist("rad_m1", "time2_start")) {
    std::string st = pin->GetString("rad_m1", "time2_start");
    if (st.compare("be") == 0) {
      t2_start = 0;
    } else if (st.compare("slope") == 0) {
      t2_start = 1;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/time2_start = '" << st
                << "' is not a choice (be | slope)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  if (pin->DoesParameterExist("rad_m1", "time2_start_frac")) {
    t2_start_frac = pin->GetReal("rad_m1", "time2_start_frac");
    if (!(t2_start_frac > 0.0 && t2_start_frac <= 1.0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/time2_start_frac must be in (0, 1]"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  t2_dbg_fail = -1;
  if (pin->DoesParameterExist("rad_m1", "time2_dbg_fail")) {
    t2_dbg_fail = pin->GetInteger("rad_m1", "time2_dbg_fail");
  }
  if (pin->DoesParameterExist("rad_m1", "dbg_nonconv")) {
    dbg_nonconv = pin->GetInteger("rad_m1", "dbg_nonconv");
  }
  // dbg_t2_admiss = N (DEBUG, default 0 = off, read only when named so the restart echo
  // and every result stay bitwise): for the first N non-admissible hesdirk2 stages print
  // which quantity went <= 0 (solved E, solved T, written-back gas eint) and, for each,
  // the cell of its minimum with the stage OLD vector (E, eint) and the stage start E
  // time2_inc_energy = gas | rad (rad_m1.hpp; read only when named)
  if (pin->DoesParameterExist("rad_m1", "time2_inc_energy")) {
    const std::string ie = pin->GetString("rad_m1", "time2_inc_energy");
    if (ie != "gas" && ie != "rad") {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/time2_inc_energy = '" << ie
                << "' (gas | rad)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    t2_ike = (ie == "rad");
    if (global_variable::my_rank == 0 && t2_ike) {
      std::cout << "<rad_m1> time2_inc_energy=rad: the stage old vector advances the gas "
                << "eint linearly (slope K_EN - K_KE); KE follows the momentum and the "
                << "nonlinear KE remainder is exchanged with the cell's E (total energy "
                << "exact, no cap)" << std::endl;
    }
  }
  // time2_exch_split = none | psi (rad_m1.hpp; read only when named)
  if (pin->DoesParameterExist("rad_m1", "time2_exch_split")) {
    const std::string xs = pin->GetString("rad_m1", "time2_exch_split");
    if (xs != "none" && xs != "psi") {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/time2_exch_split = '" << xs
                << "' (none | psi)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    t2_xsp = (xs == "psi");
    if (global_variable::my_rank == 0 && t2_xsp) {
      std::cout << "<rad_m1> time2_exch_split=psi: the carried gas-radiation exchange "
                << "terms c K_X take the weight 1/sqrt(1 + (lambda c)^2) (lambda = the "
                << "stage solve's linearised exchange rate); E and gas EN opposite, "
                << "energy exact" << std::endl;
    }
  }
  // time2_inc_exch = none | relax (rad_m1.hpp; read only when named)
  if (pin->DoesParameterExist("rad_m1", "time2_inc_exch")) {
    const std::string ie = pin->GetString("rad_m1", "time2_inc_exch");
    if (ie != "none" && ie != "relax") {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/time2_inc_exch = '" << ie
                << "' (none | relax)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    t2_iex = (ie == "relax");
    if (pin->DoesParameterExist("rad_m1", "time2_inc_exch_y0")) {
      t2_iex_y0 = pin->GetReal("rad_m1", "time2_inc_exch_y0");
    }
    if (global_variable::my_rank == 0 && t2_iex) {
      std::cout << "<rad_m1> time2_inc_exch=relax (y0=" << t2_iex_y0 << "): a falling "
                << "old-vector eint takes the relaxation form eint_h expm1(d/eint_h), "
                << "blended in by y^2/(y^2+y0^2); E pays the difference" << std::endl;
    }
  }
  // <rad_m1>/time2_reject = be | retry (read only when named; default be =
  // the backward-Euler redo, bitwise).  retry needs <time>/dt_retry = true (checked at
  // the first rejection: without it the step falls back to backward Euler as before)
  t2_reject = 0;
  if (pin->DoesParameterExist("rad_m1", "time2_reject")) {
    const std::string tr = pin->GetString("rad_m1", "time2_reject");
    if (tr == "retry") {
      t2_reject = 1;
    } else if (tr != "be") {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/time2_reject = '" << tr << "' is not a choice "
                << "(be | retry)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // <rad_m1>/time2_k1_reset (see rad_m1.hpp)
  t2_k1reset = false;
  rst_ext = true;   // USER 11-07 make some defaults: <rad_m1>/rst_state_ext default true
  if (pin->DoesParameterExist("rad_m1", "rst_state_ext")) {
    rst_ext = pin->GetBoolean("rad_m1", "rst_state_ext");
  }
  if (pin->DoesParameterExist("rad_m1", "dbg_rst_fp")) {
    dbg_rst_fp = pin->GetInteger("rad_m1", "dbg_rst_fp");
  }
  // dbg_rst_fp_cycle = C (read only when named): the first clean solve of cycle C also
  // prints fingerprints at its setup and around each linear solve; dbg_rst_fp_dump: every
  // fingerprint also writes vgd_cs to vgdcs_<tag>_r<rank>.bin (5 int extents, the data)
  if (pin->DoesParameterExist("rad_m1", "dbg_rst_fp_cycle")) {
    dbg_rst_fp_cycle = pin->GetInteger("rad_m1", "dbg_rst_fp_cycle");
  }
  if (pin->DoesParameterExist("rad_m1", "dbg_rst_fp_dump")) {
    dbg_rst_fp_dump = pin->GetBoolean("rad_m1", "dbg_rst_fp_dump");
  }
  if (pin->DoesParameterExist("rad_m1", "time2_k1_reset")) {
    const std::string kr = pin->GetString("rad_m1", "time2_k1_reset");
    if (kr == "pos") {
      t2_k1reset = true;
      t2_k1pos = true;
    } else {
      t2_k1reset = pin->GetBoolean("rad_m1", "time2_k1_reset");
    }
  }
  if (pin->DoesParameterExist("rad_m1", "time2_k1_xlim")) {
    t2_k1xlim = pin->GetReal("rad_m1", "time2_k1_xlim");
  }
  // <rad_m1>/time2_k1_xpat (default false, read only when named): Time2K1XLim's theta in
  // Patankar form, positive by construction with no threshold (rad_m1.hpp); = pos: only
  // where the unlimited loss exceeds 0.99 of the loser's holding (Time2K1XLim)
  if (pin->DoesParameterExist("rad_m1", "time2_k1_xpat")) {
    t2_k1xpos = pin->GetString("rad_m1", "time2_k1_xpat") == "pos";
    t2_k1xpat = t2_k1xpos || pin->GetBoolean("rad_m1", "time2_k1_xpat");
  }
  // <rad_m1>/time2_k1_xpat_work (default false, read only when named; needs xpat = pos):
  // E's loss L of the E side is the whole local pair, the exchange and the work that
  // pays the gas KE change, (ch/c)(de + w); the part taken back is the gas eint's
  // (momentum unchanged), at most what keeps 0.01 of its eint
  if (pin->DoesParameterExist("rad_m1", "time2_k1_xpat_work")) {
    t2_k1xwork = pin->GetBoolean("rad_m1", "time2_k1_xpat_work");
    if (t2_k1xwork && !t2_k1xpos) {
      std::cout << "### FATAL ERROR <rad_m1>/time2_k1_xpat_work needs time2_k1_xpat = pos"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // <rad_m1>/time2_k1_stiff = s (default 0 = off, bitwise)
  t2_k1stiff = 0.0;
  if (pin->DoesParameterExist("rad_m1", "time2_k1_stiff")) {
    t2_k1stiff = pin->GetReal("rad_m1", "time2_k1_stiff");
  }
  // (USER 11-08 "make some defaults", r8 key set): DEFAULT 0.1 where it is
  // legal (one radial MeshBlock; column-local), else 0 = the old default, logged; recipe
  // none keeps 0.  Read only when named as before.
  // (USER 11-13 "make some defaults", the R136/HG/RSG key set): where
  // time2_face_fct is on by default (cubed sphere / spherical polar, multi-D implicit
  // transport, hesdirk2; below) the face limiter replaces the redistribution and the
  // DEFAULT is 0; elsewhere (Cartesian, 1-D) the prod-1108 rule above stays.  The clean
  // scheme does not redistribute (a column-local compensation; the answer would depend
  // on the MeshBlocks along x1).
  const bool ffct14 = recipe::GetMode(pin) != recipe::Mode::none && sph_geom &&
                      trans_on && time_scheme == M1_TIME_HESDIRK2;
  {
    const bool onecol = pmy_pack->pmesh->mesh_indcs.nx1 == pmy_pack->pmesh->mb_indcs.nx1;
    const bool rnone = recipe::GetMode(pin) == recipe::Mode::none;
    t2_k1rd = (onecol && !rnone && !ffct14 && !clean) ? 0.1 : 0.0;
    if (!pin->DoesParameterExist("rad_m1", "time2_k1_redist")) {
      recipe::Log(std::string("<rad_m1>/time2_k1_redist = ") + (t2_k1rd > 0.0 ? "0.1" :
                  "0") + (ffct14 ? " (default, defaults-1114: time2_face_fct)" :
                  std::string(" (default, prod-1108") + (onecol ? (rnone ?
                  "; recipe none)" : ")") : "; 0: more than one radial MeshBlock)")));
    }
  }
  if (pin->DoesParameterExist("rad_m1", "time2_k1_redist")) {   // seamE-1106
    t2_k1rd = pin->GetReal("rad_m1", "time2_k1_redist");
  }
  if (t2_k1rd > 0.0 && pmy_pack->pmesh->mesh_indcs.nx1 != pmy_pack->pmesh->mb_indcs.nx1) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<rad_m1>/time2_k1_redist is column-local: one radial block only"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // step 1 <rad_m1>/time2_face_fct (rad_m1.hpp; read only when named)
  // (USER 11-13): DEFAULT true on the cubed sphere / spherical polar with
  // multi-D implicit transport and hesdirk2 (ffct14 above); false elsewhere and with
  // recipe none, as before.  An explicit key wins.  The clean scheme: on for a fresh
  // run (v20.4, recorded; a restart file without the key keeps off); its t2inc F
  // channels are face F as here, and its solve stores the energy flux Phi that moved E
  // (CleanWriteBack)
  const bool cfd = clean && !global_variable::restart_run;
  t2_ffct = (ffct14 && !clean) || cfd;
  if (pin->DoesParameterExist("rad_m1", "time2_face_fct") || cfd) {
    t2_ffct = pin->GetOrAddBoolean("rad_m1", "time2_face_fct", t2_ffct);
  }
  t2_pfst = t2_ffct || clean;
  t2_ffct_emin = clean ? 0.1 : 0.0;   // (rad_m1.hpp)
  if (pin->DoesParameterExist("rad_m1", "time2_face_fct_emin")) {
    t2_ffct_emin = pin->GetReal("rad_m1", "time2_face_fct_emin");
  }
  if (pin->DoesParameterExist("rad_m1", "time2_face_fct_sweep")) {
    t2_ffct_swp = pin->GetBoolean("rad_m1", "time2_face_fct_sweep");
  }
  if (pin->DoesParameterExist("rad_m1", "time2_face_fct_rounds")) {
    t2_ffct_rk = pin->GetInteger("rad_m1", "time2_face_fct_rounds");
  }
  // (the one-pass donor form was removed: it fails the accuracy gate)
  if (pin->DoesParameterExist("rad_m1", "time2_face_fct_form") &&
      pin->GetString("rad_m1", "time2_face_fct_form") != "zalesak") {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<rad_m1>/time2_face_fct_form: only zalesak (the default) "
              << "remains; the donor form was removed" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (pin->DoesParameterExist("rad_m1", "time2_k1_redist_r")) {
    t2_k1rd_r = pin->GetInteger("rad_m1", "time2_k1_redist_r");
  }
  if (pin->DoesParameterExist("rad_m1", "time2_k1_pos_after_stiff")) {   // seamE-1106
    t2_k1pos_as = pin->GetBoolean("rad_m1", "time2_k1_pos_after_stiff");
  }
  if (t2_k1pos && pin->DoesParameterExist("rad_m1", "time2_k1_pos_theta")) {
    t2_k1pos_th = pin->GetReal("rad_m1", "time2_k1_pos_theta");
  }
  // time2_k1_reset_grow = 0 | 1 (read only when named): the face-neighbour layer of the
  // reset cells also drops its K1 (cell and face channels) in the stage-1 solve
  t2_k1grow = 0;
  if (pin->DoesParameterExist("rad_m1", "time2_k1_reset_grow")) {
    t2_k1grow = pin->GetInteger("rad_m1", "time2_k1_reset_grow");
    if (t2_k1grow < 0 || t2_k1grow > 1) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/time2_k1_reset_grow must be 0 or 1"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  t2_ys_ok = false;
  // time2_accept_fallback is DEFAULT true under hesdirk2 (this
  // reader returns above for be, where the key is unused); false only when named so
  t2_accfb = true;
  if (pin->DoesParameterExist("rad_m1", "time2_accept_fallback")) {
    t2_accfb = pin->GetBoolean("rad_m1", "time2_accept_fallback");
  }
  // <rad_m1>/time2_floor_eint (rad_m1.hpp; read only when named)
  if (pin->DoesParameterExist("rad_m1", "time2_floor_eint")) {
    t2_flr_rho = pin->GetReal("rad_m1", "time2_floor_eint");
  }
  // <rad_m1>/time2_adm_poslim (rad_m1.hpp; read only when named)
  if (pin->DoesParameterExist("rad_m1", "time2_adm_poslim")) {
    t2_adm_pl = pin->GetBoolean("rad_m1", "time2_adm_poslim");
  }
  // <rad_m1>/time2_auto_be = N (rad_m1.hpp; read only when named)
  if (pin->DoesParameterExist("rad_m1", "time2_auto_be")) {
    t2_abe = pin->GetInteger("rad_m1", "time2_auto_be");
    if (t2_abe < 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/time2_auto_be must be >= 0" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  if (global_variable::my_rank == 0 && t2_abe > 0) {
    std::cout << "<rad_m1> time2_auto_be=" << t2_abe << ": after a stage fallback the "
              << "next 1, 2, 4, ... " << t2_abe << " steps take the backward-Euler redo "
              << "path at once" << std::endl;
  }
  if (pin->DoesParameterExist("rad_m1", "restart_pack_per_mb")) {
    rst_pack_mb = pin->GetBoolean("rad_m1", "restart_pack_per_mb");
  }

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int nmb = std::max((pmy_pack->nmb_thispack), (pmy_pack->pmesh->nmb_maxperrank));
  int n1 = indcs.nx1 + 2*(indcs.ng);
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  if (t2_k1reset) {
    Kokkos::realloc(t2ys, nmb, n3, n2, n1);
    Kokkos::deep_copy(t2ys, 0.0);
  }
  Kokkos::realloc(t2k1, nmb, M1_T2_NK, n3, n2, n1);
  Kokkos::realloc(t2k2, nmb, M1_T2_NK, n3, n2, n1);
  Kokkos::realloc(t2inc, nmb, M1_T2_NK, n3, n2, n1);
  Kokkos::deep_copy(t2k1, 0.0);
  Kokkos::deep_copy(t2k2, 0.0);
  Kokkos::deep_copy(t2inc, 0.0);
  if (t2_xsp) {   // time2_exch_split
    Kokkos::realloc(t2kx1, nmb, n3, n2, n1);
    Kokkos::realloc(t2kx2, nmb, n3, n2, n1);
    Kokkos::realloc(t2lx1, nmb, n3, n2, n1);
    Kokkos::realloc(t2lx2, nmb, n3, n2, n1);
    Kokkos::deep_copy(t2kx1, 0.0);
    Kokkos::deep_copy(t2kx2, 0.0);
    Kokkos::deep_copy(t2lx1, 0.0);
    Kokkos::deep_copy(t2lx2, 0.0);
  }
  if (t2_iex) {   // time2_inc_exch
    Kokkos::realloc(t2ixd, nmb, n3, n2, n1);
    Kokkos::deep_copy(t2ixd, 0.0);
  }
  if (t2_pfst) {   // step 1: zero = no face flux known (restart, start)
    Kokkos::realloc(t2pf1, nmb, 3, n3, n2, n1);
    Kokkos::realloc(t2pf2, nmb, 3, n3, n2, n1);
    Kokkos::deep_copy(t2pf1, 0.0);
    Kokkos::deep_copy(t2pf2, 0.0);
    t2pf1_ok = false;
    t2pf2_ok = false;
  }
  if (t2_ike) {
    Kokkos::realloc(t2ke1, nmb, n3, n2, n1);
    Kokkos::realloc(t2ke2, nmb, n3, n2, n1);
    Kokkos::realloc(t2kei, nmb, n3, n2, n1);
    Kokkos::realloc(t2ikd, nmb, n3, n2, n1);
    Kokkos::realloc(t2mo, nmb, 3, n3, n2, n1);
    Kokkos::deep_copy(t2ke1, 0.0);
    Kokkos::deep_copy(t2ke2, 0.0);
    Kokkos::deep_copy(t2kei, 0.0);
    Kokkos::deep_copy(t2ikd, 0.0);
    Kokkos::deep_copy(t2mo, 0.0);
  }
  Kokkos::realloc(t2f1, nmb, n3, n2, n1+1);
  if (trans_on) {
    Kokkos::realloc(t2f2, nmb, n3, n2+1, n1);
    if (trans_x3) {Kokkos::realloc(t2f3, nmb, n3+1, n2, n1);}
  }
  const bool noh = clean && cln_not2h;   // clean_skip_time2_hist
  if (impl_pred && !noh) {
    Kokkos::realloc(ipred2, nmb, (impl_pord == 2) ? 5 : 3, n3, n2, n1);
    Kokkos::deep_copy(ipred2, 0.0);
  }
  if (vet_sc) {
    Kokkos::realloc(vet_prev, nmb, M1_T2_NVET, n3, n2, n1);
    Kokkos::realloc(vet_now, nmb, M1_T2_NVET, n3, n2, n1);
    Kokkos::realloc(vet_opac, nmb, vscat ? (M1_NOPAC + 1) : M1_NOPAC, n3, n2,
                    n1);
    Kokkos::deep_copy(vet_prev, 0.0);
  }
  if (vet_col) {
    // Time2VetStart's save slots (E, T) and the saved stage-start opacities;
    // time2_vet_col = predict | rebuild: also KT, T(Y1), E(Y1) and the inner face flux
    // (scheme = clean: only ImplicitSolve's stage chain reads them; not made)
    if (!clean) {
      Kokkos::realloc(vet_now, nmb, (t2_vcmode == 0 || t2_vcmode >= 3) ? 2 : 6, n3, n2,
                      n1);
      Kokkos::realloc(vet_opac, nmb, vscat ? (M1_NOPAC + 1) : M1_NOPAC, n3, n2, n1);
    }
    if (t2_vcmode == 4) {
      if (!noh) {
        Kokkos::realloc(x2_tt, nmb, M1_X2_NTT, n3, n2, n1);
        Kokkos::realloc(x2_hr, nmb, M1_X2_NHR, n3, n2, n1);
        Kokkos::realloc(x2_q, nmb, n3, n2);
        Kokkos::deep_copy(x2_tt, 0.0);
        Kokkos::deep_copy(x2_hr, 0.0);
        Kokkos::deep_copy(x2_q, 0.0);
      }
      x2_ok = false;
      x2_t = 0.0;
      x2_alpha = -2.0;
    }
    if (t2_vcmode == 3) {
      Kokkos::realloc(vcol_prev, nmb, n3, n2, n1);
      Kokkos::realloc(vcol_qprev, nmb, n3, n2);
      t2_vcprev = false;
    }
  }
  if (t2_pcap > 0.0 || t2_dbg_pred > 0) {
    Kokkos::realloc(t2pth, nmb, n3, n2, n1);
    Kokkos::deep_copy(t2pth, 1.0);
  }
  if (t2_dbg_pred > 0) {
    Kokkos::realloc(vet_dbgclip, 1);
    Kokkos::deep_copy(vet_dbgclip, 0.0);
  }
  t2_ok = false;
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1> time_scheme=hesdirk2 (g=" << kT2G << "; first step, restarts "
              << "without the slope and failed stages take backward Euler)" << std::endl;
    std::cout << "         time2_enth_be=" << (t2_enth_be ? "true" : "false")
              << " implicit_cs2_vimp_npass=" << cs2_vnp << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationM1::Time2Active
//! \brief true when the next step runs the two stage solves

bool RadiationM1::Time2Active() {
  if (time_scheme != M1_TIME_HESDIRK2) return false;
  if (static_cast<int>(t2k1.extent(0)) < pmy_pack->nmb_thispack) {t2_ok = false;}
  return t2_ok;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::Time2FormStage
//! \brief called after the hydro's "after_stagen" of Heun stage `stage`: builds the stage
//! START state of the radiation (stage 1: U^n; stage 2: the Heun average of U^n and Y2,
//! which the hydro has already done for the gas) and the increment t2inc that turns it
//! into the old vector of the stage solve.  ImplicitSolve applies t2inc.

//----------------------------------------------------------------------------------------
//! \fn Real RadiationM1::InnerFluxAt
//! \brief the inner face flux of the time hooks at time t: the table value
//! (ftab_on) or the constant F0, times the sine / tanh factors when
//! named (the arithmetic of the raycheap2-1109 hooks, unchanged)

Real RadiationM1::InnerFluxAt(const Real t) const {
  Real f = ftab_on ? ftab_scale*ftab.At(0, t - ftab_t0) : fmod_f0;
  if (fmod_per > 0.0) {f *= 1.0 + fmod_amp*std::sin(2.0*M_PI*t/fmod_per);}
  if (fmod_rw > 0.0) {f *= 0.5*(1.0 + std::tanh((t - fmod_rt0)/fmod_rw));}
  return f;
}

//! \brief the exact integral of the TABLE flux over [t0, t1] (the
//! energy budget's int F dt; table only, no sine / tanh factors)
Real RadiationM1::InnerFluxIntegral(const Real t0, const Real t1) const {
  return ftab_scale*ftab.Integral(0, t0 - ftab_t0, t1 - ftab_t0);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::Time2GasTime
//! \brief time2_tableau = sdirk2 with the Heun hydro (c_E = (0, 1), b_E = (1/2, 1/2)):
//! the IMEX coupling conditions b_E.c_I = b_I.c_E = 1/2 hold when each stage has one
//! time for both parts (c = (0, g, 1, 1), the hydro's stage-2 flux an explicit stage).
//! which = 1, before the stage-1 solve: the gas U* of Heun stage 1 (t^{n+1}) is kept in
//! t2hs and the solve takes S = U^n + g (U* - U^n), the gas at t^n + g dt (convex, so
//! positive; all cells).  which = 2, first in its HydroConToPrim (gas Y = S + d): the
//! gas becomes U* + d (the explicit U* with the implicit g dt K2, so the stage-2
//! increment and the exchange, booked once in d, are unchanged), before the ring average,
//! the floor cap and the refresh act on it as they did on Y.  which = 3, after the stage
//! tasks: the hydro's stage-2 flux state U* + f d (t^{n+1}, the exchange over the full
//! dt: f = 1/g non-stiff, -> 1 stiff, below) is refreshed into w0 (its floors stay in
//! u0), then u0 is U* + d again.

void RadiationM1::Time2GasTime(Driver *pdrive, int which) {
  hydro::Hydro *ph = pmy_pack->phydro;
  if (t2_tab != 2 || !coupling || ph == nullptr || (which > 1 && t2_gt != which - 1)) return;
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
  // the flux state's factor on d: 1/g + O(k dt) (the ARK conditions), 1 where the
  // exchange is stiff (k dt >> 1: the backward-Euler ratio of a full dt to g dt of a
  // linear relaxation at rate k = ch rho kappa_E + 4 cl rho kappa_P a T^3/(rho c_v))
  auto op_ = opac;
  auto eos = ph->peos->eos_data;
  auto phc = ph->phicc0;
  auto ccl = pmy_pack->pcoord->cos_cell;
  const bool etg = ph->use_etotgrav, csk = cs_geom;
  const Real g = kT2G, cl = c_light, ch = chat, ar = arad, dt = pmy_pack->pmesh->dt;
  par_for("m1_t2_gtime", DevExeSpace(), 0, nmb-1, kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real e = 0.0;
    if (which == 2) {
      const Real dd = u0_(m,IDN,k,j,i);
      const Real eg = u0_(m,IEN,k,j,i) - CsKeMom(u0_(m,IM1,k,j,i), u0_(m,IM2,k,j,i),
                      u0_(m,IM3,k,j,i), csk ? ccl(m,k,j) : 0.0)/dd -
                      (etg ? dd*phc(m,k,j,i) : 0.0);
      const Real tk = (eg > 0.0) ? eos.Temperature(dd, eg) : 0.0;
      const Real z = dt*(ch*op_(m,M1_OP_E,k,j,i) + ((eg > 0.0) ?
                     4.0*cl*op_(m,M1_OP_P,k,j,i)*ar*SQR(SQR(tk))/eg : 0.0));
      e = (1.0 + g*z)/(g*(1.0 + z)) - 1.0;
    }
    for (int n = 0; n < nv; ++n) {
      const Real un = u1_(m,n,k,j,i), us = hs_(m,n,k,j,i), u = u0_(m,n,k,j,i);
      if (which == 1) {
        hs_(m,n,k,j,i) = u;
        u0_(m,n,k,j,i) = un + g*(u - un);
      } else if (which == 2) {
        const Real d = u - (un + g*(us - un));   // (0 where the solve left it)
        u0_(m,n,k,j,i) = us + d;
        hs_(m,n,k,j,i) = e*d;
      } else {
        u0_(m,n,k,j,i) = u + us;
      }
    }
  });
  if (which != 3) return;
  FluidRefresh(pdrive, false);   // ghosts + w0 of the flux state
  par_for("m1_t2_gtime3", DevExeSpace(), 0, nmb-1, kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int n = 0; n < nv; ++n) {u0_(m,n,k,j,i) -= hs_(m,n,k,j,i);}
  });
}

void RadiationM1::Time2FormStage(int stage, Real dt) {
  if (fmod_per > 0.0 || fmod_rw > 0.0 || ftab_on) {   // inner flux
    // hesdirk2 / sdirk2: both stage solves end at t_n + dt (c2 = a21 + g = 1, c3 = 1);
    // trbdf2: the stage-1 solve ends at t_n + 2 g dt (c2 = 2g)
    const Real cst = (stage == 1 && t2_tab == 1) ? 2.0*kT2G : 1.0;
    const Real ts = pmy_pack->pmesh->time + cst*dt;
    iflux_x1min = InnerFluxAt(ts);
  }
  if (dbg_rst_fp > 0 && stage == 1) {dbg_rst_fp -= 1; RstFingerprint("formstage1");}
  const Real g = kT2G;
  const Real w0 = t2_w0, w1 = t2_w1;
  // hesdirk2 and sdirk2 both use the Heun stage start (tab 0's original expression,
  // bitwise); only trbdf2 (tab 1) takes the w0/w1 average
  const bool heun = (t2_tab == 0 || t2_tab == 2);
  Real c1, c2;
  auto u0_ = u0;
  auto u1_ = u1;
  const int nf1 = static_cast<int>(f0x1.extent(3));
  if (stage == 1) {
    Kokkos::deep_copy(DevExeSpace(), u1, u0);
    Kokkos::deep_copy(DevExeSpace(), t2f1, f0x1);
    if (impl_hr_cmem && fax1.extent_int(0) > 0) {   // (booking only)
      if (t2fa1.extent(0) != fax1.extent(0) || t2fa1.extent(3) != fax1.extent(3)) {
        Kokkos::realloc(t2fa1, fax1.extent(0), fax1.extent(1), fax1.extent(2),
                        fax1.extent(3));
      }
      Kokkos::deep_copy(DevExeSpace(), t2fa1, fax1);
    }
    if (trans_on) {
      Kokkos::deep_copy(DevExeSpace(), t2f2, f0x2);
      if (trans_x3) {Kokkos::deep_copy(DevExeSpace(), t2f3, f0x3);}
    }
    c1 = t2_a21*dt;
    c2 = 0.0;
    // time2_k1_reset: no K1 where E was overwritten since the stage-loop start
    if (t2_k1reset && t2_ys_ok && t2ys.extent_int(0) >= pmy_pack->nmb_thispack) {
      auto &id = pmy_pack->pmesh->mb_indcs;
      const int is = id.is, ie = id.ie, js = id.js, je = id.je, ks = id.ks, ke = id.ke;
      const int nmbr = pmy_pack->nmb_thispack - 1;
      auto ys_ = t2ys;
      auto k1r_ = t2k1;
      const int ncell = M1_T2_F1;   // E, gas momentum (3), gas energy
      // Coordinates::volume is allocated (and filled) only on the cubed
      // sphere / spherical polar; elsewhere it is 1x1x1x1 and vol(m,k,j,i) read out of
      // bounds.  Cartesian: the uniform cell volume of the MeshBlock.
      auto &vol = pmy_pack->pcoord->volume;
      const bool usev = (vol.extent_int(0) >= nmbr + 1) && (vol.extent_int(3) > 1);
      auto &mbsz = pmy_pack->pmb->mb_size;
      Real nr = 0.0, de = 0.0;
      // time2_k1_reset = pos: the stage-1 old vector of E and of the gas eint
      if (t2_k1pos) {
        FluidRef flp = FluidRef::Get(pmy_pack);
        const bool hg = flp.on;
        auto ug = hg ? flp.u0 : u0;
        const bool etg = hg && flp.use_etotgrav;
        auto phc = hg ? flp.phicc0 : DvceArray4D<Real>();
        const bool emh = hg && flp.mhd;
        auto emg = emag0;
        const Real th = t2_k1pos_th, ca = c1;
        Real np = 0.0;
        Kokkos::parallel_reduce("m1_t2_k1pos", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
            DevExeSpace(), {0, ks, js, is}, {nmbr + 1, ke + 1, je + 1, ie + 1}),
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &sum) {
          if (u0_(m,M1_E,k,j,i) != ys_(m,k,j,i)) return;   // reset by the E rule anyway
          const Real e0 = u0_(m,M1_E,k,j,i);
          bool bad = !(e0 + ca*k1r_(m,M1_T2_E,k,j,i) > th*e0);
          if (hg && !bad) {
            const Real dd = ug(m,IDN,k,j,i), idd = 1.0/fmax(dd, 1.0e-300);
            Real eg0 = ug(m,IEN,k,j,i) - 0.5*(SQR(ug(m,IM1,k,j,i)) +
                                              SQR(ug(m,IM2,k,j,i)) +
                                              SQR(ug(m,IM3,k,j,i)))*idd;
            const Real m1 = ug(m,IM1,k,j,i) + ca*k1r_(m,M1_T2_M1,k,j,i);
            const Real m2 = ug(m,IM2,k,j,i) + ca*k1r_(m,M1_T2_M1+1,k,j,i);
            const Real m3 = ug(m,IM3,k,j,i) + ca*k1r_(m,M1_T2_M1+2,k,j,i);
            Real eg1 = ug(m,IEN,k,j,i) + ca*k1r_(m,M1_T2_EN,k,j,i)
                       - 0.5*(m1*m1 + m2*m2 + m3*m3)*idd;
            Real egr = etg ? dd*phc(m,k,j,i) : 0.0;
            if (emh) egr += emg(m,k,j,i);
            eg0 -= egr;
            eg1 -= egr;
            bad = (eg0 > 0.0) && !(eg1 > th*eg0);
          }
          if (bad) {
            ys_(m,k,j,i) = -1.0;   // marks the cell for the E rule below and the faces
            sum += 1.0;
          }
        }, Kokkos::Sum<Real>(np));
        t2_k1pos_n += np;
      }
      Kokkos::parallel_reduce("m1_t2_k1reset", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
          DevExeSpace(), {0, ks, js, is}, {nmbr + 1, ke + 1, je + 1, ie + 1}),
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &sum,
                    Real &sde) {
        if (u0_(m,M1_E,k,j,i) != ys_(m,k,j,i)) {
          const Real cv = usev ? vol(m,k,j,i) : (mbsz.d_view(m).dx1*mbsz.d_view(m).dx2*
                                                 mbsz.d_view(m).dx3);
          sde += (k1r_(m,M1_T2_E,k,j,i) + k1r_(m,M1_T2_EN,k,j,i))*cv;
          for (int n=0; n<ncell; ++n) {k1r_(m,n,k,j,i) = 0.0;}
          sum += 1.0;
        }
      }, Kokkos::Sum<Real>(nr), Kokkos::Sum<Real>(de));
      t2_k1reset_n += nr;
      t2_k1reset_de = de;
      // the face slopes of the faces of a reset cell are dropped by the stage-1 solve,
      // after a halo exchange of the reset mask (both blocks of a shared face agree)
      t2_ys_face = true;
    } else {
      t2_k1reset_de = 0.0;
      t2_ys_face = false;
    }
    // time2_k1_stiff = s (> 0): in cells whose gas-radiation relaxation
    // rate k = ch rho kappa_E + 4 cl rho kappa_P a T^3/(rho c_v) (rho c_v ~ eint/T) has
    // k dt > s, the EXCHANGE part of the FSAL slope K1 is moved from the gas energy to
    // the radiation E (ch/cl weighted, so the solve's total energy is unchanged): the
    // gas K1_EN keeps only the kinetic-energy change of its momentum slope (eint of the
    // old vector = eint of the stage start), K1_E keeps transport + the moved exchange.
    // The stiff exchange is then taken by the implicit stages only (L-stable there);
    // the E/gas partition of the old vector is forgotten at O(1/(k dt)).  Transport,
    // momentum and face slopes are untouched.  The stored K1 is changed (stage 2 too).
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
        Kokkos::parallel_reduce("m1_t2_k1stiff", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
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
          const Real rkp = op_(m,M1_OP_P,k,j,i), rke = op_(m,M1_OP_E,k,j,i);
          const Real rate = ch*rke + 4.0*cl*rkp*ar*tk*tk*tk*tk/eg;
          if (!(rate*dtl > sth)) return;
          const Real n1 = m1 + ca*k1s_(m,M1_T2_M1,k,j,i);
          const Real n2 = m2 + ca*k1s_(m,M1_T2_M1+1,k,j,i);
          const Real n3 = m3 + ca*k1s_(m,M1_T2_M1+2,k,j,i);
          const Real ek1 = CsKeMom(n1, n2, n3, cc)*idd;
          const Real kw = (ca != 0.0) ? (ek1 - ek0)/ca : 0.0;
          const Real dex = k1s_(m,M1_T2_EN,k,j,i) - kw;
          k1s_(m,M1_T2_EN,k,j,i) = kw;
          k1s_(m,M1_T2_E,k,j,i) += (ch/cl)*dex;
          sum += 1.0;
        }, Kokkos::Sum<Real>(ns));
        t2_k1stiff_n += ns;
        // seamE-1106: time2_k1_reset = pos tested K1_E before the exchange moved into it
        // above, so a cell the move takes to E + ca K1_E <= th E was missed (measured:
        // E 1990, increment -6100).  The same test again, same treatment.  Only with
        // time2_k1_reset = pos, time2_k1_stiff and key time2_k1_pos_after_stiff.
        if (t2_k1pos_as && t2_k1pos && t2_k1reset && t2_ys_ok &&
            t2ys.extent_int(0) >= pmy_pack->nmb_thispack) {
          auto ys2_ = t2ys;
          auto k1p_ = t2k1;
          auto ue_ = u0;
          const Real th2 = t2_k1pos_th;
          const int nc2 = M1_T2_F1;
          Real np2 = 0.0;
          Kokkos::parallel_reduce("m1_t2_k1pos2", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
              DevExeSpace(), {0, ks, js, is}, {nmbr + 1, ke + 1, je + 1, ie + 1}),
          KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &sum) {
            const Real e0 = ue_(m,M1_E,k,j,i);
            if (!(e0 + ca*k1p_(m,M1_T2_E,k,j,i) > th2*e0)) {
              for (int n = 0; n < nc2; ++n) {k1p_(m,n,k,j,i) = 0.0;}
              ys2_(m,k,j,i) = -1.0;   // the faces of the cell, as the first test
              sum += 1.0;
            }
          }, Kokkos::Sum<Real>(np2));
          t2_k1pos_n += np2;
        }
      }
    }
    // seamE-1106 time2_k1_redist = th (> 0, read only when named): where the explicit
    // increment takes the stage-1 old vector below th E, E + a21 dt K1_E < th E, the
    // missing V K1_E, V (th E - E - a21 dt K1_E)/(a21 dt), is taken from the cells of
    // the same x1 (radial) column at distance d = 1, 2, .. time2_k1_redist_r, the two
    // cells at each d in proportion to their room V (E + a21 dt K1_E - th E)/(a21 dt)
    // (> 0 only; donors keep E + a21 dt K1_E >= th E).  The moved quantity is V K1_E, so
    // sum V K1_E is exact, E channel only; what rr cells cannot pay is given back
    // (unpaid, left to the posfix).  One thread serves the receivers of a column in
    // ascending i, so when two compete for one donor the inner one is served first.  K1
    // is the stored slope, so stage 2 reads the same field.  Column-local: Time2Init
    // FATALs unless x1 is one block (mesh nx1 = block nx1), so no neighbour block or
    // rank is read.  A no-op under time2_tableau = sdirk2 (a21 = 0).
    if (t2_k1rd > 0.0) {
      auto &id = pmy_pack->pmesh->mb_indcs;
      const int is = id.is, ie = id.ie, js = id.js, je = id.je, ks = id.ks, ke = id.ke;
      const int nmbr = pmy_pack->nmb_thispack - 1;
      auto k1d_ = t2k1;
      auto ue_ = u0;
      // FIX: pcoord->volume exists only on cs / sp (1x1x1x1 elsewhere: the
      // old vol(m,k,j,i) read out of bounds on Cartesian decks, so needs and donor rooms
      // were garbage: beams drained or blown up, results changed run to run).  Cartesian:
      // the uniform cell volume of the MeshBlock (a column is inside one block).
      auto &vola = pmy_pack->pcoord->volume;
      const bool usev = (vola.extent_int(0) >= nmbr + 1) && (vola.extent_int(3) > 1);
      auto &mbsz = pmy_pack->pmb->mb_size;
      auto vol = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i)
                 -> Real {
        return usev ? vola(m,k,j,i) : (mbsz.d_view(m).dx1*mbsz.d_view(m).dx2*
                                       mbsz.d_view(m).dx3);
      };
      const Real th = t2_k1rd, ca = c1;
      const int rr = t2_k1rd_r;
      Real nb = 0.0, nun = 0.0, mvd = 0.0, mvx = 0.0;
      Kokkos::parallel_reduce("m1_t2_k1rd", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
          DevExeSpace(), {0, ks, js}, {nmbr + 1, ke + 1, je + 1}),
      KOKKOS_LAMBDA(const int m, const int k, const int j, Real &sb, Real &su, Real &sm,
                    Real &sx) {
        for (int i = is; i <= ie; ++i) {
          const Real e0 = ue_(m,M1_E,k,j,i);
          const Real ev = e0 + ca*k1d_(m,M1_T2_E,k,j,i);
          if (ev >= th*e0) continue;
          sb += 1.0;
          // the energy (rate x volume) this cell lacks
          Real need = (th*e0 - ev)/ca*vol(m,k,j,i);
          k1d_(m,M1_T2_E,k,j,i) += need/vol(m,k,j,i);
          for (int d = 1; d <= rr && need > 0.0; ++d) {
            // both sides at distance d, in proportion to their room (no inward bias)
            const int ql = i - d, qh = i + d;
            Real rl = 0.0, rh = 0.0;
            if (ql >= is) {
              const Real eq = ue_(m,M1_E,k,j,ql);
              rl = fmax((eq + ca*k1d_(m,M1_T2_E,k,j,ql) - th*eq)/ca*vol(m,k,j,ql), 0.0);
            }
            if (qh <= ie) {
              const Real eq = ue_(m,M1_E,k,j,qh);
              rh = fmax((eq + ca*k1d_(m,M1_T2_E,k,j,qh) - th*eq)/ca*vol(m,k,j,qh), 0.0);
            }
            const Real rs = rl + rh;
            if (!(rs > 0.0)) continue;
            const Real f = (need < rs) ? need/rs : 1.0;
            if (ql >= is) {k1d_(m,M1_T2_E,k,j,ql) -= f*rl/vol(m,k,j,ql);}
            if (qh <= ie) {k1d_(m,M1_T2_E,k,j,qh) -= f*rh/vol(m,k,j,qh);}
            need -= f*rs;
            sm += f*rs;
            sx += f*rs*static_cast<Real>(d);
          }
          if (need > 0.0) {
            // not enough within rr cells: undo the unpaid part, posfix decides
            k1d_(m,M1_T2_E,k,j,i) -= need/vol(m,k,j,i);
            su += 1.0;
          }
        }
      }, Kokkos::Sum<Real>(nb), Kokkos::Sum<Real>(nun), Kokkos::Sum<Real>(mvd),
         Kokkos::Sum<Real>(mvx));
      {
        Real vv[4] = {nb, nun, mvd, mvx};
#if MPI_PARALLEL_ENABLED
        MPI_Allreduce(MPI_IN_PLACE, vv, 4, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
        t2_k1rd_n += vv[0];
        t2_k1rd_un += vv[1];
        t2_k1rd_mv += vv[2];
        t2_k1rd_mx += vv[3];
        t2_k1rd_calls += 1;
        if (global_variable::my_rank == 0 && t2_k1rd_calls % 50 == 0) {
          std::cout << "<rad_m1> K1RD calls " << t2_k1rd_calls << " cells/call "
                    << t2_k1rd_n/t2_k1rd_calls << " unpaid/call "
                    << t2_k1rd_un/t2_k1rd_calls << " moved rate "
                    << t2_k1rd_mv/t2_k1rd_calls << " mean cells moved "
                    << t2_k1rd_mx/fmax(t2_k1rd_mv, 1.0e-300) << std::endl;
        }
      }
    }
    t2_ys_ok = false;
  } else {
    int nmb1 = static_cast<int>(u0.extent(0)) - 1;
    int n3 = static_cast<int>(u0.extent(2)), n2 = static_cast<int>(u0.extent(3));
    int n1 = static_cast<int>(u0.extent(4));
    par_for("m1_t2_avg", DevExeSpace(), 0, nmb1, 0, M1_NVAR-1, 0, n3-1, 0, n2-1,
            0, n1-1,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
      u0_(m,n,k,j,i) = heun ? 0.5*(u1_(m,n,k,j,i) + u0_(m,n,k,j,i))
                            : (w1*u1_(m,n,k,j,i) + w0*u0_(m,n,k,j,i));
    });
    auto avgf = [&](DvceArray4D<Real> &a, DvceArray4D<Real> &b) {
      auto a_ = a;
      auto b_ = b;
      par_for("m1_t2_avgf", DevExeSpace(), 0, static_cast<int>(a.extent(0))-1,
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
  par_for("m1_t2_inc", DevExeSpace(), 0, nmb1, 0, M1_T2_NK-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
    Real v = c1*k1_(m,n,k,j,i);
    if (two) {v += c2*k2_(m,n,k,j,i);}
    in_(m,n,k,j,i) = v;
  });
  cln_lzc = c1*cln_lzs1 + (two ? c2*cln_lzs2 : 0.0);   // (clean_cs_lz_abs = 2)
  CleanDInc(two, c1, c2, heun ? 0.5 : w1, heun ? 0.5 : w0);
  cln_kc1 = c1;   // (clean_stag_kpos: the face part of the increment, CleanStagKPos)
  cln_kc2 = two ? c2 : 0.0;
  const bool ike = t2_ike && (t2kei.extent_int(0) > nmb1);
  if (ike) {
    auto e1_ = t2ke1, e2_ = t2ke2, ei_ = t2kei;
    par_for("m1_t2_ike_inc", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real v = c1*e1_(m,k,j,i);
      if (two) {v += c2*e2_(m,k,j,i);}
      ei_(m,k,j,i) = v;
    });
  }
  // time2_exch_split = psi: the carried exchange terms in bounded form
  if (t2_xsp && t2kx1.extent_int(0) > nmb1 && t2inc.extent_int(1) > M1_T2_EN) {
    auto x1_ = t2kx1, x2_ = t2kx2, l1_ = t2lx1, l2_ = t2lx2;
    auto ix_ = t2inc;
    const Real cc1 = c1, cc2 = c2, chcl = chat/c_light;
    auto vlx = pmy_pack->pcoord->volume;
    const bool hvx = (vlx.extent_int(0) > nmb1);
    Real sde = 0.0, snc = 0.0;
    Kokkos::parallel_reduce("m1_t2_xsp", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
        DevExeSpace(), {0, 0, 0, 0}, {nmb1 + 1, n3, n2, n1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &sd,
                  Real &sn) {
      const Real a1 = cc1;
      Real z = l1_(m,k,j,i)*a1;
      Real d = a1*(1.0/sqrt(1.0 + z*z) - 1.0)*x1_(m,k,j,i);
      if (two) {
        z = l2_(m,k,j,i)*cc2;
        d += cc2*(1.0/sqrt(1.0 + z*z) - 1.0)*x2_(m,k,j,i);
      }
      if (d != 0.0) {
        ix_(m,M1_T2_EN,k,j,i) += d;
        ix_(m,M1_T2_E,k,j,i) -= chcl*d;
        sd += fabs(d)*(hvx ? vlx(m,k,j,i) : 1.0);
        sn += 1.0;
      }
    }, Kokkos::Sum<Real>(sde), Kokkos::Sum<Real>(snc));
    t2_xsp_de += sde;
    t2_xsp_cells += snc;
    t2_xsp_calls += 1;
  }
  // <hydro>/am_vector: the spin reservoir's increment from its slopes
  if (am_vec && aminc.extent_int(0) > 0) {
    auto a1_ = amk1;
    auto a2_ = amk2;
    auto ai_ = aminc;
    par_for("m1_t2_aminc", DevExeSpace(), 0, static_cast<int>(aminc.extent(0)) - 1, 0,
            static_cast<int>(aminc.extent(1)) - 1, 0,
            static_cast<int>(aminc.extent(2)) - 1,
            0, static_cast<int>(aminc.extent(3)) - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real v = c1*a1_(m,k,j,i);
      if (two) {v += c2*a2_(m,k,j,i);}
      ai_(m,k,j,i) = v;
    });
  }
  // (implicit_gas_src = conservative): the same increment of the stored
  // cell F, from its slopes (the old vector of the write-back's F^old)
  csinc_ok = false;
  if (impl_gas_src_cons && csk1_ok && csk1.extent(0) == t2k1.extent(0)) {
    auto ck1_ = csk1;
    auto ck2_ = csk2;
    auto ci_ = csinc;
    par_for("m1_t2_csinc", DevExeSpace(), 0, nmb1, 0, 2, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
      Real v = c1*ck1_(m,n,k,j,i);
      if (two) {v += c2*ck2_(m,n,k,j,i);}
      ci_(m,n,k,j,i) = v;
    });
    csinc_ok = true;
  }
  // time2_floor_eint (rad_m1.hpp): floor-class cells keep the stage start's
  // eint in the old vector; the explicit exchange part moves to E
  if (t2_flr_rho > 0.0) {
    FluidRef flf = FluidRef::Get(pmy_pack);
    if (flf.on) {
      auto &id = pmy_pack->pmesh->mb_indcs;
      const int is = id.is, ie = id.ie, js = id.js, je = id.je, ks = id.ks, ke = id.ke;
      auto ug = flf.u0;
      const Real cl = c_light, ch = chat, rc = t2_flr_rho;
      const bool ikf = ike;
      auto eif_ = t2kei;
      auto ai_ = aminc;
      const bool amz = am_vec && aminc.extent_int(0) == t2inc.extent_int(0) &&
                       aminc.extent_int(1) == t2inc.extent_int(2) &&
                       aminc.extent_int(2) == t2inc.extent_int(3) &&
                       aminc.extent_int(3) == t2inc.extent_int(4);
      Real nf = 0.0;
      Kokkos::parallel_reduce("m1_t2_flreint", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
          DevExeSpace(), {0, ks, js, is}, {nmb1 + 1, ke + 1, je + 1, ie + 1}),
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &sum) {
        const Real dd = ug(m,IDN,k,j,i);
        if (!(dd > 0.0) || !(dd < rc)) return;
        // floor gas: the old vector is the stage start's gas (momentum, KE, eint); the
        // gas-energy increment goes to the radiation E (total energy kept)
        const Real dex = in_(m,M1_T2_EN,k,j,i);
        in_(m,M1_T2_EN,k,j,i) = 0.0;
        in_(m,M1_T2_M1,k,j,i) = 0.0;
        in_(m,M1_T2_M1+1,k,j,i) = 0.0;
        in_(m,M1_T2_M1+2,k,j,i) = 0.0;
        in_(m,M1_T2_E,k,j,i) += (ch/cl)*dex;
        if (amz) {ai_(m,k,j,i) = 0.0;}
        if (ikf) {eif_(m,k,j,i) = 0.0;}
        sum += 1.0;
      }, Kokkos::Sum<Real>(nf));
      t2_flr_n += nf;
    }
  }
  const bool xlm = (t2_k1xlim > 0.0) || t2_k1xpat;
  if (xlm) {Time2K1XLim(stage, c1, c2, true);}   // (gas side)
  if (t2_ffct) {Time2FaceFct(stage, c1, c2);}   // step 1
  if (xlm) {Time2K1XLim(stage, c1, c2, false);}   // (E side)
  CleanBndAcc(stage, c1, c2, (stage == 1) ? 0.0 : (heun ? 0.5 : w0), dt);
  t2_solve = (stage == 1) ? M1_T2S_STAGE1 : M1_T2S_STAGE2;
  t2_fail = false;
  nsub = 1;
  dt_sub = g*dt;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::Time2Stash
//! \brief pgrst-1021 time2_k1_reset: the radiation E at the hesdirk2 stage-loop start
//! (the state the FSAL slope K1 belongs to; the Driver calls it before the Heun stage 1)

void RadiationM1::Time2Stash() {
  t2_ys_ok = false;
  CleanDMem(true);
  CleanRayHistory(true);
  if (!t2_k1reset || t2ys.extent_int(0) < pmy_pack->nmb_thispack) return;
  auto ys_ = t2ys;
  auto u_ = u0;
  par_for("m1_t2_stash", DevExeSpace(), 0, static_cast<int>(t2ys.extent(0)) - 1, 0,
          static_cast<int>(t2ys.extent(1)) - 1, 0, static_cast<int>(t2ys.extent(2)) - 1,
          0, static_cast<int>(t2ys.extent(3)) - 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    ys_(m,k,j,i) = u_(m,M1_E,k,j,i);
  });
  t2_ys_ok = true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanRayHistory
//! \brief save: keep the rays' history (CleanRayKeep) at the stage-loop start; else put
//! it back (a stage failed)

void RadiationM1::CleanRayHistory(const bool save) {
  if (!clean) return;
  DvceArray5D<Real> *a[7] = {&cln_fvi, &cln_fh, &cln_top, &cln_rr, &cln_rpd, &cln_rnew,
                             &vgd_xf_d};
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
  for (int q = 0; q < 7; ++q) save ? cp(k.v[q], *a[q]) : cp(*a[q], k.v[q]);
  save ? cp(k.fvd, cln_fvd) : cp(cln_fvd, k.fvd);
  auto sw = [&](auto &x, auto &y) {if (save) {y = x;} else {x = y;}};
  sw(cln_fvblk, k.fvblk); sw(cln_rb, k.rb); sw(cln_rpb, k.rpb); sw(cln_rbld, k.rbld);
  sw(cln_rpt, k.rpt); sw(cln_rs, k.rs); sw(cln_rguard_sum, k.rg); sw(cln_roff_sum, k.ro);
  sw(cln_rns, k.rns); sw(cln_rok, k.rok); sw(cln_rfix, k.rfix); sw(cln_rlast, k.rlast);
  sw(cln_rjump, k.rjump);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::Time2Restore
//! \brief a stage failed: put U^n back (radiation from u1 and the face copies, gas from
//! the hydro's u1), refresh the hydro ghosts and primitives.  The Driver then reruns the
//! Heun stages and takes the step with backward Euler.

void RadiationM1::Time2Restore(Driver *pdrive) {
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
    // m1-mhd: the field back to U^n too (b1 is the stage-1 copy of b0 incl. ghosts)
    if (fl.mhd) {
      mhd::MHD *pm = pmy_pack->pmhd;
      Kokkos::deep_copy(DevExeSpace(), pm->b0.x1f, pm->b1.x1f);
      Kokkos::deep_copy(DevExeSpace(), pm->b0.x2f, pm->b1.x2f);
      Kokkos::deep_copy(DevExeSpace(), pm->b0.x3f, pm->b1.x3f);
    }
    FluidRefresh(pdrive);
    if (fl.mhd) EmagBuild(false);
  }
  // <hydro>/am_vector: the spin reservoir back to U^n too (am_spin1,
  // taken by the hydro's stage 1 and untouched since); the redo's hydro stage 1 copies
  // am_spin into am_spin1 again
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
  // (rule 9: the backward-Euler redo that carries the step is counted, not silent)
  if (global_variable::my_rank == 0) {
    std::printf("CLEAN WARNING cycle %d t2s %d: stage solve failed, the step is redone "
                "with backward Euler (stage fallback %.0f)\n", pmy_pack->pmesh->ncycle,
                t2_solve, t2_nfall);
  }
  // m1-h2div (tests_m1/runs_5j_h2div): the predictor increments were stored by the
  // stage solve that just FAILED (the store precedes the admissibility test), and the
  // backward-Euler redo would start its Picard loop from U^n + (dt/(g dt)) x that
  // increment (+ the order-2 rate term).  Where the redo does not converge either, its
  // result stays near that start, and the failed increment of the next step extrapolates
  // it again: the diffuse-wall shadow grew x1.5 per step to E ~ 1e12.  The redo starts
  // cold from U^n instead, and the next stage solves rebuild the history.
  pred_ok = false;
  pred2_ok = false;
  t2_solve = M1_T2S_NONE;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::Time2K1XLim
//! \brief xlim-1120 <rad_m1>/time2_k1_xlim = th (read only when named; default 0 = off,
//! bitwise): the exchange pair of the stage increment (the gas eint change de and E's
//! -(ch/cl) de) is scaled by theta in [0, 1] where it would take the old vector's E or
//! gas eint to <= th x its stage start.  The FSAL exchange slope is the last solve's
//! relaxation rate; extrapolated over a21 dt (with K2) it overshoots where radiation
//! flowing in is absorbed (HG prod KSPLIT: in 98 % of the stage-1 cells whose increment
//! leaves E <= 0 without its inflow the exchange pair dominates), and once face FCT
//! limits the donors the inflow no longer pays for it.  gas = true, before
//! time2_face_fct, where the gas loses (de < 0): theta the largest with eint + theta de
//! >= th eint, but no smaller than what keeps E + its local part (the increment without
//! the faces' divergence, all FCT can fall back to) >= th E, so that FCT, which reads the
//! result, can still keep E positive. gas = false, after FCT (the face part final), where
//! E loses (de > 0): theta the largest with E + dE >= th E.  The scaled-away exchange is
//! taken by the implicit stage (L-stable).  Per cell: E + gas energy exact ((c/ch) dE +
//! dEgas = 0), momentum and faces unchanged, no new state (restart and dt_retry bitwise);
//! untouched cells bitwise.  time2_k1_xpat: theta = P/(P + L) instead, P what the loser
//! holds without the exchange pair (E: E + its increment's other parts after FCT; gas:
//! eint), L its loss: the old vector's P^2/(P + L) > 0 for every P > 0, theta = 1 -
//! O(L/P) where the exchange is not stiff (Patankar-Euler weight, MPRK; no threshold).
//! time2_k1_xpat = pos: theta = 1 where L <= a P, a = 0.99 (the unlimited loser keeps
//! >= 0.01 P: untouched, bitwise as xpat off), else the loser keeps (1 - a)^2 P^2/(L -
//! (2a - 1) P) > 0, C1 in L at L = a P (theta = x - (1 - a)^2 x^2/(1 - (2a - 1) x),
//! x = P/L); positive by construction where the unlimited pair would leave <= 0
//! time2_k1_xpat_work (with pos): on the E side L is E's whole local loss to the gas,
//! the exchange and the work w that pays the gas KE change, (ch/c)(de + w) (t9h: every
//! stage-start E < 0 left by xpat = pos is the work's); the part taken back comes from
//! the gas eint (its gain, then at most 0.99 of what it holds), the momentum unchanged

void RadiationM1::Time2K1XLim(int stage, Real c1, Real c2, bool gas) {
  FluidRef flx = FluidRef::Get(pmy_pack);
  if (!flx.on) return;
  const bool two = (stage != 1);
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int is = id.is, ie = id.ie, js = id.js, je = id.je, ks = id.ks, ke = id.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto in_ = t2inc;
  auto u0_ = u0;
  auto ug = flx.u0;
  const bool etg = flx.use_etotgrav;
  auto phc = etg ? flx.phicc0 : DvceArray4D<Real>();
  const bool emh = flx.mhd;
  auto emg = emag0;
  const bool csk = cs_geom;
  auto ccl = pmy_pack->pcoord->cos_cell;
  // (time2_inc_energy = rad: t2kei is the KE increment, unchanged here)
  // (the faces' divergence as Time2FaceFct reads it: none stored yet = none taken)
  const bool fdo = gas && t2pf1_ok && (!two || t2pf2_ok) &&
                   t2pf1.extent_int(0) >= pmy_pack->nmb_thispack;
  auto p1_ = t2pf1;
  auto p2_ = t2pf2;
  const bool tr2 = trans_on, tr3 = trans_on && trans_x3;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto vol = pmy_pack->pcoord->volume;
  auto ar1 = pmy_pack->pcoord->area.x1f;
  auto ar2 = pmy_pack->pcoord->area.x2f;
  auto ar3 = pmy_pack->pcoord->area.x3f;
  const bool curv = (cs_geom || sph_geom) && vol.extent_int(0) > nmb1 &&
                    vol.extent_int(3) > 1 && ar1.extent_int(0) > nmb1;
  const Real cc1 = c1, cc2 = two ? c2 : 0.0;
  const Real chcl = chat/c_light, th = t2_k1xpat ? 0.0 : t2_k1xlim;
  const bool pat = t2_k1xpat, pos = t2_k1xpos, xw = t2_k1xwork && !gas;
  const Real pa = 0.99, pb2 = SQR(1.0 - pa), pc = 2.0*pa - 1.0;   // (xpat = pos)
  Real nx = 0.0;
  Kokkos::parallel_reduce("m1_t2_k1xlim", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
      DevExeSpace(), {0, ks, js, is}, {nmb1 + 1, ke + 1, je + 1, ie + 1}),
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &sum) {
    const Real dd = ug(m,IDN,k,j,i), idd = 1.0/fmax(dd, 1.0e-300);
    const Real q1 = ug(m,IM1,k,j,i), q2 = ug(m,IM2,k,j,i), q3 = ug(m,IM3,k,j,i);
    const Real cc = csk ? ccl(m,k,j) : 0.0;
    const Real ek0 = CsKeMom(q1, q2, q3, cc)*idd;
    const Real w = CsKeMom(q1 + in_(m,M1_T2_M1,k,j,i), q2 + in_(m,M1_T2_M1+1,k,j,i),
                           q3 + in_(m,M1_T2_M1+2,k,j,i), cc)*idd - ek0;
    const Real de = in_(m,M1_T2_EN,k,j,i) - w;
    const Real e0 = u0_(m,M1_E,k,j,i);
    Real tq = 1.0;
    const Real ls = xw ? de + w : de;   // (E's loss to the gas: the pair, + the work)
    if (!gas) {   // E loses: E + dE >= th E
      if (!(ls > 0.0)) return;
      const Real need = th*e0 - e0 - in_(m,M1_T2_E,k,j,i);
      if (pat) {   // (P = E + the increment's parts but the exchange)
        const Real pe = -need + chcl*ls, x = pe/(chcl*ls);
        if (pos && !(chcl*ls > pa*pe)) return;
        tq = (pe > 0.0) ? (pos ? x - pb2*x*x/(1.0 - pc*x) : pe/(pe + chcl*ls)) : 0.0;
      } else if (!(e0 > 0.0)) {
        tq = 0.0;
      } else if (need > 0.0) {
        tq = fmax(1.0 - need/(chcl*de), 0.0);
      }
    } else {   // the gas loses: eint + theta de >= th eint
      if (!(de < 0.0)) return;
      Real egr = etg ? dd*phc(m,k,j,i) : 0.0;
      if (emh) egr += emg(m,k,j,i);
      const Real eg0 = ug(m,IEN,k,j,i) - ek0 - egr;
      if (!(eg0 > 0.0) || (!pat && eg0 + de > th*eg0) || (pos && !(-de > pa*eg0))) {
        return;
      }
      const Real x = eg0/(-de);
      tq = pat ? (pos ? x - pb2*x*x/(1.0 - pc*x) : eg0/(eg0 - de))
               : fmax((1.0 - th)*eg0/(-de), 0.0);
      Real fd = 0.0;
      if (fdo) {
        Real g = 0.0;
        for (int d = 0; d < 3; ++d) {
          if ((d == 1 && !tr2) || (d == 2 && !tr3)) {continue;}
          const int dk = (d == 2), dj = (d == 1), di = (d == 0);
          for (int u = 0; u < 2; ++u) {
            const int kk = k + u*dk, jj = j + u*dj, ii = i + u*di;
            Real a;
            if (curv) {
              a = (d == 0) ? ar1(m,kk,jj,ii) : ((d == 1) ? ar2(m,kk,jj,ii) :
                                                          ar3(m,kk,jj,ii));
            } else {
              const auto &s = mbsize.d_view(m);
              a = (d == 0) ? s.dx2*s.dx3 : ((d == 1) ? s.dx1*s.dx3 : s.dx1*s.dx2);
            }
            const Real gf = a*(cc1*p1_(m,d,kk,jj,ii) + cc2*p2_(m,d,kk,jj,ii));
            g += u ? -gf : gf;
          }
        }
        const auto &s = mbsize.d_view(m);
        fd = g/(curv ? vol(m,k,j,i) : s.dx1*s.dx2*s.dx3);
      }
      // (E keeps the gain its local part needs: E + local part >= th E)
      const Real el = e0 + in_(m,M1_T2_E,k,j,i) - fd;
      const Real tqe = (el > th*e0) ? 1.0 - (el - th*e0)/(chcl*(-de)) : 1.0;
      // (xpat: where the gas's weight leaves E's local part at or below 0, the middle of
      // the open range (tqe, eg0/L) in which E and eint both stay > 0)
      const Real thi = fmin(eg0/(-de), 1.0);
      tq = (pat && !(tq > tqe) && tqe < thi) ? 0.5*(tqe + thi) : fmin(fmax(tq, tqe), 1.0);
    }
    if (!(tq < 1.0)) return;
    Real mv = (1.0 - tq)*ls;   // (the exchange taken back: the stage solve's)
    if (xw) {   // (the gas eint pays at most all its gain and 0.99 of what it holds)
      Real egr = etg ? dd*phc(m,k,j,i) : 0.0;
      if (emh) egr += emg(m,k,j,i);
      const Real eg0 = ug(m,IEN,k,j,i) - ek0 - egr;
      mv = fmin(mv, fmax(de, 0.0) + pa*fmax(eg0 + fmin(de, 0.0), 0.0));
    }
    in_(m,M1_T2_EN,k,j,i) -= mv;
    in_(m,M1_T2_E,k,j,i) += chcl*mv;
    sum += 1.0;
  }, Kokkos::Sum<Real>(nx));
  t2_k1xlim_n += nx;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::Time2FaceFct
//! \brief facebound-1112 step 1 <rad_m1>/time2_face_fct (rad_m1.hpp): Zalesak positivity
//! limiter on the face fluxes of the explicit stage increment t2inc (E and face F
//! channels).  G_f = A_f (c1 t2pf1 + c2 t2pf2) is the energy the increment carries across
//! face f (+ towards +x_d); V_i dE_i = sum_lower G - sum_upper G + V_i s_i with s_i the
//! cell-local remainder.  Only cells whose old vector E + dE is not > 0 are limited, by
//! R_i = (1 - 1e-12) Q_i/P_i, Q_i = V (E + dE) + P_i - I_i + I_i(R): what it holds w/o
//! face transport plus the inflow its donors' current R still let through; a face is
//! scaled by the R of its donor.  A limited donor cuts its receivers' inflow, so rounds
//! run until no cell is newly limited (R only decreases; a receiver that loses inflow
//! is limited again by exactly the deficit, which travels down the beam one cell per
//! round; cap 4096).  Bitwise wherever no face is limited.

//! \brief the limiter R's exchange (iw slot c): with the clean scheme's direct lists, the
//! first ghost layer only (the 7-point stencil reads the face ghosts; raw copies as
//! pbval_kr's, cs_noresample), one message per neighbour rank, without the clean
//! exchanges deferred to their next fill

void RadiationM1::Time2FctHalo(const int c) {
  if (cln_hdl.state == 1 && iw.extent_int(4) == cln_hdl.n1 &&
      iw.extent_int(3) == cln_hdl.n2 && iw.extent_int(2) == cln_hdl.n3) {
    std::vector<CleanGSeg<Real>> hp;
    hp.swap(cln_hpend);
    CleanGhostFill(iw, c, 1, false, 0, 1);
    cln_hpend.swap(hp);
  } else {
    ImplicitKrylovHalo(c);
  }
}

void RadiationM1::Time2FaceFct(int stage, Real c1, Real c2) {
  const bool two = (stage != 1);
  auto rhalo = [&](const int c) {Time2FctHalo(c);};
  // no stored face flux for K1 (first stage step, restart): K1 is the slope of the last
  // stage-2 solve, whose solution is U^n, so its face flux is the stage-start face flux
  // t2f* = f0x* (the restart file's F; the cmem applied flux is not used here)
  if (!two && !t2pf1_ok && t2pf1.extent_int(0) >= pmy_pack->nmb_thispack &&
      t2f1.extent_int(0) >= pmy_pack->nmb_thispack) {
    auto pf_ = t2pf1;
    auto g1_ = t2f1, g2_ = t2f2, g3_ = t2f3;
    const bool q2 = trans_on && t2f2.extent_int(0) > 0;
    const bool q3 = trans_on && trans_x3 && t2f3.extent_int(0) > 0;
    const Real crp = chat/c_light;
    par_for("m1_t2_ffct_fill", DevExeSpace(), 0, pf_.extent_int(0) - 1,
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
    t2_ffct_nop += 1.0;   // no stored face flux for this slope (start, restart)
    return;
  }
  if (iw.extent_int(1) <= Time2FctSlot(0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<rad_m1>/time2_face_fct needs the Krylov work slots "
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
  auto un_ = u1;   // (E^n: time2_face_fct_emin keeps th V E^n in a limited cell)
  const Real eth = t2_ffct_emin;
  // V (E + dE) unlimited (w), the full outflow po and inflow pi, the inflow pl that the
  // current R of the donors lets through, and the change dw of
  // V E_new by the current R (slot ca, ghosts valid); lim: a face next to it is limited
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
    w = cvol(m,k,j,i)*(u0_(m,M1_E,k,j,i) + in_(m,M1_T2_E,k,j,i));
  };
  // R = 1 everywhere, ghosts included (physical ghosts are never filled: unlimited);
  // the rounds read slot sa and write sb, then swap them (no copy; both slots keep the
  // same ghosts: after an exchange of sa they are copied to sb)
  for (const int c : {ca, cb}) {
    Kokkos::deep_copy(DevExeSpace(), Kokkos::subview(iw, Kokkos::ALL, c, Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL), 1.0);
  }
  int sa = ca, sb = cb;
  Real ncell = 0.0;
  int nround = 0;
  // (rad_m1.hpp zalesak_halo_skip): exchange R only after a round that
  // changed a cell within 2 ng of a MeshBlock face; the first exchange always runs
  const int zd = 2*id.ng;
  const bool zx2 = pmy_pack->pmesh->multi_d, zx3 = pmy_pack->pmesh->three_d;
  const bool zsk = zl_hskip;
  bool zfirst = true;
  std::chrono::steady_clock::time_point zt0;
  if (zl_prof) {Kokkos::fence(); zt0 = std::chrono::steady_clock::now();}
  const int zlr = zl_local;   // (rad_m1.hpp zalesak_local_rounds)
  // zalesak_active: a round reads slot sa and writes sb, an R that fell with its sign
  // bit set (-0 for 0); the next round evaluates only the cells with such a face
  // neighbour (ghosts included: the exchange copies the sign), the others keep |R|
  const bool zac = zl_active;   // (rad_m1.hpp zalesak_active)
  const bool swp = t2_ffct_swp;
  if (swp) {   // (time2_face_fct_sweep: x1 columns in place, slot cb the last round's R)
    // (one cell: the Zalesak test with x1 neighbours from ca, in place, and lateral ones
    // from cb, the last round's; true if its R fell)
    auto sweepq = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (hoi && hm_(m,k,j,i) > 0.0) return false;
      const Real rc = iw_(m,ca,k,j,i);
      Real po = 0.0, pi = 0.0, pl = 0.0;
      for (int d = 0; d < 3; ++d) {
        if ((d == 1 && !tr2) || (d == 2 && !tr3)) {continue;}
        const int dk = (d == 2) ? 1 : 0, dj = (d == 1) ? 1 : 0, di = (d == 0) ? 1 : 0;
        const int sn = (d == 0) ? ca : cb;
        Real g = gf(m, d, k, j, i);
        if (g > 0.0) {pi += g; pl += iw_(m,sn,k-dk,j-dj,i-di)*g;} else {po -= g;}
        g = gf(m, d, k+dk, j+dj, i+di);
        if (g > 0.0) {
          po += g;
        } else if (g < 0.0) {
          pi -= g;
          pl -= iw_(m,sn,k+dk,j+dj,i+di)*g;
        }
      }
      const Real v = cvol(m,k,j,i);
      const Real c = v*(u0_(m,M1_E,k,j,i) + in_(m,M1_T2_E,k,j,i)) - pi + po;
      const Real wf = eth*v*un_(m,M1_E,k,j,i);
      if (c + pl - rc*po > wf) return false;
      const Real q = c + pl - wf;
      const Real rq = (q > 0.0 && po > 0.0) ? fmin(1.0, (1.0 - 1.0e-12)*q/po) : 0.0;
      if (!(rq < rc)) return false;
      iw_(m,ca,k,j,i) = rq;
      return true;
    };
    // (a round evaluates only the cells whose inputs changed since their last evaluation,
    // kept per column as a range [lo, hi]: a lateral neighbour or a ghost next to them
    // changed in the last round, or, from the last round's downward pass, the cell below
    // them; the upward pass goes on past hi while the cell below changed, the downward
    // pass starts below the highest upward change and goes on while the cell above
    // changed.  Any other cell would find what it found before, so the R are those of
    // the full sweep)
    const int n3a = iw.extent_int(2), n2a = iw.extent_int(3);
    Kokkos::View<int****> rng("m1_t2_ffct_rng", nmb1 + 1, 4, n3a, n2a);
    auto rlo = Kokkos::subview(rng, Kokkos::ALL, 0, Kokkos::ALL, Kokkos::ALL);
    auto rhi = Kokkos::subview(rng, Kokkos::ALL, 1, Kokkos::ALL, Kokkos::ALL);
    auto nlo = Kokkos::subview(rng, Kokkos::ALL, 2, Kokkos::ALL, Kokkos::ALL);
    auto nhi = Kokkos::subview(rng, Kokkos::ALL, 3, Kokkos::ALL, Kokkos::ALL);
    Kokkos::deep_copy(DevExeSpace(), rlo, is);
    Kokkos::deep_copy(DevExeSpace(), rhi, ie);
    const int n1 = ie - is + 1, n2 = je - js + 1, n3 = ke - ks + 1;
    const int f1 = n3*n2, f2 = tr2 ? n3*n1 : 0, f3 = tr3 ? n2*n1 : 0;
    const int nf = 2*(f1 + f2 + f3);
    auto mark = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Kokkos::atomic_min(&nlo(m,k,j), i);
      Kokkos::atomic_max(&nhi(m,k,j), i);
    };
    for (int r = 0; r < 4096; ++r) {
      if (r > 0) {   // (the face ghosts the exchange changed mark the cell next to them)
        Kokkos::parallel_for("m1_t2_ffct_sg", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                             (nmb1 + 1)*nf),
        KOKKOS_LAMBDA(const int t) {
          const int m = t/nf;
          int u = t - m*nf, k, j, i, kc, jc, ic;
          if (u < 2*f1) {
            i = (u < f1) ? is - 1 : ie + 1;
            ic = (u < f1) ? is : ie;
            u %= f1;
            k = kc = ks + u/n2;
            j = jc = js + u % n2;
          } else if (u < 2*(f1 + f2)) {
            u -= 2*f1;
            j = (u < f2) ? js - 1 : je + 1;
            jc = (u < f2) ? js : je;
            u %= f2;
            k = kc = ks + u/n1;
            i = ic = is + u % n1;
          } else {
            u -= 2*(f1 + f2);
            k = (u < f3) ? ks - 1 : ke + 1;
            kc = (u < f3) ? ks : ke;
            u %= f3;
            j = jc = js + u/n1;
            i = ic = is + u % n1;
          }
          if (iw_(m,ca,k,j,i) != iw_(m,cb,k,j,i)) {mark(m, kc, jc, ic);}
        });
        Kokkos::deep_copy(DevExeSpace(), rlo, nlo);
        Kokkos::deep_copy(DevExeSpace(), rhi, nhi);
      }
      Kokkos::deep_copy(DevExeSpace(), nlo, ie + 1);
      Kokkos::deep_copy(DevExeSpace(), nhi, is - 1);
      Kokkos::deep_copy(DevExeSpace(), Kokkos::subview(iw, Kokkos::ALL, cb, Kokkos::ALL,
                        Kokkos::ALL, Kokkos::ALL), Kokkos::subview(iw, Kokkos::ALL, ca,
                        Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
      Real nv = 0.0;
      Kokkos::parallel_reduce("m1_t2_ffct_sw", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
          DevExeSpace(), {0, ks, js}, {nmb1 + 1, ke + 1, je + 1}),
      KOKKOS_LAMBDA(const int m, const int k, const int j, Real &a) {
        const int lo = rlo(m,k,j), hi = rhi(m,k,j);
        if (lo > hi) return;
        auto lat = [&](const int i) {   // (a change: the lateral neighbours' cell i)
          if (tr2) {
            if (j > js) {mark(m, k, j-1, i);}
            if (j < je) {mark(m, k, j+1, i);}
          }
          if (tr3) {
            if (k > ks) {mark(m, k-1, j, i);}
            if (k < ke) {mark(m, k+1, j, i);}
          }
        };
        int cl = ie + 1, ch = is - 1;   // (the upward pass's lowest / highest change)
        bool prev = false;
        for (int i = lo; i <= ie && (i <= hi || prev); ++i) {
          prev = sweepq(m, k, j, i);
          if (prev) {a += 1.0; lat(i); cl = (i < cl) ? i : cl; ch = i;}
        }
        prev = false;
        for (int i = ch - 1; i >= is && (i >= cl - 1 || prev); --i) {
          prev = sweepq(m, k, j, i);
          if (prev) {   // (the cell above saw this one only before: it runs next round)
            a += 1.0;
            lat(i);
            if (i < ie) {mark(m, k, j, i+1);}
          }
        }
      }, Kokkos::Sum<Real>(nv));
#if MPI_PARALLEL_ENABLED
      MPI_Allreduce(MPI_IN_PLACE, &nv, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
      if (!(nv > 0.0)) {break;}
      ncell += nv;
      nround = r + 1;
      rhalo(ca);
      zl_halos[0] += 1.0;
    }
  }
  // (zalesak_active: the face ghosts of the slot read next, |R|, where no exchange
  // runs: no R near a MeshBlock face fell, so they are those of the other slot)
  auto fghost = [&](const int c, const int cf) {
    const int n1 = ie - is + 1, n2 = je - js + 1, n3 = ke - ks + 1;
    const int f1 = n3*n2, f2 = tr2 ? n3*n1 : 0, f3 = tr3 ? n2*n1 : 0;
    const int nf = 2*(f1 + f2 + f3);
    Kokkos::parallel_for("m1_t2_ffct_fg", Kokkos::RangePolicy<>(DevExeSpace(), 0,
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
  for (int r = 0; r < (swp ? 0 : 4096); ++r) {   // until no cell is newly limited
    Real nv[2] = {0.0, 0.0};
    Real nl[2] = {0.0, 0.0};
    for (int s = 0; s < (zlr > 0 ? zlr : 1); ++s) {
    nv[0] = 0.0;
    nv[1] = 0.0;
    // (time2_face_fct_rounds: the slack trigger from round K, no inflow from round 2K)
    const int rk = t2_ffct_rk, ra = (rk < 0) ? -rk : rk;   // (K < 0: the slack only)
    const Real wtr = (ra > 0 && r >= ra) ? 0.5 : 1.0;
    const Real pin_ = (rk > 0 && r >= 2*rk) ? 0.0 : 1.0;
    // (zalesak_active: every cell in round 0 and in the rounds that change wtr / pin_)
    const bool zsel = zac && r > 0 && !(ra > 0 && r == ra) && !(rk > 0 && r == 2*rk);
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
        const Real wf = eth*cvol(m,k,j,i)*un_(m,M1_E,k,j,i);
        if (!(w + dw > wtr*wf)) {
          // what the cell holds without face transport + the inflow it still gets,
          // above the floor th V E^n
          const Real q = w + po - pi + pin_*pl - wf;
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
    par_reduce_flat2("m1_t2_ffct_r", 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &a,
                  Real &ab) {
      rnd(m, k, j, i, a, ab);
    }, Kokkos::Sum<Real>(nv[0]), Kokkos::Sum<Real>(nv[1]));
    if (zlr > 0) {   // local sub-round: keep the changes, no MPI, ghosts held
      if (!(nv[0] > 0.0)) {break;}
      nl[0] += nv[0];
      nl[1] += nv[1];
      zl_sub[0] += 1.0;
      par_for("m1_t2_ffct_cl", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        iw_(m,ca,k,j,i) = iw_(m,cb,k,j,i);
      });
    }
    }   // zalesak_local_rounds sub-round loop (one pass when off)
    if (zlr > 0) {nv[0] = nl[0]; nv[1] = nl[1];}
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, nv, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
    const Real nn = nv[0];
    if (!(nn > 0.0)) {break;}
    ncell += nn;
    nround = r + 1;
    if (zlr == 0) std::swap(sa, sb);   // (the round's R is read next)
    if (!zsk || zfirst || nv[1] > 0.0) {
      rhalo(sa);
      if (!zac && zlr == 0) {   // (the other slot keeps the same ghosts)
        Kokkos::deep_copy(DevExeSpace(), Kokkos::subview(iw, Kokkos::ALL, sb,
                          Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), Kokkos::subview(iw,
                          Kokkos::ALL, sa, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
      }
      zl_halos[0] += 1.0;
      zfirst = false;
    } else if (zac) {
      fghost(sa, sb);
    }
  }
  if (zac) {   // (the last round changed nothing: |R| of its read slot, ghosts too)
    auto r0 = Kokkos::subview(iw, Kokkos::ALL, sa, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
    auto r1 = Kokkos::subview(iw, Kokkos::ALL, ca, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
    par_for("m1_t2_ffct_abs", DevExeSpace(), 0, nmb1, 0, r0.extent_int(1) - 1, 0,
            r0.extent_int(2) - 1, 0, r0.extent_int(3) - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      r1(m,k,j,i) = fabs(r0(m,k,j,i));
    });
  } else if (sa != ca) {   // (the last round changed nothing: the R of both slots)
    Kokkos::deep_copy(DevExeSpace(), Kokkos::subview(iw, Kokkos::ALL, ca, Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL), Kokkos::subview(iw, Kokkos::ALL, sa,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  zl_rounds[0] += static_cast<Real>(nround);
  if (zl_prof) {
    Kokkos::fence();
    zl_wall[0] += std::chrono::duration<double>(std::chrono::steady_clock::now()
                                               - zt0).count();
  }
  t2_ffct_cells += ncell;
  t2_ffct_rmax = std::max(t2_ffct_rmax, static_cast<Real>(nround));
  if (nround == 0) {return;}   // nothing limited: t2inc untouched (bitwise)
  // apply: the E channel of every cell next to a limited face, then the face F channels
  Real lv[3] = {0.0, 0.0, 0.0}, qmin = 1.0e300;
  // (a flat range, as par_reduce_flat: the 4-D MDRange reduction took 8 ms on the A100)
  const int n1 = ie - is + 1, nji = (je - js + 1)*n1, nkji = (ke - ks + 1)*nji;
  Kokkos::parallel_reduce("m1_t2_ffct_e", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int idx, Real &sl, Real &sw, Real &se, Real &sq) {
    const int m = idx/nkji, k = ks + (idx%nkji)/nji, j = js + (idx%nji)/n1;
    const int i = is + idx%n1;
    Real w, po, pi, pl, dw;
    bool lim;
    cell(m, k, j, i, w, po, pi, pl, dw, lim, ca);
    const Real rc = iw_(m,ca,k,j,i);
    if (lim) {in_(m,M1_T2_E,k,j,i) += dw/cvol(m,k,j,i);}
    if (rc < 1.0) {sw += (1.0 - rc)*po;}
    if (!(hoi && hm_(m,k,j,i) > 0.0) && !(w + dw > 0.0)) {sl += 1.0;}
    if (lim) {   // (the limited old vector against E^n)
      const Real q = (u0_(m,M1_E,k,j,i) + in_(m,M1_T2_E,k,j,i))/un_(m,M1_E,k,j,i);
      if (q < eth) {se += 1.0;}
      sq = fmin(sq, q);
    }
  }, Kokkos::Sum<Real>(lv[0]), Kokkos::Sum<Real>(lv[1]), Kokkos::Sum<Real>(lv[2]),
     Kokkos::Min<Real>(qmin));
  // (clean staggered: the face F is the radiation momentum, which the limiter leaves
  // alone, so momentum and L_z stay exact; the solve keeps |F| <= cE)
  for (int d = 0; d < 3 && !clean; ++d) {
    if ((d == 1 && !tr2) || (d == 2 && !tr3)) {continue;}
    const int dk = (d == 2) ? 1 : 0, dj = (d == 1) ? 1 : 0, di = (d == 0) ? 1 : 0;
    par_for("m1_t2_ffct_f", DevExeSpace(), 0, nmb1, ks, ke + dk, js, je + dj, is, ie + di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real g = gf(m, d, k, j, i);
      Real rd = 1.0;
      if (g > 0.0) {rd = iw_(m,ca,k-dk,j-dj,i-di);}
      if (g < 0.0) {rd = iw_(m,ca,k,j,i);}
      if (rd < 1.0) {in_(m,M1_T2_F1+d,k,j,i) *= rd;}
    });
  }
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
//! \fn void RadiationM1::CleanBndAcc
//! \brief the history's x1 boundary energy (rad_m1.hpp cln_lacc), production's rule of
//! booking the flux that moved E: stage 0 starts a backward-Euler step (its solves add
//! h A Phi, CleanWriteBack); stage 1 / 2 the old
//! vector of a hesdirk2 stage, whose face part is w x (the energy of Y2, in the Heun
//! start) + A (c1 t2pf1 + c2 t2pf2), each face at its donor's limiter R
//! (time2_face_fct, as Time2FaceFct applies it); the stage solve adds g dt A Phi

void RadiationM1::CleanBndAcc(int stage, Real c1, Real c2, Real w, Real dt) {
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  if (!clean || !(sph_geom || cs_geom)) return;
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

//----------------------------------------------------------------------------------------
//! \fn Real RadiationM1::EbndFluxX1
//! \brief mhdfix-1025 Ebnd: L = sum A F1 over the physical inner x1 faces minus the outer
//! ones (the radiation energy flux into the domain through x1, all ranks).
//! side = 1 only the inner faces, side = 2 only the outer faces with
//! the sign flipped (outward positive); side = 0 is the original sum (bitwise).

Real RadiationM1::EbndFluxX1(const int side) {
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int is = id.is, ie = id.ie, js = id.js, je = id.je, ks = id.ks, ke = id.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto ff = f0x1;
  auto &a1 = pmy_pack->pcoord->area.x1f;
  auto &bcs = pmy_pack->pmb->mb_bcs;
  // ebnd_adv (problem/am_hist_ebnd): + the advective end-face flux
  const bool adv = ebnd_adv && impl_bc_advect;
  const int bclo = ibc_x1min, bchi = ibc_x1max;
  const bool mfl = impl_mface_lin, sqf = vcol_sqf;
  const bool vqs = (vet_col && vcol_sq) || hr_q;
  const Real elo = iebath_x1min, ehi = iebath_x1max;
  const bool ebc = iebath_hi_col_on;
  auto ebc_ = iebath_hi_col;
  auto iw_ = iw;
  auto cx1v = pmy_pack->pcoord->x1v;
  auto cx1f = pmy_pack->pcoord->xx1f;
  Real l = 0.0;
  Kokkos::parallel_reduce("m1_ebnd", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
      DevExeSpace(), {0, ks, js}, {nmb1 + 1, ke + 1, je + 1}),
  KOKKOS_LAMBDA(const int m, const int k, const int j, Real &sum) {
    auto phys = [](const BoundaryFlag f) {
      return !(f == BoundaryFlag::block || f == BoundaryFlag::periodic ||
               f == BoundaryFlag::shear_periodic || f == BoundaryFlag::panel ||
               f == BoundaryFlag::polar || f == BoundaryFlag::undef);
    };
    const bool plo = (side != 2) && phys(bcs.d_view(m, BoundaryFace::inner_x1));
    const bool phi = (side != 1) && phys(bcs.d_view(m, BoundaryFace::outer_x1));
    if (plo) {sum += a1(m,k,j,is)*ff(m,k,j,is);}
    if (phi) {sum -= a1(m,k,j,ie+1)*ff(m,k,j,ie+1);}
    // ebnd_adv: the advective enthalpy flux A E of implicit_bc_advect through
    // a Marshak end face, with the E the row takes (ImplicitSolve's x1 rows): outflow the
    // face E (implicit_marshak_face = linear) or the end cell's, inflow the bath
    if (adv && plo && bclo == M1_IBC_MARSHAK) {
      const Real ai = iw_(m,M1_IW_ADV,k,j,is), vi = iw_(m,M1_IW_V1,k,j,is);
      Real ee = elo;
      if (!(vi > 0.0)) {
        ee = iw_(m,M1_IW_EP,k,j,is);
        if (mfl && ie > is) {
          ee = M1SphMarshakFaceE(iw_(m,M1_IW_EP,k,j,is), iw_(m,M1_IW_EP,k,j,is+1),
                                 cx1v(m,is), cx1v(m,is+1), cx1f(m,is));
        }
      }
      sum += a1(m,k,j,is)*ai*ee;
    }
    if (adv && phi && bchi == M1_IBC_MARSHAK) {
      const Real ai = iw_(m,M1_IW_ADV,k,j,ie), vi = iw_(m,M1_IW_V1,k,j,ie);
      Real ee = ebc ? ebc_(m,k,j) : ehi;
      if (vi > 0.0) {
        ee = iw_(m,M1_IW_EP,k,j,ie);
        if (mfl && (!vqs || sqf) && ie > is) {
          ee = M1SphMarshakFaceE(iw_(m,M1_IW_EP,k,j,ie), iw_(m,M1_IW_EP,k,j,ie-1),
                                 cx1v(m,ie), cx1v(m,ie-1), cx1f(m,ie+1));
        }
      }
      sum -= a1(m,k,j,ie+1)*ai*ee;
    }
  }, Kokkos::Sum<Real>(l));
#if MPI_PARALLEL_ENABLED
  Real g;
  MPI_Allreduce(&l, &g, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  l = g;
#endif
  return (side == 2) ? -l : l;
}

//! \fn Real RadiationM1::EbndEsrc
//! \brief hg262-drive-1113 (ebnd_split): (chat/c) sum esrc V over the active cells, all
//! ranks (0 without esrc)

Real RadiationM1::EbndEsrc() {
  if (!esrc_on) return 0.0;
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int is = id.is, js = id.js, ks = id.ks;
  const int ni = id.nx1, nj = id.nx2, nk = id.nx3;
  const int nmb = pmy_pack->nmb_thispack;
  auto es_ = esrc;
  auto &vol = pmy_pack->pcoord->volume;
  Real l = 0.0;
  Kokkos::parallel_reduce("m1_ebnd_src", Kokkos::RangePolicy<>(DevExeSpace(), 0,
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

//! \fn void RadiationM1::EbndStart / EbndStage1 / EbndAccept
//! \brief mhdfix-1025 Ebnd (rad_m1.hpp): L at the step start (the face flux of the slope
//! K1), of the accepted stage-1 solve (K2), and the step's applied energy at acceptance
//! (a retried cycle accumulates nothing: the Driver accepts after the dt_retry check)

void RadiationM1::EbndStart() {
  if (!ebnd_on) return;
  if (transport != M1_TRANSPORT_IMPLICIT || MRActive() || csplit != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<rad_m1> Ebnd (energy-budget diagnostic) needs the "
              << "implicit transport without multi-rate / coupling_split" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (ebnd_split) {
    ebnd_l1i = EbndFluxX1(1);
    ebnd_l1o = EbndFluxX1(2);
    ebnd_s1 = ebnd_s3_ok ? ebnd_s3 : EbndEsrc();
    return;
  }
  ebnd_l1 = EbndFluxX1();
}

void RadiationM1::EbndStage1() {
  if (!ebnd_on) return;
  if (ebnd_split) {
    ebnd_l2i = EbndFluxX1(1);
    ebnd_l2o = EbndFluxX1(2);
    ebnd_s2 = EbndEsrc();
    return;
  }
  ebnd_l2 = EbndFluxX1();
}

void RadiationM1::EbndAccept(Real dt, bool t2step) {
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
  // (clean: f0x1 is the face momentum and the iw slots are not filled; the energy its
  // solves put through the physical x1 faces is CleanBndAcc's, inner and outer)
  Real lc[2] = {0.0, 0.0};
  if (clean && !(cln_ldt > 0.0) && global_variable::my_rank == 0) {
    std::printf("CLEAN WARNING cycle %d: Ebnd without the step's boundary energy "
                "(CleanBndAcc), its faces take 0\n", pmy_pack->pmesh->ncycle);
  }
  if (clean && cln_ldt > 0.0) {
    auto la = cln_lacc;
    auto &bcs = pmy_pack->pmb->mb_bcs;
    auto &id = pmy_pack->pmesh->mb_indcs;
    Kokkos::parallel_reduce("m1_ebnd_cln", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
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
    const Real l3i = clean ? 0.0 : EbndFluxX1(1), l3o = clean ? 0.0 : EbndFluxX1(2);
    ebnd_acc_in += clean ? lc[0] : applied(ebnd_l1i, ebnd_l2i, l3i);
    ebnd_acc_out += clean ? lc[1] : applied(ebnd_l1o, ebnd_l2o, l3o);
    const Real s3 = EbndEsrc();
    ebnd_acc_src += applied(ebnd_s1, ebnd_s2, s3);
    ebnd_s3 = s3;
    ebnd_s3_ok = true;
    ebnd_acc = ebnd_acc_in - ebnd_acc_out;
    return;
  }
  ebnd_acc += clean ? lc[0] - lc[1] : applied(ebnd_l1, ebnd_l2, EbndFluxX1());
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationM1::Time2StartProbe
//! \brief time2_start = slope: on a hesdirk2 step without a stored slope,
//! a backward-Euler probe solve of dt_p = time2_start_frac dt from U^n (radiation, faces
//! and gas; no hydro) stores K1 = (Y_p - U^n)/dt_p = f_I(U^n) + O(dt_p); U^n is then put
//! back and the step runs the two stages, so no step of the run is Lie split.  Returns
//! true when a probe was taken.  The predictor history the probe leaves is dropped (it
//! only moves Picard starting points).

bool RadiationM1::Time2StartProbe(Driver *pdrive, Real dt) {
  // VetGdRotStep asks for the probe at a rotation (any time2_start)
  if (time_scheme != M1_TIME_HESDIRK2 || (t2_start != 1 && !t2_rotprobe)) return false;
  // time2_auto_be: a skipped step takes the backward-Euler path, which stores
  // a fresh slope (of the new rotation set too); no probe
  if (t2_abe_now) {
    t2_rotprobe = false;
    return false;
  }
  t2_rotprobe = false;
  if (Time2Active()) return false;
  FluidRef fl = FluidRef::Get(pmy_pack);
  // U^n into the stage-start copies (Time2FormStage(1) overwrites them again)
  Kokkos::deep_copy(DevExeSpace(), u1, u0);
  Kokkos::deep_copy(DevExeSpace(), t2f1, f0x1);
  if (trans_on) {
    Kokkos::deep_copy(DevExeSpace(), t2f2, f0x2);
    if (trans_x3) {Kokkos::deep_copy(DevExeSpace(), t2f3, f0x3);}
  }
  if (fl.on) {Kokkos::deep_copy(DevExeSpace(), fl.u1, fl.u0);}
  // <mhd>: U^n = (u1, b1) for the solve (emag1, Time2VetStart); b1 still
  // holds the last step's start here (the MHD stage 1 copies b0 into it again)
  if (fl.on && fl.mhd) {
    Kokkos::deep_copy(DevExeSpace(), fl.pm->b1.x1f, fl.pm->b0.x1f);
    Kokkos::deep_copy(DevExeSpace(), fl.pm->b1.x2f, fl.pm->b0.x2f);
    Kokkos::deep_copy(DevExeSpace(), fl.pm->b1.x3f, fl.pm->b0.x3f);
  }
  // the cold-start face init (implicit_f0_init) of the probe moves the copy too
  // (ImplicitF0Init, t2_inprobe), so it survives the put-back
  Mesh *pmesh = pmy_pack->pmesh;
  t2_inprobe = true;
  t2_solve = M1_T2S_BESTORE;
  t2_fail = false;
  nsub = 1;
  dt_sub = t2_start_frac*dt;
  pdrive->ExecuteTaskList(pmesh, "m1_before_stagen", 1);
  pdrive->ExecuteTaskList(pmesh, "m1_stagen", 1);
  pdrive->ExecuteTaskList(pmesh, "m1_after_stagen", 1);
  t2_inprobe = false;
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
  t2_solve = M1_T2S_NONE;
  t2_fail = false;
  pred_ok = false;
  pred2_ok = false;
  t2_ok = true;
  t2_dtprev = dt;
  t2_nprobe += 1.0;
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationM1::Time2RotCorrect
//! \brief vet_gd_rotate_slope = correct (asked for by VetGdRotStep on the first step of a
//! new vet_gd rotation block, with a valid FSAL slope K1 = f_old(U^n) of the OLD set):
//! one stage-type solve of dt_c = vet_gd_rotate_slope_frac dt from U^n (radiation, faces,
//! gas; no hydro) with the NEW set and the old vector U^n - dt_c K1 (cell F: - dt_c
//! csk1).  Its slope (Y - old)/dt_c = K1 + z/dt_c with, to first order in z = Y - U^n,  z
//! = dt_c (I - dt_c A_new)^{-1} (f_new - f_old)(U^n): the stored slope is K1 + (I - dt_c
//! A)^{-1} Delta f -- the closure jump of the slope where dt_c A is small (the thin
//! regime, where keep leaves an O(dt Delta f) error per jump = first order), and keep's
//! K1 where the cell is stiff (there the L-stable stages forget a stage-1 slope error and
//! keep is second order; replacing all of K1 by a BE probe slope is what broke the stiff
//! ladders).  Only the jump steps change, so an O(dt) error of the correction itself
//! stays O(dt^2) globally.  The solve runs as a STAGE1 solve (its slope lands in t2k2 /
//! csk2 and is copied to t2k1 / csk1); the predictor history is kept (implicit_predictor
//! is off for this solve); U^n is put back.  With vet_gd_rotate_slope_reuse the stage-1
//! solve of the step keeps this solve's closure build (the same set, built by
//! Time2VetColAt(1) at U^n + theta dt K1_old with the gas at U^n instead of after the
//! Heun stage: an O(dt) state shift on the jump steps only). A failed correction solve
//! keeps K1 (keep).  Returns true when a solve was taken.

bool RadiationM1::Time2RotCorrect(Driver *pdrive, Real dt) {
  if (!t2_rotcorr) return false;
  t2_rotcorr = false;
  if (time_scheme != M1_TIME_HESDIRK2 || !Time2Active()) return false;
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
  // <mhd>: U^n = (u1, b1) for the solve (emag1, Time2VetStart); b1 still
  // holds the last step's start here (the MHD stage 1 copies b0 into it again)
  if (fl.on && fl.mhd) {
    Kokkos::deep_copy(DevExeSpace(), fl.pm->b1.x1f, fl.pm->b0.x1f);
    Kokkos::deep_copy(DevExeSpace(), fl.pm->b1.x2f, fl.pm->b0.x2f);
    Kokkos::deep_copy(DevExeSpace(), fl.pm->b1.x3f, fl.pm->b0.x3f);
  }
  // the old vector U^n - dt_c K1 (every channel; the cell F of conservative2 too)
  {
    auto k1_ = t2k1;
    auto in_ = t2inc;
    const Real c = -dtc;
    int nmb1 = static_cast<int>(t2inc.extent(0)) - 1;
    int n3 = static_cast<int>(t2inc.extent(2)), n2 = static_cast<int>(t2inc.extent(3));
    int n1 = static_cast<int>(t2inc.extent(4));
    par_for("m1_t2_rcinc", DevExeSpace(), 0, nmb1, 0, M1_T2_NK-1, 0, n3-1, 0, n2-1,
            0, n1-1,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
      in_(m,n,k,j,i) = c*k1_(m,n,k,j,i);
    });
    cln_lzc = c*cln_lzs1;   // (clean_cs_lz_abs = 2)
    if (t2_ike && t2kei.extent_int(0) > nmb1) {
      auto e1_ = t2ke1, ei_ = t2kei;
      par_for("m1_t2_ike_rcinc", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        ei_(m,k,j,i) = c*e1_(m,k,j,i);
      });
    }
    csinc_ok = false;
    if (impl_gas_src_cons && csk1_ok && csk1.extent(0) == t2k1.extent(0)) {
      auto ck1_ = csk1;
      auto ci_ = csinc;
      par_for("m1_t2_rccs", DevExeSpace(), 0, nmb1, 0, 2, 0, n3-1, 0, n2-1, 0, n1-1,
      KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
        ci_(m,n,k,j,i) = c*ck1_(m,n,k,j,i);
      });
      csinc_ok = true;
    }
  }
  const bool csok0 = csk1_ok;
  const bool pred_save = impl_pred;
  impl_pred = false;
  Mesh *pmesh = pmy_pack->pmesh;
  t2_inprobe = true;
  t2_solve = M1_T2S_STAGE1;
  t2_fail = false;
  nsub = 1;
  dt_sub = dtc;
  const Real its0 = impl_itsum;
  pdrive->ExecuteTaskList(pmesh, "m1_before_stagen", 1);
  pdrive->ExecuteTaskList(pmesh, "m1_stagen", 1);
  pdrive->ExecuteTaskList(pmesh, "m1_after_stagen", 1);
  vgd_rot_ncpass += static_cast<double>(impl_itsum - its0);
  t2_inprobe = false;
  impl_pred = pred_save;
  const bool ok = !t2_fail;
  if (ok) {
    Kokkos::deep_copy(DevExeSpace(), t2k1, t2k2);
    cln_lzs1 = cln_lzs2;   // (clean_cs_lz_abs = 2)
    if (t2_pfst && t2pf1.extent_int(0) > 0) {   // step 1
      Kokkos::deep_copy(DevExeSpace(), t2pf1, t2pf2);
      t2pf1_ok = t2pf2_ok;
    }
    if (t2_ike && t2ke1.extent_int(0) > 0) {
      Kokkos::deep_copy(DevExeSpace(), t2ke1, t2ke2);
    }
    if (t2_xsp && t2kx1.extent_int(0) > 0) {   // time2_exch_split
      Kokkos::deep_copy(DevExeSpace(), t2kx1, t2kx2);
      Kokkos::deep_copy(DevExeSpace(), t2lx1, t2lx2);
    }
    if (csok0 && impl_gas_src_cons && csk2.extent(0) == csk1.extent(0)) {
      Kokkos::deep_copy(DevExeSpace(), csk1, csk2);
    }
    vgd_rot_ncorr += 1.0;
  } else {
    vgd_rot_ncfail += 1.0;
  }
  csk1_ok = csok0;
  csinc_ok = false;
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
  t2_solve = M1_T2S_NONE;
  t2_fail = false;
  // the stage-1 solve keeps this build (vet_col predict | rebuild only)
  t2_vcreuse = ok && vgd_rot_creuse && vet_col && !vet_sc &&
               (t2_vcmode == 1 || t2_vcmode == 2);
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn M1T2VcpLaunch
//! \brief the m1_t2_vcp kernel of RadiationM1::Time2VetColAt, instantiated on the
//! ideal-gas/analytic-opacity tag.
template <typename Ctx, typename Idl>
void M1T2VcpLaunch(const Ctx &ctx_, Idl) {
  auto aa = std::get<0>(ctx_);
  auto bb = std::get<1>(ctx_);
  auto dt = std::get<2>(ctx_);
  auto efl = std::get<3>(ctx_);
  auto eos = std::get<4>(ctx_);
  auto etg = std::get<5>(ctx_);
  auto gam = std::get<6>(ctx_);
  auto gdt = std::get<7>(ctx_);
  auto iw_ = std::get<8>(ctx_);
  auto k1_ = std::get<9>(ctx_);
  auto kev = std::get<10>(ctx_);
  auto kf = std::get<11>(ctx_);
  auto kp = std::get<12>(ctx_);
  auto kss = std::get<13>(ctx_);
  auto opac_ = std::get<14>(ctx_);
  auto ot = std::get<15>(ctx_);
  auto otype = std::get<16>(ctx_);
  auto phicc = std::get<17>(ctx_);
  auto rref = std::get<18>(ctx_);
  auto tref = std::get<19>(ctx_);
  auto two = std::get<20>(ctx_);
  auto u0_ = std::get<21>(ctx_);
  auto uh = std::get<22>(ctx_);
  auto vn_ = std::get<23>(ctx_);
  auto nmb1 = std::get<24>(ctx_);
  auto ks = std::get<25>(ctx_);
  auto ke = std::get<26>(ctx_);
  auto js = std::get<27>(ctx_);
  auto je = std::get<28>(ctx_);
  auto is = std::get<29>(ctx_);
  auto ie = std::get<30>(ctx_);
  auto mhd = std::get<31>(ctx_);
  auto emag_ = std::get<32>(ctx_);
  auto osc = std::get<33>(ctx_);   // vet_scatter: rho kappa_e -> opac(M1_OP_S)
  auto lim = std::get<34>(ctx_);   // time2_pred_cap: theta per cell
  auto th_ = std::get<35>(ctx_);
  auto vcb = [=] KOKKOS_FUNCTION (Idl idl, const int m, const int k, const int j,
                                  const int i) M1_INL {
#if defined(KOKKOS_ENABLE_CUDA)
    // nvcc: an extended lambda may not capture a variable for the first time inside
    // an `if constexpr` branch, so name every capture up front (no code is
    // generated; other backends capture exactly what the instantiation reads).
    (void)aa; (void)bb; (void)dt; (void)efl; (void)eos; (void)etg; (void)gam; (void)gdt;
    (void)iw_; (void)k1_; (void)kev; (void)kf; (void)kp; (void)kss; (void)opac_; (void)ot;
    (void)otype; (void)phicc; (void)rref; (void)tref; (void)two; (void)u0_; (void)uh;
    (void)vn_; (void)mhd; (void)emag_; (void)osc; (void)lim; (void)th_;
#endif
    const Real d = uh(m,IDN,k,j,i);
    Real t, e;
    if (two) {
      t = vn_(m,3,k,j,i);
      e = vn_(m,4,k,j,i);
    } else {
      // the hydro u0 is already the stage-1 OLD vector here (ImplicitSolve added
      // t2inc = (1-g) dt K1 to it): g dt K1 completes dt K1.  time2_pred_cap: theta
      // dt K1 in all, i.e. g dt - (1 - theta) dt more (= g dt, bitwise, at theta = 1)
      const Real th = lim ? th_(m,k,j,i) : 1.0;
      const Real gdtt = lim ? (gdt - (1.0 - th)*dt) : gdt;
      const Real m1 = uh(m,IM1,k,j,i) + gdtt*k1_(m,M1_T2_M1,k,j,i);
      const Real m2 = uh(m,IM2,k,j,i) + gdtt*k1_(m,M1_T2_M1+1,k,j,i);
      const Real m3 = uh(m,IM3,k,j,i) + gdtt*k1_(m,M1_T2_M1+2,k,j,i);
      Real eint = uh(m,IEN,k,j,i) + gdtt*k1_(m,M1_T2_EN,k,j,i)
                  - 0.5*(m1*m1 + m2*m2 + m3*m3)/fmax(d, 1.0e-300);
      if (etg) eint -= d*phicc(m,k,j,i);
      if (mhd) eint -= emag_(m,k,j,i);   // m1-mhd: B is not changed by K1
      if constexpr (decltype(idl)::value) {
        t = ((gam-1.0)*fmax(eint, 0.0)/d);
      } else {
        t = eos.Temperature(d, fmax(eint, 0.0));
      }
      e = u0_(m,M1_E,k,j,i) + (th*dt)*k1_(m,M1_T2_E,k,j,i);
    }
    Real op, oe, of, os;
    if constexpr (decltype(idl)::value) {
      M1Opacities(otype, d, t, kp, kev, kf, kss, rref, tref, aa, bb, op, oe, of, os);
    } else if (otype == M1_OPAC_TABLE) {
      M1TableOpacities(ot, d, t, op, oe, of, os);
    } else {
      M1Opacities(otype, d, t, kp, kev, kf, kss, rref, tref, aa, bb, op, oe, of, os);
    }
    opac_(m,M1_OP_P,k,j,i) = d*op;
    opac_(m,M1_OP_E,k,j,i) = d*oe;
    if (osc && otype == M1_OPAC_TABLE) {M1ScatterEos(ot, eos, d, t, of, os);}
    opac_(m,M1_OP_T,k,j,i) = d*(of + os);
    if (osc) {opac_(m,M1_OP_S,k,j,i) = d*os;}
    vn_(m,0,k,j,i) = iw_(m,M1_IW_EN,k,j,i);
    vn_(m,1,k,j,i) = iw_(m,M1_IW_TP,k,j,i);
    vn_(m,2,k,j,i) = iw_(m,M1_IW_KT,k,j,i);
    iw_(m,M1_IW_EN,k,j,i) = fmax(e, efl);
    iw_(m,M1_IW_TP,k,j,i) = t;
    iw_(m,M1_IW_KT,k,j,i) = d*(of + os);
  };
  par_for_lb("m1_t2_vcp", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) M1_INL {
    vcb(Idl{}, m, k, j, i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::Time2VetSourceInit
//! \brief time2_vet_source = relaxed | state | blend (was
//! time2_vet_col_source, its alias on vet_col): the VET source of the predict / rebuild
//! build, for every formal solution (#3).  vet_source = relaxed relaxes
//! the build state's gas and E by a local backward-Euler exchange over dt, i.e. it moves
//! the source one step AHEAD of the state it is given.  At U^n (lag, backward Euler)
//! that is the step's end; at the predicted P = U^n + dt K1, already at t^{n+1} +
//! O(dt^2), it lands near t^{n+2}: an O(dt) source error, FIRST order in time wherever
//! the exchange is resolved and out of equilibrium (vcpred_1010 neq1_u_T: E 1.71 1.57
//! 1.40, T 1.67 1.47 1.23).  state = the source of P itself (eps B + (1 - eps) E, or
//! VetScatterMix with J = E_P under vet_scatter): 1.92 1.87 1.91 / 1.96 1.98 1.99.
//! blend = S*(0) + w (S*(dt) - S*(0)), w = x^2/(x0^2 + x^2) with x the
//! cell's exchange stiffness c dt rho kappa (radiation + gas temperature): the state's
//! own source where the exchange is resolved (x << x0, O(dt^3) from it), the relaxed
//! source where it is stiff (VetRelaxedSourceB); time2_vet_source_x0 (alias
//! time2_vet_col_source_x0; default 1).  lag / extrap keep vet_source.
//! DEFAULT on vet_col / vet_gd: blend on a fresh run with predict | rebuild (where it
//! acts), relaxed on a restart whose file lacks the key; vet_sc: state (its old source;
//! OPEN, see VetInit); the resolved value is echoed; explicit input overrides.

void RadiationM1::Time2VetSourceInit(ParameterInput *pin, const bool act,
                                     const bool col) {
  const std::string key = col ? M1KeyName(pin, "time2_vet_source",
                                          {"time2_vet_col_source"})
                              : std::string("time2_vet_source");
  if (!act && !pin->DoesParameterExist("rad_m1", key)) {return;}
  // #3 (OPEN, rad_m1_vet.cpp VetInit): vet_sc keeps state (fresh and
  // restart), vet_col / vet_gd blend (restart lacking the key: relaxed)
  const char *dflt = col ? (global_variable::restart_run ? "relaxed" : "blend") : "state";
  std::string vs = act ?
      pin->GetOrAddString("rad_m1", key, dflt) :
      pin->GetString("rad_m1", key);
  if (vs.compare("relaxed") == 0) {
    t2_vcsrc = false;
  } else if (vs.compare("state") == 0) {
    t2_vcsrc = true;
  } else if (vs.compare("blend") == 0) {
    t2_vcsbl = true;
    const std::string kx = col ? M1KeyName(pin, "time2_vet_source_x0",
                                           {"time2_vet_col_source_x0"})
                               : std::string("time2_vet_source_x0");
    vcol_sx0 = pin->GetOrAddReal("rad_m1", kx, 1.0);
    if (!(vcol_sx0 > 0.0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/" << kx << " must be > 0" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<rad_m1>/" << key << " = '" << vs
              << "' is not a choice (relaxed | state | blend)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::Time2VetColExtrap2
//! \brief rbfast-1021 time2_vet_col = extrap2 (rad_m1.hpp): the hand-over fields of the
//! build at U^n (just made by Time2VetStart) extrapolated to t^{n+1} with the previous
//! step's build; the raw build becomes the history

void RadiationM1::Time2VetColExtrap2() {
  const Real tn = pmy_pack->pmesh->time;
  const Real dtn = pmy_pack->pmesh->dt;
  const Real r = (x2_ok && tn > x2_t) ? dtn/(tn - x2_t) : 0.0;
  // the vet_gd fields only within one rotated direction set
  // (vet_gd_rotate_xfade: the cross-faded closure is continuous across the
  // set changes, so its fields are extrapolated across them too)
  const Real rg = (vgd_on && (vgd_alpha == x2_alpha || vgd_xf)) ? r : 0.0;
  const int xf = x2_fields;
  const Real rk = (xf & 1) ? r : 0.0, rl = (xf & 2) ? rg : 0.0, rh = (xf & 4) ? rg : 0.0;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n3 = static_cast<int>(tau_ten.extent(2));
  const int n2 = static_cast<int>(tau_ten.extent(3));
  const int n1 = static_cast<int>(tau_ten.extent(4));
  const int ntt = std::min(static_cast<int>(tau_ten.extent(1)), M1_X2_NTT);
  auto tt_ = tau_ten;
  auto pt_ = x2_tt;
  const Real fkm = vcol_axis_flux ? (1.0/3.0) : vcol_fkmin;
  const int lat0 = M1_TT_LAT0;
  par_for("m1_t2_x2tt", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int c = 0; c < ntt; ++c) {
      const Real gn = tt_(m,c,k,j,i);
      Real v = gn;
      if (c == 0) {
        v = fmin(fmax(gn + rk*(gn - pt_(m,0,k,j,i)), fkm), 1.0);
      } else if (c >= lat0) {
        v = gn + rl*(gn - pt_(m,c,k,j,i));
      }
      pt_(m,c,k,j,i) = gn;
      tt_(m,c,k,j,i) = v;
    }
  });
  const int nhr = std::min(static_cast<int>(vgd_hr.extent(1)), M1_X2_NHR);
  if (vgd_on && vgd_hr.extent_int(0) > nmb1 && nhr > 0 &&
      static_cast<int>(vgd_hr.extent(4)) == n1) {
    auto hr_ = vgd_hr;
    auto ph_ = x2_hr;
    const int nh3 = static_cast<int>(vgd_hr.extent(2));
    const int nh2 = static_cast<int>(vgd_hr.extent(3));
    par_for("m1_t2_x2hr", DevExeSpace(), 0, nmb1, 0, nh3-1, 0, nh2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      // r+_r < 0 flags a cell without a formal solution: it keeps its build
      const bool ok = (hr_(m,0,k,j,i) >= 0.0) && (ph_(m,0,k,j,i) >= 0.0);
      for (int c = 0; c < nhr; ++c) {
        const Real gn = hr_(m,c,k,j,i);
        Real v = gn;
        if (ok && c < M1_VGD_HRNX) {
          v = gn + rh*(gn - ph_(m,c,k,j,i));
          if (c < 6) {
            v = (c % 2 == 0) ? fmin(fmax(v, 0.0), 1.0) : fmin(fmax(v, -1.0), 0.0);
          } else if (c == M1_VGD_HREP) {
            v = fmin(fmax(v, 0.0), 1.0);
          } else {
            v = fmax(v, 0.0);
          }
        }
        ph_(m,c,k,j,i) = gn;
        hr_(m,c,k,j,i) = v;
      }
    });
  }
  if (vcol_sq && vcol_q.extent_int(0) > nmb1) {
    auto q_ = vcol_q;
    auto pq_ = x2_q;
    const Real qlo = vcol_qmin, qhi = vcol_qmax;
    const int nq3 = static_cast<int>(vcol_q.extent(1));
    const int nq2 = static_cast<int>(vcol_q.extent(2));
    // the gd top cell's q under vet_gd_replace
    const Real rq = (xf & 8) ? (vgd_replace ? rg : r) : 0.0;
    par_for("m1_t2_x2q", DevExeSpace(), 0, nmb1, 0, nq3-1, 0, nq2-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      const Real qn = q_(m,k,j);
      q_(m,k,j) = fmin(fmax(qn + rq*(qn - pq_(m,k,j)), qlo), qhi);
      pq_(m,k,j) = qn;
    });
  }
  x2_t = tn;
  x2_alpha = vgd_alpha;
  x2_ok = true;
}

int RadiationM1::Time2FctSlot(const int q) const {
  return clean ? (q ? M1_IW_S2 : M1_IW_S1) : (q ? M1_IW_KS : M1_IW_KZ);
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::CleanRstNq
//! \brief closure components per record in the restart files of the clean extrap history
//! and clean_rot_ring, in push 16's order (CleanRstSlot): all CLN_N (the last, CLN_JR,
//! the rays' J below clean_ray_min)

int RadiationM1::CleanRstNq() const {
  return CLN_N;
}

// (clean_closure = iterate: the closure in place, the last iterate's, is state; none in
// a file of another closure, Time2RstAccept)
int RadiationM1::CleanRstNcl() const {
  if (!clean) return 0;
  if (t2_rcl >= 0) return t2_rcl;   // (a file's, Time2RstAccept)
  return cln_closit ? CleanNs(pmy_pack->pmesh->use_cubed_sphere, cln_rfh) : 0;
}

int RadiationM1::CleanRstNfv() const {
  if (!clean) return 0;
  // (the cubed sphere's sc rays; before the first sweep those a restart file holds)
  // (the sc rays' cln_fvi: their slots of the union on a ring, CleanGdFaceSums; version
  // 5: dense, CleanFvRstCopy)
  // (version 6, staggered: their face half ranges cln_fh instead, none where every
  // solve rebuilds them first, CleanFhBuilt)
  if (cln_cssc && t2_rstv >= 6) return CleanFhBuilt() ? 0 : 12;
  if (cln_cssc && t2_rstv >= 5) {
    return (cln_fvi.extent_int(0) > 0) ? CleanFvRstN() : t2_nfv;
  }
  if (cln_cssc && cln_fvi.extent_int(0) > 0) return cln_fvi.extent_int(1);
  if (cln_cssc) return (cln_fvn > 0) ? cln_fvn : t2_nfv;
  return cln_scmem ? cln_top.extent_int(1) : 0;
}

//----------------------------------------------------------------------------------------
//! \fn RadiationM1::Time2Rst*
//! \brief the restart channels: K1 (M1_T2_NK), then ipred2 (3, implicit_predictor), then
//! vet_prev (M1_T2_NVET, closure = vet_sc)

int RadiationM1::Time2RstNchWant() {
  // the hesdirk2 slope channels, then (any time scheme) the clean scheme's state.
  // (version 4, the clean scheme: no channels for what clean_skip_time2_hist leaves
  // zero, and the ring's offsets dense, CleanRingRstCopy; version 5: the cubed sphere's
  // sc rays' intensities dense too, CleanFvRstCopy)
  const bool sk = clean && cln_not2h && t2_rstv >= 4;
  const int np = (impl_pord == 2) ? 5 : 3;
  const int nt = (time_scheme != M1_TIME_HESDIRK2) ? 0 :
         M1_T2_NK + ((impl_pred && !sk) ? np : 0) + (vet_sc ? M1_T2_NVET : 0)
         + ((vet_sc && vsc_every > 1) ? 2*M1_T2_NVET : 0)    // vet_sc_every: D0, D1
         + ((vet_col && t2_vcmode == 4 && !sk) ? M1_X2_NCH : 0);   // extrap2
  return nt + ((force_ref == M1_FREF_CLEAN) ? 4 + cln_lzt : 0)   // the reference, set
         + ((force_ref == M1_FREF_WB_MD) ? 10 : 0)           // wb_arad_md: 9, set
         + ((clean && cln_extrap > 0) ? 2*CleanRstNq() + 1 : 0)    // extrap history
         + CleanRstNcl()                                     // iterate: the closure
         + ((clean && cln_ring) ? (t2_rstv < 4 ? vra_k*CleanRstNq() : CleanRingRstN()) + 1
                                : 0)                          // clean_rot_ring
         + ((clean && cln_fsolve && cln_fkin) ? 6 : 0)       // stag fkin: d and its K1
         + CleanRstNfv();                                    // the rays' intensities
}

bool RadiationM1::Time2RstAccept(const int nfile, const int ver) {
  // (a version-1 file, before 8349d48ff, holds d without its slope K1_d: read with K1_d =
  // 0 and a backward-Euler first step)
  t2_rstd3 = (ver < 2) && clean && cln_fsolve && cln_fkin;
  t2_rstv = clean ? ver : 5;
  t2_rcl = -1;
  int nb = Time2RstNchWant() - CleanRstNfv() - (t2_rstd3 ? 3 : 0);
  t2_nfv = (clean && (cln_cssc || cln_scmem)) ? nfile - nb : 0;
  // (a file of another clean_closure: the closure in place of an iterate run, CleanRstNcl,
  // is skipped by another closure; without it, an iterate run's first solve takes the
  // predicted closure, CleanSolve)
  const int ns = CleanNs(pmy_pack->pmesh->use_cubed_sphere, cln_rfh);
  auto fits = [&](const int n) {
    return cln_cssc ? (nfile - n == 0 || nfile - n == 12) : (nfile == n);
  };
  if (clean && !cln_scmem && !fits(nb) && fits(cln_closit ? nb - ns : nb + ns)) {
    t2_rcl = cln_closit ? 0 : ns;
    nb += cln_closit ? -ns : ns;
    t2_nfv = cln_cssc ? nfile - nb : 0;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> NOTE: restart from a file of another clean_closure: "
                << (cln_closit ? "no closure in it, the first solve takes the predicted "
                    "one" : "its closure (clean_closure = iterate) is not read")
                << std::endl;
    }
  }
  // (clean_ring_fh: its half ranges are ring fields, no cln_fh channels; until the
  // layout records the key, only a file of a run with it fits exactly, and a file of
  // one is refused without it: 1 .. 11 channels over)
  const bool fho = cln_rfh ? t2_nfv != 0 : (t2_nfv > 0 && t2_nfv < 12);
  if (clean && cln_cssc && ver >= 6 && fho) {
    std::cout << "### FATAL ERROR <rad_m1> restart: a file " << (cln_rfh ? "not " : "")
              << "written with clean_ring_fh (" << t2_nfv << " channels over)"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // rotwb2-1121: on the cubed sphere (version 6) the face half ranges are the last 0 or
  // 12 channels (cln_fh).  A surplus beyond them is a file written with clean_rot_ring:
  // read with the ring off its ring channels are skipped (they were read into cln_fh
  // beyond its 12 channels), anything else is refused
  t2_rskip = 0;
  if (clean && cln_cssc && ver >= 6 && t2_nfv > 12) {
    t2_rskip = cln_ring ? -1 : CleanRingRstN() + 1;
    const int nh = t2_nfv - t2_rskip;
    if (t2_rskip < 0 || (nh != 0 && nh != 12)) {
      std::cout << "### FATAL ERROR <rad_m1> restart: the clean block has " << nfile
                << " channels, not this run's layout (" << nb << " + 0 or 12) nor it "
                << "with the ring's offsets" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    t2_nfv = nh;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> NOTE: restart with clean_rot_ring off: the file's ring "
                << "offsets (" << t2_rskip << " channels) are not read" << std::endl;
    }
  }
  return (nfile == nb + t2_nfv + t2_rskip) && (t2_nfv >= 0);
}
int RadiationM1::Time2RstNch() {
  t2_rstv = 7;   // (the writer: version 7, after a restart from an older file too)
  t2_rcl = -1;
  // (the clean scheme always: its state, and K1 with cln_k1ok, also in backward-Euler
  // steps; restart.cpp stores t2_ok)
  return (t2_ok || clean || time_scheme != M1_TIME_HESDIRK2) ? Time2RstNchWant() : 0;
}

// m0 > 0 packs MeshBlocks m0 .. m0+nmb-1 into a(0 .. nmb-1) (per-MB staging)
void RadiationM1::Time2RstPack(DvceArray5D<Real> &a, int nmb, int m0, int w0) {
  auto mb = std::make_pair(0, nmb);
  auto ms = std::make_pair(m0, m0 + nmb);
  auto AL = Kokkos::ALL;
  const bool h2 = (time_scheme == M1_TIME_HESDIRK2);
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
    put(c, M1_T2_NK, t2k1, 0);
    c += M1_T2_NK;
  }
  // (clean_skip_time2_hist: no predictor and extrap2 channels, version 4; the file is
  // written in this run's layout, t2_rstv = 5)
  const bool noh = clean && cln_not2h;
  if (h2 && impl_pred && !noh) {
    const int np = (impl_pord == 2) ? 5 : 3;
    put(c, np, ipred2, 0);
    c += np;
  }
  if (h2 && vet_sc) {
    put(c, M1_T2_NVET, vet_prev, 0);
    c += M1_T2_NVET;
    if (vsc_every > 1) {
      const int nv = M1_T2_NVET;
      put(c, nv, vsc_d0, 0);
      put(c + nv, nv, vsc_d1, 0);
      c += 2*nv;
    }
  }
  // force_reference = clean (wb_arad_md): the reference force (and potentials) and
  // whether it is set
  if (force_ref == M1_FREF_CLEAN || force_ref == M1_FREF_WB_MD) {
    const int nr = (force_ref == M1_FREF_WB_MD) ? 9 : 3 + cln_lzt;
    const bool have = (cln_fref.extent_int(0) >= m0 + nmb);
    if (have) put(c, nr, cln_fref, 0);
    one(c + nr, (have && cln_fref_set) ? 1.0 : 0.0);
    c += nr + 1;
  }
  // clean_closure = extrap*: the ray moments at t^n and at the end of the step, then
  // t^n, that end time and the flags at i = 0..4
  const int nq = CleanRstNq();
  if (clean && cln_extrap > 0) {
    const bool have = (cln_mn.extent_int(0) >= m0 + nmb);
    for (int n = 0; n < nq && have; ++n) {
      put(c + n, 1, cln_mn, CleanRstSlot(n));
      put(c + nq + n, 1, cln_mc, CleanRstSlot(n));
    }
    auto a_ = a;
    const int cq = c + 2*nq - w0;
    const Real v0 = cln_tn, v1 = cln_tm, v2 = (have && cln_have_n) ? 1.0 : 0.0;
    const Real v3 = (have && cln_have_c) ? 1.0 : 0.0, v4 = cln_cand_be ? 1.0 : 0.0;
    if (cq >= 0 && cq < nw) {
      par_for("m1_t2_clnrst", DevExeSpace(), 0, nmb - 1, 0, na3 - 1, 0, na2 - 1, 0,
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
  if (clean && cln_ring) {
    const int nd = CleanRingRstN();
    const bool have = cln_rok && (cln_rr.extent_int(0) >= m0 + nmb);
    CleanRingRstCopy(a, c - w0, nmb, m0, true, have);
    auto a_ = a;
    const int cq = c + nd - w0;
    const Real v0 = have ? 1.0 : 0.0, v1 = static_cast<Real>(cln_rb);
    const Real v2 = static_cast<Real>(cln_rpb), v3 = cln_rpt;
    if (cq >= 0 && cq < nw) {
      par_for("m1_t2_ringrst", DevExeSpace(), 0, nmb - 1, 0, na3 - 1, 0, na2 - 1, 0,
              na1 - 1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        a_(m,cq,k,j,i) = (i == 0) ? v0 : ((i == 1) ? v1 : ((i == 2) ? v2 :
                         ((i == 3) ? v3 : 0.0)));
      });
    }
    c += nd + 1;
  }
  // clean_stag_fkin: the energy flux's memory d of the three face families (CLF_N + 3),
  // then its slope K1 (CLF_KD); before the first solve (CleanInit) those of the file read
  if (clean && cln_fsolve && cln_fkin) {
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
  // extrap2: the raw builds of the last step, q at i = 0, time, alpha
  if (h2 && vet_col && t2_vcmode == 4 && !noh) {
    put(c, M1_X2_NTT, x2_tt, 0);
    put(c + M1_X2_NTT, M1_X2_NHR, x2_hr, 0);
    c += M1_X2_NTT + M1_X2_NHR;
    auto a_ = a;
    auto q_ = x2_q;
    const Real t = x2_ok ? x2_t : -1.0e300, al = x2_alpha, ds = x2_dts;
    const int cq = c - w0;
    const int mq = m0;
    par_for("m1_t2_x2rst", DevExeSpace(), 0, nmb - 1, 0, na3-1, 0, na2-1, 0, na1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (cq >= 0 && cq < nw) a_(m,cq,k,j,i) = (i == 0) ? q_(m+mq,k,j) : 0.0;
      if (cq+1 >= 0 && cq+1 < nw) a_(m,cq+1,k,j,i) = t;
      if (cq+2 >= 0 && cq+2 < nw) a_(m,cq+2,k,j,i) = al;
      if (cq+3 >= 0 && cq+3 < nw) a_(m,cq+3,k,j,i) = ds;
    });
    c += 4;
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
  } else if (cln_cssc && CleanRstNfv() > 0) {
    CleanFvRstCopy(a, c, nmb, m0, true);
  } else if (CleanRstNfv() > 0 && CleanRstFv().extent_int(0) >= m0 + nmb) {
    Kokkos::deep_copy(Kokkos::subview(a, mb, std::make_pair(c, c + CleanRstNfv()), AL, AL,
                                      AL),
                      Kokkos::subview(CleanRstFv(), ms, AL, AL, AL, AL));
  }
}

void RadiationM1::Time2RstSet(int ch, const HostArray4D<Real> &w, int nmb) {
  auto mb = std::make_pair(0, nmb);
  auto AL = Kokkos::ALL;
  const bool h2 = (time_scheme == M1_TIME_HESDIRK2);
  // clean_skip_time2_hist: a file's predictor / extrap2 channels are read past (a NOTE
  // once where they hold values: files of the old scheme)
  auto skip = [&w]() {
    static bool said = false;
    for (size_t q = 0; q < w.span() && !said && global_variable::my_rank == 0; ++q) {
      if (w.data()[q] == 0.0) continue;
      std::cout << "<rad_m1> NOTE: clean_skip_time2_hist: the restart file's implicit "
                << "predictor / extrap2 channels are skipped" << std::endl;
      said = true;
    }
  };
  // the last t2_nfv channels: the sc rays' intensities (their number is not
  // known before the first sweep, Time2RstAccept)
  const int nfv0 = Time2RstNchWant() - CleanRstNfv() - (t2_rstd3 ? 3 : 0) + t2_rskip;
  if (clean && cln_cssc && t2_rstv >= 6 && ch >= nfv0) {   // (cln_fh)
    if (cln_fh.extent_int(0) != nmb) {
      Kokkos::realloc(cln_fh, nmb, 12, w.extent_int(1), w.extent_int(2), w.extent_int(3));
    }
    DeepCopyAcross(Kokkos::subview(cln_fh, mb, ch - nfv0, AL, AL, AL), w);
    return;
  }
  if (clean && cln_cssc && t2_rstv >= 5 && ch >= nfv0) {   // (version 5, CleanFvRstCopy)
    DvceArray5D<Real> t("m1_t2_fvrd", nmb, 1, w.extent_int(1), w.extent_int(2),
                        w.extent_int(3));
    DeepCopyAcross(Kokkos::subview(t, AL, 0, AL, AL, AL), w);
    CleanFvRstCopy(t, nfv0 - ch, nmb, 0, false);
    return;
  }
  if (clean && (cln_cssc || cln_scmem) && ch >= nfv0) {
    auto &fv = CleanRstFv();
    if (fv.extent_int(0) != nmb || fv.extent_int(1) != t2_nfv) {
      Kokkos::realloc(fv, nmb, t2_nfv, w.extent_int(1), w.extent_int(2), w.extent_int(3));
    }
    DeepCopyAcross(Kokkos::subview(fv, mb, ch - nfv0, AL, AL, AL), w);
    return;
  }
  if (h2 && ch < M1_T2_NK) {
    DeepCopyAcross(Kokkos::subview(t2k1, mb, ch, AL, AL, AL), w);
    return;
  }
  if (h2) ch -= M1_T2_NK;
  if (h2 && impl_pred && !(clean && cln_not2h && t2_rstv >= 4)) {
    const int np = (impl_pord == 2) ? 5 : 3;
    if (ch < np) {
      if (!(clean && cln_not2h)) {
        DeepCopyAcross(Kokkos::subview(ipred2, mb, ch, AL, AL, AL), w);
      } else {
        skip();
      }
      return;
    }
    ch -= np;
  }
  if (h2 && vet_sc) {
    if (ch < M1_T2_NVET) {
      DeepCopyAcross(Kokkos::subview(vet_prev, mb, ch, AL, AL, AL), w);
      return;
    }
    ch -= M1_T2_NVET;
    if (vsc_every > 1) {
      if (ch < 2*M1_T2_NVET) {
        if (ch < M1_T2_NVET) {
          DeepCopyAcross(Kokkos::subview(vsc_d0, mb, ch, AL, AL, AL), w);
        } else {
          DeepCopyAcross(Kokkos::subview(vsc_d1, mb, ch - M1_T2_NVET, AL, AL, AL), w);
        }
        return;
      }
      ch -= 2*M1_T2_NVET;
    }
  }
  if (force_ref == M1_FREF_CLEAN || force_ref == M1_FREF_WB_MD) {
    const int nr = (force_ref == M1_FREF_WB_MD) ? 9 : 3 + cln_lzt;
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
  if (clean && cln_extrap > 0) {
    if (ch < 2*nq) {
      auto &mx = (ch < nq) ? cln_mn : cln_mc;
      if (mx.extent_int(0) != nmb) {
        Kokkos::realloc(mx, nmb, CLN_N, w.extent_int(1), w.extent_int(2),
                        w.extent_int(3));
      }
      DeepCopyAcross(Kokkos::subview(mx, mb, CleanRstSlot(ch%nq), AL, AL, AL), w);
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
  if (t2_rskip > 0) {   // (rotwb2-1121: the file's ring with the ring off, not read)
    if (ch < t2_rskip) return;
    ch -= t2_rskip;
  }
  if (clean && cln_ring && t2_rstv >= 4) {   // (version 4, CleanRingRstCopy)
    if (ch < CleanRingRstN()) {
      DvceArray5D<Real> t("m1_t2_ringrd", nmb, 1, w.extent_int(1), w.extent_int(2),
                          w.extent_int(3));
      DeepCopyAcross(Kokkos::subview(t, AL, 0, AL, AL, AL), w);
      CleanRingRstCopy(t, -ch, nmb, 0, false, true);
      return;
    }
    ch -= CleanRingRstN() - vra_k*nq;   // (the flag slab, as in the old layout)
  }
  if (clean && cln_ring) {
    const int nr = (vra_k - 1)*nq;
    if (ch < nr + nq) {
      // (the interior cells of the carried fields, CleanRingMoments)
      auto &rx = (ch < nr) ? cln_rr : cln_rpd;
      vra_gen += 1;
      auto &id = pmy_pack->pmesh->mb_indcs;
      const int is = id.is, js = id.js, ks = id.ks;
      Kokkos::Array<int, CLN_N> fm;
      const int nf = CleanRingFields(fm.data(), true);
      const int f = fm[CleanRstSlot((ch < nr) ? ch%nq : ch - nr)];
      const int nx = (ch < nr) ? (vra_k - 1)*nf : nf;
      if (rx.extent_int(0) != nmb) {
        Kokkos::realloc(rx, nmb, nx, id.nx3, id.nx2, id.nx1);
      }
      if (f < 0) return;
      const int q = (ch < nr) ? (ch/nq)*nf + f : f;
      DvceArray4D<Real> t("m1_t2_ringrd", nmb, w.extent_int(1), w.extent_int(2),
                          w.extent_int(3));
      DeepCopyAcross(t, w);
      auto rx_ = rx;
      par_for("m1_t2_ringset", DevExeSpace(), 0, nmb - 1, 0, id.nx3 - 1, 0, id.nx2 - 1, 0,
              id.nx1 - 1,
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
  if (clean && cln_fsolve && cln_fkin) {   // (CleanInit copies it to the face arrays)
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
  // extrap2 (Time2RstPack)
  // (clean_skip_time2_hist: no extrap2 history kept, x2_ok stays false)
  if (h2 && vet_col && t2_vcmode == 4 && clean && cln_not2h) skip();
  if (h2 && vet_col && t2_vcmode == 4 && !(clean && cln_not2h)) {
    if (ch < M1_X2_NTT) {
      DeepCopyAcross(Kokkos::subview(x2_tt, mb, ch, AL, AL, AL), w);
    } else if (ch < M1_X2_NTT + M1_X2_NHR) {
      DeepCopyAcross(Kokkos::subview(x2_hr, mb, ch - M1_X2_NTT, AL, AL, AL), w);
    } else if (ch == M1_X2_NTT + M1_X2_NHR) {
      for (int m = 0; m < nmb; ++m) {   // the q column sits at i = 0
        auto hq = Kokkos::create_mirror_view(Kokkos::subview(x2_q, m, AL, AL));
        for (int k = 0; k < static_cast<int>(hq.extent(0)); ++k) {
          for (int j = 0; j < static_cast<int>(hq.extent(1)); ++j) {hq(k,j) = w(m,k,j,0);}
        }
        Kokkos::deep_copy(Kokkos::subview(x2_q, m, AL, AL), hq);
      }
    } else if (ch == M1_X2_NTT + M1_X2_NHR + 1) {
      x2_t = w(0,0,0,0);
      x2_ok = (x2_t > -1.0e299);
    } else if (ch == M1_X2_NTT + M1_X2_NHR + 2) {
      x2_alpha = w(0,0,0,0);
    } else {
      x2_dts = w(0,0,0,0);
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationM1::Time2AutoBeBegin
//! \brief autobe-1029 <rad_m1>/time2_auto_be (rad_m1.hpp): called by the Driver at the
//! cycle start (after VetGdRotStep).  In a skip run the step is put in exactly the state
//! the stage fallback leaves for its redo (Time2Restore): no valid slope (t2_ok = false,
//! so the Driver runs the Heun hydro alone and then the BESTORE backward-Euler solve, and
//! FrefWaccOn sees a BE step), cold predictor, and the dt_retry arm of a time2_reject
//! redo.  A retried cycle calls it again with the same (unchanged) run state.

bool RadiationM1::Time2AutoBeBegin() {
  t2_abe_now = false;
  // t2_beonly: every step (the first too) is a skip step
  const bool bo = t2_beonly && (time_scheme == M1_TIME_HESDIRK2);
  if (!bo && (t2_abe <= 0 || time_scheme != M1_TIME_HESDIRK2 || t2_abe_left <= 0)) {
    return false;
  }
  if (!bo && !Time2Active()) return false;   // no slope anyway (the BE path follows)
  t2_abe_now = true;
  t2_ok = false;
  t2_fail = false;
  t2_solve = M1_T2S_NONE;
  pred_ok = false;
  pred2_ok = false;
  // as the redo of a rejected stage under time2_reject = retry (arms dt_retry_t2_be)
  if (t2_reject == 1 && dtr_on) {dtr_t2hit = true;}
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::Time2AutoBeEnd
//! \brief an ACCEPTED cycle of a hesdirk2 run.  kind 0 = stage step, 1 =
//! stage fallback + BE redo, 2 = time2_auto_be step, 3 = other BE step (no slope),
//! 4 = a discarded dt_retry attempt (timing only).  tw_a / tw_b: fenced wall seconds of
//! (the attempt | the step) and (the redo).

void RadiationM1::Time2AutoBeEnd(int kind, Real tw_a, Real tw_b) {
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
//! \fn void RadiationM1::Time2Report

void RadiationM1::Time2Report() {
  if (time_scheme != M1_TIME_HESDIRK2) return;
  if (global_variable::my_rank != 0) return;
  std::cout << "<rad_m1> time_scheme=" << (t2_beonly ? "hesdirk2(be_only)" : "hesdirk2")
            << ": stage steps=" << t2_nstep
            << " backward-Euler steps=" << t2_nbe << " stage fallbacks=" << t2_nfall
            << " vet clips=" << t2_nclip;
  if (t2_k1reset) {
    std::cout << " k1_reset cells(rank 0)=" << t2_k1reset_n << " faces=" << t2_k1reset_nf;
    if (t2_k1pos) {std::cout << " (pos rule cells " << t2_k1pos_n << ")";}
  }
  if (t2_k1rd > 0.0) {
    std::cout << " k1_redist(th=" << t2_k1rd << ") cells(rank 0)=" << t2_k1rd_n
              << " unpaid=" << t2_k1rd_un;
  }
  if (t2_k1stiff > 0.0) {
    std::cout << " k1_stiff(s=" << t2_k1stiff << ") cells(rank 0)=" << t2_k1stiff_n;
  }
  if (t2_k1xlim > 0.0 || t2_k1xpat) {
    std::cout << " k1_xlim(th=" << t2_k1xlim
              << (t2_k1xpos ? " xpat=pos" : (t2_k1xpat ? " xpat" : ""))
              << ") cell-stages(rank 0)=" << t2_k1xlim_n;
  }
  if (t2_ffct) {   // step 1 (global sums)
    std::cout << " face_fct calls=" << t2_ffct_calls << " skipped(no face flux)="
              << t2_ffct_nop << " K1 flux from U^n=" << t2_ffct_fill << " limited cells="
              << t2_ffct_cells << " max rounds="
              << t2_ffct_rmax << " withheld E x vol=" << t2_ffct_de << " left (Q<=0)="
              << t2_ffct_left << " emptied (E < " << t2_ffct_emin << " E^n) cells="
              << t2_ffct_nemp << " min E/E^n="
              << t2_ffct_qmin << " (rank 0: rounds=" << zl_rounds[0] << " R exchanges="
              << zl_halos[0];
    if (zl_local > 0) {std::cout << " local sub-rounds=" << zl_sub[0];}
    if (zl_prof) {std::cout << " round-loop wall=" << zl_wall[0] << " s";}
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
  if (t2_start == 1) {
    std::cout << " time2_start=slope probes=" << t2_nprobe << " (frac " << t2_start_frac
              << ")";
  }
  if (vgd_rot_slope == 3) {
    std::cout << " vet_gd_rotate_slope=correct solves=" << vgd_rot_ncorr << " failed="
              << vgd_rot_ncfail << " Picard passes=" << vgd_rot_ncpass << " (frac "
              << vgd_rot_cfrac << ", reuse "
              << (vgd_rot_creuse ? "on" : "off") << ")";
  }
  if (t2_accfb) {std::cout << " accepted-with-fallback stages=" << t2_accfb_n;}
  if (t2_adm_pl) {std::cout << " time2_adm_poslim admitted stages=" << t2_adm_pl_n;}
  if (t2_flr_rho > 0.0) {
    std::cout << " time2_floor_eint(rho<" << t2_flr_rho << ") cell-stages(rank 0)="
              << t2_flr_n;
  }
  // VetGdAlignDt calls that shortened dt (NewTimeStep calls, incl. the start-up ones)
  if (vgd_rot_align) {std::cout << " rotate_align shortened dt=" << vgd_rot_nalign;}
  if (cln_ring) {
    std::cout << " clean_rot_ring guard cells=" << cln_rguard_sum
              << " off-boundary shifts=" << cln_roff_sum;
  }
  std::cout << std::endl;
  if (t2_xsp) {   // time2_exch_split (rank 0)
    std::cout << "<rad_m1> time2_exch_split=psi (rank 0): stages " << t2_xsp_calls
              << " cells corrected " << t2_xsp_cells << " sum |dEN| x V moved to E "
              << t2_xsp_de << std::endl;
  }
  if (t2_iex) {   // time2_inc_exch (rank 0)
    std::cout << "<rad_m1> time2_inc_exch=relax (rank 0): sum (eint given x V, E paid) "
              << t2_iex_de << " cell-solves " << t2_iex_cells
              << " cells with E_old <= 0 after paying " << t2_iex_neg << std::endl;
  }
  if (t2_ike) {   // (Time2Report runs on rank 0 only: rank-0 values)
    std::cout << "<rad_m1> time2_inc_energy=rad (rank 0): solves " << t2_ike_calls
              << " sum (KE remainder x V) to E " << t2_ike_de << " sum |.| x V "
              << t2_ike_dp << " cells with E_old <= 0 and remainder > 0 " << t2_ike_neg
              << " cells with a nonzero remainder " << t2_ike_cells << std::endl;
  }
  if (t2_abe > 0) {
    std::cout << "<rad_m1> time2_auto_be=" << t2_abe << ": skipped-stage BE steps="
              << t2_abe_n << " stage retries=" << t2_abe_ntry << " (succeeded "
              << t2_abe_nok << ")" << std::endl;
  }
  // the wall-time split (rank 0, fenced; accepted cycles by kind)
  std::cout << "<rad_m1> hesdirk2 wall split (rank 0 s): clean stage steps " << t2w_clean
            << " | failed stage attempts " << t2w_fail << " | their BE redo " << t2w_redo
            << " | auto_be steps " << t2w_abe << " | other BE steps " << t2w_be
            << " | discarded dt_retry attempts " << t2w_retry << std::endl;
}

}  // namespace radm1

namespace radm1 {
//----------------------------------------------------------------------------------------
//! \fn RadiationM1::RstFingerprint (DIAGNOSTIC <rad_m1>/dbg_rst_fp = N)
//! Sum and sum of squares of every persistent array (per component of the 5-D ones), all
//! ranks, printed by rank 0: run A (straight) and run B (restarted) are diffed line by
//! line at the same cycle to find the state a restart does not carry.  Moves no state.
//! Also the vet_gd source and the clean scheme's arrays and values; dbg_rst_fp_cycle adds
//! fingerprints inside the first clean solve of a cycle, dbg_rst_fp_dump writes vgd_cs.
namespace {
template <class V>
void M1Fp(const char *nm, const V &v, int cycle) {
  if (v.span() == 0) return;
  auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);
  const int nc = (V::rank == 5) ? static_cast<int>(v.extent(1)) : 1;
  const std::size_t per = h.span()/(V::rank == 5 ? v.extent(1)*v.extent(0) : 1);
  for (int c = 0; c < nc; ++c) {
    double s1 = 0.0, s2 = 0.0;
    if (V::rank == 5) {
      for (std::size_t m = 0; m < v.extent(0); ++m) {
        const Real *p = h.data() + (m*v.extent(1) + c)*per;
        for (std::size_t q = 0; q < per; ++q) {s1 += p[q]; s2 += p[q]*p[q];}
      }
    } else {
      const Real *p = h.data();
      for (std::size_t q = 0; q < h.span(); ++q) {s1 += p[q]; s2 += p[q]*p[q];}
    }
    double g[2] = {s1, s2};
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, g, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (global_variable::my_rank == 0) {
      std::printf("RSTFP cyc %d %s[%d] %.17e %.17e\n", cycle, nm, c, g[0], g[1]);
    }
  }
}
}  // namespace
void RadiationM1::RstFingerprint(const char *tag) {
  const int cy = pmy_pack->pmesh->ncycle;
  if (global_variable::my_rank == 0) {std::printf("RSTFP cyc %d tag %s\n", cy, tag);}
  M1Fp("u0", u0, cy); M1Fp("opac", opac, cy); M1Fp("iw", iw, cy);
  M1Fp("t2k1", t2k1, cy); M1Fp("t2k2", t2k2, cy); M1Fp("ipred", ipred, cy);
  M1Fp("ipred2", ipred2, cy); M1Fp("f0x1", f0x1, cy); M1Fp("f0x2", f0x2, cy);
  M1Fp("f0x3", f0x3, cy); M1Fp("tau_ten", tau_ten, cy); M1Fp("vet_cell", vet_cell, cy);
  M1Fp("x2_tt", x2_tt, cy); M1Fp("x2_hr", x2_hr, cy); M1Fp("x2_q", x2_q, cy);
  M1Fp("ctr_mem", ctr_mem, cy); M1Fp("csk1", csk1, cy); M1Fp("ifw", ifw, cy);
  M1Fp("ifw2", ifw2, cy); M1Fp("ifw3", ifw3, cy); M1Fp("thx2", thx2, cy);
  M1Fp("thx3", thx3, cy); M1Fp("ktd", ktd, cy); M1Fp("vgd_hr", vgd_hr, cy);
  M1Fp("vgd_hr_c", vgd_hr_c, cy); M1Fp("t2f1", t2f1, cy); M1Fp("t2f2", t2f2, cy);
  M1Fp("t2f3", t2f3, cy); M1Fp("vet_prev", vet_prev, cy);
  // the vet_gd source and the clean scheme's state
  M1Fp("vgd_cs", vgd_cs, cy); M1Fp("vlat_cs", vlat_cs, cy);
  FluidRef fl = FluidRef::Get(pmy_pack);
  if (fl.on) {M1Fp("fl_u0", fl.u0, cy); M1Fp("fl_w0", fl.w0, cy);}
  M1Fp("cln_s", cln_s, cy); M1Fp("cln_v", cln_v, cy); M1Fp("cln_w", cln_w.s, cy);
  M1Fp("cln_wo", cln_w.o, cy);
  M1Fp("cln_fc1", cln_fc1, cy); M1Fp("cln_fc2", cln_fc2, cy);
  M1Fp("cln_fc3", cln_fc3, cy);
  M1Fp("cln_fref", cln_fref, cy); M1Fp("cln_mn", cln_mn, cy); M1Fp("cln_mp", cln_mp, cy);
  M1Fp("cln_mc", cln_mc, cy); M1Fp("cln_m1", cln_m1, cy); M1Fp("cln_fvi", cln_fvi, cy);
  M1Fp("cln_fvd", cln_fvd, cy);
  M1Fp("cln_opb", cln_opb, cy); M1Fp("cln_gn", cln_gn.r, cy);
  M1Fp("cln_gn.a", cln_gn.a, cy);
  M1Fp("cln_rl", cln_rl, cy); M1Fp("cln_rr", cln_rr, cy); M1Fp("cln_rpd", cln_rpd, cy);
  M1Fp("cln_rnew", cln_rnew, cy);
  if (dbg_rst_fp_dump && vgd_cs.span() > 0) {
    auto hc = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), vgd_cs);
    char fn[256];
    std::snprintf(fn, sizeof(fn), "vgdcs_%s_r%d.bin", tag, global_variable::my_rank);
    FILE *fo = std::fopen(fn, "wb");
    int ext[5];
    for (int r = 0; r < 5; ++r) ext[r] = static_cast<int>(hc.extent(r));
    std::fwrite(ext, sizeof(int), 5, fo);
    std::fwrite(hc.data(), sizeof(Real), hc.span(), fo);
    std::fclose(fo);
  }
  if (global_variable::my_rank == 0) {
    std::printf("RSTFP cyc %d scalars pred_ok %d pred_dt %.17e pred2_ok %d "
                "pred2_dt %.17e t2_ok %d t2_dtprev %.17e t2_vprev %d x2_ok %d x2_t %.17e "
                "abe %d %d\n", cy, pred_ok ? 1 : 0, pred_dt, pred2_ok ? 1 : 0, pred2_dt,
                t2_ok ? 1 : 0, t2_dtprev, t2_vprev ? 1 : 0, x2_ok ? 1 : 0, x2_t,
                t2_abe_len, t2_abe_left);
    std::fflush(stdout);
  }  if (global_variable::my_rank == 0) {
    std::printf("RSTFP cyc %d clean init %d faces %d havec %d tc %.17e "
                "dtb %.17e fe %.17e fv %.17e\n", cy, cln_init ? 1 : 0,
                cln_faces_set ? 1 : 0, cln_havec ? 1 : 0, cln_tc,
                cln_dtb, cln_fe, cln_fv);
    std::printf("RSTFP cyc %d extrap tn %.17e tp %.17e tm %.17e t1 %.17e "
                "have %d %d %d %d cand_be %d\n", cy, cln_tn, cln_tp, cln_tm, cln_t1,
                cln_have_n ? 1 : 0, cln_have_p ? 1 : 0, cln_have_c ? 1 : 0,
                cln_have_1 ? 1 : 0, cln_cand_be ? 1 : 0);
    std::printf("RSTFP cyc %d fv n %d icut %d scut %d\n", cy, cln_fvn, vlat_icut,
                vgd_scut);
  }

}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationM1::VsrcNoes
//! \brief vet_source_noesrc in force for the current solve: every backward-Euler solve
//! (time_scheme = be; under hesdirk2 the BESTORE / fallback / ic_rad solves unless
//! vet_source_noesrc_h2be = false); never a hesdirk2 STAGE solve (its build reads iw EN,
//! where time2_vet_col = predict puts the predicted state)

bool RadiationM1::VsrcNoes() const {
  if (!vsrc_noes) return false;
  if (time_scheme == M1_TIME_BE) return true;
  return vsrc_noes_h2be && (t2_solve != M1_T2S_STAGE1) && (t2_solve != M1_T2S_STAGE2);
}

}  // namespace radm1
