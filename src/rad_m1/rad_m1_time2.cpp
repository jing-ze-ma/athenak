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
#include "coordinates/coordinates.hpp"
#include "driver/driver.hpp"
#include "hydro/hydro.hpp"
#include "hydro/hydro_am.hpp"
#include "eos/eos.hpp"
#include "rad_m1/m1_fluid.hpp"
#include "utils/deep_copy_across.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_parfor.hpp"
#include "rad_m1/rad_m1_implicit.hpp"
#include "rad_m1/rad_m1_opacity.hpp"

namespace radm1 {

using deep_copy_across::DeepCopyAcross;

namespace {
const Real kT2G = 1.0 - 1.0/std::sqrt(2.0);
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
  if (time_scheme == M1_TIME_BE) return;

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
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<rad_m1>/time2_tableau = '" << tb
                  << "' is not a choice (hesdirk2 | trbdf2)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    if (t2_tab == 0) {
      const Real b2 = g/(2.0*(1.0 - g));
      t2_a21 = 1.0 - g;
      t2_w0 = 0.5;
      t2_w1 = 0.5;
      t2_c31 = 0.5*g;
      t2_c32 = b2 - 0.5*g;
    } else {
      const Real w = 0.25*std::sqrt(2.0);
      t2_a21 = g;
      t2_w0 = 0.25;
      t2_w1 = 0.75;
      t2_c31 = w - 0.25*g;
      t2_c32 = w - 0.25*g;
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
  // open-ot-1013 time2_enth_be: with time2_enth_vel = central the stage solves build a_f
  // as the mean of the two cells, but the backward-Euler solves of the same run (the start
  // step, the time2_start slope probe, the rotation probe, the step after a failed stage)
  // take implicit_enthalpy's own a_f (plm: the 4-point plm mean).  That changes the
  // operator for one step of length dt, an O(dt) term of the T ladders (rw_p10t100 E /
  // eint 1.78 / 1.84 at l 16-32).  true: those solves take the stage rule too.
  // fin4-1014: DEFAULT true on a fresh run; a restart whose file lacks the key keeps false
  // (FIN3, bitwise) -- the same rule as implicit_cs2_vimp_npass (fresh 3, restart 0).
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
  // time2_vet_sc (torder-1010, read only when named; default lag): when the vet_sc
  // tensor of the two stage solves is built.  Both stages sit at t^{n+1} (c = 1), so
  // D(U^n) is an O(dt) lag and the step is FIRST order in time wherever the tensor
  // changes (torder_1010: pulse_k12.8 T 1.05 1.03 1.02; eddington or the extrapolated
  // D* give 2.00).  The modes are those of time2_vet_col (Time2VetColAt):
  //   lag     : D(U^n) for both stages (time2_vet_extrap = true: D* extrapolated)
  //   predict : ONE formal solution per step at P = (stage-1 start) + dt K1, i.e. at
  //             t^{n+1} + O(dt^2); same cost as lag
  //   rebuild : predict for stage 1, a second formal solution at Y1 for stage 2
  // unify-1012: DEFAULT predict on a fresh run (as time2_vet_col on the wedge), unless
  // time2_vet_extrap = true or vet_sc_every > 1 (which predict excludes: lag, silently);
  // a restart whose file lacks the key keeps lag (echoed).  torder-1010: lag is first
  // order in time wherever the tensor changes (pulse_k12.8 T 1.05), predict 2.0.
  // unify-sc-gd-1012 (#13): ONE key time2_vet for both formal solutions (lag | predict |
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
                                         global_variable::restart_run ? "lag" : "predict");
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
  // unify-1012: default true on every mesh (was the wedge only): the Cartesian vet_sc
  // ladders need it as much (torder-1010 cenPvs: rw_t10a F1 1.58 2.32 1.97 vs 1.07 1.04
  // 1.02 without); a restart whose file lacks the key keeps false (echoed)
  t2_fvnew = pin->GetOrAddBoolean("rad_m1", "time2_vstage",
                                  !global_variable::restart_run);
  t2_vcmode = 0;
  if (vet_col) {
    const bool vdef = sph_geom && !global_variable::restart_run;
    const std::string kvc = M1KeyName(pin, "time2_vet", {"time2_vet_col"});
    std::string vm = (sph_geom || pin->DoesParameterExist("rad_m1", kvc)) ?
        pin->GetOrAddString("rad_m1", kvc, vdef ? "predict" : "lag") :
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
      t2_vcmode = 4;   // rbfast-1021 (rad_m1.hpp Time2VetColExtrap2)
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
  // unify-sc-gd-1012 (#3): vet_sc reads the same time2_vet_source (its builds take the
  // VET source rule of vet_col / vet_gd, rad_m1_vet.cpp).  A restart whose file lacks the
  // key keeps vet_sc's old source of the predicted state (state); the old spelling
  // time2_vet_col_source is vet_col's own (a vet_sc file may name it, inert there).
  if (vet_sc) {
    Time2VetSourceInit(pin, t2_vsmode == 1 || t2_vsmode == 2, false);
  }
  // time2_pred_cap (predfix-1011): the predicted state of time2_vet_sc / time2_vet_col =
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
  t2_lin_tol = -1.0;
  if (pin->DoesParameterExist("rad_m1", "time2_lin_tol")) {
    t2_lin_tol = pin->GetReal("rad_m1", "time2_lin_tol");
  }
  t2_lin_fac = pin->GetOrAddReal("rad_m1", "time2_lin_tol_fac", 1.0);
  if (!(t2_lin_fac > 0.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<rad_m1>/time2_lin_tol_fac must be > 0" << std::endl;
    std::exit(EXIT_FAILURE);
  }
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
  // time2_start / time2_start_frac (h2order-1011, rad_m1.hpp; read only when named):
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
  // dbg_t2_admiss = N (DEBUG, default 0 = off, read only when named so the restart echo
  // and every result stay bitwise): for the first N non-admissible hesdirk2 stages print
  // which quantity went <= 0 (solved E, solved T, written-back gas eint) and, for each,
  // the cell of its minimum with the stage OLD vector (E, eint) and the stage start E
  t2_dbg_adm = 0;
  t2_dbg_adm_n = 0;
  if (pin->DoesParameterExist("rad_m1", "dbg_t2_admiss")) {
    t2_dbg_adm = pin->GetInteger("rad_m1", "dbg_t2_admiss");
  }
  // t2retry-1020: <rad_m1>/time2_reject = be | retry (read only when named; default be =
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

  // pgrst-1021: <rad_m1>/time2_k1_reset (see rad_m1.hpp)
  t2_k1reset = false;
  if (pin->DoesParameterExist("rad_m1", "time2_k1_reset")) {
    const std::string kr = pin->GetString("rad_m1", "time2_k1_reset");
    if (kr == "pos") {
      t2_k1reset = true;
      t2_k1pos = true;
    } else {
      t2_k1reset = pin->GetBoolean("rad_m1", "time2_k1_reset");
    }
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
                << std::endl << "<rad_m1>/time2_k1_reset_grow must be 0 or 1" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  t2_ys_ok = false;
  // t2pos-1020: <rad_m1>/time2_posfix = N, time2_posfix_grow (read only when named;
  // absent = off, every result bitwise): the local stage-1 positivity re-solve
  // (Time2PosFix).  Not wired for MHD (the field put-back); <hydro>/am_vector puts the
  // spin reservoir back with the gas (posfix-amvec-1022).
  t2_pfix = 0;
  if (pin->DoesParameterExist("rad_m1", "time2_posfix")) {
    t2_pfix = pin->GetInteger("rad_m1", "time2_posfix");
    if (t2_pfix < 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/time2_posfix must be >= 0" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // accfb-default-1022: time2_accept_fallback is DEFAULT true under hesdirk2 (this
  // reader returns above for be, where the key is unused); false only when named so
  t2_accfb = true;
  if (pin->DoesParameterExist("rad_m1", "time2_accept_fallback")) {
    t2_accfb = pin->GetBoolean("rad_m1", "time2_accept_fallback");
  }
  if (pin->DoesParameterExist("rad_m1", "restart_pack_per_mb")) {
    rst_pack_mb = pin->GetBoolean("rad_m1", "restart_pack_per_mb");
  }
  if (pin->DoesParameterExist("rad_m1", "time2_posfix_host")) {
    t2_pfix_host = pin->GetBoolean("rad_m1", "time2_posfix_host");
  }
  if (pin->DoesParameterExist("rad_m1", "time2_posfix_grow")) {
    t2_pfix_grow = pin->GetInteger("rad_m1", "time2_posfix_grow");
    if (t2_pfix_grow < 0 || t2_pfix_grow > 2) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/time2_posfix_grow must be 0, 1 or 2 (the ghost "
                << "depth)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  if (t2_pfix > 0) {
    FluidRef flp = FluidRef::Get(pmy_pack);
    if (flp.on && flp.mhd) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/time2_posfix is not wired for MHD" << std::endl;
      std::exit(EXIT_FAILURE);
    }
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
  Kokkos::realloc(t2f1, nmb, n3, n2, n1+1);
  if (trans_on) {
    Kokkos::realloc(t2f2, nmb, n3, n2+1, n1);
    if (trans_x3) {Kokkos::realloc(t2f3, nmb, n3+1, n2, n1);}
  }
  if (impl_pred) {
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
    Kokkos::realloc(vet_now, nmb, (t2_vcmode == 0 || t2_vcmode >= 3) ? 2 : 6, n3, n2, n1);
    Kokkos::realloc(vet_opac, nmb, vscat ? (M1_NOPAC + 1) : M1_NOPAC, n3, n2,
                    n1);
    if (t2_vcmode == 4) {
      Kokkos::realloc(x2_tt, nmb, M1_X2_NTT, n3, n2, n1);
      Kokkos::realloc(x2_hr, nmb, M1_X2_NHR, n3, n2, n1);
      Kokkos::realloc(x2_q, nmb, n3, n2);
      Kokkos::deep_copy(x2_tt, 0.0);
      Kokkos::deep_copy(x2_hr, 0.0);
      Kokkos::deep_copy(x2_q, 0.0);
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
  if (t2_pfix > 0) {
    Kokkos::realloc(t2pm, nmb, n3, n2, n1);
    Kokkos::deep_copy(t2pm, 0.0);
  }
  t2_ok = false;
  if (global_variable::my_rank == 0) {
    if (t2_pfix > 0) {
      std::cout << "<rad_m1> time2_posfix=" << t2_pfix << " time2_posfix_grow="
                << t2_pfix_grow << ": local stage-1 positivity re-solve"
                << (t2_pfix_host ? " (start copies in pinned host memory)" : "")
                << std::endl;
    }
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

void RadiationM1::Time2FormStage(int stage, Real dt) {
  const Real g = kT2G;
  const Real w0 = t2_w0, w1 = t2_w1;
  const bool heun = (t2_tab == 0);   // hesdirk2: the original expression, bitwise
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
    // pgrst-1021 time2_k1_reset: no K1 where E was overwritten since the stage-loop start
    if (t2_k1reset && t2_ys_ok && t2ys.extent_int(0) >= pmy_pack->nmb_thispack) {
      auto &id = pmy_pack->pmesh->mb_indcs;
      const int is = id.is, ie = id.ie, js = id.js, je = id.je, ks = id.ks, ke = id.ke;
      const int nmbr = pmy_pack->nmb_thispack - 1;
      auto ys_ = t2ys;
      auto k1r_ = t2k1;
      const int ncell = M1_T2_F1;   // E, gas momentum (3), gas energy
      auto &vol = pmy_pack->pcoord->volume;
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
          sde += (k1r_(m,M1_T2_E,k,j,i) + k1r_(m,M1_T2_EN,k,j,i))*vol(m,k,j,i);
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
    t2_ys_ok = false;
    // t2pos-1020 time2_posfix: the stage-1 start of the gas (and of the closure memory)
    // for a re-solve, and a clean mask
    if (t2_pfix > 0) {
      FluidRef flp = FluidRef::Get(pmy_pack);
      if (flp.on && t2_pfix_host) {
        // gpumem-1022 time2_posfix_host: the copy in pinned host memory (stream-ordered)
        if (t2gsv_h.extent(0) != flp.u0.extent(0) || t2gsv_h.extent(1) != flp.u0.extent(1)
            || t2gsv_h.extent(2) != flp.u0.extent(2)
            || t2gsv_h.extent(3) != flp.u0.extent(3)
            || t2gsv_h.extent(4) != flp.u0.extent(4)) {
          t2gsv_h = decltype(t2gsv_h)(Kokkos::view_alloc(std::string("t2gsv_h"),
                                      Kokkos::WithoutInitializing), flp.u0.extent(0),
                                      flp.u0.extent(1), flp.u0.extent(2),
                                      flp.u0.extent(3), flp.u0.extent(4));
        }
        Kokkos::deep_copy(DevExeSpace(), t2gsv_h, flp.u0);
      } else if (flp.on) {
        if (t2gsv.extent(0) != flp.u0.extent(0) || t2gsv.extent(1) != flp.u0.extent(1) ||
            t2gsv.extent(2) != flp.u0.extent(2) || t2gsv.extent(3) != flp.u0.extent(3) ||
            t2gsv.extent(4) != flp.u0.extent(4)) {
          Kokkos::realloc(t2gsv, flp.u0.extent(0), flp.u0.extent(1), flp.u0.extent(2),
                          flp.u0.extent(3), flp.u0.extent(4));
        }
        Kokkos::deep_copy(DevExeSpace(), t2gsv, flp.u0);
      }
      // posfix-amvec-1022 <hydro>/am_vector: the spin reservoir at the stage-1 start (the
      // solve adds its increment aminc and the write-back's radial L to it)
      if (am_vec) {
        auto &spn = hydro::AmModuleOf(pmy_pack)->am_spin;
        if (t2ssv.extent(0) != spn.extent(0) || t2ssv.extent(1) != spn.extent(1) ||
            t2ssv.extent(2) != spn.extent(2) || t2ssv.extent(3) != spn.extent(3)) {
          Kokkos::realloc(t2ssv, spn.extent(0), spn.extent(1), spn.extent(2),
                          spn.extent(3));
        }
        Kokkos::deep_copy(DevExeSpace(), t2ssv, spn);
      }
      if (ctr_mem.extent(0) > 0 && t2_pfix_host) {
        if (t2csv_h.extent(0) != ctr_mem.extent(0)) {
          t2csv_h = decltype(t2csv_h)(Kokkos::view_alloc(std::string("t2csv_h"),
                                      Kokkos::WithoutInitializing), ctr_mem.extent(0),
                                      ctr_mem.extent(1), ctr_mem.extent(2),
                                      ctr_mem.extent(3), ctr_mem.extent(4));
        }
        Kokkos::deep_copy(DevExeSpace(), t2csv_h, ctr_mem);
      } else if (ctr_mem.extent(0) > 0) {
        if (t2csv.extent(0) != ctr_mem.extent(0)) {
          Kokkos::realloc(t2csv, ctr_mem.extent(0), ctr_mem.extent(1), ctr_mem.extent(2),
                          ctr_mem.extent(3), ctr_mem.extent(4));
        }
        Kokkos::deep_copy(DevExeSpace(), t2csv, ctr_mem);
      }
      Kokkos::deep_copy(DevExeSpace(), t2pm, 0.0);
      t2_pfix_n = 0;
      t2_pfix_any = false;
      t2_pfix_want = false;
      t2_pfix_redo = false;
    }
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
  // t2pos-1020 time2_posfix: where stage 1 dropped its increment t2_a21 dt K1, stage 2
  // takes it back with the weight of Y1 in its start state (heun: 1/2), so that the step
  // is U^n + dt (b1 K1 + b2 K2 + g K3) with the tableau's b everywhere
  const bool pcomp = two && (t2_pfix > 0) && t2_pfix_any;
  const Real ccomp = t2_a21*dt*(heun ? 0.5 : w0);
  if (pcomp) {
    auto pm_ = t2pm;
    const bool d2 = (n2 > 1), d3 = (n3 > 1);
    par_for("m1_t2_pcomp", DevExeSpace(), 0, nmb1, 0, M1_T2_NK-1, 0, n3-1, 0, n2-1,
            0, n1-1,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
      Real w = pm_(m,k,j,i);
      if (n == M1_T2_F1 && i > 0) {w = fmax(w, pm_(m,k,j,i-1));}
      if (n == M1_T2_F1+1) {w = d2 ? ((j > 0) ? fmax(w, pm_(m,k,j-1,i)) : w) : 0.0;}
      if (n == M1_T2_F1+2) {w = d3 ? ((k > 0) ? fmax(w, pm_(m,k-1,j,i)) : w) : 0.0;}
      if (w > 0.0) {in_(m,n,k,j,i) += ccomp*k1_(m,n,k,j,i);}
    });
  }
  // cs-open-1018 <hydro>/am_vector: the spin reservoir's increment from its slopes
  if (am_vec && aminc.extent_int(0) > 0) {
    auto a1_ = amk1;
    auto a2_ = amk2;
    auto ai_ = aminc;
    // posfix-amvec-1022: the spin's share of the stage-2 compensation (cell slot, as the
    // gas momentum of t2inc above)
    auto pma_ = t2pm;
    const bool apc = pcomp && (t2pm.extent(0) >= aminc.extent(0)) &&
                     (t2pm.extent(1) == aminc.extent(1)) &&
                     (t2pm.extent(2) == aminc.extent(2)) &&
                     (t2pm.extent(3) == aminc.extent(3));
    if (pcomp && !apc) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "time2_posfix: t2pm and the am_vector spin differ"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    par_for("m1_t2_aminc", DevExeSpace(), 0, static_cast<int>(aminc.extent(0)) - 1, 0,
            static_cast<int>(aminc.extent(1)) - 1, 0, static_cast<int>(aminc.extent(2)) - 1,
            0, static_cast<int>(aminc.extent(3)) - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real v = c1*a1_(m,k,j,i);
      if (two) {v += c2*a2_(m,k,j,i);}
      if (apc && pma_(m,k,j,i) > 0.0) {v += ccomp*a1_(m,k,j,i);}
      ai_(m,k,j,i) = v;
    });
  }
  // consrc-1011 (implicit_gas_src = conservative): the same increment of the stored
  // cell F, from its slopes (the old vector of the write-back's F^old)
  csinc_ok = false;
  if (impl_gas_src_cons && csk1_ok && csk1.extent(0) == t2k1.extent(0)) {
    auto ck1_ = csk1;
    auto ck2_ = csk2;
    auto ci_ = csinc;
    auto pm_ = t2pm;
    par_for("m1_t2_csinc", DevExeSpace(), 0, nmb1, 0, 2, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
      Real v = c1*ck1_(m,n,k,j,i);
      if (two) {v += c2*ck2_(m,n,k,j,i);}
      if (pcomp && pm_(m,k,j,i) > 0.0) {v += ccomp*ck1_(m,n,k,j,i);}
      ci_(m,n,k,j,i) = v;
    });
    csinc_ok = true;
  }
  t2_solve = (stage == 1) ? M1_T2S_STAGE1 : M1_T2S_STAGE2;
  t2_fail = false;
  nsub = 1;
  dt_sub = g*dt;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::Time2Stash
//! \brief pgrst-1021 time2_k1_reset: the radiation E at the hesdirk2 stage-loop start (the
//! state the FSAL slope K1 belongs to; the Driver calls it before the Heun stage 1)

void RadiationM1::Time2Stash() {
  t2_ys_ok = false;
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
//! \fn void RadiationM1::Time2Restore
//! \brief a stage failed: put U^n back (radiation from u1 and the face copies, gas from
//! the hydro's u1), refresh the hydro ghosts and primitives.  The Driver then reruns the
//! Heun stages and takes the step with backward Euler.

void RadiationM1::Time2Restore(Driver *pdrive) {
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
  // posfix-amvec-1022 <hydro>/am_vector: the spin reservoir back to U^n too (am_spin1,
  // taken by the hydro's stage 1 and untouched since); the redo's hydro stage 1 copies
  // am_spin into am_spin1 again
  if (am_vec) {
    hydro::AmModule *pam = hydro::AmModuleOf(pmy_pack);
    Kokkos::deep_copy(DevExeSpace(), pam->am_spin, pam->am_spin1);
  }
  t2_fail = false;
  t2_ok = false;
  t2_nfall += 1.0;
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
//! \fn bool RadiationM1::Time2PosFix
//! \brief t2pos-1020 <rad_m1>/time2_posfix: the stage-1 solve (converged or not) left E,
//! T or the written-back gas eint <= 0 in the cells it marked in t2pm (ImplicitSolve).
//! In the thin, free-streaming floor shell (c dt/dr >> 1) the explicit part (1-g) dt K1
//! of the old vector of H-ESDIRK2 stage 1 is a stiff extrapolation that the L-stable
//! stage solve cannot undo there.
//! The stage-1 start (radiation, faces, gas, closure memory) is
//! put back, the marked set is grown by time2_posfix_grow face-neighbour layers, t2inc is
//! dropped in it (cells; faces with a marked side), and stage 1 is solved again (at most
//! time2_posfix times per step): in the set the stage is backward Euler of g dt from U^n,
//! which is positive wherever the operator is an M-matrix.  Time2FormStage(2) puts the
//! dropped increment back into the stage-2 old vector (weight of Y1 in its start), so the
//! step keeps the weights b of the tableau (consistent, conservative; stage 1 is first
//! order in the set only).  A step whose stage 1 is admissible is untouched.  Returns
//! false when no re-solve is wanted (the caller then takes the usual rejection path).

bool RadiationM1::Time2PosFix(Driver *pdrive) {
  if (!t2_pfix_want) return false;
  t2_pfix_want = false;
  t2_pfix_n += 1;
  Mesh *pmesh = pmy_pack->pmesh;
  // the stage-1 start back: radiation U^n (u1, face copies), the gas after Heun stage 1
  Kokkos::deep_copy(DevExeSpace(), u0, u1);
  Kokkos::deep_copy(DevExeSpace(), f0x1, t2f1);
  if (trans_on) {
    Kokkos::deep_copy(DevExeSpace(), f0x2, t2f2);
    if (trans_x3) {Kokkos::deep_copy(DevExeSpace(), f0x3, t2f3);}
  }
  if (t2_pfix_host) {
    // gpumem-1022 time2_posfix_host: the host copies (same bytes)
    if (ctr_mem.extent(0) > 0 && t2csv_h.extent(0) == ctr_mem.extent(0)) {
      Kokkos::deep_copy(DevExeSpace(), ctr_mem, t2csv_h);
    }
  } else if (ctr_mem.extent(0) > 0 && t2csv.extent(0) == ctr_mem.extent(0)) {
    Kokkos::deep_copy(DevExeSpace(), ctr_mem, t2csv);
  }
  FluidRef fl = FluidRef::Get(pmy_pack);
  if (fl.on) {
    if (t2_pfix_host) {
      Kokkos::deep_copy(DevExeSpace(), fl.u0, t2gsv_h);
    } else {
      Kokkos::deep_copy(DevExeSpace(), fl.u0, t2gsv);
    }
    FluidRefresh(pdrive);
  }
  // posfix-amvec-1022 <hydro>/am_vector: the spin reservoir back to the stage-1 start
  const bool amp = am_vec && (t2ssv.extent_int(0) > 0);
  if (amp) {
    Kokkos::deep_copy(DevExeSpace(), hydro::AmModuleOf(pmy_pack)->am_spin, t2ssv);
  }
  // grow the marked set (face neighbours) and fill its ghosts: iw slot KZ (Krylov
  // scratch, free between solves) carries it through the width-1 exchange; the ghosts of
  // a physical boundary stay 0
  auto &indcs = pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n1 = indcs.nx1 + 2*indcs.ng;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
  const bool d2 = (indcs.nx2 > 1), d3 = (indcs.nx3 > 1);
  auto iw_ = iw;
  auto pm_ = t2pm;
  const int sl = M1_IW_KZ;
  auto to_slot = [&]() {
    par_for("m1_t2p_ts", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const bool in = (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke);
      iw_(m,sl,k,j,i) = in ? pm_(m,k,j,i) : 0.0;
    });
    ImplicitHaloExchange(1, sl);
  };
  for (int r = 0; r < t2_pfix_grow; ++r) {
    to_slot();
    par_for("m1_t2p_grow", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real w = fmax(iw_(m,sl,k,j,i), fmax(iw_(m,sl,k,j,i-1), iw_(m,sl,k,j,i+1)));
      if (d2) {w = fmax(w, fmax(iw_(m,sl,k,j-1,i), iw_(m,sl,k,j+1,i)));}
      if (d3) {w = fmax(w, fmax(iw_(m,sl,k-1,j,i), iw_(m,sl,k+1,j,i)));}
      pm_(m,k,j,i) = w;
    });
  }
  to_slot();
  par_for("m1_t2p_back", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    pm_(m,k,j,i) = iw_(m,sl,k,j,i);
  });
  Real nc = 0.0;
  par_reduce_flat("m1_t2p_cnt", 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &a) {
    if (pm_(m,k,j,i) > 0.0) a += 1.0;
  }, Kokkos::Sum<Real>(nc));
#if MPI_PARALLEL_ENABLED
  {Real g;
  MPI_Allreduce(&nc, &g, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  nc = g;}
#endif
  // drop the stage-1 increment in the set (cells; faces with a marked side)
  auto in_ = t2inc;
  {
    const int q3 = static_cast<int>(t2inc.extent(2));
    const int q2 = static_cast<int>(t2inc.extent(3));
    const int q1 = static_cast<int>(t2inc.extent(4));
    par_for("m1_t2p_drop", DevExeSpace(), 0, nmb1, 0, M1_T2_NK-1, 0, q3-1, 0, q2-1,
            0, q1-1,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
      Real w = pm_(m,k,j,i);
      if (n == M1_T2_F1 && i > 0) {w = fmax(w, pm_(m,k,j,i-1));}
      if (n == M1_T2_F1+1) {w = d2 ? ((j > 0) ? fmax(w, pm_(m,k,j-1,i)) : w) : 0.0;}
      if (n == M1_T2_F1+2) {w = d3 ? ((k > 0) ? fmax(w, pm_(m,k-1,j,i)) : w) : 0.0;}
      if (w > 0.0) {in_(m,n,k,j,i) = 0.0;}
    });
    if (csinc_ok) {
      auto ci_ = csinc;
      par_for("m1_t2p_dropc", DevExeSpace(), 0, nmb1, 0, 2, 0, q3-1, 0, q2-1, 0, q1-1,
      KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
        if (pm_(m,k,j,i) > 0.0) {ci_(m,n,k,j,i) = 0.0;}
      });
    }
    // posfix-amvec-1022: and the spin reservoir's increment (cells)
    if (am_vec && aminc.extent_int(0) > 0) {
      auto ai_ = aminc;
      par_for("m1_t2p_dropa", DevExeSpace(), 0, static_cast<int>(aminc.extent(0)) - 1,
              0, static_cast<int>(aminc.extent(1)) - 1,
              0, static_cast<int>(aminc.extent(2)) - 1,
              0, static_cast<int>(aminc.extent(3)) - 1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        if (pm_(m,k,j,i) > 0.0) {ai_(m,k,j,i) = 0.0;}
      });
    }
  }
  t2_pfix_any = true;
  t2_pfix_nsolve += 1.0;
  t2_pfix_ncell += nc;
  if (t2_pfix_n == 1) {t2_pfix_nstep += 1.0;}
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1> time2_posfix: stage 1 re-solved at cycle " << pmesh->ncycle
              << " (round " << t2_pfix_n << ", " << nc << " cells without the stage-1 "
              << "increment)" << std::endl;
  }
  // the re-solve: a cold Picard start (the failed solve stored its increment), the
  // stage-1 closure build of the first attempt is kept (same U^n; history-free)
  pred_ok = false;
  t2_solve = M1_T2S_STAGE1;
  t2_fail = false;
  nsub = 1;
  dt_sub = kT2G*pmesh->dt;
  t2_pfix_redo = true;
  pdrive->ExecuteTaskList(pmesh, "m1_before_stagen", 1);
  pdrive->ExecuteTaskList(pmesh, "m1_stagen", 1);
  pdrive->ExecuteTaskList(pmesh, "m1_after_stagen", 1);
  t2_pfix_redo = false;
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationM1::T2EnthBe
//! \brief open-ot-1013 time2_enth_be: this (non-stage, backward-Euler) solve of a hesdirk2
//! run builds the enthalpy face a_f with the stage rule (time2_enth_vel = central)

bool RadiationM1::T2EnthBe() const {
  return t2_enth_be && (time_scheme == M1_TIME_HESDIRK2) && (t2_afmode == 2) &&
         (t2_solve != M1_T2S_STAGE1) && (t2_solve != M1_T2S_STAGE2);
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationM1::Time2StartProbe
//! \brief time2_start = slope (h2order-1011): on a hesdirk2 step without a stored slope,
//! a backward-Euler probe solve of dt_p = time2_start_frac dt from U^n (radiation, faces
//! and gas; no hydro) stores K1 = (Y_p - U^n)/dt_p = f_I(U^n) + O(dt_p); U^n is then put
//! back and the step runs the two stages, so no step of the run is Lie split.  Returns
//! true when a probe was taken.  The predictor history the probe leaves is dropped (it
//! only moves Picard starting points).

bool RadiationM1::Time2StartProbe(Driver *pdrive, Real dt) {
  // vgdrot2-1012: VetGdRotStep asks for the probe at a rotation (any time2_start)
  if (time_scheme != M1_TIME_HESDIRK2 || (t2_start != 1 && !t2_rotprobe)) return false;
  t2_rotprobe = false;
  if (Time2Active()) return false;
  FluidRef fl = FluidRef::Get(pmy_pack);
  if (fl.on && fl.mhd) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<rad_m1>/time2_start = slope is not wired for MHD"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // U^n into the stage-start copies (Time2FormStage(1) overwrites them again)
  Kokkos::deep_copy(DevExeSpace(), u1, u0);
  Kokkos::deep_copy(DevExeSpace(), t2f1, f0x1);
  if (trans_on) {
    Kokkos::deep_copy(DevExeSpace(), t2f2, f0x2);
    if (trans_x3) {Kokkos::deep_copy(DevExeSpace(), t2f3, f0x3);}
  }
  if (fl.on) {Kokkos::deep_copy(DevExeSpace(), fl.u1, fl.u0);}
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
//! \brief rotslope-1012, vet_gd_rotate_slope = correct (asked for by VetGdRotStep on the
//! first step of a new vet_gd rotation block, with a valid FSAL slope K1 = f_old(U^n) of
//! the OLD set): one stage-type solve of dt_c = vet_gd_rotate_slope_frac dt from U^n
//! (radiation, faces, gas; no hydro) with the NEW set and the old vector U^n - dt_c K1
//! (cell F: - dt_c csk1).  Its slope (Y - old)/dt_c = K1 + z/dt_c with, to first order in
//! z = Y - U^n,  z = dt_c (I - dt_c A_new)^{-1} (f_new - f_old)(U^n): the stored slope is
//! K1 + (I - dt_c A)^{-1} Delta f -- the closure jump of the slope where dt_c A is small
//! (the thin regime, where keep leaves an O(dt Delta f) error per jump = first order), and
//! keep's K1 where the cell is stiff (there the L-stable stages forget a stage-1 slope
//! error and keep is second order; replacing all of K1 by a BE probe slope is what broke
//! the stiff ladders).  Only the jump steps change, so an O(dt) error of the correction
//! itself stays O(dt^2) globally.  The solve runs as a STAGE1 solve (its slope lands in
//! t2k2 / csk2 and is copied to t2k1 / csk1); the predictor history is kept
//! (implicit_predictor is off for this solve); U^n is put back.  With
//! vet_gd_rotate_slope_reuse the stage-1 solve of the step keeps this solve's closure
//! build (the same set, built by Time2VetColAt(1) at U^n + theta dt K1_old with the gas
//! at U^n instead of after the Heun stage: an O(dt) state shift on the jump steps only).
//! A failed correction solve keeps K1 (keep).  Returns true when a solve was taken.

bool RadiationM1::Time2RotCorrect(Driver *pdrive, Real dt) {
  if (!t2_rotcorr) return false;
  t2_rotcorr = false;
  if (time_scheme != M1_TIME_HESDIRK2 || !Time2Active()) return false;
  FluidRef fl = FluidRef::Get(pmy_pack);
  if (fl.on && fl.mhd) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<rad_m1>/vet_gd_rotate_slope = correct is not wired for "
              << "MHD" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const Real dtc = vgd_rot_cfrac*dt;
  // U^n into the stage-start copies (Time2FormStage(1) overwrites them again)
  Kokkos::deep_copy(DevExeSpace(), u1, u0);
  Kokkos::deep_copy(DevExeSpace(), t2f1, f0x1);
  if (trans_on) {
    Kokkos::deep_copy(DevExeSpace(), t2f2, f0x2);
    if (trans_x3) {Kokkos::deep_copy(DevExeSpace(), t2f3, f0x3);}
  }
  if (fl.on) {Kokkos::deep_copy(DevExeSpace(), fl.u1, fl.u0);}
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
  }
  t2_solve = M1_T2S_NONE;
  t2_fail = false;
  // the stage-1 solve keeps this build (vet_col predict | rebuild only)
  t2_vcreuse = ok && vgd_rot_creuse && vet_col && !vet_sc &&
               (t2_vcmode == 1 || t2_vcmode == 2);
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::Time2VetStart
//! \brief closure = vet_sc, stage 1: the formal solution at U^n -- E^n, T^n and the
//! opacities of U^n (the hydro's u1) -- and then the extrapolated tensor
//! D* = D^n + (dt/dt_prev)(D^n - D^{n-1}), which both stage solves read.  Called from
//! ImplicitSolve at the point where the backward-Euler step calls VetShortChar; the
//! solve's own EN (old vector), TP (start T) and opacities are put back afterwards.
//! closure = vet_col (m1-sph2, tests_m1/runs_5h_sph2): the same, with VetColBuild (the
//! tensor tau_ten and the surface q at U^n) and NO extrapolation: both stages use D^n
//! (an extrapolated chi and q were measured: same order, 3e-11 from D^n).

void RadiationM1::Time2VetStart() {
  // vet_sc_every (rad_m1_mr.cpp): no formal solution on this build, the tensor is
  // extrapolated from the last two
  if (vet_sc && vsc_every > 1 && VscSkip()) {
    Time2VetExtrapolate();
    return;
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto vn_ = vet_now;
  FluidRef fl = FluidRef::Get(pmy_pack);
  const bool hh = fl.on;
  auto uh1 = hh ? fl.u1 : u1;
  const bool etg = hh ? fl.use_etotgrav : false;
  auto phicc = hh ? fl.phicc0 : arad_ref;
  // m1-mhd: U^n = (u1, b1); its |B|^2/2 goes to emag1
  const bool mhd = fl_mhd;
  if (mhd) EmagBuild(true);
  auto emag_ = emag1;
  const Real efl = e_floor;
  auto u0_ = u0;
  // cs-hydro-energy (10-05): the metric kinetic energy on the cubed sphere (covariant
  // momentum), as RadiationM1::Opacity and m1_impl_i1
  const bool csk = pmy_pack->pmesh->use_cubed_sphere;
  auto ccl = pmy_pack->pcoord->cos_cell;
  // the fast form (default): ONE fused pass over U^n (the hydro u1) gives T^n and the
  // opacities of U^n, with exactly the arithmetic of RadiationM1::Opacity, and the
  // stage-start opacities are saved and copied back afterwards instead of re-evaluated.
  // The debug opacity options (opac_freeze, dbg_opac_patch) take the general form.
  const bool fast = hh && !opac_zero && !opac_freeze && (dbg_opac_patch == 1.0) &&
                    !(opacity_type == M1_OPAC_TABLE && otab.nT <= 0);
  if (fast) {
    Kokkos::deep_copy(DevExeSpace(), vet_opac, opac);
    auto eos = fl.eos;
    auto opac_ = opac;
    const int n1 = static_cast<int>(opac.extent(4));
    const int n2 = static_cast<int>(opac.extent(3));
    const int n3 = static_cast<int>(opac.extent(2));
    const int otype = opacity_type;
    const Real kp = kappa_p, kev = kappa_e, kf = kappa_f, kss = kappa_s;
    const Real rref = opac_rho_ref, tref = opac_t_ref, aa = opac_a, bb = opac_b;
    M1OpacTab ot = otab;
    const bool osc = vscat;   // vet_scatter: rho kappa_e -> opac(M1_OP_S)
    par_for_lb("m1_t2_vsf", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real d = uh1(m,IDN,k,j,i);
      Real ke_dens = 0.5*(SQR(uh1(m,IM1,k,j,i)) + SQR(uh1(m,IM2,k,j,i)) +
                          SQR(uh1(m,IM3,k,j,i)))/fmax(d, 1.0e-300);
      if (csk) {
        const Real c = ccl(m,k,j);
        const Real m2 = uh1(m,IM2,k,j,i), m3 = uh1(m,IM3,k,j,i);
        ke_dens = 0.5*(SQR(uh1(m,IM1,k,j,i)) + (m2*m2 + m3*m3 - 2.0*c*m2*m3)/(1.0 - c*c))
                  /fmax(d, 1.0e-300);
      }
      Real eint = uh1(m,IEN,k,j,i) - ke_dens;
      if (etg) eint -= d*phicc(m,k,j,i);
      if (mhd) eint -= emag_(m,k,j,i);
      Real t = eos.Temperature(d, fmax(eint, 0.0));
      Real op, oe, of, os;
      if (otype == M1_OPAC_TABLE) {
        M1TableOpacities(ot, d, t, op, oe, of, os);
      } else {
        M1Opacities(otype, d, t, kp, kev, kf, kss, rref, tref, aa, bb, op, oe, of, os);
      }
      if (osc && otype == M1_OPAC_TABLE) {M1ScatterEos(ot, eos, d, t, of, os);}
      opac_(m,M1_OP_P,k,j,i) = d*op;
      opac_(m,M1_OP_E,k,j,i) = d*oe;
      opac_(m,M1_OP_T,k,j,i) = d*(of + os);
      if (osc) {opac_(m,M1_OP_S,k,j,i) = d*os;}
      if (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke) {
        vn_(m,0,k,j,i) = iw_(m,M1_IW_EN,k,j,i);
        vn_(m,1,k,j,i) = iw_(m,M1_IW_TP,k,j,i);
        iw_(m,M1_IW_EN,k,j,i) = fmax(u0_(m,M1_E,k,j,i), efl);
        iw_(m,M1_IW_TP,k,j,i) = t;
      }
    });
  } else if (hh) {
    auto eos = fl.eos;
    par_for("m1_t2_vs0", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      vn_(m,0,k,j,i) = iw_(m,M1_IW_EN,k,j,i);
      vn_(m,1,k,j,i) = iw_(m,M1_IW_TP,k,j,i);
      iw_(m,M1_IW_EN,k,j,i) = fmax(u0_(m,M1_E,k,j,i), efl);
      Real dd = uh1(m,IDN,k,j,i);
      Real idd = 1.0/fmax(dd, 1.0e-300);
      Real ekin = 0.5*(SQR(uh1(m,IM1,k,j,i)) + SQR(uh1(m,IM2,k,j,i)) +
                       SQR(uh1(m,IM3,k,j,i)))*idd;
      if (csk) {
        const Real c = ccl(m,k,j);
        const Real m2 = uh1(m,IM2,k,j,i), m3 = uh1(m,IM3,k,j,i);
        ekin = 0.5*(SQR(uh1(m,IM1,k,j,i)) + (m2*m2 + m3*m3 - 2.0*c*m2*m3)/(1.0 - c*c))
               *idd;
      }
      Real egrv = etg ? (dd*phicc(m,k,j,i)) : 0.0;
      if (mhd) egrv += emag_(m,k,j,i);
      Real eg = uh1(m,IEN,k,j,i) - ekin - egrv;
      iw_(m,M1_IW_TP,k,j,i) = eos.Temperature(dd, fmax(eg, 1.0e-300));
    });
    // the opacities of U^n: Opacity reads the fluid's u0 and emag0, so both are swapped
    // with their U^n twins for the call (emag_hold: Opacity must not rebuild emag0)
    fl.SwapU01();
    std::swap(emag0, emag1);
    emag_hold = true;
    (void) Opacity(nullptr, 1);
    emag_hold = false;
    std::swap(emag0, emag1);
    fl.SwapU01();
  } else {
    par_for("m1_t2_vs0r", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      vn_(m,0,k,j,i) = iw_(m,M1_IW_EN,k,j,i);
      iw_(m,M1_IW_EN,k,j,i) = fmax(u0_(m,M1_E,k,j,i), efl);
    });
  }
  // closure = vet_col (m1-sph2): the per-column formal solution of the same state;
  // its tensor (tau_ten) and surface q are then kept for both stages (no extrapolation)
  if (vet_col) {
    VetColBuild();
  } else {
    VetShortChar();
    if (vsc_every > 1) {VscStore();}
  }
  if (fast) {
    Kokkos::deep_copy(DevExeSpace(), opac, vet_opac);
  } else if (hh) {
    (void) Opacity(nullptr, 1);
  }
  par_for("m1_t2_vs1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    iw_(m,M1_IW_EN,k,j,i) = vn_(m,0,k,j,i);
    if (hh) {iw_(m,M1_IW_TP,k,j,i) = vn_(m,1,k,j,i);}
  });
  if (vet_sc) {Time2VetExtrapolate();}
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::Time2VetExtrapolate
//! \brief vet_cell <- D* = D^n + r (D^n - D^{n-1}), r = dt/dt_prev, for the uniaxial
//! (chi, n) and the full D_ab; a cell whose D* is not realizable keeps D^n (counted).
//! vet_prev <- D^n for the next step.

void RadiationM1::Time2VetExtrapolate() {
  auto vc_ = vet_cell;
  auto vp_ = vet_prev;
  // implicit_mr_every: the step is the multi-rate window Delta
  const Real dt = mr_on ? mr_dt : pmy_pack->pmesh->dt;
  const Real r = (t2_vprev && t2_dtprev > 0.0 && t2_vext) ? (dt/t2_dtprev) : 0.0;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const bool full = vet_full;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  int n3 = static_cast<int>(vet_cell.extent(2));
  int n2 = static_cast<int>(vet_cell.extent(3));
  int n1 = static_cast<int>(vet_cell.extent(4));
  Real nclip = 0.0;
  if (r == 0.0 && impl_fastk) {
    // implicit_fast_kernels: nothing is extrapolated, only D^n is saved (bitwise the
    // same values as the reduction below writes, which adds no clip when r = 0)
    const int nc = M1_T2_NVET;
    const int nt = (nmb1 + 1)*n3*n2*n1;
    par_for("m1_t2_vcopy", DevExeSpace(), 0, nc-1, 0, nt-1,
    KOKKOS_LAMBDA(const int c, const int idx) {
      int m = idx/(n3*n2*n1);
      int r3 = idx - m*(n3*n2*n1);
      int k = r3/(n2*n1);
      r3 -= k*(n2*n1);
      int j = r3/n1;
      int i = r3 - j*n1;
      vp_(m,c,k,j,i) = vc_(m,M1_VET_CHI+c,k,j,i);
    });
    t2_vprev = true;
    return;
  }
  Kokkos::parallel_reduce("m1_t2_vext",
  Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,0,0,0},
                                         {nmb1+1,n3,n2,n1}),
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lc) {
    Real dn[M1_T2_NVET], ds[M1_T2_NVET];
    for (int c = 0; c < M1_T2_NVET; ++c) {
      dn[c] = vc_(m,M1_VET_CHI+c,k,j,i);
      ds[c] = dn[c] + r*(dn[c] - vp_(m,c,k,j,i));
      vp_(m,c,k,j,i) = dn[c];
    }
    if (r == 0.0) return;
    bool ok = (ds[0] >= 0.0) && (ds[0] <= 1.0);
    Real nn = sqrt(ds[1]*ds[1] + ds[2]*ds[2] + ds[3]*ds[3]);
    ok = ok && (nn > 0.5);
    if (ok) {
      ds[1] /= nn;
      ds[2] /= nn;
      ds[3] /= nn;
    }
    if (full) {
      // 11 22 33 12 13 23 at 4..9: symmetric positive semi-definite, diagonal <= 1
      Real a = ds[4], b = ds[5], c = ds[6], d = ds[7], e = ds[8], f = ds[9];
      Real m2 = a*b - d*d;
      Real det = a*(b*c - f*f) - d*(d*c - f*e) + e*(d*f - b*e);
      ok = ok && (a >= 0.0) && (b >= 0.0) && (c >= 0.0) && (a <= 1.0) && (b <= 1.0) &&
           (c <= 1.0) && (m2 >= 0.0) && (det >= 0.0);
    }
    if (ok) {
      for (int q = 0; q < M1_T2_NVET; ++q) {vc_(m,M1_VET_CHI+q,k,j,i) = ds[q];}
    } else if (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke) {
      lc += 1.0;   // counted over the ACTIVE cells only
    }
  }, Kokkos::Sum<Real>(nclip));
  t2_nclip += nclip;
  t2_vprev = true;
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
  auto lim = std::get<34>(ctx_);   // time2_pred_cap: theta per cell (predfix-1011)
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
//! \brief time2_vet_source = relaxed | state | blend (vcpred-1010, sporder-1012; was
//! time2_vet_col_source, its alias on vet_col): the VET source of the predict / rebuild
//! build, for every formal solution (unify-sc-gd-1012 #3).  vet_source = relaxed relaxes
//! the build state's gas and E by a local backward-Euler exchange over dt, i.e. it moves
//! the source one step AHEAD of the state it is given.  At U^n (lag, backward Euler)
//! that is the step's end; at the predicted P = U^n + dt K1, already at t^{n+1} +
//! O(dt^2), it lands near t^{n+2}: an O(dt) source error, FIRST order in time wherever
//! the exchange is resolved and out of equilibrium (vcpred_1010 neq1_u_T: E 1.71 1.57
//! 1.40, T 1.67 1.47 1.23).  state = the source of P itself (eps B + (1 - eps) E, or
//! VetScatterMix with J = E_P under vet_scatter): 1.92 1.87 1.91 / 1.96 1.98 1.99.
//! blend (sporder-1012) = S*(0) + w (S*(dt) - S*(0)), w = x^2/(x0^2 + x^2) with x the
//! cell's exchange stiffness c dt rho kappa (radiation + gas temperature): the state's
//! own source where the exchange is resolved (x << x0, O(dt^3) from it), the relaxed
//! source where it is stiff (VetRelaxedSourceB); time2_vet_source_x0 (alias
//! time2_vet_col_source_x0; default 1).  lag / extrap keep vet_source.
//! DEFAULT on vet_col / vet_gd: blend on a fresh run with predict | rebuild (where it
//! acts), relaxed on a restart whose file lacks the key; vet_sc: state (its old source;
//! OPEN, see VetInit); the resolved value is echoed; explicit input overrides.

void RadiationM1::Time2VetSourceInit(ParameterInput *pin, const bool act, const bool col) {
  const std::string key = col ? M1KeyName(pin, "time2_vet_source", {"time2_vet_col_source"})
                              : std::string("time2_vet_source");
  if (!act && !pin->DoesParameterExist("rad_m1", key)) {return;}
  // unify-sc-gd-1012 #3 (OPEN, rad_m1_vet.cpp VetInit): vet_sc keeps state (fresh and
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
//! \fn void RadiationM1::Time2VetColAt
//! \brief closure = vet_col, time2_vet_col = predict | rebuild (m1-sp-order2b), and
//! closure = vet_sc, time2_vet_sc = predict | rebuild (torder-1010, no inner face flux):
//! the formal solution of a state at t^{n+1} accurate to O(dt^2), in place of U^n.
//!   which = 1 (stage 1): P = the stage-1 start + dt K1: radiation E^n + dt K1_E, the
//!             inner face flux F^n + dt K1_F, gas = the Heun predictor + dt K1 of the
//!             momentum and total energy (the implicit coupling's rate; the hydro u0 is
//!             already the old vector, Heun + (1-g) dt K1, so g dt K1 is added);
//!   which = 2 (stage 2, rebuild): the stage-1 solution Y1: E and T saved at the end of
//!             the stage-1 solve (Time2VetColSaveY1), rho of the stage-2 start, the inner
//!             face flux 2 F_avg - F^n.
//! The source (M1_IW_EN, M1_IW_TP), the extinction (M1_IW_KT), the opacities and f0x1 at
//! the inner face are swapped in for the build and put back afterwards; the solve itself
//! sees its own stage-start state as before.  Opacity options without the fused form
//! (opac_freeze, dbg_opac_patch, an empty table) fall back to D(U^n) (counted).

void RadiationM1::Time2VetColAt(int which) {
  FluidRef fl = FluidRef::Get(pmy_pack);
  const bool hh = fl.on;
  const bool fast = hh && !opac_zero && !opac_freeze && (dbg_opac_patch == 1.0) &&
                    !(opacity_type == M1_OPAC_TABLE && otab.nT <= 0);
  if (hh && !fast) {
    t2_vcnfb += 1.0;
    if (which == 1) {
      Time2VetStart();
    }
    return;
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  // implicit_mr_every (rad_m1_mr.cpp): the step is the window Delta, and the stage-A old
  // vector is Y_0 itself (t2inc = 0), so the gas takes the whole Delta K1
  const Real dt = mr_on ? mr_dt : pmy_pack->pmesh->dt;
  const Real gdt = mr_on ? dt : kT2G*dt;
  const bool two = (which == 2);
  auto iw_ = iw;
  auto vn_ = vet_now;
  auto u0_ = u0;
  auto k1_ = t2k1;
  auto f0_ = f0x1;
  auto t2f1_ = t2f1;
  const Real efl = e_floor;
  // time2_pred_cap (predfix-1011, Time2Init): theta per cell, P = U^n + theta dt K1
  const bool lim = !two && (t2_pcap > 0.0);
  const bool dbg = !two && (t2_dbg_pred_n < t2_dbg_pred);
  if (lim || dbg) {Time2PredTheta(dt, gdt, dbg);}
  auto th_ = t2pth;
  // the inner face flux (vet_col_order2 core rays): saved in slot 5 at i = is.
  // closure = vet_sc (time2_vet_sc, torder-1010) does not read it
  const bool vcol = vet_col;
  if (vcol) {
    par_for("m1_t2_vcf", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      const Real f = f0_(m,k,j,is);
      const Real tdt = lim ? th_(m,k,j,is)*dt : dt;
      vn_(m,5,k,j,is) = f;
      f0_(m,k,j,is) = two ? (2.0*f - t2f1_(m,k,j,is))
                          : (f + tdt*k1_(m,M1_T2_F1,k,j,is));
    });
  }
  if (hh) {
    Kokkos::deep_copy(DevExeSpace(), vet_opac, opac);
    auto uh = fl.u0;
    const bool etg = fl.use_etotgrav;
    auto phicc = fl.phicc0;
    auto eos = fl.eos;
    auto opac_ = opac;
    const bool mhd = fl_mhd;      // m1-mhd: emag0 = the current field's |B|^2/2
    auto emag_ = emag0;
    const int otype = opacity_type;
    const Real kp = kappa_p, kev = kappa_e, kf = kappa_f, kss = kappa_s;
    const Real rref = opac_rho_ref, tref = opac_t_ref, aa = opac_a, bb = opac_b;
    M1OpacTab ot = otab;
    // m1-fast5-sp: an ideal-gas EOS with analytic opacities gets its own kernel (T and
    // the opacities by the same expressions, bitwise): the generic one carries the EOS
    // and opacity table code, 392 B of scratch per lane, even when both are off
    const Real gam = eos.gamma;
    const bool vcid = !eos.tbl.active && (otype != M1_OPAC_TABLE);
    // nvcc forbids generic (auto) extended lambdas, so the tag-dependent
    // helper and its kernel are the function template M1T2VcpLaunch
    // (above); vcb_ctx is what the helper captured, by value.
    auto vcb_ctx = std::make_tuple(aa, bb, dt, efl, eos, etg, gam, gdt, iw_, k1_, kev, kf,
                                   kp, kss, opac_, ot, otype, phicc, rref, tref, two, u0_,
                                   uh, vn_, nmb1, ks, ke, js, je, is, ie, mhd, emag_,
                                   vscat, lim, th_);
    if (vcid) {
      M1T2VcpLaunch(vcb_ctx, std::true_type{});
    } else {
      M1T2VcpLaunch(vcb_ctx, std::false_type{});
    }
  } else {
    par_for("m1_t2_vcr", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      vn_(m,0,k,j,i) = iw_(m,M1_IW_EN,k,j,i);
      const Real tdt = lim ? th_(m,k,j,i)*dt : dt;
      const Real e = two ? vn_(m,4,k,j,i)
                         : (u0_(m,M1_E,k,j,i) + tdt*k1_(m,M1_T2_E,k,j,i));
      iw_(m,M1_IW_EN,k,j,i) = fmax(e, efl);
    });
  }
  if (vcol) {
    // time2_vet_col_source = state (vcpred-1010): the source of P itself, not relaxed
    // by a further dt (vet_col, vet_col_lat and vet_gd all read vcol_srelax)
    const bool srs = vcol_srelax;
    if (t2_vcsrc) {vcol_srelax = false;}
    vcol_sbl = t2_vcsbl;
    VetColBuild();
    vcol_sbl = false;
    vcol_srelax = srs;
  } else {
    // unify-sc-gd-1012 (#3): time2_vet_source for vet_sc too, the same flag swap
    const bool srs = vcol_srelax;
    if (t2_vcsrc) {vcol_srelax = false;}
    vcol_sbl = t2_vcsbl;
    VetShortChar();
    vcol_sbl = false;
    vcol_srelax = srs;
    if (dbg) {Time2PredDebugSc();}
  }
  if (hh) {
    Kokkos::deep_copy(DevExeSpace(), opac, vet_opac);
  }
  par_for("m1_t2_vcb", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    iw_(m,M1_IW_EN,k,j,i) = vn_(m,0,k,j,i);
    if (hh) {
      iw_(m,M1_IW_TP,k,j,i) = vn_(m,1,k,j,i);
      iw_(m,M1_IW_KT,k,j,i) = vn_(m,2,k,j,i);
    }
    if (vcol && i == is) {
      f0_(m,k,j,is) = vn_(m,5,k,j,is);
    }
  });
  if (dbg) {t2_dbg_pred_n += 1;}
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::Time2PredTheta
//! \brief time2_pred_cap (predfix-1011): theta = min(1, cap/r) of every active cell into
//! t2pth, r = max(|dt K1_E|/E^n, |eint(P) - eint_b|/eint_b), eint_b the gas internal
//! energy of the stage-1 start without the coupling's K1 (u0 of the hydro holds Heun +
//! (dt - gdt) K1).  Counts the limited cells; dbg prints the predicted state.

void RadiationM1::Time2PredTheta(const Real dt, const Real gdt, const bool dbg) {
  FluidRef fl = FluidRef::Get(pmy_pack);
  const bool hh = fl.on;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto th_ = t2pth;
  auto u0_ = u0;
  auto k1_ = t2k1;
  auto uh = hh ? fl.u0 : u0;
  const bool etg = hh ? fl.use_etotgrav : false;
  auto phicc = hh ? fl.phicc0 : arad_ref;
  const bool mhd = hh && fl_mhd;
  auto emag_ = emag0;
  const Real cap = t2_pcap, efl = e_floor;
  const Real bdt = gdt - dt;   // the base: the stage-1 start without K1
  const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
  Real nlim = 0.0, emin = 1.0e300, nneg = 0.0, ngneg = 0.0, rmax = 0.0, thmin = 1.0;
  Kokkos::parallel_reduce("m1_t2_pth",
  Kokkos::RangePolicy<DevExeSpace>(DevExeSpace(), 0, (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int idx, Real &ln, Real &le, Real &lneg, Real &lg, Real &lr,
                Real &lt) {
    const int m = idx/nkji;
    int q = idx - m*nkji;
    const int k = q/nji + ks;
    q -= (k - ks)*nji;
    const int j = q/ni + js;
    const int i = q - (j - js)*ni + is;
    const Real e0 = u0_(m,M1_E,k,j,i);
    const Real ep = e0 + dt*k1_(m,M1_T2_E,k,j,i);
    Real r = fabs(ep - e0)/fmax(e0, efl);
    if (hh) {
      const Real d = fmax(uh(m,IDN,k,j,i), 1.0e-300);
      Real eb = 0.0, e1 = 0.0;
      for (int w = 0; w < 2; ++w) {
        const Real c = (w == 0) ? bdt : gdt;
        const Real m1 = uh(m,IM1,k,j,i) + c*k1_(m,M1_T2_M1,k,j,i);
        const Real m2 = uh(m,IM2,k,j,i) + c*k1_(m,M1_T2_M1+1,k,j,i);
        const Real m3 = uh(m,IM3,k,j,i) + c*k1_(m,M1_T2_M1+2,k,j,i);
        Real ei = uh(m,IEN,k,j,i) + c*k1_(m,M1_T2_EN,k,j,i)
                  - 0.5*(m1*m1 + m2*m2 + m3*m3)/d;
        if (etg) ei -= uh(m,IDN,k,j,i)*phicc(m,k,j,i);
        if (mhd) ei -= emag_(m,k,j,i);
        if (w == 0) {eb = ei;} else {e1 = ei;}
      }
      if (eb > 0.0) {r = fmax(r, fabs(e1 - eb)/eb);}
      if (!(e1 > 0.0)) {lg += 1.0;}
    }
    if (!(r < 1.0e300)) {r = 1.0e300;}   // NaN / inf: fully lagged
    const Real th = (cap > 0.0 && r > cap) ? cap/r : 1.0;
    th_(m,k,j,i) = th;
    if (th < 1.0) {ln += 1.0;}
    le = fmin(le, ep);
    if (!(ep > 0.0)) {lneg += 1.0;}
    lr = fmax(lr, r);
    lt = fmin(lt, th);
  }, Kokkos::Sum<Real>(nlim), Kokkos::Min<Real>(emin), Kokkos::Sum<Real>(nneg),
     Kokkos::Sum<Real>(ngneg), Kokkos::Max<Real>(rmax), Kokkos::Min<Real>(thmin));
  Real ncell = static_cast<Real>((nmb1 + 1)*nkji);
#if MPI_PARALLEL_ENABLED
  {
    Real a[4] = {nlim, nneg, ngneg, ncell}, b[4];
    MPI_Allreduce(a, b, 4, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    nlim = b[0]; nneg = b[1]; ngneg = b[2]; ncell = b[3];
    Real c[2] = {rmax, -fmin(emin, 1.0e300)}, e[2];
    MPI_Allreduce(c, e, 2, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
    rmax = e[0]; emin = -e[1];
    Real g;
    MPI_Allreduce(&thmin, &g, 1, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
    thmin = g;
  }
#endif
  if (cap > 0.0) {
    t2_pnlim += nlim;
    t2_pncell += ncell;
  }
  if (dbg && global_variable::my_rank == 0) {
    std::cout << "<rad_m1> DEBUG dbg_t2_pred cycle " << pmy_pack->pmesh->ncycle
              << ": P = U^n + dt K1: min E_P=" << emin << " cells E_P<=0: " << nneg
              << " eint_P<=0: " << ngneg << " max r=" << rmax << " | cap=" << cap
              << " limited cells " << nlim << "/" << ncell << " min theta=" << thmin
              << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::Time2PredDebugSc
//! \brief DEBUG dbg_t2_pred: the formal solution of the predict build (vet_sc): mean and
//! extremes of H_1/J and H_2/J over the active cells, mean J, and the sweep's clip count

void RadiationM1::Time2PredDebugSc() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto vc_ = vet_cell;
  const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
  Real sj = 0.0, sh1 = 0.0, sh2 = 0.0, mx1 = -1.0e300, mn1 = 1.0e300, mxa = 0.0;
  Kokkos::parallel_reduce("m1_t2_pdbg",
  Kokkos::RangePolicy<DevExeSpace>(DevExeSpace(), 0, (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int idx, Real &lj, Real &l1, Real &l2, Real &lx, Real &ln,
                Real &la) {
    const int m = idx/nkji;
    int q = idx - m*nkji;
    const int k = q/nji + ks;
    q -= (k - ks)*nji;
    const int j = q/ni + js;
    const int i = q - (j - js)*ni + is;
    const Real jj = vc_(m,M1_VET_J,k,j,i);
    const Real h1 = (jj > 0.0) ? vc_(m,M1_VET_H1,k,j,i)/jj : 0.0;
    const Real h2 = (jj > 0.0) ? vc_(m,M1_VET_H1+1,k,j,i)/jj : 0.0;
    lj += jj;
    l1 += h1;
    l2 += h2;
    lx = fmax(lx, h1);
    ln = fmin(ln, h1);
    la = fmax(la, sqrt(h1*h1 + h2*h2));
  }, Kokkos::Sum<Real>(sj), Kokkos::Sum<Real>(sh1), Kokkos::Sum<Real>(sh2),
     Kokkos::Max<Real>(mx1), Kokkos::Min<Real>(mn1), Kokkos::Max<Real>(mxa));
  Real ncl = 0.0;
  if (vet_dbgclip.extent_int(0) > 0) {
    auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), vet_dbgclip);
    ncl = h(0);
    Kokkos::deep_copy(vet_dbgclip, 0.0);
  }
  const Real nc = static_cast<Real>((nmb1 + 1)*nkji);
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1> DEBUG dbg_t2_pred cycle " << pmy_pack->pmesh->ncycle
              << ": SC of P (rank 0): mean J=" << sj/nc << " mean H1/J=" << sh1/nc
              << " [" << mn1 << ", " << mx1 << "] mean H2/J=" << sh2/nc
              << " max |H|/J=" << mxa << " clipped I<0: " << ncl << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::Time2VetColSaveY1
//! \brief time2_vet_col = rebuild: at the end of the stage-1 solve, E and T of Y1

void RadiationM1::Time2VetColSaveY1() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto vn_ = vet_now;
  auto u0_ = u0;
  par_for("m1_t2_vcy", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    vn_(m,3,k,j,i) = iw_(m,M1_IW_TP,k,j,i);
    vn_(m,4,k,j,i) = u0_(m,M1_E,k,j,i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::Time2VetColExtrap
//! \brief time2_vet_col = extrap (DIAGNOSTIC): f_K and q of D(U^n) extrapolated to
//! t^{n+1} with the previous step's, f_K clipped to [vet_col_fk_min, 1], q to its range

void RadiationM1::Time2VetColExtrap() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const Real dtn = mr_on ? mr_dt : pmy_pack->pmesh->dt;   // implicit_mr_every: Delta
  const Real r = (t2_vcprev && t2_dtprev > 0.0) ? (dtn/t2_dtprev) : 0.0;
  auto tt_ = tau_ten;
  auto vp_ = vcol_prev;
  const Real fkm = vcol_axis_flux ? (1.0/3.0) : vcol_fkmin;
  par_for("m1_t2_vcx", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real dn = tt_(m,0,k,j,i);
    tt_(m,0,k,j,i) = fmin(fmax(dn + r*(dn - vp_(m,k,j,i)), fkm), 1.0);
    vp_(m,k,j,i) = dn;
  });
  if (vcol_sq) {
    auto vq_ = vcol_q;
    auto qp_ = vcol_qprev;
    const Real qlo = vcol_qmin, qhi = vcol_qmax;
    par_for("m1_t2_vcq", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      const Real qn = vq_(m,k,j);
      vq_(m,k,j) = fmin(fmax(qn + r*(qn - qp_(m,k,j)), qlo), qhi);
      qp_(m,k,j) = qn;
    });
  }
  t2_vcprev = true;
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
  const Real rg = (vgd_on && vgd_alpha == x2_alpha) ? r : 0.0;
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

//----------------------------------------------------------------------------------------
//! \fn RadiationM1::Time2Rst*
//! \brief the restart channels: K1 (M1_T2_NK), then ipred2 (3, implicit_predictor), then
//! vet_prev (M1_T2_NVET, closure = vet_sc)

int RadiationM1::Time2RstNchWant() {
  if (time_scheme != M1_TIME_HESDIRK2) return 0;
  const int np = (impl_pord == 2) ? 5 : 3;
  return M1_T2_NK + (impl_pred ? np : 0) + (vet_sc ? M1_T2_NVET : 0)
         + ((vet_sc && vsc_every > 1) ? 2*M1_T2_NVET : 0)    // vet_sc_every: D0, D1
         + ((vet_col && t2_vcmode == 4) ? M1_X2_NCH : 0);    // rbfast-1021 extrap2
}

int RadiationM1::Time2RstNch() {
  return t2_ok ? Time2RstNchWant() : 0;
}

// gpumem-1022: m0 > 0 packs MeshBlocks m0 .. m0+nmb-1 into a(0 .. nmb-1) (per-MB staging)
void RadiationM1::Time2RstPack(DvceArray5D<Real> &a, int nmb, int m0) {
  auto mb = std::make_pair(0, nmb);
  auto ms = std::make_pair(m0, m0 + nmb);
  auto AL = Kokkos::ALL;
  int c = 0;
  Kokkos::deep_copy(Kokkos::subview(a, mb, std::make_pair(c, c+M1_T2_NK), AL, AL, AL),
                    Kokkos::subview(t2k1, ms, AL, AL, AL, AL));
  c += M1_T2_NK;
  if (impl_pred) {
    const int np = (impl_pord == 2) ? 5 : 3;
    Kokkos::deep_copy(Kokkos::subview(a, mb, std::make_pair(c, c+np), AL, AL, AL),
                      Kokkos::subview(ipred2, ms, AL, AL, AL, AL));
    c += np;
  }
  if (vet_sc) {
    Kokkos::deep_copy(Kokkos::subview(a, mb, std::make_pair(c, c+M1_T2_NVET), AL, AL, AL),
                      Kokkos::subview(vet_prev, ms, AL, AL, AL, AL));
    c += M1_T2_NVET;
    if (vsc_every > 1) {
      const int nv = M1_T2_NVET;
      Kokkos::deep_copy(Kokkos::subview(a, mb, std::make_pair(c, c+nv), AL, AL, AL),
                        Kokkos::subview(vsc_d0, ms, AL, AL, AL, AL));
      Kokkos::deep_copy(Kokkos::subview(a, mb, std::make_pair(c+nv, c+2*nv), AL, AL, AL),
                        Kokkos::subview(vsc_d1, ms, AL, AL, AL, AL));
      c += 2*nv;
    }
  }
  // rbfast-1021 extrap2: the raw builds of the last step, q at i = 0, time, alpha
  if (vet_col && t2_vcmode == 4) {
    Kokkos::deep_copy(Kokkos::subview(a, mb, std::make_pair(c, c+M1_X2_NTT), AL, AL, AL),
                      Kokkos::subview(x2_tt, ms, AL, AL, AL, AL));
    c += M1_X2_NTT;
    Kokkos::deep_copy(Kokkos::subview(a, mb, std::make_pair(c, c+M1_X2_NHR), AL, AL, AL),
                      Kokkos::subview(x2_hr, ms, AL, AL, AL, AL));
    c += M1_X2_NHR;
    auto a_ = a;
    auto q_ = x2_q;
    const Real t = x2_ok ? x2_t : -1.0e300, al = x2_alpha, ds = x2_dts;
    const int cq = c;
    const int mq = m0;
    const int na3 = static_cast<int>(a.extent(2)), na2 = static_cast<int>(a.extent(3));
    const int na1 = static_cast<int>(a.extent(4));
    par_for("m1_t2_x2rst", DevExeSpace(), 0, nmb - 1, 0, na3-1, 0, na2-1, 0, na1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      a_(m,cq,k,j,i) = (i == 0) ? q_(m+mq,k,j) : 0.0;
      a_(m,cq+1,k,j,i) = t;
      a_(m,cq+2,k,j,i) = al;
      a_(m,cq+3,k,j,i) = ds;
    });
  }
}

void RadiationM1::Time2RstSet(int ch, const HostArray4D<Real> &w, int nmb) {
  auto mb = std::make_pair(0, nmb);
  auto AL = Kokkos::ALL;
  if (ch < M1_T2_NK) {
    DeepCopyAcross(Kokkos::subview(t2k1, mb, ch, AL, AL, AL), w);
    return;
  }
  ch -= M1_T2_NK;
  if (impl_pred) {
    const int np = (impl_pord == 2) ? 5 : 3;
    if (ch < np) {
      DeepCopyAcross(Kokkos::subview(ipred2, mb, ch, AL, AL, AL), w);
      return;
    }
    ch -= np;
  }
  if (vet_sc) {
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
  // rbfast-1021 extrap2 (Time2RstPack)
  if (vet_col && t2_vcmode == 4) {
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
//! \fn void RadiationM1::Time2Report

void RadiationM1::Time2Report() {
  if (time_scheme != M1_TIME_HESDIRK2) return;
  if (global_variable::my_rank != 0) return;
  std::cout << "<rad_m1> time_scheme=hesdirk2: stage steps=" << t2_nstep
            << " backward-Euler steps=" << t2_nbe << " stage fallbacks=" << t2_nfall
            << " vet clips=" << t2_nclip;
  if (t2_k1reset) {
    std::cout << " k1_reset cells(rank 0)=" << t2_k1reset_n << " faces=" << t2_k1reset_nf;
    if (t2_k1pos) {std::cout << " (pos rule cells " << t2_k1pos_n << ")";}
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
  if (t2_pfix > 0) {
    std::cout << " time2_posfix: steps=" << t2_pfix_nstep << " re-solves="
              << t2_pfix_nsolve << " cells=" << t2_pfix_ncell;
  }
  if (t2_accfb) {std::cout << " accepted-with-fallback stages=" << t2_accfb_n;}
  // VetGdAlignDt calls that shortened dt (NewTimeStep calls, incl. the start-up ones)
  if (vgd_rot_align) {std::cout << " rotate_align shortened dt=" << vgd_rot_nalign;}
  std::cout << std::endl;
}

}  // namespace radm1
