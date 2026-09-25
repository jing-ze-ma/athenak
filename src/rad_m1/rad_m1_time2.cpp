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
#include <type_traits>
#include <utility>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "hydro/hydro.hpp"
#include "eos/eos.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_parfor.hpp"
#include "rad_m1/rad_m1_implicit.hpp"
#include "rad_m1/rad_m1_opacity.hpp"

namespace radm1 {

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
  t2_fvnew = false;
  if (sph_geom) {
    t2_fvnew = pin->GetOrAddBoolean("rad_m1", "time2_vstage",
                                    !global_variable::restart_run);
  } else if (pin->DoesParameterExist("rad_m1", "time2_vstage")) {
    t2_fvnew = pin->GetBoolean("rad_m1", "time2_vstage");
  }
  t2_vcmode = 0;
  if (vet_col) {
    const bool vdef = sph_geom && !global_variable::restart_run;
    std::string vm = (sph_geom || pin->DoesParameterExist("rad_m1", "time2_vet_col")) ?
        pin->GetOrAddString("rad_m1", "time2_vet_col", vdef ? "predict" : "lag") :
        std::string("lag");
    if (vm.compare("lag") == 0) {
      t2_vcmode = 0;
    } else if (vm.compare("predict") == 0) {
      t2_vcmode = 1;
    } else if (vm.compare("rebuild") == 0) {
      t2_vcmode = 2;
    } else if (vm.compare("extrap") == 0) {
      t2_vcmode = 3;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/time2_vet_col = '" << vm
                << "' is not a choice (lag | predict | rebuild | extrap)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (t2_vcmode != 0 && vcol_every != 1) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/time2_vet_col = " << vm
                << " needs vet_col_every = 1" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // time2_lin_tol / time2_lin_tol_fac (tests_m1/runs_5f_h2fast): the linear (Krylov)
  // tolerance of the two stage solves: time2_lin_tol when named, else implicit_lin_tol
  // x time2_lin_tol_fac.  Default 10 (the input files set implicit_lin_tol =
  // implicit_tol/100; the stage solves then take implicit_tol/10); 3-D He box: 21 -> 15
  // Krylov iterations per stage solve; the radwave G1 order is kept together with
  // time2_one_pass_safety = 30.  As for the runs_4a_accel levers, a restart whose file
  // does not carry the key keeps the old behaviour (1).
  const bool rs = global_variable::restart_run;
  t2_lin_tol = -1.0;
  if (pin->DoesParameterExist("rad_m1", "time2_lin_tol")) {
    t2_lin_tol = pin->GetReal("rad_m1", "time2_lin_tol");
  }
  t2_lin_fac = pin->GetOrAddReal("rad_m1", "time2_lin_tol_fac", rs ? 1.0 : 10.0);
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
  t2_dbg_fail = -1;
  if (pin->DoesParameterExist("rad_m1", "time2_dbg_fail")) {
    t2_dbg_fail = pin->GetInteger("rad_m1", "time2_dbg_fail");
  }

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int nmb = std::max((pmy_pack->nmb_thispack), (pmy_pack->pmesh->nmb_maxperrank));
  int n1 = indcs.nx1 + 2*(indcs.ng);
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
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
    Kokkos::realloc(vet_opac, nmb, M1_NOPAC, n3, n2, n1);
    Kokkos::deep_copy(vet_prev, 0.0);
  }
  if (vet_col) {
    // Time2VetStart's save slots (E, T) and the saved stage-start opacities;
    // time2_vet_col = predict | rebuild: also KT, T(Y1), E(Y1) and the inner face flux
    Kokkos::realloc(vet_now, nmb, (t2_vcmode == 0 || t2_vcmode == 3) ? 2 : 6, n3, n2, n1);
    Kokkos::realloc(vet_opac, nmb, M1_NOPAC, n3, n2, n1);
    if (t2_vcmode == 3) {
      Kokkos::realloc(vcol_prev, nmb, n3, n2, n1);
      Kokkos::realloc(vcol_qprev, nmb, n3, n2);
      t2_vcprev = false;
    }
  }
  t2_ok = false;
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1> time_scheme=hesdirk2 (g=" << kT2G << "; first step, restarts "
              << "without the slope and failed stages take backward Euler)" << std::endl;
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
  const Real b2 = g/(2.0*(1.0 - g));
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
    c1 = (1.0 - g)*dt;
    c2 = 0.0;
  } else {
    int nmb1 = static_cast<int>(u0.extent(0)) - 1;
    int n3 = static_cast<int>(u0.extent(2)), n2 = static_cast<int>(u0.extent(3));
    int n1 = static_cast<int>(u0.extent(4));
    par_for("m1_t2_avg", DevExeSpace(), 0, nmb1, 0, M1_NVAR-1, 0, n3-1, 0, n2-1,
            0, n1-1,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
      u0_(m,n,k,j,i) = 0.5*(u1_(m,n,k,j,i) + u0_(m,n,k,j,i));
    });
    auto avgf = [&](DvceArray4D<Real> &a, DvceArray4D<Real> &b) {
      auto a_ = a;
      auto b_ = b;
      par_for("m1_t2_avgf", DevExeSpace(), 0, static_cast<int>(a.extent(0))-1,
              0, static_cast<int>(a.extent(1))-1, 0, static_cast<int>(a.extent(2))-1,
              0, static_cast<int>(a.extent(3))-1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        a_(m,k,j,i) = 0.5*(b_(m,k,j,i) + a_(m,k,j,i));
      });
    };
    avgf(f0x1, t2f1);
    if (trans_on) {
      avgf(f0x2, t2f2);
      if (trans_x3) {avgf(f0x3, t2f3);}
    }
    c1 = 0.5*g*dt;
    c2 = (b2 - 0.5*g)*dt;
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
  t2_solve = (stage == 1) ? M1_T2S_STAGE1 : M1_T2S_STAGE2;
  t2_fail = false;
  nsub = 1;
  dt_sub = g*dt;
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
  hydro::Hydro *ph = pmy_pack->phydro;
  if (ph != nullptr) {
    Kokkos::deep_copy(DevExeSpace(), ph->u0, ph->u1);
    (void) ph->RestrictU(pdrive, 0);
    (void) ph->InitRecv(pdrive, -1);
    (void) ph->SendU(pdrive, 0);
    (void) ph->ClearSend(pdrive, -1);
    (void) ph->ClearRecv(pdrive, -1);
    (void) ph->RecvU(pdrive, 0);
    (void) ph->ApplyPhysicalBCs(pdrive, 0);
    (void) ph->Prolongate(pdrive, 0);
    (void) ph->ConToPrim(pdrive, 0);
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
  hydro::Hydro *ph = pmy_pack->phydro;
  const bool hh = (ph != nullptr);
  auto uh1 = hh ? ph->u1 : u1;
  const bool etg = hh ? ph->use_etotgrav : false;
  auto phicc = hh ? ph->phicc0 : arad_ref;
  const Real efl = e_floor;
  auto u0_ = u0;
  // the fast form (default): ONE fused pass over U^n (the hydro u1) gives T^n and the
  // opacities of U^n, with exactly the arithmetic of RadiationM1::Opacity, and the
  // stage-start opacities are saved and copied back afterwards instead of re-evaluated.
  // The debug opacity options (opac_freeze, dbg_opac_patch) take the general form.
  const bool fast = hh && !opac_zero && !opac_freeze && (dbg_opac_patch == 1.0) &&
                    !(opacity_type == M1_OPAC_TABLE && otab.nT <= 0);
  if (fast) {
    Kokkos::deep_copy(DevExeSpace(), vet_opac, opac);
    auto eos = ph->peos->eos_data;
    auto opac_ = opac;
    const int n1 = static_cast<int>(opac.extent(4));
    const int n2 = static_cast<int>(opac.extent(3));
    const int n3 = static_cast<int>(opac.extent(2));
    const int otype = opacity_type;
    const Real kp = kappa_p, kev = kappa_e, kf = kappa_f, kss = kappa_s;
    const Real rref = opac_rho_ref, tref = opac_t_ref, aa = opac_a, bb = opac_b;
    M1OpacTab ot = otab;
    par_for_lb("m1_t2_vsf", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real d = uh1(m,IDN,k,j,i);
      Real ke_dens = 0.5*(SQR(uh1(m,IM1,k,j,i)) + SQR(uh1(m,IM2,k,j,i)) +
                          SQR(uh1(m,IM3,k,j,i)))/fmax(d, 1.0e-300);
      Real eint = uh1(m,IEN,k,j,i) - ke_dens;
      if (etg) eint -= d*phicc(m,k,j,i);
      Real t = eos.Temperature(d, fmax(eint, 0.0));
      Real op, oe, of, os;
      if (otype == M1_OPAC_TABLE) {
        M1TableOpacities(ot, d, t, op, oe, of, os);
      } else {
        M1Opacities(otype, d, t, kp, kev, kf, kss, rref, tref, aa, bb, op, oe, of, os);
      }
      opac_(m,M1_OP_P,k,j,i) = d*op;
      opac_(m,M1_OP_E,k,j,i) = d*oe;
      opac_(m,M1_OP_T,k,j,i) = d*(of + os);
      if (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke) {
        vn_(m,0,k,j,i) = iw_(m,M1_IW_EN,k,j,i);
        vn_(m,1,k,j,i) = iw_(m,M1_IW_TP,k,j,i);
        iw_(m,M1_IW_EN,k,j,i) = fmax(u0_(m,M1_E,k,j,i), efl);
        iw_(m,M1_IW_TP,k,j,i) = t;
      }
    });
  } else if (hh) {
    auto eos = ph->peos->eos_data;
    par_for("m1_t2_vs0", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      vn_(m,0,k,j,i) = iw_(m,M1_IW_EN,k,j,i);
      vn_(m,1,k,j,i) = iw_(m,M1_IW_TP,k,j,i);
      iw_(m,M1_IW_EN,k,j,i) = fmax(u0_(m,M1_E,k,j,i), efl);
      Real dd = uh1(m,IDN,k,j,i);
      Real idd = 1.0/fmax(dd, 1.0e-300);
      Real ekin = 0.5*(SQR(uh1(m,IM1,k,j,i)) + SQR(uh1(m,IM2,k,j,i)) +
                       SQR(uh1(m,IM3,k,j,i)))*idd;
      Real egrv = etg ? (dd*phicc(m,k,j,i)) : 0.0;
      Real eg = uh1(m,IEN,k,j,i) - ekin - egrv;
      iw_(m,M1_IW_TP,k,j,i) = eos.Temperature(dd, fmax(eg, 1.0e-300));
    });
    std::swap(ph->u0, ph->u1);
    (void) Opacity(nullptr, 1);
    std::swap(ph->u0, ph->u1);
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
//! \fn void RadiationM1::Time2VetColAt
//! \brief closure = vet_col, time2_vet_col = predict | rebuild (m1-sp-order2b): the
//! formal solution of a state at t^{n+1} accurate to O(dt^2), in place of U^n.
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
  hydro::Hydro *ph = pmy_pack->phydro;
  const bool hh = (ph != nullptr);
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
  // the inner face flux (vet_col_order2 core rays): saved in slot 5 at i = is
  par_for("m1_t2_vcf", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    const Real f = f0_(m,k,j,is);
    vn_(m,5,k,j,is) = f;
    f0_(m,k,j,is) = two ? (2.0*f - t2f1_(m,k,j,is)) : (f + dt*k1_(m,M1_T2_F1,k,j,is));
  });
  if (hh) {
    Kokkos::deep_copy(DevExeSpace(), vet_opac, opac);
    auto uh = ph->u0;
    const bool etg = ph->use_etotgrav;
    auto phicc = ph->phicc0;
    auto eos = ph->peos->eos_data;
    auto opac_ = opac;
    const int otype = opacity_type;
    const Real kp = kappa_p, kev = kappa_e, kf = kappa_f, kss = kappa_s;
    const Real rref = opac_rho_ref, tref = opac_t_ref, aa = opac_a, bb = opac_b;
    M1OpacTab ot = otab;
    // m1-fast5-sp: an ideal-gas EOS with analytic opacities gets its own kernel (T and
    // the opacities by the same expressions, bitwise): the generic one carries the EOS
    // and opacity table code, 392 B of scratch per lane, even when both are off
    const Real gam = eos.gamma;
    const bool vcid = !eos.tbl.active && (otype != M1_OPAC_TABLE);
    auto vcb = [=] KOKKOS_FUNCTION (auto idl, const int m, const int k, const int j,
                                    const int i) M1_INL {
      const Real d = uh(m,IDN,k,j,i);
      Real t, e;
      if (two) {
        t = vn_(m,3,k,j,i);
        e = vn_(m,4,k,j,i);
      } else {
        // the hydro u0 is already the stage-1 OLD vector here (ImplicitSolve added
        // t2inc = (1-g) dt K1 to it): g dt K1 completes dt K1
        const Real m1 = uh(m,IM1,k,j,i) + gdt*k1_(m,M1_T2_M1,k,j,i);
        const Real m2 = uh(m,IM2,k,j,i) + gdt*k1_(m,M1_T2_M1+1,k,j,i);
        const Real m3 = uh(m,IM3,k,j,i) + gdt*k1_(m,M1_T2_M1+2,k,j,i);
        Real eint = uh(m,IEN,k,j,i) + gdt*k1_(m,M1_T2_EN,k,j,i)
                    - 0.5*(m1*m1 + m2*m2 + m3*m3)/fmax(d, 1.0e-300);
        if (etg) eint -= d*phicc(m,k,j,i);
        if constexpr (decltype(idl)::value) {
          t = ((gam-1.0)*fmax(eint, 0.0)/d);
        } else {
          t = eos.Temperature(d, fmax(eint, 0.0));
        }
        e = u0_(m,M1_E,k,j,i) + dt*k1_(m,M1_T2_E,k,j,i);
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
      opac_(m,M1_OP_T,k,j,i) = d*(of + os);
      vn_(m,0,k,j,i) = iw_(m,M1_IW_EN,k,j,i);
      vn_(m,1,k,j,i) = iw_(m,M1_IW_TP,k,j,i);
      vn_(m,2,k,j,i) = iw_(m,M1_IW_KT,k,j,i);
      iw_(m,M1_IW_EN,k,j,i) = fmax(e, efl);
      iw_(m,M1_IW_TP,k,j,i) = t;
      iw_(m,M1_IW_KT,k,j,i) = d*(of + os);
    };
    if (vcid) {
      par_for_lb("m1_t2_vcp", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) M1_INL {
        vcb(std::true_type{}, m, k, j, i);
      });
    } else {
      par_for_lb("m1_t2_vcp", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) M1_INL {
        vcb(std::false_type{}, m, k, j, i);
      });
    }
  } else {
    par_for("m1_t2_vcr", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      vn_(m,0,k,j,i) = iw_(m,M1_IW_EN,k,j,i);
      const Real e = two ? vn_(m,4,k,j,i)
                         : (u0_(m,M1_E,k,j,i) + dt*k1_(m,M1_T2_E,k,j,i));
      iw_(m,M1_IW_EN,k,j,i) = fmax(e, efl);
    });
  }
  VetColBuild();
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
    if (i == is) {
      f0_(m,k,j,is) = vn_(m,5,k,j,is);
    }
  });
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
//! \fn RadiationM1::Time2Rst*
//! \brief the restart channels: K1 (M1_T2_NK), then ipred2 (3, implicit_predictor), then
//! vet_prev (M1_T2_NVET, closure = vet_sc)

int RadiationM1::Time2RstNchWant() {
  if (time_scheme != M1_TIME_HESDIRK2) return 0;
  const int np = (impl_pord == 2) ? 5 : 3;
  return M1_T2_NK + (impl_pred ? np : 0) + (vet_sc ? M1_T2_NVET : 0)
         + ((vet_sc && vsc_every > 1) ? 2*M1_T2_NVET : 0);   // vet_sc_every: D0, D1
}

int RadiationM1::Time2RstNch() {
  return t2_ok ? Time2RstNchWant() : 0;
}

void RadiationM1::Time2RstPack(DvceArray5D<Real> &a, int nmb) {
  auto mb = std::make_pair(0, nmb);
  auto AL = Kokkos::ALL;
  int c = 0;
  Kokkos::deep_copy(Kokkos::subview(a, mb, std::make_pair(c, c+M1_T2_NK), AL, AL, AL),
                    Kokkos::subview(t2k1, mb, AL, AL, AL, AL));
  c += M1_T2_NK;
  if (impl_pred) {
    const int np = (impl_pord == 2) ? 5 : 3;
    Kokkos::deep_copy(Kokkos::subview(a, mb, std::make_pair(c, c+np), AL, AL, AL),
                      Kokkos::subview(ipred2, mb, AL, AL, AL, AL));
    c += np;
  }
  if (vet_sc) {
    Kokkos::deep_copy(Kokkos::subview(a, mb, std::make_pair(c, c+M1_T2_NVET), AL, AL, AL),
                      Kokkos::subview(vet_prev, mb, AL, AL, AL, AL));
    c += M1_T2_NVET;
    if (vsc_every > 1) {
      const int nv = M1_T2_NVET;
      Kokkos::deep_copy(Kokkos::subview(a, mb, std::make_pair(c, c+nv), AL, AL, AL),
                        Kokkos::subview(vsc_d0, mb, AL, AL, AL, AL));
      Kokkos::deep_copy(Kokkos::subview(a, mb, std::make_pair(c+nv, c+2*nv), AL, AL, AL),
                        Kokkos::subview(vsc_d1, mb, AL, AL, AL, AL));
    }
  }
}

void RadiationM1::Time2RstSet(int ch, const HostArray4D<Real> &w, int nmb) {
  auto mb = std::make_pair(0, nmb);
  auto AL = Kokkos::ALL;
  if (ch < M1_T2_NK) {
    Kokkos::deep_copy(Kokkos::subview(t2k1, mb, ch, AL, AL, AL), w);
    return;
  }
  ch -= M1_T2_NK;
  if (impl_pred) {
    const int np = (impl_pord == 2) ? 5 : 3;
    if (ch < np) {
      Kokkos::deep_copy(Kokkos::subview(ipred2, mb, ch, AL, AL, AL), w);
      return;
    }
    ch -= np;
  }
  if (vet_sc) {
    if (ch < M1_T2_NVET) {
      Kokkos::deep_copy(Kokkos::subview(vet_prev, mb, ch, AL, AL, AL), w);
      return;
    }
    ch -= M1_T2_NVET;
    if (vsc_every > 1) {
      if (ch < M1_T2_NVET) {
        Kokkos::deep_copy(Kokkos::subview(vsc_d0, mb, ch, AL, AL, AL), w);
      } else {
        Kokkos::deep_copy(Kokkos::subview(vsc_d1, mb, ch - M1_T2_NVET, AL, AL, AL), w);
      }
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
  if (vet_col && t2_vcmode != 0) {
    std::cout << " time2_vet_col=" << ((t2_vcmode == 1) ? "predict" :
                                       ((t2_vcmode == 2) ? "rebuild" : "extrap"))
              << " fallbacks=" << t2_vcnfb;
  }
  std::cout << std::endl;
}

}  // namespace radm1
