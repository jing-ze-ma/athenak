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
#include <utility>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "hydro/hydro.hpp"
#include "eos/eos.hpp"
#include "rad_m1/rad_m1.hpp"
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
  if (tau_closure) {
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
  t2_solve = M1_T2S_NONE;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::Time2VetStart
//! \brief closure = vet_sc, stage 1: the formal solution at U^n -- E^n, T^n and the
//! opacities of U^n (the hydro's u1) -- and then the extrapolated tensor
//! D* = D^n + (dt/dt_prev)(D^n - D^{n-1}), which both stage solves read.  Called from
//! ImplicitSolve at the point where the backward-Euler step calls VetShortChar; the
//! solve's own EN (old vector), TP (start T) and opacities are put back afterwards.

void RadiationM1::Time2VetStart() {
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
    par_for("m1_t2_vsf", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
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
  VetShortChar();
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
  Time2VetExtrapolate();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::Time2VetExtrapolate
//! \brief vet_cell <- D* = D^n + r (D^n - D^{n-1}), r = dt/dt_prev, for the uniaxial
//! (chi, n) and the full D_ab; a cell whose D* is not realizable keeps D^n (counted).
//! vet_prev <- D^n for the next step.

void RadiationM1::Time2VetExtrapolate() {
  auto vc_ = vet_cell;
  auto vp_ = vet_prev;
  const Real dt = pmy_pack->pmesh->dt;
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
//! \fn RadiationM1::Time2Rst*
//! \brief the restart channels: K1 (M1_T2_NK), then ipred2 (3, implicit_predictor), then
//! vet_prev (M1_T2_NVET, closure = vet_sc)

int RadiationM1::Time2RstNchWant() {
  if (time_scheme != M1_TIME_HESDIRK2) return 0;
  const int np = (impl_pord == 2) ? 5 : 3;
  return M1_T2_NK + (impl_pred ? np : 0) + (vet_sc ? M1_T2_NVET : 0);
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
    Kokkos::deep_copy(Kokkos::subview(vet_prev, mb, ch, AL, AL, AL), w);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::Time2Report

void RadiationM1::Time2Report() {
  if (time_scheme != M1_TIME_HESDIRK2) return;
  if (global_variable::my_rank != 0) return;
  std::cout << "<rad_m1> time_scheme=hesdirk2: stage steps=" << t2_nstep
            << " backward-Euler steps=" << t2_nbe << " stage fallbacks=" << t2_nfall
            << " vet clips=" << t2_nclip << std::endl;
}

}  // namespace radm1
