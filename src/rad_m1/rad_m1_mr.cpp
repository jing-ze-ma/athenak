//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_mr.cpp
//! \brief <rad_m1>/implicit_mr_every = k: MULTI-RATE implicit radiation (m1-fast4,
//! tests_m1/runs_5r_fast4).
//!
//! The implicit radiation (transport + gas coupling) is unconditionally stable, and where
//! its state follows the hydro quasi-statically nothing requires it to be solved on
//! every hydro step.  With k > 1 the hydro (Heun) takes every step ALONE, and the
//! radiation R advances the gas + radiation over a whole window Delta of k hydro steps,
//! operator split:
//!   H^{k/2} R(Delta) H^k R(Delta) H^k ...   =   (H(Delta/2) R(Delta) H(Delta/2))^N,
//! i.e. R sits at the window CENTRES (the first one after k/2 steps, with Delta twice
//! the elapsed time), which is Strang splitting: second order in time when R is.  At a
//! time that is a window END (t = N Delta from the start) the state is the Strang one;
//! in between, it is off by a non-accumulating half-window splitting term.
//!
//! R(Delta) is the two-stage, L-stable, stiffly accurate SDIRK2 (g = 1 - 1/sqrt 2):
//!   Y_A = Y_0 + g Delta K_A,  Y_B = Y_0 + (1-g) Delta K_A + g Delta K_B,  Y_1 = Y_B,
//! built from the hesdirk2 stage machinery (rad_m1_time2.cpp): stage A is a STAGE1
//! solve (dt = g Delta, t2inc = 0, K_A stored in t2k2), stage B a STAGE2 solve from
//! Y_A with t2inc = (1 - 2g) Delta K_A.  The closure tensor is built once per window at
//! Y_0, as the hesdirk2 step builds it at U^n.  A stage that fails puts Y_0 back and
//! takes R(Delta) as backward Euler (counted).
//!
//! Needs time_scheme = hesdirk2 (its stage arrays).  Default k = 1: nothing here runs,
//! and every other path is bitwise unchanged.

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
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

namespace radm1 {

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::MRInit
//! \brief read <rad_m1>/implicit_mr_every (read only when named; 1 = off)

void RadiationM1::MRInit(ParameterInput *pin) {
  // vet_sc_every (read only when named): see VscSkip below
  vsc_every = 1;
  if (vet_sc && pin->DoesParameterExist("rad_m1", "vet_sc_every")) {
    vsc_every = std::max(1, pin->GetInteger("rad_m1", "vet_sc_every"));
    if (vsc_every > 1) {
      if (time_scheme != M1_TIME_HESDIRK2) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<rad_m1>/vet_sc_every > 1 needs time_scheme = hesdirk2"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      const int nmb = pmy_pack->nmb_thispack;
      Kokkos::realloc(vsc_d0, nmb, M1_T2_NVET, static_cast<int>(vet_prev.extent(2)),
                      static_cast<int>(vet_prev.extent(3)),
                      static_cast<int>(vet_prev.extent(4)));
      Kokkos::realloc(vsc_d1, nmb, M1_T2_NVET, static_cast<int>(vet_prev.extent(2)),
                      static_cast<int>(vet_prev.extent(3)),
                      static_cast<int>(vet_prev.extent(4)));
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_m1> vet_sc_every = " << vsc_every << ": formal solution on "
                  << "every " << vsc_every << "-th tensor build, linear extrapolation "
                  << "in time in between" << std::endl;
      }
    }
  }
  mr_every = 1;
  if (pin->DoesParameterExist("rad_m1", "implicit_mr_every")) {
    mr_every = pin->GetInteger("rad_m1", "implicit_mr_every");
  }
  if (mr_every <= 1) {
    mr_every = 1;
    return;
  }
  if (time_scheme != M1_TIME_HESDIRK2 || (mr_every % 2) != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<rad_m1>/implicit_mr_every = " << mr_every << " needs time_scheme = "
              << "hesdirk2 and an EVEN k (Strang placement)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // implicit_mr_theta > 0: the ADAPTIVE window.  theta = the largest relative change of
  // the gas internal energy that a radiation step makes (active cells); above theta the
  // next window is halved (down to 2 steps), below theta/4 doubled (up to k).  Read
  // only when named.
  if (pin->DoesParameterExist("rad_m1", "implicit_mr_theta")) {
    mr_theta = pin->GetReal("rad_m1", "implicit_mr_theta");
  }
  mr_first = true;
  mr_cnt = 0;
  mr_acc = 0.0;
  mr_kc = 0;
  mr_kn = mr_every;
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1> implicit_mr_every = " << mr_every << ": multi-rate radiation "
              << "(operator-split SDIRK2 over " << mr_every << " hydro steps, Strang)"
              << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::MRStep
//! \brief called by the Driver after the Heun stages of every step: count the step, and
//! take R(Delta) when the window centre is reached

void RadiationM1::MRStep(Driver *pdrive, Real dt) {
  mr_acc += dt;
  mr_cnt += 1;
  // R sits at the window centres: the distance between two of them is half the last
  // window plus half the next one (k steps each for a fixed k)
  const int need = (mr_kc + mr_kn)/2;
  if (mr_cnt < need) return;
  // Strang: the first window is half long; R(Delta) covers the whole window
  const Real dlt = mr_first ? (2.0*mr_acc) : mr_acc;
  mr_first = false;
  mr_cnt = 0;
  mr_acc = 0.0;
  mr_kc = mr_kn;
  mr_ksum += static_cast<Real>(mr_kc);
  MRSolve(pdrive, dlt);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::MRSolve
//! \brief R(dlt): the SDIRK2 radiation step (see the file header)

void RadiationM1::MRSolve(Driver *pdrive, Real dlt) {
  Mesh *pmesh = pmy_pack->pmesh;
  hydro::Hydro *ph = pmy_pack->phydro;
  const Real g = 1.0 - 1.0/std::sqrt(2.0);
  mr_on = true;
  mr_dt = dlt;
  // Y_0, kept for a failed stage: radiation in u1 and the t2f* face copies, gas in the
  // hydro u1 (free between steps; stage A's formal solution also reads U^n there)
  Kokkos::deep_copy(DevExeSpace(), u1, u0);
  Kokkos::deep_copy(DevExeSpace(), t2f1, f0x1);
  if (trans_on) {
    Kokkos::deep_copy(DevExeSpace(), t2f2, f0x2);
    if (trans_x3) {Kokkos::deep_copy(DevExeSpace(), t2f3, f0x3);}
  }
  if (ph != nullptr) {Kokkos::deep_copy(DevExeSpace(), ph->u1, ph->u0);}
  auto run_lists = [&]() {
    pdrive->ExecuteTaskList(pmesh, "m1_before_stagen", 1);
    pdrive->ExecuteTaskList(pmesh, "m1_stagen", 1);
    pdrive->ExecuteTaskList(pmesh, "m1_after_stagen", 1);
  };
  // stage A: old vector = Y_0, dt = g Delta, K_A -> t2k2
  Kokkos::deep_copy(DevExeSpace(), t2inc, 0.0);
  t2_solve = M1_T2S_STAGE1;
  t2_fail = false;
  nsub = 1;
  dt_sub = g*dlt;
  run_lists();
  if (!t2_fail) {
    // stage B from Y_A: old vector = Y_0 + (1-g) Delta K_A = Y_A + (1-2g) Delta K_A
    const Real c = (1.0 - 2.0*g)*dlt;
    auto k2_ = t2k2;
    auto in_ = t2inc;
    const int nmb1 = static_cast<int>(t2inc.extent(0)) - 1;
    const int n3 = static_cast<int>(t2inc.extent(2));
    const int n2 = static_cast<int>(t2inc.extent(3));
    const int n1 = static_cast<int>(t2inc.extent(4));
    par_for("m1_mr_inc", DevExeSpace(), 0, nmb1, 0, M1_T2_NK-1, 0, n3-1, 0, n2-1,
            0, n1-1,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
      in_(m,n,k,j,i) = c*k2_(m,n,k,j,i);
    });
    t2_solve = M1_T2S_STAGE2;
    t2_fail = false;
    nsub = 1;
    dt_sub = g*dlt;
    run_lists();
  }
  if (t2_fail) {
    // a stage was not admissible: Y_0 back, R(Delta) as backward Euler, cold predictor
    Kokkos::deep_copy(DevExeSpace(), u0, u1);
    Kokkos::deep_copy(DevExeSpace(), f0x1, t2f1);
    if (trans_on) {
      Kokkos::deep_copy(DevExeSpace(), f0x2, t2f2);
      if (trans_x3) {Kokkos::deep_copy(DevExeSpace(), f0x3, t2f3);}
    }
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
    pred_ok = false;
    pred2_ok = false;
    t2_vprev = false;
    mr_nfall += 1.0;
    t2_solve = M1_T2S_NONE;
    int ns = SetSubsteps(dlt);
    for (int n = 0; n < ns; ++n) {
      run_lists();
      (void) ApplyClosureLimits(pdrive, 1);
    }
  }
  t2_solve = M1_T2S_NONE;
  (void) ApplyClosureLimits(pdrive, 1);
  (void) NewTimeStep(pdrive, 1);
  // the guard: theta = max |e_int(Y_1) - e_int(Y_0)|/e_int(Y_0) over the active cells
  if (mr_theta > 0.0 && ph != nullptr) {
    auto &indcs = pmesh->mb_indcs;
    const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
    const int ks = indcs.ks, ke = indcs.ke;
    const int nmb1 = pmy_pack->nmb_thispack - 1;
    auto ua = ph->u0;
    auto ub = ph->u1;
    const bool etg = ph->use_etotgrav;
    auto phi = ph->phicc0;
    const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
    Real th = 0.0;
    Kokkos::parallel_reduce("m1_mr_theta",
    Kokkos::RangePolicy<DevExeSpace>(DevExeSpace(), 0, (nmb1 + 1)*nkji),
    KOKKOS_LAMBDA(const int idx, Real &lmx) {
      int m = idx/nkji;
      int r = idx - m*nkji;
      int k = r/nji;
      r -= k*nji;
      int j = r/ni;
      int i = r - j*ni;
      k += ks; j += js; i += is;
      Real ei[2];
      for (int q = 0; q < 2; ++q) {
        auto &u = (q == 0) ? ua : ub;
        Real d = u(m,IDN,k,j,i);
        Real ek = 0.5*(SQR(u(m,IM1,k,j,i)) + SQR(u(m,IM2,k,j,i)) +
                       SQR(u(m,IM3,k,j,i)))/fmax(d, 1.0e-300);
        ei[q] = u(m,IEN,k,j,i) - ek - (etg ? d*phi(m,k,j,i) : 0.0);
      }
      Real v = fabs(ei[0] - ei[1])/fmax(fabs(ei[1]), 1.0e-300);
      lmx = (v > lmx) ? v : lmx;
    }, Kokkos::Max<Real>(th));
#if MPI_PARALLEL_ENABLED
    {Real g;
    MPI_Allreduce(&th, &g, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
    th = g;}
#endif
    mr_thmax = std::max(mr_thmax, th);
    if (th > mr_theta) {
      mr_kn = std::max(2, 2*((mr_kc/2 + 1)/2));   // half the window, even, >= 2
    } else if (th < 0.25*mr_theta) {
      mr_kn = std::min(mr_every, 2*mr_kc);
    } else {
      mr_kn = mr_kc;
    }
  }
  t2_dtprev = dlt;
  // the stage state (K1, ipred2, vet_prev) is then written to restart files
  t2_ok = true;
  mr_nr += 1.0;
  mr_on = false;
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationM1::VscSkip
//! \brief vet_sc_every = N: called by Time2VetStart for every tensor build (once per
//! step, or per multi-rate radiation step).  Every N-th build, and while fewer than two
//! formal solutions are stored, runs the formal solution (returns false;
//! Time2VetStart then stores it with VscStore).  The others put D* = D1 + r (D1 - D0),
//! r = (t - t1)/(t1 - t0), into vet_cell (a cell whose D* is not realizable keeps D1,
//! counted) and return true.  D is a smooth function of the state, so the error is
//! O((N dt)^2) (second order); the tensor was lagged by a step anyway.

bool RadiationM1::VscSkip() {
  Mesh *pmesh = pmy_pack->pmesh;
  vsc_tn = pmesh->time + (mr_on ? pmesh->dt : 0.0);
  const bool skip = (vsc_nb >= 2) && ((vsc_cnt % vsc_every) != 0) && (vsc_t1 > vsc_t0);
  vsc_cnt += 1;
  if (!skip) return false;
  const Real r = (vsc_tn - vsc_t1)/(vsc_t1 - vsc_t0);
  auto vc_ = vet_cell;
  auto d0_ = vsc_d0;
  auto d1_ = vsc_d1;
  const bool full = vet_full;
  auto &indcs = pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n3 = static_cast<int>(vsc_d0.extent(2));
  const int n2 = static_cast<int>(vsc_d0.extent(3));
  const int n1 = static_cast<int>(vsc_d0.extent(4));
  const int n21 = n2*n1, n321 = n3*n21;
  Real nclip = 0.0;
  Kokkos::parallel_reduce("m1_vsc_ext",
  Kokkos::RangePolicy<DevExeSpace>(DevExeSpace(), 0, (nmb1 + 1)*n321),
  KOKKOS_LAMBDA(const int idx, Real &lc) {
    int m = idx/n321;
    int q = idx - m*n321;
    int k = q/n21;
    q -= k*n21;
    int j = q/n1;
    int i = q - j*n1;
    Real ds[M1_T2_NVET];
    for (int c = 0; c < M1_T2_NVET; ++c) {
      ds[c] = d1_(m,c,k,j,i) + r*(d1_(m,c,k,j,i) - d0_(m,c,k,j,i));
    }
    // the realizability test of Time2VetExtrapolate
    bool ok = (ds[0] >= 0.0) && (ds[0] <= 1.0);
    Real nn = sqrt(ds[1]*ds[1] + ds[2]*ds[2] + ds[3]*ds[3]);
    ok = ok && (nn > 0.5);
    if (ok) {
      ds[1] /= nn;
      ds[2] /= nn;
      ds[3] /= nn;
    }
    if (full) {
      Real a = ds[4], b = ds[5], c = ds[6], d = ds[7], e = ds[8], f = ds[9];
      Real m2 = a*b - d*d;
      Real det = a*(b*c - f*f) - d*(d*c - f*e) + e*(d*f - b*e);
      ok = ok && (a >= 0.0) && (b >= 0.0) && (c >= 0.0) && (a <= 1.0) && (b <= 1.0) &&
           (c <= 1.0) && (m2 >= 0.0) && (det >= 0.0);
    }
    for (int c = 0; c < M1_T2_NVET; ++c) {
      vc_(m,M1_VET_CHI+c,k,j,i) = ok ? ds[c] : d1_(m,c,k,j,i);
    }
    if (!ok && i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke) {
      lc += 1.0;
    }
  }, Kokkos::Sum<Real>(nclip));
  vsc_nclip += nclip;
  vsc_nskip += 1.0;
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VscStore
//! \brief vet_sc_every: D0 <- D1, D1 <- the formal solution just made (vet_cell), times

void RadiationM1::VscStore() {
  std::swap(vsc_d0, vsc_d1);
  auto AL = Kokkos::ALL;
  Kokkos::deep_copy(DevExeSpace(), vsc_d1,
                    Kokkos::subview(vet_cell, AL, std::make_pair(M1_VET_CHI,
                                    M1_VET_CHI + M1_T2_NVET), AL, AL, AL));
  vsc_t0 = vsc_t1;
  vsc_t1 = vsc_tn;
  vsc_nb = std::min(vsc_nb + 1, 2);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::MRReport
//! \brief one line at the end of the run (rank 0; called from ImplicitReport)

void RadiationM1::MRReport() {
  if (vsc_every > 1) {
    std::cout << "<rad_m1> vet_sc_every = " << vsc_every << ": builds=" << vsc_cnt
              << " extrapolated=" << vsc_nskip << " clipped cells=" << vsc_nclip
              << std::endl;
  }
  if (mr_every <= 1) return;
  std::cout << "<rad_m1> implicit_mr_every = " << mr_every << ": radiation steps="
            << mr_nr << " mean window=" << ((mr_nr > 0.0) ? mr_ksum/mr_nr : 0.0)
            << " backward-Euler fallbacks=" << mr_nfall;
  if (mr_theta > 0.0) {
    std::cout << " theta max=" << mr_thmax << " (guard " << mr_theta << ")";
  }
  std::cout << std::endl;
}

} // namespace radm1
