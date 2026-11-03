#ifndef UTILS_CK_2ND_HPP_
#define UTILS_CK_2ND_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file ck_2nd.hpp
//! \brief ck2nd-1102: the pieces that make the correlated-k heating SECOND ORDER in time
//! under <problem>/ck_time_scheme = hesdirk2.  Every key is read only when named; absent,
//! nothing here is allocated or launched and the code is bitwise the code without it.
//! Included by utils/two_stream_rt_driver.cpp and the deep_hot_jupiter_rt pgen only (not
//! a wide header).
//!
//! (1) problem/ck_impl_tol_inc = tau > 0: THE NEWTON TEST RELATIVE TO THE STEP.  The
//!     residual test max_i |R_i|/(e_i + eps e_max) <= ck_impl_tol leaves an error of up
//!     to tol e per call, i.e. a per-unit-time error that grows like 1/dt (measured: the
//!     frozen ladder does not converge at tol 1e-8).  With the key, cell i passes when
//!       |R_i| <= tau |e_i - e*_i| + ck_impl_tol_floor |e_i|
//!                  + ck_impl_tol_xfloor bdt (|S_i| + 2 E_i),
//!     so the error per call is tau times the increment the call applies (a fixed small
//!     fraction of the heating, independent of dt) plus round-off floors: 1e-15 of the
//!     cell's OWN energy and 1e-12 of its gross exchange bdt (A + E) over the call (the
//!     scale of the round-off of S = A - E where A ~ E; without it measured: cells in
//!     radiative balance never pass, every call runs to maxit).  The eps e_max share of today's norm is dropped
//!     on purpose: in a thin top cell it is 1e5-1e9 x e, so tol 1e-13 there still allows
//!     ~1e-8 e per call, which accumulates over the 1/dt calls into a dt-independent
//!     error floor (measured: the frozen ladder's finest members stop converging).  Implemented as the per-cell norm s = tol |R| / (that bound) in the fused
//!     T4 kernel, so every other use of tol (the KKT test, the column skip) is unchanged.
//! (2) problem/ck_t2_seed = true (hesdirk2): THE OPACITY AT A SECOND-ORDER STATE.  The
//!     column solve freezes kappa (and the beam) at its first iterate.  Without a seed
//!     that is the stage's old vector, O(dt) away from the stage value, so the stage is
//!     solved with an O(dt) opacity error and the step is first order.  The seed starts
//!     the Newton at old + g dt K (stage 1: K1, stage 2: K2), O(dt^2) away; e* is still
//!     recorded net of it (the ck_impl_warm mechanism), so the balance solved is the same
//!     and only the frozen opacity moves.  Full (non-partial) calls only.
//! (3) problem/ck_cad_drift = true (ck_impl_every > 1): THE LINEARISED STEP WITH A RATE
//!     DRIFT.  CkCadLin's rate Q0 + D (T - T0) holds the rate the full call measured, plus
//!     the cell's own diagonal response; what the other cells (and rho) do to it since is
//!     lagged, an O(N dt) rate error -> first order.  Each full store now also measures
//!     the drift G = (Q0_new - Q0_old - D_old (T0_new - T0_old)) / (t_new - t_old) of
//!     the part the diagonal model missed, and the linearised step uses
//!     Q0 + G (t - t0) + D (T - T0): an O((N dt)^2) rate error -> second order.  The
//!     evaluation time t is pm->time + pm->dt for every call (both hesdirk2 stages have
//!     c = 1).  Not in the restart file: after a restart G = 0 until the second full call.
//! (4b) problem/ck_vef_dt = D > 0 [s]: the VEF refresh counted in SIMULATED TIME: a storing
//!     pass refreshes when t + dt >= t_last + D (instead of every ck_vef_every cycles),
//!     so the lag of the factors no longer scales with the number of cycles.  t_last is
//!     not in the restart file: a restarted run refreshes on its first storing pass.
//! (5) problem/ck_impl_xpass = n > 0: AN EXTRA NEWTON PASS.  Once a call has converged
//!     by the usual test, n more passes are run on every column with the step forced
//!     (no column is accepted as it stands), each multiplying the remaining error by the
//!     Newton contraction (<< 1), so the error a call leaves is ~ contraction^n x tol
//!     instead of tol, at the cost of n sweeps per call.  The usual stop test then ends
//!     the call.
//! (3b) problem/ck_cad_tri = true (with ck_cad_drift): THE LINEARISED STEP IMPLICIT IN
//!     THE NEIGHBOURS TOO.  The diagonal step treats the coupling to the adjacent cells
//!     explicitly (through the drift), which in the stiff top (radiative exchange time <<
//!     dt) is first order whatever the extrapolation (measured: drift alone 10x smaller
//!     error, ratio 5-10).  The full store also keeps the Newton tridiagonal couplings
//!     J0, J2 (dS_i/dT_{i-1,i+1}, thick cells), and each linearised step solves the
//!     linearised BE balance with them (one Thomas sweep per column).
//! (6) problem/ck_top_live = true: THE TOP SLOT FROM THE CURRENT ITERATE.  The band
//!     solver's slot i = ie+1 (the unresolved column above, rt_top_clamp = false) reads
//!     the hydro GHOST, which ApplyPhysicalBCs filled at the end of the previous stage:
//!     inside an implicit stage the top boundary is lagged by one stage, an O(dt) error in
//!     the top cell's heating that the hesdirk2 tableau does not cancel (first order).
//!     With the key the user boundary function is re-run before every Newton pass, so the
//!     ghost follows the iterate (same boundary model, no lag).
//! (4) problem/ck_vef_tint = true (ck_sph_face 5/6): THE VEF FACTORS EXTRAPOLATED IN
//!     TIME.  The formal solution refreshes s = sqrt(3f) and the face factors every
//!     ck_vef_every cycles and holds them in between: an O(N dt) lag -> first order.  With
//!     the key the last two refreshes (t_P, t_K) are kept and every storing pass in
//!     between uses x_K + (x_K - x_P)(t - t_K)/(t_K - t_P) per element (x_K where the two
//!     differ in sign or one is 0): an O((N dt)^2) error.  Not in the restart file: after
//!     a restart the factors are held until the second refresh.

#include <cmath>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"

class Mesh;

namespace two_stream_rt {

inline Real ck_tol_inc = 0.0;          // (1) 0 = off
inline Real ck_tol_floor = 1.0e-15;
inline Real ck_tol_xfloor = 1.0e-12;
inline bool ck_t2_seed = false;        // (2)
inline bool ck_seed_ext_on = false;    //     set around a seeded solve
inline DvceArray4D<Real> *ck_seed_ext_ptr = nullptr;   // the seed [energy/vol]
inline bool ck_cad_drift = false;      // (3)
inline bool ck_cad_tri = false;        // (3b) the tridiagonal linearised step
// (3c) ck2ndt-1102 problem/ck_cad_drift2 = true (with ck_cad_drift, not ck_cad_tri): the
//     drift G of (3) is the MEAN slope over the last store interval, i.e. the slope at its
//     midpoint, an O(N dt) error in G -> an O((N dt)^2) rate error.  With the key the
//     last two mean slopes give the curvature H, G is moved to the store time
//     (G + H tl/2) and the model rate is Q0 + G tau + H tau^2/2: O((N dt)^3).
inline bool ck_cad_drift2 = false;
// (3d) ck2ndt-1102 problem/ck_cad_lin_iter = n (default 0): the linearised (model) step is
//     today ONE linearisation of the backward-Euler model step in T (e - e_old = cv dT).
//     With n > 0 the model step e - e_old = dt (Q + D (T(e) - T0)) is solved by n more
//     Newton iterations in e with the EOS T(e); not with ck_cad_tri.
inline int ck_cad_lin_iter = 0;
inline DvceArray5D<Real> *ck_cad2_ptr = nullptr;       // (m,2,k,j,i): G, t0
inline Real ck_cad_teval = 0.0;        //     the evaluation time of this call
inline int64_t ck_cad_ndrift = 0;      //     stores that measured a drift
inline bool ck_top_live = false;       // (6)
inline int ck_impl_xpass = 0;          // (5) extra forced Newton passes after convergence
inline bool ck_force_step = false;     //     set for such a pass
inline int64_t ck_impl_nxpass = 0;
inline bool ck_vef_tint = false;       // (4)
inline Real ck_vef_dt = 0.0;           // (4b) VEF refresh every ck_vef_dt SIMULATED seconds
inline Real ck_vef_tlast = -1.0e300;   //      time of the last refresh
inline int64_t ck_vef_ntint = 0;       //     storing passes that extrapolated
// (7) ck2ndt-1102 problem/ck_cad_dt = P [s] (ck_impl_every > 1): THE CADENCE IN SIMULATED
//     TIME.  The full call is made on the step whose end t + dt passes a multiple of P
//     (absolute time, so a restart keeps the schedule: stateless), linearised steps in
//     between; with a dt that divides P every ladder member refreshes at the same instants,
//     so the cadence's model error is the same function of time in every member.
//     problem/ck_xstep_dt = X [s] (a multiple of P; with ck_impl_xstep > 0): the stored
//     operator is re-formed on the full call that passes a multiple of X, re-applied on
//     the others.  With ck_cad_dt, problem/ck_vef_dt = V (a multiple of X, or of P when
//     ck_xstep_dt is absent) refreshes the VEF on the storing pass that passes a multiple
//     of V (same stateless rule).
inline Real ck_cad_dt = 0.0;
inline Real ck_xstep_dt = 0.0;
inline Real ck_cad_t0 = 0.0;   // problem/ck_cad_t0: the schedule is t0 + n P (default 0)
// a point t0 + n P in (t, t + dt] (to 1e-3 dt)
inline bool CkSchedCross(Real t, const Real dt, const Real p) {
  t -= ck_cad_t0;
  const Real eps = 1.0e-3*dt;
  return std::floor((t + dt + eps)/p) > std::floor((t + eps)/p);
}
// (4) called by the pass after the refresh block: refreshed = this pass re-formed the
// factors, storing = this pass forms the stored operator (defined in the rt driver)
void CkVefTint(Mesh *pm, const bool refreshed, const bool storing);

// (8) ck2ndt-1102 problem/ck_cad_probe = true (hesdirk2, ck_impl_every > 1,
//     ck_t2_start = probe): THE CADENCE REFRESH AS A PROBE AT A STEP BOUNDARY.  Today the
//     cycle's full call is the stage-1 solve itself: inside that step the RHS is the true
//     heating, in the others the linearised model, so every refresh puts a model jump of
//     one dt into the integral (measured: first order in dt even with the refreshes at
//     fixed simulated times, ck_cad_dt).  With the key, at the START of a due step
//     (ncycle % ck_impl_every == 0, or with ck_cad_dt a multiple of P at t) the model is
//     re-measured on U^n by a backward-Euler probe of ck_cad_probe_fac x dt (u0 = u1 for
//     it, then put back; Q0 = the probe's rate, D, T0 at U^n, model time t), and BOTH
//     stages of every step are the linearised step.  The ODE integrated is then exactly
//     the piecewise model, refreshed at step boundaries (hesdirk2 is second order on it).
inline bool ck_cad_probe = false;
// (8b) problem/ck_cad_cont = true (with ck_cad_probe and ck_cad_drift): A CONTINUOUS
//     MODEL.  A refresh that replaces the model by the measured rate makes the RHS jump at
//     every refresh; each jump starts a relaxation transient in the cells whose exchange
//     time is ~dt, which the step integrates with order reduction (measured: first order
//     in dt at a FIXED refresh schedule, error proportional to the jump).  With the key the
//     new model starts from the OLD model's rate at U^n (continuous in time) and ramps to
//     the measured rate extrapolated one interval ahead:
//       Q0 = M_old(t_n, T_n),  G = G_meas + (S_n - Q0)/P,  T0 = T_n,  D = D_new,
//     S_n the probe's rate at U^n, G_meas the drift of S between the last two refreshes,
//     P = ck_cad_dt (or ck_impl_every dt).  The model error stays O(P^2) and only kinks.
inline bool ck_cad_cont = false;
// problem/ck_cad_cont_len = f (default 1): the inherited offset S_n - M_old is removed
// linearly over f P instead of P (model rate S_n + G_meas tau + (M_old - S_n)
// max(0, 1 - tau/(f P))); continuous for any f > 0.
inline Real ck_cad_cont_len = 1.0;
inline Real ck_cad_cont_l = 1.0;     // L [s] of the last cont store
// (8d) ckdrift-1102 problem/ck_cad_drift_lim = f > 0 (default 0 = off, read only when
//     named): A LIMITER ON THE MEASURED DRIFT G (ck_cad_drift, both store paths).  The
//     drift is an explicit linear extrapolation of the rate.  In cells whose rate is
//     driven by a stiff coupling (the neighbours' response, a heating that grows with T,
//     which D <= 0 cannot hold) the measured rate alternates or grows between refreshes,
//     and extrapolating it amplifies the alternation by up to 3x per refresh (2 S_n -
//     S_n-1): measured on the fresh r1m (WASP-121b MHD) as a runaway to the EOS ceiling
//     1e6 K in a few columns at ~1e-4 bar by t ~1300-1900 s.  With the key a cell's new G
//     is set to 0 when it has the opposite sign of its previous G (an alternating rate is
//     not a drift), and |G| is capped so that the drift alone moves T by at most f T over
//     one refresh interval P (|G| P^2 <= f c_v T).  A smooth drift (same sign, small) is
//     unchanged, so the order is kept where the model is meaningful.
inline Real ck_cad_drift_lim = 0.0;
inline int64_t ck_cad_nlim_flip = 0, ck_cad_nlim_cap = 0;
// problem/ck_cad_drift_diag = true (TEST, read only when named): every store reports the
//     cell with the largest |G_meas| P^2/(c_v T) (location, T, rho, rates, sign flips).
inline bool ck_cad_drift_diag = false;
inline DvceArray5D<Real> *ck_cad_dg_ptr = nullptr;   // (m,3,k,j,i): x, flip, capped
// (9) ckdrift-1102 problem/ck_t2_pos_cap = f > 0 (default 0 = off): the hesdirk2 stage
//     increment inc = c1 K1 + c2 K2 is clamped to +f e_int on the heating side too
//     (the cooling side is clamped to -(1 - ck_t2_pos_th) e_int today).
inline Real ck_t2_pos_cap = 0.0;
inline DvceArray5D<Real> ck_cad_un;   // U^n during a probe store (cont)
inline Real ck_cad_probe_fac = 1.0e-4;
// problem/ck_cad_probe_dt = s > 0 [s]: the probe length in seconds instead of fac x dt.
// The probe moves every cell by ~pdt S, and the rate it measures includes the
// neighbours' response to that move (the diagonal correction to U^n cannot remove it):
// with pdt proportional to dt that is an O(dt) error in the model (first order).
inline Real ck_cad_probe_dt = 0.0;
// problem/ck_cad_probe_tolinc = tau (default 1e-6; 0 = the call's own test): the probe's
// Newton test is (1) of this file with tau, i.e. relative to the probe's own increment
// pdt S.  The production test |R| <= tol (e + eps e_max) is ~1e-8 e, larger than the
// whole probe increment of a slow cell (measured: garbage rates, dt -15 %).
inline Real ck_cad_probe_tolinc = 1.0e-6;
inline int64_t ck_cad_nprobe = 0;
inline bool ck_cad_probe_on = false;   // set around a probe call
// (8c) with ck_cad_probe and problem/ck_vef_dt = V: the VEF factors are refreshed by the
//     probe whose state U^n (time t) is the first at or after a point ck_cad_t0 + n V,
//     i.e. a point in (t - ck_impl_every dt, t] (stateless, so a restart keeps the
//     schedule).  With dt dividing V / ck_impl_every every ladder member refreshes at the
//     same instants; the factors' lag is then a function of time only.
// at the step start: is a probe refresh due (cycle or simulated-time schedule)
inline bool CkCadProbeDue(const Real t, const Real dt, const int64_t ncycle,
                          const int every) {
  if (ck_cad_dt > 0.0) {
    const Real eps = 1.0e-3*dt;
    const Real tt = t - ck_cad_t0;
    return std::floor((tt + dt - eps)/ck_cad_dt) > std::floor((tt - eps)/ck_cad_dt);
  }
  return (ncycle % static_cast<int64_t>(every) == 0);
}
// the store of a probe call of pdt made at model time tev (rt driver)
void CkCadProbeStore(Mesh *pm, const Real pdt, const Real tev);

// (4b) the VEF refresh test: ck_vef_dt > 0 replaces the cycle count by simulated time
inline bool CkVefDue2(Mesh *pm, const bool ckfop, const bool ckfst);

// read the keys (only when named); fused/glob are known by now
inline void Ck2ndInit(ParameterInput *pin, const bool fused, const int glob,
                      const int every) {
  std::string why;
  if (pin->DoesParameterExist("problem", "ck_impl_tol_inc")) {
    ck_tol_inc = pin->GetReal("problem", "ck_impl_tol_inc");
    if (pin->DoesParameterExist("problem", "ck_impl_tol_floor")) {
      ck_tol_floor = pin->GetReal("problem", "ck_impl_tol_floor");
    }
    if (pin->DoesParameterExist("problem", "ck_impl_tol_xfloor")) {
      ck_tol_xfloor = pin->GetReal("problem", "ck_impl_tol_xfloor");
    }
    if (ck_tol_inc > 0.0 && (!fused || glob != 0)) {
      why += " ck_impl_tol_inc needs ck_impl_fuse, ck_impl_debug <= 0, ck_impl_glob none;";
    }
  }
  if (pin->DoesParameterExist("problem", "ck_impl_xpass")) {
    ck_impl_xpass = pin->GetInteger("problem", "ck_impl_xpass");
    if (ck_impl_xpass > 0 && (!fused || glob != 0)) {
      why += " ck_impl_xpass needs ck_impl_fuse, ck_impl_debug <= 0, ck_impl_glob none;";
    }
  }
  if (pin->DoesParameterExist("problem", "ck_top_live")) {
    ck_top_live = pin->GetBoolean("problem", "ck_top_live");
  }
  if (pin->DoesParameterExist("problem", "ck_vef_dt")) {
    ck_vef_dt = pin->GetReal("problem", "ck_vef_dt");
  }
  if (pin->DoesParameterExist("problem", "ck_vef_tint")) {
    ck_vef_tint = pin->GetBoolean("problem", "ck_vef_tint");
  }
  if (pin->DoesParameterExist("problem", "ck_cad_drift")) {
    ck_cad_drift = pin->GetBoolean("problem", "ck_cad_drift");
    (void) every;   // (with ck_impl_every = 1 it acts under ck_t2_lin2 only)
  }
  if (pin->DoesParameterExist("problem", "ck_cad_probe")) {
    ck_cad_probe = pin->GetBoolean("problem", "ck_cad_probe");
    if (pin->DoesParameterExist("problem", "ck_cad_probe_tolinc")) {
      ck_cad_probe_tolinc = pin->GetReal("problem", "ck_cad_probe_tolinc");
    }
    if (ck_cad_probe && ck_cad_probe_tolinc > 0.0 && (!fused || glob != 0)) {
      why += " ck_cad_probe_tolinc needs ck_impl_fuse, ck_impl_debug <= 0, glob none;";
    }
    if (pin->DoesParameterExist("problem", "ck_cad_probe_dt")) {
      ck_cad_probe_dt = pin->GetReal("problem", "ck_cad_probe_dt");
    }
    if (pin->DoesParameterExist("problem", "ck_cad_probe_fac")) {
      ck_cad_probe_fac = pin->GetReal("problem", "ck_cad_probe_fac");
    }
    if (ck_cad_probe && every <= 1) why += " ck_cad_probe needs ck_impl_every > 1;";
    if (ck_cad_probe && (!pin->DoesParameterExist("problem", "ck_t2_start") ||
        pin->GetString("problem", "ck_t2_start").compare("probe") != 0)) {
      why += " ck_cad_probe needs ck_t2_start = probe;";
    }
  }
  if (pin->DoesParameterExist("problem", "ck_cad_lin_iter")) {
    ck_cad_lin_iter = pin->GetInteger("problem", "ck_cad_lin_iter");
  }
  if (pin->DoesParameterExist("problem", "ck_cad_cont")) {
    ck_cad_cont = pin->GetBoolean("problem", "ck_cad_cont");
    if (pin->DoesParameterExist("problem", "ck_cad_cont_len")) {
      ck_cad_cont_len = pin->GetReal("problem", "ck_cad_cont_len");
    }
    if (ck_cad_cont && (!ck_cad_probe || !ck_cad_drift)) {
      why += " ck_cad_cont needs ck_cad_probe and ck_cad_drift;";
    }
  }
  if (pin->DoesParameterExist("problem", "ck_cad_drift_lim")) {
    ck_cad_drift_lim = pin->GetReal("problem", "ck_cad_drift_lim");
    if (ck_cad_drift_lim > 0.0 && !ck_cad_drift) {
      why += " ck_cad_drift_lim needs ck_cad_drift;";
    }
  }
  if (pin->DoesParameterExist("problem", "ck_cad_drift_diag")) {
    ck_cad_drift_diag = pin->GetBoolean("problem", "ck_cad_drift_diag");
  }
  if (pin->DoesParameterExist("problem", "ck_t2_pos_cap")) {
    ck_t2_pos_cap = pin->GetReal("problem", "ck_t2_pos_cap");
  }
  if (pin->DoesParameterExist("problem", "ck_cad_drift2")) {
    ck_cad_drift2 = pin->GetBoolean("problem", "ck_cad_drift2");
    if (ck_cad_drift2 && !ck_cad_drift) why += " ck_cad_drift2 needs ck_cad_drift;";
  }
  if (pin->DoesParameterExist("problem", "ck_cad_tri")) {
    ck_cad_tri = pin->GetBoolean("problem", "ck_cad_tri");
    if (ck_cad_tri && !ck_cad_drift) why += " ck_cad_tri needs ck_cad_drift;";
  }
  if (pin->DoesParameterExist("problem", "ck_cad_dt")) {
    ck_cad_dt = pin->GetReal("problem", "ck_cad_dt");
    if (ck_cad_dt > 0.0 && every <= 1) why += " ck_cad_dt needs ck_impl_every > 1;";
  }
  if (pin->DoesParameterExist("problem", "ck_cad_t0")) {
    ck_cad_t0 = pin->GetReal("problem", "ck_cad_t0");
  }
  if (pin->DoesParameterExist("problem", "ck_xstep_dt")) {
    ck_xstep_dt = pin->GetReal("problem", "ck_xstep_dt");
    if (ck_xstep_dt > 0.0) {
      const Real r = ck_xstep_dt/ck_cad_dt;
      if (!(ck_cad_dt > 0.0) || std::fabs(r - std::round(r)) > 1.0e-9) {
        why += " ck_xstep_dt needs ck_cad_dt > 0 and a multiple of it;";
      }
      if (pin->GetOrAddInteger("problem", "ck_impl_xstep", 0) <= 0) {
        why += " ck_xstep_dt needs ck_impl_xstep > 0;";
      }
    }
  }
  if (ck_cad_dt > 0.0 && ck_vef_dt > 0.0) {
    const Real b = (ck_xstep_dt > 0.0) ? ck_xstep_dt : ck_cad_dt;
    const Real r = ck_vef_dt/b;
    if (std::fabs(r - std::round(r)) > 1.0e-9 || r < 0.5) {
      why += " ck_vef_dt must be a multiple of ck_xstep_dt (or ck_cad_dt);";
    }
  }
  if (!why.empty()) {
    std::cout << "### FATAL ERROR (ck2nd):" << why << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

}  // namespace two_stream_rt

#endif  // UTILS_CK_2ND_HPP_
