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
// (4) called by the pass after the refresh block: refreshed = this pass re-formed the
// factors, storing = this pass forms the stored operator (defined in the rt driver)
void CkVefTint(Mesh *pm, const bool refreshed, const bool storing);

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
  if (pin->DoesParameterExist("problem", "ck_cad_tri")) {
    ck_cad_tri = pin->GetBoolean("problem", "ck_cad_tri");
    if (ck_cad_tri && !ck_cad_drift) why += " ck_cad_tri needs ck_cad_drift;";
  }
  if (!why.empty()) {
    std::cout << "### FATAL ERROR (ck2nd):" << why << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

}  // namespace two_stream_rt

#endif  // UTILS_CK_2ND_HPP_
