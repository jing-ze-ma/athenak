#ifndef UTILS_TWO_STREAM_COLUMN_CK_HPP_
#define UTILS_TWO_STREAM_COLUMN_CK_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file two_stream_column_ck.hpp
//! \brief problem/ck_implicit: a BACKWARD-EULER column solve for the CORRELATED-K
//! two-stream (deep_hot_jupiter_rt).  Default OFF, in which case nothing in this file is
//! allocated, launched or read and the answer is bitwise today's.
//!
//! WHY.  What the gas receives on the correlated-k path today is NOT the flux divergence
//! the sweep computed.  The apply block relaxes each cell separately toward the
//! equilibrium of ONE sweep -- absorbed field frozen -- and then clips the step with
//! rt_de_max and the floors.  The existing rt_desum diagnostic measures the gap:
//! tests_ck_sph/README.md 0.3 reports -2.9 .. -5.9 % with ck_spherical on and
//! -36.8 .. -40.5 % with it off, on the production radial grid.  A flux-form scheme whose
//! face fluxes telescope exactly (that README 0.1) still loses that fraction of the
//! energy in transients, because the applied de is not the divergence of any flux.
//!
//! WHAT THIS DOES.  Per column, unknown = the cell internal energy e_i, the balance is
//!     R_i(e) = e_i - e*_i - bdt * S_i(e) = 0,
//! S_i the FULL correlated-k source the apply block already forms: the (1-w)-weighted
//! two-stream flux divergence, the tau-blend handover, and the direct beam.  R is
//! evaluated by RE-RUNNING THE ORDINARY SWEEP at the current iterate -- one ck sweep per
//! Newton pass -- so the fixed point is the exact backward-Euler balance whatever the
//! Jacobian gets wrong.  At convergence the gas is handed EXACTLY bdt*S(e^{n+1}) with no
//! relaxation factor, no rt_de_max clip and no per-cell damping: flux-form conservation
//! to the Newton tolerance.
//!
//! THE JACOBIAN.  At frozen opacity the thermal exchange is LINEAR in the band Planck
//! functions B_b(T), so dS_i/dT_j = sum_{bands,g} w_g (dSrc_i/dB_b,j)(dB_b/dT)_j.  The
//! NEAREST-NEIGHBOUR part of dSrc_i/dB is accumulated inside the ck sweep itself out of
//! the layer coefficients already in registers (see "ck_implicit" in the rt_chain_ck
//! kernel); everything two cells away or further, and the direct beam -- which at frozen
//! opacity is temperature-INDEPENDENT and therefore has no Jacobian term at all -- stay
//! in the residual.  Rows:
//!     J_{i,i-1} = -bdt jac0_i/cv_{i-1},  J_{ii} = 1 - bdt jac1_i/cv_i,
//!     J_{i,i+1} = -bdt jac2_i/cv_{i+1}
//! with cv = de/dT.  jac1 < 0 (the cell's own emission) and jac0, jac2 >= 0 (what the
//! neighbours send it), so J is an M-matrix for any bdt and the Thomas sweep is stable
//! without pivoting.
//!
//! WHAT IS NOT IN THE JACOBIAN, deliberately, and why it costs nothing at the fixed
//! point: (i) the direct beam, exactly temperature-independent here; (ii) the
//! ck_spherical
//! PROBE passes -- the two extra passes that build the face-mixing constants Cmx are
//! re-run every Newton pass, so they are exact in R, but their B-dependence is dropped
//! from J (the face-mixing FACTORS (1+beta), (1-beta) and the up ray's area frame change
//! ARE carried, since they are free); (iii) the tau-blend handover d[wF]/dr; (iv) the
//! opacity's own temperature dependence, which is frozen over the step by construction.
//! All four make this a quasi-Newton, not a Newton: the iteration converges linearly
//! rather than quadratically where they matter, and lands on the same root.
//!
//! WHAT IS NOT FOLDED IN.  The radial radiative CONDUCTION below/around ck_pcut_bar
//! (rad_tau_lo/hi) is NOT folded into this tridiagonal -- unlike the grey
//! rt_implicit_column = 3, which merges it into Conduction::ImplicitRadialUpdate.  It
//! stays the separate operator it is today, so the 10 bar hand-over remains split: the
//! blend's handover flux enters S_i explicitly (lagged by one operator), exactly as it
//! does with ck_implicit off.  Folding it in needs the conduction tridiagonal to accept a
//! correlated-k row, which is a second deliverable.
//!
//! HEAT CAPACITY.  cv = e/T, i.e. e ~ T -- the SAME approximation the semi-implicit
//! apply makes (see rt_semi_lin and the e_eq = e (A/E)^(1/4) form).  It is poor where H2
//! is partly dissociated, but it appears ONLY in the Jacobian: it changes how many passes
//! the Newton takes and not what it converges to.
//!
//! CADENCE.  This solves the balance over the interval bdt it is handed, which on this
//! pgen is one RK STAGE (the source term is called per stage; with rk2 that is two ck
//! sweeps per hydro step today).  ck_implicit does not change that: it makes each stage's
//! source implicit, at ck_impl_maxit sweeps per stage instead of one.

#include <math.h>
#include <cstdint>
#include <cstdio>

#include <iostream>

#include "athena.hpp"
#include "mesh/mesh.hpp"

namespace two_stream_rt {

//----------------------------------------------------------------------------------------
// problem/ck_implicit: the switch.  false = bitwise the semi-implicit per-cell apply.
inline bool ck_implicit = false;
// problem/ck_impl_tol: the RESIDUAL tolerance, max_i |R_i|/(e_i + eps e_max) over the
// column set, with eps = ck_impl_norm_eps.  Measuring a thin top cell against its own
// energy demands an absolute precision nothing else in the scheme has (the same reason
// rt_impl_norm exists on the grey path), so it is measured against its share instead.
inline Real ck_impl_tol = 1.0e-8;
inline Real ck_impl_norm_eps = 1.0e-3;
// problem/ck_impl_dtol: also stop when the Newton STEP is small, max_i |de_i|/e_i.  Both
// tests must pass.
inline Real ck_impl_dtol = 1.0e-8;
// problem/ck_impl_maxit: the pass cap.  Each pass is one full ck sweep.
inline int ck_impl_maxit = 8;
// problem/ck_impl_dtmax: cap on |de|/e per Newton pass.  A belt, not a limiter: it is
// applied to the STEP, so a capped pass is simply a shorter one and the fixed point is
// unchanged.  A column that needed it is counted.
inline Real ck_impl_dtmax = 0.5;
// problem/ck_impl_demax: bound on the TOTAL excursion from e^n over the step,
// |e - e*| <= ck_impl_demax e*.  This is the same bound rt_de_max (default 0.5) imposes
// on the semi-implicit apply today, and it exists for one measured reason: with
// problem/ck_spherical = FALSE the plane-parallel correlated-k source has NO
// positive-energy backward-Euler root in the thin top cells of the production grid (the
// undiluted A F manufactures a cooling rate ~13x the cell's own internal energy per
// stage that is almost independent of T, so the Newton runs the cell to e < 0).  The
// Newton then does not converge, and this bound keeps the answer no worse than today's
// while the fallback counter says so.  With ck_spherical = true it is never reached --
// measured, capped = 0 on every gate run.  0 = no bound.
inline Real ck_impl_demax = 0.5;
// problem/ck_impl_verbose: print the per-call pass count and residual.
inline bool ck_impl_verbose = false;
// problem/ck_impl_debug: print the worst cell of every pass with its whole row.
inline int ck_impl_debug = 0;
// problem/ck_impl_refresh_kappa: re-look-up the correlated-k opacity (and hence the beam
// transmission tau_ray) at every Newton pass instead of freezing it over the step.  OFF
// by default.  Cost: the pre-opacity kernel is a ck_continuum call plus two table index
// look-ups per cell and the beam's tau_ray is rebuilt inside the sweep anyway, so a
// refresh pass costs one extra rt_pre_opac launch, ~5 % of a sweep -- but it makes the
// Jacobian inconsistent with the residual (dkappa/dT is nowhere in J), so it converges
// more slowly.  Provided because "frozen opacity" is an assumption, not a law.
inline bool ck_impl_refresh_kappa = false;
// problem/ck_impl_tau_min: THE TWO-LEVEL SPLIT, ported from the grey mode-3 column's
// rt_impl_tau_min (deleted in 912ef43c; the design note is in that commit's parent at
// two_stream_rt.hpp:763).  Only a cell whose OWN optical depth kappa rho dr reaches this
// goes into the tridiagonal.  Linearising the emission about the current state gives a
// row whose diagonal is 1 - bdt dE/de, and as the cell's own emission E -> 0 that
// diagonal degenerates to 1: the Newton step is then the bare Picard step
// de = e* + bdt S - e, which at large bdt asks a thin top cell for many times its own
// energy, drives it through zero and diverges -- measured before this switch existed at
// 10x and 100x the production dt (README_phase2.md gate e).  In a cell with dtau >~ 1 the
// field is within a factor of its own B, the linearisation is excellent, and those are
// exactly the cells the semi-implicit (1-e^-x)/x damping mistreats.
//
// A thin cell is NOT dropped from the balance: it is solved by its own BRACKETED,
// positivity-preserving implicit step (CkThinSolve below) with the absorbed field frozen
// at the current iterate, and it enters the thick cells' rows only through that field.
// The fixed point of the two-level iteration is therefore still the exact backward-Euler
// balance, and the residual test measures thin and thick cells alike.  Set 0 to put every
// cell in the tridiagonal (the pre-split behaviour).
//
// MEASURED, and this is a DEVIATION from the grey template: on this grid the cell's own
// dtau is NOT the useful discriminant, because the only optical depth available before
// the sweep is the CONTINUUM one (the line opacity is a function of the g index, not of
// the cell) and it under-estimates dtau by enough that tau_min = 1 calls 63 % of the
// production column thin.  The iteration is then per-cell Picard almost everywhere and
// the residual after 20 passes at the production dt is 2.6e-03, against 5.7e-05 with no
// split at all.  So this gate is kept, and available, but its DEFAULT IS 0 (off); what
// does the work is ck_impl_arat below.  Both must pass for a cell to be called thick.
inline Real ck_impl_tau_min = 0.0;
// problem/ck_impl_arat: THE DISCRIMINANT THAT ACTUALLY WORKS.  A cell is thick only
// while the field it absorbs is within this factor of its own emission, A <= arat E.
// That is the ratio the linearisation's error is governed by, and the apply block's own
// note derives it: the linear step asks for de = e (A - E)/(4 E) where the true
// backward-Euler root is e ((A/E)^(1/4) - 1), so the overshoot is a function of A/E
// alone -- 1.3x at A/E = 2, 3.8x at 16, unbounded as E -> 0.  dtau is only ever a proxy
// for it (a thermalised cell has A/E = O(1)), and here it is a bad one.  A COOLING cell
// (A < E) is always safe -- the linearisation under-steps -- so only A > arat E is split
// out.  Evaluated from the stored S and E, i.e. exactly, with a one-pass lag: the apply
// of pass p classifies, rt_pre_opac of pass p+1 acts on it, so the decoupling of the
// matrix and the rows the solve builds always refer to the SAME classification.
inline Real ck_impl_arat = 2.0;
// problem/ck_impl_colskip: once a COLUMN's residual is below ck_impl_tol, stop sweeping
// it.  The columns of this pgen are independent -- the ck sweep, the beam and the apply
// are all per-column -- so a converged column contributes nothing further, and its
// arrays already hold the converged sweep.  The pass count is a max over columns, so a
// handful of stiff columns otherwise pay for all of them.
inline bool ck_impl_colskip = true;
// problem/ck_impl_once: run the whole correlated-k radiation ONCE per hydro step, with
// the full cycle dt, after the last RK stage (ProblemGenerator::user_split_func with
// user_split_once), instead of once per stage.  Halves the sweeps under rk2.  Enrolled by
// the pgen; refused unless ck_implicit is on.
inline bool ck_impl_once = false;
// problem/ck_impl_frozen_op: FREEZE THE EXCHANGE OPERATOR over the Newton passes.
// At frozen opacity the thermal two-stream is LINEAR in the band Planck functions: the
// sweep is the application of a fixed operator M to B, and everything in M -- the layer
// optical depths, the exponential coefficients, the BFace weights, the face mixing --
// depends on the opacity and the geometry alone.  Pass 0 therefore STORES what the sweep
// computed from the opacity (kappa rho per cell and chain, and -- see ck_impl_frozen_cof
// -- the three half-layer coefficients of `step`), and every later pass re-applies the
// operator to the new B_b(T) by the same recurrences with the table look-ups and the
// exponentials replaced by loads.  The DIRECT BEAM is exactly temperature-independent at
// frozen opacity, so its whole deposit Qb_g is frozen too and its pseudo-spherical ray
// integration (the O(N^2) chord walk of ck_beam_sph) is not re-run at all.
//
// The arithmetic of a frozen pass is the SAME arithmetic: the stored quantities are the
// bits pass 0 computed, the optical depths are re-formed from the stored kappa rho by the
// same expressions, so a frozen pass is BITWISE the pass it replaces and the converged
// state is bitwise the phase-2 one.  Refused (silently ignored) with
// ck_impl_refresh_kappa, which by construction wants the opacity rebuilt every pass.
inline bool ck_impl_frozen_op = false;
// problem/ck_impl_frozen_cof: store the three half-layer coefficients (e0, cin, cout) as
// well as kappa rho.  ON costs 4 Reals per (cell, chain) and removes every expm1 from a
// frozen pass; OFF costs 1 Real and pays one expm1 per half layer, keeping only the table
// look-ups.  THE MEMORY IS THE WHOLE TRADE: on the production mesh (6 x 32 x 32 columns,
// ~128 radial cells, 88 chains at nquad = 1) the four arrays are 69 M entries, i.e.
// 2.2 GB, against 0.55 GB for kappa rho alone.  The column solve is a CPU-side
// deliverable today; on a GPU the default should be reconsidered per mesh, and the
// natural third option -- storing per band block and looping the blocks -- is not built.
inline bool ck_impl_frozen_cof = true;
// problem/ck_impl_warm: start the Newton from the PREVIOUS call's converged increment.
// The fixed point does not move -- only the initial iterate does -- so a warm start can
// change the pass count and nothing else.  The seed is the total de the previous call
// applied, clipped to ck_impl_dtmax of the cell's energy; e^n is recorded net of it, so
// the balance being solved is unchanged.  NOT written to the restart file: after a
// restart (and on the first call of a run) the seed array is zero and the call is a cold
// start, which is correct, only slower for one step.
inline bool ck_impl_warm = false;

// the Newton pass index inside one RT call; read by the pass function to decide whether
// the opacity is rebuilt.  -1 = ck_implicit is off.
inline int ck_impl_pass = -1;

// dSrc_i/dT_{i-1,i,i+1}, (m, 3, k, j, i), summed over bands and g-points
inline DvceArray5D<Real> *ck_jac_ptr = nullptr;
// dB_b/dT at the cell, (m, CK_NB, i, k, j) -- the Bb_g layout exactly
inline DvceArray5D<Real> *ck_dbdt_ptr = nullptr;
// the full source S_i the apply block formed, and the internal energy it formed it at
inline DvceArray4D<Real> *ck_src_ptr = nullptr;
inline DvceArray4D<Real> *ck_ei_ptr = nullptr;
// e*_i, the internal energy the step started from
inline DvceArray4D<Real> *ck_estar_ptr = nullptr;
// the cell's OWN emission, with exactly the weights the semi-implicit apply gives it
// (the (1-w) tau-blend factor, zero below the band cut).  Needed by the thin-cell solve
// and by the fallback, which both have to separate S = A - E.
inline DvceArray4D<Real> *ck_em_ptr = nullptr;
// 1 = this cell is in the tridiagonal, 0 = it takes CkThinSolve.  ck_thk is written by
// the apply block of pass p; ck_thu is the copy rt_pre_opac of pass p+1 froze and acted
// on, and is the one the matrix and the solve of that pass both read, so the decoupling
// and the rows can never disagree.  Zero-initialised = every cell thin, which is the
// safe state the very first pass of a run starts from.
inline DvceArray4D<Real> *ck_thk_ptr = nullptr;
inline DvceArray4D<Real> *ck_thu_ptr = nullptr;
// per column, 1 = converged and no longer swept.  See ck_impl_colskip.
inline DvceArray3D<Real> *ck_done_ptr = nullptr;
// the Thomas sweep's two work rows, (m,k,j,i)
inline DvceArray4D<Real> *ck_cp_ptr = nullptr;
inline DvceArray4D<Real> *ck_dp_ptr = nullptr;
// ---- problem/ck_impl_frozen_op: the STORED OPERATOR, laid out (m, chain, i, k, j) like
// Bb_g and kc_g, so that the lanes of a wave -- which vary in j -- read adjacent Reals.
// kappa rho per cell and chain, and the three coefficients of that cell's HALF layer
// (every `step` in the sweep is the half layer of one cell, taken at that chain's mu, so
// one triple per (cell, chain) serves the two probe passes, the down-sweep and the
// up-sweep alike).
inline DvceArray5D<Real> *ck_kro_ptr = nullptr;
inline DvceArray5D<Real> *ck_c0_ptr = nullptr;
inline DvceArray5D<Real> *ck_ci_ptr = nullptr;
inline DvceArray5D<Real> *ck_co_ptr = nullptr;
// the top boundary layer's transmission factor (1 - e^-dtau) per column and chain, which
// multiplies the ghost cell's Planck function: the one piece of the operator that is not
// per cell
inline DvceArray4D<Real> *ck_tpf_ptr = nullptr;
// ---- problem/ck_impl_warm: the total increment this call has applied so far, carried
// over to seed the next one, and the seed actually applied (so that e^n can be recorded
// net of it)
inline DvceArray4D<Real> *ck_dep_ptr = nullptr;
inline DvceArray4D<Real> *ck_seed_ptr = nullptr;
// 0 = max |R|/(e + eps e_max), 1 = max |de|/e, 2 = capped-step count, 3 = fallback
// count, 4 = sum (e - e*) dx, 5 = sum bdt S dx, 6 = columns still active, 7 = thin
// cells.  Slots 4 and 5 are the DIRECT analogue
// of the rt_desum line -- what the gas received against what the sweep asked for -- and
// their ratio is the number tests_ck_sph/README.md 0.3 reports for the semi-implicit
// apply.  Here it is 1 to the Newton tolerance by construction.
inline DvceArray1D<Real> *ck_conv_ptr = nullptr;
inline Real ck_impl_last_gap = 0.0;
// host-side record of the last call, for the report line
inline int ck_impl_last_it = 0;
inline Real ck_impl_last_res = 0.0;
inline Real ck_impl_last_dst = 0.0;
inline int ck_impl_ncap = 0;
inline int ck_impl_nfall = 0;
inline int ck_impl_nonconv = 0;
inline int ck_impl_nthin = 0;
inline int ck_impl_nactive = 0;
// cumulative sweep count over the run, for the cost report
inline std::int64_t ck_impl_nsweep = 0;

//----------------------------------------------------------------------------------------
//! \fn Real CkThinSolve
//! \brief the OPTICALLY THIN cell's own implicit step.  Solves, exactly and by
//! bisection on a guaranteed bracket,
//!     x - e* - bdt (A - Em (x/e)^4) = 0,   A = S + Em frozen at the current iterate,
//! i.e. the backward-Euler balance of the cell with the absorbed field lagged (the outer
//! Newton pass is what un-lags it).  Em ~ T^4 and, under the e ~ T the whole scheme
//! assumes for c_v, E(x) = Em (x/e)^4.  f is strictly increasing with
//! f(0) = -(e* + bdt A)
//! < 0 and f(max(e*+bdt A, e)) > 0, so the root is unique, bracketed, and POSITIVE: the
//! cell can never be driven through zero, which is the whole failure mode at large bdt.
//! Where bdt is small the root is e* + bdt S to O(bdt^2); where it is large the root is
//! radiative equilibrium E = A.  Both limits are exact, unlike the semi-implicit
//! (1-e^-x)/x step, and unlike that step the fixed point over the Newton passes IS the
//! backward-Euler balance, so the residual test converges on thin cells too.

KOKKOS_INLINE_FUNCTION
Real CkThinSolve(const Real ei, const Real est, const Real src, const Real em,
                 const Real bdt) {
  const Real absn = src + em;
  const Real r0 = est + bdt*absn;
  if (!(em > 0.0)) {
    // no emission to make implicit: the balance is linear and exact
    return (r0 > 0.0) ? r0 : 0.5*ei;
  }
  if (!(r0 > 0.0)) return 0.5*ei;   // no positive root at all: halve, retry next pass
  Real lo = 0.0;
  Real hi = (r0 > ei) ? r0 : ei;
  for (int q=0; q<60; ++q) {
    const Real x = 0.5*(lo + hi);
    const Real u = x/ei;
    const Real f = x - r0 + bdt*em*(u*u)*(u*u);
    if (f > 0.0) {
      hi = x;
    } else {
      lo = x;
    }
    if (hi - lo <= 1.0e-15*hi) break;
  }
  return 0.5*(lo + hi);
}

//----------------------------------------------------------------------------------------
//! \fn void CkImplAlloc
//! \brief allocate the column-solve scratch on the first call.  Never reached with
//! ck_implicit off.

inline void CkImplAlloc(const int nmb, const int nb, const int nch, const int n1,
                        const int n2, const int n3) {
  if (ck_jac_ptr != nullptr &&
      ck_jac_ptr->extent(0) == static_cast<size_t>(nmb) &&
      ck_jac_ptr->extent(4) == static_cast<size_t>(n1)) {
    return;
  }
  if (ck_jac_ptr != nullptr) {
    delete ck_jac_ptr;
    delete ck_dbdt_ptr;
    delete ck_src_ptr;
    delete ck_ei_ptr;
    delete ck_estar_ptr;
    delete ck_cp_ptr;
    delete ck_dp_ptr;
    delete ck_em_ptr;
    delete ck_thk_ptr;
    delete ck_thu_ptr;
    delete ck_done_ptr;
    delete ck_dep_ptr;
    delete ck_seed_ptr;
    if (ck_kro_ptr != nullptr) {
      delete ck_kro_ptr;
      delete ck_c0_ptr;
      delete ck_ci_ptr;
      delete ck_co_ptr;
      delete ck_tpf_ptr;
      ck_kro_ptr = nullptr;
    }
  }
  ck_jac_ptr = new DvceArray5D<Real>("ck_jac", nmb, 3, n3, n2, n1);
  ck_dbdt_ptr = new DvceArray5D<Real>("ck_dbdt", nmb, nb, n1, n3, n2);
  ck_src_ptr = new DvceArray4D<Real>("ck_src", nmb, n3, n2, n1);
  ck_ei_ptr = new DvceArray4D<Real>("ck_ei", nmb, n3, n2, n1);
  ck_estar_ptr = new DvceArray4D<Real>("ck_estar", nmb, n3, n2, n1);
  ck_cp_ptr = new DvceArray4D<Real>("ck_cp", nmb, n3, n2, n1);
  ck_dp_ptr = new DvceArray4D<Real>("ck_dp", nmb, n3, n2, n1);
  ck_em_ptr = new DvceArray4D<Real>("ck_em", nmb, n3, n2, n1);
  ck_thk_ptr = new DvceArray4D<Real>("ck_thk", nmb, n3, n2, n1);
  ck_thu_ptr = new DvceArray4D<Real>("ck_thu", nmb, n3, n2, n1);
  ck_done_ptr = new DvceArray3D<Real>("ck_done", nmb, n3, n2);
  ck_dep_ptr = new DvceArray4D<Real>("ck_dep", nmb, n3, n2, n1);
  ck_seed_ptr = new DvceArray4D<Real>("ck_seed", nmb, n3, n2, n1);
  // ---- problem/ck_impl_frozen_op: the stored operator.  The coefficient triple is the
  // memory-hungry half and is allocated only when ck_impl_frozen_cof asks for it; the
  // arrays are 1-element dummies otherwise, captured and never read.
  if (ck_impl_frozen_op && !ck_impl_refresh_kappa) {
    const int nc3 = ck_impl_frozen_cof ? nch : 1;
    const int n13 = ck_impl_frozen_cof ? n1 : 1;
    const int n33 = ck_impl_frozen_cof ? n3 : 1;
    const int n23 = ck_impl_frozen_cof ? n2 : 1;
    ck_kro_ptr = new DvceArray5D<Real>("ck_kro", nmb, nch, n1, n3, n2);
    ck_c0_ptr = new DvceArray5D<Real>("ck_c0", nmb, nc3, n13, n33, n23);
    ck_ci_ptr = new DvceArray5D<Real>("ck_ci", nmb, nc3, n13, n33, n23);
    ck_co_ptr = new DvceArray5D<Real>("ck_co", nmb, nc3, n13, n33, n23);
    ck_tpf_ptr = new DvceArray4D<Real>("ck_tpf", nmb, nch, n3, n2);
  }
  if (ck_conv_ptr == nullptr) {
    ck_conv_ptr = new DvceArray1D<Real>("ck_conv", 8);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CkWarmSeed
//! \brief problem/ck_impl_warm.  Called ONCE at the top of every RT call, before the
//! first sweep, and it does two things: it RESETS the per-call increment accumulator, and
//! -- when the switch is on and the arrays carry a previous call -- it seeds the gas with
//! the increment that call converged on, clipped to ck_impl_dtmax of the cell's energy.
//!
//! The seed is recorded in ck_seed so that the apply block of pass 0 can record e^n NET
//! of it: the balance solved is then exactly the one the cold start solves, and the only
//! thing that changed is where the Newton starts.  With the switch off the seed is 0 and
//! the array is a no-op the apply subtracts, i.e. bitwise today's.
//!
//! RESTART.  Nothing here is written to the restart file; a fresh run and a restart both
//! start from zeros, which is a cold start and is correct.

inline void CkWarmSeed(Mesh *pm, DvceArray5D<Real> u0) {
  if (ck_dep_ptr == nullptr) return;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pm->pmb_pack->nmb_thispack - 1;
  auto dep_ = *ck_dep_ptr;
  auto sd_ = *ck_seed_ptr;
  auto ei_ = *ck_ei_ptr;
  const bool warm = ck_impl_warm;
  const Real cap = ck_impl_dtmax;
  par_for("ck_warm_seed", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real s = 0.0;
    if (warm) {
      const Real e = ei_(m,k,j,i);          // the previous call's last iterate
      s = dep_(m,k,j,i);
      if (e > 0.0) {
        const Real lim = cap*e;
        if (s > lim) s = lim;
        if (s < -lim) s = -lim;
      } else {
        s = 0.0;
      }
      if (!(s == s)) s = 0.0;
      u0(m,IEN,k,j,i) += s;
    }
    sd_(m,k,j,i) = s;
    dep_(m,k,j,i) = s;     // the accumulator starts from what has already been applied
  });
}

//----------------------------------------------------------------------------------------
//! \fn int CkImplStep
//! \brief form the residual of the CURRENT iterate, test it, and -- if it is not yet
//! converged -- solve the tridiagonal and apply the step to u0.  Returns 1 if a step was
//! taken, 0 if the iterate was already converged (in which case NOTHING is applied and
//! the state handed to the gas satisfies |R| <= tol).
//!
//! The Thomas sweep runs over the cells the correlated-k domain owns, icut .. ie.  Cells
//! below the cut receive no radiative source at all and are not in the system.
//!
//! FALLBACK.  A row whose diagonal is not positive, or a solve that produces a
//! non-finite step, drops that CELL back to the bounded relaxation the semi-implicit
//! apply would have taken, de = deq (1 - exp(-x)) with deq = e((A/E)^(1/4) - 1) -- here
//! with A and E built from the same source, so it degenerates to the explicit step where
//! the emission is unavailable.  It is counted in slot 3.

inline int CkImplStep(Mesh *pm, DvceArray5D<Real> u0, DvceArray3D<int> icut_,
                      DvceArray4D<Real> T_, DvceArray4D<Real> dx1_, const Real bdt) {
  auto &indcs = pm->mb_indcs;
  const int ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pm->pmb_pack->nmb_thispack - 1;
  auto jac_ = *ck_jac_ptr;
  auto src_ = *ck_src_ptr;
  auto ei_ = *ck_ei_ptr;
  auto est_ = *ck_estar_ptr;
  auto cnv_ = *ck_conv_ptr;
  auto cp_ = *ck_cp_ptr;
  auto dp_ = *ck_dp_ptr;
  auto em_ = *ck_em_ptr;
  auto thk_ = *ck_thu_ptr;
  auto done_ = *ck_done_ptr;
  // problem/ck_impl_warm: the running total this call has applied, carried to the next
  auto dep_ = *ck_dep_ptr;
  const bool skipcol = ck_impl_colskip;
  Kokkos::deep_copy(cnv_, 0.0);
  const Real tol = ck_impl_tol;
  const Real dtol = ck_impl_dtol;
  const Real eps = ck_impl_norm_eps;
  const Real dcap = ck_impl_dtmax;
  const Real detot = ck_impl_demax;

  // ---- pass 1: the residual norm of the CURRENT iterate --------------------------
  par_for("ck_impl_res", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    const int ic = icut_(m,k,j);
    Real emax = 0.0;
    for (int i=ic; i<ie+1; ++i) {
      const Real e = ei_(m,k,j,i);
      if (e > emax) emax = e;
    }
    if (!(emax > 0.0)) {
      done_(m,k,j) = 1.0;      // nothing usable here; it is not part of the system
      return;
    }
    Real rn = 0.0;
    Real sg = 0.0, ss = 0.0;
    for (int i=ic; i<ie+1; ++i) {
      const Real r = ei_(m,k,j,i) - est_(m,k,j,i) - bdt*src_(m,k,j,i);
      const Real s = fabs(r)/(ei_(m,k,j,i) + eps*emax);
      if (s > rn) rn = s;
      const Real dx = dx1_(m,k,j,i);
      sg += (ei_(m,k,j,i) - est_(m,k,j,i))*dx;
      ss += bdt*src_(m,k,j,i)*dx;
    }
    Kokkos::atomic_max(&cnv_(0), rn);
    Kokkos::atomic_add(&cnv_(4), sg);
    Kokkos::atomic_add(&cnv_(5), ss);
    // ---- ck_impl_colskip: this COLUMN's own verdict.  Once set it is never cleared:
    // the sweep no longer touches the column, so neither its source nor its energy can
    // change and the residual it was retired on is the one it keeps.
    if (rn <= tol) {
      done_(m,k,j) = 1.0;
    } else {
      Kokkos::atomic_add(&cnv_(6), 1.0);
    }
  });
  auto hc = Kokkos::create_mirror_view(cnv_);
  Kokkos::deep_copy(hc, cnv_);
  ck_impl_last_res = hc(0);
  ck_impl_last_gap = (hc(5) != 0.0) ? (hc(4)/hc(5) - 1.0) : 0.0;
  ck_impl_nactive = static_cast<int>(hc(6));
  // ---- problem/ck_impl_debug: print the WORST cell of the pass, with everything the
  // row is built from.  Diagnosis only; nothing reads it.
  if (ck_impl_debug > 0) {
    const Real rmax = hc(0);
    const Real emref = ck_impl_norm_eps;
    par_for("ck_impl_dbg", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      const int ic = icut_(m,k,j);
      Real emax = 0.0;
      for (int i=ic; i<ie+1; ++i) {
        if (ei_(m,k,j,i) > emax) emax = ei_(m,k,j,i);
      }
      if (!(emax > 0.0)) return;
      for (int i=ic; i<ie+1; ++i) {
        const Real r = ei_(m,k,j,i) - est_(m,k,j,i) - bdt*src_(m,k,j,i);
        const Real s = fabs(r)/(ei_(m,k,j,i) + emref*emax);
        if (s >= 0.999*rmax) {
          const Real emc = em_(m,k,j,i);
          Kokkos::printf("### ckdbg m=%d k=%d j=%d i=%d r=%.4e e=%.5e est=%.5e "
                         "src=%.5e em=%.5e AoverE=%.4e thick=%d T=%.4e bdt=%.3e\n",
                         m, k, j, i, s, ei_(m,k,j,i), est_(m,k,j,i), src_(m,k,j,i),
                         emc, (emc > 0.0) ? (src_(m,k,j,i) + emc)/emc : -1.0,
                         (thk_(m,k,j,i) > 0.0) ? 1 : 0, T_(m,k,j,i), bdt);
          return;
        }
      }
    });
  }
  if (hc(6) == 0.0) return 0;

  // ---- pass 2: assemble, solve, apply --------------------------------------------
  par_for("ck_impl_tri", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    const int ic = icut_(m,k,j);
    if (ic > ie) return;
    if (done_(m,k,j) > 0.0) return;      // converged; see ck_impl_colskip
    // the Thomas sweep.  One thread owns the column, exactly as the sweeps do; the two
    // work rows live in device arrays rather than on the stack so the kernel carries no
    // per-thread footprint that would keep it off an accelerator.
    const int n = ie - ic + 1;
    if (n < 1) return;
    bool bad = false;
    for (int q=0; q<n; ++q) {
      const int i = ic + q;
      const Real ei = ei_(m,k,j,i);
      const Real Ti = T_(m,k,j,i);
      if (!(ei > 0.0) || !(Ti > 0.0)) {
        bad = true;
        break;
      }
      const Real cvi = ei/Ti;
      const Real eim = (q > 0) ? ei_(m,k,j,i-1) : 0.0;
      const Real Tim = (q > 0) ? T_(m,k,j,i-1) : 1.0;
      const Real eip = (q < n-1) ? ei_(m,k,j,i+1) : 0.0;
      const Real Tip = (q < n-1) ? T_(m,k,j,i+1) : 1.0;
      const Real cvm = (q > 0 && Tim > 0.0 && eim > 0.0) ? (eim/Tim) : 1.0;
      const Real cvp = (q < n-1 && Tip > 0.0 && eip > 0.0) ? (eip/Tip) : 1.0;
      // ---- THE TWO-LEVEL SPLIT, see ck_impl_tau_min ------------------------------
      // A THIN row is an identity row whose right-hand side is the cell's own
      // bracketed implicit step.  It is decoupled from the tridiagonal in BOTH
      // directions: from the thin cell outwards because its row is the identity, and
      // from a thick neighbour inwards because the sweep zeroed its dB/dT, which is
      // the factor every off-diagonal entry pointing at it carries (see the jck block
      // in rt_chain_ck).  So the matrix stays the same M-matrix, one size smaller.
      const bool thin = !(thk_(m,k,j,i) > 0.0);
      Real a, b, c, d;
      if (thin) {
        a = 0.0;
        b = 1.0;
        c = 0.0;
        d = CkThinSolve(ei, est_(m,k,j,i), src_(m,k,j,i), em_(m,k,j,i), bdt) - ei;
        Kokkos::atomic_add(&cnv_(7), 1.0);
      } else {
        a = (q > 0) ? (-bdt*jac_(m,0,k,j,i)/cvm) : 0.0;
        b = 1.0 - bdt*jac_(m,1,k,j,i)/cvi;
        c = (q < n-1) ? (-bdt*jac_(m,2,k,j,i)/cvp) : 0.0;
        d = -(ei - est_(m,k,j,i) - bdt*src_(m,k,j,i));
      }
      if (!(b > 0.0)) {
        bad = true;
        break;
      }
      const Real den = b - a*((q > 0) ? cp_(m,k,j,i-1) : 0.0);
      if (!(fabs(den) > 0.0)) {
        bad = true;
        break;
      }
      cp_(m,k,j,i) = c/den;
      dp_(m,k,j,i) = (d - a*((q > 0) ? dp_(m,k,j,i-1) : 0.0))/den;
    }
    if (bad) {
      Kokkos::atomic_add(&cnv_(3), 1.0);
      // the bounded fallback for the whole column: every cell takes the per-cell
      // bracketed implicit step of the thin branch, which is positivity-preserving and
      // has the same fixed point as the tridiagonal, only with the neighbours lagged.
      for (int i=ic; i<ie+1; ++i) {
        const Real ei = ei_(m,k,j,i);
        if (!(ei > 0.0)) continue;
        Real de = CkThinSolve(ei, est_(m,k,j,i), src_(m,k,j,i), em_(m,k,j,i), bdt) - ei;
        const Real lim = dcap*ei;
        if (de > lim) de = lim;
        if (de < -lim) de = -lim;
        u0(m,IEN,k,j,i) += de;
        dep_(m,k,j,i) += de;
      }
      return;
    }
    Real dn = 0.0;
    Real prev = 0.0;
    bool cap = false;
    for (int q=n-1; q>=0; --q) {
      const int i = ic + q;
      Real de = dp_(m,k,j,i) - cp_(m,k,j,i)*prev;
      // NaN guard; <cmath> is not available on the device side
      if (!(de == de)) de = 0.0;
      prev = de;
      const Real ei = ei_(m,k,j,i);
      const Real lim = dcap*ei;
      if (de > lim) {
        de = lim;
        cap = true;
      }
      if (de < -lim) {
        de = -lim;
        cap = true;
      }
      // the total-excursion bound, see ck_impl_demax
      const Real es = est_(m,k,j,i);
      if (detot > 0.0 && es > 0.0) {
        Real en = ei + de;
        const Real ehi = es*(1.0 + detot);
        const Real elo = es*(1.0 - detot);
        if (en > ehi) {
          en = ehi;
          cap = true;
        }
        if (en < elo) {
          en = elo;
          cap = true;
        }
        de = en - ei;
      }
      u0(m,IEN,k,j,i) += de;
      dep_(m,k,j,i) += de;
      const Real s = (ei > 0.0) ? fabs(de)/ei : 0.0;
      if (s > dn) dn = s;
    }
    Kokkos::atomic_max(&cnv_(1), dn);
    if (cap) Kokkos::atomic_add(&cnv_(2), 1.0);
  });
  Kokkos::deep_copy(hc, cnv_);
  ck_impl_last_dst = hc(1);
  ck_impl_ncap = static_cast<int>(hc(2));
  ck_impl_nfall = static_cast<int>(hc(3));
  ck_impl_nthin = static_cast<int>(hc(7));
  // the step test.  The residual of the NEW iterate is measured at the top of the next
  // pass, which is the only honest place to measure it.
  if (hc(1) <= dtol && hc(0) <= tol) return 0;
  return 1;
}

}  // namespace two_stream_rt

#endif  // UTILS_TWO_STREAM_COLUMN_CK_HPP_
