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
#include <string>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"

namespace two_stream_rt {

//----------------------------------------------------------------------------------------
// problem/ck_implicit: the switch.  false = bitwise the semi-implicit per-cell apply.
// The variable starts false; deep_hot_jupiter_rt defaults the KEY to true on fresh rt_ck
// runs since defaults2-0927 (false on restarts; see the pgen).
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
// A DEAD cell: its conserved state gives no positive internal energy, so the RT read
// (eiN in two_stream_rt.hpp) clamped e to e(rho, tfloor) or, where even that is not
// positive, to 1e-300.  Such a cell cannot take a Newton step -- its "heat capacity"
// e/T is ~1e-302, so any residual asks for a temperature change of ~1e300 K -- and
// before this guard it was coupled into its neighbours' rows through that capacity,
// which stalled the column (measured: the wp transient, nquad = 2, a thick cell next to
// a dead one kept its residual at 1.98e-6 for 30 passes).  It is given an identity row
// with no step, and its neighbours' rows drop their coupling to it.  Healthy columns
// never reach the threshold, so they are unchanged bit for bit.
constexpr Real CkDeadE = 1.0e-250;
// problem/ck_impl_maxit: the pass cap.  Each pass is one full ck sweep.
inline int ck_impl_maxit = 24;  // default changed 09-30 (user), was 8
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
// problem/ck_impl_floorbound: THE NEWTON IS A BOUND-CONSTRAINED SOLVE, e >= e_floor.
// false = the unconstrained Newton, bit for bit.  deep_hot_jupiter_rt defaults the key
// (and ck_impl_kkt_demax) to true on fresh fused ck_implicit runs since defaults2-0927.
//
// e_floor = max(e(rho, tfloor), e(rho, pfloor)) is the state ConsToPrim restores a cell
// to anyway (problem/rt_floor_consistent applies the same rule to the semi-implicit
// apply).  In the night-side top slab of a fresh WASP-121b start (SMOKE_DIAG.md, 09-25)
// the backward-Euler root of the coldest top cells lies below that state: the Newton
// drives them down, the caps clip the step, the column never converges, and the floors
// put the energy straight back (eos_efloor ~1e7 events per 400 cycles).  With this on:
//   (a) a step that would take a cell below e_floor stops AT e_floor (a cell already at
//       or below it is not cooled further; heating steps are untouched);
//   (b) a cell ON the bound (e <= e_floor (1 + 1e-10)) whose residual points through it
//       (R = e - e* - h S > 0, i.e. the Newton would cool it further) is at its KKT
//       point: its residual is not counted in the convergence test.
// Away from the floors both are the identity, so the arithmetic of every other cell is
// unchanged; but a cell the old code pushed below e_floor now stops at it, so results
// change wherever the bound is active.  Counted in slot 16 (cells stopped at the floor)
// and slot 17 (cells excluded from the test as KKT).  Fused, ck_impl_glob = 0 and
// ck_impl_aa = 0 path only (refused otherwise).
inline bool ck_impl_floorbound = false;
// problem/ck_impl_kkt_demax: the same KKT rule (b) for the TOTAL-excursion bound
// ck_impl_demax: a cell sitting on e* (1 -+ demax) whose residual points through that
// bound is converged FOR THE CONSTRAINED PROBLEM.  The bound itself is unchanged (it is
// the limiter that keeps a dt collapse away; removing it, demax = 0 with dtmax 3,
// collapsed dt at cycle 22 in smokediag/dtm3d0); only the pass count and the
// NOT-CONVERGED flag change: the passes a capped cell cannot use are no longer spent.
// false = bitwise (default true on fresh fused ck_implicit dhj runs since
// defaults2-0927).  Counted in slot 17 with (b).
inline bool ck_impl_kkt_demax = false;
// problem/ck_impl_stalldbg: diagnosis of columns that stall at a fixed residual.  0 =
// off (bitwise).  N > 0: on every fused pass with index >= N, each column still active
// prints its worst cell, that cell's neighbours and every cell whose step was clipped
// (dtmax / demax / floor), with the unclipped Newton step, the applied step, the KKT
// flag and the magnitudes that set the round-off floor of the residual.  Nothing reads
// it.
inline int ck_impl_stalldbg = 0;
// problem/ck_impl_negdbg = N > 0 (diagnostic): on every fused pass print up to N Newton
// rows whose diagonal b = 1 - h J_ii/cv is not positive (they send the whole column to
// the bounded per-cell fallback), with the three Jacobian entries.  Print only.
inline int ck_impl_negdbg = 0;
// problem/ck_impl_rowfb (fused step only, default false = bitwise off): a Newton row
// whose diagonal b = 1 - h J_ii/cv is not positive (a cell whose net heating grows with
// its own temperature faster than cv/h: the frozen ck_sph_face = 5 operator has such
// cells in a strongly driven top) takes the bounded per-cell step of a thin cell; the
// rest of the column keeps its Newton step.  Off, such a row sends the WHOLE column to
// the per-cell fallback.  The residual and the tolerance are unchanged.
inline bool ck_impl_rowfb = false;
// problem/ck_impl_negpiv (fused step only, default false = bitwise off): keep the exact
// Newton row of a cell with a non-positive diagonal (b <= 0: its net heating grows
// with its temperature faster than cv/h) instead of treating it as bad; only a pivot
// |b| <= 1e-12 sends the column to the fallback.  Takes precedence over ck_impl_rowfb.
inline bool ck_impl_negpiv = false;
// problem/ck_impl_jnet (default false = bitwise off): keep a NEGATIVE net off-diagonal
// of the ck_impl_jac_lin tridiagonal (ck_jlin_sum) instead of dropping it; with
// ck_impl_jneg (the per-chain parts) the rows are then the full derivative.
inline bool ck_impl_jnet = false;
// problem/ck_impl_tol_ptop = p [bar] > 0 (default 0 = bitwise off): in the residual test
// of the fused step a cell with pressure below p is measured against e + eps_top e_max
// (eps_top = problem/ck_impl_norm_eps_top) instead of e + ck_impl_norm_eps e_max, i.e.
// a looser tolerance in the top, where the accuracy does not matter (user rule: p <
// 1e-6 bar); the cells below p keep the full test.  ck_pbar_ptr is the RT's pressure.
inline Real ck_impl_tol_ptop = 0.0;
inline Real ck_impl_norm_eps_top = 0.1;
inline DvceArray4D<Real> *ck_pbar_ptr = nullptr;
// problem/ck_sph_face = 5 (two_stream_rt.hpp): the cycle of the last refresh of the VEF
// factors (-1 = none yet) and the restart's rebuild of that refresh
// (two_stream_ck_rst_state.hpp: the factors are re-formed from the state the refresh
// read, so a restart is bitwise); ck_vef_on = (ck_sph_face == 5), set by the pgen
inline int ck_vef_last = -1;
inline bool ck_vef_force = false;
inline bool ck_vef_on = false;
// problem/ck_impl_jfd (diagnostic, > 0 on): on every pass that builds the tridiagonal
// from the factorisation, check its diagonal against a finite difference of the linear
// re-apply (two_stream_rt.hpp).  Print only; costs n1 extra linear passes.
inline int ck_impl_jfd = 0;
// problem/ck_impl_kkt_row: the ACTIVE-SET row for the KKT cells of ck_impl_floorbound /
// ck_impl_kkt_demax.  A cell on a bound whose residual points through it cannot move,
// but its row still asks the tridiagonal for a step, and that unrealised step feeds its
// neighbours' solutions through the off-diagonals: the neighbours then converge to the
// fixed point of the clipped iteration, not to R = 0, and stall there.  With this on the
// KKT cell gets the identity row with zero right-hand side (a = c = 0, b = 1, d = 0), so
// its neighbours solve the reduced system with that cell held.  The key defaults to
// TRUE since default-flips (deep_hot_jupiter_rt reads it; false = the old rows,
// bitwise).  Fused, glob = 0 path only.
inline bool ck_impl_kkt_row = false;
// problem/ck_impl_conserve: 0 = off (default, bitwise).  1 = the face fluxes a call
// reports are made consistent with the energy it applied in the cells it leaves on a
// bound (ck_impl_floorbound / ck_impl_kkt_demax KKT cells); 2 = in every cell.  Only the
// reported arrays change, never the solution.  See CkConserveFlux (two_stream_rt.hpp).
inline int ck_impl_conserve = 1;  // default changed 09-30 (user), was 0
inline DvceArray1D<Real> *ck_cons_stat_ptr = nullptr;
inline DvceArray5D<Real> *ck_cons_buf_ptr = nullptr;   // (m, 2, k, j, i) scratch
inline int64_t ck_cons_ncell = 0;          // cells corrected, last call (verbose only)
inline Real ck_cons_de = 0.0;              // their sum r V [erg], last call, this rank
// problem/ck_impl_osc = N > 0: PER-CELL AITKEN DAMPING of oscillating cells, from pass N
// of a call on.  0 = off (default, bitwise).  Measured on the 10x WASP-121b production
// (cknewton_0928/RESULTS.md): the calls that do not converge even with 24 or 60 passes
// sit in a period-2/4 limit cycle of the top cells (p ~ 1e-10 bar, T swinging 1700 <->
// 5300 K): the band Planck emission there rises far faster than its linearisation, so
// the Newton step of the hot iterate and of the cold one overshoot each other, and the
// neighbours, heated and cooled by that cell, swing with it at the dtmax cap.  With this
// on, a thick cell whose Newton step x reverses the sign of the step it applied on the
// previous pass, s, is a cell with a negative local multiplier lambda = x/s of the
// iteration, and it takes the step that removes that mode, w x with w = 1/(1 - lambda)
// = s/(s - x) in (0, 1) (Aitken's delta^2 on the cell).  A monotone (same-sign) step is
// untouched, so the iterate of a column that converges monotonically is bitwise; the
// residual, and therefore the fixed point, is untouched everywhere.  Counted in slot 18.
// Fused, glob = none, aa = 0 path only.  NOT NEEDED with the recommended 10x keys: a
// per-pass cap ck_impl_dtmax = 0.25 with ck_impl_maxit = 24 removes the limit cycles by
// itself (0/1500 NOT-CONVERGED over 6000 cycles, osc on or off); osc 12 is what makes
// dtmax 0.5 + maxit 24 converge (0/500 over 2000 cycles, 5/1500 over 6000).
inline int ck_impl_osc = 0;
inline DvceArray4D<Real> *ck_osc_ptr = nullptr;   // the step each cell applied last pass
inline int ck_impl_nosc = 0;                       // cells damped, last pass
inline int ck_impl_nfloor = 0;             // cells stopped at e_floor, last pass
inline int ck_impl_nkkt = 0;               // cells excluded from the test as KKT
// problem/ck_impl_verbose: print the per-call pass count and residual.
inline bool ck_impl_verbose = false;
// problem/ck_impl_debug: print the worst cell of every pass with its whole row.
inline int ck_impl_debug = 0;
// problem/ck_impl_ncloc = N > 0: on a NOT-CONVERGED call, print the location of up to
// N still-active columns (gid, m, k, j), their worst cell i with its residual, pressure
// and T.  Print only (works under ck_impl_every); 0 (default) = off, bitwise.
inline int ck_impl_ncloc = 0;
// problem/ck_impl_ncref = N > 0 (diagnostic, for the accuracy of a NOT-CONVERGED call):
// after pass N of a call that has not converged, the running increment (ck_dep) is
// latched; the call then goes on to ck_impl_maxit, and at its end the latched state is
// compared with the final one ("### ck_ncref" line: max |de|/e over the cells, and the
// column heating error sum|de| dx / sum|e_final - e*| dx, max over columns).  Run it
// with a large ck_impl_maxit so that the final state is the converged reference.  0
// (default) = off, bitwise.  Print only.
inline int ck_impl_ncref = 0;
inline DvceArray4D<Real> *ck_ncref_ptr = nullptr;
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
// ---- ck-fast levers (tests_ck_implicit/README_fast.md).  All default off, and off each
// is bitwise the code without it.
// problem/ck_impl_xstep = k > 0 (lever 2): keep the STORED OPERATOR of ck_impl_frozen_op
// (kappa rho, the coefficient triple, the tm factorisation, the cut and the frozen beam
// deposit Qb) ACROSS calls.  A call re-applies the operator an earlier call stored, i.e.
// its pass 0 is a frozen (linear) pass instead of the storing sweep, unless (a) no
// operator is stored yet, (b) k or more cycles have passed since the store, or (c)
// problem/ck_impl_xstep_thr > 0 and some cell of the ck domain has moved by more than
// that relative amount in T or rho since the store (checked on the current state, after
// rt_pre_tp of pass 0).  Needs frozen_op + lin + jac_lin.  The KEY defaults to 8 since
// default-flips where ck_impl_every > 1 and the stored operator exists, 0 otherwise
// (deep_hot_jupiter_rt reads it).
inline int ck_impl_xstep = 0;
inline Real ck_impl_xstep_thr = 0.0;
// set per call by the pass function (pass 0): this call re-applies a stored operator
inline bool ck_impl_reuse_op = false;
inline int ck_impl_xs_cyc = -1;           // ncycle of the last store, -1 = none yet
inline Real ck_impl_xs_dmax = 0.0;        // the last measured max relative change
inline int64_t ck_impl_nstore = 0;        // calls that stored / re-applied
inline int64_t ck_impl_nreuse = 0;
inline DvceArray4D<Real> *ck_xsT_ptr = nullptr;   // T and rho at the store
inline DvceArray4D<Real> *ck_xsD_ptr = nullptr;
// problem/ck_impl_jreuse = rho > 0 (lever 3): the Jacobian is built on pass 0 of a call
// only and reused (chord) until the max residual contracts by less than rho per pass,
// after which the next pass rebuilds it.
inline Real ck_impl_jreuse = 0.0;
// problem/ck_impl_jreuse_xc: with ck_impl_xstep, a call that re-applies the stored
// operator also starts from the previous call's Jacobian instead of building one
inline bool ck_impl_jreuse_xc = false;
// problem/ck_impl_jreuse_act: with ck_impl_jreuse, also rebuild on any pass that follows
// one on which fewer than this fraction of the rank's columns were still active
inline Real ck_impl_jreuse_act = 0.25;
inline bool ck_impl_jac_built = false;
inline Real ck_impl_prev_res = -1.0;
// problem/ck_impl_pred (lever 4): no confirmation pass.  A column whose residual r_p
// (after the step of pass p) is predicted by the observed contraction to be below tol,
// r_p^2 / r_{p-1} <= ck_impl_pred_fac tol, is marked converged right after the step, so
// the call ends without a sweep whose only job is to verify it.  Needs ck_impl_fuse.
inline bool ck_impl_pred = false;
inline Real ck_impl_pred_fac = 0.5;
inline DvceArray3D<Real> *ck_rprev_ptr = nullptr;
inline int ck_impl_npred = 0;              // columns ended by the prediction, last pass
// problem/ck_impl_pred_chk (diagnostic, for the accuracy gates only): after a call the
// prediction ended, run ONE more sweep over every column and evaluate (no step)
// the residual and the energy gap of the state actually handed to the gas, so that the
// reported res / ckdesum (and ck_src, which the test hooks read) are the final state's.
inline bool ck_impl_pred_chk = false;
inline bool ck_impl_evalonly = false;
// problem/ck_impl_nosync (lever 5): no device allocation and no blocking scalar
// deep_copy in the RT pass path (cached dummies and host mirrors, stream-ordered
// fills), and the apply's clip count reduced into a device View.
inline bool ck_impl_nosync = false;
// problem/ck_impl_warm_step (lever 4, with ck_impl_warm): pass 0 of a warm-started call
// always takes a Newton step, even where the seeded state already meets tol, so that
// the seed is never accepted as it is (the residual norm is loose in the thin top,
// where e << ck_impl_norm_eps e_max).
inline bool ck_impl_warm_step = false;
// problem/ck_impl_cvkeep (lever 4, with ck_impl_cvsec): pass 0 of a call builds its rows
// with the secant heat capacity the previous call ended on (bounded to 1/50 .. 50 of
// e/T, else e/T as before) instead of e/T, which README_T1_stall measured 5x off in the
// dissociating top.  Only the Jacobian changes, not the balance solved.
inline bool ck_impl_cvkeep = false;
// ---- ck-cadence (tests_ck_implicit/README_cadence.md).  Default off (every = 1), and
// off it is bitwise the code without it.
// problem/ck_impl_every = N > 1: the full implicit call (with whatever levers are on)
// runs on the cycles with ncycle % N == 0 only, each with THAT STEP'S OWN dt.  After it,
// per cell, Q0 = (e - e*)/dt (the heating rate the gas actually received), D = dQ/dT at
// fixed rho (the Newton diagonal jac1, or -4 E/T for a thin cell; clamped <= 0), T0 and
// rho0 are stored.  On every other cycle each cell takes the linearised backward-Euler
// step
//     dT = (Q0 + D (T - T0)) dt / (c_v - D dt),   de = c_v dT,
// c_v = rho c_v(EOS) per kelvin, |de| capped at min(ck_impl_dtmax, ck_impl_demax) e.
// Each step's RT therefore covers exactly its own dt (nothing is lumped, lagged or
// counted twice); see the README for why this form and not a summed dt_rad.
// problem/ck_impl_every_thr = thr > 0: before the linearised step, a column any of
// whose RT cells has |T/T0 - 1| or |rho/rho0 - 1| > thr takes a full solve instead
// (only those columns: the others are masked out of the Newton by ck_done = 2).
// RESTART: Q0, D, T0, rho0 travel in the restart file (utils/two_stream_ck_rst.hpp),
// so a restarted run follows the schedule bitwise; a file without them (older) forces a
// full call first.
inline int ck_impl_every = 1;
inline Real ck_impl_every_thr = 0.0;
inline bool ck_cad_partial = false;          // this call solves the unmasked columns only
inline DvceArray3D<Real> *ck_cad_mask_ptr = nullptr;   // (m,k,j): 0 solve, 2 masked
inline DvceArray5D<Real> *ck_cad_ptr = nullptr;        // (m,5,k,j,i): Q0 D T0 rho0 Tn
inline DvceArray1D<Real> *ck_cad_stat_ptr = nullptr;   // per-step sums, see CkCadStep
inline int64_t ck_cad_nfull = 0;             // scheduled full calls
inline int64_t ck_cad_nlin = 0;              // linearised steps
inline int64_t ck_cad_nguard = 0;            // linearised steps with a guard call
inline double ck_cad_fref = 0.0;             // sum over linearised steps of nref/ncol
// ckhesdirk-1102 problem/ck_time_scheme = hesdirk2 (DhjCkT2Stage, deep_hot_jupiter_rt):
// ck_cad_force_lin = the cadence step of THIS call is linearised even on a full cycle
// (the second SDIRK stage); ck_t2_k1_ptr = the FSAL slope K1 (m,k,j,i) [energy/vol/s]
// and ck_t2_ok its validity, both carried by the restart block (slabs T2K1, T2OK).
inline bool ck_cad_force_lin = false;
inline bool ck_t2_on = false;
inline bool ck_t2_ok = false;
inline DvceArray4D<Real> *ck_t2_k1_ptr = nullptr;
// ---- ck-fast2 lever 1: problem/ck_dif_dtau = thr > 0, THE DEEP DIFFUSION HANDOVER.
// Default 0 (off), and off nothing below is allocated, launched or read (bitwise).
// Route B only (the ck two-stream owns the whole column down to the inner wall).  Per
// column, on every pass that builds the opacity, the chain cut ich is moved up to the
// shallowest face below which EVERY centre-to-centre layer is optically thick in EVERY
// chain: min over (band, g) of dtau_bg(f) = kr_{f-1} dz_{f-1}/2 + kr_f dz_f/2 >= thr for
// all faces f = is+1 .. ich (kr = kappa rho, the sweep's own; minus ck_dif_margin cells).
// Below ich the ck chains are replaced by the EXACT THICK LIMIT of the same discrete
// two-stream, band resolved and from the same face quantities:
//     F(f) = sum_b G_b(f) (B_b(f-1) - B_b(f)),
//     G_b(f) = sum_{g,q in b} wfc 2 mu_q (w_l + w_u - 1)/dtau_bg(f),
// wfc/mu_q the chains' flux weights and angles and w_l, w_u the BFace weights of the
// layer (1 wherever kappa rho varies smoothly), i.e. (4 pi/3) sum_g gw dB_b/dtau_bg for
// nquad = 2 and the 0.904 of nquad = 1.  The wall face keeps route B's datum, F(is) =
// sum_c wfc I_int,b.  The cells is .. ich-1 get the flux-form divergence of F with the
// chains' area/volume factors; the jacobian rows are exact (F is linear in B_b at frozen
// opacity).  THE HANDOVER IS ALGEBRAIC: the chains start at face ich with a flux
// boundary condition, R = 1 and Sc = g_c (B_b(ich-1) - B_b(ich)), g_c = 2 mu_q (w_l + w_u
// - 1)/dtau_c(ich), so the chains' own face flux at ich is sum_c wfc g_c dB = F(ich) to
// round-off, and the intensities above it start from the diffusion-limit field.
inline Real ck_dif_dtau = 0.0;
inline int ck_dif_margin = 0;
inline DvceArray3D<int> *ck_ich_ptr = nullptr;    // (m,k,j) chain cut
inline DvceArray5D<Real> *ck_difG_ptr = nullptr;  // (m,b,i,k,j) G_b at face i
inline DvceArray5D<Real> *ck_difE_ptr = nullptr;  // (m,b,i,k,j) emission weight, cell i
inline DvceArray4D<Real> *ck_difg_ptr = nullptr;  // (m,c,k,j) chain datum g_c; -1 = wall
inline DvceArray4D<int> *ck_diff_ptr = nullptr;   // (m,b,k,j) first thin face of band b
inline int64_t ck_dif_ncol = 0;                   // last store: columns with a handover
inline double ck_dif_fsum = 0.0;                  // ... and the sum of (ich - icut)

//! \fn CkDum
//! \brief a 1-element View for a capture that is never read.  CkDum<V>(label) builds a
//! fresh one (the old behaviour: an allocation, a zero fill and a free, each with a
//! device synchronisation) unless problem/ck_impl_nosync, in which case one cached View
//! per type is returned (never freed: it is leaked on purpose so that no View outlives
//! Kokkos::finalize).
template <typename V>
inline V CkMakeDummy(const std::string &lab) {
  if constexpr (V::rank == 1) {
    return V(lab, 1);
  } else if constexpr (V::rank == 2) {
    return V(lab, 1, 1);
  } else if constexpr (V::rank == 3) {
    return V(lab, 1, 1, 1);
  } else if constexpr (V::rank == 4) {
    return V(lab, 1, 1, 1, 1);
  } else {
    return V(lab, 1, 1, 1, 1, 1);
  }
}
template <typename V>
inline V CkDum(const std::string &lab) {
  if (!ck_impl_nosync) return CkMakeDummy<V>(lab);
  static V *p = nullptr;
  if (p == nullptr) p = new V(CkMakeDummy<V>("ck_cached_dummy"));
  return *p;
}
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
// problem/ck_impl_lin: THE LINEAR RE-APPLY KERNEL (tm sweep only, tests_ck_implicit/
// DESIGN_tm.md phase T3).  With the opacity and the beam frozen for the call
// (ck_impl_frozen_op), the tm sweep is AFFINE in the band Planck functions, Src = M B +
// s0, and its reflectance R does not depend on B at all: pass 1's R recurrence IS the
// factorisation of the column.  So after the storing pass a small kernel (ck_lin_build)
// parks R and 1/(1 + R beta) at every face and the three BFace/face-interpolation weights
// of every layer, and every later Newton pass that does not assemble the Jacobian runs
// rt_chain_ck_lin instead of the chain kernel: the Sc forward substitution and the ray
// back substitution on the new B_b, with no k-table look-up, no exponential, no divide
// and no beam.  The answer is the chain kernel's to round-off (the divides become
// multiplies by stored reciprocals).  Needs ck_impl_frozen_op, ck_impl_frozen_cof and
// ck_sweep_form = 1; the pass that builds the Jacobian still runs the chain kernel.
inline bool ck_impl_lin = false;
// problem/ck_impl_lin_check: on every pass the linear kernel runs, run the chain kernel
// on the same B first and print max|lin - chain|/max|chain| of Src, Fb, Em.  Gate only.
inline int ck_impl_lin_check = 0;
// problem/ck_impl_lin_thr: the linear kernel's threading.  4 = one thread per block of
// RT_NB chains, writing Src/Fb/Em directly (the chain kernel's layout); 1 = one thread
// per CHAIN, 4x the threads and a quarter of the registers, into per-chain partials that
// a second small kernel sums over the block (rt_chain_ck_lin_sum).  Same answer to
// round-off; the choice is a GPU occupancy one (tests_ck_implicit/README_T3.md).
inline int ck_impl_lin_thr = 1;
// problem/ck_impl_warm: start the Newton from the PREVIOUS call's converged increment.
// The fixed point does not move -- only the initial iterate does -- so a warm start can
// change the pass count and nothing else.  The seed is the total de the previous call
// applied, clipped to ck_impl_dtmax of the cell's energy; e^n is recorded net of it, so
// the balance being solved is unchanged.  NOT written to the restart file: after a
// restart (and on the first call of a run) the seed array is zero and the call is a cold
// start, which is correct, only slower for one step.
inline bool ck_impl_warm = false;

// problem/ck_impl_reuse_jac: QUASI-NEWTON.  The tridiagonal is rebuilt from scratch by
// every pass today, at the cost of ~10 atomic_add per (cell, chain) inside the sweep plus
// the CK_NB pair of Planck-fraction look-ups per cell that dB_b/dT needs.  With the
// exchange operator frozen (ck_impl_frozen_op) the ONLY thing in J that moves between
// passes is dB_b/dT, which is a smooth O(T^3) factor, so the matrix of the FIRST
// Newton pass is a good approximation to every later one.  This is the chord (modified
// Newton) method: the fixed point is untouched -- the residual is still re-formed by a
// full sweep every pass -- and only the convergence RATE can suffer.
//   0 = off, rebuild every pass (today's).
//   1 = CHORD: build J on the first pass that takes a Newton step, reuse the bits.
//   2 = SCALED CHORD: as 1, but each entry is rescaled by (T_col/T_col^0)^3, T_col being
//       the temperature of the CELL THE COLUMN OF J POINTS AT (entry 0 -> i-1,
//       1 -> i, 2 -> i+1).  dB_b/dT = d(sigma T^4 f_b(T)/pi)/dT is 4 sigma T^3 f_b/pi to
//       within the band fraction's own (weak) T dependence, so this recovers the dominant
//       part of the update for three multiplies per row and no sweep work at all.
// The M-matrix property is preserved by both (the rescaling is positive), so the Thomas
// sweep still needs no pivoting.
inline int ck_impl_reuse_jac = 0;
// problem/ck_impl_seed: THE INITIAL GUESS.  ck_impl_warm starts the Newton from the
// PREVIOUS call's increment, which perturbs calls that do not converge (README_phase3
// section 4).  This instead starts it from a guess formed inside THIS call, from this
// call's own first sweep, and so carries none of that history:
//   0 = off: pass 0 takes the ordinary Newton step (today's).
//   1 = SEMI-IMPLICIT: pass 0 takes the step the DEFAULT (ck_implicit = false) scheme
//       would take, de = (S/lambda)(1 - e^{-lambda bdt}), lambda = 4 E/e -- the
//       rt_semi_lin form of the apply block, per cell, absorbed field lagged.
//   2 = EXACT LAGGED: pass 0 takes the per-cell BRACKETED backward-Euler step
//       CkThinSolve, i.e. the exact root of x - e* - bdt(A - E (x/e)^4) with A lagged.
//       Strictly better than 1 (it is the same balance without the (1-e^-x) damping and
//       without the linearisation), positivity-preserving, and it costs 60 bisections in
//       a kernel that is already memory bound.
// Both are applied through the SAME Thomas sweep as an identity row (a = c = 0, b = 1),
// so the per-pass caps ck_impl_dtmax / ck_impl_demax act on them exactly as on a Newton
// step and nothing downstream knows the difference.  The residual of the seeded iterate
// is measured at the top of pass 1, which is where every other pass is measured.
inline int ck_impl_seed = 0;
// ---- phase T4 (tests_ck_implicit/README_T4.md): the cost of a call.  All default off,
// in which case the call is T3's to the bit.
// problem/ck_impl_fuse: ONE kernel for the residual test and the tridiagonal step, one
// thread TEAM per column.  The rows are assembled and the step is capped and applied in
// parallel over the cells (coalesced loads, CkThinSolve in parallel), and only the
// Thomas recurrences run on one lane, out of team scratch.  Same arithmetic in the same
// order as ck_impl_res + ck_impl_tri, so the state is bitwise theirs; only the order of
// the two diagnostic sums (slots 4 and 5) changes.  Refused with ck_impl_debug > 0.
inline bool ck_impl_fuse = false;
// problem/ck_impl_jac_lin: build the tridiagonal from the STORED FACTORISATION of the
// linear kernel instead of in the JAC instantiation of the chain kernel.  The JAC
// recurrences of the tm sweep (dSc/dB carried up, dd^+/dB and the reflection scalar W
// carried down; see the JAC blocks in rt_chain_ck) read only the half-layer triple, R,
// 1/(1 + R beta), the BFace and face-interpolation weights, beta, the area ratio and 1/dz
// -- all of which ck_lin_build already parks in lP/lG.  So a separate one-chain-per-
// thread kernel re-runs them on the stored data, writes each chain's three entries to
// per-chain partials, and a small kernel sums them in chain order (no atomics).  Same
// entries as the JAC pass to round-off (the divides by 1 + R beta become multiplies by
// the stored reciprocal), negative per-chain off-diagonal parts dropped exactly as there.
// The residual of that pass comes from the linear kernel.  Needs ck_impl_lin with
// ck_impl_lin_thr = 1.
inline bool ck_impl_jac_lin = false;
// problem/ck_impl_jneg: under ck_impl_jac_lin, keep the NEGATIVE per-chain off-diagonal
// parts the JAC assembly drops, and drop only a negative NET entry (the M-matrix).
inline bool ck_impl_jneg = false;
// problem/ck_impl_cvsec: the heat capacity of the Jacobian rows.  cv = e/T (the header
// note) is off by up to ~5x where H2 dissociates (README_T1_stall.md: d ln T/d ln e = 5.4
// in the day-side top), and it enters every row.  With this on, each cell's cv is the
// SECANT of its own last two iterates, (e_p - e_{p-1})/(T_p - T_{p-1}), taken from pass 1
// on (pass 0 -> 1 is the seed step, a large and clean secant) and kept from the previous
// pass wherever the step is too small to resolve it (|dT| < 1e-7 T) or the slope is not
// within 1/50 .. 50 of e/T.  Jacobian only: the fixed point does not move.  Needs
// ck_impl_fuse.
inline bool ck_impl_cvsec = false;
// problem/ck_impl_rsec: a SECANT BOUND on each thick row's diagonal.  Default 0 = off
// (bitwise).  The Jacobian drops the opacity's own temperature dependence (header note,
// item iv).  In the 10x WASP-121b day-side upper layers (T 5000-6500 K, t_rad 1-4 % of
// the stage) that term makes the true dR_i/de_i larger than the modelled one, so the
// Newton overshoots: a period-2 limit cycle with the e/T Jacobian, ~1e-6 residual stalls
// with cvsec (SMOKE_DIAG.md section 6).  With rsec > 0, from pass 1 of a call on, each
// thick row compares its diagonal b with the secant of its OWN residual over the last
// pass, (R_p - R_{p-1})/(e_p - e_{p-1}), and takes the secant where it is LARGER, capped
// at rsec*b.  A larger diagonal keeps the M-matrix and only shortens a step; the
// residual, and therefore the fixed point, is untouched.  The secant also carries what
// the neighbours' steps did to R_i, which is why it is a bound (max) and capped, not a
// replacement.  Needs ck_impl_fuse, ck_impl_glob = none, ck_impl_aa = 0.
inline Real ck_impl_rsec = 0.0;
// problem/ck_impl_jac0: with ck_impl_seed > 0, build the Jacobian on pass 0 (at e^n,
// together with the storing pass) instead of on pass 1 (at the seeded state).  Merges
// the two full sweeps of T3 into one; the chord matrix is then one seed step stale.
inline bool ck_impl_jac0 = false;
// problem/ck_impl_lw: launch the implicit-only kernels (the fused step, the linear
// kernel, its block sum and the Jacobian kernels) as LIGHT-WEIGHT Kokkos kernels.  A
// functor of 512 B .. 32 kB otherwise takes the HIP constant-memory path, which waits for
// the previous such kernel and copies the functor before every launch.  Under 4 kB the
// light-weight path passes it as a kernel argument instead.  Same kernel body.
inline bool ck_impl_lw = false;
// problem/ck_impl_glob: NEWTON GLOBALISATION of the fused column step (tests_ck_implicit/
// README_glob.md).  0 = none (default, bitwise T4); 1 = ls; 2 = ls_sub.
//  ls: a per-column BACKTRACKING LINE SEARCH on the merit phi = ||R_i/(e*_i + eps
//    e*_max)||_2 (e* = e^n, fixed over the call, so phi is one function of the iterate).
//    The first trial is the capped Newton step (ck_impl_dtmax / ck_impl_demax, i.e. T4's
//    step exactly); the NEXT pass's sweep evaluates it, and a column whose merit did not
//    drop to (1 - c alpha) phi_old is pulled back to half its trial (no Newton step that
//    pass, no new Jacobian: the stored factorisation is reused), up to ck_impl_ls_ntry
//    halvings, after which the last trial is accepted.  Accepted columns take their next
//    step in the same pass, so a rejection costs no extra sweep for anyone else.
//  ls_sub: ls, plus per-column SUB-STEPPING: a column still above tol after
//    ck_impl_maxit passes of its (sub-)step is put back to e^n and redone as 2, 4, ...
//    (at most ck_impl_sub_max) backward-Euler sub-steps of bdt/2^L, each with its own
//    maxit budget; the source of the finished sub-steps is carried in ck_sacc, so the
//    ckdesum gap still compares sum (e - e^n) dx with the source the gas received.  A
//    column that fails at the finest level is left at its last iterate and counted.
// Needs ck_impl_fuse (the other step path is not globalised).
inline int ck_impl_glob = 0;
inline int ck_impl_ls_ntry = 4;
inline Real ck_impl_ls_c = 1.0e-4;
inline int ck_impl_sub_max = 32;
// per-column state (m, q, k, j), q = 0 merit at the accepted iterate, 1 pending trial's
// alpha (0 = none), 2 halvings so far, 3 sub-step level L, 4 sub-step index, 5 residual
// evaluations in this sub-step
inline DvceArray4D<Real> *ck_lsc_ptr = nullptr;
// the increment of the pending trial, and (ls_sub) the source of the finished sub-steps
inline DvceArray4D<Real> *ck_lsd_ptr = nullptr;
inline DvceArray4D<Real> *ck_sacc_ptr = nullptr;
// host counters over one call: rejected trials, sub-step restarts, failed columns
inline int ck_impl_nrej = 0;
inline int ck_impl_nsub = 0;
inline int ck_impl_nsubfail = 0;
// ls_sub: the chord Jacobian (ck_impl_reuse_jac) was built at the whole step's seeded
// state; a sub-step starts with its own seed step (ck_impl_seed = 2), and the pass after
// one rebuilds the Jacobian at the seeded state (for every live column)
inline bool ck_impl_jac_again = false;
// problem/ck_impl_esc (tests_ck_implicit/README_jac.md): CHEAPER ESCALATION for ls_sub.
// 0 = off (default, bitwise ls_sub).  1 = (i) a (sub-)step that has used its maxit
// residual evaluations but whose merit still falls by ck_impl_esc_rho or better per pass,
// and whose projected passes to tol at that rate fit in ck_impl_esc_extra more, keeps
// iterating instead of being split; (ii) a split restarts ONLY the failing sub-step (from
// its own start state e^n + ck_sacc, redone as two sub-steps of half its length), not the
// whole column from e^n; (iii) after a sub-step that converged in <= maxit/2 steps (its
// first residual evaluation is not a step) the column returns to the coarser level where
// the sub-step grid allows it.
// Needs ck_impl_glob = ls_sub.
inline int ck_impl_esc = 0;
inline Real ck_impl_esc_rho = 0.9;
inline int ck_impl_esc_extra = 8;
// problem/ck_impl_aa (README_jac.md): ANDERSON (Pulay / DIIS) ACCELERATION of the chord
// iteration, depth ck_impl_aa (0 = off, default, bitwise; at most CK_AA_MAX).  Each
// column keeps the last aa differences of its iterates dX and of its backward-Euler
// residuals dR
// (the residual every pass's ordinary sweep already evaluates, so no extra sweep).  The
// step is taken from the combination of the past iterates that minimises the weighted
// linearised residual, gamma = argmin || W (r - dR gamma) ||_2 (W = the merit weights
// 1/(e*_i + eps e*_max)), and is the tridiagonal solve of that minimised residual:
//     e_new = e - dX gamma + P^{-1} (-(r - dR gamma)),  P the chord tridiagonal.
// On a linear problem with a fixed P this is GMRES preconditioned by the tridiagonal
// (Walker & Ni 2011): the non-local two-stream coupling and the dropped negative
// per-chain parts that the tridiagonal leaves out are picked up from the secant
// information, so the iteration is superlinear instead of linear.  Thin rows (whose step
// is the per-cell thin solve) are left out of the least squares and are not mixed.  The
// history restarts with every new backward-Euler residual (pass 0, every sub-step), after
// a fallback step and whenever the weighted residual norm grew over the last step (the
// secants then span a region where the column is far from linear); a seed step's pair is
// never stored (the seed is not a step of this iteration).  Needs ck_impl_fuse.
constexpr int CK_AA_MAX = 8;
inline int ck_impl_aa = 0;
// problem/ck_impl_aa_rst: restart the history when the weighted residual norm grew
inline bool ck_impl_aa_rst = true;
// ck_aah: (m, s, k, j, i), s = 0..aa-1 dX, aa..2aa-1 dR, 2aa the last iterate, 2aa+1 its
// residual; ck_aac: (m, q, k, j), q = 0 stored pairs, 1 ring head, 2 last iterate valid,
// 3 its weighted residual norm^2
inline DvceArray5D<Real> *ck_aah_ptr = nullptr;
inline DvceArray4D<Real> *ck_aac_ptr = nullptr;
// host counters over one call: extra passes granted, coarsenings, accelerated steps
inline int ck_impl_nesc = 0;
inline int ck_impl_ncoarse = 0;
inline int ck_impl_naa = 0;

// the Newton pass index inside one RT call; read by the pass function to decide whether
// the opacity is rebuilt.  -1 = ck_implicit is off.
inline int ck_impl_pass = -1;
// the pass that BUILDS the Jacobian under ck_impl_reuse_jac: the first one that takes a
// Newton step, i.e. pass 1 when pass 0 is a seed step and pass 0 otherwise.
// problem/ck_impl_jac0 moves it to pass 0 in any case.
inline int CkImplJacPass() { return (ck_impl_seed > 0 && !ck_impl_jac0) ? 1 : 0; }

//----------------------------------------------------------------------------------------
//! \fn void CkParFor4
//! \brief the 4-D par_for of athena.hpp, launched LIGHT-WEIGHT when lw (ck_impl_lw) and
//! by par_for itself otherwise.  For the implicit-only kernels.
template <typename Function>
inline void CkParFor4(const std::string &name, const bool lw, const int nl, const int nu,
                      const int kl, const int ku, const int jl, const int ju,
                      const int il, const int iu, const Function &function) {
  if (!lw) {
    par_for(name, DevExeSpace(), nl, nu, kl, ku, jl, ju, il, iu, function);
    return;
  }
  const int nn = nu - nl + 1;
  const int nk = ku - kl + 1;
  const int nj = ju - jl + 1;
  const int ni = iu - il + 1;
  const int nkji = nk*nj*ni;
  const int nji  = nj*ni;
  const int nnkji = nn*nk*nj*ni;
  auto pol = Kokkos::Experimental::require(
      Kokkos::RangePolicy<DevExeSpace>(DevExeSpace(), 0, nnkji),
      Kokkos::Experimental::WorkItemProperty::HintLightWeight);
  Kokkos::parallel_for(name, pol, KOKKOS_LAMBDA(const int &idx) {
    int n = (idx)/nkji;
    int k = (idx - n*nkji)/nji;
    int j = (idx - n*nkji - k*nji)/ni;
    int i = (idx - n*nkji - k*nji - j*ni) + il;
    n += nl;
    k += kl;
    j += jl;
    function(n, k, j, i);
  });
}

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
// problem/ck_impl_reuse_jac = 2: the temperature the reused Jacobian was built at
inline DvceArray4D<Real> *ck_t0_ptr = nullptr;
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
// ---- problem/ck_impl_lin: the factorisation, packed so that the linear kernel
// captures few Views.  lP (m, q*nch + chain, i, k, j), q = 0..5: R on the below side
// of face i (i = icut .. ie+1) and, slots 2-4, the layer joining cells i and i+1
// (i = icut .. ie-1): dslv/dB_i, dsuu/dB_{i+1} and the face interpolation dt_l/dtc.
// ck-next: slots 1 (1/(1 + R beta_i)) and 5 (2 (wfc/mu) kappa rho, the emission) are
// no longer stored: the readers form them (CkIdn, CkEmW in two_stream_rt.hpp).
// lG (m, q, i, k, j), q = 0..4: beta at face i, A_{i-1}/A_i, the flux frame factor,
// 1/dz and dz (ck-next).  lC (q, chain): the flux weight, the internal-flux datum
// at the cut and (ck-next) the emission weight 2 (wfc/mu).
// The half-layer triple (e0, cin, cout; slots 0-2 before ck-store) is not packed: the
// kernels read it from ck_c0, ck_ci, ck_co, which hold the same numbers for as long as
// the factorisation lives.
inline DvceArray5D<Real> *ck_linP_ptr = nullptr;
inline DvceArray5D<Real> *ck_linG_ptr = nullptr;
inline DvceArray2D<Real> *ck_linC_ptr = nullptr;
// ck_impl_lin_thr = 1: the per-chain partial Src and Fb, (m, chain, i, k, j)
inline DvceArray5D<Real> *ck_lps_ptr = nullptr;
inline DvceArray5D<Real> *ck_lpf_ptr = nullptr;
// ck_impl_lin_check: the chain kernel's Src, Fb, Em on the checked pass, (m, 3*nblk, ...)
inline DvceArray5D<Real> *ck_lchk_ptr = nullptr;
// ck_impl_cvsec: the previous iterate's e and T, and the cv the rows use, (m,k,j,i)
inline DvceArray4D<Real> *ck_ep_ptr = nullptr;
// problem/ck_impl_rsec: each thick row's residual and energy at the previous pass
inline DvceArray4D<Real> *ck_rsr_ptr = nullptr;
inline DvceArray4D<Real> *ck_rse_ptr = nullptr;
inline DvceArray4D<Real> *ck_tp_ptr = nullptr;
inline DvceArray4D<Real> *ck_cv_ptr = nullptr;
// ck_impl_jac_lin: the per-chain partial of the third Jacobian entry (the first two go to
// ck_lpf and ck_lps, free by then), (m, chain, i, k, j)
inline DvceArray5D<Real> *ck_lpj_ptr = nullptr;
// ck-jlin: the per-chain Jacobian coefficients of rt_chain_ck_jlin BEFORE the dB_b/dT
// factor, (m, chain, i, k, j): dSrc_i/dB_{i-1} (masked to 0 where the row drops it) and
// dSrc_i/dB_i; dSrc_i/dB_{i+1} (masked) lives in ck_lpj.  They depend on the stored
// factorisation alone, so rt_chain_ck_jlin runs once per stored operator (ck_jlc_ok)
// and every pass that builds the tridiagonal only multiplies them by dB/dT and sums.
inline DvceArray5D<Real> *ck_ljm_ptr = nullptr;
inline DvceArray5D<Real> *ck_lj0_ptr = nullptr;
inline bool ck_jlc_ok = false;
//! \struct CkJlP1
//! \brief the state of pass 1 of rt_chain_ck_jlin: dSc/dB in the window (B_{i-2},
//! B_{i-1}, B_i) below face i, and the face weights the layer below hands up
struct CkJlP1 {
  Real s0, s1, s2, fc0, fc1, uc0, uc1;
};
//! \fn void CkJlP1Step
//! \brief one layer (cell i) of pass 1, shared by rt_chain_ck_jlin and ck_lin_build
//! (which fuses pass 1 on the pass that builds the tridiagonal): the same expressions
//! in one place.  up = (i < ie); wl, wu, ffj are then lP slots 5-7 of the layer i, i+1;
//! face i maps R to r1 and Sc to c1 Sc (CkJlP1Step: the Eddington face; CkVefMap).
KOKKOS_INLINE_FUNCTION
void CkJlP1StepRC(CkJlP1 &w, const bool up, const Real wl, const Real wu, const Real ffj,
                  const Real ci, const Real co, const Real e0, const Real r1,
                  const Real c1) {
  // slots (B_{i-1}, B_i, B_{i+1}); the first is 0 for all three
  Real dlv1 = 1.0, dlv2 = 0.0, duu1 = 1.0, duu2 = 0.0;
  Real dfw1 = 1.0, dfw2 = 0.0;
  if (up) {
    const Real dfl = (1.0 - ffj)*wl + ffj*(1.0 - wu);
    dlv1 = wl;
    dlv2 = 1.0 - wl;
    duu1 = 1.0 - wu;
    duu2 = wu;
    dfw1 = dfl;
    dfw2 = 1.0 - dfl;
  }
  const Real tj = 1.0 - e0;
  const Real r2 = tj*tj*r1;
  auto upd = [&](const Real dsf, const Real dsu, const Real dlv, const Real dfw,
                 const Real so) {
    const Real dpl = ci*dsf + co*dsu;
    const Real dql = ci*dsu + co*dsf;
    const Real dpu = ci*dlv + co*dfw;
    const Real dqu = ci*dfw + co*dlv;
    const Real s2 = tj*(r1*dql + c1*so) + dpl;
    return tj*(r2*dqu + s2) + dpu;
  };
  const Real n0 = upd(w.fc0, w.uc0, 0.0, 0.0, w.s1);
  const Real n1 = upd(w.fc1, w.uc1, dlv1, dfw1, w.s2);
  const Real n2 = upd(0.0, 0.0, dlv2, dfw2, 0.0);
  w.s0 = n0;
  w.s1 = n1;
  w.s2 = n2;
  w.fc0 = dfw1;
  w.fc1 = dfw2;
  w.uc0 = duu1;
  w.uc1 = duu2;
}
//! \fn void CkJlP1Step
//! \brief CkJlP1StepRC with the Eddington face's map of face i, r1 = (R + beta)/(1 +
//! R beta) and c1 = (1 - beta)/(1 + R beta) (idn = 1/(1 + R beta)); kept as its own
//! copy so that the default kernels are the code they always were
KOKKOS_INLINE_FUNCTION
void CkJlP1Step(CkJlP1 &w, const bool up, const Real wl, const Real wu, const Real ffj,
                const Real ci, const Real co, const Real e0, const Real rj,
                const Real idn, const Real bt) {
  // slots (B_{i-1}, B_i, B_{i+1}); the first is 0 for all three
  Real dlv1 = 1.0, dlv2 = 0.0, duu1 = 1.0, duu2 = 0.0;
  Real dfw1 = 1.0, dfw2 = 0.0;
  if (up) {
    const Real dfl = (1.0 - ffj)*wl + ffj*(1.0 - wu);
    dlv1 = wl;
    dlv2 = 1.0 - wl;
    duu1 = 1.0 - wu;
    duu2 = wu;
    dfw1 = dfl;
    dfw2 = 1.0 - dfl;
  }
  const Real tj = 1.0 - e0;
  const Real r1 = (rj + bt)*idn;
  const Real r2 = tj*tj*r1;
  const Real c1 = (1.0 - bt)*idn;
  auto upd = [&](const Real dsf, const Real dsu, const Real dlv, const Real dfw,
                 const Real so) {
    const Real dpl = ci*dsf + co*dsu;
    const Real dql = ci*dsu + co*dsf;
    const Real dpu = ci*dlv + co*dfw;
    const Real dqu = ci*dfw + co*dlv;
    const Real s2 = tj*(r1*dql + c1*so) + dpl;
    return tj*(r2*dqu + s2) + dpu;
  };
  const Real n0 = upd(w.fc0, w.uc0, 0.0, 0.0, w.s1);
  const Real n1 = upd(w.fc1, w.uc1, dlv1, dfw1, w.s2);
  const Real n2 = upd(0.0, 0.0, dlv2, dfw2, 0.0);
  w.s0 = n0;
  w.s1 = n1;
  w.s2 = n2;
  w.fc0 = dfw1;
  w.fc1 = dfw2;
  w.uc0 = duu1;
  w.uc1 = duu2;
}
//! \fn void CkVefMap
//! \brief problem/ck_sph_face = 5 (two_stream_rt::ck_sph_face): the face map of the VEF
//! jump conditions on the relation u = R d + Sc below a face, rp = rho' = (A_b/A_a)
//! s_b/s_a and gm = gamma: above it R' = r1, Sc' = c1 Sc (the chain kernel's tm map)
KOKKOS_INLINE_FUNCTION
void CkVefMap(const Real rj, const Real rp, const Real gm, Real &r1, Real &c1) {
  const Real ev = (1.0 + rj) + gm*rp*(1.0 - rj);
  r1 = 1.0 - 2.0*gm*rp*(1.0 - rj)/ev;
  c1 = 2.0*rp/ev;
}
//! \fn CkJlP1 CkJlP1Init
//! \brief the window at the cut; gdf >= 0: the flux datum of a ck_dif_dtau handover
KOKKOS_INLINE_FUNCTION
CkJlP1 CkJlP1Init(const Real gdf) {
  CkJlP1 w{0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
  if (gdf >= 0.0) {
    w.s1 = gdf;
    w.s2 = -gdf;
  }
  return w;
}
//! \fn Real CkMulRn
//! \brief a*b rounded on its own, never contracted into an FMA with a following add:
//! the product ck_jlin_sum adds is then bitwise the one rt_chain_ck_jlin used to store
KOKKOS_INLINE_FUNCTION
double CkMulRn(const double a, const double b) {
#if defined(__CUDA_ARCH__)
  return __dmul_rn(a, b);
#else
  return a*b;
#endif
}
KOKKOS_INLINE_FUNCTION
float CkMulRn(const float a, const float b) {
#if defined(__CUDA_ARCH__)
  return __fmul_rn(a, b);
#else
  return a*b;
#endif
}
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
                        const int n2, const int n3, const int nbk = 1) {
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
    delete ck_t0_ptr;
    if (ck_kro_ptr != nullptr) {
      delete ck_kro_ptr;
      delete ck_c0_ptr;
      delete ck_ci_ptr;
      delete ck_co_ptr;
      delete ck_tpf_ptr;
      ck_kro_ptr = nullptr;
    }
    if (ck_linP_ptr != nullptr) {
      delete ck_linP_ptr;
      delete ck_linG_ptr;
      delete ck_linC_ptr;
      ck_linP_ptr = nullptr;
    }
    if (ck_lps_ptr != nullptr) {
      delete ck_lps_ptr;
      delete ck_lpf_ptr;
      ck_lps_ptr = nullptr;
    }
    if (ck_lchk_ptr != nullptr) {
      delete ck_lchk_ptr;
      ck_lchk_ptr = nullptr;
    }
    if (ck_lpj_ptr != nullptr) {
      delete ck_lpj_ptr;
      delete ck_ljm_ptr;
      delete ck_lj0_ptr;
      ck_lpj_ptr = nullptr;
      ck_ljm_ptr = nullptr;
      ck_lj0_ptr = nullptr;
    }
    ck_jlc_ok = false;
    if (ck_osc_ptr != nullptr) {
      delete ck_osc_ptr;
      ck_osc_ptr = nullptr;
    }
    if (ck_rsr_ptr != nullptr) {
      delete ck_rsr_ptr;
      delete ck_rse_ptr;
      ck_rsr_ptr = nullptr;
      ck_rse_ptr = nullptr;
    }
    if (ck_ep_ptr != nullptr) {
      delete ck_ep_ptr;
      delete ck_tp_ptr;
      delete ck_cv_ptr;
      ck_ep_ptr = nullptr;
    }
    if (ck_lsc_ptr != nullptr) {
      delete ck_lsc_ptr;
      delete ck_lsd_ptr;
      ck_lsc_ptr = nullptr;
    }
    if (ck_sacc_ptr != nullptr) {
      delete ck_sacc_ptr;
      ck_sacc_ptr = nullptr;
    }
    if (ck_aah_ptr != nullptr) {
      delete ck_aah_ptr;
      delete ck_aac_ptr;
      ck_aah_ptr = nullptr;
      ck_aac_ptr = nullptr;
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
  // problem/ck_impl_reuse_jac = 2 only; a 1-element dummy otherwise
  ck_t0_ptr = new DvceArray4D<Real>("ck_t0", (ck_impl_reuse_jac == 2) ? nmb : 1,
                                    (ck_impl_reuse_jac == 2) ? n3 : 1,
                                    (ck_impl_reuse_jac == 2) ? n2 : 1,
                                    (ck_impl_reuse_jac == 2) ? n1 : 1);
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
    // problem/ck_impl_lin: six more Reals per (cell, chain), four per cell
    if (ck_impl_lin) {
      ck_linP_ptr = new DvceArray5D<Real>("ck_lP", nmb, 6*nch, n1, n3, n2);
      ck_linG_ptr = new DvceArray5D<Real>("ck_lG", nmb, 5, n1, n3, n2);
      ck_linC_ptr = new DvceArray2D<Real>("ck_lC", 3, nch);
      if (ck_impl_lin_thr == 1) {
        ck_lps_ptr = new DvceArray5D<Real>("ck_lps", nmb, nch, n1, n3, n2);
        ck_lpf_ptr = new DvceArray5D<Real>("ck_lpf", nmb, nch, n1, n3, n2);
      }
      if (ck_impl_lin_check > 0) {
        ck_lchk_ptr = new DvceArray5D<Real>("ck_lchk", nmb, 3*nbk, n1, n3, n2);
      }
      if (ck_impl_jac_lin && ck_impl_lin_thr == 1) {
        ck_lpj_ptr = new DvceArray5D<Real>("ck_lpj", nmb, nch, n1, n3, n2);
        ck_ljm_ptr = new DvceArray5D<Real>("ck_ljm", nmb, nch, n1, n3, n2);
        ck_lj0_ptr = new DvceArray5D<Real>("ck_lj0", nmb, nch, n1, n3, n2);
        ck_jlc_ok = false;
      }
    }
  }
  if (ck_impl_osc > 0) {
    ck_osc_ptr = new DvceArray4D<Real>("ck_osc", nmb, n3, n2, n1);
  }
  if (ck_impl_rsec > 0.0) {
    ck_rsr_ptr = new DvceArray4D<Real>("ck_rsr", nmb, n3, n2, n1);
    ck_rse_ptr = new DvceArray4D<Real>("ck_rse", nmb, n3, n2, n1);
  }
  if (ck_impl_cvsec) {
    ck_ep_ptr = new DvceArray4D<Real>("ck_ep", nmb, n3, n2, n1);
    ck_tp_ptr = new DvceArray4D<Real>("ck_tp", nmb, n3, n2, n1);
    ck_cv_ptr = new DvceArray4D<Real>("ck_cv", nmb, n3, n2, n1);
  }
  // problem/ck_impl_glob: the per-column line-search / sub-step state
  if (ck_impl_glob > 0) {
    ck_lsc_ptr = new DvceArray4D<Real>("ck_lsc", nmb, 6, n3, n2);
    ck_lsd_ptr = new DvceArray4D<Real>("ck_lsd", nmb, n3, n2, n1);
    if (ck_impl_glob == 2) {
      ck_sacc_ptr = new DvceArray4D<Real>("ck_sacc", nmb, n3, n2, n1);
    }
  }
  // problem/ck_impl_pred: the previous pass's column residual
  if (ck_impl_pred) {
    if (ck_rprev_ptr != nullptr) delete ck_rprev_ptr;
    ck_rprev_ptr = new DvceArray3D<Real>("ck_rprev", nmb, n3, n2);
  }
  // problem/ck_impl_aa: the per-column Anderson history
  if (ck_impl_aa > 0) {
    ck_aah_ptr = new DvceArray5D<Real>("ck_aah", nmb, 2*ck_impl_aa + 2, n3, n2, n1);
    ck_aac_ptr = new DvceArray4D<Real>("ck_aac", nmb, 4, n3, n2);
  }
  if (ck_conv_ptr == nullptr) {
    // slots 8-11: ck_impl_glob's rejected trials, sub-step restarts, failed columns,
    // sub-step seed steps; 12-14: ck_impl_esc's extra passes and coarsenings, and
    // ck_impl_aa's accelerated steps
    // 16-17: ck_impl_floorbound / ck_impl_kkt_demax (cells at the floor, KKT cells)
    // 18: ck_impl_osc (cells damped)
    ck_conv_ptr = new DvceArray1D<Real>("ck_conv", 19);
  }
  // ck-fast2 lever 1 (problem/ck_dif_dtau)
  if (ck_dif_dtau > 0.0) {
    if (ck_ich_ptr != nullptr) {
      delete ck_ich_ptr;
      delete ck_difG_ptr;
      delete ck_difE_ptr;
      delete ck_difg_ptr;
      delete ck_diff_ptr;
    }
    ck_ich_ptr = new DvceArray3D<int>("ck_ich", nmb, n3, n2);
    ck_difG_ptr = new DvceArray5D<Real>("ck_difG", nmb, nb, n1, n3, n2);
    ck_difE_ptr = new DvceArray5D<Real>("ck_difE", nmb, nb, n1, n3, n2);
    ck_difg_ptr = new DvceArray4D<Real>("ck_difg", nmb, nch, n3, n2);
    ck_diff_ptr = new DvceArray4D<int>("ck_diff", nmb, nb, n3, n2);
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

void CkWarmSeed(Mesh *pm, DvceArray5D<Real> u0);
// (defined in utils/two_stream_rt_driver.cpp, tusplit-1102)

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

//----------------------------------------------------------------------------------------
//! \fn bool CkKktCell
//! \brief the KKT test of ck_impl_floorbound / ck_impl_kkt_demax for one cell (the rule
//! of the fused kernel's kktcell), for ck_impl_kkt_row and ck_impl_stalldbg.

KOKKOS_INLINE_FUNCTION
bool CkKktCell(const EOS_Data &eos, const bool fb, const bool kd, const Real detot,
               const Real rho, const Real e, const Real es, const Real r) {
  bool kkt = false;
  if (fb && r > 0.0) {
    Real efl = eos.EnergyFromTemperature(rho, eos.tfloor);
    if (e < eos.EnergyFloorBound(rho)) {
      const Real ep = eos.EnergyFromPressure(rho, eos.pfloor);
      if (ep > efl) efl = ep;
    }
    kkt = (e <= efl*(1.0 + 1.0e-10));
  }
  if (kd && !kkt && es > 0.0) {
    kkt = (r > 0.0 && e <= es*(1.0 - detot)*(1.0 + 1.0e-10))
          || (r < 0.0 && e >= es*(1.0 + detot)*(1.0 - 1.0e-10));
  }
  return kkt;
}

int CkImplStep(Mesh *pm, DvceArray5D<Real> u0, DvceArray3D<int> icut_,
               DvceArray4D<Real> T_, DvceArray4D<Real> dx1_, const Real bdt);
// (defined in utils/two_stream_rt_driver.cpp, tusplit-1102)

}  // namespace two_stream_rt

#endif  // UTILS_TWO_STREAM_COLUMN_CK_HPP_
