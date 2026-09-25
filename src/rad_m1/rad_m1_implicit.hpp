#ifndef RAD_M1_RAD_M1_IMPLICIT_HPP_
#define RAD_M1_RAD_M1_IMPLICIT_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_implicit.hpp
//! \brief constants and device-inline helpers of the IMPLICIT x1 column transport,
//! <rad_m1>/transport = implicit_x1 (docs/dev/rad_m1_implicit_design.md, milestone 3a).
//!
//! Nothing in here is reachable from the explicit path: the file is included only by
//! rad_m1_implicit.cpp and by the small delimited hooks of rad_m1.cpp / rad_m1_tasks.cpp
//! / rad_m1_newdt.cpp, so `transport = explicit` stays bitwise what it was.

#include <math.h>

#include "athena.hpp"
#include "eos/eos_table.hpp"

namespace radm1 {

// implicit_predictor = step as RESTART STATE: the marker of the little header block
// (int32 have, int32 ncomp, Real pred_dt) written behind the two-stream warm-start
// header, and read by the same eight-byte peek (src/outputs/restart.cpp, pgen.cpp).
// Present only when a stored increment exists; its ncomp per-cell slabs of ipred then
// follow the warm-start slabs at the end of every MeshBlock record.
constexpr char kM1PredRstMagic[8] = {'M', '1', 'P', 'R', 'E', 'D', '0', '1'};

// <rad_m1>/time_scheme (docs/dev/rad_m1_time2_design.md, tests_m1/runs_3x_hesdirk2)
constexpr int M1_TIME_BE = 0;         // Lie splitting: Heun hydro, then one BE solve
constexpr int M1_TIME_HESDIRK2 = 1;   // H-ESDIRK2: one stage solve inside each Heun stage
// what the next ImplicitSolve does (RadiationM1::t2_solve)
constexpr int M1_T2S_NONE = 0;        // plain backward Euler (time_scheme = be)
constexpr int M1_T2S_BESTORE = 1;     // backward Euler, and K1 = (Y - rhs)/dt is stored
constexpr int M1_T2S_STAGE1 = 2;      // stage solve 1 (gamma dt), K2 stored
constexpr int M1_T2S_STAGE2 = 3;      // stage solve 2 (gamma dt), K1 (FSAL) stored
// channels of t2k1 / t2k2 / t2inc: cell E, gas momentum and total energy, and the face
// fluxes F0 of x1/x2/x3, each face stored at its own (lower) index
constexpr int M1_T2_E = 0;
constexpr int M1_T2_M1 = 1;
constexpr int M1_T2_EN = 4;
constexpr int M1_T2_F1 = 5;
constexpr int M1_T2_NK = 8;
// vet_sc: the extrapolated components of vet_cell, M1_VET_CHI .. M1_VET_D11+5
constexpr int M1_T2_NVET = 10;
// the restart block of the stored slope: int32 have, int32 nch, Real pred2_dt,
// Real dt_prev, Real vprev; then nch cell slabs (K1, ipred2, vet_prev) behind the
// predictor slabs.  Absent when no slope is stored: old files stay byte-identical.
constexpr char kM1Time2RstMagic[8] = {'M', '1', 'T', 'I', 'M', 'E', '2', 'A'};
// <rad_m1>/implicit_one_pass (tests_m1/runs_4a_accel): a marked header behind the
// hesdirk2 one -- 9 Reals: the last two measured contractions and the solve counters of
// the three kinds of solve.  Written only when implicit_one_pass > 0.
constexpr char kM1OnePassRstMagic[8] = {'M', '1', 'O', 'N', 'E', 'P', '0', '1'};
// implicit_mr_every > 1 or vet_sc_every > 1 (rad_m1_mr.cpp): the multi-rate window and
// tensor-cadence state, behind the one-pass header: 10 Reals, no slabs
constexpr char kM1MRRstMagic[8] = {'M', '1', 'M', 'R', 'W', 'I', 'N', '1'};

// <rad_m1>/transport
constexpr int M1_TRANSPORT_EXPLICIT   = 0;   // stages 1-2: PD-ARS, sub-cycled
constexpr int M1_TRANSPORT_IMPLICIT_X1 = 1;  // stage 3a: backward Euler on x1 columns
constexpr int M1_TRANSPORT_IMPLICIT = 2;     // stage 3b phase B: the same backward-Euler
                                             // solve with the x2/x3 couplings added, the
                                             // x1 direction still solved exactly per
                                             // column and the transverse ones LAGGED
                                             // (line Jacobi) or wrapped in a Krylov
                                             // iteration.

// <rad_m1>/implicit_solver: how the full 7-point system is solved (transport = implicit)
constexpr int M1_ISOLV_LINE_JACOBI = 0;  // outer line-Jacobi: the x1 tridiagonal
                                         // system is solved exactly, the x2/x3
                                         // couplings are taken
                                         // from the PREVIOUS pass.  Their DIAGONAL part
                                         // stays on the matrix diagonal, so the full
                                         // 7-point M-matrix structure (and E' > 0) is
                                         // preserved at any dt.
constexpr int M1_ISOLV_BICGSTAB = 1;     // milestone 3b phase C: matrix-free BiCGStab on
                                         // the SAME frozen 7-point system, RIGHT-
                                         // preconditioned by the exact x1 line solve
                                         // (the tridiagonal part plus the full diagonal,
                                         // with the gather partition when blocks are
                                         // stacked along x1).  Its fixed point is the
                                         // fixed point line Jacobi converges to, so the
                                         // two solvers answer the same question; only
                                         // the number of passes differs.

// <rad_m1>/implicit_offdiag: what is done with the OFF-DIAGONAL Eddington terms
// sum_{e != d} d_e P_de of the face-flux equations (milestone 3b phase D).  They are the
// only part of the operator that is not in the 7-point stencil, and in an optically THIN
// cell they are multiplied by (c dt/dx)^2: lagging them is then a fixed-point iteration
// with no contraction and the outer loop diverges (the seeded 2-D He slab of phase C).
constexpr int M1_OD_LAGGED   = 0;  // phase B/C: evaluated at the previous Picard iterate
                                   // and put on the right-hand side of the face equation.
constexpr int M1_OD_OPERATOR = 1;  // phase D: with the closure (chi, n, hence D_ab)
                                   // FROZEN inside a Picard pass, d_e(D_de E) is LINEAR
                                   // in E', so it belongs in the matrix.  The right-hand
                                   // side is then built WITHOUT it and the matrix-free
                                   // operator application adds it back -- a 9-point
                                   // stencil in 2-D, 19-point in 3-D, through the halo
                                   // the Krylov vector is exchanged in anyway, so at no
                                   // extra communication.  The system is no longer an
                                   // M-matrix (the cross-derivative coefficients have
                                   // mixed signs) and E' > 0 is no longer guaranteed by
                                   // construction, so min E is monitored and a step that
                                   // produces a non-positive cell falls back to `none`.
                                   // Requires implicit_solver = bicgstab: the line-Jacobi
                                   // preconditioner alone cannot carry it.
constexpr int M1_OD_NONE     = 2;  // drop them: P is taken DIAGONAL in the grid frame.
                                   // Exact wherever the flux is along a grid axis or the
                                   // field is isotropic (chi = 1/3), wrong for an oblique
                                   // beam; the size of that error is gate G-oblique.

// <rad_m1>/implicit_trans_limit: the REALIZABILITY limiter of the transverse (x2/x3)
// face fluxes.  The face-eliminated flux F_f' = theta_f [...] has no free-streaming
// bound: in steady state it is the unlimited diffusive flux F = -c grad P/(rho kappa),
// which is super-luminal wherever rho kappa is tiny (the optically thin top of the 2-D
// He slab, where |F_2|/(c E) saturates at the post-solve clip of 1).
constexpr int M1_TLIM_NONE = 0;  // no limiter; the arithmetic of phase D, bit for bit.
constexpr int M1_TLIM_LP   = 1;  // add a lagged "limiter opacity" klim_f = |G_f|/(phi_f
                                 // E_f) to the face transport opacity, with G_f the face
                                 // pressure divergence the flux kernel already forms and
                                 // phi_f = fmax sqrt(max(1 - f1_f^2, 0.01)).  In steady
                                 // state |F_f| = c|G|/(rho kappa + |G|/(phi E)) <= c phi
                                 // E, and where R = |G|/(rho kappa E) << 1 the change is
                                 // O(R).  theta_f only SHRINKS, so the 7-point M-matrix
                                 // property of phase B is untouched.

// x1 boundary conditions of the implicit solve, in FACE-FLUX form (design sect. 3).
// The face value used is the TOTAL normal flux F at the boundary face.
constexpr int M1_IBC_MARSHAK  = 0;   // free surface, F_f = +- c*marshak_q*E_boundary
constexpr int M1_IBC_FLUX     = 1;   // imposed flux, F_f = implicit_flux_x1min/max
constexpr int M1_IBC_REFLECT  = 2;   // F_f = 0
constexpr int M1_IBC_PERIODIC = 3;   // cyclic tridiagonal (Sherman-Morrison)
constexpr int M1_IBC_EFIX     = 4;   // Dirichlet: the boundary CELL keeps the E it has
                                     // at the start of every step, i.e. it is frozen at
                                     // its initial value.  This is the face-flux form of
                                     // the fixed-E ghost the explicit tests T3b and T6
                                     // use, and the only boundary that anchors the LEVEL
                                     // of E in a pure-scattering column (imposed flux on
                                     // both ends leaves the operator singular).

// <rad_m1>/implicit_accel: how the OUTER (Picard) fixed-point iteration on the lagged
// closure is accelerated (milestone 3e).  The Picard map of one pass is
//   x = (E, F_1, F_2, F_3)  ->  G(x) = the state the pass leaves behind,
// the closure (chi, n) being rebuilt from x at the top of the pass.  In the optically
// thin top of the seeded 2-D He slab the map has gain >> 1 (the transverse pressure
// divergence acts as an advection at c dt/dx ~ 7e3 treated explicitly), so plain Picard
// -- and plain under-relaxation -- diverge.
constexpr int M1_IACC_NONE     = 0;  // the bare Picard map; bitwise the pre-3e code.
constexpr int M1_IACC_ANDERSON = 1;  // Anderson acceleration (Walker & Ni 2011): the next
                                     // iterate is the beta-mixed G over the affine span
                                     // of the last m iterates that minimises the
                                     // fixed-point residual g = G(x) - x in the 2-norm.
constexpr Real M1_AND_REG = 1.0e-10;  // the Tikhonov weight of the normal equations,
                                      // relative to tr(dG^T dG)/m
constexpr int M1_AND_MMAX = 10;      // the hard cap on implicit_anderson_m: the m x m
                                     // normal equations are solved on the HOST and the
                                     // per-column dot products ride one GlobalSum of
                                     // m + 1 <= NREDUCTION_VARIABLES entries.

// <rad_m1>/implicit_flux: which spatial form the implicit E-flux takes at a face
// (milestone 3a2, LIMIT 1 of the 3a findings).
constexpr int M1_IFLUX_CENTRAL = 0;  // 3a: the face-eliminated central (diffusion) form
                                     // at EVERY optical depth.  Exact in the thick limit,
                                     // heavily damped in the thin one (a non-upwinded
                                     // backward-Euler discretisation of the wave system).
constexpr int M1_IFLUX_APHLL   = 1;  // the same asymptotic-preserving blend the EXPLICIT
                                     // scheme uses, linearised in E':
                                     //   G_f = A_up + alpha F_HLL + (1-alpha) F_diff
                                     // with F_HLL built from the LAGGED reduced flux and
                                     // the closed-form M1 wave speeds.  alpha -> 1 in the
                                     // thin limit (upwind transport is recovered) and
                                     // alpha -> 0 in the thick one (3a is recovered).
                                     // The HLL DISSIPATION carries alpha^2 here, not
                                     // alpha, or the physical diffusion is counted twice
                                     // (see rad_m1_implicit.cpp).
constexpr int M1_IFLUX_BERTHON = 2;  // G_f = A_up + alpha F_HLL, with NO F_diff at all:
                                     // Berthon's asymptotic-preserving flux, which is
                                     // what the explicit scheme uses with piecewise-
                                     // constant states.  alpha is built so that
                                     // alpha (F_adv + F_dis) IS the physical diffusion at
                                     // every tau, so no second diffusive flux may be
                                     // added; the price is that the face diffusivity is
                                     // the arithmetic tau_face and not the exact harmonic
                                     // mean the face-eliminated form carries.
constexpr int M1_IFLUX_BLEND   = 3;  // milestone 3c: a SMOOTH per-face convex blend of
                                     // the two forms above,
                                     //   F_f = (1 - w_f) F_central + w_f F_berthon,
                                     // with w_f in [0,1] built from LAGGED face
                                     // quantities (the previous Picard iterate), so the
                                     // row stays linear.  Both forms contribute to the
                                     // SAME two unknowns with the sign structure of the
                                     // column-wise M-matrix argument of sect. 7, and a
                                     // convex combination of two column-sum-preserving
                                     // contributions preserves the column sums, so
                                     // E' > 0 survives at any dt.  w_f = 1 is bitwise
                                     // `berthon` and w_f = 0 is bitwise `central`.

// <rad_m1>/implicit_blend: how the per-face weight w_f of M1_IFLUX_BLEND is built.
constexpr int M1_IBLEND_TAU  = 0;    // w = exp(-(tau_face/tau0)^2): upwind where the face
                                     // is optically THIN, central where it is thick
                                     // (Jiang 2021 uses 1 - exp(-tau^2) for the opposite
                                     // sense).  <rad_m1>/implicit_blend_tau0.
constexpr int M1_IBLEND_F    = 1;    // w = smoothstep((f_face - f_lo)/(f_hi - f_lo)) on
                                     // the LAGGED comoving reduced flux at the face: a
                                     // DIFFUSE field (f -> 1/2, which is what a grey
                                     // surface with marshak_q = 1/2 has) stays central, a
                                     // beam or a front (f -> 1) goes upwind.
constexpr int M1_IBLEND_TAUF = 2;    // the product of the two: upwind only where the face
                                     // is thin AND the field is beamed.

// <rad_m1>/implicit_blend_fmode: which of the two cells' lagged reduced fluxes sets
// f_face.
constexpr int M1_IBFM_MAX  = 0;
constexpr int M1_IBFM_MEAN = 1;

// <rad_m1>/implicit_blend_mode: WHAT is blended.
constexpr int M1_IBMODE_FLUX   = 0;  // the whole face flux (the convex blend above)
constexpr int M1_IBMODE_DISSIP = 1;  // the central flux in FULL plus w_f times the HLL
                                     // DISSIPATION alone: F = F_central + w alpha dk
                                     // (E_L - E_R).  The extra term contributes +w dk to
                                     // the diagonal of the left cell and -w dk to the
                                     // upper off-diagonal (and the mirror image to the
                                     // right cell), so the column sums are untouched and
                                     // the M-matrix survives as well.

// <rad_m1>/implicit_recon: the states the HLL part of the implicit flux is built from.
constexpr int M1_IRECON_DC    = 0;   // piecewise constant; the matrix IS the operator
constexpr int M1_IRECON_PLMDC = 1;   // plm as a DEFERRED CORRECTION: the difference
                                     // between the plm and the dc flux, evaluated at the
                                     // PREVIOUS Picard iterate, goes into the right-hand
                                     // side, so the matrix stays the low-order M-matrix.

// <rad_m1>/implicit_enthalpy: the face value of the advective enthalpy flux A_d = a_d E
// (a_d = v_d + (v.D)_d) of every face of the implicit operator.  The MATRIX always
// carries the donor-cell form a_up E_up (upwind on the face velocity), which keeps the
// M-matrix; the other two add the difference between their face flux and the donor-cell
// one, evaluated at the PREVIOUS Picard iterate, to the right-hand side (a deferred
// correction, like implicit_recon = plm_dc).  At the Picard fixed point the scheme is the
// high-order one.  The donor-cell form is FIRST order even for a linear wave: with v = 0
// in the background its error (4/3) E0 (dx/2) sign(v) dv/dx is independent of the
// amplitude (tests_m1/runs_3r_radwave, tests_m1/runs_3s_space2).
constexpr int M1_IENTH_UPWIND  = 0;  // donor cell (old scheme; default to d0c59f7c)
constexpr int M1_IENTH_CENTRAL = 1;  // a_f (E_L + E_R)/2, a_f = (a_L + a_R)/2
constexpr int M1_IENTH_PLM     = 2;  // DEFAULT: a_f E_f, E_f the van Leer plm value from
                                     // the side upwind of a_f; central where the 4-cell
                                     // stencil is not available on both sides of the face

// <rad_m1>/implicit_vimp: the gas velocity v' of the enthalpy flux a(v') E' taken
// IMPLICIT through the radiative force of the same solve (tests_m1/runs_3v_vimplicit).
// Per face f of axis d the flux gains
//   E_f^k [ (1 + D_dd) dv_d^k + sum_{e != d} D_de dv_e^k ]_f
//     + E_f^k (1 + D_dd)_f [ (P dv_d)(E') - (P dv_d)(E^k) ]_f,
// dv^k the velocity change of the ITERATE's face fluxes (the write-back rule) and P the
// Jacobian of dv_d in E' (the face-normal eliminated flux, od/g0 terms lagged).  At the
// Picard fixed point the flux is E_f a_f(v') whatever the Jacobian; the operator couples
// cells two apart along every axis.  The components below are APPENDED to iw (offset
// RadiationM1::iw_vimp); the first M1_NVIMP_X are exchanged once per pass.
constexpr int M1_NVIMP_X  = 12;   // P1m,P10,P1p, P2m,P20,P2p, P3m,P30,P3p, DV1,DV2,DV3
constexpr int M1_IV_DV    = 9;
constexpr int M1_IV_X1M2  = 12;   // operator coefficients not in the 7-point row:
constexpr int M1_IV_X1P2  = 13;   // x1 +-2; x2 -2,-1,+1,+2; x3 -2,-1,+1,+2
constexpr int M1_IV_X2M2  = 14;
constexpr int M1_IV_X3M2  = 18;
constexpr int M1_IV_JD    = 22;   // added to the row: diagonal, x1 -1, x1 +1
constexpr int M1_IV_J1M   = 23;
constexpr int M1_IV_J1P   = 24;
constexpr int M1_IV_JRHS  = 25;   // added to the right-hand side
constexpr int M1_NIW_VIMP = 26;
// time_scheme = hesdirk2 (time2_enth_vel = start): 3 more components behind the block,
// DA_d = a_d(old vector) - a_d(stage start) = dv_d + (D dv)_d, dv = t2inc momentum/rho
constexpr int M1_IV_DA    = 26;
constexpr int M1_NIW_VIMP_T2 = 3;

// <rad_m1>/implicit_partition: how a column that spans several MeshBlocks (and ranks) is
// solved (milestone 3a2, LIMIT 4).
constexpr int M1_IPART_NONE   = 0;   // one MeshBlock per column; fatal otherwise (3a)
constexpr int M1_IPART_GATHER = 1;   // the column's tridiagonal rows are GATHERED onto
                                     // the rank that owns its lowest block, solved by the
                                     // very same serial Thomas sweep, and scattered back.
                                     // Bitwise independent of the partition by
                                     // construction: the arithmetic order is identical.

// components of the implicit FACE work array RadiationM1::ifw (m,n,k,j,i), i the x1 FACE
// index (face i is the one between cells i-1 and i), filled once per Picard iteration by
// the implicit_flux = ap_hll path and identically zero under `central`.
constexpr int M1_IFW_AL  = 0;   // alpha, the asymptotic-preserving weight of F_HLL
constexpr int M1_IFW_HCL = 1;   // alpha * (coefficient of E'_L in F_HLL), >= 0
constexpr int M1_IFW_HCR = 2;   // alpha * (coefficient of E'_R in F_HLL), <= 0
constexpr int M1_IFW_DG  = 3;   // the plm deferred correction, alpha (F^plm - F^dc)
constexpr int M1_NIFW = 4;

// components of the implicit work array RadiationM1::iw (m,n,k,j,i).  All are per-cell
// and live only for the duration of one solve, except EN/EGN which carry the
// start-of-step state through the Picard loop.
constexpr int M1_IW_EN   = 0;   // E^n, the radiation energy at the start of the step
constexpr int M1_IW_EP   = 1;   // the current Picard iterate of E'
constexpr int M1_IW_TP   = 2;   // the current Picard iterate of T'
constexpr int M1_IW_EGN  = 3;   // rho e^n, the gas INTERNAL energy density at the start
constexpr int M1_IW_WCHI = 4;   // w = P_11/E of the lagged closure (= chi in 1-D)
constexpr int M1_IW_ADV  = 5;   // a = v1 (1 + w), so that A = v E + v.P = a E
constexpr int M1_IW_DE0  = 6;   // E0 - E, the O(beta) comoving correction (explicit)
constexpr int M1_IW_G0   = 7;   // g0 = rho(kappa_E E0 - kappa_P a T^4), lagged
constexpr int M1_IW_V1   = 8;   // the x1 gas velocity (0 without <hydro>)
constexpr int M1_IW_SRCB = 9;   // the emission/absorption addition to the DIAGONAL
constexpr int M1_IW_SRCR = 10;  // the emission/absorption addition to the RHS
constexpr int M1_IW_TA   = 11;  // tridiagonal lower diagonal
constexpr int M1_IW_TB   = 12;  // tridiagonal diagonal
constexpr int M1_IW_TC   = 13;  // tridiagonal upper diagonal
constexpr int M1_IW_TR   = 14;  // tridiagonal right-hand side
constexpr int M1_IW_S1   = 15;  // Thomas scratch (c')
constexpr int M1_IW_S2   = 16;  // Thomas solution
constexpr int M1_IW_S3   = 17;  // second Thomas solution (cyclic / Sherman-Morrison)
constexpr int M1_IW_RES  = 18;  // per-cell Picard residual
constexpr int M1_IW_F1   = 19;  // the derived cell-centred x1 flux of the iterate
constexpr int M1_IW_RF0  = 20;  // the COMOVING reduced flux F0_cell/(c E) of the iterate,
                                // clipped to [-1,1]: what the HLL part of the ap_hll
                                // implicit flux lags (the closure chi keeps using the LAB
                                // reduced flux, as the explicit scheme does)
constexpr int M1_IW_KT   = 21;  // a verbatim copy of opac(M1_OP_T) = rho (kappa_F +
                                // kappa_s).  It exists only so that EVERY quantity the
                                // assembly and the face update read at a neighbouring
                                // cell lives in ONE array, which is what the x1 ghost
                                // halo of the partitioned solve exchanges (milestone 3b,
                                // LIMIT 4).  Single-block runs are bitwise unaffected:
                                // it is a copy, read where opac was read before.
// ---- the components above are ALL the x1-only solve (transport = implicit_x1) needs.
constexpr int M1_NIW_X1 = 22;
// ---- MILESTONE 3b phase B, transport = implicit: the transverse (x2/x3) couplings.
constexpr int M1_IW_V2   = 22;  // the x2 gas velocity (0 without <hydro>)
constexpr int M1_IW_V3   = 23;  // the x3 gas velocity
constexpr int M1_IW_N1   = 24;  // the LAGGED flux direction n_1 = F_1/|F| (unit vector)
constexpr int M1_IW_N2   = 25;  // n_2
constexpr int M1_IW_N3   = 26;  // n_3
constexpr int M1_IW_A2   = 27;  // a_2 = v_2 + (v.D)_2, so that A_2 = a_2 E
constexpr int M1_IW_A3   = 28;  // a_3 = v_3 + (v.D)_3
constexpr int M1_IW_TDIA = 29;  // the DIAGONAL part of the transverse couplings: what the
                                // line-Jacobi pass keeps on the matrix diagonal
constexpr int M1_IW_TRHS = 30;  // -(T(E^k) - TDIA E^k), the lagged transverse term that
                                // goes to the right-hand side
constexpr int M1_IW_LRES = 31;  // |U(E^{k+1}) - U(E^k)|, the TRUE residual of the full
                                // 7-point linear system (see the note in ImplicitSolve)
constexpr int M1_IW_F2   = 32;  // the derived cell-centred x2 flux of the iterate
constexpr int M1_IW_F3   = 33;  // ...and the x3 one; both drive the lagged closure
constexpr int M1_NIW = 34;

//----------------------------------------------------------------------------------------
//! \fn M1AccComp
//! \brief MILESTONE 3e: the iw component that carries component c of the Anderson
//! fixed-point vector.  c = 0 is E and c = 1..3 the cell-centred fluxes F_1, F_2, F_3 --
//! that is, the state the top of a Picard pass rebuilds the closure (chi, n) from, so
//! that one pass IS the map x -> G(x) whose fixed point is being sought.  The x1-only
//! solve has no F_2/F_3 slot in iw and uses c = 0,1 alone.

KOKKOS_INLINE_FUNCTION
int M1AccComp(const int c) {
  switch (c) {
    case 0: return M1_IW_EP;
    case 1: return M1_IW_F1;
    case 2: return M1_IW_F2;
    default: return M1_IW_F3;
  }
}

// ---- MILESTONE 3b phase C, implicit_solver = bicgstab.  The four (six in 3-D)
// TRANSVERSE
// OFF-DIAGONAL coefficients of the frozen 7-point row, and the ten Krylov vectors.  They
// are allocated only when the solver is bicgstab, so line_jacobi keeps the array it had.
//
// The row of cell c is
//   TA E_{i-1} + TB E_i + TC E_{i+1} + CJM E_{j-1} + CJP E_{j+1} + CKM E_{k-1}
//     + CKP E_{k+1}  =  b_c ,
// where TB ALREADY carries the transverse diagonal M1_IW_TDIA and b_c = TR + sum_nb C_nb
// E^k_nb, i.e. the right-hand side line Jacobi would use PLUS the lagged off-diagonal
// term it moved there.  Subtracting it back is what makes the two solvers solve the very
// same linear system: one line-Jacobi pass is exactly x <- M^{-1}(b - sum C x).
constexpr int M1_IW_CJM  = 34;  // coefficient of E_{j-1} (0 at a physical x2 face)
constexpr int M1_IW_CJP  = 35;  // coefficient of E_{j+1}
constexpr int M1_IW_CKM  = 36;  // coefficient of E_{k-1} (0 in 2-D)
constexpr int M1_IW_CKP  = 37;  // coefficient of E_{k+1}
constexpr int M1_IW_KB   = 38;  // b, the right-hand side of the frozen 7-point system
constexpr int M1_IW_KX   = 39;  // the BiCGStab iterate (starts at the Picard iterate)
constexpr int M1_IW_KR   = 40;  // the recursive residual r
constexpr int M1_IW_KRH  = 41;  // the shadow residual rhat (fixed between restarts)
constexpr int M1_IW_KP   = 42;  // the search direction p
constexpr int M1_IW_KV   = 43;  // v = A M^{-1} p
constexpr int M1_IW_KS   = 44;  // s = r - alpha v
constexpr int M1_IW_KTT  = 45;  // t = A M^{-1} s
constexpr int M1_IW_KY   = 46;  // y = M^{-1} p   (the preconditioned direction)
constexpr int M1_IW_KZ   = 47;  // z = M^{-1} s
constexpr int M1_NIW_K = 48;

// ---- MILESTONE 3g: the GAS-RADIATION energy coupling, <rad_m1>/implicit_gas_newton and
// <rad_m1>/implicit_eos_cache.  Both default to false, neither allocates anything then
// and neither is referenced, so every earlier configuration is bitwise unchanged.
//
// The five components below are APPENDED after whichever base component count the solve
// uses (M1_NIW_X1 / M1_NIW / M1_NIW_K), so the base indices above never move; the offset
// of the block is RadiationM1::iw_gas, a runtime int the kernels capture.
constexpr int M1_NIW_GAS = 5;
constexpr int M1_IWG_BK = 0;   // B_k = rho c_v + 4 c dt rho kappa_P a T_k^3, the
                               // derivative d/dT of the LOCAL gas energy equation
constexpr int M1_IWG_RK = 1;   // R_k = rho e^n - rho e(T_k) - c dt rho kappa_P a T_k^4
                               //       + c dt rho kappa_E de0, its residual at E' = 0
constexpr int M1_IWG_YR = 2;   // |R_k + c dt rho kappa_E E'| of the pass that produced
                               // the current T iterate: the nonlinear gas residual the
                               // NEXT pass checks the Newton step against (0 = the step
                               // was an exact root find, no check needed)
constexpr int M1_IWG_FB = 3;   // per-cell count of Newton fallbacks in this step
constexpr int M1_IWG_MS = 4;   // per-cell count of EOS-cache misses in this step

// The TRUST REGION of the Newton temperature update: a step larger than this fraction of
// T_k is refused and the bracketed root find runs instead.  y(T) = rho e(T) + A T^4 -
// target is increasing and (for a table whose c_v increases with T, which is what
// ionization does) convex, so Newton from either side converges monotonically; the guard
// is there for the cells where that is not true -- a recombination shoulder in c_v, or a
// first pass that moves T by a factor.  It is a heuristic, not a proof: the PROOF that
// the converged state solves the exact nonlinear equation is that at the fixed point
// dT -> 0 forces R_k + c dt rho kappa_E E' -> 0, which IS that equation, with e(T) and
// T^4 evaluated exactly (never linearised) at the final T.
constexpr Real M1_NEWT_TRUST = 0.5;

// BiCGStab breakdown thresholds: |rho| and |rhat.v| below these times the scale of the
// right-hand side mean the shadow residual has become orthogonal to the Krylov space.
constexpr Real M1_BCG_EPS = 1.0e-300;
// implicit_krylov_dev (rad_m1_launch.cpp): the size of its device scalar array
constexpr int M1_KD_SIZE = 20;

// the LAGGED quantities the transverse halo exchanges once per Picard pass.  Everything
// the x2/x3 face fluxes, the lagged off-diagonal Eddington terms and the x1 assembly read
// at a NEIGHBOURING cell is in this list, so the two blocks that share a face build it
// from bit-identical numbers.
// M1_IW_F1 is in the list ONLY for the transverse realizability limiter, which needs the
// face mean of the lagged x1 reduced flux and must build it from numbers both blocks of
// a shared face agree on; it is not read anywhere else, so carrying it changes no
// arithmetic under implicit_trans_limit = none.
//
// The list is ordered so that a PREFIX of it is a legal exchange on its own.  The first
// M1_NHALO_Q entries are the only ones a Picard pass can move once the closure is frozen
// (implicit_closure_lag = step): E and T' of the iterate feed G0, the pass re-derives the
// cell flux F1, and KT moves only under implicit_opac_update.  Everything after them --
// the closure (chi, n), the advective coefficients built from it and the gas velocities
// -- is recomputed by the pass from FROZEN inputs with the same expressions, so it comes
// out bit-identical and its ghost layer is still the one the previous exchange left.
// The last three (the velocities) never move inside a step at all.
constexpr int M1_NHALO_T = 14;
constexpr int M1_NHALO_Q = 4;
KOKKOS_INLINE_FUNCTION
int M1HaloCompT(const int n) {
  switch (n) {
    case 0: return M1_IW_EP;
    case 1: return M1_IW_G0;
    case 2: return M1_IW_F1;
    case 3: return M1_IW_KT;
    case 4: return M1_IW_WCHI;
    case 5: return M1_IW_N1;
    case 6: return M1_IW_N2;
    case 7: return M1_IW_N3;
    case 8: return M1_IW_ADV;
    case 9: return M1_IW_A2;
    case 10: return M1_IW_A3;
    case 11: return M1_IW_V1;
    case 12: return M1_IW_V2;
    default: return M1_IW_V3;
  }
}

//----------------------------------------------------------------------------------------
//! \fn M1EddDiag
//! \brief the DIAGONAL Eddington-tensor component D_dd = P_dd/E of direction d (0,1,2)
//! from the lagged closure: D_ab = (1-chi)/2 delta_ab + (3 chi - 1)/2 n_a n_b.

KOKKOS_INLINE_FUNCTION
Real M1EddDiag(const Real chi, const Real nd) {
  return 0.5*(1.0 - chi) + 0.5*(3.0*chi - 1.0)*nd*nd;
}

//----------------------------------------------------------------------------------------
//! \fn M1EddOff
//! \brief the OFF-diagonal component D_ab = (3 chi - 1)/2 n_a n_b, a != b.

KOKKOS_INLINE_FUNCTION
Real M1EddOff(const Real chi, const Real na, const Real nb) {
  return 0.5*(3.0*chi - 1.0)*na*nb;
}

//----------------------------------------------------------------------------------------
//! \fn M1POff
//! \brief the OFF-diagonal radiation pressure P_ab = D_ab E at one cell, from the LAGGED
//! closure of the previous Picard pass.  a, b are 0-based directions.  `ec` names the
//! component of the work array the energy is taken from: M1_IW_EP for the lagged form,
//! and a KRYLOV VECTOR when the term is applied as part of the linear operator
//! (implicit_offdiag = operator), where the closure is frozen and the term is linear.

//!
//! FULL TENSOR (<rad_m1>/closure = vet_sc, vet_tensor = full): with `full` the component
//! is read from the guarded D_ab = K_ab/J of the formal solution, RadiationM1::vet_cell
//! (vd), instead of being built from the uniaxial (chi, n).  With full = false vd is
//! never touched and the arithmetic is that of the uniaxial form, term for term.

template <class V>
KOKKOS_INLINE_FUNCTION
Real M1POff(const V &iw, const int m, const int a, const int b,
            const int k, const int j, const int i, const int ec,
            const V &vd, const bool full) {
  if (full) {return vd(m,M1_VET_D11+2+a+b,k,j,i)*iw(m,ec,k,j,i);}
  return M1EddOff(iw(m,M1_IW_WCHI,k,j,i), iw(m,M1_IW_N1+a,k,j,i),
                  iw(m,M1_IW_N1+b,k,j,i))*iw(m,ec,k,j,i);
}

//----------------------------------------------------------------------------------------
//! \fn M1DDiag
//! \brief the DIAGONAL Eddington component D_dd of direction d (0,1,2) at one cell of the
//! MULTI-D solve: M1EddDiag of the lagged uniaxial closure (WCHI, N1+d), or with `full`
//! the guarded D_dd of the formal solution (see M1POff).

template <class V>
KOKKOS_INLINE_FUNCTION
Real M1DDiag(const V &iw, const V &vd, const bool full, const int m, const int d,
             const int k, const int j, const int i) {
  if (full) {return vd(m,M1_VET_D11+d,k,j,i);}
  return M1EddDiag(iw(m,M1_IW_WCHI,k,j,i), iw(m,M1_IW_N1+d,k,j,i));
}

//----------------------------------------------------------------------------------------
//! \fn M1OffDiv
//! \brief sum_{e != d} d_e P_de at one cell: the OFF-DIAGONAL part of the divergence of
//! the radiation pressure that drives the face-normal flux of direction d.  It is fully
//! LAGGED (previous-pass E and closure) and enters the right-hand side of the face-flux
//! equation; the diagonal part d_d P_dd is what the matrix carries.
//!
//! Centred differences, with the index range [il,iu] x [jl,ju] x [kl,ku] the caller may
//! read -- one ghost layer wide wherever a neighbouring MeshBlock (or a periodic wrap)
//! has filled it, and the active range at a PHYSICAL boundary, where the difference
//! silently becomes one-sided (the divisor counts the cells actually used).

template <class V>
KOKKOS_INLINE_FUNCTION
Real M1OffDiv(const V &iw, const int m, const int d, const int k, const int j,
              const int i, const Real dx1, const Real dx2, const Real dx3,
              const bool thrd,
              const int il, const int iu, const int jl, const int ju,
              const int kl, const int ku, const int ec,
              const V &vd, const bool full) {
  Real s = 0.0;
  if (d != 0) {
    int ia = (i+1 <= iu) ? (i+1) : i;
    int ib = (i-1 >= il) ? (i-1) : i;
    if (ia != ib) {
      s += (M1POff(iw,m,d,0,k,j,ia,ec,vd,full)
            - M1POff(iw,m,d,0,k,j,ib,ec,vd,full))/((ia - ib)*dx1);
    }
  }
  if (d != 1) {
    int ja = (j+1 <= ju) ? (j+1) : j;
    int jb = (j-1 >= jl) ? (j-1) : j;
    if (ja != jb) {
      s += (M1POff(iw,m,d,1,k,ja,i,ec,vd,full)
            - M1POff(iw,m,d,1,k,jb,i,ec,vd,full))/((ja - jb)*dx2);
    }
  }
  if (thrd && d != 2) {
    int ka = (k+1 <= ku) ? (k+1) : k;
    int kb = (k-1 >= kl) ? (k-1) : k;
    if (ka != kb) {
      s += (M1POff(iw,m,d,2,ka,j,i,ec,vd,full)
            - M1POff(iw,m,d,2,kb,j,i,ec,vd,full))/((ka - kb)*dx3);
    }
  }
  return s;
}

//----------------------------------------------------------------------------------------
//! \fn M1SphDrr
//! \brief STAGE S2 (spherical-polar wedge, closure m1/minerbo/kershaw): the radial face
//! weight of one cell with the INTEGRATING FACTOR of design sect. 1.4 (branch
//! m1-curv-design).  With P = p I + q n n, p = (1-chi) E/2, q = (3 chi-1) E/2,
//!   (div P)_r = d_r p + (1/r^2) d_r(r^2 q n_r^2) - q (1 - n_r^2)/r + [tangential],
//! and the face gradient at r_f is
//!   (p_R - p_L)/dr_f + (r_R^2 Q_R E_R - r_L^2 Q_L E_L)/(r_f^2 dr_f),  Q = q n_r^2/E,
//! i.e. the Cartesian two-point difference of D E with D_rr replaced, per cell and per
//! face, by (1-chi)/2 + (3 chi-1)/2 n_r^2 (r_c/r_f)^2.  s2 = (r_c/r_f)^2.  Both parts are
//! >= 0, so the row stays an M-matrix; chi = 1, n = r gives r^2 E = const exactly.

KOKKOS_INLINE_FUNCTION
Real M1SphDrr(const Real chi, const Real n1, const Real s2) {
  return 0.5*(1.0 - chi) + 0.5*(3.0*chi - 1.0)*n1*n1*s2;
}

//----------------------------------------------------------------------------------------
//! \fn M1SphMarshakFaceE
//! \brief <rad_m1>/implicit_marshak_face = linear (m1-sp-order2, spherical-polar wedge,
//! tests_m1/runs_5o_sporder2): the E of a Marshak boundary FACE at r_f, from the end cell
//! (E ec at its centroid rc) and its interior neighbour (en at rn), by the linear
//! extrapolation of G = r^2 E to the face: second order, and exact for the free-
//! streaming state r^2 E = const.  E_lin = ca ec - cb en with
//!   alpha = (r_f - rc)/(rc - rn) > 0 (either end),  ca = (1 + alpha) (rc/r_f)^2,
//!   cb = alpha (rn/r_f)^2,  ca > cb >= 0.
//! The implicit row carries E_lin itself (ca on the diagonal, -cb on the neighbour: a
//! non-positive off-diagonal and a diagonally dominant row, so the M-matrix sign pattern
//! stays).  The returned value is E_lin LIMITED to [0, 2 ca ec / (1 + alpha)] (i.e. G_f
//! in [0, 2 G_end], the bound of the van Leer face value of implicit_enthalpy = plm);
//! the difference limited - linear at the lagged iterate is a deferred correction that
//! vanishes wherever the limiter is inactive (every smooth profile).

KOKKOS_INLINE_FUNCTION
void M1SphMarshakCoef(const Real rc, const Real rn, const Real rf, Real &ca, Real &cb) {
  const Real al = (rf - rc)/(rc - rn);
  ca = (1.0 + al)*SQR(rc/rf);
  cb = al*SQR(rn/rf);
}

KOKKOS_INLINE_FUNCTION
Real M1SphMarshakFaceE(const Real ec, const Real en, const Real rc, const Real rn,
                       const Real rf) {
  Real ca, cb;
  M1SphMarshakCoef(rc, rn, rf, ca, cb);
  const Real ecap = 2.0*SQR(rc/rf)*ec;
  return fmin(fmax(ca*ec - cb*en, 0.0), ecap);
}

//----------------------------------------------------------------------------------------
//! \fn M1SphCurv
//! \brief STAGE S2: the LAGGED part of (div P)_d at one cell centre of the
//! spherical-polar wedge (physical orthonormal components r, theta, phi), i.e. what
//! the face equation of direction d does not carry implicitly.  The implicit parts are
//! d_r of the integrating-factor form (M1SphDrr, radial faces) and (1/r) d_th P_thth,
//! (1/(r sin)) d_ph P_phph on the transverse faces (the S1 two-point differences over
//! dxface).
//! With P = p I + q n n (uniaxial closure, components from WCHI and N1..N3):
//!   d = 0:  - q (1 - n_r^2)/r
//!           + [od] (1/r) d_th P_rt + cot P_rt / r + (1/(r sin)) d_ph P_rp
//!   d = 1:  cot (P_tt - P_pp)/r
//!           + [od] (1/r^3) d_r (r^3 P_rt) + (1/(r sin)) d_ph P_tp
//!   d = 2:  [od] (1/r^3) d_r (r^3 P_rp) + (1/r) d_th P_tp + 2 cot P_tp / r
//! [od] = the off-diagonal (tangential-derivative) terms, present when `od` is set
//! (implicit_offdiag = lagged; the sp default is none); the curvature of the diagonal
//! components is always on.
//! Centred differences over the coordinate distances of the neighbours, one-sided at a
//! physical boundary exactly as M1OffDiv ([il,iu] x [jl,ju] x [kl,ku] the readable
//! range).  E is read from component `ec` (M1_IW_EP: the lagged form).  The wedge stays
//! clear of the poles, so cot and 1/sin are finite at every cell centre.

template <class V>
KOKKOS_INLINE_FUNCTION
Real M1SphPab(const V &iw, const int m, const int a, const int b, const int k,
              const int j, const int i, const int ec) {
  return 0.5*(3.0*iw(m,M1_IW_WCHI,k,j,i) - 1.0)*iw(m,M1_IW_N1+a,k,j,i)
         *iw(m,M1_IW_N1+b,k,j,i)*iw(m,ec,k,j,i);
}

template <class V, class A>
KOKKOS_INLINE_FUNCTION
Real M1SphCurv(const V &iw, const A &x1v, const A &x2v, const A &x3v, const int m,
               const int d, const int k, const int j, const int i, const bool od,
               const bool thrd, const int il, const int iu, const int jl, const int ju,
               const int kl, const int ku, const int ec) {
  const Real r = x1v(m,i);
  const Real sn = sin(x2v(m,j));
  const Real ct = cos(x2v(m,j))/sn;
  const Real e = iw(m,ec,k,j,i);
  const Real chi = iw(m,M1_IW_WCHI,k,j,i);
  const Real n1 = iw(m,M1_IW_N1,k,j,i);
  const Real n2 = iw(m,M1_IW_N2,k,j,i);
  const Real n3 = iw(m,M1_IW_N3,k,j,i);
  const Real q = 0.5*(3.0*chi - 1.0)*e;
  Real s = 0.0;
  if (d == 0) {
    s = -q*(1.0 - n1*n1)/r;
  } else if (d == 1) {
    s = ct*q*(n2*n2 - n3*n3)/r;
  }
  if (!od) {return s;}
  const int ia = (i+1 <= iu) ? (i+1) : i;
  const int ib = (i-1 >= il) ? (i-1) : i;
  const int ja = (j+1 <= ju) ? (j+1) : j;
  const int jb = (j-1 >= jl) ? (j-1) : j;
  const int ka = (thrd && k+1 <= ku) ? (k+1) : k;
  const int kb = (thrd && k-1 >= kl) ? (k-1) : k;
  if (d == 0) {
    if (ja != jb) {
      s += (M1SphPab(iw,m,0,1,k,ja,i,ec) - M1SphPab(iw,m,0,1,k,jb,i,ec))
           /(r*(x2v(m,ja) - x2v(m,jb)));
    }
    s += ct*q*n1*n2/r;
    if (ka != kb) {
      s += (M1SphPab(iw,m,0,2,ka,j,i,ec) - M1SphPab(iw,m,0,2,kb,j,i,ec))
           /(r*sn*(x3v(m,ka) - x3v(m,kb)));
    }
  } else if (d == 1) {
    if (ia != ib) {
      const Real ra = x1v(m,ia), rb = x1v(m,ib);
      s += (ra*ra*ra*M1SphPab(iw,m,0,1,k,j,ia,ec) - rb*rb*rb*M1SphPab(iw,m,0,1,k,j,ib,ec))
           /(r*r*r*(ra - rb));
    }
    if (ka != kb) {
      s += (M1SphPab(iw,m,1,2,ka,j,i,ec) - M1SphPab(iw,m,1,2,kb,j,i,ec))
           /(r*sn*(x3v(m,ka) - x3v(m,kb)));
    }
  } else {
    if (ia != ib) {
      const Real ra = x1v(m,ia), rb = x1v(m,ib);
      s += (ra*ra*ra*M1SphPab(iw,m,0,2,k,j,ia,ec) - rb*rb*rb*M1SphPab(iw,m,0,2,k,j,ib,ec))
           /(r*r*r*(ra - rb));
    }
    if (ja != jb) {
      s += (M1SphPab(iw,m,1,2,k,ja,i,ec) - M1SphPab(iw,m,1,2,k,jb,i,ec))
           /(r*(x2v(m,ja) - x2v(m,jb)));
    }
    s += 2.0*ct*q*n2*n3/r;
  }
  return s;
}

// the six LAGGED quantities the x1 halo of the partitioned solve exchanges once per
// Picard iteration (halo "A"), in the order the pack kernel uses.  The seventh exchange
// (halo "B") carries M1_IW_EP alone, after the line solve has accepted the new iterate.
constexpr int M1_NHALO_A = 6;
KOKKOS_INLINE_FUNCTION
int M1HaloCompA(const int n) {
  switch (n) {
    case 0: return M1_IW_WCHI;
    case 1: return M1_IW_ADV;
    case 2: return M1_IW_G0;
    case 3: return M1_IW_V1;
    case 4: return M1_IW_RF0;
    default: return M1_IW_KT;
  }
}

//----------------------------------------------------------------------------------------
//! \fn M1BlendWeight
//! \brief the per-face blend weight w_f in [0,1] of implicit_flux = blend, from the
//! LAGGED face optical depth and the LAGGED comoving reduced fluxes of the two cells.
//! Every quantity it reads lives in the iw work array, so the x1 halo of the partitioned
//! solve already carries what a face on a block boundary needs (milestone 3b, LIMIT 4).

KOKKOS_INLINE_FUNCTION
Real M1BlendWeight(const int kind, const int fmode, const Real tauf, const Real tau0,
                   const Real rfl, const Real rfr, const Real flo, const Real fhi) {
  Real wt = 1.0;
  if (kind != M1_IBLEND_F) {
    Real x = tauf/tau0;
    wt = exp(-x*x);
  }
  Real wf = 1.0;
  if (kind != M1_IBLEND_TAU) {
    Real afl = fabs(rfl), afr = fabs(rfr);
    Real ff = (fmode == M1_IBFM_MEAN) ? (0.5*(afl + afr)) : fmax(afl, afr);
    Real s = (ff - flo)/fmax(fhi - flo, 1.0e-300);
    if (s < 0.0) {s = 0.0;}
    if (s > 1.0) {s = 1.0;}
    wf = s*s*(3.0 - 2.0*s);
  }
  return wt*wf;
}

//----------------------------------------------------------------------------------------
//! MILESTONE 3g, the FROZEN-DENSITY EOS CACHE, <rad_m1>/implicit_eos_cache.
//!
//! The density is frozen over the radiation step, so every e(rho,T) the gas solve asks
//! for during the Picard loop lies on ONE line of the table.  The tabulated surface is a
//! bicubic Hermite patch in (x,y) = (log10 rho, log10 T),
//!
//!   f(u,v) = sum_{a,b} c_ab hu_a(u) hv_b(v),
//!
//! with u frozen; collapsing the x direction ONCE per step,
//!
//!   C_b = sum_a hu_a(u) c_ab   ->   f(v) = sum_b C_b hv_b(v),
//!
//! leaves a 1-D cubic Hermite in ln T whose four coefficients are stored per cell.  It is
//! the table's OWN spline, on the table's OWN temperature nodes -- not a fit -- so a
//! cached evaluation differs from the direct one only by the summation order, i.e. by
//! floating-point round-off (measured in ImplicitReport as `eos_cache: max |de|/e`).
//!
//! One set of coefficients is valid on one table cell in T.  The cache therefore holds
//! 2*nt+1 consecutive cells centred on the cell of T^n (<rad_m1>/implicit_eos_cache_nt),
//! and a temperature that leaves that window falls back to the real table for that ONE
//! evaluation and is counted as a miss.
//!
//! NOT cached, and falling back to the table unconditionally: a density off the table, a
//! sub-table temperature (the <block>/eos_floor_consistent continuation), and a tapered
//! radiation term (<block>/eos_rad_taper), whose weight depends on BOTH x and y.
//!
//! Layout of the per-cell record ec(m,n,k,j,i): n = 0 is iy0, the lowest table T node of
//! the window, or -1 when the cell has no cache; n = 1 + 4c + b is C_b of cell c.

constexpr int M1_EC_NTMAX = 8;   // cap on implicit_eos_cache_nt

KOKKOS_INLINE_FUNCTION
int M1EosCacheNComp(const int nt) { return 4*(2*nt + 1) + 1; }

//----------------------------------------------------------------------------------------
//! \fn M1EosCacheBuild
//! \brief collapse the density direction of the tabulated energy surface at (rho, T^n)
//! and store the 1-D Hermite coefficients of the window.  Called once per cell per step.

template <class EosT, class V>
KOKKOS_INLINE_FUNCTION
void M1EosCacheBuild(const EosT &eos, const V &ec, const int m, const int k,
                     const int j, const int i, const int nt,
                     const Real dd, const Real tt) {
  ec(m,0,k,j,i) = -1.0;
  const auto &tb = eos.tbl;
  if (!tb.active || tb.rad_taper) {return;}
  const Real rho = dd*eos.dens_cgs;
  const Real tk = tt*eos.temp_cgs;
  if (!(rho > 0.0) || !(tk > 0.0)) {return;}
  const Real gx = (log10(rho) - tb.xmin)*tb.dxi;
  const int ix = static_cast<int>(floor(gx));
  if (ix < 0 || ix > tb.nx-2) {return;}
  const Real u = gx - static_cast<Real>(ix);
  if (u < 0.0 || u > 1.0) {return;}
  const int nc = 2*nt + 1;
  if (tb.ny-1 < nc) {return;}
  const Real gy = (log10(tk) - tb.ymin)*tb.dyi;
  int iy0 = static_cast<int>(floor(gy)) - nt;
  if (iy0 < 0) {iy0 = 0;}
  if (iy0 > tb.ny-1-nc) {iy0 = tb.ny-1-nc;}
  const Real u2 = u*u, u3 = u2*u;
  const Real hu[4] = {2.0*u3 - 3.0*u2 + 1.0, u3 - 2.0*u2 + u,
                      -2.0*u3 + 3.0*u2, u3 - u2};
  for (int c=0; c<nc; ++c) {
    const int iy = iy0 + c;
    for (int b=0; b<4; ++b) {
      // the Hermite coefficient matrix of HermitePatch(): row r -> node ix + (r>>1) and
      // the dx-scaled slope surface when r is odd, column b -> node iy + (b>>1) and the
      // dy-scaled slope surface when b is odd.
      const int j0 = iy + (b >> 1);
      const Real sy = (b & 1) ? tb.dy : 1.0;
      Real s = 0.0;
      for (int a=0; a<4; ++a) {
        const int i0 = ix + (a >> 1);
        const Real sx = (a & 1) ? tb.dx : 1.0;
        s += hu[a]*tb.tbl(j0, i0, ITE + (a & 1) + 2*(b & 1))*sx*sy;
      }
      ec(m,1+4*c+b,k,j,i) = s;
    }
  }
  ec(m,0,k,j,i) = static_cast<Real>(iy0);
}

//----------------------------------------------------------------------------------------
//! \fn M1EosCacheEval
//! \brief e(T) and c_v = de/dT at the frozen density from the cached coefficients, in
//! CODE units and with exactly the arithmetic EOS_Data::ThermoAt() would apply to the
//! same interpolated surface.  Returns false (leaving ee, cv untouched) on a miss.

template <class EosT, class V>
KOKKOS_INLINE_FUNCTION
bool M1EosCacheEval(const EosT &eos, const V &ec, const int m, const int k, const int j,
                    const int i, const int nt, const Real dd, const Real tt,
                    Real &ee, Real &cv) {
  const Real fy0 = ec(m,0,k,j,i);
  if (fy0 < 0.0) {return false;}
  const auto &tb = eos.tbl;
  const Real tk = tt*eos.temp_cgs;
  if (!(tk > 0.0)) {return false;}
  const Real gy = (log10(tk) - tb.ymin)*tb.dyi;
  const int iq = static_cast<int>(floor(gy));
  const int c = iq - static_cast<int>(fy0);
  if (c < 0 || c > 2*nt) {return false;}
  const Real v = gy - static_cast<Real>(iq);
  if (v < 0.0 || v > 1.0) {return false;}
  const Real v2 = v*v, v3 = v2*v;
  const Real hv[4] = {2.0*v3 - 3.0*v2 + 1.0, v3 - 2.0*v2 + v,
                      -2.0*v3 + 3.0*v2, v3 - v2};
  const Real dv[4] = {6.0*v2 - 6.0*v, 3.0*v2 - 4.0*v + 1.0,
                      -6.0*v2 + 6.0*v, 3.0*v2 - 2.0*v};
  Real f = 0.0, fy = 0.0;
  for (int b=0; b<4; ++b) {
    const Real q = ec(m,1+4*c+b,k,j,i);
    f += q*hv[b];
    fy += q*dv[b];
  }
  fy *= tb.dyi;
  const Real rho = dd*eos.dens_cgs;
  const Real egas = rho*EOSTable::Pow10(f);
  Real e = egas;
  Real cvg = (egas/rho)*fy/tk;
  if (tb.radiation) {
    const Real erad = tb.arad*tk*tk*tk*tk;
    e += erad;
    cvg += 4.0*erad/(rho*tk);
  }
  ee = e/eos.pres_cgs;
  cv = cvg*eos.dens_cgs*eos.temp_cgs/eos.pres_cgs;
  return true;
}

//----------------------------------------------------------------------------------------
//! \struct M1EosDirect
//! \brief the (e, c_v) accessor the gas solve uses when the cache is OFF: one
//! EOS_Data::ThermoAt() call, in the order the pre-3g code made it.

template <class EosT>
struct M1EosDirect {
  EosT eos;
  KOKKOS_INLINE_FUNCTION
  void operator()(const Real dd, const Real t, Real &ee, Real &cv) const {
    Real pp, cr, ct;
    eos.ThermoAt(dd, t, ee, pp, cr, ct, cv);
  }
};

//----------------------------------------------------------------------------------------
//! \struct M1EosCached
//! \brief the same accessor served from the per-cell cache, with the table as fallback.
//! `nm` points at a thread-local miss counter (may be null).

template <class EosT, class V>
struct M1EosCached {
  EosT eos;
  V ec;
  int m, k, j, i, nt;
  Real *nm;
  KOKKOS_INLINE_FUNCTION
  void operator()(const Real dd, const Real t, Real &ee, Real &cv) const {
    if (M1EosCacheEval(eos, ec, m, k, j, i, nt, dd, t, ee, cv)) {return;}
    if (nm != nullptr) {*nm += 1.0;}
    Real pp, cr, ct;
    eos.ThermoAt(dd, t, ee, pp, cr, ct, cv);
  }
};

// safeguarded root find for T' inside the Picard loop
constexpr int  M1_IMPL_TMAXIT = 100;
constexpr Real M1_IMPL_TRTOL  = 1.0e-12;

//----------------------------------------------------------------------------------------
//! \fn M1ImplTemperature
//! \brief solve the LOCAL gas energy equation of the implicit transport scheme,
//!
//!   y(T) = rho e(rho,T) + c dt rho kappa_P a T^4 - [ rho e^n + c dt rho kappa_E E0' ]
//!        = 0,
//!
//! for T at fixed E' (E0' = E' + de0), by the same safeguarded Newton (`rtsafe`) that
//! the explicit coupling uses -- bracket by walking, then a Newton step accepted only
//! when it stays inside the bracket AND shrinks it at least as fast as bisection.  y is
//! strictly increasing in T, so the bracket is unique.
//!
//! This is NOT the explicit path's function (there E'(T) is eliminated locally and the
//! invariant is the bracket variable), so it cannot call rad_m1_coupling.cpp's copy;
//! the ALGORITHM is the same and the explicit arithmetic is left untouched.
//!
//! Returns the number of function evaluations; `ok` is false if no bracket was found.

template <class ThermoT>
KOKKOS_INLINE_FUNCTION
int M1ImplTemperatureT(const ThermoT &th, const Real dd, const Real tguess,
                       const Real egn, const Real cdtkp_a, const Real rhs_abs,
                       Real &tout, bool &ok) {
  // y(T) = ee(T) + cdtkp_a*T^4 - (egn + rhs_abs)
  Real target = egn + rhs_abs;
  Real tlo = tguess, thi = tguess;
  Real ylo, yhi;
  int nev = 0;
  {
    Real ee, cv;
    th(dd, tguess, ee, cv);
    Real t2 = tguess*tguess;
    Real y0 = ee + cdtkp_a*t2*t2 - target;
    ylo = y0;
    yhi = y0;
    ++nev;
  }
  ok = true;
  if (ylo > 0.0) {
    for (int it=0; it<80 && ylo > 0.0; ++it) {
      tlo *= 0.5;
      Real ee, cv;
      th(dd, tlo, ee, cv);
      Real t2 = tlo*tlo;
      ylo = ee + cdtkp_a*t2*t2 - target;
      ++nev;
    }
    ok = (ylo <= 0.0);
  } else if (yhi < 0.0) {
    for (int it=0; it<80 && yhi < 0.0; ++it) {
      thi *= 2.0;
      Real ee, cv;
      th(dd, thi, ee, cv);
      Real t2 = thi*thi;
      yhi = ee + cdtkp_a*t2*t2 - target;
      ++nev;
    }
    ok = (yhi >= 0.0);
  }
  if (!ok) {
    tout = tguess;
    return nev;
  }
  Real tp = 0.5*(tlo + thi);
  if (tlo == thi) {tp = tlo;}
  Real dxold = fabs(thi - tlo);
  Real dx = dxold;
  bool conv = (tlo == thi);
  for (int it=0; it<M1_IMPL_TMAXIT && !conv; ++it) {
    Real ee, cv;
    th(dd, tp, ee, cv);
    ++nev;
    Real t3 = tp*tp*tp;
    Real y = ee + cdtkp_a*tp*t3 - target;
    if (y > 0.0) {thi = tp;} else {tlo = tp;}
    Real dy = dd*cv + 4.0*cdtkp_a*t3;
    dxold = dx;
    Real tn = (dy > 0.0) ? (tp - y/dy) : tp;
    if (!(tn > tlo && tn < thi) || (fabs(2.0*y) > fabs(dxold*dy))) {
      tn = 0.5*(tlo + thi);
    }
    dx = fabs(tn - tp);
    conv = (dx <= M1_IMPL_TRTOL*fabs(tn)) || (fabs(thi - tlo) <= M1_IMPL_TRTOL*fabs(tp));
    tp = tn;
  }
  tout = tp;
  return nev;
}

//----------------------------------------------------------------------------------------
//! \fn M1ImplTemperature
//! \brief the pre-3g entry point: the same safeguarded root find, reading the EOS table
//! directly.  Kept so that every call site that does not use the cache makes exactly the
//! sequence of EOS_Data::ThermoAt() calls it made before, bit for bit.

template <class EosT>
KOKKOS_INLINE_FUNCTION
int M1ImplTemperature(const EosT &eos, const Real dd, const Real tguess,
                      const Real egn, const Real cdtkp_a, const Real rhs_abs,
                      Real &tout, bool &ok) {
  M1EosDirect<EosT> th{eos};
  return M1ImplTemperatureT(th, dd, tguess, egn, cdtkp_a, rhs_abs, tout, ok);
}

} // namespace radm1
#endif // RAD_M1_RAD_M1_IMPLICIT_HPP_
