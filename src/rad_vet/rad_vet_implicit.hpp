#ifndef RAD_VET_RAD_VET_IMPLICIT_HPP_
#define RAD_VET_RAD_VET_IMPLICIT_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_implicit.hpp
//! \brief constants and device-inline helpers of the IMPLICIT x1 column transport,
//! <rad_vet>/transport = implicit_x1 (docs/dev/rad_vet_implicit_design.md, milestone 3a).
//!
//! Nothing in here is reachable from the explicit path: the file is included only by
//! rad_vet_implicit.cpp and by the small delimited hooks of rad_vet.cpp /
//! rad_vet_tasks.cpp / rad_vet_newdt.cpp, so `transport = explicit` stays bitwise what it
//! was.

#include <math.h>

#include "athena.hpp"
#include "eos/eos_table.hpp"

namespace radvet {

// implicit_predictor = step as RESTART STATE: the marker of the little header block
// (int32 have, int32 ncomp, Real pred_dt) written behind the two-stream warm-start
// header, and read by the same eight-byte peek (src/outputs/restart.cpp, pgen.cpp).
// Present only when a stored increment exists; its ncomp per-cell slabs of ipred then
// follow the warm-start slabs at the end of every MeshBlock record.
constexpr char kM1PredRstMagic[8] = {'M', '1', 'P', 'R', 'E', 'D', '0', '1'};

// <rad_vet>/time_scheme (docs/dev/rad_vet_time2_design.md, tests_m1/runs_3x_hesdirk2)
constexpr int M1_TIME_BE = 0;         // Lie splitting: Heun hydro, then one BE solve
constexpr int M1_TIME_HESDIRK2 = 1;   // H-ESDIRK2: one stage solve inside each Heun stage
// what the next ImplicitSolve does (RadiationVet::t2_solve)
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
// <rad_vet>/implicit_one_pass (tests_m1/runs_4a_accel): a marked header behind the
// hesdirk2 one -- 9 Reals: the last two measured contractions and the solve counters of
// the three kinds of solve.  Written only when implicit_one_pass > 0.
constexpr char kM1OnePassRstMagic[8] = {'M', '1', 'O', 'N', 'E', 'P', '0', '1'};
// implicit_mr_every > 1 or vet_sc_every > 1 (rad_vet_mr.cpp): the multi-rate window and
// tensor-cadence state, behind the one-pass header: 10 Reals, no slabs
constexpr char kM1MRRstMagic[8] = {'M', '1', 'M', 'R', 'W', 'I', 'N', '1'};
// the hesdirk2 time2_auto_be back-off state (int32 len, int32 left), written
// only when a skip run is active (len > 0 or left > 0), behind the multi-rate header; a
// file without it restarts in stage mode, as before (and stays readable by older builds)
constexpr char kM1AutoBeRstMagic[8] = {'M', '1', 'A', 'U', 'T', 'O', 'B', 'E'};
// implicit_closure_thin_relax (runs_5c_thinstab): the relaxed closure of the last step
// (ctr_mem: chi, n1, n2, n3), header behind EVERY other one (after kCkRstMagic): int32
// have, int32 nch = 4; its 4 slabs are the LAST of each MeshBlock record.  Written only
// when the relaxation is active and has a stored step (ctr_init); a file without it
// restarts with an unrelaxed first step (the behaviour before 2026-09-27).
constexpr char kM1CtrRstMagic[8] = {'M', '1', 'C', 'T', 'R', 'L', 'X', '1'};
// <hydro>/rad_signal_speed (RadiationVet::RssPrime): the signal-speed inputs of the
// last M1 step, header behind the thin_relax one: int32 have = 1, nch, mode, 0; its nch
// slabs (opac(M1_OP_T), then the closure tensor channels of RssSlot) are the LAST of
// each MeshBlock record.  Written only when rad_signal_speed is on; a file without it
// restarts with opac filled from the restored state (no longer bitwise).
constexpr char kM1RssRstMagic[8] = {'M', '1', 'R', 'S', 'S', 'I', 'G', '1'};
// <rad_vet>/vet_gd_async (GD_ASYNC.md): the source (ln chi, ln S = vlat_cs) of the
// pending gd build, header behind the rad_signal_speed one: int32 have = 1, nch = 2,
// the build's ncycle, 0; its 2 slabs are the LAST of each MeshBlock record.  Written
// only when vet_gd_async is on with a pending build; a file without it restarts with a
// synchronous first build (D of the restart cycle unlagged: not bitwise).
constexpr char kM1VgdRstMagic[8] = {'M', '1', 'V', 'G', 'D', 'A', 'S', '1'};
// <rad_vet>/vet_rotate_avg: the averaging history (the ring of offsets,
// the pending Delta and the scalars), header behind the vet_gd_async one: int64 x
// RadiationVet::M1_RA_RST_NH (RotAvgRstHdr); its nslab slabs are the LAST of each
// MeshBlock record.  Written only with a complete history; a file without it restarts
// re-measuring the history at the first build (not bitwise).
// <rad_vet>/implicit_gas_src = conservative(2) with hesdirk2: the stored
// cell-F slope csk1 of the write-back's F^old increment and the first-solve flag, header
// behind the vet_gd_async one: int32 have = 1, nch (3 or 0), csk1_ok, 0; its nch slabs
// follow the vet_gd_async ones.  Written once a conservative solve has run; a file
// without it restarts with the cold-start F^old (gcs_first) and the face/E increments
// (not bitwise).
constexpr char kM1CskRstMagic[8] = {'M', '1', 'C', 'S', 'K', '1', 'S', '1'};
// time2_face_fct: the face fluxes t2pf1 (RadiationVet::FaxRst*; old files: also fax*)
constexpr char kM1FaxRstMagic[8] = {'M', '1', 'F', 'A', 'X', 'R', 'S', '1'};
constexpr char kM1RaRstMagic[8] = {'M', '1', 'R', 'O', 'T', 'A', 'V', '1'};
// time2_stage_local_be: the carried stage-2 drop of old files (not read)
constexpr char kM1ThcRstMagic[8] = {'M', '1', 'T', 'H', 'C', 'R', 'S', '1'};
// vet_gd_rotate_xfade: the old scheme's cross-fade state, behind the thc slabs.  Read
// past; no longer written
constexpr char kM1XfdRstMagic[8] = {'M', '1', 'X', 'F', 'D', 'R', 'S', '1'};
constexpr int kM1XfdRstNi = 24;   // int64 header words
constexpr int kM1XfdRstNr = 8;    // Real header words

                                             // (line Jacobi) or wrapped in a Krylov
                                             // iteration.

// <rad_vet>/implicit_solver: how the full 7-point system is solved (transport = implicit)
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

// <rad_vet>/implicit_offdiag: what is done with the OFF-DIAGONAL Eddington terms
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

// <rad_vet>/implicit_trans_limit: the REALIZABILITY limiter of the transverse (x2/x3)
// face fluxes.  The face-eliminated flux F_f' = theta_f [...] has no free-streaming
// bound: in steady state it is the unlimited diffusive flux F = -c grad P/(rho kappa),
// which is super-luminal wherever rho kappa is tiny (the optically thin top of the 2-D
// He slab, where |F_2|/(c E) saturates at the post-solve clip of 1).
constexpr int M1_TLIM_NONE = 0;  // no limiter; the arithmetic of phase D, bit for bit.
// realizability bound max(kt, |G|/(fmax E)), tapered by min(1, kt dx)
constexpr int M1_TLIM_LPR  = 2;
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
// transverse (x2) OUTFLOW face, read only when named: the face carries the edge cell's
// own outgoing x2 flux (zero gradient of F, outflow only; lagged in the Picard pass), F_f
// = max(+-F2_cell, 0); never the default (auto stays reflect)
constexpr int M1_IBC_OUTFLOW  = 5;
// transverse (x2) IMPOSED-flux face (read only when named): F_f = ifx2_lo/hi(m,k,i)
// from the problem generator (RadiationVet::SetX2FaceFlux)
constexpr int M1_IBC_X2FLUX   = 6;
                                     // at the start of every step, i.e. it is frozen at
                                     // its initial value.  This is the face-flux form of
                                     // the fixed-E ghost the explicit tests T3b and T6
                                     // use, and the only boundary that anchors the LEVEL
                                     // of E in a pure-scattering column (imposed flux on
                                     // both ends leaves the operator singular).

// <rad_vet>/implicit_accel: how the OUTER (Picard) fixed-point iteration on the lagged
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

// <rad_vet>/implicit_flux: which spatial form the implicit E-flux takes at a face
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
                                     // (see rad_vet_implicit.cpp).
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

// <rad_vet>/implicit_blend: how the per-face weight w_f of M1_IFLUX_BLEND is built.
constexpr int M1_IBLEND_TAU  = 0;    // w = exp(-(tau_face/tau0)^2): upwind where the face
                                     // is optically THIN, central where it is thick
                                     // (Jiang 2021 uses 1 - exp(-tau^2) for the opposite
                                     // sense).  <rad_vet>/implicit_blend_tau0.
constexpr int M1_IBLEND_F    = 1;    // w = smoothstep((f_face - f_lo)/(f_hi - f_lo)) on
                                     // the LAGGED comoving reduced flux at the face: a
                                     // DIFFUSE field (f -> 1/2, which is what a grey
                                     // surface with marshak_q = 1/2 has) stays central, a
                                     // beam or a front (f -> 1) goes upwind.
constexpr int M1_IBLEND_IDF = 5;     // idort x smoothstep(|H|/J)
constexpr int M1_IBLEND_IDA = 6;     // idort x smoothstep(h asymmetry)
constexpr int M1_IBLEND_KN = 4;      // (half-range only): Knudsen x tau
constexpr int M1_IBLEND_IDORT = 3;   // (half-range only): w = 1/(1 + x + x^2/
                                     // tau0), x = implicit_blend_alpha tau_face
constexpr int M1_IBLEND_TAUF = 2;    // the product of the two: upwind only where the face
                                     // is thin AND the field is beamed.

// <rad_vet>/implicit_thin_switch, the dt-free per-face weight of the
// half-range blend (M1ThinFace, rad_vet_implicit.cpp); none keeps implicit_blend/xthin.
constexpr int M1_THSW_NONE = 0;
constexpr int M1_THSW_HJ   = 1;   // smoothstep of the formal-solution |H|/J
constexpr int M1_THSW_KN   = 2;   // the face Knudsen number of E^n
constexpr int M1_THSW_FS   = 3;   // smoothstep of max(|H|/J, C1/1.5, ...), C1 = <|I-S|>/J
constexpr int M1_THSW_HJKN = 4;   // max(hj, kn)
constexpr int M1_THSW_FSKN = 5;   // max(fs, kn)
constexpr int M1_THSW_FSI  = 7;   // fs weight; h+- -> isotropic +-1/4 where the rays
                                  // do not describe E^n (M1ThinHr)
constexpr int M1_THSW_FSE  = 8;   // max(|H|/J, eps C1/1.5), eps the thermal fraction
constexpr int M1_THSW_FSQ  = 6;   // max(|H|/J, q_c C1/1.5), q_c = min(E^n/J, J/E^n)^qp
// <rad_vet>/implicit_thin_corr: the lagged formal-solution face flux (C4)
constexpr int M1_THCR_NONE = 0;
constexpr int M1_THCR_LAG  = 1;   // + w (F_sc^n - F_hr(E^n)), a fixed source
constexpr int M1_THCR_SC   = 2;   // w F_sc^n replaces the half-range part
constexpr int M1_THCR_HCS  = 4;   // hc weighted by the scattering fraction 1 - eps
constexpr int M1_THCR_HC   = 3;   // + w (F_sc,central - F_sc,upwind)^n (no E^n term)

// <rad_vet>/implicit_blend_fmode: which of the two cells' lagged reduced fluxes sets
// f_face.
constexpr int M1_IBFM_MAX  = 0;
constexpr int M1_IBFM_MEAN = 1;

// <rad_vet>/implicit_blend_mode: WHAT is blended.
constexpr int M1_IBMODE_FLUX   = 0;  // the whole face flux (the convex blend above)
constexpr int M1_IBMODE_DISSIP = 1;  // the central flux in FULL plus w_f times the HLL
                                     // DISSIPATION alone: F = F_central + w alpha dk
                                     // (E_L - E_R).  The extra term contributes +w dk to
                                     // the diagonal of the left cell and -w dk to the
                                     // upper off-diagonal (and the mirror image to the
                                     // right cell), so the column sums are untouched and
                                     // the M-matrix survives as well.

// <rad_vet>/implicit_recon: the states the HLL part of the implicit flux is built from.
constexpr int M1_IRECON_DC    = 0;   // piecewise constant; the matrix IS the operator
constexpr int M1_IRECON_PLMDC = 1;   // plm as a DEFERRED CORRECTION: the difference
                                     // between the plm and the dc flux, evaluated at the
                                     // PREVIOUS Picard iterate, goes into the right-hand
                                     // side, so the matrix stays the low-order M-matrix.

// <rad_vet>/implicit_enthalpy: the face value of the advective enthalpy flux A_d = a_d E
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
// HIDDEN TERM diagnostics (read only when named; every mode >= PLM takes the
// plm a_f and stencil rule, only E_f differs):
constexpr int M1_IENTH_PLMC    = 3;  // E_f = mean of the two van Leer plm face values (no
                                     // upwind-side switch, limiter kept)
// E_f the UNLIMITED linear (Fromm) value from the side
constexpr int M1_IENTH_LIN     = 4;
                                     // upwind of a_f (side switch kept, no limiter)

// <rad_vet>/implicit_vimp: the gas velocity v' of the enthalpy flux a(v') E' taken
// IMPLICIT through the radiative force of the same solve (tests_m1/runs_3v_vimplicit).
// Per face f of axis d the flux gains
//   E_f^k [ (1 + D_dd) dv_d^k + sum_{e != d} D_de dv_e^k ]_f
//     + E_f^k (1 + D_dd)_f [ (P dv_d)(E') - (P dv_d)(E^k) ]_f,
// dv^k the velocity change of the ITERATE's face fluxes (the write-back rule) and P the
// Jacobian of dv_d in E' (the face-normal eliminated flux, od/g0 terms lagged).  At the
// Picard fixed point the flux is E_f a_f(v') whatever the Jacobian; the operator couples
// cells two apart along every axis.  The components below are APPENDED to iw (offset
// RadiationVet::iw_vimp); the first M1_NVIMP_X are exchanged once per pass.
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
// <rad_vet>/implicit_hr_recon = plm (Fix B): the half-range face values
// E_L + s_L/2, E_R - s_R/2 with the limited slope s_c = sig_c (E_c+1 - E_c-1)/2, sig_c
// in [0, 1] the van Leer limiter of the Picard-lagged E, frozen for the pass, applied
// as a linear operator outside the 7-point row (M1MusclRow): sig per direction, then
// the row coefficients diag, x1 -2 -1 +1 +2, x2 -2 -1 +1 +2, x3 -2 -1 +1 +2.
constexpr int M1_IM_SIG   = 0;    // sig_1, sig_2, sig_3
constexpr int M1_IM_D     = 3;
constexpr int M1_IM_X1    = 4;    // -2, -1, +1, +2
constexpr int M1_IM_X2    = 8;
constexpr int M1_IM_X3    = 12;
constexpr int M1_IM_KILL  = 16;   // 1: dc for the rest of the step (positivity)
constexpr int M1_IM_B     = 17;   // led: b_1, b_2, b_3 (SIG holds a_d)
constexpr int M1_IM_EPREV = 20;   // qs_mode = step: E^n of the previous step
constexpr int M1_NIW_MUSCL = 21;

// <rad_vet>/implicit_partition: how a column that spans several MeshBlocks (and ranks) is
// solved (milestone 3a2, LIMIT 4).
constexpr int M1_IPART_NONE   = 0;   // one MeshBlock per column; fatal otherwise (3a)
constexpr int M1_IPART_GATHER = 1;   // the column's tridiagonal rows are GATHERED onto
                                     // the rank that owns its lowest block, solved by the
                                     // very same serial Thomas sweep, and scattered back.
                                     // Bitwise independent of the partition by
                                     // construction: the arithmetic order is identical.

// components of the implicit FACE work array RadiationVet::ifw (m,n,k,j,i), i the x1 FACE
// index (face i is the one between cells i-1 and i), filled once per Picard iteration by
// the implicit_flux = ap_hll path and identically zero under `central`.
// vet-positivity counters (RadiationVet::pos_cnt_d / pos_cnt)
constexpr int M1_POS_G0     = 0;   // cell-passes whose lagged g0 was clipped
constexpr int M1_POS_GAS    = 1;   // cell-solves where the gas eint limiter fired
constexpr int M1_POS_GAS_DE = 2;   // energy (erg, volume-weighted) moved gas <- radiation
constexpr int M1_POS_FLR    = 3;   // cell-solves where E was raised to e_floor
constexpr int M1_POS_FLR_DE = 4;   // floor energy taken from the gas (erg)
constexpr int M1_POS_FLR_UN = 5;   // floor energy NOT covered by the gas (created, erg)
constexpr int M1_POS_FCLIP  = 6;   // cell-solves whose |F| > c E was scaled back
constexpr int M1_POS_FCLIPM = 7;   // sum of (|F|/(c E) - 1) over those cell-solves
constexpr int M1_POS_FLRS   = 8;   // F3: residual raise below 1e-6 e_floor (x volume)
constexpr int M1_POS_FLR_MV = 9;   // implicit_pos_floor_cons: floor energy taken
                                   // from the column's radiation (erg, gas units x
                                   // volume)
constexpr int M1_POS_N      = 10;
constexpr int M1_IFW_AL  = 0;   // alpha, the asymptotic-preserving weight of F_HLL
constexpr int M1_IFW_HCL = 1;   // alpha * (coefficient of E'_L in F_HLL), >= 0
constexpr int M1_IFW_HCR = 2;   // alpha * (coefficient of E'_R in F_HLL), <= 0
constexpr int M1_IFW_DG  = 3;   // the plm deferred correction, alpha (F^plm - F^dc)
constexpr int M1_NIFW = 4;

// components of the implicit work array RadiationVet::iw (m,n,k,j,i).  All are per-cell
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

// ---- MILESTONE 3g: the GAS-RADIATION energy coupling, <rad_vet>/implicit_gas_newton and
// <rad_vet>/implicit_eos_cache.  Both default to false, neither allocates anything then
// and neither is referenced, so every earlier configuration is bitwise unchanged.
//
// The five components below are APPENDED after whichever base component count the solve
// uses (M1_NIW_X1 / M1_NIW / M1_NIW_K), so the base indices above never move; the offset
// of the block is RadiationVet::iw_gas, a runtime int the kernels capture.
// (implicit_hr_face = sf): per cell, the formal-solution half-range fluxes
// H+_d, H-_d (d = 1..3), J and the donor scale Jh (M1HrSf), appended like the gas block
constexpr int M1_NIW_HSF = 8;
constexpr int M1_HSF_HP = 0;   // + 2 d: H+_d = sum w I max(mu_d, 0) (E units)
constexpr int M1_HSF_HM = 1;   // + 2 d: H-_d = H_d - H+_d <= 0
constexpr int M1_HSF_J  = 6;   // J of the formal solution
// Jh = max(J, max over the faces of the outgoing face SC flux)
constexpr int M1_HSF_JH = 7;
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
// implicit_opac_newton: the relative temperature step of the opacity derivative
constexpr Real M1_OPN_H = 1.0e-5;
// implicit_opac_newton: one cell's share of the Newton term of a face.  q = -chat dt th
// G(E^k)/2 of the face, kt = d(rho kappa_T)/dT, bk/rk = the gas Newton step's B_k/R_k
// and ke = rho kappa_E of that cell: its E' coefficient goes into coef (the row entry of
// that cell), the rest to the right-hand side.  Cartesian and sp rows share it.
KOKKOS_INLINE_FUNCTION
void M1OpnCell(const Real q, const Real kt, const Real bk, const Real rk, const Real ke,
               const Real cl, const Real dt, Real &coef, Real &rr) {
  if (bk > 0.0) {
    const Real s = q*kt/bk;
    coef += s*cl*dt*ke;
    rr -= s*rk;
  }
}
// BiCGStab breakdown thresholds: |rho| and |rhat.v| below these times the scale of the
// right-hand side mean the shadow residual has become orthogonal to the Krylov space.
constexpr Real M1_BCG_EPS = 1.0e-300;
// implicit_krylov_dev (rad_vet_launch.cpp): the size of its device scalar array
constexpr int M1_KD_SIZE = 20;

//----------------------------------------------------------------------------------------
//! \fn M1EddDiag
//! \brief the DIAGONAL Eddington-tensor component D_dd = P_dd/E of direction d (0,1,2)
//! from the lagged closure: D_ab = (1-chi)/2 delta_ab + (3 chi - 1)/2 n_a n_b.

KOKKOS_INLINE_FUNCTION
Real M1EddDiag(const Real chi, const Real nd) {
  return 0.5*(1.0 - chi) + 0.5*(3.0*chi - 1.0)*nd*nd;
}

//----------------------------------------------------------------------------------------
//! \fn M1POff
//! \brief the OFF-diagonal radiation pressure P_ab = D_ab E at one cell, from the LAGGED
//! closure of the previous Picard pass.  a, b are 0-based directions.  `ec` names the
//! component of the work array the energy is taken from: M1_IW_EP for the lagged form,
//! and a KRYLOV VECTOR when the term is applied as part of the linear operator
//! (implicit_offdiag = operator), where the closure is frozen and the term is linear.

//----------------------------------------------------------------------------------------
//! \fn M1CsSeamGeom
//! \brief STAGE CS1 (cubed sphere): the CANONICAL geometry of the panel-seam face next to
//! along-seam cell q (0..n-1, n cells per panel edge, equiangular) of a seam, as angles
//! on the unit sphere: dth = the seam arc of the face (area 0.5 (r_r^2 - r_l^2) dth),
//! angm = the angle between the cell centre next to the seam and its MIRROR image across
//! the seam plane (the neighbour panel's cell next to the seam, which is what the
//! no-resample halo puts in the first ghost; the pair chord is normal to the seam),
//! sfoot = the along-seam position of the foot of that normal, smid = the along-seam
//! position of the face midpoint (both from the seam midpoint, increasing with q).
//!
//! Everything is a function of the integers (q, n) alone, evaluated for min(q, n-1-q)
//! with the along-seam positions negated for the reflected half: the two panels sharing
//! the face get BITWISE the same dth and angm (and opposite-oriented sfoot, smid), so the
//! two-point seam coefficient is single valued and the seam flux telescopes.

KOKKOS_INLINE_FUNCTION
void M1CsSeamGeom(const int q, const int n, Real &dth, Real &angm, Real &sfoot,
                  Real &smid) {
  const bool flip = (q > n - 1 - q);
  const int qq = flip ? (n - 1 - q) : q;
  const Real dl = 0.5*M_PI/static_cast<Real>(n);
  const Real el = -0.25*M_PI + static_cast<Real>(qq)*dl;
  const Real er = -0.25*M_PI + static_cast<Real>(qq + 1)*dl;
  const Real ec = -0.25*M_PI + (static_cast<Real>(qq) + 0.5)*dl;
  const Real xc = tan(0.25*M_PI - 0.5*dl);
  const Real yl = tan(el), yr = tan(er), yc = tan(ec);
  const Real r2 = sqrt(2.0);
  const Real sl = atan(yl/r2), sr = atan(yr/r2);
  dth = sr - sl;
  smid = 0.5*(sl + sr);
  sfoot = atan(yc*r2/(1.0 + xc));
  angm = 2.0*asin((1.0 - xc)/(r2*sqrt(1.0 + xc*xc + yc*yc)));
  if (flip) {
    smid = -smid;
    sfoot = -sfoot;
  }
}

//----------------------------------------------------------------------------------------
//! \fn M1CsSeamMid
//! \brief the value at the face MIDPOINT (smid, slot 3 of the seam geometry
//! g = csg2 / csg3) of a seam quantity known at the FEET of the mirror-pair chords
//! (sfoot, slot 2) of the seam faces c-1, c, c+1 along the seam: the quadratic through
//! the three foot values (one-sided c, c+1, c+2 / c-2, c-1, c at a cube vertex, where
//! up / dn is false).  cs-s4-1017 moved the foot value by (smid - sfoot) d/ds, which
//! drops the (smid - sfoot)^2/2 d2/ds2 term: an O(h^2) face error that the O(h) seam
//! cell turns into a first-order seam truncation (the cs MMS); the quadratic is exact
//! to O(h^3) (second-order seam and vertex cells).  fv(c') is the foot value of face c'.

//----------------------------------------------------------------------------------------
//! \fn M1CsSeamDepth
//! \brief depth of index c (active lo..hi, ghosts allowed) below the nearest block face
//! along that axis that is a PANEL SEAM (slo / shi): 0 for the cell next to the seam and
//! for the first ghost beyond it (the neighbour panel's seam cell), 1 one cell further,
//! ...; a large number when neither face is a seam.

KOKKOS_INLINE_FUNCTION
int M1CsSeamDepth(const int c, const int lo, const int hi, const bool slo,
                  const bool shi) {
  int d = 1 << 20;
  if (slo) {d = (c < lo) ? (lo - 1 - c) : (c - lo);}
  if (shi) {
    const int e = (c > hi) ? (c - hi - 1) : (hi - c);
    d = (e < d) ? e : d;
  }
  return d;
}

//! \brief Chebyshev distance of cell (k,j) to the cube vertex of its block
//! (both tangential depths finite only in a block holding a panel corner)

KOKKOS_INLINE_FUNCTION
int M1CsVtxDist(const int k, const int j, const int ks, const int ke, const int js,
                const int je, const bool s2lo, const bool s2hi, const bool s3lo,
                const bool s3hi) {
  const int dj = M1CsSeamDepth(j, js, je, s2lo, s2hi);
  const int dk = M1CsSeamDepth(k, ks, ke, s3lo, s3hi);
  return (dj > dk) ? dj : dk;
}

template <class F, class G>
KOKKOS_INLINE_FUNCTION
Real M1CsSeamMid(const F &fv, const G &g, const int m, const int c, const int sd,
                 const bool up, const bool dn) {
  int a, b, e;
  if (up && dn) {
    a = c - 1; b = c; e = c + 1;
  } else if (up) {
    a = c; b = c + 1; e = c + 2;
  } else if (dn) {
    a = c - 2; b = c - 1; e = c;
  } else {
    return fv(c);
  }
  const Real sa = g(m,2,a,sd), sb = g(m,2,b,sd), se = g(m,2,e,sd), s = g(m,3,c,sd);
  return fv(a)*((s - sb)*(s - se))/((sa - sb)*(sa - se))
         + fv(b)*((s - sa)*(s - se))/((sb - sa)*(sb - se))
         + fv(e)*((s - sa)*(s - sb))/((se - sa)*(se - sb));
}


//----------------------------------------------------------------------------------------
//! \fn M1CsTanDer
//! \brief STAGE CS1: the derivative along the OTHER tangential axis (e = 1: x2, e = 2:
//! x3) of D_ee E at cell (k,j,i), per unit arc length, for the lagged cross term of the
//! skewed face gradient.  Centred over the two centre-to-centre arcs dxface; ONE-SIDED
//! (second order, three cells inside the panel) at a cell next to a panel seam along e,
//! whose ghost there is a mirror cell off this cell's coordinate line.  slo/shi: the
//! block's low/high face along e is a panel seam.

// seamE-1106 cs_cross_mono: the one-sided derivative d at the cell next to a seam, kept
// between 0 and 2 s1 (s1 = the one-sided first difference quotient to the first
// neighbour), but only where |q0 - 2 q1 + q2| > k (|q0| + |q1| + |q2|)/3 with
// k = cs_cross_mono_k.  The three-point formula extrapolates across a sharp front
// (E = 20, 174, 3390 gives the wrong sign and 5 times the size), and through the skew
// cross term that reaches the face flux and the RHS.  k = 0, the DEFAULT when
// cs_cross_mono_k is not named, is the first form: the clamp acts wherever the second
// difference is nonzero, also at smooth extrema, where Linf falls to first order (audit
// 11-07).  Only k > 0 leaves smooth data alone: a sinusoid with (wavenumber h)^2 < k
// does not reach the clamp, so the formula stays second order there.
KOKKOS_INLINE_FUNCTION
Real M1CsTanMono(const Real d, const Real s1, const Real q0, const Real q1, const Real q2,
                 const Real fk) {
  const Real a0 = (q0 < 0.0) ? -q0 : q0, a1 = (q1 < 0.0) ? -q1 : q1,
             a2 = (q2 < 0.0) ? -q2 : q2;
  const Real c2 = q0 - 2.0*q1 + q2;
  if ((c2 < 0.0 ? -c2 : c2) <= (fk/3.0)*(a0 + a1 + a2)) {return d;}
  const Real lo = (s1 > 0.0) ? 0.0 : 2.0*s1, hi = (s1 > 0.0) ? 2.0*s1 : 0.0;
  return (d < lo) ? lo : ((d > hi) ? hi : d);
}

// the derivative itself, of any cell field q(o) at offset o along the axis,
// with h(o) the arc between offsets o-1 and o (shared by M1CsTanDer and M1CsTanFaceLag)
template <class Q, class H>
KOKKOS_INLINE_FUNCTION
Real M1CsAxisDer(const Q &q, const H &h, const int c, const int lo, const int hi,
                 const bool slo, const bool shi, const Real mono = -1.0) {
  if (c == lo && slo) {
    const Real h1 = h(1), h2 = h(2);
    const Real d = (-(2.0*h1 + h2)/(h1*(h1 + h2)))*q(0) + ((h1 + h2)/(h1*h2))*q(1)
                   - (h1/(h2*(h1 + h2)))*q(2);
    return (mono >= 0.0) ? M1CsTanMono(d, (q(1) - q(0))/h1, q(0), q(1), q(2), mono) : d;
  }
  if (c == hi && shi) {
    const Real h1 = h(0), h2 = h(-1);
    const Real d = ((2.0*h1 + h2)/(h1*(h1 + h2)))*q(0) - ((h1 + h2)/(h1*h2))*q(-1)
                   + (h1/(h2*(h1 + h2)))*q(-2);
    return (mono >= 0.0) ? M1CsTanMono(d, (q(0) - q(-1))/h1, q(0), q(-1), q(-2), mono)
                         : d;
  }
  return (q(1) - q(-1))/(h(0) + h(1));
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

//! <mesh>/cs_vertex_3cell (vx.on): at a CUBE VERTEX the ng x ng corner
//! ghost block has no physical counterpart.  A seam GHOST cell (column j outside
//! [js,je], or row k outside [ks,ke]) next to a vertex corner takes the derivative ALONG
//! its seam one-sided, from the cell itself and its in-range neighbour, instead of the
//! centred difference across the corner block.  Off (the default argument): bitwise.
struct M1Vtx {
  bool on = false;
  int js = 0, je = 0, ks = 0, ke = 0;
  bool jl = false, jh = false, kl = false, kh = false;   // panel edges of the block
};
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

//! \fn M1SphMarshakCoef
//! \brief <rad_vet>/implicit_marshak_face = linear (m1-sp-order2, spherical-polar wedge,
//! tests_m1/runs_5o_sporder2): the E of a Marshak boundary FACE at r_f, from the end cell
//! (E ec at its centroid rc) and its interior neighbour (en at rn), by the linear
//! extrapolation of G = r^2 E to the face: second order, and exact for the free-
//! streaming state r^2 E = const.  E_lin = ca ec - cb en with
//!   alpha = (r_f - rc)/(rc - rn) > 0 (either end),  ca = (1 + alpha) (rc/r_f)^2,
//!   cb = alpha (rn/r_f)^2,  ca > cb >= 0.
//! The implicit row carries E_lin itself (ca on the diagonal, -cb on the neighbour: a
//! non-positive off-diagonal and a diagonally dominant row, so the M-matrix sign pattern
//! stays).

KOKKOS_INLINE_FUNCTION
void M1SphMarshakCoef(const Real rc, const Real rn, const Real rf, Real &ca, Real &cb) {
  const Real al = (rf - rc)/(rc - rn);
  ca = (1.0 + al)*SQR(rc/rf);
  cb = al*SQR(rn/rf);
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

//! the overload with x2t takes sin / cot of the cell from x2t (the theta
//! of the data the row holds, RadiationVet::VtanX2t; see M1SphTan): at a theta-PERIODIC
//! seam the d = 1 curvature cot (P_tt - P_pp)/r of a ghost row is then the value of its
//! source cell, and the seam face of the two blocks carries one flux.  Coordinate
//! distances stay x2v.  The overload without x2t is x2t = x2v (bitwise the old form).
template <class V, class A>
KOKKOS_INLINE_FUNCTION
Real M1SphCurv(const V &iw, const A &x1v, const A &x2v, const A &x2t, const A &x3v,
               const int m, const int d, const int k, const int j, const int i,
               const bool od, const bool thrd, const int il, const int iu, const int jl,
               const int ju, const int kl, const int ku, const int ec) {
  const Real r = x1v(m,i);
  const Real sn = sin(x2t(m,j));
  const Real ct = cos(x2t(m,j))/sn;
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

template <class V, class A>
KOKKOS_INLINE_FUNCTION
Real M1SphCurv(const V &iw, const A &x1v, const A &x2v, const A &x3v, const int m,
               const int d, const int k, const int j, const int i, const bool od,
               const bool thrd, const int il, const int iu, const int jl, const int ju,
               const int kl, const int ku, const int ec) {
  return M1SphCurv(iw, x1v, x2v, x2v, x3v, m, d, k, j, i, od, thrd, il, iu, jl, ju, kl,
                   ku, ec);
}

//----------------------------------------------------------------------------------------
//! \fn M1SphLat
//! \brief <rad_vet>/vet_col_lat (m1-vetcol-lat): the LAGGED divergence of the lateral
//! off-diagonal pressure P_ra = D_ra E (a = theta, phi) at one cell centre of the sp
//! wedge.  D_rt, D_rp are slots c0+1, c0+2 of tau_ten (c0 = M1_TT_LAT0; the short-
//! characteristics correction of rad_vet_vetlat.cpp, ghosts exchanged), E is component
//! `ec` of iw:
//!   d = 0:  (1/r) d_th P_rt + cot P_rt / r + (1/(r sin)) d_ph P_rp
//!   d = 1:  (1/r^3) d_r (r^3 P_rt)
//!   d = 2:  (1/r^3) d_r (r^3 P_rp)
//! (the [od] terms of M1SphCurv with P_tp = 0), centred over the coordinate distances,
//! one-sided at a physical boundary exactly as M1SphCurv.

template <class V, class A>
KOKKOS_INLINE_FUNCTION
Real M1SphLat(const V &iw, const V &tt, const int c0, const A &x1v, const A &x2v,
              const A &x3v, const int m, const int d, const int k, const int j,
              const int i, const bool thrd, const int il, const int iu, const int jl,
              const int ju, const int kl, const int ku, const int ec) {
  const Real r = x1v(m,i);
  const int ia = (i+1 <= iu) ? (i+1) : i;
  const int ib = (i-1 >= il) ? (i-1) : i;
  const int ja = (j+1 <= ju) ? (j+1) : j;
  const int jb = (j-1 >= jl) ? (j-1) : j;
  const int ka = (thrd && k+1 <= ku) ? (k+1) : k;
  const int kb = (thrd && k-1 >= kl) ? (k-1) : k;
  Real s = 0.0;
  if (d == 0) {
    const Real sn = sin(x2v(m,j));
    const Real ct = cos(x2v(m,j))/sn;
    if (ja != jb) {
      s += (tt(m,c0+1,k,ja,i)*iw(m,ec,k,ja,i) - tt(m,c0+1,k,jb,i)*iw(m,ec,k,jb,i))
           /(r*(x2v(m,ja) - x2v(m,jb)));
    }
    s += ct*tt(m,c0+1,k,j,i)*iw(m,ec,k,j,i)/r;
    if (ka != kb) {
      s += (tt(m,c0+2,ka,j,i)*iw(m,ec,ka,j,i) - tt(m,c0+2,kb,j,i)*iw(m,ec,kb,j,i))
           /(r*sn*(x3v(m,ka) - x3v(m,kb)));
    }
  } else if (ia != ib) {
    const int c = c0 + d;
    const Real ra = x1v(m,ia), rb = x1v(m,ib);
    s = (ra*ra*ra*tt(m,c,k,j,ia)*iw(m,ec,k,j,ia)
         - rb*rb*rb*tt(m,c,k,j,ib)*iw(m,ec,k,j,ib))
        /(r*r*r*(ra - rb));
  }
  return s;
}

//----------------------------------------------------------------------------------------
//! \fn M1SphTan
//! \brief <rad_vet>/vet_gd (m1-vet-gd): the LAGGED divergence of the TANGENTIAL
//! anisotropy of the Eddington tensor at one cell centre of the sp wedge.  The operator
//! carries the isotropic tangential pressure (1 - D_rr)/2 E; the deviation
//!   dP_tt = a E, dP_pp = -a E, P_tp = D_tp E,   a = tt(c0+3) = D_tt - (1 - D_rr)/2,
//! D_tp = tt(c0+4) (tau_ten slots of rad_vet_vetgd.cpp, ghosts exchanged), E = iw(ec),
//! adds (the trace-free part has no radial component)
//!   d = 1:  (1/(r sin)) d_th (sin dP_tt) + (1/(r sin)) d_ph P_tp - cot dP_pp / r
//!   d = 2:  (1/(r sin^2)) d_th (sin^2 P_tp) + (1/(r sin)) d_ph dP_pp
//! The NORMAL-gradient parts (1/r) d_th dP_tt (d = 1) and (1/(r sin)) d_ph dP_pp (d = 2)
//! are NOT here: the caller adds +a / -a to the face diagonal D_dd of the compact face
//! gradient (ImplicitTransverseTerms), exactly as the isotropic part.  A lagged wide
//! centred difference of them was unstable in the thin top (He 128^2: L_top grew e-fold
//! per ~90 s).  What remains here are zeroth-derivative and cross terms:
//!   d = 1:  2 cot a E / r + (1/(r sin)) d_ph P_tp
//!   d = 2:  (1/(r sin^2)) d_th (sin^2 P_tp)
//! centred over the coordinate distances, one-sided at a physical boundary as M1SphLat.
//! the trigonometric factors (sin, cot, sin^2) take the theta of the cell
//! whose DATA the slot holds, x2t (RadiationVet::VtanX2t: = x2v, except in the ghost rows
//! of a theta-PERIODIC wedge's x2min/x2max blocks, where it is the theta of the periodic
//! source cell).  With the plain x2v the two blocks at the theta seam evaluated the same
//! seam face with cot(x2min) ~ +1 and cot(x2max) ~ -1 for the same cells: the face flux
//! F2 of the seam differed between its two sides and the solve created energy (He 128^2
//! wedge +0.35 L_in, EAUDIT RSOLVE).  The coordinate distances stay x2v.

template <class V, class A>
KOKKOS_INLINE_FUNCTION
Real M1SphTan(const V &iw, const V &tt, const int c0, const A &x1v, const A &x2v,
              const A &x2t, const A &x3v, const int m, const int d, const int k,
              const int j, const int i, const bool thrd, const int jl, const int ju,
              const int kl, const int ku, const int ec) {
  const Real r = x1v(m,i);
  const int ja = (j+1 <= ju) ? (j+1) : j;
  const int jb = (j-1 >= jl) ? (j-1) : j;
  const int ka = (thrd && k+1 <= ku) ? (k+1) : k;
  const int kb = (thrd && k-1 >= kl) ? (k-1) : k;
  const Real sn = sin(x2t(m,j));
  const Real ct = cos(x2t(m,j))/sn;
  Real s = 0.0;
  if (d == 1) {
    if (ka != kb) {
      s += (tt(m,c0+4,ka,j,i)*iw(m,ec,ka,j,i) - tt(m,c0+4,kb,j,i)*iw(m,ec,kb,j,i))
           /(r*sn*(x3v(m,ka) - x3v(m,kb)));
    }
    s += 2.0*ct*tt(m,c0+3,k,j,i)*iw(m,ec,k,j,i)/r;
  } else if (d == 2) {
    if (ja != jb) {
      const Real sa = sin(x2t(m,ja)), sb = sin(x2t(m,jb));
      s += (sa*sa*tt(m,c0+4,k,ja,i)*iw(m,ec,k,ja,i)
            - sb*sb*tt(m,c0+4,k,jb,i)*iw(m,ec,k,jb,i))
           /(r*sn*sn*(x2v(m,ja) - x2v(m,jb)));
    }
  }
  return s;
}

//----------------------------------------------------------------------------------------
//! \fn M1CsLat0
//! \brief cs-s4-1017 (CS_VETGD.md sect. 4C): the cubed-sphere counterpart of M1SphLat
//! (d = 0) -- the lateral divergence of P_r,lat = v E at one cell, with v = D.rhat -
//! D_rr rhat the tangential part of the r row of the vet_gd tensor held as a CARTESIAN
//! vector in tau_ten slots c0+1..c0+3 (VetGdCsMoments; ghosts exchanged as scalars,
//! no resample: a seam ghost is the mirror cell).  Finite-volume face sum over the four
//! lateral faces, (1/V) sum_f +- A_f nhat_f . (w_c + w_nb)/2, w = v E (E = component ec
//! of iw), with the planar face normals nhat_f (fn2/fn3: Cartesian, constant along the
//! radial line) and the cs face areas (the canonical areas at a seam).  The x1 faces
//! carry nothing (w is tangential).  Every lateral face of a cs block has a neighbour
//! (block or seam), so there is no one-sided case.
//! corr (the lagged assembly; false in the operator/stencil, whose 19 points cannot hold
//! it): at a PANEL-SEAM face the mirror pair's mean is the value at the FOOT of the pair
//! chord, not at the face midpoint; it is moved along the seam by (smid - sfoot) d/ds
//! over the neighbouring seam faces (one-sided at a cube vertex), exactly as the seam
//! face gradient of ImplicitTransverseTerms (csg2/csg3 of M1CsSeamGeom).  Without it the
//! seam cells carry an O(1) error (the O(h) foot offset over the O(h) cell: the cs MMS).

template <class V, class VV, class VF, class VN, class VB, class VG>
KOKKOS_INLINE_FUNCTION
Real M1CsLat0(const V &iw, const V &tt, const int c0, const VV &cvol, const VF &car,
              const VN &fn2, const VN &fn3, const int m, const int k, const int j,
              const int i, const bool thrd, const int ec, const bool corr,
              const VB &cse, const VG &g2, const VG &g3, const int js, const int je,
              const int ks, const int ke) {
  auto wn = [&](const int kq, const int jq, const Real *nf) {
    return (tt(m,c0+1,kq,jq,i)*nf[0] + tt(m,c0+2,kq,jq,i)*nf[1]
            + tt(m,c0+3,kq,jq,i)*nf[2])*iw(m,ec,kq,jq,i);
  };
  const bool s2lo = corr && (cse(m,0) != 0) && (j == js);
  const bool s2hi = corr && (cse(m,1) != 0) && (j == je);
  const bool s3lo = corr && thrd && (cse(m,2) != 0) && (k == ks);
  const bool s3hi = corr && thrd && (cse(m,3) != 0) && (k == ke);
  const bool b3lo = (cse(m,2) != 0), b3hi = (cse(m,3) != 0);
  const bool b2lo = (cse(m,0) != 0), b2hi = (cse(m,1) != 0);
  // the seam face value of the x2 face jf (pair jf-1, jf) moved to the face midpoint
  auto seam2 = [&](const int jf, const int sd, const Real *nf) {
    auto fv = [&](const int kq) {return 0.5*(wn(kq,jf-1,nf) + wn(kq,jf,nf));};
    const bool up = !(k == ke && b3hi), dn = !(k == ks && b3lo);
    // the quadratic through the foot values, at the face midpoint
    return M1CsSeamMid(fv, g2, m, k, sd, up, dn);
  };
  auto seam3 = [&](const int kf, const int sd, const Real *nf) {
    auto fv = [&](const int jq) {return 0.5*(wn(kf-1,jq,nf) + wn(kf,jq,nf));};
    const bool up = !(j == je && b2hi), dn = !(j == js && b2lo);
    // the quadratic through the foot values, at the face midpoint
    return M1CsSeamMid(fv, g3, m, j, sd, up, dn);
  };
  Real s = 0.0;
  {
    const Real nh[3] = {fn2(m,j+1,0), fn2(m,j+1,1), fn2(m,j+1,2)};
    const Real nl[3] = {fn2(m,j,0), fn2(m,j,1), fn2(m,j,2)};
    s += car.x2f(m,k,j+1,i)*(s2hi ? seam2(j+1, 1, nh)
                                   : 0.5*(wn(k,j,nh) + wn(k,j+1,nh)));
    s -= car.x2f(m,k,j,i)*(s2lo ? seam2(j, 0, nl) : 0.5*(wn(k,j-1,nl) + wn(k,j,nl)));
  }
  if (thrd) {
    const Real nh[3] = {fn3(m,k+1,0), fn3(m,k+1,1), fn3(m,k+1,2)};
    const Real nl[3] = {fn3(m,k,0), fn3(m,k,1), fn3(m,k,2)};
    s += car.x3f(m,k+1,j,i)*(s3hi ? seam3(k+1, 1, nh)
                                   : 0.5*(wn(k,j,nh) + wn(k+1,j,nh)));
    s -= car.x3f(m,k,j,i)*(s3lo ? seam3(k, 0, nl) : 0.5*(wn(k-1,j,nl) + wn(k,j,nl)));
  }
  return s/cvol(m,k,j,i);
}

//----------------------------------------------------------------------------------------
//! \fn M1CsLatF
//! \brief the counterpart of M1SphLat (d = 1, 2) for the lateral face with
//! the planar unit normal nf: (1/r^3) d_r (r^3 (nf . v) E) at cell (k,j) of the column,
//! centred over i -+ 1 and one-sided at a physical x1 boundary exactly as M1SphLat.
//! nhat_f is a fixed Cartesian vector along the radial line (the face plane holds the
//! origin), so (div P).nhat_f = div(P nhat_f), and the r rhat-row part of it is this
//! term (the derivation in CS_VETGD.md sect. 4C: d_r(r^2 w_n)/r^2 + w_n/r).

template <class V, class A>
KOKKOS_INLINE_FUNCTION
Real M1CsLatF(const V &iw, const V &tt, const int c0, const A &x1v, const Real *nf,
              const int m, const int k, const int j, const int i, const int il,
              const int iu, const int ec) {
  const int ia = (i+1 <= iu) ? (i+1) : i;
  const int ib = (i-1 >= il) ? (i-1) : i;
  if (ia == ib) {return 0.0;}
  auto q = [&](const int iq) {
    const Real r = x1v(m,iq);
    return r*r*r*(tt(m,c0+1,k,j,iq)*nf[0] + tt(m,c0+2,k,j,iq)*nf[1]
                  + tt(m,c0+3,k,j,iq)*nf[2])*iw(m,ec,k,j,iq);
  };
  const Real r = x1v(m,i);
  return (q(ia) - q(ib))/(r*r*r*(x1v(m,ia) - x1v(m,ib)));
}

//----------------------------------------------------------------------------------------
//! \brief cs-tan-1106 (vet_gd_tangential on the cubed sphere): the trace-free tangential
//! part Q of the vet_gd tensor, the cs counterpart of M1SphTan.  Q . r_hat = 0, tr Q = 0,
//! so a cell holds two numbers in tau_ten slots c0+4, c0+5 about an orthonormal tangent
//! frame (e1, e2) of THAT cell (qf, Cartesian, (m, 6, k, j), ghosts exchanged like the
//! slots, so a seam ghost holds its mirror cell's frame next to its mirror cell's
//! values):
//!   Q = s1 (e1 e1 - e2 e2) + s2 (e1 e2 + e2 e1).
//! Contractions with fixed Cartesian vectors need no rotation at a seam or a vertex.  At
//! the face midpoint, with the face's planar normal n (Cartesian, constant on the face)
//! and the unit tangent t along the face, w = Q n E is a tangential field: w . r_hat = 0
//! on the radial line through the midpoint, so d_r(w . r_mid) = 0 there and
//!   n . div(Q E) = div(w) = d_n (n.Q.n E) + d_t (t.Q.n E),
//! with n and t FROZEN at their midpoint values over the stencil: the curvature terms of
//! M1SphTan (2 cot a E / r, ...) are in the variation of the Cartesian contractions.  The
//! normal derivative is the compact two-point one (the wide centred difference is what
//! drove the sp thin-top L_top runaway, fc693cbe1), corrected for the skew like the
//! isotropic gradient, d_n = (d_a - cos d_t)/sin; at a seam the mirror pair, moved from
//! the chord's foot to the face midpoint by M1CsSeamMid.

template <class VQ, class VF, class VS, class VG, class VB, class VN, class VX>
struct M1CsTanGeo {
  VQ qf;
  VF cdxf;
  VS sn2, cs2, sn3, cs3;
  VG g2, g3;
  VB cse;
  VN fn2, fn3, cn2, cn3;
  VX x1v;
  int js, je, ks, ke;
};

template <class VQ, class VF, class VS, class VG, class VB, class VN, class VX>
KOKKOS_INLINE_FUNCTION
M1CsTanGeo<VQ,VF,VS,VG,VB,VN,VX> M1MakeCsTanGeo(
    const VQ &qf, const VF &cdxf, const VS &sn2, const VS &cs2, const VS &sn3,
    const VS &cs3, const VG &g2, const VG &g3, const VB &cse, const VN &fn2,
    const VN &fn3, const VN &cn2, const VN &cn3, const VX &x1v, const int js,
    const int je, const int ks, const int ke) {
  return M1CsTanGeo<VQ,VF,VS,VG,VB,VN,VX>{qf, cdxf, sn2, cs2, sn3, cs3, g2, g3, cse, fn2,
                                          fn3, cn2, cn3, x1v, js, je, ks, ke};
}

// a . Q . b at cell (k,j,i)
template <class VT, class VQ>
KOKKOS_INLINE_FUNCTION
Real M1CsQab(const VT &tt, const int c0, const VQ &qf, const int m, const int k,
             const int j, const int i, const Real *a, const Real *b) {
  const Real a1 = a[0]*qf(m,0,k,j) + a[1]*qf(m,1,k,j) + a[2]*qf(m,2,k,j);
  const Real a2 = a[0]*qf(m,3,k,j) + a[1]*qf(m,4,k,j) + a[2]*qf(m,5,k,j);
  const Real b1 = b[0]*qf(m,0,k,j) + b[1]*qf(m,1,k,j) + b[2]*qf(m,2,k,j);
  const Real b2 = b[0]*qf(m,3,k,j) + b[1]*qf(m,4,k,j) + b[2]*qf(m,5,k,j);
  return tt(m,c0+4,k,j,i)*(a1*b1 - a2*b2) + tt(m,c0+5,k,j,i)*(a1*b2 + a2*b1);
}

// the face's n and t: an x2 face (e = 1) t along +x3 from the x3 plane through the face
// midpoint (its normal cn3(k) is tangential there), an x3 face (e = 2) t along +x2
template <class G>
KOKKOS_INLINE_FUNCTION
void M1CsFaceNT(const G &tg, const int e, const int m, const int k, const int j, Real *n,
                Real *t) {
  for (int c = 0; c < 3; ++c) {
    n[c] = (e == 1) ? tg.fn2(m,j,c) : tg.fn3(m,k,c);
    t[c] = (e == 1) ? tg.cn3(m,k,c) : tg.cn2(m,j,c);
  }
  const Real d = n[0]*t[0] + n[1]*t[1] + n[2]*t[2];
  for (int c = 0; c < 3; ++c) {t[c] -= d*n[c];}
  const Real nn = 1.0/sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);
  for (int c = 0; c < 3; ++c) {t[c] *= nn;}
}

// d/ds of the M1CsSeamMid quadratic at the face midpoint (per unit angle); 0 at a foot-
// value face (a vertex face under cs_vertex_mono, where Q is zero anyway)
template <class F, class G>
KOKKOS_INLINE_FUNCTION
Real M1CsSeamDs(const F &fv, const G &g, const int m, const int c, const int sd,
                const bool up, const bool dn) {
  int a, b, e;
  if (up && dn) {
    a = c - 1; b = c; e = c + 1;
  } else if (up) {
    a = c; b = c + 1; e = c + 2;
  } else if (dn) {
    a = c - 2; b = c - 1; e = c;
  } else {
    return 0.0;
  }
  const Real sa = g(m,2,a,sd), sb = g(m,2,b,sd), se = g(m,2,e,sd), s = g(m,3,c,sd);
  return fv(a)*(2.0*s - sb - se)/((sa - sb)*(sa - se))
         + fv(b)*(2.0*s - sa - se)/((sb - sa)*(sb - se))
         + fv(e)*(2.0*s - sa - sb)/((se - sa)*(se - sb));
}

//! \fn M1CsTanFaceLag: n . div(Q E) at face (e, k, j, i), the second-order lagged form
//! (E = component ec of iw).  vmo, vmn: cs_vertex_mono, as the isotropic gradient.
template <class V, class VT, class G>
KOKKOS_INLINE_FUNCTION
Real M1CsTanFaceLag(const V &iw, const VT &tt, const int c0, const G &tg, const int e,
                    const int m, const int k, const int j, const int i, const int ec,
                    const bool vmo, const int vmn, const Real xmono = -1.0) {
  Real n[3], t[3];
  M1CsFaceNT(tg, e, m, k, j, n, t);
  // cell of side s (0 = low, 1 = high) at position q along the face
  auto ck = [&](const int s, const int q) {return (e == 1) ? q : (k - 1 + s);};
  auto cj = [&](const int s, const int q) {return (e == 1) ? (j - 1 + s) : q;};
  auto pq = [&](const int s, const int q, const Real *a) {
    return M1CsQab(tt, c0, tg.qf, m, ck(s,q), cj(s,q), i, a, n)*
           iw(m,ec,ck(s,q),cj(s,q),i);
  };
  const int c = (e == 1) ? k : j;
  const int lo = (e == 1) ? tg.ks : tg.js, hi = (e == 1) ? tg.ke : tg.je;
  const bool slo = (e == 1) ? (tg.cse(m,2) != 0) : (tg.cse(m,0) != 0);
  const bool shi = (e == 1) ? (tg.cse(m,3) != 0) : (tg.cse(m,1) != 0);
  const int f = (e == 1) ? j : k;
  const int flo = (e == 1) ? tg.js : tg.ks, fhi = ((e == 1) ? tg.je : tg.ke) + 1;
  const int cl = (e == 1) ? 0 : 2;
  const bool seam = (f == flo && tg.cse(m,cl) != 0) || (f == fhi && tg.cse(m,cl+1) != 0);
  if (seam) {
    const int sd = (f == flo) ? 0 : 1;
    const Real r = tg.x1v(m,i);
    bool up = !(c == hi && shi), dn = !(c == lo && slo);
    if (vmo && M1CsSeamDepth(c, lo, hi, slo, shi) <= vmn) {up = false; dn = false;}
    if (e == 1) {
      auto dpr = [&](const int q) {
        return (pq(1,q,n) - pq(0,q,n))/(r*tg.g2(m,1,q,sd));
      };
      auto fv = [&](const int q) {return 0.5*(pq(1,q,t) + pq(0,q,t));};
      return M1CsSeamMid(dpr, tg.g2, m, c, sd, up, dn)
             + M1CsSeamDs(fv, tg.g2, m, c, sd, up, dn)/r;
    }
    auto dpr = [&](const int q) {
      return (pq(1,q,n) - pq(0,q,n))/(r*tg.g3(m,1,q,sd));
    };
    auto fv = [&](const int q) {return 0.5*(pq(1,q,t) + pq(0,q,t));};
    return M1CsSeamMid(dpr, tg.g3, m, c, sd, up, dn)
           + M1CsSeamDs(fv, tg.g3, m, c, sd, up, dn)/r;
  }
  auto der = [&](const int s, const Real *a) {
    auto q = [&](const int o) {return pq(s, c + o, a);};
    auto h = [&](const int o) {
      return (e == 1) ? tg.cdxf.x3f(m,k+o,cj(s,c),i) : tg.cdxf.x2f(m,ck(s,c),j+o,i);
    };
    return M1CsAxisDer(q, h, c, lo, hi, slo, shi, xmono);   // cs_cross_mono (seamE-1106)
  };
  const Real sn = (e == 1) ? tg.sn2(m,k,j) : tg.sn3(m,k,j);
  const Real cs = (e == 1) ? tg.cs2(m,k,j) : tg.cs3(m,k,j);
  const Real al = (e == 1) ? tg.cdxf.x2f(m,k,j,i) : tg.cdxf.x3f(m,k,j,i);
  return ((pq(1,c,n) - pq(0,c,n))/al - cs*0.5*(der(0,n) + der(1,n)))/sn
         + 0.5*(der(0,t) + der(1,t));
}

//! \fn M1CsTanFaceOp: the same term as frozen coefficients on the 19-point stencil,
//! vis(kq, jq, c) per cell of the column i (one body for VetLatOpCs and its stencil).  No
//! foot-to-midpoint move at a seam and a first-order one-sided d_t next to a seam
//! (neither fits the stencil): the operator is the implicit part only, the fixed point is
//! the lagged form above (the right-hand side adds it and removes this one at E^k).
template <class VT, class G, class F>
KOKKOS_INLINE_FUNCTION
void M1CsTanFaceOp(const VT &tt, const int c0, const G &tg, const int e, const int m,
                   const int k, const int j, const int i, const F &vis) {
  Real n[3], t[3];
  M1CsFaceNT(tg, e, m, k, j, n, t);
  auto ck = [&](const int s, const int q) {return (e == 1) ? q : (k - 1 + s);};
  auto cj = [&](const int s, const int q) {return (e == 1) ? (j - 1 + s) : q;};
  auto qq = [&](const int s, const int q, const Real *a) {
    return M1CsQab(tt, c0, tg.qf, m, ck(s,q), cj(s,q), i, a, n);
  };
  const int c = (e == 1) ? k : j;
  const int lo = (e == 1) ? tg.ks : tg.js, hi = (e == 1) ? tg.ke : tg.je;
  const bool slo = (e == 1) ? (tg.cse(m,2) != 0) : (tg.cse(m,0) != 0);
  const bool shi = (e == 1) ? (tg.cse(m,3) != 0) : (tg.cse(m,1) != 0);
  const int f = (e == 1) ? j : k;
  const int flo = (e == 1) ? tg.js : tg.ks, fhi = ((e == 1) ? tg.je : tg.ke) + 1;
  const int cl = (e == 1) ? 0 : 2;
  const bool seam = (f == flo && tg.cse(m,cl) != 0) || (f == fhi && tg.cse(m,cl+1) != 0);
  const int a = (c == hi && shi) ? c : c + 1, b = (c == lo && slo) ? c : c - 1;
  const Real r = tg.x1v(m,i);
  if (seam) {
    const int sd = (f == flo) ? 0 : 1;
    const Real pair = r*((e == 1) ? tg.g2(m,1,c,sd) : tg.g3(m,1,c,sd));
    vis(ck(1,c), cj(1,c), qq(1,c,n)/pair);
    vis(ck(0,c), cj(0,c), -qq(0,c,n)/pair);
    const Real ds = r*((e == 1) ? (tg.g2(m,2,a,sd) - tg.g2(m,2,b,sd))
                                : (tg.g3(m,2,a,sd) - tg.g3(m,2,b,sd)));
    for (int s = 0; s < 2; ++s) {
      vis(ck(s,a), cj(s,a), 0.5*qq(s,a,t)/ds);
      vis(ck(s,b), cj(s,b), -0.5*qq(s,b,t)/ds);
    }
    return;
  }
  const Real sn = (e == 1) ? tg.sn2(m,k,j) : tg.sn3(m,k,j);
  const Real cs = (e == 1) ? tg.cs2(m,k,j) : tg.cs3(m,k,j);
  const Real al = (e == 1) ? tg.cdxf.x2f(m,k,j,i) : tg.cdxf.x3f(m,k,j,i);
  vis(ck(1,c), cj(1,c), qq(1,c,n)/(al*sn));
  vis(ck(0,c), cj(0,c), -qq(0,c,n)/(al*sn));
  const Real cot = cs/sn;
  for (int s = 0; s < 2; ++s) {
    Real den = 0.0;
    if (a != c) {
      den += (e == 1) ? tg.cdxf.x3f(m,k+1,cj(s,c),i) : tg.cdxf.x2f(m,ck(s,c),j+1,i);
    }
    if (b != c) {
      den += (e == 1) ? tg.cdxf.x3f(m,k,cj(s,c),i) : tg.cdxf.x2f(m,ck(s,c),j,i);
    }
    vis(ck(s,a), cj(s,a), 0.5*(qq(s,a,t) - cot*qq(s,a,n))/den);
    vis(ck(s,b), cj(s,b), -0.5*(qq(s,b,t) - cot*qq(s,b,n))/den);
  }
}

//----------------------------------------------------------------------------------------
//! MILESTONE 3g, the FROZEN-DENSITY EOS CACHE, <rad_vet>/implicit_eos_cache.
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
//! 2*nt+1 consecutive cells centred on the cell of T^n (<rad_vet>/implicit_eos_cache_nt),
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
//! \struct M1EosIdeal
//! \brief the accessor for an ideal-gas EOS_Data (table off): e and c_v by the same
//! expressions as the ideal branch of EOS_Data::ThermoAt (bitwise), without the table
//! code in the kernel (m1-fast5-sp).

struct M1EosIdeal {
  Real gamma;
  KOKKOS_INLINE_FUNCTION
  void operator()(const Real dd, const Real t, Real &ee, Real &cv) const {
    ee = dd*t/(gamma-1.0);
    cv = 1.0/(gamma-1.0);
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
//! \fn M1ImplTemperatureOpacNear
//! \brief <rad_vet>/implicit_tsolve_opac_mode bit 1 (tsolve_root_fix_1003): the same
//! equation as M1ImplTemperatureOpac, but the root NEAREST the guess in log T.  y(T) need
//! not be monotone where kappa_P(T) T^4 falls with T, and the halve/double bracket of
//! M1ImplTemperatureOpac is grown in the direction a MONOTONE y would point to, so it can
//! step over the root next to the iterate and bisect onto a distant one.  Here the scan
//! steps outward from the guess by the factor M1_TSO_NEAR_Q alternately up and down, so
//! the first sign change met is the closest one (to within one step), then bisects that
//! cell in log T (either orientation of the crossing).  No sign change within a factor
//! 2^(M1_TSO_NEAR_NS/16) of the guess: ok = false, tout = guess.

constexpr int  M1_TSO_NEAR_NS = 320;                 // steps per side (factor 2^20)
constexpr Real M1_TSO_NEAR_Q  = 1.0442737824274138;  // 2^(1/16)

template <class ThermoT, class OpF>
KOKKOS_INLINE_FUNCTION
int M1ImplTemperatureOpacNear(const ThermoT &th, const OpF &opf, const Real dd,
                              const Real tguess, const Real egn, const Real cdta,
                              const Real cdt, const Real erhs, Real &tout, bool &ok) {
  auto yfn = [&](const Real t) {
    Real ee, cv, rkp, rke;
    th(dd, t, ee, cv);
    opf(dd, t, rkp, rke);
    const Real t2 = t*t;
    return ee + cdta*rkp*t2*t2 - egn - cdt*rke*erhs;
  };
  const Real y0 = yfn(tguess);
  int nev = 1;
  ok = true;
  tout = tguess;
  if (y0 == 0.0) {return nev;}
  const bool s0 = (y0 > 0.0);
  Real tu = tguess, td = tguess;      // the scan fronts above and below the guess
  Real ta = tguess, tb = tguess;      // the bracket: sign(y(ta)) = s0 != sign(y(tb))
  bool found = false;
  for (int s=0; s<M1_TSO_NEAR_NS && !found; ++s) {
    const Real t1 = tu*M1_TSO_NEAR_Q;
    const Real y1 = yfn(t1);
    ++nev;
    if (y1 == 0.0) {
      tout = t1;
      return nev;
    }
    if ((y1 > 0.0) != s0) {
      ta = tu;
      tb = t1;
      found = true;
    } else {
      tu = t1;
      const Real t2 = td/M1_TSO_NEAR_Q;
      const Real y2 = yfn(t2);
      ++nev;
      if (y2 == 0.0) {
        tout = t2;
        return nev;
      }
      if ((y2 > 0.0) != s0) {
        ta = td;
        tb = t2;
        found = true;
      } else {
        td = t2;
      }
    }
  }
  if (!found) {
    ok = false;
    return nev;
  }
  // bisection in log T between ta (sign s0) and tb (the other sign), either order
  for (int it=0; it<M1_IMPL_TMAXIT && fabs(tb - ta) > M1_IMPL_TRTOL*fmax(ta, tb); ++it) {
    const Real tm = sqrt(ta*tb);
    const Real y = yfn(tm);
    ++nev;
    if ((y > 0.0) == s0) {ta = tm;} else {tb = tm;}
  }
  tout = sqrt(ta*tb);
  return nev;
}

//----------------------------------------------------------------------------------------
//! \fn M1ImplTemperatureOpacStable / M1ImplTemperatureOpacOde
//! \brief <rad_vet>/implicit_tsolve_opac_mode bit 4 (value 8, r2b_shift_1003): the same
//! equation as M1ImplTemperatureOpac, but the STABLE root (up-crossing of y in T) nearest
//! the END POINT of the exact exchange ODE de/dt = c (rho kE E' - rho kP(T) a T^4) at
//! fixed E' over dt, started from the no-exchange T_gn (e(T_gn) = e_gn).  Why: every
//! root of the backward-Euler gas equation lies between T_gn and T_rad; where kappa_P
//! rises steeply with T (H recombination, absorbing gas) it has several.  Neither the
//! root nearest the lagged T (bit 1: lag-dependent, 2-cycles) nor the one nearest T_rad
//! (tried as value 4: wrong when the exchange is not stiff, BSG R2N diverged) is the
//! physical one; the ODE end point is, and it is lag-independent and continuous in E'.
//! BSG cell rank 2 (0,38,44,181), cycle 14572: ODE ends at T_rad = 11241 K, root
//! 11232 K (not 9423.6 K); at dt/10 it ends at 8897 K, root 8020 K (value 4: 11119 K).
//! T_gn: safeguarded Newton on e(T) = e_gn from the lagged T.  End point: elapsed time
//! int de/rate on M1_TSO_ODE_NS (32) log-T points from T_gn toward T_rad (trapezoid),
//! stopped where it reaches c dt.  The root: the tightened outward log scan of
//! tsolve-near-tight 09b12e4e started there (first step ln2/1024, x1.25 per step, up to
//! 2^20, nearer of the brackets of the first bracketing step, bisected in log T).
//! No stable root, E' + de0 <= 0 or no T_gn: ok = false, tout = the lagged T.

constexpr int  M1_TSO_RAD_NS = 64;                  // max steps per side
constexpr Real M1_TSO_RAD_H0 = 6.7690286301289e-4;  // first step in ln T (ln 2 / 1024)
constexpr Real M1_TSO_RAD_G  = 1.25;                // step growth per step
constexpr Real M1_TSO_RAD_LMAX = 13.862943611198906;  // ln(2^20)

template <class ThermoT, class OpF>
KOKKOS_INLINE_FUNCTION
int M1ImplTemperatureOpacStable(const ThermoT &th, const OpF &opf, const Real dd,
                              const Real tguess, const Real egn, const Real cdta,
                              const Real cdt, const Real erhs, Real &tout, bool &ok) {
  auto yfn = [&](const Real t) {
    Real ee, cv, rkp, rke;
    th(dd, t, ee, cv);
    opf(dd, t, rkp, rke);
    const Real t2 = t*t;
    return ee + cdta*rkp*t2*t2 - egn - cdt*rke*erhs;
  };
  const Real y0 = yfn(tguess);
  int nev = 1;
  ok = true;
  tout = tguess;
  if (y0 == 0.0) {return nev;}
  // scan fronts (T and y) above and below the guess; an up-crossing bracket [lo, hi] has
  // y(lo) < 0 <= y(hi)
  Real tu = tguess, yu = y0, td = tguess, yd = y0;
  Real lu = 0.0, hu = 0.0, ld = 0.0, hd = 0.0;
  bool fu = false, fd = false;
  Real dl = 0.0, h = M1_TSO_RAD_H0;
  for (int s=0; s<M1_TSO_RAD_NS && !fu && !fd && dl < M1_TSO_RAD_LMAX; ++s) {
    dl += h;
    h *= M1_TSO_RAD_G;
    const Real t1 = tguess*exp(dl);
    const Real y1 = yfn(t1);
    const Real t2 = tguess*exp(-dl);
    const Real y2 = yfn(t2);
    nev += 2;
    if (yu < 0.0 && y1 >= 0.0) {
      fu = true;
      lu = tu;
      hu = t1;
    }
    if (y2 < 0.0 && yd >= 0.0) {
      fd = true;
      ld = t2;
      hd = td;
    }
    tu = t1;
    yu = y1;
    td = t2;
    yd = y2;
  }
  if (!fu && !fd) {
    ok = false;
    return nev;
  }
  // bisection in log T of each bracket found (y(lo) < 0 <= y(hi))
  Real ru = tguess, rd = tguess;
  if (fu) {
    for (int it=0; it<M1_IMPL_TMAXIT && (hu - lu) > M1_IMPL_TRTOL*hu; ++it) {
      const Real tm = sqrt(lu*hu);
      const Real y = yfn(tm);
      ++nev;
      if (y < 0.0) {lu = tm;} else {hu = tm;}
    }
    ru = sqrt(lu*hu);
  }
  if (fd) {
    for (int it=0; it<M1_IMPL_TMAXIT && (hd - ld) > M1_IMPL_TRTOL*hd; ++it) {
      const Real tm = sqrt(ld*hd);
      const Real y = yfn(tm);
      ++nev;
      if (y < 0.0) {ld = tm;} else {hd = tm;}
    }
    rd = sqrt(ld*hd);
  }
  if (fu && fd) {
    tout = (fabs(log(ru/tguess)) <= fabs(log(rd/tguess))) ? ru : rd;
  } else {
    tout = fu ? ru : rd;
  }
  return nev;
}

constexpr int M1_TSO_ODE_NS = 32;    // log-T points of the exchange-ODE quadrature

} // namespace radvet
#endif // RAD_VET_RAD_VET_IMPLICIT_HPP_
