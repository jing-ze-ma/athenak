#ifndef RAD_VET_RAD_VET_IMPLICIT_HPP_
#define RAD_VET_RAD_VET_IMPLICIT_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_implicit.hpp
//! \brief constants (key values, work-array layouts, restart block markers) and
//! device-inline helpers of the implicit radiation transport.

#include <math.h>

#include "athena.hpp"
#include "eos/eos_table.hpp"

namespace radvet {

// Restart block markers (8 bytes).  Each block is written behind the previous one only
// when its state exists; a file without it restarts from the cold state.
// implicit_predictor = step: int32 have, int32 ncomp, Real pred_dt; ncomp slabs of ipred
constexpr char kM1PredRstMagic[8] = {'M', '1', 'P', 'R', 'E', 'D', '0', '1'};

// time_scheme
constexpr int M1_TIME_BE = 0;        // Lie splitting: Heun hydro, then one BE solve
constexpr int M1_TIME_HESDIRK2 = 1;  // H-ESDIRK2: one stage solve inside each Heun stage
// t2_solve: what the next ImplicitSolve does
constexpr int M1_T2S_NONE = 0;     // plain backward Euler
constexpr int M1_T2S_BESTORE = 1;  // backward Euler, K1 = (Y - rhs)/dt stored
constexpr int M1_T2S_STAGE1 = 2;   // stage solve 1 (gamma dt), K2 stored
constexpr int M1_T2S_STAGE2 = 3;   // stage solve 2 (gamma dt), K1 (FSAL) stored
// channels of t2k1 / t2k2 / t2inc: E, gas momentum, total energy, face fluxes F0 x1..x3
constexpr int M1_T2_E = 0;
constexpr int M1_T2_M1 = 1;
constexpr int M1_T2_EN = 4;
constexpr int M1_T2_F1 = 5;
constexpr int M1_T2_NK = 8;
// vet_sc: extrapolated components of vet_cell (M1_VET_CHI .. M1_VET_D11+5)
constexpr int M1_T2_NVET = 10;
// hesdirk2 slope: int32 have, int32 nch, Real pred2_dt, dt_prev, vprev; nch slabs
// (K1, ipred2, vet_prev) behind the predictor slabs
constexpr char kM1Time2RstMagic[8] = {'M', '1', 'T', 'I', 'M', 'E', '2', 'A'};
// implicit_one_pass > 0: 9 Reals (last two contractions, solve counters)
constexpr char kM1OnePassRstMagic[8] = {'M', '1', 'O', 'N', 'E', 'P', '0', '1'};
// implicit_mr_every / vet_sc_every > 1: multi-rate window and tensor cadence, 10 Reals
constexpr char kM1MRRstMagic[8] = {'M', '1', 'M', 'R', 'W', 'I', 'N', '1'};
// time2_auto_be back-off while active: int32 len, int32 left
constexpr char kM1AutoBeRstMagic[8] = {'M', '1', 'A', 'U', 'T', 'O', 'B', 'E'};
// implicit_closure_thin_relax: ctr_mem (chi, n1, n2, n3): int32 have, nch = 4; slabs last
constexpr char kM1CtrRstMagic[8] = {'M', '1', 'C', 'T', 'R', 'L', 'X', '1'};
// rad_signal_speed inputs: int32 have = 1, nch, mode, 0; nch slabs (opac T, tensor)
constexpr char kM1RssRstMagic[8] = {'M', '1', 'R', 'S', 'S', 'I', 'G', '1'};
// vet_gd_async pending source (ln chi, ln S): int32 have = 1, nch = 2, ncycle, 0
constexpr char kM1VgdRstMagic[8] = {'M', '1', 'V', 'G', 'D', 'A', 'S', '1'};
// conservative gas source with hesdirk2: int32 have = 1, nch (3 | 0), csk1_ok, 0; nch
// slabs of csk1
constexpr char kM1CskRstMagic[8] = {'M', '1', 'C', 'S', 'K', '1', 'S', '1'};
// time2_face_fct: face fluxes t2pf1 (FaxRst*)
constexpr char kM1FaxRstMagic[8] = {'M', '1', 'F', 'A', 'X', 'R', 'S', '1'};
// vet_rotate_avg history: M1_RA_RST_NH int64 words (RotAvgRstHdr), then nslab slabs
constexpr char kM1RaRstMagic[8] = {'M', '1', 'R', 'O', 'T', 'A', 'V', '1'};
// old files: stage-2 drop slabs (skipped)
constexpr char kM1ThcRstMagic[8] = {'M', '1', 'T', 'H', 'C', 'R', 'S', '1'};
// old files: cross-fade state behind the thc slabs (skipped)
constexpr char kM1XfdRstMagic[8] = {'M', '1', 'X', 'F', 'D', 'R', 'S', '1'};
constexpr int kM1XfdRstNi = 24;  // int64 header words
constexpr int kM1XfdRstNr = 8;   // Real header words

// implicit_solver (both converge to the same fixed point)
constexpr int M1_ISOLV_LINE_JACOBI = 0;  // x1 tridiagonal exact, x2/x3 lagged (M-matrix)
constexpr int M1_ISOLV_BICGSTAB = 1;  // BiCGStab, x1 line solve as right preconditioner

// implicit_offdiag: the off-diagonal Eddington terms sum_{e != d} d_e P_de (lagging them
// diverges in thin cells, where they scale as (c dt/dx)^2)
constexpr int M1_OD_LAGGED   = 0;  // previous Picard iterate, on the right-hand side
constexpr int M1_OD_OPERATOR = 1;  // in the operator (bicgstab; not an M-matrix)
constexpr int M1_OD_NONE     = 2;  // dropped: P diagonal in the grid frame

// implicit_trans_limit: transverse realizability limiter.  lp adds the lagged opacity
// |G_f|/(phi_f E_f), phi_f = fmax sqrt(max(1 - f1^2, 0.01)), so |F_f| <= c phi_f E_f
constexpr int M1_TLIM_NONE = 0;
constexpr int M1_TLIM_LPR  = 2;  // bound max(kt, |G|/(fmax E)), tapered by min(1, kt dx)
constexpr int M1_TLIM_LP   = 1;  // lagged limiter opacity (above)

// boundary types of the implicit solve, as the total normal flux F_f at the face
constexpr int M1_IBC_MARSHAK  = 0;  // free surface, F_f = +- c marshak_q E
constexpr int M1_IBC_FLUX     = 1;  // imposed F_f = implicit_flux_x1min/max
constexpr int M1_IBC_REFLECT  = 2;  // F_f = 0
constexpr int M1_IBC_PERIODIC = 3;  // cyclic tridiagonal (Sherman-Morrison)
constexpr int M1_IBC_EFIX     = 4;  // Dirichlet: boundary cell E frozen at its start
constexpr int M1_IBC_OUTFLOW  = 5;  // x2: F_f = max(+-F2_cell, 0), lagged
constexpr int M1_IBC_X2FLUX   = 6;  // x2: F_f = ifx2_lo/hi from the pgen (SetX2FaceFlux)

// implicit_accel: acceleration of the outer Picard iteration
constexpr int M1_IACC_NONE     = 0;   // bare Picard map
constexpr int M1_IACC_ANDERSON = 1;   // Anderson acceleration (Walker & Ni 2011)
constexpr Real M1_AND_REG = 1.0e-10;  // Tikhonov weight, relative to tr(dG^T dG)/m
constexpr int M1_AND_MMAX = 10;  // cap on implicit_anderson_m (host solve, one GlobalSum)

// implicit_flux: spatial form of the implicit face E-flux
constexpr int M1_IFLUX_CENTRAL = 0;  // face-eliminated central (diffusion) form
constexpr int M1_IFLUX_APHLL   = 1;  // A_up + alpha F_HLL + (1-alpha) F_diff, HLL lagged
constexpr int M1_IFLUX_BERTHON = 2;  // A_up + alpha F_HLL (Berthon), no F_diff
constexpr int M1_IFLUX_BLEND   = 3;  // (1 - w_f) central + w_f berthon, w_f lagged

// implicit_blend: the face weight w_f of M1_IFLUX_BLEND
constexpr int M1_IBLEND_TAU  = 0;   // exp(-(tau_face/tau0)^2): upwind where thin
constexpr int M1_IBLEND_F    = 1;   // smoothstep of the lagged face reduced flux
constexpr int M1_IBLEND_IDF = 5;    // idort x smoothstep(|H|/J)
constexpr int M1_IBLEND_IDA = 6;    // idort x smoothstep(h asymmetry)
constexpr int M1_IBLEND_KN = 4;     // half-range only: Knudsen x tau
constexpr int M1_IBLEND_IDORT = 3;  // half-range: 1/(1 + x + x^2/tau0), x = alpha tau
constexpr int M1_IBLEND_TAUF = 2;   // TAU x F

// implicit_thin_switch: dt-free face weight of the half-range blend (M1ThinFace)
constexpr int M1_THSW_NONE = 0;
constexpr int M1_THSW_HJ   = 1;  // smoothstep of the formal-solution |H|/J
constexpr int M1_THSW_KN   = 2;  // face Knudsen number of E^n
constexpr int M1_THSW_FS   = 3;  // smoothstep of max(|H|/J, C1/1.5), C1 = <|I-S|>/J
constexpr int M1_THSW_HJKN = 4;  // max(hj, kn)
constexpr int M1_THSW_FSKN = 5;  // max(fs, kn)
constexpr int M1_THSW_FSI  = 7;  // fs; h+- isotropic where the rays miss E^n (M1ThinHr)
constexpr int M1_THSW_FSE  = 8;  // max(|H|/J, eps C1/1.5), eps the thermal fraction
constexpr int M1_THSW_FSQ  = 6;  // max(|H|/J, q_c C1/1.5), q_c = min(E^n/J, J/E^n)^qp
// implicit_thin_corr: lagged formal-solution face flux
constexpr int M1_THCR_NONE = 0;
constexpr int M1_THCR_LAG  = 1;  // + w (F_sc^n - F_hr(E^n)), a fixed source
constexpr int M1_THCR_SC   = 2;  // w F_sc^n replaces the half-range part
constexpr int M1_THCR_HCS  = 4;  // hc weighted by the scattering fraction 1 - eps
constexpr int M1_THCR_HC   = 3;  // + w (F_sc,central - F_sc,upwind)^n

// implicit_blend_fmode: which cell's lagged reduced flux sets f_face
constexpr int M1_IBFM_MAX  = 0;
constexpr int M1_IBFM_MEAN = 1;

// implicit_blend_mode: what is blended
constexpr int M1_IBMODE_FLUX   = 0;  // the whole face flux
constexpr int M1_IBMODE_DISSIP = 1;  // central + w alpha dk (E_L - E_R)

// implicit_recon: states of the HLL part
constexpr int M1_IRECON_DC    = 0;  // piecewise constant
constexpr int M1_IRECON_PLMDC = 1;  // plm as a lagged deferred correction (rhs)

// implicit_enthalpy: face value of the enthalpy flux a_d E, a_d = v_d + (v.D)_d.  The
// matrix carries the donor form; the others add their lagged difference to the rhs.
constexpr int M1_IENTH_UPWIND  = 0;  // donor cell (first order)
constexpr int M1_IENTH_CENTRAL = 1;  // a_f (E_L + E_R)/2, a_f = (a_L + a_R)/2
constexpr int M1_IENTH_PLM     = 2;  // a_f E_f, E_f the upwind van Leer plm value
constexpr int M1_IENTH_PLMC    = 3;  // mean of the two plm face values
constexpr int M1_IENTH_LIN     = 4;  // unlimited upwind linear (Fromm) value

// implicit_vimp: gas velocity implicit in the enthalpy flux through the radiative force
// of the same solve (Jacobian P of dv in E'; couples cells two apart).  Components are
// appended to iw at iw_vimp; the first M1_NVIMP_X are exchanged once per pass.
constexpr int M1_NVIMP_X  = 12;  // P1m,P10,P1p, P2m,P20,P2p, P3m,P30,P3p, DV1,DV2,DV3
constexpr int M1_IV_DV    = 9;
constexpr int M1_IV_X1M2  = 12;  // coefficients outside the 7-point row:
constexpr int M1_IV_X1P2  = 13;  // x1 +-2; x2 -2,-1,+1,+2; x3 -2,-1,+1,+2
constexpr int M1_IV_X2M2  = 14;
constexpr int M1_IV_X3M2  = 18;
constexpr int M1_IV_JD    = 22;  // added to the row: diagonal, x1 -1, x1 +1
constexpr int M1_IV_J1M   = 23;
constexpr int M1_IV_J1P   = 24;
constexpr int M1_IV_JRHS  = 25;  // added to the right-hand side
constexpr int M1_NIW_VIMP = 26;
// hesdirk2: DA_d = a_d(old) - a_d(stage start) = dv_d + (D dv)_d, behind the block
constexpr int M1_IV_DA    = 26;
constexpr int M1_NIW_VIMP_T2 = 3;
// implicit_hr_recon = plm: half-range faces E_L + s_L/2, E_R - s_R/2, s_c = sig_c
// (E_c+1 - E_c-1)/2 with the lagged van Leer limiter sig_c (M1MusclRow)
constexpr int M1_IM_SIG   = 0;  // sig_1, sig_2, sig_3
constexpr int M1_IM_D     = 3;  // row coefficients: diagonal
constexpr int M1_IM_X1    = 4;  // -2, -1, +1, +2
constexpr int M1_IM_X2    = 8;
constexpr int M1_IM_X3    = 12;
constexpr int M1_IM_KILL  = 16;  // 1: dc for the rest of the step (positivity)
constexpr int M1_IM_B     = 17;  // led: b_1, b_2, b_3 (SIG holds a_d)
constexpr int M1_IM_EPREV = 20;  // qs_mode = step: E^n of the previous step
constexpr int M1_NIW_MUSCL = 21;

// implicit_partition: columns that span several MeshBlocks
constexpr int M1_IPART_NONE   = 0;  // one MeshBlock per column
constexpr int M1_IPART_GATHER = 1;  // gathered to one rank, one Thomas sweep (bitwise)

// vet-positivity counters (pos_cnt_d / pos_cnt); energies volume-weighted
constexpr int M1_POS_G0     = 0;  // cell-passes whose lagged g0 was clipped
constexpr int M1_POS_GAS    = 1;  // cell-solves where the gas eint limiter fired
constexpr int M1_POS_GAS_DE = 2;  // energy moved gas <- radiation
constexpr int M1_POS_FLR    = 3;  // cell-solves where E was raised to e_floor
constexpr int M1_POS_FLR_DE = 4;  // floor energy taken from the gas
constexpr int M1_POS_FLR_UN = 5;  // floor energy not covered by the gas (created)
constexpr int M1_POS_FCLIP  = 6;  // cell-solves whose |F| > c E was scaled back
constexpr int M1_POS_FCLIPM = 7;  // sum of (|F|/(c E) - 1) over those
constexpr int M1_POS_FLRS   = 8;  // residual raise below 1e-6 e_floor
constexpr int M1_POS_FLR_MV = 9;  // floor energy taken from the column's radiation
constexpr int M1_POS_N      = 10;

// implicit face work array ifw (m,n,k,j,i), i the x1 face between cells i-1 and i
constexpr int M1_IFW_AL  = 0;  // alpha, the AP weight of F_HLL
constexpr int M1_IFW_HCL = 1;  // alpha x coefficient of E'_L in F_HLL, >= 0
constexpr int M1_IFW_HCR = 2;  // alpha x coefficient of E'_R in F_HLL, <= 0
constexpr int M1_IFW_DG  = 3;  // plm deferred correction alpha (F^plm - F^dc)
constexpr int M1_NIFW = 4;

// implicit work array iw (m,n,k,j,i): per solve, except the start-of-step EN / EGN
constexpr int M1_IW_EN   = 0;   // E^n
constexpr int M1_IW_EP   = 1;   // Picard iterate of E'
constexpr int M1_IW_TP   = 2;   // Picard iterate of T'
constexpr int M1_IW_EGN  = 3;   // rho e^n, gas internal energy density
constexpr int M1_IW_WCHI = 4;   // w = P_11/E of the lagged closure (chi in 1-D)
constexpr int M1_IW_ADV  = 5;   // a = v1 (1 + w), A = a E
constexpr int M1_IW_DE0  = 6;   // E0 - E, the O(beta) comoving correction
constexpr int M1_IW_G0   = 7;   // g0 = rho (kappa_E E0 - kappa_P a T^4), lagged
constexpr int M1_IW_V1   = 8;   // x1 gas velocity
constexpr int M1_IW_SRCB = 9;   // emission / absorption: diagonal
constexpr int M1_IW_SRCR = 10;  // emission / absorption: rhs
constexpr int M1_IW_TA   = 11;  // tridiagonal: lower
constexpr int M1_IW_TB   = 12;  // diagonal
constexpr int M1_IW_TC   = 13;  // upper
constexpr int M1_IW_TR   = 14;  // rhs
constexpr int M1_IW_S1   = 15;  // Thomas scratch (c')
constexpr int M1_IW_S2   = 16;  // Thomas solution
constexpr int M1_IW_S3   = 17;  // second solution (Sherman-Morrison)
constexpr int M1_IW_RES  = 18;  // per-cell Picard residual
constexpr int M1_IW_F1   = 19;  // cell-centred x1 flux of the iterate
constexpr int M1_IW_RF0  = 20;  // comoving F0/(c E) in [-1,1] (lagged by HLL)
constexpr int M1_IW_KT   = 21;  // copy of opac(M1_OP_T), for the x1 halo exchange
constexpr int M1_NIW_X1 = 22;   // x1-only solve
constexpr int M1_IW_V2   = 22;  // transverse: x2 gas velocity
constexpr int M1_IW_V3   = 23;  // x3 gas velocity
constexpr int M1_IW_N1   = 24;  // lagged flux direction n = F/|F|
constexpr int M1_IW_N2   = 25;
constexpr int M1_IW_N3   = 26;
constexpr int M1_IW_A2   = 27;  // a_2 = v_2 + (v.D)_2
constexpr int M1_IW_A3   = 28;  // a_3 = v_3 + (v.D)_3
constexpr int M1_IW_TDIA = 29;  // diagonal part of the transverse couplings
constexpr int M1_IW_TRHS = 30;  // -(T(E^k) - TDIA E^k), lagged, to the rhs
constexpr int M1_IW_LRES = 31;  // |U(E^{k+1}) - U(E^k)|, 7-point system residual
constexpr int M1_IW_F2   = 32;  // cell-centred x2 flux of the iterate
constexpr int M1_IW_F3   = 33;  // cell-centred x3 flux
constexpr int M1_NIW = 34;

//----------------------------------------------------------------------------------------
//! \fn M1AccComp
//! \brief iw component of entry c of the Anderson vector (E, F_1, F_2, F_3)

KOKKOS_INLINE_FUNCTION
int M1AccComp(const int c) {
  switch (c) {
    case 0: return M1_IW_EP;
    case 1: return M1_IW_F1;
    case 2: return M1_IW_F2;
    default: return M1_IW_F3;
  }
}

// implicit_solver = bicgstab: transverse off-diagonals of the frozen 7-point row
//   TA E_{i-1} + TB E_i + TC E_{i+1} + CJM E_{j-1} + CJP E_{j+1} + CKM E_{k-1}
//     + CKP E_{k+1} = b,
// TB including TDIA and b = TR + sum C E^k (the line-Jacobi system), and Krylov vectors
constexpr int M1_IW_CJM  = 34;  // E_{j-1} (0 at a physical x2 face)
constexpr int M1_IW_CJP  = 35;  // E_{j+1}
constexpr int M1_IW_CKM  = 36;  // E_{k-1} (0 in 2-D)
constexpr int M1_IW_CKP  = 37;  // E_{k+1}
constexpr int M1_IW_KB   = 38;  // b
constexpr int M1_IW_KX   = 39;  // iterate x
constexpr int M1_IW_KR   = 40;  // residual r
constexpr int M1_IW_KRH  = 41;  // shadow residual rhat
constexpr int M1_IW_KP   = 42;  // direction p
constexpr int M1_IW_KV   = 43;  // v = A M^{-1} p
constexpr int M1_IW_KS   = 44;  // s = r - alpha v
constexpr int M1_IW_KTT  = 45;  // t = A M^{-1} s
constexpr int M1_IW_KY   = 46;  // y = M^{-1} p
constexpr int M1_IW_KZ   = 47;  // z = M^{-1} s
constexpr int M1_NIW_K = 48;

// implicit_hr_face = sf: per-cell formal-solution half-range fluxes (M1HrSf)
constexpr int M1_NIW_HSF = 8;
constexpr int M1_HSF_HP = 0;  // + 2 d: H+_d = sum w I max(mu_d, 0)
constexpr int M1_HSF_HM = 1;  // + 2 d: H-_d = H_d - H+_d <= 0
constexpr int M1_HSF_J  = 6;  // J of the formal solution
constexpr int M1_HSF_JH = 7;  // Jh = max(J, outgoing face SC fluxes)
// implicit_gas_newton / implicit_eos_cache: gas block, appended at offset iw_gas
constexpr int M1_NIW_GAS = 5;
constexpr int M1_IWG_BK = 0;  // B_k = rho c_v + 4 c dt rho kappa_P a T_k^3
constexpr int M1_IWG_RK = 1;  // R_k = rho (e^n - e(T_k)) - c dt rho (kP a T_k^4 - kE de0)
constexpr int M1_IWG_YR = 2;  // |R_k + c dt rho kappa_E E'| of the last pass (0: exact)
constexpr int M1_IWG_FB = 3;  // per-cell Newton fallbacks this step
constexpr int M1_IWG_MS = 4;  // per-cell EOS-cache misses this step

// Newton T step trust region (fraction of T_k); larger steps use the bracketed root find
constexpr Real M1_NEWT_TRUST = 0.5;
// implicit_opac_newton: relative T step of the opacity derivative
constexpr Real M1_OPN_H = 1.0e-5;

//----------------------------------------------------------------------------------------
//! \fn M1OpnCell
//! \brief implicit_opac_newton: one cell's share of a face's Newton term (q = -chat dt
//! th G(E^k)/2, kt = d(rho kappa_T)/dT); E' coefficient to coef, the rest to rr

KOKKOS_INLINE_FUNCTION
void M1OpnCell(const Real q, const Real kt, const Real bk, const Real rk, const Real ke,
               const Real cl, const Real dt, Real &coef, Real &rr) {
  if (bk > 0.0) {
    const Real s = q*kt/bk;
    coef += s*cl*dt*ke;
    rr -= s*rk;
  }
}
// BiCGStab breakdown: |rho|, |rhat.v| below this times the rhs scale
constexpr Real M1_BCG_EPS = 1.0e-300;
constexpr int M1_KD_SIZE = 20;  // implicit_krylov_dev: size of its device scalar array

//----------------------------------------------------------------------------------------
//! \fn M1EddDiag
//! \brief diagonal Eddington component D_dd of the lagged closure,
//! D_ab = (1-chi)/2 delta_ab + (3 chi - 1)/2 n_a n_b

KOKKOS_INLINE_FUNCTION
Real M1EddDiag(const Real chi, const Real nd) {
  return 0.5*(1.0 - chi) + 0.5*(3.0*chi - 1.0)*nd*nd;
}

//----------------------------------------------------------------------------------------
//! \fn M1CsSeamGeom
//! \brief cubed sphere: canonical geometry of the seam face next to along-seam cell q
//! (0..n-1): dth its seam arc, angm the angle between the cell centre and its mirror
//! across the seam, sfoot / smid the along-seam positions of the pair chord's foot and
//! of the face midpoint.  A function of (q, n) only: bitwise equal on both panels.

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
//! \fn M1CsSeamDepth
//! \brief depth of index c below the nearest panel-seam block face (slo / shi): 0 next to
//! the seam and in the first ghost beyond it; large when neither face is a seam

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

//----------------------------------------------------------------------------------------
//! \fn M1CsVtxDist
//! \brief Chebyshev distance of cell (k,j) to a cube vertex of its block

KOKKOS_INLINE_FUNCTION
int M1CsVtxDist(const int k, const int j, const int ks, const int ke, const int js,
                const int je, const bool s2lo, const bool s2hi, const bool s3lo,
                const bool s3hi) {
  const int dj = M1CsSeamDepth(j, js, je, s2lo, s2hi);
  const int dk = M1CsSeamDepth(k, ks, ke, s3lo, s3hi);
  return (dj > dk) ? dj : dk;
}

//----------------------------------------------------------------------------------------
//! \fn M1CsSeamMid
//! \brief seam quantity at the face midpoint (smid, slot 3 of g) from its values fv at
//! the chord feet (sfoot, slot 2) of faces c-1, c, c+1: the quadratic through them,
//! one-sided at a cube vertex (up / dn false); O(h^3)

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
//! \fn M1CsTanMono
//! \brief cs_cross_mono: clamps the one-sided seam derivative d to [0, 2 s1] where
//! |q0 - 2 q1 + q2| > fk (|q0| + |q1| + |q2|)/3; fk > 0 leaves smooth extrema alone

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

//----------------------------------------------------------------------------------------
//! \fn M1CsAxisDer
//! \brief derivative of q(o) along a tangential axis (h(o): arc between o-1 and o):
//! centred, second-order one-sided next to a panel seam (M1CsTanMono if mono >= 0)

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
//! \struct M1Vtx
//! \brief cs_vertex_3cell: a seam ghost next to a cube-vertex corner differences along
//! its seam one-sided (the corner ghost block has no physical counterpart)

struct M1Vtx {
  bool on = false;
  int js = 0, je = 0, ks = 0, ke = 0;
  bool jl = false, jh = false, kl = false, kh = false;  // panel edges of the block
};

//----------------------------------------------------------------------------------------
//! \fn M1SphDrr
//! \brief spherical polar: radial face weight with the integrating factor,
//! (1-chi)/2 + (3 chi-1)/2 n_r^2 s2, s2 = (r_c/r_f)^2, so the two-point difference of
//! D E gives d_r p + (1/r^2) d_r(r^2 q n_r^2).  Both parts >= 0 (M-matrix); chi = 1,
//! n = r keeps r^2 E = const exactly.

KOKKOS_INLINE_FUNCTION
Real M1SphDrr(const Real chi, const Real n1, const Real s2) {
  return 0.5*(1.0 - chi) + 0.5*(3.0*chi - 1.0)*n1*n1*s2;
}

//----------------------------------------------------------------------------------------
//! \fn M1SphMarshakCoef
//! \brief implicit_marshak_face = linear: Marshak face E_lin = ca ec - cb en from linear
//! extrapolation of r^2 E (end cell rc, neighbour rn, face rf), al = (rf - rc)/(rc - rn),
//! ca = (1 + al)(rc/rf)^2, cb = al (rn/rf)^2 > 0 (M-matrix sign pattern kept)

KOKKOS_INLINE_FUNCTION
void M1SphMarshakCoef(const Real rc, const Real rn, const Real rf, Real &ca, Real &cb) {
  const Real al = (rf - rc)/(rc - rn);
  ca = (1.0 + al)*SQR(rc/rf);
  cb = al*SQR(rn/rf);
}

//----------------------------------------------------------------------------------------
//! \fn M1SphCurv
//! \brief spherical polar: the lagged part of (div P)_d at a cell centre (orthonormal
//! r, theta, phi), P = p I + q n n:
//!   d = 0:  - q (1 - n_r^2)/r
//!           + [od] (1/r) d_th P_rt + cot P_rt / r + (1/(r sin)) d_ph P_rp
//!   d = 1:  cot (P_tt - P_pp)/r
//!           + [od] (1/r^3) d_r (r^3 P_rt) + (1/(r sin)) d_ph P_tp
//!   d = 2:  [od] (1/r^3) d_r (r^3 P_rp) + (1/r) d_th P_tp + 2 cot P_tp / r
//! [od] only when od is set.  Centred differences, one-sided at a physical boundary
//! ([il,iu] x [jl,ju] x [kl,ku] readable); E = iw(ec).

// P_ab - p delta_ab = q n_a n_b at cell (k,j,i)
template <class V>
KOKKOS_INLINE_FUNCTION
Real M1SphPab(const V &iw, const int m, const int a, const int b, const int k,
              const int j, const int i, const int ec) {
  return 0.5*(3.0*iw(m,M1_IW_WCHI,k,j,i) - 1.0)*iw(m,M1_IW_N1+a,k,j,i)
         *iw(m,M1_IW_N1+b,k,j,i)*iw(m,ec,k,j,i);
}

// x2t: theta for sin / cot (VtanX2t: the source cell's theta in the ghost rows of a
// theta-periodic seam, so the seam face carries one flux); distances use x2v
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

// x2t = x2v
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
//! \brief vet_col_lat: lagged divergence of the lateral pressure P_ra = D_ra E at a cell
//! of the sp wedge (D_rt, D_rp in tau_ten slots c0+1, c0+2; E = iw(ec)):
//!   d = 0:  (1/r) d_th P_rt + cot P_rt / r + (1/(r sin)) d_ph P_rp
//!   d = 1:  (1/r^3) d_r (r^3 P_rt)
//!   d = 2:  (1/r^3) d_r (r^3 P_rp)
//! differenced as M1SphCurv

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
//! \brief vet_gd: lagged part of the divergence of the tangential anisotropy
//! dP_tt = a E, dP_pp = -a E, P_tp = D_tp E (a = tt(c0+3), D_tp = tt(c0+4), E = iw(ec)):
//!   d = 1:  2 cot a E / r + (1/(r sin)) d_ph P_tp
//!   d = 2:  (1/(r sin^2)) d_th (sin^2 P_tp)
//! The normal-gradient parts sit on the compact face diagonal (ImplicitTransverseTerms):
//! a lagged wide difference of them is unstable in thin regions.  Trig factors use x2t.

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
//! \brief cubed sphere counterpart of M1SphLat d = 0: lateral divergence of v E, v the
//! tangential part of the vet_gd r row as a Cartesian vector (tt slots c0+1..c0+3), as
//! the face sum (1/V) sum_f +- A_f nhat_f . (w_c + w_nb)/2.  corr: at a panel-seam face
//! the mirror-pair mean moves from the chord foot to the face midpoint (M1CsSeamMid).

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
  // value of the seam x2 face jf (pair jf-1, jf), moved to the face midpoint
  auto seam2 = [&](const int jf, const int sd, const Real *nf) {
    auto fv = [&](const int kq) {return 0.5*(wn(kq,jf-1,nf) + wn(kq,jf,nf));};
    const bool up = !(k == ke && b3hi), dn = !(k == ks && b3lo);
    return M1CsSeamMid(fv, g2, m, k, sd, up, dn);
  };
  auto seam3 = [&](const int kf, const int sd, const Real *nf) {
    auto fv = [&](const int jq) {return 0.5*(wn(kf-1,jq,nf) + wn(kf,jq,nf));};
    const bool up = !(j == je && b2hi), dn = !(j == js && b2lo);
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
//! \brief cubed sphere counterpart of M1SphLat d = 1, 2 for the lateral face with the
//! planar unit normal nf: (1/r^3) d_r (r^3 (nf . v) E), differenced as M1SphLat

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
//! \struct M1CsTanGeo
//! \brief cubed sphere vet_gd_tangential: geometry of the trace-free tangential part
//! Q = s1 (e1 e1 - e2 e2) + s2 (e1 e2 + e2 e1) of the vet_gd tensor (tt slots c0+4,
//! c0+5, per-cell tangent frame qf).  Face term, n and t frozen at the face midpoint:
//!   n . div(Q E) = d_n (n.Q.n E) + d_t (t.Q.n E),  d_n = (d_a - cos d_t)/sin compact.

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

// builds an M1CsTanGeo
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

// face normal n and unit tangent t: x2 face (e = 1) t along +x3, x3 face t along +x2
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

// d/ds of the M1CsSeamMid quadratic at the face midpoint; 0 at a foot-value face
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

//! \fn M1CsTanFaceLag
//! \brief n . div(Q E) at face (e, k, j, i), lagged, second order (E = iw(ec)); vmo,
//! vmn: cs_vertex_mono
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
    return M1CsAxisDer(q, h, c, lo, hi, slo, shi, xmono);
  };
  const Real sn = (e == 1) ? tg.sn2(m,k,j) : tg.sn3(m,k,j);
  const Real cs = (e == 1) ? tg.cs2(m,k,j) : tg.cs3(m,k,j);
  const Real al = (e == 1) ? tg.cdxf.x2f(m,k,j,i) : tg.cdxf.x3f(m,k,j,i);
  return ((pq(1,c,n) - pq(0,c,n))/al - cs*0.5*(der(0,n) + der(1,n)))/sn
         + 0.5*(der(0,t) + der(1,t));
}

//! \fn M1CsTanFaceOp
//! \brief the same term as frozen 19-point coefficients vis(kq, jq, c); no foot move and
//! first-order d_t next to a seam (the rhs carries the lagged form)
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
//! frozen-density EOS cache (implicit_eos_cache): per cell, the table's bicubic Hermite
//! patch in (log10 rho, log10 T) collapsed at the frozen rho to 1-D cubics in T on the
//! 2*nt+1 table cells around T^n (implicit_eos_cache_nt); equal to the table up to
//! round-off.  Outside the window, off-table rho, sub-table T, eos_rad_taper: the table.
//! ec(m,n,k,j,i): n = 0 the lowest T node iy0 (-1: no cache); n = 1 + 4c + b: C_b of
//! cell c

constexpr int M1_EC_NTMAX = 8;  // cap on implicit_eos_cache_nt

KOKKOS_INLINE_FUNCTION
int M1EosCacheNComp(const int nt) { return 4*(2*nt + 1) + 1; }

//----------------------------------------------------------------------------------------
//! \fn M1EosCacheEval
//! \brief e(T) and c_v at the frozen density from the cache, in code units, with the
//! arithmetic of EOS_Data::ThermoAt; false (ee, cv untouched) on a miss

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
//! \brief ideal-gas accessor: e and c_v as the ideal branch of EOS_Data::ThermoAt

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
//! \brief the accessor served from the cache, table fallback; nm: miss counter or null

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
//! \brief implicit_tsolve_opac_mode bit 1: the root of y(T) nearest the guess in log T
//! (y need not be monotone where kappa_P T^4 falls with T): outward scan by
//! M1_TSO_NEAR_Q alternately up and down, bisection of the first sign change.  None
//! within 2^(M1_TSO_NEAR_NS/16): ok = false, tout = guess.

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
  Real tu = tguess, td = tguess;  // scan fronts above and below the guess
  Real ta = tguess, tb = tguess;  // bracket: sign(y(ta)) = s0 != sign(y(tb))
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
  // bisection in log T between ta and tb
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
//! \fn M1ImplTemperatureOpacStable
//! \brief implicit_tsolve_opac_mode bit 4 (value 8): the stable root (up-crossing of y)
//! nearest tguess, the end point of the exact exchange ODE at fixed E' (lag-independent
//! where a steep kappa_P(T) gives several roots).  Outward log scan, the nearer bracket
//! bisected in log T; no stable root: ok = false, tout = tguess.

constexpr int  M1_TSO_RAD_NS = 64;                    // max steps per side
constexpr Real M1_TSO_RAD_H0 = 6.7690286301289e-4;    // first step in ln T (ln 2 / 1024)
constexpr Real M1_TSO_RAD_G  = 1.25;                  // step growth per step
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
  // scan fronts above / below the guess; up-crossing bracket: y(lo) < 0 <= y(hi)
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
  // bisection in log T of each bracket found
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

constexpr int M1_TSO_ODE_NS = 32;  // log-T points of the exchange-ODE quadrature

} // namespace radvet
#endif // RAD_VET_RAD_VET_IMPLICIT_HPP_
