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
constexpr char kVetPredRstMagic[8] = {'M', '1', 'P', 'R', 'E', 'D', '0', '1'};

// time_scheme
constexpr int VET_TIME_BE = 0;        // Lie splitting: Heun hydro, then one BE solve
constexpr int VET_TIME_HESDIRK2 = 1;  // H-ESDIRK2: one stage solve inside each Heun stage
// t2_solve: what the next ImplicitSolve does
constexpr int VET_T2S_NONE = 0;     // plain backward Euler
constexpr int VET_T2S_BESTORE = 1;  // backward Euler, K1 = (Y - rhs)/dt stored
constexpr int VET_T2S_STAGE1 = 2;   // stage solve 1 (gamma dt), K2 stored
constexpr int VET_T2S_STAGE2 = 3;   // stage solve 2 (gamma dt), K1 (FSAL) stored
// channels of t2k1 / t2k2 / t2inc: E, gas momentum, total energy, face fluxes F0 x1..x3
constexpr int VET_T2_E = 0;
constexpr int VET_T2_M1 = 1;
constexpr int VET_T2_EN = 4;
constexpr int VET_T2_F1 = 5;
constexpr int VET_T2_NK = 8;
// vet_sc: extrapolated components of vet_cell (VETRAY_CHI .. VETRAY_D11+5)
constexpr int VET_T2_NVET = 10;
// hesdirk2 slope: int32 have, int32 nch, Real pred2_dt, dt_prev, vprev; nch slabs
// (K1, ipred2, vet_prev) behind the predictor slabs
constexpr char kVetTime2RstMagic[8] = {'M', '1', 'T', 'I', 'M', 'E', '2', 'A'};
// implicit_one_pass > 0: 9 Reals (last two contractions, solve counters)
constexpr char kVetOnePassRstMagic[8] = {'M', '1', 'O', 'N', 'E', 'P', '0', '1'};
// implicit_mr_every / vet_sc_every > 1: multi-rate window and tensor cadence, 10 Reals
constexpr char kVetMRRstMagic[8] = {'M', '1', 'M', 'R', 'W', 'I', 'N', '1'};
// time2_auto_be back-off while active: int32 len, int32 left
constexpr char kVetAutoBeRstMagic[8] = {'M', '1', 'A', 'U', 'T', 'O', 'B', 'E'};
// removed closure memory of old files (read past): int32 have, nch = 4; slabs last
constexpr char kVetCtrRstMagic[8] = {'M', '1', 'C', 'T', 'R', 'L', 'X', '1'};
// rad_signal_speed inputs: int32 have = 1, nch, mode, 0; nch slabs (opac T, tensor)
constexpr char kVetRssRstMagic[8] = {'M', '1', 'R', 'S', 'S', 'I', 'G', '1'};
// vet_gd_async pending source (ln chi, ln S): int32 have = 1, nch = 2, ncycle, 0
constexpr char kVetVgdRstMagic[8] = {'M', '1', 'V', 'G', 'D', 'A', 'S', '1'};
// conservative gas source with hesdirk2: int32 have = 1, nch (3 | 0), csk1_ok, 0; nch
// slabs of csk1
constexpr char kVetCskRstMagic[8] = {'M', '1', 'C', 'S', 'K', '1', 'S', '1'};
// time2_face_fct: face fluxes t2pf1 (FaxRst*)
constexpr char kVetFaxRstMagic[8] = {'M', '1', 'F', 'A', 'X', 'R', 'S', '1'};
// vet_rotate_avg history: VET_RA_RST_NH int64 words (RotAvgRstHdr), then nslab slabs
constexpr char kVetRaRstMagic[8] = {'M', '1', 'R', 'O', 'T', 'A', 'V', '1'};
// old files: stage-2 drop slabs (skipped)
constexpr char kVetThcRstMagic[8] = {'M', '1', 'T', 'H', 'C', 'R', 'S', '1'};
// old files: cross-fade state behind the thc slabs (skipped)
constexpr char kVetXfdRstMagic[8] = {'M', '1', 'X', 'F', 'D', 'R', 'S', '1'};
constexpr int kVetXfdRstNi = 24;  // int64 header words
constexpr int kVetXfdRstNr = 8;   // Real header words
// vet_x1_periodic (adaptive): the last build's exit plane vet_x1ent, ray r of row (k, j)
// at slab r / n1, cell i = r % n1 (RadiationVet::X1eRst*), behind the fax slabs
constexpr char kVetX1eRstMagic[8] = {'M', '1', 'X', '1', 'E', 'R', 'S', '1'};

// implicit_solver (both converge to the same fixed point)
constexpr int VET_ISOLV_LINE_JACOBI = 0;  // x1 tridiagonal exact, x2/x3 lagged (M-matrix)
constexpr int VET_ISOLV_BICGSTAB = 1;  // BiCGStab, x1 line solve as right preconditioner

// implicit_offdiag: the off-diagonal Eddington terms sum_{e != d} d_e P_de (lagging them
// diverges in thin cells, where they scale as (c dt/dx)^2)
constexpr int VET_OD_LAGGED   = 0;  // previous Picard iterate, on the right-hand side
constexpr int VET_OD_OPERATOR = 1;  // in the operator (bicgstab; not an M-matrix)
constexpr int VET_OD_NONE     = 2;  // dropped: P diagonal in the grid frame

// implicit_trans_limit: transverse realizability limiter.  lp adds the lagged opacity
// |G_f|/(phi_f E_f), phi_f = fmax sqrt(max(1 - f1^2, 0.01)), so |F_f| <= c phi_f E_f
constexpr int VET_TLIM_NONE = 0;
constexpr int VET_TLIM_LPR  = 2;  // bound max(kt, |G|/(fmax E)), tapered by min(1, kt dx)
constexpr int VET_TLIM_LP   = 1;  // lagged limiter opacity (above)

// boundary types of the implicit solve, as the total normal flux F_f at the face
constexpr int VET_IBC_MARSHAK  = 0;  // free surface, F_f = +- c marshak_q E
constexpr int VET_IBC_FLUX     = 1;  // imposed F_f = implicit_flux_x1min/max
constexpr int VET_IBC_REFLECT  = 2;  // F_f = 0
constexpr int VET_IBC_PERIODIC = 3;  // cyclic tridiagonal (Sherman-Morrison)
constexpr int VET_IBC_EFIX     = 4;  // Dirichlet: boundary cell E frozen at its start
constexpr int VET_IBC_OUTFLOW  = 5;  // x2: F_f = max(+-F2_cell, 0), lagged
constexpr int VET_IBC_X2FLUX   = 6;  // x2: F_f = ifx2_lo/hi from the pgen (SetX2FaceFlux)

// implicit_flux: spatial form of the implicit face E-flux
constexpr int VET_IFLUX_CENTRAL = 0;  // face-eliminated central (diffusion) form
constexpr int VET_IFLUX_BERTHON = 2;  // A_up + alpha F_HLL (Berthon), no F_diff
constexpr int VET_IFLUX_BLEND   = 3;  // (1 - w_f) central + w_f berthon, w_f lagged

// implicit_thin_switch: dt-free face weight of the half-range blend (VetThinFace)
constexpr int VET_THSW_FS   = 3;  // smoothstep of max(|H|/J, C1/1.5), C1 = <|I-S|>/J
constexpr int VET_THSW_FSKN = 5;  // max(fs, kn)
constexpr int VET_THSW_FSI  = 7;  // fs; h+- isotropic where the rays miss E^n (VetThinHr)
constexpr int VET_THSW_FSE  = 8;  // max(|H|/J, eps C1/1.5), eps the thermal fraction
constexpr int VET_THSW_FSQ  = 6;  // max(|H|/J, q_c C1/1.5), q_c = min(E^n/J, J/E^n)^qp

// implicit_recon: states of the HLL part
constexpr int VET_IRECON_DC    = 0;  // piecewise constant
constexpr int VET_IRECON_PLMDC = 1;  // plm as a lagged deferred correction (rhs)

// implicit_enthalpy: face value of the enthalpy flux a_d E, a_d = v_d + (v.D)_d.  The
// matrix carries the donor form; the others add their lagged difference to the rhs.
constexpr int VET_IENTH_UPWIND  = 0;  // donor cell (first order)
constexpr int VET_IENTH_CENTRAL = 1;  // a_f (E_L + E_R)/2, a_f = (a_L + a_R)/2
constexpr int VET_IENTH_PLM     = 2;  // a_f E_f, E_f the upwind van Leer plm value
constexpr int VET_IENTH_PLMC    = 3;  // mean of the two plm face values
constexpr int VET_IENTH_LIN     = 4;  // unlimited upwind linear (Fromm) value

// vet-positivity counters (pos_cnt_d / pos_cnt); energies volume-weighted
constexpr int VET_POS_GAS    = 1;  // cell-solves where the gas eint limiter fired
constexpr int VET_POS_GAS_DE = 2;  // energy moved gas <- radiation
constexpr int VET_POS_FLR    = 3;  // cell-solves where E was raised to e_floor
constexpr int VET_POS_FLR_DE = 4;  // floor energy taken from the gas
constexpr int VET_POS_FLR_UN = 5;  // floor energy not covered by the gas (created)
constexpr int VET_POS_FLRS   = 8;  // residual raise below 1e-6 e_floor
constexpr int VET_POS_FLR_MV = 9;  // floor energy taken from the column's radiation
constexpr int VET_POS_N      = 10;

// implicit work array iw (m,n,k,j,i): per solve, except the start-of-step EN / EGN
constexpr int VET_IW_EN   = 0;   // E^n
constexpr int VET_IW_TP   = 1;   // T of the state the rays read
constexpr int VET_IW_F1   = 2;   // cell-centred flux of that state (3)
constexpr int VET_IW_F2   = 3;
constexpr int VET_IW_F3   = 4;
constexpr int VET_IW_KT   = 5;   // copy of opac(VET_OP_T), for the x1 halo exchange
constexpr int VET_IW_S1   = 6;   // FCT scratch (Time2FctSlot)
constexpr int VET_IW_S2   = 7;
constexpr int VET_NIW = 8;

KOKKOS_INLINE_FUNCTION
void VetCsSeamGeom(const int q, const int n, Real &dth, Real &angm, Real &sfoot,
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

KOKKOS_INLINE_FUNCTION
void VetSphMarshakCoef(const Real rc, const Real rn, const Real rf, Real &ca, Real &cb) {
  const Real al = (rf - rc)/(rc - rn);
  ca = (1.0 + al)*SQR(rc/rf);
  cb = al*SQR(rn/rf);
}

} // namespace radvet
#endif // RAD_VET_RAD_VET_IMPLICIT_HPP_
