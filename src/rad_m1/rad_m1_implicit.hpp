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

namespace radm1 {

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

// BiCGStab breakdown thresholds: |rho| and |rhat.v| below these times the scale of the
// right-hand side mean the shadow residual has become orthogonal to the Krylov space.
constexpr Real M1_BCG_EPS = 1.0e-300;

// the LAGGED quantities the transverse halo exchanges once per Picard pass.  Everything
// the x2/x3 face fluxes, the lagged off-diagonal Eddington terms and the x1 assembly read
// at a NEIGHBOURING cell is in this list, so the two blocks that share a face build it
// from bit-identical numbers.
// M1_IW_F1 is in the list ONLY for the transverse realizability limiter, which needs the
// face mean of the lagged x1 reduced flux and must build it from numbers both blocks of
// a shared face agree on; it is not read anywhere else, so carrying it changes no
// arithmetic under implicit_trans_limit = none.
constexpr int M1_NHALO_T = 14;
KOKKOS_INLINE_FUNCTION
int M1HaloCompT(const int n) {
  switch (n) {
    case 0: return M1_IW_EP;
    case 1: return M1_IW_WCHI;
    case 2: return M1_IW_N1;
    case 3: return M1_IW_N2;
    case 4: return M1_IW_N3;
    case 5: return M1_IW_KT;
    case 6: return M1_IW_V1;
    case 7: return M1_IW_V2;
    case 8: return M1_IW_V3;
    case 9: return M1_IW_ADV;
    case 10: return M1_IW_A2;
    case 11: return M1_IW_A3;
    case 12: return M1_IW_G0;
    case 13: return M1_IW_F1;
    default: return M1_IW_EP;
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

template <class V>
KOKKOS_INLINE_FUNCTION
Real M1POff(const V &iw, const int m, const int a, const int b,
            const int k, const int j, const int i, const int ec) {
  return M1EddOff(iw(m,M1_IW_WCHI,k,j,i), iw(m,M1_IW_N1+a,k,j,i),
                  iw(m,M1_IW_N1+b,k,j,i))*iw(m,ec,k,j,i);
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
              const int kl, const int ku, const int ec) {
  Real s = 0.0;
  if (d != 0) {
    int ia = (i+1 <= iu) ? (i+1) : i;
    int ib = (i-1 >= il) ? (i-1) : i;
    if (ia != ib) {
      s += (M1POff(iw,m,d,0,k,j,ia,ec) - M1POff(iw,m,d,0,k,j,ib,ec))/((ia - ib)*dx1);
    }
  }
  if (d != 1) {
    int ja = (j+1 <= ju) ? (j+1) : j;
    int jb = (j-1 >= jl) ? (j-1) : j;
    if (ja != jb) {
      s += (M1POff(iw,m,d,1,k,ja,i,ec) - M1POff(iw,m,d,1,k,jb,i,ec))/((ja - jb)*dx2);
    }
  }
  if (thrd && d != 2) {
    int ka = (k+1 <= ku) ? (k+1) : k;
    int kb = (k-1 >= kl) ? (k-1) : k;
    if (ka != kb) {
      s += (M1POff(iw,m,d,2,ka,j,i,ec) - M1POff(iw,m,d,2,kb,j,i,ec))/((ka - kb)*dx3);
    }
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

template <class EosT>
KOKKOS_INLINE_FUNCTION
int M1ImplTemperature(const EosT &eos, const Real dd, const Real tguess,
                      const Real egn, const Real cdtkp_a, const Real rhs_abs,
                      Real &tout, bool &ok) {
  // y(T) = ee(T) + cdtkp_a*T^4 - (egn + rhs_abs)
  Real target = egn + rhs_abs;
  Real tlo = tguess, thi = tguess;
  Real ylo, yhi;
  int nev = 0;
  {
    Real ee, pp, cr, ct, cv;
    eos.ThermoAt(dd, tguess, ee, pp, cr, ct, cv);
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
      Real ee, pp, cr, ct, cv;
      eos.ThermoAt(dd, tlo, ee, pp, cr, ct, cv);
      Real t2 = tlo*tlo;
      ylo = ee + cdtkp_a*t2*t2 - target;
      ++nev;
    }
    ok = (ylo <= 0.0);
  } else if (yhi < 0.0) {
    for (int it=0; it<80 && yhi < 0.0; ++it) {
      thi *= 2.0;
      Real ee, pp, cr, ct, cv;
      eos.ThermoAt(dd, thi, ee, pp, cr, ct, cv);
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
    Real ee, pp, cr, ct, cv;
    eos.ThermoAt(dd, tp, ee, pp, cr, ct, cv);
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

} // namespace radm1
#endif // RAD_M1_RAD_M1_IMPLICIT_HPP_
