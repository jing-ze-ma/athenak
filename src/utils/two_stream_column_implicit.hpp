#ifndef UTILS_TWO_STREAM_COLUMN_IMPLICIT_HPP_
#define UTILS_TWO_STREAM_COLUMN_IMPLICIT_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file two_stream_column_implicit.hpp
//! \brief THE EXACT IMPLICIT COLUMN SOLVE of the grey two-stream transfer coupled to the
//! gas energy (<problem>/rt_implicit_column = 3).  Jiang-2021 in form -- the specific
//! intensities are unknowns, so the coupled system is block-tridiagonal and needs no
//! Jacobi/Gauss-Seidel sweep -- but with OUR exact exponential layer coefficients as the
//! flux function rather than a tau-throttled HLLE.
//!
//! WHY.  The operator-split semi-implicit source damps each cell's NON-LOCAL exchange by
//! its own per-cell factor (1 - e^-x)/x.  Neighbours are damped differently, the O(1e3)
//! absorption/emission cancellation of an optically thick interior no longer cancels, and
//! the residual is a dt-LINEAR numerical forcing of the box's acoustic modes (memory note
//! rt-source-dt-forcing; He box v_rms linear in dt over a 6x range).  Every cheaper cure
//! failed: sub-cycling the local relaxation, re-sweeping the absorbed field, and the
//! single-pass linearised tridiagonal of rt_implicit_column = 1 (which is an INCONSISTENT
//! solve -- it linearises the sweep about the old state and never re-sweeps).  The cure
//! has to be a consistent backward-Euler solve of the whole column, which is this.
//!
//! THE UNKNOWNS, per cell i between the band cut and the top active cell:
//!   D_q(i)  the DOWNWARD intensity at the cell's LOWER face, angle q = 0,1
//!   U_q(i)  the UPWARD   intensity at the cell's UPPER face, angle q = 0,1
//!   b(i)    the Planck function sigma T^4/pi at the cell centre
//! This placement is what makes the transport rows exactly nearest-neighbour: the ray
//! enters cell i at one of its own faces and leaves at the other, crossing the two
//! half-layers the centre-to-centre construction already uses (see rt_layer in
//! two_stream_rt.hpp), and the intermediate CENTRE intensity is substituted out
//! analytically.  5 unknowns per cell, 5x5 blocks, block Thomas, one thread per column.
//!
//! THE ROWS.  With E = 1 - e^{-dtau_half/mu}, t = 1 - E and c_in/c_out the exact
//! linear-in-tau coefficients of RTLayerCoef (the SAME function the explicit sweep uses,
//! so at zero source change the solve reproduces the sweep to round-off):
//!   (1) D_q(i) = t^2 D_q(i+1) + t*(c_in S_f^up + c_out S_l) + (c_in S_u + c_out S_f^dn)
//!   (2) U_q(i) = t^2 U_q(i-1) + t*(c_in S_f^dn + c_out S_u) + (c_in S_l + c_out S_f^up)
//!   (3) e(rho_i, T(b_i)) - e*_i = a dt [ (1-w_i) Src_i + src_ex_i + Q_i ]
//!       Src_i = sum_q W_q [ E(1+t)(D_q(i+1) + U_q(i-1))
//!                           - (t c_in + c_out)(S_f^up + S_f^dn)
//!                           - (t c_out + c_in)(S_l + S_u) ]
//! Every S is a convex combination of b(i-1), b(i), b(i+1) (the BFace emissivity
//! weighting and the dt_l/(dt_l+dt_u) face interpolation), so (1) and (2) are LINEAR in
//! the unknowns and (3)'s only nonlinearity is the gas energy e(rho, T).
//!
//! THE ENERGY ROW IS THE EXACT EOS.  Under the tabulated EOS e is not c_v T: hydrogen and
//! helium ionisation and the LTE radiation term make de/dT swing by more than an order of
//! magnitude across the FeCZ.  So the residual carries e(rho_i, T_i(b_i)) from the same
//! table inversion the apply block uses, and the Jacobian diagonal carries
//! rho c_v(rho,T) dT/db at the CURRENT Newton iterate -- not a frozen heat capacity times
//! (T - T*).  The converged state therefore satisfies the exact energy balance, and where
//! the source vanishes it returns e = e* to the accuracy of the EOS round trip (which
//! this file measures and reports; see RTCOL3_ROUNDTRIP).
//!
//! M-MATRIX.  dSrc/db_i < 0 (the cell's own emission) and dSrc/db_{i+-1} >= 0, dSrc/dD,
//! dSrc/dU >= 0, so row (3) has a positive diagonal and non-positive off-diagonals for
//! ANY dt; rows (1) and (2) have unit diagonals and non-positive off-diagonals because
//! t^2 <= 1 and every source weight is non-negative.  Newton converges in 2-4 iterations
//! and the count of diagonal-sign violations is reported.
//!
//! WHAT IS FROZEN.  The opacity (kappa_R(rho, T*) from the sweep's own cache), the
//! deep-limit gradient dB/dtau at the cut, the top boundary intensity, the tau-blend
//! handover src_ex and the stellar heating Q.  All of these are frozen in the EXPLICIT
//! path too, so mode 3 is not a weaker statement than mode 0 anywhere.
//!
//! NOT IMPLEMENTED: the I/(c dt) diagonal regularisation.  Here the intensities are not
//! evolved variables -- there is no stored I^n -- and the block system is already a
//! strict M-matrix without it, so it would buy nothing but would need the entry-state
//! intensity field stored at every half-layer endpoint.

#include <math.h>
#include <cstdio>

#include "athena.hpp"
#include "eos/eos.hpp"
#include "utils/eint_from_cons.hpp"

namespace two_stream_rt {

// workspace slots per cell, see RTCol3::Solve
#define RTCOL3_NW 33

//----------------------------------------------------------------------------------------
//! \fn bool RTCol3Inv5
//! \brief in-place 5x5 inverse by Gauss-Jordan with partial pivoting.  Returns false on a
//! singular block, which the caller treats as "leave this column alone".

KOKKOS_INLINE_FUNCTION
bool RTCol3Inv5(const Real a[5][5], Real inv[5][5]) {
  Real m[5][10];
  for (int r=0; r<5; ++r) {
    for (int c=0; c<5; ++c) {
      m[r][c] = a[r][c];
      m[r][5+c] = (r == c) ? 1.0 : 0.0;
    }
  }
  for (int c=0; c<5; ++c) {
    int p = c;
    Real best = fabs(m[c][c]);
    for (int r=c+1; r<5; ++r) {
      const Real v = fabs(m[r][c]);
      if (v > best) {
        best = v;
        p = r;
      }
    }
    if (!(best > 0.0)) return false;
    if (p != c) {
      for (int q=0; q<10; ++q) {
        const Real tmp = m[c][q];
        m[c][q] = m[p][q];
        m[p][q] = tmp;
      }
    }
    const Real iv = 1.0/m[c][c];
    for (int q=0; q<10; ++q) m[c][q] *= iv;
    for (int r=0; r<5; ++r) {
      if (r == c) continue;
      const Real f = m[r][c];
      if (f != 0.0) {
        for (int q=0; q<10; ++q) m[r][q] -= f*m[c][q];
      }
    }
  }
  for (int r=0; r<5; ++r) {
    for (int c=0; c<5; ++c) inv[r][c] = m[r][5+c];
  }
  return true;
}

//----------------------------------------------------------------------------------------
//! \struct RTCol3
//! \brief every input the column solve needs, as a POD the launching lambda copies by
//! value.  Nothing here dereferences a host pointer on the device.

struct RTCol3 {
  // ---- state and geometry ----------------------------------------------------------
  DvceArray5D<Real> u0;           // conserved state; IEN is updated in place
  DvceArray5D<Real> bcc;          // cell-centred field (MHD only)
  DvceArray4D<Real> phicc;        // gravitational potential (use_etotgrav only)
  DvceArray3D<Real> cosc;         // cubed-sphere cell cosine
  DvceArray5D<Real> kc;           // (m,0,i,k,j) kappa_R [cm^2/g], frozen
  DvceArray5D<Real> Bb;           // (m,0,i,k,j) entry-state Planck function
  DvceArray5D<Real> Fb;           // (m,0,i,k,j) entry-state net face flux (read only)
  DvceArray5D<Real> Src;          // (m,0,i,k,j) entry-state per-cell source (read only)
  DvceArray5D<Real> Qb;           // (m,0,i,k,j) stellar heating (read only)
  DvceArray4D<Real> Tg;           // (m,k,j,i) T* in Kelvin
  DvceArray4D<Real> wblend;       // (m,k,j,i) tau-blend face weight
  DvceArray3D<int>  icut;         // (m,k,j) the band cut
  DvceArray5D<Real> wk;           // (m, RTCOL3_NW, i, k, j) the per-column workspace
  DvceArray4D<Real> dtop;         // (m,q,k,j) the frozen top-face downward intensity
  DvceArray1D<Real> stat;         // 12 reduction slots, see the launcher
  DualArray1D<RegionSize> size;
  DvceArray4D<Real> dx1;
  EOS_Data eos;
  // ---- switches and scalars --------------------------------------------------------
  Real bdt = 0.0;
  Real sigma = 5.6704e-5;         // the two-stream's own sigma_SB, bit for bit
  Real Iint = 0.0;                // the internal-flux intensity at the cut
  Real mu[2] = {0.0, 0.0};
  Real wf[2] = {0.0, 0.0};
  Real tol = 1.0e-6;
  Real dfloor = 0.0;
  Real rgas = 1.0;                // ideal branch: p = rho Rgas T, see PresTempFromEint
  Real gm1 = 0.6666666666666666;  // ideal branch: gamma - 1
  int nq = 2;
  int maxit = 6;
  int is = 0, ie = 0;
  int is_pp = 0, nx1_pp = 1;
  bool pp = false;
  bool mhd = false;
  bool etg = false;
  bool cs = false;
  bool bface = false;
  bool taublend = false;
  bool int_at_cut = false;
  bool cut_legacy = false;
  bool direct = true;
  bool dump = false;              // one-shot per-cell assembly dump of column (0,ks,js)

  // ---- the two state accessors, the rt_use_cons forms (mode 3 requires it) ----------
  KOKKOS_INLINE_FUNCTION
  Real Rho(const int m, const int k, const int j, const int i) const {
    const Real d = u0(m,IDN,k,j,i);
    return (d > 0.0) ? d : dfloor;
  }
  KOKKOS_INLINE_FUNCTION
  Real Ei(const int m, const int k, const int j, const int i) const {
    const Real e = EintFromCons(u0, m, k, j, i, cs ? cosc(m,k,j) : 0.0, cs,
                                etg, etg ? phicc(m,k,j,i) : 0.0,
                                mhd ? MagEnergyCC(bcc,m,k,j,i) : 0.0);
    if (e > 0.0) return e;
    Real ef = eos.EnergyFromTemperature(u0(m,IDN,k,j,i), eos.tfloor);
    if (!(ef > 0.0)) ef = 1.0e-300;
    return ef;
  }
  //! \brief e(rho, T) with T in KELVIN, in code energy units.  THE EXACT EOS: under the
  //! table this is the same inversion PresTempFromEint uses, so ionisation and the LTE
  //! radiation term are carried in full; under an ideal gas it is the problem's own
  //! p/(Rgas rho) relation, which is NOT eos.EnergyFromTemperature's code-temperature
  //! argument (an ideal EOS_Data carries temp_cgs = 1 and no Rgas).
  KOKKOS_INLINE_FUNCTION
  Real EFromT(const Real rho, const Real tk) const {
    return eos.IsGeneral() ? eos.EnergyFromTemperature(rho, tk/eos.temp_cgs)
                           : rho*rgas*tk/gm1;
  }
  //! \brief de/dT at that state, per KELVIN: rho c_v(rho,T) from the table, or the ideal
  //! constant.  The Jacobian diagonal needs it at the CURRENT iterate, not frozen.
  KOKKOS_INLINE_FUNCTION
  Real dEdT(const Real rho, const Real e, const Real tk) const {
    return eos.IsGeneral()
        ? rho*eos.SpecificHeatCv(rho, e, tk/eos.temp_cgs)/eos.temp_cgs
        : rho*rgas/gm1;
  }
  KOKKOS_INLINE_FUNCTION
  Real Dx(const int m, const int k, const int j, const int i) const {
    return pp ? size.d_view(m).dx1 : dx1(m,k,j,i);
  }
  //! kappa_R * rho, the emissivity the layers are built from
  KOKKOS_INLINE_FUNCTION
  Real Kr(const int m, const int k, const int j, const int i) const {
    return kc(m,0,i,k,j)*Rho(m,k,j,i);
  }
  //! half-layer optical thickness of cell i
  KOKKOS_INLINE_FUNCTION
  Real Ht(const int m, const int k, const int j, const int i) const {
    return 0.5*Kr(m,k,j,i)*Dx(m,k,j,i);
  }

  //! \brief BFaceW of two_stream_rt.hpp: d(BFace)/d(b_far).  Duplicated here rather than
  //! included so that this header stays free of the solver's configuration namespace;
  //! the two are checked against each other by the bitwise gate (rt_implicit_column = 0).
  KOKKOS_INLINE_FUNCTION
  Real FaceW(const Real k_own, const Real k_far) const {
    if (!bface) return 1.0;
    const Real kt = 0.1*k_own;               // RT_BFACE_R
    if (k_far >= kt) return 1.0;
    if (!(k_own > 0.0)) return 1.0;
    if (!(k_far > 0.0)) return 0.0;
    const Real w = k_far/kt;
    return w + (1.0 - w)*k_far/(k_own + k_far);
  }

  KOKKOS_INLINE_FUNCTION
  void SourceCoef(const int m, const int k, const int j, const int i, const int ic,
                  Real sl[3], Real su[3], Real sfu[3], Real sfd[3]) const;
  KOKKOS_INLINE_FUNCTION
  void SourceVals(const int m, const int k, const int j, const int i, const int ic,
                  const Real cutc, Real &sl, Real &su, Real &sfu, Real &sfd) const;
  KOKKOS_INLINE_FUNCTION
  void BuildRow(const int m, const int k, const int j, const int i, const int ic,
                const Real cutc, Real A3[5][3], Real Bm[5][5], Real C3[5][3],
                Real rv[5], Real &rsc) const;
  KOKKOS_INLINE_FUNCTION
  void Solve(const int m, const int k, const int j) const;
};

//----------------------------------------------------------------------------------------
//! \fn void RTCol3::SourceCoef
//! \brief the four source functions cell i needs, as their coefficients on
//! (b(i-1), b(i), b(i+1)).  They are, in the notation of rt_layer:
//!   sl  = s_l of layer i     -- the ray's endpoint at centre i coming from ABOVE
//!   su  = s_u of layer i-1   -- the ray's endpoint at centre i coming from BELOW
//!   sfu = s_f of layer i     -- the source at the cell's UPPER face
//!   sfd = s_f of layer i-1   -- the source at the cell's LOWER face
//! At the top cell every layer-i quantity collapses to b(ie) (the sweep enters the top
//! half layer with the cell's own Planck function at both ends); at the cut cell every
//! layer-(i-1) quantity collapses to b(icut), the face one carrying in addition the
//! frozen deep-limit offset dbdtau*dt_cut, which SourceVals adds.

KOKKOS_INLINE_FUNCTION
void RTCol3::SourceCoef(const int m, const int k, const int j, const int i, const int ic,
                        Real sl[3], Real su[3], Real sfu[3], Real sfd[3]) const {
  for (int c=0; c<3; ++c) {
    sl[c] = 0.0;
    su[c] = 0.0;
    sfu[c] = 0.0;
    sfd[c] = 0.0;
  }
  const Real k0 = Kr(m,k,j,i), h0 = Ht(m,k,j,i);
  if (i < ie) {
    const Real kp = Kr(m,k,j,i+1), hp = Ht(m,k,j,i+1);
    const Real pl = FaceW(kp, k0);            // ds_l(layer i)/db(i)
    const Real pu = FaceW(k0, kp);            // ds_u(layer i)/db(i+1)
    const Real dtc = h0 + hp;
    const Real f = (dtc > 0.0) ? (h0/dtc) : 0.5;
    sl[1] = pl;
    sl[2] = 1.0 - pl;
    sfu[1] = (1.0 - f)*pl + f*(1.0 - pu);
    sfu[2] = (1.0 - f)*(1.0 - pl) + f*pu;
  } else {
    sl[1] = 1.0;
    sfu[1] = 1.0;
  }
  if (i > ic) {
    const Real km = Kr(m,k,j,i-1), hm = Ht(m,k,j,i-1);
    const Real plm = FaceW(k0, km);           // ds_l(layer i-1)/db(i-1)
    const Real pum = FaceW(km, k0);           // ds_u(layer i-1)/db(i)
    const Real dtc = hm + h0;
    const Real f = (dtc > 0.0) ? (hm/dtc) : 0.5;
    su[0] = 1.0 - pum;
    su[1] = pum;
    sfd[0] = (1.0 - f)*plm + f*(1.0 - pum);
    sfd[1] = (1.0 - f)*(1.0 - plm) + f*pum;
  } else {
    su[1] = 1.0;
    sfd[1] = 1.0;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RTCol3::SourceVals
//! \brief the same four sources, evaluated at the CURRENT b.

KOKKOS_INLINE_FUNCTION
void RTCol3::SourceVals(const int m, const int k, const int j, const int i, const int ic,
                        const Real cutc, Real &sl, Real &su, Real &sfu,
                        Real &sfd) const {
  Real cl[3], cu[3], cfu[3], cfd[3];
  SourceCoef(m, k, j, i, ic, cl, cu, cfu, cfd);
  const int BBs = 26;
  const Real bm = (i > ic) ? wk(m,BBs,i-1,k,j) : 0.0;
  const Real b0 = wk(m,BBs,i,k,j);
  const Real bp = (i < ie) ? wk(m,BBs,i+1,k,j) : 0.0;
  sl  = cl[0]*bm + cl[1]*b0 + cl[2]*bp;
  su  = cu[0]*bm + cu[1]*b0 + cu[2]*bp;
  sfu = cfu[0]*bm + cfu[1]*b0 + cfu[2]*bp;
  sfd = cfd[0]*bm + cfd[1]*b0 + cfd[2]*bp + ((i == ic) ? cutc : 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn void RTCol3::BuildRow
//! \brief the 5x5 blocks of one cell.  A3's three columns are the (U_0, U_1, b) slots of
//! cell i-1, C3's the (D_0, D_1, b) slots of cell i+1 -- the only ones that can be
//! non-zero.  rv is the NEGATIVE residual, so the solve returns the Newton increment.
//! rsc is the scale the energy residual is measured against (the cell's internal energy).

KOKKOS_INLINE_FUNCTION
void RTCol3::BuildRow(const int m, const int k, const int j, const int i, const int ic,
                      const Real cutc, Real A3[5][3], Real Bm[5][5], Real C3[5][3],
                      Real rv[5], Real &rsc) const {
  for (int r=0; r<5; ++r) {
    rv[r] = 0.0;
    for (int c=0; c<3; ++c) {
      A3[r][c] = 0.0;
      C3[r][c] = 0.0;
    }
    for (int c=0; c<5; ++c) Bm[r][c] = 0.0;
  }
  const int EEs = 20, CIs = 22, COs = 24, BBs = 26, EXs = 31, SAs = 32;
  Real cl[3], cu[3], cfu[3], cfd[3];
  SourceCoef(m, k, j, i, ic, cl, cu, cfu, cfd);
  const Real W = 1.0/Dx(m,k,j,i);
  Real dsdb[3] = {0.0, 0.0, 0.0};
  Real cdu[2] = {0.0, 0.0};
  for (int q=0; q<nq; ++q) {
    const Real E = wk(m,EEs+q,i,k,j), t = 1.0 - E;
    const Real ci = wk(m,CIs+q,i,k,j), co = wk(m,COs+q,i,k,j);
    // ---- the two transport rows -----------------------------------------------------
    Bm[q][q] = 1.0;
    Bm[2+q][2+q] = 1.0;
    if (i < ie) C3[q][q] = -t*t;
    if (i > ic) {
      A3[2+q][q] = -t*t;
    } else {
      Bm[2+q][4] -= t*t;              // U enters at the cut as b(icut) + a constant
    }
    for (int c=0; c<3; ++c) {
      const Real demd = t*ci*cfu[c] + t*co*cl[c] + ci*cu[c] + co*cfd[c];
      const Real demu = t*ci*cfd[c] + t*co*cu[c] + ci*cl[c] + co*cfu[c];
      if (c == 0) {
        A3[q][2] -= demd;
        A3[2+q][2] -= demu;
      } else if (c == 1) {
        Bm[q][4] -= demd;
        Bm[2+q][4] -= demu;
      } else {
        C3[q][2] -= demd;
        C3[2+q][2] -= demu;
      }
    }
    // ---- the energy row's radiative derivatives --------------------------------------
    const Real g1 = t*ci + co, g2 = t*co + ci;
    for (int c=0; c<3; ++c) {
      dsdb[c] -= wf[q]*W*(g1*(cfu[c] + cfd[c]) + g2*(cl[c] + cu[c]));
    }
    cdu[q] = wf[q]*W*E*(1.0 + t);
  }
  // with the hemispheric mean (nq = 1) the second angle is not an unknown: give its two
  // rows a unit diagonal so the block stays non-singular
  for (int q=nq; q<2; ++q) {
    Bm[q][q] = 1.0;
    Bm[2+q][2+q] = 1.0;
  }
  // ---- the energy row ---------------------------------------------------------------
  const Real tk = Tg(m,k,j,i);
  const Real rho = Rho(m,k,j,i);
  const Real es = Ei(m,k,j,i);
  rsc = (es > 0.0) ? es : 1.0;
  const Real b = wk(m,BBs,i,k,j);
  if (!(tk > 0.0) || !(b > 0.0) || !(es > 0.0)) {
    // an unusable state: the cell takes no part, exactly as it does in the sweep
    Bm[4][4] = 1.0;
    rv[4] = 0.0;
    return;
  }
  const Real wb = taublend ? (1.0 - 0.5*(wblend(m,k,j,i) + wblend(m,k,j,i+1))) : 1.0;
  const Real fj = -bdt*wb;
  const Real tnew = sqrt(sqrt(b*M_PI/sigma));
  const Real enew = EFromT(rho, tnew);
  // de/db = (de/dT) dT/db with b = sigma T^4/pi, i.e. dT_K/db = pi/(4 sigma T^3)
  const Real dedb = dEdT(rho, enew, tnew)*M_PI/(4.0*sigma*tnew*tnew*tnew);
  Bm[4][4] = dedb + fj*dsdb[1];
  A3[4][2] = fj*dsdb[0];
  C3[4][2] = fj*dsdb[2];
  for (int q=0; q<nq; ++q) {
    if (i < ie) C3[4][q] = fj*cdu[q];
    if (i > ic) {
      A3[4][q] = fj*cdu[q];
    } else {
      Bm[4][4] += fj*cdu[q];
    }
  }
  rv[4] = -(enew - es - bdt*(wb*wk(m,SAs,i,k,j) + wk(m,EXs,i,k,j)));
  if (!(Bm[4][4] > 0.0)) Kokkos::atomic_add(&stat(11), 1.0);
}

//----------------------------------------------------------------------------------------
//! \fn void RTCol3::Solve
//! \brief one column: assemble, Newton, block Thomas, apply.

KOKKOS_INLINE_FUNCTION
void RTCol3::Solve(const int m, const int k, const int j) const {
  const int ic = icut(m,k,j);
  if (ic > ie) return;                        // nothing radiative in this column

  const Real sopi = sigma/M_PI;
  // slot layout
  const int G0 = 0, DP = 15, EE = 20, CI = 22, CO = 24, BB = 26, DD = 27, UU = 29,
            EX = 31, SA = 32;

  // ---- 1. the frozen per-cell layer coefficients, and b^0 --------------------------
  for (int i=ic; i<=ie; ++i) {
    const Real h = Ht(m,k,j,i);
    for (int q=0; q<nq; ++q) {
      const Real x = h/mu[q];
      const Real e0 = -expm1(-x);
      wk(m,EE+q,i,k,j) = e0;
      wk(m,CI+q,i,k,j) = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0 - SQR(x)/3.0);
      wk(m,CO+q,i,k,j) = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0 - SQR(x)/6.0);
    }
    wk(m,BB,i,k,j) = Bb(m,0,i,k,j);
  }

  // ---- 2. the frozen boundary data --------------------------------------------------
  // the deep-limit gradient at the cut, exactly as the sweep forms it
  Real dbdtau = 0.0;
  if (!cut_legacy && ic + 1 <= ie) {
    const Real dtc = Ht(m,k,j,ic) + Ht(m,k,j,ic+1);
    if (dtc > 0.0) dbdtau = (Bb(m,0,ic,k,j) - Bb(m,0,ic+1,k,j))/dtc;
  }
  const Real hcut = Ht(m,k,j,ic);
  const Real cutc = dbdtau*hcut;              // b_cutf - b(ic), a frozen offset
  Real Dtop[2], Ucut[2];
  for (int q=0; q<nq; ++q) {
    Dtop[q] = 0.0;                            // filled by the caller's top model below
    Ucut[q] = cutc + (int_at_cut ? Iint : 0.0) + mu[q]*dbdtau;
  }
  // the downward intensity entering the top face, frozen.  The launcher fills it with
  // exactly the sweep's unresolved-column model (RTTopDtau of the top slot's opacity and
  // pressure), so mode 3 and mode 0 see the same boundary.
  for (int q=0; q<nq; ++q) Dtop[q] = dtop(m,q,k,j);

  // ---- 3. the tau-blend handover, frozen from the entry-state sweep -----------------
  for (int i=ic; i<=ie; ++i) {
    Real ex = 0.0;
    if (taublend) {
      const Real dxi = Dx(m,k,j,i);
      const Real wb = 0.5*(wblend(m,k,j,i) + wblend(m,k,j,i+1));
      const Real ft = Fb(m,0,i+1,k,j), fb = Fb(m,0,i,k,j);
      if (direct) {
        ex = wb*Src(m,0,i,k,j)
           + (wblend(m,k,j,i+1)*ft - wblend(m,k,j,i)*fb)/dxi;
      } else {
        ex = -((1.0 - wblend(m,k,j,i+1))*ft - (1.0 - wblend(m,k,j,i))*fb)/dxi
           + (1.0 - wb)*(ft - fb)/dxi;
      }
    }
    wk(m,EX,i,k,j) = ex + Qb(m,0,i,k,j);
  }

  // ---- 4. Newton ---------------------------------------------------------------------
  int nit = 0;
  int nclamp = 0;
  Real dbmax = 0.0;
  for (int it=0; it<maxit; ++it) {
    nit = it + 1;
    // ---- 4a. the formal solution at the current b, and Src ---------------------------
    for (int i=ic; i<=ie; ++i) wk(m,SA,i,k,j) = 0.0;
    Real Din[2];
    for (int q=0; q<nq; ++q) Din[q] = Dtop[q];
    for (int i=ie; i>=ic; --i) {
      Real sl, su, sfu, sfd;
      SourceVals(m, k, j, i, ic, cutc, sl, su, sfu, sfd);
      const Real W = 1.0/Dx(m,k,j,i);
      Real acc = 0.0;
      for (int q=0; q<nq; ++q) {
        const Real E = wk(m,EE+q,i,k,j), t = 1.0 - E;
        const Real ci = wk(m,CI+q,i,k,j), co = wk(m,CO+q,i,k,j);
        const Real P = ci*sfu + co*sl;
        const Real Q = ci*su + co*sfd;
        acc += wf[q]*W*(E*(1.0 + t)*Din[q] + E*P - (P + Q));
        Din[q] = t*t*Din[q] + t*P + Q;
        wk(m,DD+q,i,k,j) = Din[q];
      }
      wk(m,SA,i,k,j) += acc;
    }
    Real Uin[2];
    for (int q=0; q<nq; ++q) Uin[q] = wk(m,BB,ic,k,j) + Ucut[q];
    for (int i=ic; i<=ie; ++i) {
      Real sl, su, sfu, sfd;
      SourceVals(m, k, j, i, ic, cutc, sl, su, sfu, sfd);
      const Real W = 1.0/Dx(m,k,j,i);
      Real acc = 0.0;
      for (int q=0; q<nq; ++q) {
        const Real E = wk(m,EE+q,i,k,j), t = 1.0 - E;
        const Real ci = wk(m,CI+q,i,k,j), co = wk(m,CO+q,i,k,j);
        const Real Pu = ci*sfd + co*su;
        const Real Qu = ci*sl + co*sfu;
        acc += wf[q]*W*(E*(1.0 + t)*Uin[q] + E*Pu - (Pu + Qu));
        Uin[q] = t*t*Uin[q] + t*Pu + Qu;
        wk(m,UU+q,i,k,j) = Uin[q];
      }
      wk(m,SA,i,k,j) += acc;
    }

    // ---- 4b. forward elimination ----------------------------------------------------
    Real Gp[5][3], dp[5];
    for (int r=0; r<5; ++r) {
      dp[r] = 0.0;
      for (int c=0; c<3; ++c) Gp[r][c] = 0.0;
    }
    Real rmax = 0.0;
    bool ok = true;
    for (int i=ic; i<=ie && ok; ++i) {
      Real A3[5][3], Bm[5][5], C3[5][3], rv[5];
      Real rsc = 1.0;
      BuildRow(m, k, j, i, ic, cutc, A3, Bm, C3, rv, rsc);
      const Real rr = (rsc > 0.0) ? fabs(rv[4])/rsc : 0.0;
      if (rr > rmax) rmax = rr;
      if (i > ic) {
        for (int r=0; r<5; ++r) {
          for (int c=0; c<3; ++c) {
            Real s = 0.0;
            for (int t=0; t<3; ++t) s += A3[r][t]*Gp[t+2][c];
            Bm[r][(c < 2) ? c : 4] -= s;
          }
          Real s2 = 0.0;
          for (int t=0; t<3; ++t) s2 += A3[r][t]*dp[t+2];
          rv[r] -= s2;
        }
      }
      Real Bi[5][5];
      if (!RTCol3Inv5(Bm, Bi)) {
        ok = false;
        break;
      }
      for (int r=0; r<5; ++r) {
        Real s = 0.0;
        for (int c=0; c<5; ++c) s += Bi[r][c]*rv[c];
        dp[r] = s;
        wk(m,DP+r,i,k,j) = s;
        for (int c=0; c<3; ++c) {
          Real g = 0.0;
          for (int t=0; t<5; ++t) g += Bi[r][t]*C3[t][c];
          Gp[r][c] = g;
          wk(m,G0+3*r+c,i,k,j) = g;
        }
      }
    }
    if (!ok) {
      Kokkos::atomic_add(&stat(6), 1.0);
      return;
    }
    if (rmax < tol) break;

    // ---- 4c. back substitution and the clamped update -------------------------------
    Real ynext[5];
    for (int r=0; r<5; ++r) ynext[r] = 0.0;
    Real dbm = 0.0;
    for (int i=ie; i>=ic; --i) {
      Real y[5];
      for (int r=0; r<5; ++r) {
        Real s = wk(m,DP+r,i,k,j);
        if (i < ie) {
          s -= wk(m,G0+3*r+0,i,k,j)*ynext[0] + wk(m,G0+3*r+1,i,k,j)*ynext[1]
             + wk(m,G0+3*r+2,i,k,j)*ynext[4];
        }
        y[r] = s;
      }
      const Real b = wk(m,BB,i,k,j);
      Real db = y[4];
      if (b > 0.0) {
        if (db > 3.0*b) {
          db = 3.0*b;
          ++nclamp;
        } else if (db < -0.75*b) {
          db = -0.75*b;
          ++nclamp;
        }
        const Real rel = fabs(db)/b;
        if (rel > dbm) dbm = rel;
        wk(m,BB,i,k,j) = b + db;
      }
      y[4] = db;
      for (int r=0; r<5; ++r) ynext[r] = y[r];
    }
    if (dbm > dbmax) dbmax = dbm;   // the MAX OVER ITERATIONS, not the last one
    if (dbm < 1.0e-14) break;
  }

  // ---- 5. apply, and the column energy budget --------------------------------------
  // Two identities the converged solve must satisfy to round-off, and which nothing in
  // it enforces by construction:
  //   (i)  sum_i de_i dx_i = a dt sum_i [ (1-w_i) Src_i + src_ex_i ] dx_i  -- the energy
  //        rows, summed; it fails if the EOS round trip or a clamp has bitten;
  //   (ii) sum_i Src_i dx_i = F(cut) - F(top)  -- the sweep telescopes, i.e. the source
  //        really is the divergence of the face flux the same intensities carry.
  Real budget = 0.0, bscale = 0.0, rtmax = 0.0, ubmax = 0.0;
  Real rhsum = 0.0, srsum = 0.0;
  for (int i=ic; i<=ie; ++i) {
    const Real dxi = Dx(m,k,j,i);
    const Real wb = taublend ? (1.0 - 0.5*(wblend(m,k,j,i) + wblend(m,k,j,i+1))) : 1.0;
    rhsum += bdt*(wb*wk(m,SA,i,k,j) + wk(m,EX,i,k,j))*dxi;
    srsum += wk(m,SA,i,k,j)*dxi;
  }
  Real fnet = 0.0;
  for (int q=0; q<nq; ++q) {
    const Real ftop = wf[q]*(wk(m,UU+q,ie,k,j) - Dtop[q]);
    const Real fcut = wf[q]*((wk(m,BB,ic,k,j) + Ucut[q]) - wk(m,DD+q,ic,k,j));
    fnet += fcut - ftop;
  }
  Kokkos::atomic_add(&stat(12), rhsum);
  Kokkos::atomic_add(&stat(13), srsum);
  Kokkos::atomic_add(&stat(14), fnet);
  Kokkos::atomic_add(&stat(15), fabs(srsum) + fabs(fnet));
  // the emergent flux of the CONVERGED field against the one the explicit sweep left in
  // Fb (which is what the flux dump, rad_f2s and rt_rad_force still see, by design: the
  // frozen handover in src_ex is paired with exactly that flux).  On a relaxed column
  // the two must agree, and how well is the statement that mode 3 has not moved the
  // radiation field away from the formal solution the sweep was validated against.
  Real ftop3 = 0.0;
  for (int q=0; q<nq; ++q) ftop3 += wf[q]*(wk(m,UU+q,ie,k,j) - Dtop[q]);
  Kokkos::atomic_add(&stat(16), ftop3);
  Kokkos::atomic_add(&stat(17), Fb(m,0,ie+1,k,j));
  Kokkos::atomic_max(&stat(18), ubmax);
  for (int i=ic; i<=ie; ++i) {
    const Real tk = Tg(m,k,j,i);
    if (!(tk > 0.0)) continue;
    const Real rho = Rho(m,k,j,i);
    const Real es = Ei(m,k,j,i);
    const Real b = wk(m,BB,i,k,j);
    if (!(b > 0.0)) continue;
    const Real tnew = sqrt(sqrt(b/sopi));
    const Real enew = EFromT(rho, tnew);
    if (!(enew > 0.0) || !isfinite(enew)) {
      Kokkos::atomic_add(&stat(7), 1.0);
      continue;
    }
    const Real eref = EFromT(rho, tk);
    if (es > 0.0 && eref > 0.0) {
      const Real rt = fabs(eref - es)/es;
      if (rt > rtmax) rtmax = rt;
    }
    Real de = enew - es;
    // how far the implicit solve moved this cell from the state the explicit sweep saw,
    // and what the explicit source would have deposited instead
    const Real b0 = Bb(m,0,i,k,j);
    if (b0 > 0.0) {
      const Real rb = fabs(b/b0 - 1.0);
      if (rb > ubmax) ubmax = rb;
    }
    if (dump && m == 0) {
      const Real wb = taublend ? (1.0 - 0.5*(wblend(m,k,j,i) + wblend(m,k,j,i+1))) : 1.0;
      Kokkos::printf("### rt_col3_cell i=%d dtau=%.4e T=%.6e rho=%.4e b0=%.6e "
                     "b=%.6e du_u=%.4e de=%.6e de_expl=%.6e e=%.6e cv=%.4e\n",
                     i, 2.0*Ht(m,k,j,i), tk, rho, b0, b,
                     (b0 > 0.0) ? (b/b0 - 1.0) : 0.0, de,
                     bdt*(wb*Src(m,0,i,k,j) + wk(m,EX,i,k,j)), es,
                     dEdT(rho, es, tk));
    }
    if (!((es + de) > 0.0)) {
      de = -0.999*es;
      Kokkos::atomic_add(&stat(8), 1.0);
    }
    u0(m,IEN,k,j,i) += de;
    const Real dxi = Dx(m,k,j,i);
    budget += de*dxi;
    bscale += fabs(de)*dxi;
  }
  Kokkos::atomic_add(&stat(0), static_cast<Real>(nit));
  Kokkos::atomic_add(&stat(1), 1.0);
  Kokkos::atomic_max(&stat(2), static_cast<Real>(nit));
  Kokkos::atomic_add(&stat(3), static_cast<Real>(nclamp));
  Kokkos::atomic_max(&stat(4), dbmax);
  Kokkos::atomic_max(&stat(5), rtmax);
  Kokkos::atomic_add(&stat(9), budget);
  Kokkos::atomic_add(&stat(10), bscale);
}

//----------------------------------------------------------------------------------------
//! \fn void RTCol3Launch
//! \brief one thread per column.  The struct is copied into the lambda by value, so
//! nothing on the host is dereferenced on the device.

inline void RTCol3Launch(const RTCol3 &c, const int nmb1, const int ks, const int ke,
                         const int js, const int je) {
  RTCol3 cc = c;
  par_for("rt_col3", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    cc.Solve(m, k, j);
  });
}

}  // namespace two_stream_rt
#endif  // UTILS_TWO_STREAM_COLUMN_IMPLICIT_HPP_
