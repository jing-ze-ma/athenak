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
#define RTCOL3_NW 35

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
//! \struct RTCol3Hyb
//! \brief THE HYBRID INTERFACE (<problem>/rt_col3_hybrid_tau), per column.
//!
//! WHY.  Deep in a stellar column the two-stream rows are an O(1e3) cancellation: U and D
//! are each ~2 pi b and their difference is the flux, a part in tau of it.  The exact
//! deep limit of the same two-stream is diffusion, F = (4 pi/3) dB/dtau, which carries the
//! same number with NO cancellation and ONE unknown per cell instead of five.  So above
//! the interface (column tau < rt_col3_hybrid_tau) the full two-stream is solved, below
//! it the diffusion equation is, in ONE Newton system -- the deep part is a scalar
//! tridiagonal, the thin part the 5x5 block system, and the two meet at a 1x5 / 5x1 pair.
//!
//! THE INTERFACE CONDITIONS, at the face between cell isp-1 (deep) and isp (thin):
//!   (a) the thin segment's LOWER boundary is the deep limit, exactly the Ucut the whole
//!       column uses at its bottom cut -- U_q = b_f + mu_q dB/dtau -- except that b_f and
//!       dB/dtau are now the IMPLICIT interpolation and difference of b(isp-1) and b(isp)
//!       rather than frozen numbers.  With g = h(isp)/(h(isp-1)+h(isp)),
//!           b_f = (1-g) b(isp) + g b(isp-1) = b(isp) + dB/dtau * h(isp),
//!       so the thin rows are the unchanged cut rows plus a column on b(isp-1), whose
//!       entries are am_q = d U_q / d b(isp-1) = g + mu_q/dtc.
//!   (b) the deep segment's TOP face flux is that same face's two-stream net flux,
//!           F_if = sum_q w_q (U_q - D_q(isp)),
//!       so the energy the deep cell loses upward is exactly what the thin cell gains at
//!       its lower face and the column telescopes across the interface to round-off.
//! The deep segment's bottom face keeps the column's own bottom condition as its flux,
//! F = (4 pi/3) dB/dtau + (int_at_cut ? 2 pi I_int : 0), which under rt_bottom_flux
//! is bit for bit the imposed internal flux.

struct RTCol3Hyb {
  int isp = 0;              // the thin segment's bottom cell; == icut when off
  bool on = false;          // a deep segment exists (isp > icut)
  Real g = 0.0;             // h(isp)/(h(isp-1) + h(isp))
  Real dtc = 1.0;           // h(isp-1) + h(isp)
  Real am[2] = {0.0, 0.0};  // d U_q(interface) / d b(isp-1)
};

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
  DvceArray5D<Real> wk;           // the per-column workspace, see Wk()
  // the PARTITIONED path only (problem/rt_impl_solver = pcr, see
  // two_stream_column_partition.hpp): the per-SEGMENT workspace of the reduced system,
  // (m, nseg*RTCOL3_NRD, k, j), and the number of segments = the Kokkos team size
  DvceArray4D<Real> rd;
  DvceArray4D<Real> dtop;         // (m,q,k,j) the frozen top-face downward intensity
  // THE WARM START (problem/rt_impl_warm).  The converged b of the PREVIOUS call, per
  // cell, and the one before it; a non-positive entry means "no history here" (the
  // Views are zero-initialised, so the first call of a run falls back everywhere).
  DvceArray4D<Real> bprev;        // (m,k,j,i) last converged b
  DvceArray4D<Real> bprev2;       // (m,k,j,i) the call before that (warm = 2)
  DvceArray1D<Real> stat;         // 12 reduction slots, see the launcher
  DualArray1D<RegionSize> size;
  DvceArray4D<Real> dx1;
  EOS_Data eos;
  // ---- switches and scalars --------------------------------------------------------
  Real bdt = 0.0;
  Real sigma = 5.6704e-5;         // the two-stream's own sigma_SB, bit for bit
  Real Iint = 0.0;                // the internal-flux intensity at the cut
  Real bot_flux = 0.0;            // problem/rt_bottom_flux, see rt_bot_flux
  Real mu[2] = {0.0, 0.0};
  Real wf[2] = {0.0, 0.0};
  Real tol = 1.0e-6;
  Real norm_eps = 1.0e-3;         // see problem/rt_impl_norm
  Real dfloor = 0.0;
  Real rgas = 1.0;                // ideal branch: p = rho Rgas T, see PresTempFromEint
  Real gm1 = 0.6666666666666666;  // ideal branch: gamma - 1
  int nq = 2;
  int nseg = 1;                   // problem/rt_impl_nseg, the partitioned path only
  // problem/rt_impl_redpar: solve the REDUCED system of the partitioned path by parallel
  // cyclic reduction across the team's lanes instead of serially on one of them.
  bool redpar = false;
  int maxit = 6;
  int norm = 1;                   // problem/rt_impl_norm
  int cvfreeze = 0;               // problem/rt_impl_cvfreeze
  int warm = 0;                   // problem/rt_impl_warm
  Real dtr = 1.0;                 // bdt/bdt_prev, the warm = 2 extrapolation ratio
  bool exjac = true;              // problem/rt_impl_exjac
  bool dstop = true;              // problem/rt_impl_dstop
  // problem/rt_impl_rescheck.  THE CONFIRMING PASS, WITHOUT THE FACTORISATION.  The
  // convergence test of step 4b reads one number out of the assembled system -- the
  // energy-row residual -- but the test sits AFTER the 5x5 inverses and the elimination,
  // so the pass that does nothing but confirm convergence (it_mean is ~2.0 on the B
  // star: pass 1 solves, pass 2 only checks) pays the whole block factorisation to learn
  // that it had nothing to do.  With this switch the residual is formed FIRST, by
  // ResidRel, which is the SAME arithmetic BuildRow uses for rv[4] and rsc and nothing
  // else; if it is already below tol the pass breaks having done only the formal
  // solution.  The test, the value compared and the iterate are identical, so the answer
  // is BITWISE the answer of the standard path -- this is a scheduling change, not an
  // approximation.  Off by default only because every switch here is.
  bool rescheck = false;
  // ---- TIMING INSTRUMENTATION (problem/rt_impl_ablate, problem/rt_impl_fixit) -------
  // To apportion the cost of a Newton pass between its steps, ablate REPEATS a step
  // rather than skipping it: steps 4a and 4b are IDEMPOTENT (4a zeroes Src before it
  // accumulates, and 4b reads only b, which it does not touch), so running one of them
  // twice leaves the answer BITWISE unchanged while adding exactly its own cost.
  // Skipping instead would change the iterate, and with it the pass count, the timestep
  // and -- with a diverging solve -- the branch mix, which is what makes a difference of
  // wall times unreadable.  The bit mask: 1 = one extra step 4a, 2 = one extra step 4b,
  // 4 = one extra residual-only loop (what a rescheck pass costs).  fixit removes the
  // early exits so that every column runs exactly maxit passes, which gives the cost of
  // a WHOLE pass as the slope in maxit; it is the only one of the two that moves the
  // answer, and only by converging further.  Both branch on RUNTIME members, so nothing
  // is dead-code eliminated.
  int ablate = 0;
  bool fixit = false;
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
  bool ex_iter = false;          // see problem/rt_col3_ex_iter
  // problem/rt_col3_skip_sweep: the entry sweep did not run, so Fb holds nothing.  WRITE
  // the converged face flux of this solve into it instead, so that rad_f2s, the
  // radiative momentum source and the flux/surface/history dumps read the flux of the
  // very field mode 3 applies.  Src/Em are left alone (zero); their consumers are
  // refused at startup.
  bool wrflux = false;
  // problem/rt_col3_hybrid_tau: the column optical depth below which the column is
  // solved as DIFFUSION (1 unknown per cell) instead of as the full two-stream.
  // 0 = off, bitwise the whole-column solve.  See RTCol3Hyb.
  Real hyb_tau = 0.0;
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
  //! \brief the per-cell workspace element.  THE PARTITIONED PATH TRANSPOSES IT.  With
  //! one thread per column the fast index must be j, so that neighbouring threads read
  //! neighbouring words; with a TEAM per column the threads of one wavefront differ in i
  //! instead, and the same layout scatters every access over as many cache lines as
  //! there are lanes.  rt_impl_solver = pcr therefore allocates the workspace as
  //! (m, slot, k, j, i) and this accessor hides which one is live.  The arithmetic is
  //! untouched either way, so the thomas path stays bit for bit what it was.
  template <bool TLAY>
  KOKKOS_INLINE_FUNCTION
  Real &Wk(const int m, const int n, const int i, const int k, const int j) const {
    if constexpr (TLAY) {
      return wk(m,n,k,j,i);
    } else {
      return wk(m,n,i,k,j);
    }
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
                  const RTCol3Hyb &hb, Real sl[3], Real su[3], Real sfu[3],
                  Real sfd[3]) const;
  template <bool TLAY>
  KOKKOS_INLINE_FUNCTION
  void SourceVals(const int m, const int k, const int j, const int i, const int ic,
                  const RTCol3Hyb &hb, const Real cutc, Real &sl, Real &su, Real &sfu,
                  Real &sfd) const;
  template <bool TLAY>
  KOKKOS_INLINE_FUNCTION
  Real ResidRel(const int m, const int k, const int j, const int i,
                const Real eoff) const;
  template <bool TLAY>
  KOKKOS_INLINE_FUNCTION
  void BuildRow(const int m, const int k, const int j, const int i, const int ic,
                const RTCol3Hyb &hb, const Real cutc, const int it, Real A3[5][3],
                Real A1[5], Real Bm[5][5], Real C3[5][3], Real rv[5], Real &rsc) const;
  //! the scalar row of one DEEP cell: diffusion, one unknown, plus the 1x5 interface
  //! coupling C5 on (D_0, D_1, b) of the thin segment's bottom cell.  See RTCol3Hyb.
  template <bool TLAY>
  KOKKOS_INLINE_FUNCTION
  void DeepRow(const int m, const int k, const int j, const int i, const int ic,
               const RTCol3Hyb &hb, const Real kflx, const int it, Real &aa, Real &dd,
               Real &cc, Real C5[3], Real &rv, Real &rsc) const;
  //! the interface cell index of this column: the first cell (from the cut up) whose
  //! column optical depth has fallen below hyb_tau, clamped to [ic+2, ie-2].  Returns ic
  //! when the hybrid is off or the column cannot carry both segments.
  KOKKOS_INLINE_FUNCTION
  int Interface(const int m, const int k, const int j, const int ic) const;
  KOKKOS_INLINE_FUNCTION
  void Solve(const int m, const int k, const int j) const;
};

//----------------------------------------------------------------------------------------
//! \fn Real RTCol3WarmStart
//! \brief THE SHARED WARM START, used by BOTH the serial (RTCol3::Solve) and the
//! partitioned (RTCol3TeamSolve) path -- the initial NEWTON ITERATE only.
//!
//! WHY.  Nothing else about the solve changes: the entry state e*, T*, the frozen layer
//! coefficients, the boundary data and the handover are all still built from the current
//! state, and the tolerance is unchanged.  Newton therefore converges to the SAME root;
//! starting it at the previous call's answer instead of at the entry Planck function
//! only removes passes.  Results are equal to round-off, not bitwise, which is why the
//! default (rt_impl_warm = 0) returns b0 unchanged.
//!
//! THE GUARD.  A warm value that is non-positive, non-finite, or further from the entry
//! Planck function than one clamped Newton step could take it (db in [-0.75 b, 3 b], so
//! the band [0.25 b0, 4 b0]) is discarded for that cell and b0 is used instead.  The
//! discards are counted in stat(21).

KOKKOS_INLINE_FUNCTION
Real RTCol3WarmStart(const RTCol3 &c, const int m, const int k, const int j, const int i,
                     const Real b0) {
  if (c.warm <= 0) return b0;
  const Real bp = c.bprev(m,k,j,i);
  if (!(bp > 0.0) || !isfinite(bp)) return b0;
  Real bw = bp;
  if (c.warm >= 2) {
    const Real bp2 = c.bprev2(m,k,j,i);
    if (bp2 > 0.0 && isfinite(bp2)) bw = bp + (bp - bp2)*c.dtr;
  }
  if (!(bw > 0.0) || !isfinite(bw) || !(b0 > 0.0) ||
      bw < 0.25*b0 || bw > 4.0*b0) {
    Kokkos::atomic_add(&c.stat(21), 1.0);
    return b0;
  }
  return bw;
}

//----------------------------------------------------------------------------------------
//! \fn void RTCol3WarmStore
//! \brief the other half of the same helper: push this call's converged b into the
//! history, shifting the previous level down when the extrapolation needs it.

KOKKOS_INLINE_FUNCTION
void RTCol3WarmStore(const RTCol3 &c, const int m, const int k, const int j, const int i,
                     const Real b) {
  if (c.warm <= 0) return;
  if (c.warm >= 2) c.bprev2(m,k,j,i) = c.bprev(m,k,j,i);
  c.bprev(m,k,j,i) = (b > 0.0 && isfinite(b)) ? b : 0.0;
}

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
                        const RTCol3Hyb &hb, Real sl[3], Real su[3], Real sfu[3],
                        Real sfd[3]) const {
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
  } else if (hb.on) {
    // THE HYBRID INTERFACE.  This cell is the bottom of the THIN segment and its lower
    // face is the interface; the source there is the deep segment's own Planck function
    // interpolated to the face, b_f = (1-g) b(i) + g b(i-1), which is the cut's
    // b(icut) + dB/dtau h_cut with the gradient carried IMPLICITLY by the deep unknown
    // b(i-1).  The c = 0 slot, dead at a real cut, therefore carries d/db(i-1) and
    // BuildRow routes those contributions into the 5x1 interface column A1.
    su[1] = 1.0;
    sfd[0] = hb.g;
    sfd[1] = 1.0 - hb.g;
  } else {
    su[1] = 1.0;
    sfd[1] = 1.0;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RTCol3::SourceVals
//! \brief the same four sources, evaluated at the CURRENT b.

template <bool TLAY>
KOKKOS_INLINE_FUNCTION
void RTCol3::SourceVals(const int m, const int k, const int j, const int i, const int ic,
                        const RTCol3Hyb &hb, const Real cutc, Real &sl, Real &su,
                        Real &sfu, Real &sfd) const {
  Real cl[3], cu[3], cfu[3], cfd[3];
  SourceCoef(m, k, j, i, ic, hb, cl, cu, cfu, cfd);
  const int BBs = 26;
  const Real bm = (i > ic || hb.on) ? Wk<TLAY>(m,BBs,i-1,k,j) : 0.0;
  const Real b0 = Wk<TLAY>(m,BBs,i,k,j);
  const Real bp = (i < ie) ? Wk<TLAY>(m,BBs,i+1,k,j) : 0.0;
  sl  = cl[0]*bm + cl[1]*b0 + cl[2]*bp;
  su  = cu[0]*bm + cu[1]*b0 + cu[2]*bp;
  sfu = cfu[0]*bm + cfu[1]*b0 + cfu[2]*bp;
  // at a hybrid interface the offset is already in cfd[0]*b(i-1), not a frozen constant
  sfd = cfd[0]*bm + cfd[1]*b0 + cfd[2]*bp + ((i == ic && !hb.on) ? cutc : 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn Real RTCol3::ResidRel
//! \brief THE ENERGY-ROW RESIDUAL ALONE, relative to the norm scale -- exactly the
//! number step 4b extracts from BuildRow as fabs(rv[4])/(rsc + eoff), formed with the
//! same operations in the same order so that the comparison against tol is bitwise the
//! same comparison.  It needs only the formal solution of step 4a (through the cached
//! Src and the handover), the entry energy e* cached in step 1, and ONE equation-of-state
//! call e(rho, T(b)); none of the Jacobian, no SourceCoef, no c_v, no 5x5 inverse.
//!
//! The unusable-state branch of BuildRow (a non-positive T*, b or e*) writes rv[4] = 0
//! and leaves the cell out of the solve, so it contributes nothing to the maximum here
//! either.

template <bool TLAY>
KOKKOS_INLINE_FUNCTION
Real RTCol3::ResidRel(const int m, const int k, const int j, const int i,
                      const Real eoff) const {
  const int BBs = 26, EXs = 31, SAs = 32, ESs = 34;
  const Real tk = Tg(m,k,j,i);
  const Real es = Wk<TLAY>(m,ESs,i,k,j);
  const Real b = Wk<TLAY>(m,BBs,i,k,j);
  if (!(tk > 0.0) || !(b > 0.0) || !(es > 0.0)) return 0.0;
  const Real wlo = taublend ? wblend(m,k,j,i) : 0.0;
  const Real whi = taublend ? wblend(m,k,j,i+1) : 0.0;
  const Real wb = 1.0 - 0.5*(wlo + whi);
  const Real tnew = sqrt(sqrt(b*M_PI/sigma));
  const Real enew = EFromT(Rho(m,k,j,i), tnew);
  const Real rv = -(enew - es - bdt*(wb*Wk<TLAY>(m,SAs,i,k,j) + Wk<TLAY>(m,EXs,i,k,j)));
  const Real den = es + eoff;
  return (den > 0.0) ? fabs(rv)/den : 0.0;
}

//----------------------------------------------------------------------------------------
//! \fn void RTCol3::BuildRow
//! \brief the 5x5 blocks of one cell.  A3's three columns are the (U_0, U_1, b) slots of
//! cell i-1, C3's the (D_0, D_1, b) slots of cell i+1 -- the only ones that can be
//! non-zero.  rv is the NEGATIVE residual, so the solve returns the Newton increment.
//! rsc is the scale the energy residual is measured against (the cell's internal energy).

template <bool TLAY>
KOKKOS_INLINE_FUNCTION
void RTCol3::BuildRow(const int m, const int k, const int j, const int i, const int ic,
                      const RTCol3Hyb &hb, const Real cutc, const int it, Real A3[5][3],
                      Real A1[5], Real Bm[5][5], Real C3[5][3], Real rv[5],
                      Real &rsc) const {
  for (int r=0; r<5; ++r) {
    rv[r] = 0.0;
    A1[r] = 0.0;
    for (int c=0; c<3; ++c) {
      A3[r][c] = 0.0;
      C3[r][c] = 0.0;
    }
    for (int c=0; c<5; ++c) Bm[r][c] = 0.0;
  }
  const int EEs = 20, CIs = 22, COs = 24, BBs = 26, EXs = 31, SAs = 32, CVs = 33,
            ESs = 34;
  Real cl[3], cu[3], cfu[3], cfd[3];
  SourceCoef(m, k, j, i, ic, hb, cl, cu, cfu, cfd);
  const Real W = 1.0/Dx(m,k,j,i);
  Real dsdb[3] = {0.0, 0.0, 0.0};
  Real cdu[2] = {0.0, 0.0};
  for (int q=0; q<nq; ++q) {
    const Real E = Wk<TLAY>(m,EEs+q,i,k,j), t = 1.0 - E;
    const Real ci = Wk<TLAY>(m,CIs+q,i,k,j), co = Wk<TLAY>(m,COs+q,i,k,j);
    // ---- the two transport rows -----------------------------------------------------
    Bm[q][q] = 1.0;
    Bm[2+q][2+q] = 1.0;
    if (i < ie) C3[q][q] = -t*t;
    if (i > ic) {
      A3[2+q][q] = -t*t;
    } else if (hb.on) {
      // U enters at the INTERFACE as the deep limit, linear in b(i) and in b(i-1)
      Bm[2+q][4] -= t*t*(1.0 - hb.am[q]);
      A1[2+q] -= t*t*hb.am[q];
    } else {
      Bm[2+q][4] -= t*t;              // U enters at the cut as b(icut) + a constant
    }
    for (int c=0; c<3; ++c) {
      const Real demd = t*ci*cfu[c] + t*co*cl[c] + ci*cu[c] + co*cfd[c];
      const Real demu = t*ci*cfd[c] + t*co*cu[c] + ci*cl[c] + co*cfu[c];
      if (c == 0) {
        if (i > ic) {
          A3[q][2] -= demd;
          A3[2+q][2] -= demu;
        } else if (hb.on) {
          A1[q] -= demd;               // the interface's b(isp-1), a 5x1 column
          A1[2+q] -= demu;
        }
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
  const Real es = Wk<TLAY>(m,ESs,i,k,j);          // e*, cached in step 1: it never changes
  rsc = (es > 0.0) ? es : 1.0;
  const Real b = Wk<TLAY>(m,BBs,i,k,j);
  if (!(tk > 0.0) || !(b > 0.0) || !(es > 0.0)) {
    // an unusable state: the cell takes no part, exactly as it does in the sweep
    Bm[4][4] = 1.0;
    rv[4] = 0.0;
    return;
  }
  const Real wlo = taublend ? wblend(m,k,j,i) : 0.0;
  const Real whi = taublend ? wblend(m,k,j,i+1) : 0.0;
  const Real wb = 1.0 - 0.5*(wlo + whi);
  // THE EXACT ex_iter JACOBIAN (problem/rt_impl_exjac).  With rt_col3_ex_iter the
  // applied source is  A_i = (1-w_i) Src_i + src_ex_i = Src_i + (w F_3)|_lo^hi / dx + Q,
  // the 0.5(w_lo+w_hi) Src of the re-formed handover cancelling the (1-w_i) exactly.  So
  // the Src coefficient is 1, not (1-w_i), and the face-flux term -- LINEAR in the same
  // D and U the block already carries as unknowns -- is differentiated here instead of
  // being left as a Picard lag.
  const bool xj = exjac && ex_iter && taublend && direct;
  const Real fj = -bdt*(xj ? 1.0 : wb);
  const Real tnew = sqrt(sqrt(b*M_PI/sigma));
  const Real enew = EFromT(rho, tnew);
  // de/db = (de/dT) dT/db with b = sigma T^4/pi, i.e. dT_K/db = pi/(4 sigma T^3)
  Real dedb;
  if (cvfreeze > 0 && it >= cvfreeze) {
    dedb = Wk<TLAY>(m,CVs,i,k,j);
  } else {
    dedb = dEdT(rho, enew, tnew)*M_PI/(4.0*sigma*tnew*tnew*tnew);
    Wk<TLAY>(m,CVs,i,k,j) = dedb;
  }
  Bm[4][4] = dedb + fj*dsdb[1];
  if (i > ic) {
    A3[4][2] = fj*dsdb[0];
  } else if (hb.on) {
    A1[4] += fj*dsdb[0];
  }
  C3[4][2] = fj*dsdb[2];
  for (int q=0; q<nq; ++q) {
    if (i < ie) C3[4][q] = fj*cdu[q];
    if (i > ic) {
      A3[4][q] = fj*cdu[q];
    } else if (hb.on) {
      Bm[4][4] += fj*cdu[q]*(1.0 - hb.am[q]);
      A1[4] += fj*cdu[q]*hb.am[q];
    } else {
      Bm[4][4] += fj*cdu[q];
    }
  }
  if (xj) {
    const Real cfl = bdt*wlo*W, cfh = bdt*whi*W;
    for (int q=0; q<nq; ++q) {
      Bm[4][2+q] -= cfh*wf[q];                    // U_q(i),  the cell's upper face
      Bm[4][q] -= cfl*wf[q];                      // D_q(i),  the cell's lower face
      if (i < ie) C3[4][q] += cfh*wf[q];          // D_q(i+1); Dtop is frozen at i = ie
      if (i > ic) {
        A3[4][q] += cfl*wf[q];                    // U_q(i-1)
      } else if (hb.on) {
        Bm[4][4] += cfl*wf[q]*(1.0 - hb.am[q]);   // U(interface), implicit in both b
        A1[4] += cfl*wf[q]*hb.am[q];
      } else {
        Bm[4][4] += cfl*wf[q];                    // U(ic-1) = b(ic) + Ucut
      }
    }
  }
  rv[4] = -(enew - es - bdt*(wb*Wk<TLAY>(m,SAs,i,k,j) + Wk<TLAY>(m,EXs,i,k,j)));
  if (!(Bm[4][4] > 0.0)) Kokkos::atomic_add(&stat(11), 1.0);
}

//----------------------------------------------------------------------------------------
//! \fn int RTCol3::Interface
//! \brief where this column stops being solved as diffusion.  The column optical depth is
//! accumulated from the top face down (the same half/full layer thicknesses the sweep
//! uses); the first cell whose depth reaches rt_col3_hybrid_tau is the topmost DEEP cell,
//! so the interface is the one above it.  Clamped so that the deep segment keeps at least
//! two cells and the two-stream segment at least three.

KOKKOS_INLINE_FUNCTION
int RTCol3::Interface(const int m, const int k, const int j, const int ic) const {
  if (!(hyb_tau > 0.0)) return ic;
  if (ie - 2 < ic + 2) return ic;              // too short to carry both segments
  Real tf = 0.0;
  int isp = ic;
  bool hit = false;
  for (int i=ie; i>=ic; --i) {
    const Real h = Ht(m,k,j,i);
    if (tf + h >= hyb_tau) {
      isp = i + 1;
      hit = true;
      break;
    }
    tf += 2.0*h;
  }
  if (!hit) return ic;                         // the whole column is thin
  if (isp < ic + 2) isp = ic + 2;
  if (isp > ie - 2) isp = ie - 2;
  return isp;
}

//----------------------------------------------------------------------------------------
//! \fn void RTCol3::DeepRow
//! \brief the scalar row of one deep cell.  The unknown is b alone and the transport is
//! the exact two-stream deep limit, F = kflx dB/dtau with kflx = 2 sum_q w_q mu_q (= 4
//! pi/3 for the two-point Gauss-Legendre pair), evaluated between cell CENTRES over the
//! Rosseland depth between them -- the same construction the cut's dB/dtau uses.  The
//! residual is the SAME energy row the thin cells carry, with Src_i = (F_lo - F_hi)/dx.
//!
//! The row is an M-matrix by inspection: the two off-diagonals are -bdt c W kflx/dtau < 0
//! and the diagonal is de/db plus their magnitudes.  At the top deep cell the upper face
//! is the INTERFACE flux sum_q w_q (U_q - D_q(isp)), which is linear in this cell's own b
//! (through U), in b(isp) and in the two D(isp) -- the 1x5 coupling C5, stored on the
//! three columns (D_0, D_1, b) that can be non-zero.

template <bool TLAY>
KOKKOS_INLINE_FUNCTION
void RTCol3::DeepRow(const int m, const int k, const int j, const int i, const int ic,
                     const RTCol3Hyb &hb, const Real kflx, const int it, Real &aa,
                     Real &dd, Real &cc, Real C5[3], Real &rv, Real &rsc) const {
  const int BBs = 26, EXs = 31, SAs = 32, CVs = 33, ESs = 34;
  aa = 0.0;
  cc = 0.0;
  rv = 0.0;
  for (int c=0; c<3; ++c) C5[c] = 0.0;
  const Real tk = Tg(m,k,j,i);
  const Real rho = Rho(m,k,j,i);
  const Real es = Wk<TLAY>(m,ESs,i,k,j);
  rsc = (es > 0.0) ? es : 1.0;
  const Real b = Wk<TLAY>(m,BBs,i,k,j);
  if (!(tk > 0.0) || !(b > 0.0) || !(es > 0.0)) {
    dd = 1.0;                                  // an unusable state takes no part
    return;
  }
  const Real W = 1.0/Dx(m,k,j,i);
  const Real wlo = taublend ? wblend(m,k,j,i) : 0.0;
  const Real whi = taublend ? wblend(m,k,j,i+1) : 0.0;
  const Real wb = 1.0 - 0.5*(wlo + whi);
  // the same exact-handover Jacobian the thin rows use: with rt_impl_exjac the applied
  // source is ((1-w_lo) F_lo - (1-w_hi) F_hi)/dx, otherwise w_b Src plus a lagged src_ex
  const bool xj = exjac && ex_iter && taublend && direct;
  const Real clo = xj ? (1.0 - wlo) : wb;
  const Real chi = xj ? (1.0 - whi) : wb;
  const Real tnew = sqrt(sqrt(b*M_PI/sigma));
  const Real enew = EFromT(rho, tnew);
  Real dedb;
  if (cvfreeze > 0 && it >= cvfreeze) {
    dedb = Wk<TLAY>(m,CVs,i,k,j);
  } else {
    dedb = dEdT(rho, enew, tnew)*M_PI/(4.0*sigma*tnew*tnew*tnew);
    Wk<TLAY>(m,CVs,i,k,j) = dedb;
  }
  dd = dedb;
  if (i > ic) {                                // the lower face: pure diffusion
    const Real dtm = Ht(m,k,j,i-1) + Ht(m,k,j,i);
    const Real kk = (dtm > 0.0) ? kflx/dtm : 0.0;
    aa -= bdt*clo*W*kk;
    dd += bdt*clo*W*kk;
  }                                            // at i = ic the flux is imposed, d/db = 0
  if (i + 1 < hb.isp) {                        // the upper face: pure diffusion
    const Real dtp = Ht(m,k,j,i) + Ht(m,k,j,i+1);
    const Real kk = (dtp > 0.0) ? kflx/dtp : 0.0;
    dd += bdt*chi*W*kk;
    cc -= bdt*chi*W*kk;
  } else {                                     // the upper face IS the interface
    Real sam = 0.0, sap = 0.0;
    for (int q=0; q<nq; ++q) {
      sam += wf[q]*hb.am[q];                   // d F_if / d b(isp-1), through U
      sap += wf[q]*(1.0 - hb.am[q]);           // d F_if / d b(isp)
      C5[q] = -bdt*chi*W*wf[q];                // d F_if / d D_q(isp) = -w_q
    }
    dd += bdt*chi*W*sam;
    C5[2] = bdt*chi*W*sap;
  }
  rv = -(enew - es - bdt*(wb*Wk<TLAY>(m,SAs,i,k,j) + Wk<TLAY>(m,EXs,i,k,j)));
  if (!(dd > 0.0)) Kokkos::atomic_add(&stat(11), 1.0);
}

//----------------------------------------------------------------------------------------
//! \fn void RTCol3::Solve
//! \brief one column: assemble, Newton, block Thomas, apply.

KOKKOS_INLINE_FUNCTION
void RTCol3::Solve(const int m, const int k, const int j) const {
  const int ic = icut(m,k,j);
  if (ic > ie) return;                        // nothing radiative in this column

  const Real sopi = sigma/M_PI;
  // slot layout.  FL aliases DD: in a DEEP cell there is no downward intensity unknown,
  // and the slot carries that cell's LOWER-face diffusion flux instead.
  const int G0 = 0, DP = 15, EE = 20, CI = 22, CO = 24, BB = 26, DD = 27, UU = 29,
            EX = 31, SA = 32, ES = 34, FL = 27;

  // ---- 0. THE HYBRID SPLIT (problem/rt_col3_hybrid_tau) ----------------------------
  // kflx = 2 sum_q w_q mu_q is the deep-limit flux constant OF THIS QUADRATURE (4 pi/3
  // for the two-point Gauss-Legendre pair), so F = kflx dB/dtau is the exact deep limit
  // of the very two-stream solved above the interface -- and reproduces rt_bottom_flux
  // bit for bit, since the sweep sets dB/dtau = 3 F_int/(4 pi) there.
  RTCol3Hyb hb;
  hb.isp = ic;
  Real kflx = 0.0, wsum = 0.0;
  for (int q=0; q<nq; ++q) {
    kflx += 2.0*wf[q]*mu[q];
    wsum += wf[q];
  }
  {
    const int isp0 = Interface(m, k, j, ic);
    if (isp0 > ic) {
      const Real hm = Ht(m,k,j,isp0-1), h0 = Ht(m,k,j,isp0);
      if (hm + h0 > 0.0) {
        hb.on = true;
        hb.isp = isp0;
        hb.dtc = hm + h0;
        hb.g = h0/hb.dtc;
        for (int q=0; q<nq; ++q) hb.am[q] = hb.g + mu[q]/hb.dtc;
      }
    }
  }
  const int ib = hb.isp;      // the bottom cell of the TWO-STREAM segment

  // ---- 1. the frozen per-cell layer coefficients, and b^0 --------------------------
  for (int i=ic; i<=ie; ++i) {
    const Real h = Ht(m,k,j,i);
    for (int q=0; q<nq; ++q) {
      const Real x = h/mu[q];
      const Real e0 = -expm1(-x);
      Wk<false>(m,EE+q,i,k,j) = e0;
      Wk<false>(m,CI+q,i,k,j) = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0 - SQR(x)/3.0);
      Wk<false>(m,CO+q,i,k,j) = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0 - SQR(x)/6.0);
    }
    Wk<false>(m,BB,i,k,j) = RTCol3WarmStart(*this, m, k, j, i, Bb(m,0,i,k,j));
    Wk<false>(m,ES,i,k,j) = Ei(m,k,j,i);     // e*, fixed for the whole solve: cache it once
  }
  // THE NORM SCALE.  problem/rt_impl_norm = 1 measures each cell's energy residual
  // against e_i + eps e_max rather than against e_i alone.  In a stellar column e spans
  // five decades, so the old norm is set entirely by the top cell, whose whole thermal
  // content is a millionth of the base's: it demanded absolute precision there and
  // reported the thick interior, which carries all the energy, as converged long before
  // it was.
  Real emax = 0.0;
  if (norm == 1) {
    for (int i=ic; i<=ie; ++i) {
      const Real e = Wk<false>(m,ES,i,k,j);
      if (e > emax) emax = e;
    }
  }
  const Real eoff = norm_eps*emax;

  // ---- 2. the frozen boundary data --------------------------------------------------
  // the deep-limit gradient at the cut, exactly as the sweep forms it
  Real dbdtau = 0.0;
  if (!cut_legacy && ic + 1 <= ie) {
    const Real dtc = Ht(m,k,j,ic) + Ht(m,k,j,ic+1);
    if (dtc > 0.0) dbdtau = (Bb(m,0,ic,k,j) - Bb(m,0,ic+1,k,j))/dtc;
  }
  // problem/rt_bottom_flux: the cut IS the bottom wall and carries the imposed internal
  // flux, so its gradient is set by that flux, not by the column's own two deepest cells.
  if (bot_flux > 0.0) dbdtau = 3.0*bot_flux/(4.0*M_PI);
  const Real hcut = Ht(m,k,j,ic);
  const Real cutc = dbdtau*hcut;              // b_cutf - b(ic), a frozen offset
  // the deep segment's BOTTOM face flux: the same number the two-stream's cut face
  // carries in the deep limit, sum_q w_q ((b + cutc + I_int + mu dB/dtau) - (b + cutc -
  // mu dB/dtau)).  Under rt_bottom_flux this is exactly the imposed internal flux.
  const Real fbot = kflx*dbdtau + (int_at_cut ? wsum*Iint : 0.0);
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
    Wk<false>(m,EX,i,k,j) = ex + Qb(m,0,i,k,j);
  }

  // ---- 4. Newton ---------------------------------------------------------------------
  int nit = 0;
  int nclamp = 0;
  // THE INTERFACE STATE, re-formed from the current b at every pass: the deep limit the
  // thin segment's bottom cell sees in place of the cut's frozen Ucut.  Equal to the cut
  // values when the hybrid is off, which is what keeps that path bitwise.
  Real cutc_i = cutc;
  Real Ucut_i[2] = {Ucut[0], Ucut[1]};
  Real fif = 0.0;                             // the interface net flux
  Real dbmax = 0.0;
  Real rfin = 0.0;   // DIAGNOSTIC: the residual norm at the LAST iterate examined
  for (int it=0; it<maxit; ++it) {
    nit = it + 1;
    // ---- 4a. the formal solution at the current b, and Src ---------------------------
    for (int rep=(ablate & 1); rep>=0; --rep) {
    for (int i=ic; i<=ie; ++i) Wk<false>(m,SA,i,k,j) = 0.0;
    if (hb.on) {
      const Real dbd = (Wk<false>(m,BB,ib-1,k,j) - Wk<false>(m,BB,ib,k,j))/hb.dtc;
      cutc_i = dbd*Ht(m,k,j,ib);
      for (int q=0; q<nq; ++q) Ucut_i[q] = cutc_i + mu[q]*dbd;
    }
    Real Din[2];
    for (int q=0; q<nq; ++q) Din[q] = Dtop[q];
    for (int i=ie; i>=ib; --i) {
      Real sl, su, sfu, sfd;
      SourceVals<false>(m, k, j, i, ib, hb, cutc, sl, su, sfu, sfd);
      const Real W = 1.0/Dx(m,k,j,i);
      Real acc = 0.0;
      for (int q=0; q<nq; ++q) {
        const Real E = Wk<false>(m,EE+q,i,k,j), t = 1.0 - E;
        const Real ci = Wk<false>(m,CI+q,i,k,j), co = Wk<false>(m,CO+q,i,k,j);
        const Real P = ci*sfu + co*sl;
        const Real Q = ci*su + co*sfd;
        acc += wf[q]*W*(E*(1.0 + t)*Din[q] + E*P - (P + Q));
        Din[q] = t*t*Din[q] + t*P + Q;
        Wk<false>(m,DD+q,i,k,j) = Din[q];
      }
      Wk<false>(m,SA,i,k,j) += acc;
    }
    Real Uin[2];
    for (int q=0; q<nq; ++q) Uin[q] = Wk<false>(m,BB,ib,k,j) + Ucut_i[q];
    for (int i=ib; i<=ie; ++i) {
      Real sl, su, sfu, sfd;
      SourceVals<false>(m, k, j, i, ib, hb, cutc, sl, su, sfu, sfd);
      const Real W = 1.0/Dx(m,k,j,i);
      Real acc = 0.0;
      for (int q=0; q<nq; ++q) {
        const Real E = Wk<false>(m,EE+q,i,k,j), t = 1.0 - E;
        const Real ci = Wk<false>(m,CI+q,i,k,j), co = Wk<false>(m,CO+q,i,k,j);
        const Real Pu = ci*sfd + co*su;
        const Real Qu = ci*sl + co*sfu;
        acc += wf[q]*W*(E*(1.0 + t)*Uin[q] + E*Pu - (Pu + Qu));
        Uin[q] = t*t*Uin[q] + t*Pu + Qu;
        Wk<false>(m,UU+q,i,k,j) = Uin[q];
      }
      Wk<false>(m,SA,i,k,j) += acc;
    }
    // ---- 4a''. THE DEEP SEGMENT: the diffusion fluxes and their divergence ---------
    // The interface flux is the two-stream's OWN net flux at that face, built from the
    // deep-limit U (implicit in b(isp-1), b(isp)) and the solved D(isp): the deep cell
    // loses upward exactly what the thin cell gains at its lower face, so the column
    // still telescopes across the interface.
    if (hb.on) {
      fif = 0.0;
      for (int q=0; q<nq; ++q) {
        fif += wf[q]*((Wk<false>(m,BB,ib,k,j) + Ucut_i[q]) - Wk<false>(m,DD+q,ib,k,j));
      }
      for (int i=ic; i<ib; ++i) {
        Real flo = fbot;
        if (i > ic) {
          const Real dtm = Ht(m,k,j,i-1) + Ht(m,k,j,i);
          flo = (dtm > 0.0)
              ? kflx*(Wk<false>(m,BB,i-1,k,j) - Wk<false>(m,BB,i,k,j))/dtm : 0.0;
        }
        Wk<false>(m,FL,i,k,j) = flo;
      }
      for (int i=ic; i<ib; ++i) {
        const Real fhi = (i + 1 < ib) ? Wk<false>(m,FL,i+1,k,j) : fif;
        Wk<false>(m,SA,i,k,j) = (Wk<false>(m,FL,i,k,j) - fhi)/Dx(m,k,j,i);
      }
    }
    // ---- 4a'. the handover, re-formed from THIS column's OWN flux ------------------
    // problem/rt_col3_ex_iter.  Frozen (the default), src_ex carries the ENTRY sweep's
    // face flux while the (1-w) term carries the converged source, and the two are no
    // longer the divergence of one field: what the cells see is
    //     div[(1-w) F_sweep] + (1-w)(Src_3 - Src_sweep),
    // whose second piece is a real deposit with no flux behind it -- the 2.6 % of F the
    // relaxed B-star column loses between the cut and the top.  With the switch on the
    // whole source is div[(1-w) F_3] of the SAME intensities the energy row absorbs,
    // and telescopes to the faces again.  The Jacobian keeps only the (1-w) part, so
    // the extra coupling enters as a Picard lag, not as a Newton term.  NOT for
    // rad_blend_use_2s > 0, where the conduction operator takes up w F_sweep by
    // construction and the two-stream must give up exactly that.
    if (ex_iter && taublend && direct) {
      if (hb.on) {
        for (int i=ic; i<ib; ++i) {
          const Real fhi = (i + 1 < ib) ? Wk<false>(m,FL,i+1,k,j) : fif;
          const Real wlo = wblend(m,k,j,i), whi = wblend(m,k,j,i+1);
          Wk<false>(m,EX,i,k,j) = 0.5*(wlo + whi)*Wk<false>(m,SA,i,k,j)
                     + (whi*fhi - wlo*Wk<false>(m,FL,i,k,j))/Dx(m,k,j,i) + Qb(m,0,i,k,j);
        }
      }
      for (int i=ib; i<=ie; ++i) {
        Real f3lo = 0.0, f3hi = 0.0;
        for (int q=0; q<nq; ++q) {
          const Real ulo = (i == ib) ? (Wk<false>(m,BB,ib,k,j) + Ucut_i[q]) : Wk<false>(m,UU+q,i-1,k,j);
          const Real dhi = (i == ie) ? Dtop[q] : Wk<false>(m,DD+q,i+1,k,j);
          f3lo += wf[q]*(ulo - Wk<false>(m,DD+q,i,k,j));
          f3hi += wf[q]*(Wk<false>(m,UU+q,i,k,j) - dhi);
        }
        const Real wlo = wblend(m,k,j,i), whi = wblend(m,k,j,i+1);
        Wk<false>(m,EX,i,k,j) = 0.5*(wlo + whi)*Wk<false>(m,SA,i,k,j)
                       + (whi*f3hi - wlo*f3lo)/Dx(m,k,j,i) + Qb(m,0,i,k,j);
      }
    }
    }

    // ---- 4a''. THE RESIDUAL-ONLY CONVERGENCE TEST (problem/rt_impl_rescheck) ---------
    // The same test step 4b makes, made before the factorisation instead of after it, so
    // that a pass which only confirms convergence never assembles or inverts a block.
    // See the switch's note on the struct: this is bitwise, not an approximation.
    if (rescheck || (ablate & 4)) {
      Real rpre = 0.0;
      for (int rep=((ablate & 4) ? 1 : 0); rep>=0; --rep) {
        rpre = 0.0;
        for (int i=ic; i<=ie; ++i) {
          const Real rr = ResidRel<false>(m, k, j, i, eoff);
          if (rr > rpre) rpre = rr;
        }
      }
      if (rescheck) {
        rfin = rpre;
        if (rpre < tol && !fixit) break;
      }
    }

    // ---- 4b. forward elimination ----------------------------------------------------
    Real Gp[5][3], dp[5];
    for (int r=0; r<5; ++r) {
      dp[r] = 0.0;
      for (int c=0; c<3; ++c) Gp[r][c] = 0.0;
    }
    Real rmax = 0.0;
    bool ok = true;
    for (int rep=(ablate & 2); rep>=0; --rep) {
    rmax = 0.0;
    ok = true;
    // ---- 4b-0. the DEEP segment: a SCALAR tridiagonal, one unknown per cell --------
    // The same forward elimination, with 1x1 blocks: no inverse, one divide per cell.
    // Its top row carries the 1x3 interface coupling instead of a scalar c.
    Real dg[3] = {0.0, 0.0, 0.0};
    Real ddp = 0.0;
    if (hb.on) {
      Real gprev = 0.0, dprev = 0.0;
      for (int i=ic; i<ib; ++i) {
        Real aa, dd, cc, C5[3], rvd, rscd;
        DeepRow<false>(m, k, j, i, ic, hb, kflx, it, aa, dd, cc, C5, rvd, rscd);
        const Real dend = rscd + eoff;
        const Real rrd = (dend > 0.0) ? fabs(rvd)/dend : 0.0;
        if (rrd > rmax) rmax = rrd;
        const Real piv = dd - aa*gprev;
        if (!(fabs(piv) > 0.0)) {
          ok = false;
          break;
        }
        dprev = (rvd - aa*dprev)/piv;
        Wk<false>(m,DP+0,i,k,j) = dprev;
        if (i + 1 < ib) {
          gprev = cc/piv;
          Wk<false>(m,G0+0,i,k,j) = gprev;
        } else {
          for (int c=0; c<3; ++c) {
            dg[c] = C5[c]/piv;
            Wk<false>(m,G0+c,i,k,j) = dg[c];
          }
        }
      }
      ddp = dprev;
    }
    for (int i=ib; i<=ie && ok; ++i) {
      Real A3[5][3], A1[5], Bm[5][5], C3[5][3], rv[5];
      Real rsc = 1.0;
      BuildRow<false>(m, k, j, i, ib, hb, cutc, it, A3, A1, Bm, C3, rv, rsc);
      const Real den = rsc + eoff;
      const Real rr = (den > 0.0) ? fabs(rv[4])/den : 0.0;
      if (rr > rmax) rmax = rr;
      if (i == ib && hb.on) {
        // the 5x1 interface column against the deep segment's eliminated top row
        for (int r=0; r<5; ++r) {
          Bm[r][0] -= A1[r]*dg[0];
          Bm[r][1] -= A1[r]*dg[1];
          Bm[r][4] -= A1[r]*dg[2];
          rv[r] -= A1[r]*ddp;
        }
      }
      if (i > ib) {
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
        Wk<false>(m,DP+r,i,k,j) = s;
        for (int c=0; c<3; ++c) {
          Real g = 0.0;
          for (int t=0; t<5; ++t) g += Bi[r][t]*C3[t][c];
          Gp[r][c] = g;
          Wk<false>(m,G0+3*r+c,i,k,j) = g;
        }
      }
    }
    }
    if (!ok) {
      Kokkos::atomic_add(&stat(6), 1.0);
      return;
    }
    rfin = rmax;
    if (rmax < tol && !fixit) break;

    // ---- 4c. back substitution and the clamped update -------------------------------
    Real dbm = 0.0;
    {
    Real ynext[5];
    for (int r=0; r<5; ++r) ynext[r] = 0.0;
    for (int i=ie; i>=ib; --i) {
      Real y[5];
      for (int r=0; r<5; ++r) {
        Real s = Wk<false>(m,DP+r,i,k,j);
        if (i < ie) {
          s -= Wk<false>(m,G0+3*r+0,i,k,j)*ynext[0] + Wk<false>(m,G0+3*r+1,i,k,j)*ynext[1]
             + Wk<false>(m,G0+3*r+2,i,k,j)*ynext[4];
        }
        y[r] = s;
      }
      const Real b = Wk<false>(m,BB,i,k,j);
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
        Wk<false>(m,BB,i,k,j) = b + db;
      }
      y[4] = db;
      for (int r=0; r<5; ++r) ynext[r] = y[r];
    }
    // the DEEP segment, back-substituted from the interface down
    if (hb.on) {
      Real dnext = 0.0;
      for (int i=ib-1; i>=ic; --i) {
        Real sdb = Wk<false>(m,DP+0,i,k,j);
        if (i + 1 < ib) {
          sdb -= Wk<false>(m,G0+0,i,k,j)*dnext;
        } else {
          sdb -= Wk<false>(m,G0+0,i,k,j)*ynext[0] + Wk<false>(m,G0+1,i,k,j)*ynext[1]
               + Wk<false>(m,G0+2,i,k,j)*ynext[4];
        }
        const Real b = Wk<false>(m,BB,i,k,j);
        Real db = sdb;
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
          Wk<false>(m,BB,i,k,j) = b + db;
        }
        dnext = db;
      }
    }
    if (dbm > dbmax) dbmax = dbm;   // the MAX OVER ITERATIONS, not the last one
    if (dump && m == 0 && k == 0 && j == 0) {
      Kokkos::printf("### rt_col3_it it=%d resid=%.6e dbmax=%.6e nclamp=%d\n",
                     it, rmax, dbm, nclamp);
    }
    }
    // the STEP-SIZE stopping rule (problem/rt_impl_dstop).  A Newton that has just taken
    // a negligible step is converged; making it assemble one more block system only to
    // read the residual back costs a whole pass and changes nothing.
    if (fixit) continue;
    if (dstop && dbm < tol) break;
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
    rhsum += bdt*(wb*Wk<false>(m,SA,i,k,j) + Wk<false>(m,EX,i,k,j))*dxi;
    srsum += Wk<false>(m,SA,i,k,j)*dxi;
    RTCol3WarmStore(*this, m, k, j, i, Wk<false>(m,BB,i,k,j));
  }
  Real fnet = 0.0, ftopn = 0.0;
  for (int q=0; q<nq; ++q) ftopn += wf[q]*(Wk<false>(m,UU+q,ie,k,j) - Dtop[q]);
  if (hb.on) {
    fnet = Wk<false>(m,FL,ic,k,j) - ftopn;      // the imposed deep bottom flux
  } else {
    for (int q=0; q<nq; ++q) {
      fnet += wf[q]*((Wk<false>(m,BB,ic,k,j) + Ucut[q]) - Wk<false>(m,DD+q,ic,k,j));
    }
    fnet -= ftopn;
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
  for (int q=0; q<nq; ++q) ftop3 += wf[q]*(Wk<false>(m,UU+q,ie,k,j) - Dtop[q]);
  Kokkos::atomic_add(&stat(16), ftop3);
  Kokkos::atomic_add(&stat(17), Fb(m,0,ie+1,k,j));

  // ---- THE FLUX-TELESCOPING DUMP (problem/rt_outer_verbose, first cycle only) -------
  // Identity (ii) cell by cell.  Src_i is the quadrature-weighted (I_in - I_out) of the
  // two half layers of cell i, so it MUST equal the divergence of the face flux built
  // from the very same converged D/U with the same weights; any difference is an
  // assembly inconsistency between the energy row and the transport rows.  Next to it
  // sit the sweep's own entry-state face flux Fb and source Src -- which is where a
  // difference lives if the two FIELDS differ rather than the algebra -- and the blend
  // weight w, because the applied heating is (1-w) Src + src_ex, not Src.
  if (dump && m == 0) {
    Kokkos::printf("### rt_col3_fsum ic=%d ie=%d srsum=%.10e fnet=%.10e rhsum=%.10e "
                   "F3top=%.10e Fbtop=%.10e Fbcut=%.10e\n",
                   ic, ie, srsum, fnet, rhsum, ftop3, Fb(m,0,ie+1,k,j),
                   Fb(m,0,ic,k,j));
    for (int i=ic; i<ib; ++i) {
      const Real fhi = (i + 1 < ib) ? Wk<false>(m,FL,i+1,k,j) : fif;
      Kokkos::printf("### rt_col3_deep i=%d dtau=%.4e srcdx=%.10e Flo=%.10e Fhi=%.10e "
                     "b=%.6e db_rel=%.4e\n",
                     i, 2.0*Ht(m,k,j,i), Wk<false>(m,SA,i,k,j)*Dx(m,k,j,i),
                     Wk<false>(m,FL,i,k,j), fhi, Wk<false>(m,BB,i,k,j),
                     (Bb(m,0,i,k,j) > 0.0) ? (Wk<false>(m,BB,i,k,j)/Bb(m,0,i,k,j) - 1.0)
                                           : 0.0);
    }
    for (int i=ib; i<=ie; ++i) {
      Real f3lo = 0.0, f3hi = 0.0;
      for (int q=0; q<nq; ++q) {
        const Real ulo = (i == ib) ? (Wk<false>(m,BB,ib,k,j) + Ucut_i[q]) : Wk<false>(m,UU+q,i-1,k,j);
        const Real dhi = (i == ie) ? Dtop[q] : Wk<false>(m,DD+q,i+1,k,j);
        f3lo += wf[q]*(ulo - Wk<false>(m,DD+q,i,k,j));
        f3hi += wf[q]*(Wk<false>(m,UU+q,i,k,j) - dhi);
      }
      const Real dxi = Dx(m,k,j,i);
      const Real wlo = taublend ? wblend(m,k,j,i) : 0.0;
      const Real whi = taublend ? wblend(m,k,j,i+1) : 0.0;
      const Real wb = 1.0 - 0.5*(wlo + whi);
      const Real srcdx = Wk<false>(m,SA,i,k,j)*dxi;
      const Real divf = f3lo - f3hi;
      Kokkos::printf("### rt_col3_flux i=%d dtau=%.4e w=%.4e srcdx=%.10e divF=%.10e "
                     "dif=%.4e sw_srcdx=%.10e F3lo=%.10e Fblo=%.10e exdx=%.10e "
                     "appdx=%.10e db_rel=%.4e\n",
                     i, 2.0*Ht(m,k,j,i), 1.0 - wb, srcdx, divf, srcdx - divf,
                     Src(m,0,i,k,j)*dxi, f3lo, Fb(m,0,i,k,j), Wk<false>(m,EX,i,k,j)*dxi,
                     (wb*Wk<false>(m,SA,i,k,j) + Wk<false>(m,EX,i,k,j))*dxi,
                     (Bb(m,0,i,k,j) > 0.0) ? (Wk<false>(m,BB,i,k,j)/Bb(m,0,i,k,j) - 1.0) : 0.0);
    }
  }
  for (int i=ic; i<=ie; ++i) {
    const Real tk = Tg(m,k,j,i);
    if (!(tk > 0.0)) continue;
    const Real rho = Rho(m,k,j,i);
    const Real es = Ei(m,k,j,i);
    const Real b = Wk<false>(m,BB,i,k,j);
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
                     bdt*(wb*Src(m,0,i,k,j) + Wk<false>(m,EX,i,k,j)), es,
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
  Kokkos::atomic_max(&stat(18), ubmax);   // AFTER the apply loop, which fills it
  Kokkos::atomic_add(&stat(9), budget);
  Kokkos::atomic_add(&stat(10), bscale);
  Kokkos::atomic_max(&stat(19), rfin);
  Kokkos::atomic_add(&stat(20), rfin);
  // problem/rt_col3_skip_sweep: Fb IS this solve's own flux.  Face i carries
  // sum_q w_q (I_up(i-) - I_down(i)), exactly the sweep's definition, with the cut face
  // taking the thermalised upward intensity and the top face the frozen Dtop.  Written
  // after stat(17) so that slot still reports what Fb held on entry (the PREVIOUS call's
  // converged Ftop here, and 0 on the first one).
  if (wrflux) {
    for (int i=ic; i<ib; ++i) Fb(m,0,i,k,j) = Wk<false>(m,FL,i,k,j);
    for (int i=ib; i<=ie; ++i) {
      Real f3lo = 0.0;
      for (int q=0; q<nq; ++q) {
        const Real ulo = (i == ib) ? (Wk<false>(m,BB,ib,k,j) + Ucut_i[q])
                                   : Wk<false>(m,UU+q,i-1,k,j);
        f3lo += wf[q]*(ulo - Wk<false>(m,DD+q,i,k,j));
      }
      Fb(m,0,i,k,j) = f3lo;
    }
    Fb(m,0,ie+1,k,j) = ftop3;
  }
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
