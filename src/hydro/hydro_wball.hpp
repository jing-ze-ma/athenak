#ifndef HYDRO_HYDRO_WBALL_HPP_
#define HYDRO_HYDRO_WBALL_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hydro_wball.hpp
//! \brief wball-1104: the dynamic (Kappeli & Mishra) well-balanced scheme in x2 and x3,
//! <hydro>/wb_dirs = all, with a per-cell, per-direction VALIDITY WEIGHT w in [0, 1].
//!
//! Each cell carries, per direction d = x2, x3, its own local hydrostatic background
//! walked from its own state with the same WBBackgroundStencil / closed forms the x1 cache
//! uses (Hydro::wbq2 / wbq3, 15 channels: 5 points {d-1, d-1/2, d, d+1/2, d+1} x (d, e,
//! p)), and a weight w_d (Hydro::wbw) from the stencil:
//!     r_d = |dp + rho dPhi| / (rho |dPhi| + |dp| + eta p),  d over the cells d-1, d+1
//!     w_d = ramp(r_d; r0 -> r1) x ramp(Mach; m0 -> m1)        (ramp: 1 below, 0 above)
//! A face takes w_f = min(w of its two cells); the reconstruction AND the gravity source
//! are blended with that same w_f:
//!     face state  q_f = w_f q_f(WB deviation PLM) + (1 - w_f) q_f(plain PLM)
//!     source      S_d = (1/dx) sum_{f = d+-1/2} sgn_f [w_f (p_bg,f - p_c)
//!                                                  + (1 - w_f)(-rho_c (Phi_f - Phi_c))]
//! The second source term is the linearisation of the first (dp = -rho dPhi), the
//! central -rho dPhi/dx in the same face-sum form, so w = 1 is the exact well-balanced
//! pair, w = 0 the plain scheme with the analytic force, and any w in between a convex
//! combination of two O(h^2)-consistent discretisations (no w-gradient term).
//! x1 keeps the existing global scheme (wbq0, w = 1); under wb_dirs = all its
//! reconstruction is this same routine (dir = 1), whose only fallback is a negative face:
//! the x1 routine's extremum fallback (sL sR <= 0 -> plain PLM) would break the balance
//! wherever the STATE has an extremum along a grid line with a zero deviation -- generic
//! for a potential with an extremum off the grid axes (point mass, tide).

#include "athena.hpp"
#include "reconstruct/plm.hpp"
#include "utils/wb_background.hpp"
#include "hydro/hydro.hpp"

namespace hydro {

//! the validity ramp: 1 for q <= a, 0 for q >= b, linear between
KOKKOS_INLINE_FUNCTION
Real WbRamp(const Real q, const Real a, const Real b) {
  return (q <= a) ? 1.0 : ((q >= b) ? 0.0 : (b - q)/(b - a));
}

//! the blended x2 (dir = 2) or x3 (dir = 3) reconstruction of one pencil: ql_p1(n,i) =
//! the left state at the face +1/2 of cell (k,j,i), qr(n,i) = the right state at its face
//! -1/2.  wbq = Hydro::wbq2 / wbq3, wbw = Hydro::wbw (index dir-1).
KOKKOS_INLINE_FUNCTION
void WbAllPLM(TeamMember_t const &member, const bool use_wb_rho, const int dir,
              const int m, const int k, const int j, const int il, const int iu,
              const DvceArray5D<Real> &q, const DvceArray5D<Real> &wbq,
              const DvceArray5D<Real> &wbw, ScrArray2D<Real> &ql_p1,
              ScrArray2D<Real> &qr) {
  const int dk = (dir == 3) ? 1 : 0, dj = (dir == 2) ? 1 : 0, di = (dir == 1) ? 1 : 0;
  const int nvar = q.extent_int(1);
  const int n2 = wbw.extent_int(3), n3 = wbw.extent_int(2);
  for (int n=0; n<nvar; ++n) {
    const bool wbv = (n == IEN) || (n == IDN && use_wb_rho);
    par_for_inner(member, il, iu, [&](const int i) {
      const Real qm = q(m,n,k-dk,j-dj,i-di), qc = q(m,n,k,j,i), qp = q(m,n,k+dk,j+dj,i+di);
      Real lp, rp;
      PLM(qm, qc, qp, lp, rp);
      if (!wbv) {
        ql_p1(n,i+di) = lp; qr(n,i) = rp;
        return;
      }
      const int var = (n == IEN) ? WBVar::wb_eint : WBVar::wb_dens;
      Real b_m, b_mh, b_c, b_ph, b_p;
      WBReadCache(wbq, var, m, k, j, i, b_m, b_mh, b_c, b_ph, b_p);
      Real lw, rw;
      PLM(qm - b_m, qc - b_c, qp - b_p, lw, rw);
      lw += b_ph; rw += b_mh;
      if (lw < 0.0 || rw < 0.0) {lw = lp; rw = rp;}
      // face weights: min over the two cells of each face (clamped to the built range)
      if (di) {   // x1: the global scheme (w = 1)
        ql_p1(n,i+1) = lw; qr(n,i) = rw;
        return;
      }
      const int km = (k - dk < 0) ? 0 : k - dk, kp = (k + dk > n3 - 1) ? n3 - 1 : k + dk;
      const int jm = (j - dj < 0) ? 0 : j - dj, jp = (j + dj > n2 - 1) ? n2 - 1 : j + dj;
      const Real wc = wbw(m,dir-1,k,j,i);
      const Real wp = fmin(wc, wbw(m,dir-1,kp,jp,i));
      const Real wm = fmin(wc, wbw(m,dir-1,km,jm,i));
      ql_p1(n,i) = (wp >= 1.0) ? lw : wp*lw + (1.0 - wp)*lp;
      qr(n,i) = (wm >= 1.0) ? rw : wm*rw + (1.0 - wm)*rp;
    });
  }
}

//! the same for the derived thermodynamic variables of a general EOS: the pressure
//! channel against the background pressure, Gamma_1 plain (as WbPiecewiseLinearDerX1)
KOKKOS_INLINE_FUNCTION
void WbAllPLMDer(TeamMember_t const &member, const int dir, const int m, const int k,
                 const int j, const int il, const int iu, const DvceArray5D<Real> &qd,
                 const DvceArray5D<Real> &wbq, const DvceArray5D<Real> &wbw,
                 ScrArray2D<Real> &dl_p1, ScrArray2D<Real> &dr) {
  const int dk = (dir == 3) ? 1 : 0, dj = (dir == 2) ? 1 : 0, di = (dir == 1) ? 1 : 0;
  const int nvar = qd.extent_int(1);
  const int n2 = wbw.extent_int(3), n3 = wbw.extent_int(2);
  for (int n=0; n<nvar; ++n) {
    par_for_inner(member, il, iu, [&](const int i) {
      const Real qm = qd(m,n,k-dk,j-dj,i-di), qc = qd(m,n,k,j,i), qp = qd(m,n,k+dk,j+dj,i+di);
      Real lp, rp;
      PLM(qm, qc, qp, lp, rp);
      if (n != IDPR) {
        dl_p1(n,i+di) = lp; dr(n,i) = rp;
        return;
      }
      Real b_m, b_mh, b_c, b_ph, b_p;
      WBReadCache(wbq, WBVar::wb_pres, m, k, j, i, b_m, b_mh, b_c, b_ph, b_p);
      Real lw, rw;
      PLM(qm - b_m, qc - b_c, qp - b_p, lw, rw);
      lw += b_ph; rw += b_mh;
      if (lw < 0.0 || rw < 0.0) {lw = lp; rw = rp;}
      if (di) {
        dl_p1(n,i+1) = lw; dr(n,i) = rw;
        return;
      }
      const int km = (k - dk < 0) ? 0 : k - dk, kp = (k + dk > n3 - 1) ? n3 - 1 : k + dk;
      const int jm = (j - dj < 0) ? 0 : j - dj, jp = (j + dj > n2 - 1) ? n2 - 1 : j + dj;
      const Real wc = wbw(m,dir-1,k,j,i);
      const Real wp = fmin(wc, wbw(m,dir-1,kp,jp,i));
      const Real wm = fmin(wc, wbw(m,dir-1,km,jm,i));
      dl_p1(n,i) = (wp >= 1.0) ? lw : wp*lw + (1.0 - wp)*lp;
      dr(n,i) = (wm >= 1.0) ? rw : wm*rw + (1.0 - wm)*rp;
    });
  }
}

//! the cs shared-seam face state of an edge cell (Hydro::CsSeamFaceStates under
//! wb_dirs = all): the WB deviation face value (lo: the -1/2 face, else +1/2) blended with
//! the plain one by the face weight w, as WbAllPLM forms it in the sweep
KOKKOS_INLINE_FUNCTION
Real WbSeamBlend(const Real qm, const Real qc, const Real qp, const DvceArray5D<Real> &c,
                 const int var, const int m, const int k, const int j, const int i,
                 const bool lo, const Real w, const Real plain) {
  Real b_m, b_mh, b_c, b_ph, b_p;
  WBReadCache(c, var, m, k, j, i, b_m, b_mh, b_c, b_ph, b_p);
  Real lw, rw;
  PLM(qm - b_m, qc - b_c, qp - b_p, lw, rw);
  lw += b_ph; rw += b_mh;
  Real qw = lo ? rw : lw;
  if (lw < 0.0 || rw < 0.0) qw = plain;
  return (w >= 1.0) ? qw : w*qw + (1.0 - w)*plain;
}

//! uniform grids (Cartesian, cs panel coordinates): as WbAllBlendGrid with the uniform PLM
//! (the plain faces of any reconstruction already in ql_p1 / qr)
KOKKOS_INLINE_FUNCTION
void WbAllBlendUni(TeamMember_t const &member, const bool use_wb_rho, const int dir,
                   const bool der, const int m, const int k, const int j, const int il,
                   const int iu, const DvceArray5D<Real> &q, const DvceArray5D<Real> &wbq,
                   const DvceArray5D<Real> &wbw, ScrArray2D<Real> &ql_p1,
                   ScrArray2D<Real> &qr) {
  const int dk = (dir == 3) ? 1 : 0, dj = (dir == 2) ? 1 : 0;
  const int nvar = q.extent_int(1);
  const int n2 = wbw.extent_int(3), n3 = wbw.extent_int(2);
  member.team_barrier();
  for (int n=0; n<nvar; ++n) {
    const bool wbv = der ? (n == IDPR) : ((n == IEN) || (n == IDN && use_wb_rho));
    if (!wbv) continue;
    const int var = der ? WBVar::wb_pres : ((n == IEN) ? WBVar::wb_eint : WBVar::wb_dens);
    par_for_inner(member, il, iu, [&](const int i) {
      const Real qm = q(m,n,k-dk,j-dj,i), qc = q(m,n,k,j,i), qp = q(m,n,k+dk,j+dj,i);
      Real b_m, b_mh, b_c, b_ph, b_p;
      WBReadCache(wbq, var, m, k, j, i, b_m, b_mh, b_c, b_ph, b_p);
      Real lw, rw;
      PLM(qm - b_m, qc - b_c, qp - b_p, lw, rw);
      lw += b_ph; rw += b_mh;
      const Real lp = ql_p1(n,i), rp = qr(n,i);
      if (lw < 0.0 || rw < 0.0) {lw = lp; rw = rp;}
      const int km = (k - dk < 0) ? 0 : k - dk, kp = (k + dk > n3 - 1) ? n3 - 1 : k + dk;
      const int jm = (j - dj < 0) ? 0 : j - dj, jp = (j + dj > n2 - 1) ? n2 - 1 : j + dj;
      const Real wc = wbw(m,dir-1,k,j,i);
      const Real wp = fmin(wc, wbw(m,dir-1,kp,jp,i));
      const Real wm = fmin(wc, wbw(m,dir-1,km,jm,i));
      ql_p1(n,i) = (wp >= 1.0) ? lw : wp*lw + (1.0 - wp)*lp;
      qr(n,i) = (wm >= 1.0) ? rw : wm*rw + (1.0 - wm)*rp;
    });
  }
  member.team_barrier();
}

//! spherical polar (x2 = theta, x3 = phi): the plain faces are ALREADY in ql_p1 / qr
//! (GridPiecewiseLinearX2/X3, the polar-row quadratic included); this forms the WB
//! deviation on the non-uniform grid (PLM_nonuniform at the reconstruction coordinates
//! xv, xf: x2v_rec / xx2f, or x3v / xx3f) and blends it in with the face weights.  nd > 0:
//! the derived channel IDPR of qd (general EOS) against the background pressure.
KOKKOS_INLINE_FUNCTION
void WbAllBlendGrid(TeamMember_t const &member, const bool use_wb_rho, const int dir,
                    const bool der, const int m, const int k, const int j, const int il,
                    const int iu, const DvceArray5D<Real> &q,
                    const DvceArray2D<Real> &xv, const DvceArray2D<Real> &xf,
                    const DvceArray5D<Real> &wbq, const DvceArray5D<Real> &wbw,
                    ScrArray2D<Real> &ql_p1, ScrArray2D<Real> &qr) {
  const int dk = (dir == 3) ? 1 : 0, dj = (dir == 2) ? 1 : 0;
  const int c = (dir == 3) ? k : j;
  const Real x_m = xv(m,c-1), x_mh = xf(m,c), x_c = xv(m,c), x_ph = xf(m,c+1),
             x_p = xv(m,c+1);
  const Real dxLh = x_c - x_mh, dxRh = x_ph - x_c, dxL = x_c - x_m, dxR = x_p - x_c;
  const int nvar = q.extent_int(1);
  const int n2 = wbw.extent_int(3), n3 = wbw.extent_int(2);
  member.team_barrier();
  for (int n=0; n<nvar; ++n) {
    const bool wbv = der ? (n == IDPR) : ((n == IEN) || (n == IDN && use_wb_rho));
    if (!wbv) continue;
    const int var = der ? WBVar::wb_pres : ((n == IEN) ? WBVar::wb_eint : WBVar::wb_dens);
    par_for_inner(member, il, iu, [&](const int i) {
      const Real qm = q(m,n,k-dk,j-dj,i), qc = q(m,n,k,j,i), qp = q(m,n,k+dk,j+dj,i);
      Real b_m, b_mh, b_c, b_ph, b_p;
      WBReadCache(wbq, var, m, k, j, i, b_m, b_mh, b_c, b_ph, b_p);
      Real lw, rw;
      Hydro::PLM_nonuniform(qm - b_m, qc - b_c, qp - b_p, dxL, dxR, dxLh, dxRh, lw, rw);
      lw += b_ph; rw += b_mh;
      const Real lp = ql_p1(n,i), rp = qr(n,i);
      if (lw < 0.0 || rw < 0.0) {lw = lp; rw = rp;}
      const int km = (k - dk < 0) ? 0 : k - dk, kp = (k + dk > n3 - 1) ? n3 - 1 : k + dk;
      const int jm = (j - dj < 0) ? 0 : j - dj, jp = (j + dj > n2 - 1) ? n2 - 1 : j + dj;
      const Real wc = wbw(m,dir-1,k,j,i);
      const Real wp = fmin(wc, wbw(m,dir-1,kp,jp,i));
      const Real wm = fmin(wc, wbw(m,dir-1,km,jm,i));
      ql_p1(n,i) = (wp >= 1.0) ? lw : wp*lw + (1.0 - wp)*lp;
      qr(n,i) = (wm >= 1.0) ? rw : wm*rw + (1.0 - wm)*rp;
    });
  }
  member.team_barrier();
}

//! the fluid-independent view of a Hydro or MHD for the wb_dirs = all routines (shared
//! code: Hydro:: and MHD:: members fill it and call the WbAll* functions below)
struct WbAllCtx {
  MeshBlockPack *pmbp = nullptr;
  EOS_Data eos;
  WBOption wbo;
  bool wgt = true, wbgi = false, wbs = false;
  Real r0 = 0.1, r1 = 0.3, m0 = 0.2, m1 = 0.5, eta = 1.0e-4;
  DvceArray5D<Real> w0, u0, wder, wbq0, wbq2, wbq3, wbw;
  DvceArray4D<Real> phicc, wtemp, pwb;
  DvceArray4D<Real> phi1f, phi2f, phi3f;   // the face potentials (phi0.x1f/x2f/x3f)
  class AmModule *pam = nullptr;
  DvceArray5D<Real> *wbam_u = nullptr, *wbam_w = nullptr;
  DvceFaceFld5D<Real> **wbam_f = nullptr;
};

void WbAllBuild(WbAllCtx &c);
void WbAllGravity(WbAllCtx &c, const Real bdt, const bool with_x1);
void WbAllAmBackground(WbAllCtx &c);

}  // namespace hydro
#endif  // HYDRO_HYDRO_WBALL_HPP_
