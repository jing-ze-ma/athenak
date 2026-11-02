#ifndef HYDRO_RSOLVERS_LHLLC_C1_HYD_HPP_
#define HYDRO_RSOLVERS_LHLLC_C1_HYD_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file lhllc_c1_hyd.hpp
//! \brief lhllc2nd-1102: the LHLLC Riemann solver (lhllc_hyd.hpp) with a flux that is C1
//! in the left/right states.  Selected by <hydro>/lhllc_c1 = eps > 0 (read only when
//! named; without it rsolver = lhllc is the legacy LHLLC, bit for bit).
//!
//! WHY: an RK2 step is second order in dt only if the semi-discrete right-hand side is
//! smooth (C1) along the trajectory.  A Lipschitz kink that the state crosses once costs
//! one O(dt^2) local error and keeps the global error O(dt^2), but a kink the state SITS
//! ON (crosses at every step) costs O(dt^2) per step, i.e. O(dt) globally: first order.
//! A low-Mach flow sits on several switches of LHLLC (and of HLLC):
//!   (1) the contact switch: the flux takes the L* or the R* state by sign(S*); with
//!       u_n ~ 0 the contact speed S* changes sign everywhere, all the time;
//!   (2) the PVRS wave-speed factor q = 1 for p* <= p, sqrt(1 + alpha (p*/p - 1)) above;
//!       at low Mach p* = p (1 + O(M)) on every face;
//!   (3) chi = min(1, max(|u_L|, |u_R|)/c): the abs at u_n = 0 and the max at
//!       |u_L| = |u_R|;
//!   (4) the odd-even sensor s (utils/cb_sensor.hpp): the clamps at J d2 = 0 and at
//!       |J| = |d2| and the sqrt(1 - s) of its chi floor (infinite slope at s = 1);
//!   (5) the x1 phi floor max(phi, x1_phi_min).
//! HLLE has none of (1)-(4) at low Mach (its min/max wave speeds switch only at
//! |u| ~ c), which is why it keeps second order where HLLC and LHLLC do not.
//!
//! WHAT: each switch becomes a C1 blend over a band of width eps (bit k-1 of
//! <hydro>/lhllc_c1_parts, default 31 = all):
//!   bit 1: S*: the L* and R* fluxes are blended with the cubic step
//!          w = 1/2 + 3x/4 - x^3/4, x = clamp(S*/(eps c_max), -1, 1) (exact HLLC choice
//!          outside |S*| < eps c_max; both star states are valid there, since
//!          S_L < -c < S* < c < S_R);
//!   bit 2: q = sqrt(1 + alpha ramp(p*/p - 1)), ramp the C1 max(0, x) over |x| < eps;
//!   bit 4: |u| and max/min in chi by their C1 versions over eps c_max (eps for the
//!          min with 1);
//!   bit 8: sensor: phi_eff = phi + (1 - phi) s with s = r^2 (r > 0, else 0),
//!          r = 2 J d2/(J^2 + d2^2): C1 at r = 0, smooth at |J| = |d2| (where r = 1, the
//!          pure checkerboard, s = 1), no chi-floor sqrt; the key value p of
//!          lhllc_cb_sensor still raises s to the power p;
//!   bit 16: phi = smax(phi, x1_phi_min) over eps.
//! The change to the dissipation is O(eps): e.g. the blended contact flux differs from
//! HLLC only where |S*| < eps c, where the upwind dissipation |S*| dU is itself below
//! eps c dU.

#include <cmath>

#include "utils/cb_sensor.hpp"

namespace hydro {

//! C1 |x|: |x| for |x| >= d, (x^2/d + d)/2 inside (value and slope match at |x| = d)
KOKKOS_INLINE_FUNCTION
Real C1Abs(const Real x, const Real d) {
  const Real ax = fabs(x);
  return (ax >= d) ? ax : 0.5*(x*x/d + d);
}
KOKKOS_INLINE_FUNCTION
Real C1Max(const Real a, const Real b, const Real d) {
  return 0.5*(a + b) + 0.5*C1Abs(a - b, d);
}
KOKKOS_INLINE_FUNCTION
Real C1Min(const Real a, const Real b, const Real d) {
  return 0.5*(a + b) - 0.5*C1Abs(a - b, d);
}

//! C1 odd-even fraction (bit 8): s = r^2 for r = 2 J d2/(J^2 + d2^2) > 0, else 0, with
//! J and d2 of the face-normal mass flux as in CBSensorChiFloor (utils/cb_sensor.hpp)
KOKKOS_INLINE_FUNCTION
Real CBSensorC1(const DvceArray5D<Real> &w0, const int m, const int k, const int j,
                const int i, const int ivx, const Real gain) {
  const int cdi = (ivx == IVX) ? 1 : 0;
  const int cdj = (ivx == IVY) ? 1 : 0;
  const int cdk = (ivx == IVZ) ? 1 : 0;
  const int nn = (cdi ? static_cast<int>(w0.extent(4)) :
                 (cdj ? static_cast<int>(w0.extent(3)) :
                        static_cast<int>(w0.extent(2)))) - 1;
  const int c0 = cdi*i + cdj*j + cdk*k;
  Real qs[4];
  for (int o = 0; o < 4; ++o) {
    int c = c0 - 2 + o;
    c = (c < 0) ? 0 : ((c > nn) ? nn : c);
    const int kk = cdk ? c : k, jj = cdj ? c : j, ii = cdi ? c : i;
    qs[o] = w0(m,IDN,kk,jj,ii)*w0(m,ivx,kk,jj,ii);
  }
  const Real d2 = qs[2] - qs[1];
  const Real jodd = 0.25*(2.0*d2 - (qs[1] - qs[0]) - (qs[3] - qs[2]));
  const Real den = jodd*jodd + d2*d2;
  const Real r = (den > 0.0) ? 2.0*jodd*d2/den : 0.0;
  Real s = (r > 0.0) ? r*r : 0.0;
  if (gain != 1.0 && s > 0.0) s = pow(s, gain);
  return s;
}

//----------------------------------------------------------------------------------------
//! \fn void LHLLC_C1
//! \brief LHLLC with C1 switches (see the file comment); same arguments as LHLLC

KOKKOS_INLINE_FUNCTION
void LHLLC_C1(TeamMember_t const &member, const EOS_Data &eos,
     const RegionIndcs &indcs,const DualArray1D<RegionSize> &size,const CoordData &coord,
     const int m, const int k, const int j, const int il, const int iu, const int ivx,
     const ScrArray2D<Real> &wl, const ScrArray2D<Real> &wr,
     const ScrArray2D<Real> &dl, const ScrArray2D<Real> &dr, DvceArray5D<Real> flx,
     const DvceArray5D<Real> &w0) {
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  Real gm1 = eos.gamma - 1.0;
  Real igm1 = 1.0/gm1;
  Real alpha = ((eos.gamma) + 1.0)/(2.0*(eos.gamma));

  const Real phi_min = (ivx == IVX) ? eos.x1_phi_min : 0.0;
  const Real cb_sensor = eos.cb_sensor;
  const bool cb_cont = eos.cb_cont;
  const Real eps = eos.c1_eps;
  const int parts = eos.c1_parts;
  const bool p_s = (parts & 1), p_q = (parts & 2), p_chi = (parts & 4),
             p_cb = (parts & 8), p_pm = (parts & 16);

  par_for_inner(member, il, iu, [&](const int i) {
    Real &wl_idn = wl(IDN,i);
    Real &wl_ivx = wl(ivx,i);
    Real &wl_ivy = wl(ivy,i);
    Real &wl_ivz = wl(ivz,i);

    Real &wr_idn = wr(IDN,i);
    Real &wr_ivx = wr(ivx,i);
    Real &wr_ivy = wr(ivy,i);
    Real &wr_ivz = wr(ivz,i);

    Real wl_ipr, wr_ipr, g1l, g1r, alpl, alpr;
    if (eos.IsGeneral()) {
      wl_ipr = dl(IDPR,i);
      wr_ipr = dr(IDPR,i);
      g1l = dl(IDG1,i);
      g1r = dr(IDG1,i);
      alpl = (g1l + 1.0)/(2.0*g1l);
      alpr = (g1r + 1.0)/(2.0*g1r);
    } else {
      wl_ipr = eos.IdealGasPressure(wl(IEN,i));
      wr_ipr = eos.IdealGasPressure(wr(IEN,i));
      g1l = eos.gamma;
      g1r = eos.gamma;
      alpl = alpha;
      alpr = alpha;
    }

    Real qa,qb,qc,qd,qe,qf;
    qa = eos.SoundSpeedFromP(wl_idn, wl_ipr, g1l);
    qb = eos.SoundSpeedFromP(wr_idn, wr_ipr, g1r);
    Real el, er;
    if (eos.IsGeneral()) {
      el = wl(IEN,i) + 0.5*wl_idn*(SQR(wl_ivx) + SQR(wl_ivy) + SQR(wl_ivz));
      er = wr(IEN,i) + 0.5*wr_idn*(SQR(wr_ivx) + SQR(wr_ivy) + SQR(wr_ivz));
    } else {
      el = wl_ipr*igm1 + 0.5*wl_idn*(SQR(wl_ivx) + SQR(wl_ivy) + SQR(wl_ivz));
      er = wr_ipr*igm1 + 0.5*wr_idn*(SQR(wr_ivx) + SQR(wr_ivy) + SQR(wr_ivz));
    }

    qc = 0.25*(wl_idn + wr_idn)*(qa + qb);
    qd = 0.5 * (wl_ipr + wr_ipr + (wl_ivx - wr_ivx) * qc);  // P_mid

    // coefficients for the low-Mach fix; cmax is a smooth function away from qa = qb
    // only through the max, which is harmless: it scales the bands, not the flux
    const Real cmax = fmax(qa,qb);
    const Real du = eps*cmax;
    Real chi;
    if (p_chi) {
      qe = C1Max(C1Abs(wl_ivx, du), C1Abs(wr_ivx, du), du);
      chi = C1Min(1.0, qe/cmax, eps);
    } else {
      qe = fmax(fabs(wl_ivx),fabs(wr_ivx));
      chi = fmin(1.0,qe/cmax);
    }
    Real phi;
    if (p_cb) {
      phi = chi*(2.0 - chi);
      if (cb_sensor > 0.0) {
        const Real s = CBSensorC1(w0, m, k, j, i, ivx, cb_sensor);
        phi = phi + (1.0 - phi)*s;
      }
    } else {
      if (cb_sensor > 0.0) chi = fmax(chi, CBSensorChiFloor(w0, m, k, j, i, ivx, chi,
                                                            cb_sensor, cb_cont));
      phi = chi*(2.0 - chi);
    }
    if (phi_min > 0.0) phi = p_pm ? C1Max(phi, phi_min, eps) : fmax(phi, phi_min);

    // PVRS wave-speed factors
    if (p_q) {
      const Real xl = qd/wl_ipr - 1.0, xr = qd/wr_ipr - 1.0;
      qe = sqrt(1.0 + alpl*0.5*(xl + C1Abs(xl, eps)));
      qf = sqrt(1.0 + alpr*0.5*(xr + C1Abs(xr, eps)));
    } else {
      qe = (qd <= wl_ipr) ? 1.0 : sqrt(1.0 + alpl * ((qd / wl_ipr) - 1.0));
      qf = (qd <= wr_ipr) ? 1.0 : sqrt(1.0 + alpr * ((qd / wr_ipr) - 1.0));
    }

    qc = wl_ivx - qa*qe;  // al
    qd = wr_ivx + qb*qf;  // ar

    qa = qd > 0.0 ? qd : 1.0e-20;   // bp
    qb = qc < 0.0 ? qc : -1.0e-20;  // bm

    qe = wl_ivx - qc; // vxl
    qf = wr_ivx - qd; // vxr

    qc = wl_ipr + qe*wl_idn*wl_ivx;  // tl
    qd = wr_ipr + qf*wr_idn*wr_ivx;  // tr

    Real ml =   wl_idn*qe;
    Real mr = -(wr_idn*qf);

    Real am = (qc - qd)/(ml + mr);
    qc = wl_ipr + phi*qe*wl_idn*wl_ivx;  // tl
    qd = wr_ipr + phi*qf*wr_idn*wr_ivx;  // tr
    Real cp = (ml*qd + mr*qc)/(ml + mr);
    cp = cp > 0.0 ? cp : 0.0;

    qe = wl_idn*(wl_ivx - qb);
    qf = wr_idn*(wr_ivx - qa);

    HydCons1D fl, fr;
    fl.d  = qe;
    fr.d  = qf;

    fl.mx = qe*wl_ivx + wl_ipr;
    fr.mx = qf*wr_ivx + wr_ipr;

    fl.my = qe*wl_ivy;
    fr.my = qf*wr_ivy;

    fl.mz = qe*wl_ivz;
    fr.mz = qf*wr_ivz;

    fl.e  = el*(wl_ivx - qb) + wl_ipr*wl_ivx;
    fr.e  = er*(wr_ivx - qa) + wr_ipr*wr_ivx;

    // flux weights: the L* branch (S* >= 0) and the R* branch, blended over |S*| < du
    Real w = (am >= 0.0) ? 1.0 : 0.0;
    if (p_s && du > 0.0 && fabs(am) < du) {
      const Real x = am/du;
      w = 0.5 + 0.75*x - 0.25*x*x*x;
    }
    qc = 0.0; qd = 0.0; qe = 0.0;
    if (w > 0.0) {
      qc += w*am/(am - qb);
      qe += w*(-qb/(am - qb));
    }
    if (w < 1.0) {
      qd += (1.0 - w)*(-am/(qa - am));
      qe += (1.0 - w)*qa/(qa - am);
    }

    flx(m,IDN,k,j,i) = qc*fl.d  + qd*fr.d;
    flx(m,ivx,k,j,i) = qc*fl.mx + qd*fr.mx + qe*cp;
    flx(m,ivy,k,j,i) = qc*fl.my + qd*fr.my;
    flx(m,ivz,k,j,i) = qc*fl.mz + qd*fr.mz;
    flx(m,IEN,k,j,i) = qc*fl.e  + qd*fr.e  + qe*cp*am;
  });
  return;
}
}  // namespace hydro
#endif  // HYDRO_RSOLVERS_LHLLC_C1_HYD_HPP_
