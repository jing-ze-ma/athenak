#ifndef HYDRO_RSOLVERS_LHLLC_HYD_HPP_
#define HYDRO_RSOLVERS_LHLLC_HYD_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file lhllc_hyd.hpp
//! \brief The HLLC Riemann solver for hydrodynamics with a low-Mach fix.  Only works for ideal gas EOS in hydrodynamics.
//!
//! REFERENCES:
//! - E.F. Toro, "Riemann Solvers and numerical methods for fluid dynamics", 2nd ed.,
//!   Springer-Verlag, Berlin, (1999) chpt. 10.
//!
//! - P. Batten, N. Clarke, C. Lambert, and D. M. Causon, "On the Choice of Wavespeeds
//!   for the HLLC Riemann Solver", SIAM J. Sci. & Stat. Comp. 18, 6, 1553-1570, (1997).
//!
//! - T. Minoshima and T. Miyoshi, "A low-dissipation HLLD approximate Riemann solver for a
//!   very wide range of Mach numbers", JCP, 446, 110639, (2021).

#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

namespace hydro {
//----------------------------------------------------------------------------------------
//! \fn void LHLLC
//! \brief The LHLLC Riemann solver for ideal gas hydrodynamics (use HLLE for isothermal)

KOKKOS_INLINE_FUNCTION
void LHLLC(TeamMember_t const &member, const EOS_Data &eos,
     const RegionIndcs &indcs,const DualArray1D<RegionSize> &size,const CoordData &coord,
     const int m, const int k, const int j, const int il, const int iu, const int ivx,
     const ScrArray2D<Real> &wl, const ScrArray2D<Real> &wr,
     const ScrArray2D<Real> &dl, const ScrArray2D<Real> &dr, DvceArray5D<Real> flx) {
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  Real gm1 = eos.gamma - 1.0;
  Real igm1 = 1.0/gm1;
  Real alpha = ((eos.gamma) + 1.0)/(2.0*(eos.gamma));

  par_for_inner(member, il, iu, [&](const int i) {
    //--- Step 1.  Create local references for L/R states (helps compiler vectorize)

    Real &wl_idn = wl(IDN,i);
    Real &wl_ivx = wl(ivx,i);
    Real &wl_ivy = wl(ivy,i);
    Real &wl_ivz = wl(ivz,i);

    Real &wr_idn = wr(IDN,i);
    Real &wr_ivx = wr(ivx,i);
    Real &wr_ivy = wr(ivy,i);
    Real &wr_ivz = wr(ivz,i);

    // general EOS: p and Gamma_1 were precomputed in ConsToPrim and
    // reconstructed here; Gamma_1 replaces gamma in the wave structure
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

    //--- Step 2.  Compute middle state estimates with PVRS (Toro 10.5.2)

    // define 6 registers used below
    Real qa,qb,qc,qd,qe,qf;
    qa = eos.SoundSpeedFromP(wl_idn, wl_ipr, g1l);
    qb = eos.SoundSpeedFromP(wr_idn, wr_ipr, g1r);
    Real el, er;
    if (eos.IsGeneral()) {
      // internal energy density is a reconstructed primitive
      el = wl(IEN,i) + 0.5*wl_idn*(SQR(wl_ivx) + SQR(wl_ivy) + SQR(wl_ivz));
      er = wr(IEN,i) + 0.5*wr_idn*(SQR(wr_ivx) + SQR(wr_ivy) + SQR(wr_ivz));
    } else {
      el = wl_ipr*igm1 + 0.5*wl_idn*(SQR(wl_ivx) + SQR(wl_ivy) + SQR(wl_ivz));
      er = wr_ipr*igm1 + 0.5*wr_idn*(SQR(wr_ivx) + SQR(wr_ivy) + SQR(wr_ivz));
    }
      
//      qc = sqrt(SQR(wl_ivx) + SQR(wl_ivy) + SQR(wl_ivz))/qa;
//      qd = sqrt(SQR(wr_ivx) + SQR(wr_ivy) + SQR(wr_ivz))/qb;
//      Real chi = fmax(qc,qd);
//      Real phi = chi*(2.0-chi);
//      
//    qc = fmin(wl_ivx,wr_ivx);
//    qd = fmax(wl_ivx,wr_ivx);
//    qe = fmax(qa,qb);
//    Real sl = qc - qe;
//    Real sr = qd - qe;
//    
//    qc = sl - wl_ivx;
//    qd = sr - wr_ivx;
//      qe = wl_idn*qc - wr_idn*qd;
//      Real ustar = (wr_ipr-wl_ipr+wl_idn*wl_ivx*qc-wr_idn*wr_ivx*qd)/qe;
//      
//      Real rhobar = 0.5*(wl_idn + wr_idn);
//      Real abar = 0.5*(qa + qb);
//      Real pstar = 0.5*(wl_ipr+wr_ipr) - 0.5*phi*(wr_ivx-wl_ivx)*rhobar*abar;
//      
//      Real ustaralpha;
//      if ((sl < 0.0) and (ustar >= 0.0)) {
//          ustaralpha = ustar*wl_idn*qc/(sl-ustar);
//      }
      
      
    qc = 0.25*(wl_idn + wr_idn)*(qa + qb);  // average density * average sound speed
    qd = 0.5 * (wl_ipr + wr_ipr + (wl_ivx - wr_ivx) * qc);  // P_mid
    
    // coefficients for the low-Mach fix
    qe = fmax(fabs(wl_ivx),fabs(wr_ivx)); // maximum speed
    qf = fmax(qa,qb); // maximum sound speed
    Real chi = fmin(1.0,qe/qf);
    Real phi = chi * (2.0 - chi);

    //--- Step 3.  Compute sound speed in L,R

    qe = (qd <= wl_ipr) ? 1.0 : sqrt(1.0 + alpl * ((qd / wl_ipr) - 1.0));  // ql
    qf = (qd <= wr_ipr) ? 1.0 : sqrt(1.0 + alpr * ((qd / wr_ipr) - 1.0));  // qr

    //--- Step 4.  Compute the max/min wave speeds based on L/R

    qc = wl_ivx - qa*qe;  // al
    qd = wr_ivx + qb*qf;  // ar

    // following min/max set to TINY_NUMBER to fix bug found in converging supersonic flow
    qa = qd > 0.0 ? qd : 1.0e-20;   // bp
    qb = qc < 0.0 ? qc : -1.0e-20;  // bm

    //--- Step 5. Compute the contact wave speed and pressure

    qe = wl_ivx - qc; // vxl
    qf = wr_ivx - qd; // vxr

    qc = wl_ipr + qe*wl_idn*wl_ivx;  // tl
    qd = wr_ipr + qf*wr_idn*wr_ivx;  // tr

    Real ml =   wl_idn*qe;
    Real mr = -(wr_idn*qf);

    // Determine the contact wave speed...
    Real am = (qc - qd)/(ml + mr);
    // ...and the pressure at the contact surface
    // Low-Mach fix (Minoshima & Miyoshi 2021)
    qc = wl_ipr + phi*qe*wl_idn*wl_ivx;  // tl
    qd = wr_ipr + phi*qf*wr_idn*wr_ivx;  // tr
    Real cp = (ml*qd + mr*qc)/(ml + mr);
    cp = cp > 0.0 ? cp : 0.0;

    //--- Step 6. Compute L/R fluxes along the line bm (qb), bp (qa)

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

    //--- Step 8. Compute flux weights or scales

    if (am >= 0.0) {
      qc =  am/(am - qb);
      qd = 0.0;
      qe = -qb/(am - qb);
    } else {
      qc =  0.0;
      qd = -am/(qa - am);
      qe =  qa/(qa - am);
    }

    //--- Step 9. Compute the LHLLC flux at interface, including weighted contribution
    // of the flux along the contact

    flx(m,IDN,k,j,i) = qc*fl.d  + qd*fr.d;
    flx(m,ivx,k,j,i) = qc*fl.mx + qd*fr.mx + qe*cp;
    flx(m,ivy,k,j,i) = qc*fl.my + qd*fr.my;
    flx(m,ivz,k,j,i) = qc*fl.mz + qd*fr.mz;
    flx(m,IEN,k,j,i) = qc*fl.e  + qd*fr.e  + qe*cp*am;
  });
  return;
}
} // namespace hydro
#endif // HYDRO_RSOLVERS_LHLLC_HYD_HPP_
