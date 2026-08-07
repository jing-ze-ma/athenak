#ifndef HYDRO_RSOLVERS_HLLCLM_HYD_HPP_
#define HYDRO_RSOLVERS_HLLCLM_HYD_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hllclm_hyd.hpp
//! \brief The HLLC Riemann solver for hydrodynamics with a low-Mach fix in the central form.  Only works for ideal gas EOS in hydrodynamics.
//!
//! REFERENCES:
//! - E.F. Toro, "Riemann Solvers and numerical methods for fluid dynamics", 2nd ed.,
//!   Springer-Verlag, Berlin, (1999) chpt. 10.
//!
//! - P. Batten, N. Clarke, C. Lambert, and D. M. Causon, "On the Choice of Wavespeeds
//!   for the HLLC Riemann Solver", SIAM J. Sci. & Stat. Comp. 18, 6, 1553-1570, (1997).
//!
//! - N. Fleischmann, S. Adami, and N. A. Adams, "A shock-stable modification of the HLLC Riemann solver with reduced numerical dissipation", JCP, 423, 109762, (2020).

#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

namespace hydro {
//----------------------------------------------------------------------------------------
//! \fn void HLLCLM
//! \brief The HLLCLM Riemann solver for ideal gas hydrodynamics (use HLLE for isothermal)

KOKKOS_INLINE_FUNCTION
void HLLCLM(TeamMember_t const &member, const EOS_Data &eos,
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
    qc = 0.25*(wl_idn + wr_idn)*(qa + qb);  // average density * average sound speed
    qd = 0.5 * (wl_ipr + wr_ipr + (wl_ivx - wr_ivx) * qc);  // P_mid
      
    // coefficients for the low-Mach fix
    qe = fmax(fabs(wl_ivx/qa),fabs(wr_ivx/qb)); // maximum local Mach number
      Real phi = 1.0; // sin(fmin(1.0,qe/0.01)*M_PI/2.0);

    //--- Step 3.  Compute sound speed in L,R

    qe = (qd <= wl_ipr) ? 1.0 : sqrt(1.0 + alpl * ((qd / wl_ipr) - 1.0));  // ql
    qf = (qd <= wr_ipr) ? 1.0 : sqrt(1.0 + alpr * ((qd / wr_ipr) - 1.0));  // qr

    //--- Step 4.  Compute the max/min wave speeds based on L/R

    qc = wl_ivx - qa*qe;  // al
    qd = wr_ivx + qb*qf;  // ar

    //--- Step 5. Compute the contact wave speed and pressure

    qe = wl_ivx - qc; // vxl
    qf = wr_ivx - qd; // vxr

    qa = wl_ipr + qe*wl_idn*wl_ivx;  // tl
    qb = wr_ipr + qf*wr_idn*wr_ivx;  // tr

    Real ml =   wl_idn*qe;
    Real mr = -(wr_idn*qf);

    // Determine the contact wave speed...
    Real am = (qa - qb)/(ml + mr);
    // ...and the pressure at the contact surface
    Real cp = (ml*qb + mr*qa)/(ml + mr);
    cp = cp > 0.0 ? cp : 0.0;
      
    //--- Step 6. Compute L/R U_star
      
    qa = wl_idn*qe/(am - qc);
    qb = wr_idn*qf/(am - qd);
        
    HydCons1D usl, usr;
      
    usl.d = qa;
    usr.d = qb;
      
    usl.mx = qa*am;
    usr.mx = qb*am;
      
    usl.my = qa*wl_ivy;
    usr.my = qb*wl_ivy;
      
    usl.mz = qa*wl_ivz;
    usr.mz = qb*wl_ivz;
      
    usl.e = qe/(am - qc)*(el + (am - wl_ivx)*(wl_idn*am - wl_ipr/qe));
    usr.e = qf/(am - qd)*(er + (am - wr_ivx)*(wr_idn*am - wr_ipr/qf));

    //--- Step 7. Compute L/R fluxes

    qe = wl_idn*wl_ivx;
    qf = wr_idn*wr_ivx;

    HydCons1D fl, fr, fs;
    fl.d  = qe;
    fr.d  = qf;

    fl.mx = qe*wl_ivx + wl_ipr;
    fr.mx = qf*wr_ivx + wr_ipr;

    fl.my = qe*wl_ivy;
    fr.my = qf*wr_ivy;

    fl.mz = qe*wl_ivz;
    fr.mz = qf*wr_ivz;

    fl.e  = (el + wl_ipr)*wl_ivx;
    fr.e  = (er + wr_ipr)*wr_ivx;

    //--- Step 8. Compute flux star
      
    // Low-Mach fix (Fleischmann et al. 2020)
    qa = phi*qc;
    qb = phi*qd;

    fs.d = 0.5*(fl.d+fr.d) + 0.5*(qa*(usl.d - wl_idn) + fabs(am)*(usl.d - usr.d) + qb*(usr.d - wr_idn));
    fs.mx = 0.5*(fl.mx+fr.mx) + 0.5*(qa*(usl.mx - wl_idn*wl_ivx) + fabs(am)*(usl.mx - usr.mx) + qb*(usr.mx - wr_idn*wr_ivx));
    fs.my = 0.5*(fl.my+fr.my) + 0.5*(qa*(usl.my - wl_idn*wl_ivy) + fabs(am)*(usl.my - usr.my) + qb*(usr.my - wr_idn*wr_ivy));
    fs.mz = 0.5*(fl.mz+fr.mz) + 0.5*(qa*(usl.mz - wl_idn*wl_ivz) + fabs(am)*(usl.mz - usr.mz) + qb*(usr.mz - wr_idn*wr_ivz));
    fs.e = 0.5*(fl.e+fr.e) + 0.5*(qa*(usl.e - el) + fabs(am)*(usl.e - usr.e) + qb*(usr.e - er));
      
    //--- Step 9. Compute flux weights or scales
      
    qa = (qc >= 0.0);
    qb = (qd <= 0.0);
    qe = 1.0 - qa - qb;

    //--- Step 10. Compute the HLLCLM flux at interface

    flx(m,IDN,k,j,i) = qa*fl.d  + qb*fr.d  + qe*fs.d;
    flx(m,ivx,k,j,i) = qa*fl.mx + qb*fr.mx + qe*fs.mx;
    flx(m,ivy,k,j,i) = qa*fl.my + qb*fr.my + qe*fs.my;
    flx(m,ivz,k,j,i) = qa*fl.mz + qb*fr.mz + qe*fs.mz;
    flx(m,IEN,k,j,i) = qa*fl.e  + qb*fr.e  + qe*fs.e;
  });
  return;
}
} // namespace hydro
#endif // HYDRO_RSOLVERS_HLLCLM_HYD_HPP_
