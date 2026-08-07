#ifndef HYDRO_RSOLVERS_AUSMPUP_HYD_HPP_
#define HYDRO_RSOLVERS_AUSMPUP_HYD_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file ausmpup_hyd.hpp
//! \brief The AUSM+-up Riemann solver for hydrodynamics.  Only works for ideal gas EOS in hydrodynamics.
//!
//! REFERENCES:
//! - M.-S. Liou, "A sequel to AUSM, Part II: AUSM +-up for all speeds", JCP, 214, 137, (2006).

#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()
#include "hydro/hydro.hpp"

namespace hydro {
//----------------------------------------------------------------------------------------
//! \fn void AUSMPUP
//! \brief The AUSMPUP Riemann solver for ideal gas hydrodynamics (use HLLE for isothermal)

KOKKOS_INLINE_FUNCTION
void AUSMPUP(TeamMember_t const &member, const EOS_Data &eos,
     const RegionIndcs &indcs,const DualArray1D<RegionSize> &size,const CoordData &coord,
     const int m, const int k, const int j, const int il, const int iu, const int ivx,
     const ScrArray2D<Real> &wl, const ScrArray2D<Real> &wr,
     const ScrArray2D<Real> &dl, const ScrArray2D<Real> &dr, DvceArray5D<Real> flx) {
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  Real gm1 = eos.gamma - 1.0;
  Real igm1 = 1.0/gm1;
  Real alpha = 2.0*gm1/(eos.gamma+1.0);
    
  const Real Mref = 1.0;

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
      alpl = 2.0*(g1l - 1.0)/(g1l + 1.0);
      alpr = 2.0*(g1r - 1.0)/(g1r + 1.0);
    } else {
      wl_ipr = eos.IdealGasPressure(wl(IEN,i));
      wr_ipr = eos.IdealGasPressure(wr(IEN,i));
      g1l = eos.gamma;
      g1r = eos.gamma;
      alpl = alpha;
      alpr = alpha;
    }

    //--- Step 2.  Dimensionless coefficients (Liou 3.3)

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
      
    // Liou 2.3
    qc = sqrt(alpl * (el+wl_ipr)/wl_idn);  // a_star left
    qd = sqrt(alpr * (er+wr_ipr)/wr_idn);  // a_star right
    qe = SQR(qc)/fmax(qc,wl_ivx);
    qf = SQR(qd)/fmax(qd,-wr_ivx);
    qc = fmin(qe,qf);  // average sound speed
      
    qd = 0.5 * (wl_idn + wr_idn);  // average density
    
    qe = sqrt((SQR(wl_ivx)+SQR(wr_ivx))/(2.0*SQR(qc))); // M_bar
    qf = fmin(1.0,fmax(qe,1.0e-13)); // M_o
    Real fa = qf * (2.0 - qf);
    qf = fmin(1.0,fmax(qe,Mref)); // M_o^p (Edelmann et al. 2021)
    Real fap = qf * (2.0 - qf);
    Real alp = 3.0/16.0*(-4.0+5.0*SQR(fa));
    Real bet = 1.0/8.0;
      
    //--- Step 3.  Compute the mass flux and pressure flux
      
    Real Mhalf = - 0.25/fap * fmax(1.0-SQR(qe),0.0) * (wr_ipr - wl_ipr) / (qd*SQR(qc));
    Real pu = - 0.75 * (wl_idn + wr_idn) * fa * qc * (wr_ivx - wl_ivx);
      
    qa = wl_ivx/qc; // left Mach
    qb = 0.5*(qa+fabs(qa));
    qd = 0.25*SQR(qa+1.0);
    qe = -0.25*SQR(qa-1.0);
    qf = (fabs(qa) >= 1.0) ? qb : (qd*(1.0-16.0*bet*qe));
    Mhalf += qf;
    qf = (fabs(qa) >= 1.0) ? (qb/qa) : (qd*(2.0-qa-16.0*alp*qa*qe));
    pu *= qf;
    Real phalf = qf*wl_ipr;
    
    qa = wr_ivx/qc; // right Mach
    qb = 0.5*(qa-fabs(qa));
    qd = -0.25*SQR(qa-1.0);
    qe = 0.25*SQR(qa+1.0);
    qf = (fabs(qa) >= 1.0) ? qb : (qd*(1.0+16.0*bet*qe));
    Mhalf += qf;
    qf = (fabs(qa) >= 1.0) ? (qb/qa) : (qd*(-2.0-qa+16.0*alp*qa*qe));
    pu *= qf;
    phalf += qf*wr_ipr + pu; // pressure flux
                                    
    qa = (Mhalf > 0.0) ? wl_idn : wr_idn;
    qa *= qc * Mhalf; // mass flux

    //--- Step 4. Compute L/R fluxes

    HydCons1D fl, fr;
    fl.d  = 1.0;
    fr.d  = 1.0;

    fl.mx = wl_ivx;
    fr.mx = wr_ivx;

    fl.my = wl_ivy;
    fr.my = wr_ivy;

    fl.mz = wl_ivz;
    fr.mz = wr_ivz;

    fl.e  = (el + wl_ipr)/wl_idn;
    fr.e  = (er + wr_ipr)/wr_idn;

    //--- Step 5. Compute flux weights or scales

    qc = (qa > 0.0);
    qd = 1.0 - qc;

    //--- Step 6. Compute the AUSM+-up flux at interface

    flx(m,IDN,k,j,i) = qa*(qc*fl.d  + qd*fr.d);
    flx(m,ivx,k,j,i) = qa*(qc*fl.mx + qd*fr.mx) + phalf;
    flx(m,ivy,k,j,i) = qa*(qc*fl.my + qd*fr.my);
    flx(m,ivz,k,j,i) = qa*(qc*fl.mz + qd*fr.mz);
    flx(m,IEN,k,j,i) = qa*(qc*fl.e  + qd*fr.e);
  });
  return;
}
} // namespace hydro
#endif // HYDRO_RSOLVERS_AUSMPUP_HYD_HPP_
