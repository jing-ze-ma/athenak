#ifndef HYDRO_HYDRO_AM_HPP_
#define HYDRO_HYDRO_AM_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hydro_am.hpp
//! \brief cs-open-1018: <hydro>/am_transport, the axial angular-momentum-conservative
//! transport correction of the hydro on the cubed sphere and on spherical polar, as a
//! GENERIC hydro option (the deep_hot_jupiter_rt pgen carries the same correction as
//! problem/am_transport since 5a119ba9 / f121c61f / 53e55b4c, AmCell there; this is its
//! hydro-only transport part, without the rotating-frame Coriolis part).
//!
//! The zonal (about z) part of every cell's update is corrected by T_cons - T_code:
//! T_code the torque the update applied (the x1/x2/x3 flux divergences incl. the x1
//! r-factor source, and the angular metric source of SrcTermsGnomonicEquiangle /
//! SrcTermsSphericalPolarHydro, re-evaluated here), T_cons the torque of the same face
//! fluxes taken as ANGULAR-MOMENTUM fluxes with exact face first moments (cs: r_f rhat x F
//! on the radial faces, (int_f x dA) x F on the planar xi / eta faces; sp: F_phi A
//! Rbar_f).  sum_c T_cons is the boundary flux alone, so L_z changes only by boundary
//! fluxes and real torques (to round-off).  Applied as a force along (-y, x, 0): no
//! radial or meridional part; the total energy is not touched.

#include "athena.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "eos/eos.hpp"

namespace hydro {

struct HydroAmV {
  DvceArray5D<Real> w0, wder, f1, f2, f3;
  DvceArray4D<Real> a1, a2, a3, vol, xovr, yovr;
  DvceArray3D<Real> ccell, scell;
  DvceArray2D<Real> x1v, x1f, x2v, x2f;
  DualArray1D<RegionSize> size;
  DualArray1D<int> panel;
  EOS_Data eos;
  int is, js, ks, nx2, nx3;
  bool cs, gen;
};

//! \brief the zonal torque correction dT = T_cons - T_code of cell (m,k,j,i) (already x
//! vol) and the coefficients that turn it into momentum density: m2 += dT c2 bdt,
//! m3 += dT c3 bdt.  Returns false where R = 0 (no zonal direction).
KOKKOS_INLINE_FUNCTION
bool HydroAmTorque(const HydroAmV &v, const int m, const int k, const int j,
                   const int i, Real &dtq, Real &c2o, Real &c3o) {
  const Real dv = v.vol(m,k,j,i);
  const Real ia = 1.0/dv;
  const Real r = v.x1v(m,i);
  const Real rl = v.x1f(m,i), rr = v.x1f(m,i+1);
  const Real fac = (rr - rl)/(rr + rl);
  const Real rho = v.w0(m,IDN,k,j,i);
  const Real v2 = v.w0(m,IVY,k,j,i), v3 = v.w0(m,IVZ,k,j,i);
  if (v.cs) {
    const Real x2l = v.size.d_view(m).x2min, x2u = v.size.d_view(m).x2max;
    const Real x3l = v.size.d_view(m).x3min, x3u = v.size.d_view(m).x3max;
    const Real xc = 0.25*M_PI*CellCenterX(j-v.js, v.nx2, x2l, x2u);
    const Real ec = 0.25*M_PI*CellCenterX(k-v.ks, v.nx3, x3l, x3u);
    const Real xl = 0.25*M_PI*LeftEdgeX(j-v.js, v.nx2, x2l, x2u);
    const Real xr = 0.25*M_PI*LeftEdgeX(j+1-v.js, v.nx2, x2l, x2u);
    const Real el = 0.25*M_PI*LeftEdgeX(k-v.ks, v.nx3, x3l, x3u);
    const Real er = 0.25*M_PI*LeftEdgeX(k+1-v.ks, v.nx3, x3l, x3u);
    const int pnl = v.panel.d_view(m);
    Real qc[3], e1[3], e2[3];
    cubed_sphere::PanelToCart(pnl, xc, ec, qc);
    cubed_sphere::PanelTangents(pnl, xc, ec, e1, e2);
    const Real qn = 1.0/sqrt(qc[0]*qc[0] + qc[1]*qc[1] + qc[2]*qc[2]);
    const Real rh0 = qc[0]*qn, rh1 = qc[1]*qn;
    const Real xx = r*rh0, yy = r*rh1;
    const Real wc = xx*xx + yy*yy;
    if (!(wc > 0.0)) return false;
    c2o = (-yy*e1[0] + xx*e1[1])/(wc*dv);
    c3o = (-yy*e2[0] + xx*e2[1])/(wc*dv);
    // T_code: the flux divergences and the angular source of SrcTermsGnomonicEquiangle
    const Real c = v.ccell(m,k,j);
    const Real s2 = 1.0 - c*c;
    Real dx1[2], dx2[2], dx3[2];
    for (int n=0; n<2; ++n) {
      const int nm = IM2 + n;
      dx1[n] = -(v.f1(m,nm,k,j,i+1)*v.a1(m,k,j,i+1) - v.f1(m,nm,k,j,i)*v.a1(m,k,j,i))*ia
               - fac*(v.f1(m,nm,k,j,i)*v.a1(m,k,j,i)
                      + v.f1(m,nm,k,j,i+1)*v.a1(m,k,j,i+1))*ia;
      dx2[n] = -(v.f2(m,nm,k,j+1,i)*v.a2(m,k,j+1,i) - v.f2(m,nm,k,j,i)*v.a2(m,k,j,i))*ia;
      dx3[n] = -(v.f3(m,nm,k+1,j,i)*v.a3(m,k+1,j,i) - v.f3(m,nm,k,j,i)*v.a3(m,k,j,i))*ia;
    }
    const Real pr = v.gen ? v.wder(m,IDPR,k,j,i) : v.eos.Pressure(rho, v.w0(m,IEN,k,j,i));
    const Real sine = v.scell(m,k,j);
    const Real sine2 = sine*sine;
    const Real src2 = v.xovr(m,k,j,i)*(pr + rho*v3*v3*sine2);
    const Real src3 = v.yovr(m,k,j,i)*(pr + rho*v2*v2*sine2);
    const Real d2[4] = {dx1[0], dx2[0], dx3[0], src2};
    const Real d3[4] = {dx1[1], dx2[1], dx3[1], src3};
    Real tcode = 0.0;
    for (int n=0; n<4; ++n) {
      const Real u2 = (d2[n] - c*d3[n])/s2, u3 = (d3[n] - c*d2[n])/s2;
      tcode += (xx*(u2*e1[1] + u3*e2[1]) - yy*(u2*e1[0] + u3*e2[0]))*dv;
    }
    // T_cons, radial faces: arm r_f rhat_c, (rhat x F)_z, F on the cell's basis
    Real t = 0.0;
    for (int n=0; n<2; ++n) {
      const int ii = i + n;
      const Real g2 = v.f1(m,IM2,k,j,ii), g3 = v.f1(m,IM3,k,j,ii);
      const Real fa = (g2 - c*g3)/s2, fb = (g3 - c*g2)/s2;
      const Real tr = v.a1(m,k,j,ii)*v.x1f(m,ii)*(rh0*(fa*e1[1] + fb*e2[1])
                                                  - rh1*(fa*e1[0] + fb*e2[0]));
      t += n ? -tr : tr;
    }
    // planar xi / eta faces: (int_f x dA) x F, int_f x dA = (rr^3 - rl^3)/3 tan(psi/2)
    // (q_a + q_b), q_a/q_b the directions of the face's two ends, psi the angle between
    const Real mom = (rr*rr*rr - rl*rl*rl)/3.0;
    for (int dir=0; dir<2; ++dir) {
      for (int n=0; n<2; ++n) {
        Real xf, ef, xa, ea, xb, eb, g2, g3;
        if (dir == 0) {
          xf = n ? xr : xl;  ef = ec;  xa = xf;  ea = el;  xb = xf;  eb = er;
          g2 = v.f2(m,IM2,k,j+n,i);  g3 = v.f2(m,IM3,k,j+n,i);
        } else {
          xf = xc;  ef = n ? er : el;  xa = xl;  ea = ef;  xb = xr;  eb = ef;
          g2 = v.f3(m,IM2,k+n,j,i);  g3 = v.f3(m,IM3,k+n,j,i);
        }
        Real qa[3], qb[3], b1[3], b2[3];
        cubed_sphere::PanelToCart(pnl, xa, ea, qa);
        cubed_sphere::PanelToCart(pnl, xb, eb, qb);
        cubed_sphere::PanelTangents(pnl, xf, ef, b1, b2);
        const Real cf = b1[0]*b2[0] + b1[1]*b2[1] + b1[2]*b2[2];
        const Real sf = 1.0 - cf*cf;
        const Real fa = (g2 - cf*g3)/sf, fb = (g3 - cf*g2)/sf;
        const Real cx = qa[1]*qb[2] - qa[2]*qb[1];
        const Real cy = qa[2]*qb[0] - qa[0]*qb[2];
        const Real cz = qa[0]*qb[1] - qa[1]*qb[0];
        const Real dot = qa[0]*qb[0] + qa[1]*qb[1] + qa[2]*qb[2];
        const Real th = sqrt(cx*cx + cy*cy + cz*cz)/(1.0 + dot);   // tan(psi/2)
        const Real mx = mom*th*(qa[0] + qb[0]), my = mom*th*(qa[1] + qb[1]);
        const Real tf = mx*(fa*b1[1] + fb*b2[1]) - my*(fa*b1[0] + fb*b2[0]);
        t += n ? -tf : tf;
      }
    }
    dtq = t - tcode;
    return true;
  }
  // spherical polar: R = r sin(theta), the phi momentum is IM3 (orthonormal)
  const Real th = v.x2v(m,j);
  const Real thl = v.x2f(m,j), thr = v.x2f(m,j+1);
  const Real rc = r*sin(th);
  if (!(rc > 0.0)) return false;
  c2o = 0.0;
  c3o = ia/rc;
  const Real sl = sin(thl), sr = sin(thr);
  const Real fl1 = v.f1(m,IM3,k,j,i)*v.a1(m,k,j,i);
  const Real fr1 = v.f1(m,IM3,k,j,i+1)*v.a1(m,k,j,i+1);
  const Real fl2 = v.f2(m,IM3,k,j,i)*v.a2(m,k,j,i);
  const Real fr2 = v.f2(m,IM3,k,j+1,i)*v.a2(m,k,j+1,i);
  const Real fl3 = v.f3(m,IM3,k,j,i)*v.a3(m,k,j,i);
  const Real fr3 = v.f3(m,IM3,k+1,j,i)*v.a3(m,k+1,j,i);
  const Real tcode = rc*(-(fr1 - fl1) - fac*(fl1 + fr1)) + rc*(-(fr2 - fl2))
                     + rc*(-(fr3 - fl3)) + rc*(-v.yovr(m,k,j,i)*(fl2 + fr2));
  const Real cl = cos(thl), cr = cos(thr);
  const Real sin2int = 0.5*(thr - thl) - 0.25*(sin(2.0*thr) - sin(2.0*thl));
  const Real dcos = cl - cr;
  const Real r3 = (rr*rr*rr - rl*rl*rl)/3.0, r2h = 0.5*(rr*rr - rl*rl);
  const Real rb1 = (dcos > 0.0) ? sin2int/dcos : 0.0;
  const Real rb2 = r3/r2h;
  const Real rb3 = (thr > thl) ? r3*dcos/(r2h*(thr - thl)) : 0.0;
  const Real tcons = -(fr1*rr*rb1 - fl1*rl*rb1) - (fr2*rb2*sr - fl2*rb2*sl)
                     - (fr3 - fl3)*rb3;
  dtq = tcons - tcode;
  return true;
}

}  // namespace hydro
#endif  // HYDRO_HYDRO_AM_HPP_
