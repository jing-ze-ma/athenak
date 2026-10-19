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
#include "coordinates/cell_locations.hpp"
#include "eos/eos.hpp"
#include <string>

class MeshBlockPack;
class ParameterInput;
class Driver;
class MeshBoundaryValuesCC;

namespace hydro {

struct HydroAmV {
  DvceArray5D<Real> w0, wder, f1, f2, f3;
  DvceArray5D<Real> bcc;   // MHD: the cell field (cs: orthonormal triple, see hydro_am.cpp)
  bool mhd = false;
  DvceArray4D<Real> a1, a2, a3, vol, xovr, yovr;
  DvceArray3D<Real> ccell, scell;
  DvceArray2D<Real> x1v, x1f, x2v, x2f, x3v, x3f;
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
    // the angular source of SrcTermsGnomonicEquiangleImpl, incl. its Maxwell part (MHD:
    // bcc the orthonormal triple, |B|^2 = b1^2 + b2^2 + b3^2, (s B^xi)^2 = (s b2 - c b3)^2,
    // (s B^eta)^2 = b3^2)
    Real ptot = pr, bb2 = 0.0, bb3 = 0.0;
    if (v.mhd) {
      const Real b1 = v.bcc(m,IBX,k,j,i), b2 = v.bcc(m,IBY,k,j,i), b3 = v.bcc(m,IBZ,k,j,i);
      bb2 = (sine*b2 - c*b3)*(sine*b2 - c*b3);
      bb3 = b3*b3;
      ptot += 0.5*(b1*b1 + b2*b2 + b3*b3);
    }
    const Real src2 = v.xovr(m,k,j,i)*(ptot + rho*v3*v3*sine2 - bb3);
    const Real src3 = v.yovr(m,k,j,i)*(ptot + rho*v2*v2*sine2 - bb2);
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


KOKKOS_INLINE_FUNCTION
void HydroAmCross(const Real *a, const Real *b, Real *c) {
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

//! \brief cs-open-1018 <hydro>/am_vector: the full VECTOR form of HydroAmTorque.  dT =
//! T_cons - T_code (3 Cartesian components, x vol), the cell-centre position rc
//! (= x1v rhat_c, the arm of L = sum V rc x m) and the cell's tangential basis (b1, b2:
//! cs e_xi, e_eta, covariant momenta raised with cc = e_xi.e_eta; sp e_theta, e_phi,
//! cc = 0).  T_code = V rc x (the tangential momentum change rate of the flux divergences
//! and of the angular metric sources); T_cons = the face fluxes taken as momentum-flux
//! vectors on the face's basis (at its centre), with the faces' exact first moments
//! int_f x dA (cs: planar xi/eta faces as HydroAmTorque, radial faces A r_f rhat_c;
//! sp: radial caps, theta cones and phi half-planes analytically).  The radial part of the
//! flux exerts no torque about the origin (x || rhat on every face point) and is left out.
//! Each face's (moment, flux) pair is formed from the face's own data only, so
//! sum_c T_cons is the boundary torque alone (to round-off).

KOKKOS_INLINE_FUNCTION
bool HydroAmVecTorque(const HydroAmV &v, const int m, const int k, const int j,
                      const int i, Real *dT, Real *rc, Real *b1, Real *b2, Real &cc) {
  const Real dv = v.vol(m,k,j,i);
  const Real ia = 1.0/dv;
  const Real r = v.x1v(m,i);
  const Real rl = v.x1f(m,i), rr = v.x1f(m,i+1);
  const Real fac = (rr - rl)/(rr + rl);
  const Real rho = v.w0(m,IDN,k,j,i);
  const Real v2 = v.w0(m,IVY,k,j,i), v3 = v.w0(m,IVZ,k,j,i);
  const Real pr = v.gen ? v.wder(m,IDPR,k,j,i) : v.eos.Pressure(rho, v.w0(m,IEN,k,j,i));
  const Real mom = (rr*rr*rr - rl*rl*rl)/3.0;
  Real tc[3] = {0.0, 0.0, 0.0}, tk[3] = {0.0, 0.0, 0.0};
  Real rh[3];
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
    Real qc[3];
    cubed_sphere::PanelToCart(pnl, xc, ec, qc);
    cubed_sphere::PanelTangents(pnl, xc, ec, b1, b2);
    const Real qn = 1.0/sqrt(qc[0]*qc[0] + qc[1]*qc[1] + qc[2]*qc[2]);
    for (int q = 0; q < 3; ++q) {rh[q] = qc[q]*qn; rc[q] = r*rh[q];}
    const Real c = v.ccell(m,k,j);
    cc = c;
    const Real s2 = 1.0 - c*c;
    // T_code
    Real d2 = 0.0, d3 = 0.0;
    for (int n=0; n<2; ++n) {
      const int nm = IM2 + n;
      Real dd = -(v.f1(m,nm,k,j,i+1)*v.a1(m,k,j,i+1) - v.f1(m,nm,k,j,i)*v.a1(m,k,j,i))*ia
                - fac*(v.f1(m,nm,k,j,i)*v.a1(m,k,j,i)
                       + v.f1(m,nm,k,j,i+1)*v.a1(m,k,j,i+1))*ia
                - (v.f2(m,nm,k,j+1,i)*v.a2(m,k,j+1,i) - v.f2(m,nm,k,j,i)*v.a2(m,k,j,i))*ia
                - (v.f3(m,nm,k+1,j,i)*v.a3(m,k+1,j,i) - v.f3(m,nm,k,j,i)*v.a3(m,k,j,i))*ia;
      if (n == 0) {d2 = dd;} else {d3 = dd;}
    }
    const Real sine = v.scell(m,k,j), sine2 = sine*sine;
    {   // the angular source incl. the Maxwell part (as HydroAmTorque)
      Real ptot = pr, bb2 = 0.0, bb3 = 0.0;
      if (v.mhd) {
        const Real b1 = v.bcc(m,IBX,k,j,i), b2 = v.bcc(m,IBY,k,j,i);
        const Real b3 = v.bcc(m,IBZ,k,j,i);
        bb2 = (sine*b2 - c*b3)*(sine*b2 - c*b3);
        bb3 = b3*b3;
        ptot += 0.5*(b1*b1 + b2*b2 + b3*b3);
      }
      d2 += v.xovr(m,k,j,i)*(ptot + rho*v3*v3*sine2 - bb3);
      d3 += v.yovr(m,k,j,i)*(ptot + rho*v2*v2*sine2 - bb2);
    }
    {
      const Real u2 = (d2 - c*d3)/s2, u3 = (d3 - c*d2)/s2;
      Real dp[3];
      for (int q = 0; q < 3; ++q) {dp[q] = (u2*b1[q] + u3*b2[q])*dv;}
      HydroAmCross(rc, dp, tk);
    }
    // T_cons, radial faces: arm r_f rhat_c, F on the cell's basis, x area
    for (int n=0; n<2; ++n) {
      const int ii = i + n;
      const Real g2 = v.f1(m,IM2,k,j,ii), g3 = v.f1(m,IM3,k,j,ii);
      const Real fa = (g2 - c*g3)/s2, fb = (g3 - c*g2)/s2;
      Real fv[3], mv[3], t[3];
      for (int q = 0; q < 3; ++q) {
        fv[q] = fa*b1[q] + fb*b2[q];
        mv[q] = v.a1(m,k,j,ii)*v.x1f(m,ii)*rh[q];
      }
      HydroAmCross(mv, fv, t);
      for (int q = 0; q < 3; ++q) {tc[q] += n ? -t[q] : t[q];}
    }
    // planar xi / eta faces: (int_f x dA) x F
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
        Real qa[3], qb[3], f1[3], f2[3];
        cubed_sphere::PanelToCart(pnl, xa, ea, qa);
        cubed_sphere::PanelToCart(pnl, xb, eb, qb);
        cubed_sphere::PanelTangents(pnl, xf, ef, f1, f2);
        const Real cf = f1[0]*f2[0] + f1[1]*f2[1] + f1[2]*f2[2];
        const Real sf = 1.0 - cf*cf;
        const Real fa = (g2 - cf*g3)/sf, fb = (g3 - cf*g2)/sf;
        Real cx[3];
        HydroAmCross(qa, qb, cx);
        const Real dot = qa[0]*qb[0] + qa[1]*qb[1] + qa[2]*qb[2];
        const Real th = sqrt(cx[0]*cx[0] + cx[1]*cx[1] + cx[2]*cx[2])/(1.0 + dot);
        Real mv[3], fv[3], t[3];
        for (int q = 0; q < 3; ++q) {
          mv[q] = mom*th*(qa[q] + qb[q]);
          fv[q] = fa*f1[q] + fb*f2[q];
        }
        HydroAmCross(mv, fv, t);
        for (int q = 0; q < 3; ++q) {tc[q] += n ? -t[q] : t[q];}
      }
    }
  } else {
    // spherical polar: orthonormal (e_r, e_theta, e_phi), IM2 = theta, IM3 = phi
    const Real th = v.x2v(m,j), ph = v.x3v(m,k);
    const Real thl = v.x2f(m,j), thr = v.x2f(m,j+1);
    const Real phl = v.x3f(m,k), phr = v.x3f(m,k+1);
    auto basis = [&](const Real t, const Real p, Real *er_, Real *et_, Real *ep_) {
      const Real st = sin(t), ct = cos(t), sp = sin(p), cp = cos(p);
      er_[0] = st*cp; er_[1] = st*sp; er_[2] = ct;
      et_[0] = ct*cp; et_[1] = ct*sp; et_[2] = -st;
      ep_[0] = -sp;   ep_[1] = cp;    ep_[2] = 0.0;
    };
    basis(th, ph, rh, b1, b2);
    for (int q = 0; q < 3; ++q) {rc[q] = r*rh[q];}
    cc = 0.0;
    // T_code: the theta / phi momentum changes of the update and SrcTermsSphericalPolarHydro
    const Real fl1t = v.f1(m,IM2,k,j,i)*v.a1(m,k,j,i);
    const Real fr1t = v.f1(m,IM2,k,j,i+1)*v.a1(m,k,j,i+1);
    const Real fl1p = v.f1(m,IM3,k,j,i)*v.a1(m,k,j,i);
    const Real fr1p = v.f1(m,IM3,k,j,i+1)*v.a1(m,k,j,i+1);
    Real d2 = -(fr1t - fl1t)*ia
              - (v.f2(m,IM2,k,j+1,i)*v.a2(m,k,j+1,i) - v.f2(m,IM2,k,j,i)*v.a2(m,k,j,i))*ia
              - (v.f3(m,IM2,k+1,j,i)*v.a3(m,k+1,j,i) - v.f3(m,IM2,k,j,i)*v.a3(m,k,j,i))*ia;
    Real d3 = -(fr1p - fl1p)*ia
              - (v.f2(m,IM3,k,j+1,i)*v.a2(m,k,j+1,i) - v.f2(m,IM3,k,j,i)*v.a2(m,k,j,i))*ia
              - (v.f3(m,IM3,k+1,j,i)*v.a3(m,k+1,j,i) - v.f3(m,IM3,k,j,i)*v.a3(m,k,j,i))*ia;
    // SrcTermsSphericalPolar{Hydro,MHD}: m_pp = p + rho v3^2 (+ (b1^2 + b2^2 - b3^2)/2)
    Real mpp = pr + rho*v3*v3;
    if (v.mhd) {
      const Real b1 = v.bcc(m,IBX,k,j,i), b2 = v.bcc(m,IBY,k,j,i), b3 = v.bcc(m,IBZ,k,j,i);
      mpp += 0.5*(b1*b1 + b2*b2 - b3*b3);
    }
    d2 += v.xovr(m,k,j,i)*mpp - fac*(fl1t + fr1t)*ia;
    d3 += -v.yovr(m,k,j,i)*(v.f2(m,IM3,k,j,i)*v.a2(m,k,j,i)
                            + v.f2(m,IM3,k,j+1,i)*v.a2(m,k,j+1,i))*ia
          - fac*(fl1p + fr1p)*ia;
    {
      Real dp[3];
      for (int q = 0; q < 3; ++q) {dp[q] = (d2*b1[q] + d3*b2[q])*dv;}
      HydroAmCross(rc, dp, tk);
    }
    // T_cons
    auto addf = [&](const Real *mv, const Real *fv, const bool hi) {
      Real t[3];
      HydroAmCross(mv, fv, t);
      for (int q = 0; q < 3; ++q) {tc[q] += hi ? -t[q] : t[q];}
    };
    const Real dph = phr - phl;
    const Real s2i = 0.5*(thr - thl) - 0.25*(sin(2.0*thr) - sin(2.0*thl));
    const Real c2i = 0.5*(sin(thr)*sin(thr) - sin(thl)*sin(thl));
    // radial caps: int x dA = r_f^3 int rhat dOmega; F tangential on the cell's basis
    for (int n=0; n<2; ++n) {
      const int ii = i + n;
      const Real rf = v.x1f(m,ii), r3 = rf*rf*rf;
      const Real mv[3] = {r3*s2i*(sin(phr) - sin(phl)), r3*s2i*(cos(phl) - cos(phr)),
                          r3*c2i*dph};
      Real fv[3];
      for (int q = 0; q < 3; ++q) {
        fv[q] = v.f1(m,IM2,k,j,ii)*b1[q] + v.f1(m,IM3,k,j,ii)*b2[q];
      }
      addf(mv, fv, n == 1);
    }
    // theta cones: int x dA = sin(t_f) (r^3/3) int (sin t_f cos p, sin t_f sin p, cos t_f) dp
    for (int n=0; n<2; ++n) {
      const int jj = j + n;
      const Real tf = v.x2f(m,jj), st = sin(tf), ct = cos(tf);
      const Real mv[3] = {st*mom*st*(sin(phr) - sin(phl)), st*mom*st*(cos(phl) - cos(phr)),
                          st*mom*ct*dph};
      Real er_[3], et_[3], ep_[3], fv[3];
      basis(tf, ph, er_, et_, ep_);
      for (int q = 0; q < 3; ++q) {
        fv[q] = v.f2(m,IM2,k,jj,i)*et_[q] + v.f2(m,IM3,k,jj,i)*ep_[q];
      }
      addf(mv, fv, n == 1);
    }
    // phi half-planes: int x dA = (r^3/3) int (sin t cos p_f, sin t sin p_f, cos t) dt
    for (int n=0; n<2; ++n) {
      const int kk = k + n;
      const Real pf = v.x3f(m,kk);
      const Real dc = cos(thl) - cos(thr), ds = sin(thr) - sin(thl);
      const Real mv[3] = {mom*dc*cos(pf), mom*dc*sin(pf), mom*ds};
      Real er_[3], et_[3], ep_[3], fv[3];
      basis(th, pf, er_, et_, ep_);
      for (int q = 0; q < 3; ++q) {
        fv[q] = v.f3(m,IM2,kk,j,i)*et_[q] + v.f3(m,IM3,kk,j,i)*ep_[q];
      }
      addf(mv, fv, n == 1);
    }
  }
  for (int q = 0; q < 3; ++q) {dT[q] = tc[q] - tk[q];}
  return true;
}


//! \brief cs-open-1018 <hydro>/frame_omega: the rotating-frame force density of cell
//! (m,k,j,i), f = -2 rho Omega x v - rho Omega x (Omega x rc) (Cartesian, cell centre), and
//! (amc) the L-conserving torque correction dT = T* - V rc x f, with
//!   T* = -Q_c - Omega x L_abs,c,  L_abs,c = V rc x rho (v + Omega x rc),
//!   Q_c = sum_f Phi_f^out (w_f - w_c),  w(x) = x x (Omega x x),
//! Phi^out the outward mass flux (x area) of the update and w_f at the face centre (the
//! generalisation of the dhj problem/coriolis_am from z to a vector Omega): the frame terms
//! then change sum_c L_abs by -Omega x L_abs plus the boundary flux of w exactly (the
//! component along Omega is conserved to round-off; |L_abs| to the time integration).
//! rc, the cell basis (b1, b2) and cc as HydroAmVecTorque; v from w0 (cs: contravariant).

KOKKOS_INLINE_FUNCTION
void HydroAmFrame(const HydroAmV &v, const Real *om, const bool amc, const int m,
                  const int k, const int j, const int i, Real *f, Real *dT, Real *rc,
                  Real *rh, Real *b1, Real *b2, Real &cc) {
  const Real dv = v.vol(m,k,j,i);
  const Real r = v.x1v(m,i);
  const Real rho = v.w0(m,IDN,k,j,i);
  const Real v1 = v.w0(m,IVX,k,j,i), v2 = v.w0(m,IVY,k,j,i), v3 = v.w0(m,IVZ,k,j,i);
  auto wv = [&](const Real *x, Real *w) {   // x x (Omega x x)
    Real ox[3];
    HydroAmCross(om, x, ox);
    HydroAmCross(x, ox, w);
  };
  Real xf[6][3];       // the six face centres
  Real phi[6];         // outward mass flux x area
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
    Real q[3];
    cubed_sphere::PanelToCart(pnl, xc, ec, q);
    cubed_sphere::PanelTangents(pnl, xc, ec, b1, b2);
    const Real qn = 1.0/sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
    for (int c = 0; c < 3; ++c) {rh[c] = q[c]*qn;}
    cc = v.ccell(m,k,j);
    const Real xs[4] = {xl, xr, xc, xc}, es[4] = {ec, ec, el, er};
    for (int n = 0; n < 4; ++n) {
      cubed_sphere::PanelToCart(pnl, xs[n], es[n], q);
      const Real nq = r/sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
      for (int c = 0; c < 3; ++c) {xf[2+n][c] = q[c]*nq;}
    }
  } else {
    const Real th = v.x2v(m,j), ph = v.x3v(m,k);
    const Real st = sin(th), ct = cos(th), sp = sin(ph), cp = cos(ph);
    rh[0] = st*cp; rh[1] = st*sp; rh[2] = ct;
    b1[0] = ct*cp; b1[1] = ct*sp; b1[2] = -st;
    b2[0] = -sp;   b2[1] = cp;    b2[2] = 0.0;
    cc = 0.0;
    const Real ts[4] = {v.x2f(m,j), v.x2f(m,j+1), th, th};
    const Real ps[4] = {ph, ph, v.x3f(m,k), v.x3f(m,k+1)};
    for (int n = 0; n < 4; ++n) {
      xf[2+n][0] = r*sin(ts[n])*cos(ps[n]);
      xf[2+n][1] = r*sin(ts[n])*sin(ps[n]);
      xf[2+n][2] = r*cos(ts[n]);
    }
  }
  for (int c = 0; c < 3; ++c) {
    rc[c] = r*rh[c];
    xf[0][c] = v.x1f(m,i)*rh[c];
    xf[1][c] = v.x1f(m,i+1)*rh[c];
  }
  Real vel[3];
  for (int c = 0; c < 3; ++c) {vel[c] = v1*rh[c] + v2*b1[c] + v3*b2[c];}
  Real ov[3], orc[3], oorc[3];
  HydroAmCross(om, vel, ov);
  HydroAmCross(om, rc, orc);
  HydroAmCross(om, orc, oorc);
  for (int c = 0; c < 3; ++c) {f[c] = -2.0*rho*ov[c] - rho*oorc[c];}
  for (int c = 0; c < 3; ++c) {dT[c] = 0.0;}
  if (!amc) return;
  phi[0] = -v.f1(m,IDN,k,j,i)*v.a1(m,k,j,i);
  phi[1] = v.f1(m,IDN,k,j,i+1)*v.a1(m,k,j,i+1);
  phi[2] = -v.f2(m,IDN,k,j,i)*v.a2(m,k,j,i);
  phi[3] = v.f2(m,IDN,k,j+1,i)*v.a2(m,k,j+1,i);
  phi[4] = -v.f3(m,IDN,k,j,i)*v.a3(m,k,j,i);
  phi[5] = v.f3(m,IDN,k+1,j,i)*v.a3(m,k+1,j,i);
  Real wc[3], qv[3] = {0.0, 0.0, 0.0};
  wv(rc, wc);
  for (int n = 0; n < 6; ++n) {
    Real wf[3];
    wv(xf[n], wf);
    for (int c = 0; c < 3; ++c) {qv[c] += phi[n]*(wf[c] - wc[c]);}
  }
  Real pa[3], la[3], ol[3], tcd[3];
  for (int c = 0; c < 3; ++c) {pa[c] = rho*(vel[c] + orc[c])*dv;}
  HydroAmCross(rc, pa, la);
  HydroAmCross(om, la, ol);
  Real fv[3];
  for (int c = 0; c < 3; ++c) {fv[c] = f[c]*dv;}
  HydroAmCross(rc, fv, tcd);
  for (int c = 0; c < 3; ++c) {dT[c] = (-qv[c] - ol[c]) - tcd[c];}
}


//! \brief cs-open-1018: the fluid arrays the angular-momentum module works on (Hydro or MHD)
struct AmFluid {
  DvceArray5D<Real> u0, u1, w0, wder, bcc;
  DvceFaceFld5D<Real> uflx;
  EOS_Data eos;
  bool mhd;
  bool wbs;   // the fluid's wellbalance_static (refused)
};

//! \class AmModule
//! \brief cs-open-1018: <hydro|mhd>/am_transport, am_vector, am_spin_return, frame_omega,
//! frame_exact (all read only when named), owned by Hydro or MHD (hydro_am.cpp)
class AmModule {
 public:
  AmModule(MeshBlockPack *ppack, ParameterInput *pin, const std::string &blk);
  ~AmModule();
  static bool Named(ParameterInput *pin, const std::string &blk);
  void Apply(Driver *pdrive, int stage, Real bdt, const AmFluid &fl);
  MeshBlockPack *pmy_pack;
  std::string blk;
  bool am_transport = false;
  bool am_checked = false;
  bool am_vector = false;
  DvceArray4D<Real> am_spin, am_spin1, am_spins;
  Real am_eta = 0.0;
  MeshBoundaryValuesCC *pbval_spin = nullptr;
  DvceArray5D<Real> spw, spw_c;
  bool frame_on = false;
  Real frame_om[3] = {0.0, 0.0, 0.0};
  bool frame_exact = false;
  DvceArray5D<Real> frame_ln;
  void AmTransport(Real bdt, Driver *pdrive, int stage, const AmFluid &fl);
  void AmSpinReturn(const AmFluid &fl);
  void FrameRotation(Real bdt, Driver *pdrive, int stage, const AmFluid &fl);
  void AmFillV(HydroAmV &v, const AmFluid &fl);
  void BasisCheck();
};

//! \brief the module of the pack's fluid (Hydro first, else MHD), nullptr if none
AmModule *AmModuleOf(MeshBlockPack *pp);
// cs-open-1018: the restart block of <hydro|mhd>/am_vector (outputs/restart.cpp, read in
// pgen.cpp): header {1, nslab, 0, 0}; slab 0 = am_spin, slab 1 (if present) = the
// radiation's hesdirk2 FSAL spin slope amk1 (RadiationM1)
constexpr char kAmSpinRstMagic[8] = {'A', 'M', 'S', 'P', 'I', 'N', '0', '1'};

}  // namespace hydro
#endif  // HYDRO_HYDRO_AM_HPP_
