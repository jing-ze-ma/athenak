#ifndef RAD_VET_RAD_VET_FACES_HPP_
#define RAD_VET_RAD_VET_FACES_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_faces.hpp
//! \brief device functions on faces: metric, seams, weights, flux, stress, sources.
//!
//! Data layout and device functions of the clean implicit radiation solve.
//! Each stage solves for E per cell, with face fluxes in terms of E:
//!   F_f = (1-w) [theta (G_f - h c^2 grad(D E)) + (1-theta) A_f]
//!         + w c (h+_f E_L + h-_f E_R)
//! theta = 1/(1 + h c chi_f); A = a E, a_d = v_d + D_db v_b; h+-_f the rays' half-range
//! fluxes over the donor's J; w a smooth function of C1 = <|I - S|>/J.  Cell equation:
//!   E (1 + h c kE) + h div F = E_old + h c kP a T^4.

#include "rad_vet/rad_vet_args.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn Real CleanDuu
//! \brief u.D.u of the Eddington tensor along u = g(m, gs..gs+2, k,j,i), for
//! <hydro>/rad_signal_speed

KOKKOS_INLINE_FUNCTION
Real CleanDuu(const DvceArray5D<Real> &s, const CleanGeo &g, const int gs,
              const int m, const int k, const int j, const int i) {
  Real d = 0.0;
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      d += g(m,gs+a,k,j,i)*g(m,gs+b,k,j,i)*s(m,CLN_D11+((a == b) ? a : 2 + a + b),k,j,i);
    }
  }
  return d;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanSmoothStep
//! \brief 3x^2 - 2x^3 on [0,1], 0 below and 1 above

KOKKOS_INLINE_FUNCTION
Real CleanSmoothStep(Real x) {
  x = fmin(fmax(x, 0.0), 1.0);
  return x*x*(3.0 - 2.0*x);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanExchangeTemperature
//! \brief gas T after an implicit local exchange: the root of e(d,T) + A T^4 = rhs
//! (A = h c kP a/(1 + h c kE)); unique, below T(d, rhs); safeguarded Newton.
//! rtol > 0 (clean_troot_rtol): stop at |g| <= rtol (e + A T^4 + rhs), the noise level of
//! a tabulated e(T); the bound T(d, rhs) is then computed only when a step leaves it

KOKKOS_INLINE_FUNCTION
Real CleanExchangeTemperature(const EOS_Data &eos, const Real d, const Real a,
                              const Real rhs, const Real tguess, const Real rtol = 0.0) {
  if (!(rhs > 0.0)) return 0.0;
  Real lo = 0.0, hi = (rtol > 0.0 && tguess > 0.0) ? -1.0 : eos.Temperature(d, rhs);
  Real t = (tguess > 0.0 && (hi < 0.0 || tguess < hi)) ? tguess : 0.5*hi;
  for (int n = 0; n < 100; ++n) {
    Real e, p, chi_rho, chi_t, cv;
    eos.ThermoAt(d, t, e, p, chi_rho, chi_t, cv);
    const Real t3 = t*t*t;
    const Real g = e + a*t3*t - rhs;
    if (rtol > 0.0 && fabs(g) <= rtol*(fabs(e) + a*t3*t + rhs)) {
      return t - g/(d*cv + 4.0*a*t3);
    }
    if (g > 0.0) {
      hi = t;
    } else {
      lo = t;
    }
    Real tn = t - g/(d*cv + 4.0*a*t3);
    if (!(tn > lo && (hi < 0.0 || tn < hi))) {
      if (hi < 0.0) hi = eos.Temperature(d, rhs);
      tn = 0.5*(lo + hi);
    }
    if (fabs(tn - t) <= 1.0e-15*t) return tn;
    t = tn;
  }
  return t;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanStagSph
//! \brief spherical-polar metric of face (k,j,i) normal to dir: q = 0 inverse centre
//! distance, q = 1 inverse radius, q = 2, 3 curvature coefficient h of cell L, R in
//! d_n(D_nn E) + transverse + (h_L E_L + h_R E_R)/(2 r):
//!   r faces:     h = 2 D_rr - D_thth - D_phph + cot(theta) D_rth
//!   theta faces: h = cot(theta_f) (D_thth - D_phph) + 3 D_rth
//!   phi faces:   none (torque form, so L_z changes only by boundary torques)
//! q = 4: s = 1/r_f (r), cot(theta_f)/(2 r) (theta), so P enters as (1/dx +- s) P_R,L,
//! the discrete (1/A) d_n(A P).  D = I/3 uniform is exactly steady

KOKKOS_INLINE_FUNCTION
Real CleanStagSph(const CleanArgs &x, const int q, const int dir, const int m,
                  const int k, const int j, const int i) {
  if (q == 0) {
    return (dir == 2) ? x.area[2](m,k,j,i)/x.vol(m,k,j,i) : 1.0/x.dxf[dir](m,k,j,i);
  }
  if (dir == 2 && q >= 2) return 0.0;   // (phi faces: torque form, CleanFaceStress)
  // (a theta-periodic seam face takes one cot in both copies: single-valued flux)
  const BoundaryFlag pf = BoundaryFlag::periodic;
  const bool sm = (dir == 1) &&
                  ((j == x.js && x.mbc(m,BoundaryFace::inner_x2) == pf) ||
                   (j == x.je + 1 && x.mbc(m,BoundaryFace::outer_x2) == pf));
  if (q == 4) {
    return (dir == 0) ? 1.0/x.x1f(m,i) :
           (sm ? 0.5*x.cseam/x.x1v(m,i) : 0.5/(x.tg(m,1,j)*x.x1v(m,i)));
  }
  if (q == 1) return 1.0/((dir == 0) ? x.x1f(m,i) : x.x1v(m,i));
  const int kc = (q == 2) ? k - (dir == 2) : k, jc = (q == 2) ? j - (dir == 1) : j;
  const int ic = (q == 2) ? i - (dir == 0) : i;
  auto d = [&](const int n) {return x.cs(m,CLN_D11+n,kc,jc,ic);};
  if (dir == 0) {
    return 2.0*d(0) - d(1) - d(2) + (x.multi_d ? d(3)/x.tg(m,0,jc) : 0.0);
  }
  if (sm) return (d(1) - d(2))*x.cseam + 3.0*d(3);
  return (d(1) - d(2))/x.tg(m,1,j) + 3.0*d(3);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanDivW
//! \brief weight of the low (side 0) or high (side 1) face flux in the divergence of
//! cell (k,j,i): 1/dx_d; A_f/V on curved meshes, on the cubed sphere with the hydro's
//! curved radial area (not the closed vector areas of CleanCsGeometry)

KOKKOS_INLINE_FUNCTION
Real CleanDivW(const CleanArgs &x, const int dir, const int m, const int k, const int j,
               const int i, const int side) {
  if (!x.sph && !x.cub) return x.idx[dir];
  const auto &a = (x.cub && dir == 0) ? x.ae0 : x.area[dir];
  return a(m,k+side*(dir == 2),j+side*(dir == 1),i+side*(dir == 0))/x.vol(m,k,j,i);
}

//! \brief divergence of face values fl, fh along dir (CleanDivW); CleanVol: cell volume

KOKKOS_INLINE_FUNCTION
Real CleanDivF(const CleanArgs &x, const int dir, const int m, const int k, const int j,
               const int i, const Real fl, const Real fh) {
  if (!x.sph && !x.cub) return x.idx[dir]*(fh - fl);
  return CleanDivW(x, dir, m, k, j, i, 1)*fh - CleanDivW(x, dir, m, k, j, i, 0)*fl;
}

KOKKOS_INLINE_FUNCTION
Real CleanVol(const CleanArgs &x, const int m, const int k, const int j, const int i) {
  return (x.sph || x.cub) ? x.vol(m,k,j,i) : 1.0/(x.idx[0]*x.idx[1]*x.idx[2]);
}

//----------------------------------------------------------------------------------------
//! \fn bool CleanCsSeamFace
//! \brief cubed sphere: true for a panel seam face normal to dir, not next to a vertex

KOKKOS_INLINE_FUNCTION
bool CleanCsSeamFace(const CleanArgs &x, const int dir, const int m, const int k,
                     const int j, const int i) {
  if (!x.cub || dir == 0) return false;
  const BoundaryFlag pf = BoundaryFlag::panel;
  const int b = 3 - dir, p = (dir == 1) ? j : k, o = (dir == 1) ? k : j;
  const int pl = (dir == 1) ? x.js : x.ks, ph = ((dir == 1) ? x.je : x.ke) + 1;
  const int ol = (dir == 1) ? x.ks : x.js, oh = (dir == 1) ? x.ke : x.je;
  return ((p == pl && x.mbc(m,2*dir) == pf) || (p == ph && x.mbc(m,2*dir+1) == pf)) &&
         !((o == ol && x.mbc(m,2*b) == pf) || (o == oh && x.mbc(m,2*b+1) == pf));
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanCsNormA
//! \brief cubed sphere: area A_n in the face momentum's A_n (P_R - P_L)/V_f: the face
//! area, but on a seam face (mirror-image cells, whose volume variation adds in V_f)
//! A_n = V_f/d, d = 2 |X| |r_hat.n_f| the chord of the mirror pair's centres
//! clean_cs_lat_wb: a radial face takes the hydro's r^2 Omega (as V_f), not its vector
//! area r^2 |S| (|S|/Omega - 1 = O(dtheta^2) varies over the panel)

KOKKOS_INLINE_FUNCTION
Real CleanCsNormA(const CleanArgs &x, const int dir, const int m, const int k,
                  const int j, const int i) {
  if (x.lwb && dir == 0) return x.ae0(m,k,j,i);
  if (!CleanCsSeamFace(x, dir, m, k, j, i)) return x.area[dir](m,k,j,i);
  const int kl = k - (dir == 2), jl = j - (dir == 1);
  const bool lo = ((dir == 1) ? j : k) == ((dir == 1) ? x.js : x.ks);
  const int ko = lo ? k : kl, jo = lo ? j : jl;
  Real xr = 0.0, rn = 0.0;
  for (int c = 0; c < 3; ++c) {
    xr += SQR(x.g.At<CLG_X>(m,c,ko,jo,i));
    rn += x.g.At<CLG_R>(m,c,ko,jo,i)*x.g.At<CLG_N>(m,3*dir+c,k,j,i);
  }
  return 0.5*(x.g(m,CLG_F+3,k,j,i) + x.g(m,CLG_F+3,kl,jl,i))/(2.0*std::sqrt(xr)*fabs(rn));
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanStagDnn
//! \brief n.D.n of cell (k,j,i), n the normal of face (kf,jf,iff): D_dir,dir off the
//! cubed sphere

KOKKOS_INLINE_FUNCTION
Real CleanStagDnn(const CleanArgs &x, const int dir, const int m, const int kf,
                  const int jf, const int iff, const int k, const int j, const int i) {
  if (!x.cub) return x.cs(m,CLN_D11+dir,k,j,i);
  Real d = 0.0;
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      d += x.g.At<CLG_N>(m,3*dir+a,kf,jf,iff)*x.g.At<CLG_N>(m,3*dir+b,kf,jf,iff)*
           x.cs(m,CLN_D11+((a == b) ? a : 2 + a + b),k,j,i);
    }
  }
  return d;
}

//----------------------------------------------------------------------------------------
//! \struct CleanFaceMet
//! \brief metric of a face momentum row: a (P_R - P_L)/v + s (P_R + P_L)
//! + hr (h_L P_L + h_R P_R), h of CleanStagSph:
//!   Cartesian:       a = 1/dx, v = 1, s = hr = 0;
//!   spherical polar: a = CleanStagSph q = 0, v = 1, s its q = 4, hr half its q = 1;
//!   cubed sphere:    a = A_n (CleanCsNormA), v = V_f (owners' half cells), s = hr = 0

struct CleanFaceMet {Real a, v, s, hr;};

KOKKOS_INLINE_FUNCTION
CleanFaceMet CleanFaceMetric(const CleanArgs &x, const int dir, const int m, const int k,
                             const int j, const int i) {
  const int kl = k - (dir == 2), jl = j - (dir == 1);
  if (x.cub) {
    return {CleanCsNormA(x, dir, m, k, j, i),
            0.5*(x.g(m,CLG_F+3,k,j,i) + x.g(m,CLG_F+3,kl,jl,i - (dir == 0))), 0.0, 0.0};
  }
  if (!x.sph) return {x.idx[dir], 1.0, 0.0, 0.0};
  return {CleanStagSph(x, 0, dir, m, k, j, i), 1.0, CleanStagSph(x, 4, dir, m, k, j, i),
          0.5*CleanStagSph(x, 1, dir, m, k, j, i)};
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanStagMemUp
//! \brief clean_stag_gflux: the flux memory (1 - w) theta (G + d) taken implicitly upwind
//! as q E_donor, q = memory/E_old_donor, |q| <= c, where its explicit outflow h A/V |mem|
//! would exceed CLN_MEMUP E_old_donor.  Returns q (> 0: from L, < 0: from R), else 0

KOKKOS_INLINE_FUNCTION
Real CleanStagMemUp(const CleanArgs &x, const int dir, const int m, const int k,
                    const int j, const int i) {
  if (!x.gf || x.hi) return 0.0;
  const int kl = k - (dir == 2), jl = j - (dir == 1), il = i - (dir == 0);
  // (clean_stag_vimp: G takes theta (1 - b), CleanFaceSource)
  const Real qd = x.drag ? x.fc[dir](m,CLF_DB,k,j,i) : 0.0;
  const Real g = (1.0 - x.fc[dir](m,CLF_W,k,j,i))*x.fc[dir](m,CLF_TH,k,j,i)*
                 ((1.0 - qd)*x.fc[dir](m,CLF_G,k,j,i) +
                  (x.fd ? x.fc[dir](m,CLF_N+3,k,j,i) : 0.0));
  if (g == 0.0) return 0.0;
  const bool up = (g > 0.0);   // (the donor: L for an outflow to R)
  const Real eo = up ? x.w(m,CLW_EO,kl,jl,il) : x.w(m,CLW_EO,k,j,i);
  // (cubed sphere: A_f over the owner's V, so both sides of a seam decide alike)
  const Real out = x.h*fabs(g)*(!x.cub ? (up ? CleanDivW(x, dir, m, kl, jl, il, 1) :
                                               CleanDivW(x, dir, m, k, j, i, 0)) :
                   ((dir == 0) ? x.ae0 : x.area[dir])(m,k,j,i)/
                   (up ? x.g(m,CLG_F+3,kl,jl,il) : x.g(m,CLG_F+3,k,j,i)));
  if (out <= CLN_MEMUP*eo) return 0.0;
  const Real q = (eo > 0.0) ? fmin(fabs(g)/eo, x.c) : x.c;
  return up ? q : -q;
}

//----------------------------------------------------------------------------------------
//! \fn int CleanCsSlot
//! \brief cubed sphere: face (kf,jf,iff) of cell (k,j,i) (far: the other face of its
//! family) in the owner's axes: returns 2 q + hi (family q, hi = owner's upper face),
//! s = +-1 the owner's normal along ours (a seam ghost's axes may be swapped/flipped)

KOKKOS_INLINE_FUNCTION
int CleanCsSlotQ(const DvceArray5D<Real> &g, const int dir, const int m, const int kf,
                 const int jf, const int iff, const int k, const int j, const int i,
                 Real &s) {
  Real t2 = 0.0, t3 = 0.0;
  for (int c = 0; c < 3; ++c) {
    t2 += g(m,CLG_F+c,k,j,i)*g(m,CLG_N+3*dir+c,kf,jf,iff);
    t3 += g(m,CLG_F+4+c,k,j,i)*g(m,CLG_N+3*dir+c,kf,jf,iff);
  }
  const int q = (fabs(t2) > fabs(t3)) ? 1 : 2;
  s = (((q == 1) ? t2 : t3) > 0.0) ? 1.0 : -1.0;
  return q;
}

KOKKOS_INLINE_FUNCTION
int CleanCsSlot(const CleanArgs &x, const int dir, const int m, const int kf,
                const int jf, const int iff, const int k, const int j, const int i,
                const bool far, Real &s) {
  s = 1.0;
  const bool up = (dir == 0) ? (iff > i) : ((dir == 1) ? (jf > j) : (kf > k));
  if (dir == 0) return (up != far) ? 1 : 0;
  // (CleanCsSlotQ of the face's cell on this side, CLG_Q)
  const int b = static_cast<int>(x.g(m,CLG_Q,kf,jf,iff)) >> (4*(dir - 1) + 2*up);
  const int q = 1 + (b & 1);
  s = (b & 2) ? 1.0 : -1.0;
  return 2*q + (((up == (s > 0.0)) != far) ? 1 : 0);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanStagCsLat
//! \brief cubed sphere: n.(Y_p + Y_p') of cell (k,j,i), n the face's normal, p, p' the
//! other two face families (owner's axes).  clean_cs_lat_wb, a lateral face: n less its
//! part along the cell's radial (centroid) direction, unit: the radial sum of a
//! spherically symmetric P gives the face no force (the equiangular cells are not mirror
//! images about it); O(dtheta^2) for a smooth P

KOKKOS_INLINE_FUNCTION
Real CleanStagCsLat(const CleanArgs &x, const int dir, const int m, const int kf,
                    const int jf, const int iff, const int k, const int j, const int i) {
  Real s, d = 0.0, n[3];   // (off the cubed sphere, clean_stag_lat_sums: q = dir, e_q)
  const int q = x.cub ? CleanCsSlot(x, dir, m, kf, jf, iff, k, j, i, false, s)/2 : dir;
  const int p = (q == 0) ? 1 : 0, pp = (q == 2) ? 1 : 2;
  for (int c = 0; c < 3; ++c) n[c] = x.cub ? x.g.At<CLG_N>(m,3*dir+c,kf,jf,iff)
                                           : (c == q);
  if (x.lwb && dir > 0) {
    Real a[3], aa = 0.0, an = 0.0, nn = 0.0;
    for (int c = 0; c < 3; ++c) {
      a[c] = x.g.At<CLG_X>(m,c,k,j,i);
      aa += a[c]*a[c];
      an += a[c]*n[c];
    }
    for (int c = 0; c < 3; ++c) {
      n[c] -= an/aa*a[c];
      nn += n[c]*n[c];
    }
    for (int c = 0; c < 3; ++c) n[c] /= sqrt(nn);
  }
  for (int c = 0; c < 3; ++c) {
    d += n[c]*(x.cub ? x.w(m,CLW_L+3*p+c,k,j,i) + x.w(m,CLW_L+3*pp+c,k,j,i)
                     : x.yl(m,3*p+c,k,j,i) + x.yl(m,3*pp+c,k,j,i));
  }
  return d;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanStagN
//! \brief normal component n.v of work vector v..v+2 of cell (k,j,i) at face (kf,jf,iff)

KOKKOS_INLINE_FUNCTION
Real CleanStagN(const CleanArgs &x, const int v, const int dir, const int m,
                const int kf, const int jf, const int iff, const int k, const int j,
                const int i) {
  if (!x.cub) return x.w(m,v+dir,k,j,i);
  Real d = 0.0;
  for (int b = 0; b < 3; ++b) d += x.g.At<CLG_N>(m,3*dir+b,kf,jf,iff)*x.w(m,v+b,k,j,i);
  return d;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStagArea
//! \brief the net streaming part s = h+ + h- in area form: adds w c s (A_d/A_f - 1) to
//! the donor's coefficient pl or pr (A_d the area at the donor's centre), so a free
//! stream keeps A_d E_d; signs (M-matrix) kept.  CleanAreaRatio: A_d/A_f

KOKKOS_INLINE_FUNCTION
Real CleanAreaRatio(const CleanArgs &x, const int dir, const int m, const int k,
                    const int j, const int i, const int kd, const int jd, const int id) {
  if (x.cub) {
    Real sd;
    return x.g.At<CLG_A>(m,CleanCsSlot(x, dir, m, k, j, i, kd, jd, id, false, sd)/2,
                         kd,jd,id)/x.area[dir](m,k,j,i);
  }
  return (dir == 0) ? x.x1v(m,id)*x.x1v(m,id)/(x.x1f(m,i)*x.x1f(m,i))
                    : x.tg(m,2,jd)/x.tg(m,3,j);
}

KOKKOS_INLINE_FUNCTION
void CleanStagArea(const CleanArgs &x, const int dir, const int m, const int k,
                   const int j, const int i, const Real w, Real &pl, Real &pr) {
  if (!x.sar || x.hi || !(x.sph || x.cub) || (!x.cub && dir == 2)) return;
  const Real s = w*x.c*(x.fc[dir](m,CLF_HP,k,j,i) + x.fc[dir](m,CLF_HM,k,j,i));
  const int d = (s > 0.0) ? 1 : 0, id = i - d*(dir == 0), jd = j - d*(dir == 1);
  const int kd = k - d*(dir == 2);
  ((d == 1) ? pl : pr) += s*(CleanAreaRatio(x, dir, m, k, j, i, kd, jd, id) - 1.0);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFaceChi
//! \brief chi of a face: mean of its two cells; an x1 face off physical boundaries takes
//! the distance-weighted harmonic mean

KOKKOS_INLINE_FUNCTION
Real CleanFaceChi(const DvceArray5D<Real> &op, const DvceArray2D<Real> &x1f,
                  const DvceArray2D<BoundaryFlag> &bc, const int dir,
                  const int m, const int k, const int j, const int i, const int is,
                  const int ie) {
  const Real cl = op(m,VET_OP_T,k-(dir == 2),j-(dir == 1),i-(dir == 0));
  const Real cr = op(m,VET_OP_T,k,j,i);
  auto phys = [&](const int s) {
    return bc(m,s) != BoundaryFlag::block && bc(m,s) != BoundaryFlag::periodic;
  };
  if (dir != 0 || (i == is && phys(0)) || (i == ie + 1 && phys(1)) ||
      !(cl > 0.0 && cr > 0.0)) {
    return 0.5*(cl + cr);
  }
  if (x1f.extent_int(1) <= i + 1) return 2.0*cl*cr/(cl + cr);   // (Cartesian: uniform)
  const Real hl = x1f(m,i) - x1f(m,i-1), hr = x1f(m,i+1) - x1f(m,i);
  return (hl + hr)/(hl/cl + hr/cr);
}

//----------------------------------------------------------------------------------------
//! \fn void CleanFaceWeights
//! \brief coefficients al, ar of E_L, E_R in the face flux (all but the cross terms and
//! G); physical boundaries take their boundary form; up: clean_stag_upwind's share phi

KOKKOS_INLINE_FUNCTION
bool CleanBoundaryFace(const CleanArgs &x, const int dir, const int m, const int k,
                       const int j, const int i);

KOKKOS_INLINE_FUNCTION
void CleanFaceWeights(const CleanArgs &x, const int dir, const int m, const int k,
                      const int j, const int i, Real &al, Real &ar, Real *up = nullptr) {
  if (x.fwc) {
    al = x.fc[dir](m,CLF_QL,k,j,i);
    ar = x.fc[dir](m,CLF_QL+1,k,j,i);
    return;
  }
  const int kl = k - (dir == 2), jl = j - (dir == 1), il = i - (dir == 0);
  const int pos = (dir == 0) ? i : ((dir == 1) ? j : k);
  const int lo = (dir == 0) ? x.is : ((dir == 1) ? x.js : x.ks);
  const int hi = (dir == 0) ? x.ie + 1 : ((dir == 1) ? x.je + 1 : x.ke + 1);
  al = 0.0;
  ar = 0.0;
  // (an open r end's net outflow h+ + h- in area form, as CleanStagArea)
  const bool ba = dir == 0 && (x.sph || x.cub);
  if (pos == lo && CleanBc(x, m, 2*dir) != VET_IBC_PERIODIC) {
    if (!CleanOpen(x, m, 2*dir)) return;   // (a wall; an imposed flux is all source)
    const bool beam = x.use_beam && dir == 0 && x.beam(m,k-x.ks,j-x.js) > 0.0;
    ar = x.c*x.cs(m,CLN_HM+dir,k,j,i);
    if (ba) ar += x.c*fmin(x.cs(m,CLN_HP,k,j,i) + x.cs(m,CLN_HM,k,j,i), 0.0)*
                  (CleanAreaRatio(x, 0, m, k, j, i, k, j, i) - 1.0);
    // (clean_stag_upwind: an open end carries a E, the end cell's on outflow)
    if (x.upw && x.feedback && !x.hi && CleanOpen(x, m, 2*dir) && !beam) {
      ar += fmin(CleanStagN(x, CLW_A1, dir, m, k, j, i, k, j, i), 0.0);
    }
    return;
  }
  if (pos == hi && CleanBc(x, m, 2*dir+1) != VET_IBC_PERIODIC) {
    if (!CleanOpen(x, m, 2*dir+1)) return;
    al = x.c*x.cs(m,CLN_HP+dir,kl,jl,il);
    if (ba) al += x.c*fmax(x.cs(m,CLN_HP,kl,jl,il) + x.cs(m,CLN_HM,kl,jl,il), 0.0)*
                  (CleanAreaRatio(x, 0, m, k, j, i, kl, jl, il) - 1.0);
    if (x.upw && x.feedback && !x.hi && CleanOpen(x, m, 2*dir+1)) {
      al += fmax(CleanStagN(x, CLW_A1, dir, m, k, j, i, kl, jl, il), 0.0);
    }
    return;
  }
  const Real w = x.hi ? 0.0 : x.fc[dir](m,CLF_W,k,j,i);
  // (clean_stag_vimp: stress theta_d = theta + b (1 - theta), advective 1 - theta_d,
  // beta term theta (1 - b))
  const Real t0 = x.fc[dir](m,CLF_TH,k,j,i);
  const Real qd = x.drag ? x.fc[dir](m,CLF_DB,k,j,i) : 0.0, th = t0 + qd*(1.0 - t0);
  // (the face-normal momentum's stress a n.D.n E/v, curved: and its curvature terms)
  const CleanFaceMet fm = CleanFaceMetric(x, dir, m, k, j, i);
  const Real dl = CleanStagDnn(x, dir, m, k, j, i, kl, jl, il);
  const Real dr = CleanStagDnn(x, dir, m, k, j, i, k, j, i);
  const Real gs = (1.0 - w)*th*x.h*x.c*x.c;
  if (x.sph) {
    al = gs*(fm.a*dl - fm.hr*CleanStagSph(x, 2, dir, m, k, j, i));
    ar = -gs*(fm.a*dr + fm.hr*CleanStagSph(x, 3, dir, m, k, j, i));
  } else {
    al = gs*fm.a/fm.v*dl;
    ar = -(gs*fm.a/fm.v)*dr;
  }
  al += w*x.c*x.fc[dir](m,CLF_HP,k,j,i);
  ar += w*x.c*x.fc[dir](m,CLF_HM,k,j,i);
  CleanStagArea(x, dir, m, k, j, i, w, al, ar);
  const Real q = CleanStagMemUp(x, dir, m, k, j, i);
  if (q > 0.0) al += q;
  if (q < 0.0) ar += q;
  if (x.feedback) {
    const Real s = 0.5*(1.0 - w)*(1.0 - th);
    const Real ul = CleanStagN(x, CLW_A1, dir, m, k, j, i, kl, jl, il);
    const Real ur = CleanStagN(x, CLW_A1, dir, m, k, j, i, k, j, i);
    const Real sl = s*ul, sr = s*ur;
    // the beta term th c beta_f (hc kP a T^4 - hc kE E)_f of eq. 57 (Ma et al. 2025)
    const Real b = -0.25*(1.0 - w)*t0*(1.0 - qd)*x.c*
                   (CleanStagN(x, CLW_BT, dir, m, k, j, i, kl, jl, il) +
                    CleanStagN(x, CLW_BT, dir, m, k, j, i, k, j, i))*x.h*x.c;
    // of (1 - r) kE E: emission follows E through the gas T (r = s/(d cv + s)); the full
    // hc kE beta would break the M-matrix of an opaque moving cell
    const Real rl = 1.0 - x.w(m,CLW_GX,kl,jl,il);
    const Real rr = 1.0 - x.w(m,CLW_GX,k,j,i);
    const Real bl = b*x.op(m,VET_OP_E,kl,jl,il)*rl, br = b*x.op(m,VET_OP_E,k,j,i)*rr;
    al += sl;
    ar += sr;
    al += bl;
    ar += br;
    // (clean_stag_upwind: where al < 0 or ar > 0, the least donor-cell share phi that
    // restores the M-matrix, continuous; CleanStagAdvectCorr adds phi x the rest, lagged)
    if (x.upw && !(al >= 0.0 && ar <= 0.0)) {
      auto need = [](const Real c, const Real d) {   // least p >= 0 with c + p d >= 0
        return (c >= 0.0) ? 0.0 : ((d > 0.0) ? -c/d : 1.0);
      };
      const Real o = (ul + ur > 0.0) ? 1.0 : -1.0;
      const Real phi = fmin(fmax(need(al, o*sl), need(-ar, o*sr)), 1.0);
      if (up) *up = phi;
      al += o*phi*sl;
      ar -= o*phi*sr;
    }
  }
  // (clean_stag_xdonor: the cross flux q E_donor, CleanStagXDonor; keeps al >= 0 >= ar)
  if (x.xq >= 0 && !x.hi) {   // (the energy flux; F^H keeps the central form)
    al += x.fc[dir](m,x.xq,k,j,i);
    ar += x.fc[dir](m,x.xq+1,k,j,i);
  }
}

//----------------------------------------------------------------------------------------
//! \fn bool CleanBoundaryFace
//! \brief true for a face on a physical (non-periodic) boundary

KOKKOS_INLINE_FUNCTION
bool CleanBoundaryFace(const CleanArgs &x, const int dir, const int m, const int k,
                       const int j, const int i) {
  const int pos = (dir == 0) ? i : ((dir == 1) ? j : k);
  const int lo = (dir == 0) ? x.is : ((dir == 1) ? x.js : x.ks);
  const int hi = (dir == 0) ? x.ie + 1 : ((dir == 1) ? x.je + 1 : x.ke + 1);
  return (pos == lo && CleanBc(x, m, 2*dir) != VET_IBC_PERIODIC) ||
         (pos == hi && CleanBc(x, m, 2*dir+1) != VET_IBC_PERIODIC);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanRayFace
//! \brief the rays' flux through face (k,j,i) normal to dir: the upper cell's (across a
//! MeshBlock edge the neighbour's, one value per face); a panel seam's high face: its last
//! cell's (CLN_FH, as the ghost's x2 or x3 may be another axis)

KOKKOS_INLINE_FUNCTION
Real CleanRayFace(const CleanArgs &x, const int dir, const int m, const int k,
                  const int j, const int i) {
  if (x.cub && dir > 0 && ((dir == 1) ? j > x.je : k > x.ke) &&
      x.mbc(m,2*dir+1) == BoundaryFlag::panel) {
    return x.cs(m,CLN_FH+4*(dir-1),k-(dir == 2),j-(dir == 1),i);
  }
  return x.cs(m,(dir == 0) ? CLN_FX : ((dir == 1) ? CLN_F2 : CLN_F3),k,j,i);
}

//----------------------------------------------------------------------------------------
//! \fn int CleanTransSide
//! \brief 1 if cell (k,j,i)'s transverse stress terms reach across its side s along b
//! (interior, MeshBlock, periodic, reflecting); 0 at open, imposed-flux, pole and seam
//! faces, where the cell itself replaces the ghost (one sided)

KOKKOS_INLINE_FUNCTION
int CleanTransSide(const CleanArgs &x, const int b, const int m, const int k, const int j,
                   const int i, const int s) {
  const int pos = (b == 0) ? i : ((b == 1) ? j : k);
  const int end = s ? ((b == 0) ? x.ie : ((b == 1) ? x.je : x.ke))
                    : ((b == 0) ? x.is : ((b == 1) ? x.js : x.ks));
  const int t = CleanBc(x, m, 2*b + s);
  const BoundaryFlag f = x.mbc(m,2*b+s);
  return (pos != end || (f != BoundaryFlag::polar && f != BoundaryFlag::panel &&
                         (t == VET_IBC_PERIODIC || t == VET_IBC_REFLECT))) ? 1 : 0;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCsSeamCoef
//! \brief cubed sphere: the seam ghost's P - P0 as c0 P0 + c1 P1 + c2 P2 (cells 0, si,
//! 2 si in panel, P = D_n E): second-order extrapolant, d = (4 P1 - 3 P0 - P2)/2 clamped
//! to [0, 2 (P1 - P0)]; the clamp case from E^n keeps the operator linear over the step

KOKKOS_INLINE_FUNCTION
void CleanCsSeamCoef(const CleanArgs &x, const int n, const int q, const int m,
                     const int k, const int j, const int i, const int si, Real c[3]) {
  const int dk = (q == 2)*si, dj = (q == 1)*si, di = (q == 0)*si;
  auto p = [&](const int o) {
    return x.cs(m,n,k+o*dk,j+o*dj,i+o*di)*x.w(m,CLW_EO,k+o*dk,j+o*dj,i+o*di);
  };
  const Real s1 = p(1) - p(0), d = 0.5*(4.0*p(1) - 3.0*p(0) - p(2));
  const Real lo = fmin(0.0, 2.0*s1), hi = fmax(0.0, 2.0*s1);
  const bool cl = (d < lo || d > hi), z = (s1 > 0.0) ? (d < lo) : (d > hi);
  c[0] = !cl ? 2.0 : (z ? -1.0 : 3.0);   // (d: P1 - 2 d; d = 0: P1; d = 2 s1: P1 - 4 s1)
  c[1] = !cl ? -3.0 : (z ? 1.0 : -3.0);
  c[2] = !cl ? 1.0 : 0.0;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFaceStress
//! \brief d_b(D_{dir,b} E) at a face (E = work vector e): the normal part from the two
//! cells (if normal), plus off physical boundaries the transverse part over four cells

CLEAN_HIP_INLINE KOKKOS_INLINE_FUNCTION
Real CleanFaceStress(const CleanArgs &x, const int e, const int dir, const int m,
                     const int k, const int j, const int i, const bool normal) {
  const int kl = k - (dir == 2), jl = j - (dir == 1), il = i - (dir == 0);
  Real div = 0.0;
  const CleanFaceMet fm = CleanFaceMetric(x, dir, m, k, j, i);
  if (normal) {   // (the normal stress of the two (half) cells, CleanFaceMet)
    div = fm.a*(CleanStagDnn(x, dir, m, k, j, i, k, j, i)*x.w(m,e,k,j,i) -
                CleanStagDnn(x, dir, m, k, j, i, kl, jl, il)*x.w(m,e,kl,jl,il));
    if (x.sph) {
      div += fm.hr*(CleanStagSph(x, 2, dir, m, k, j, i)*x.w(m,e,kl,jl,il) +
                    CleanStagSph(x, 3, dir, m, k, j, i)*x.w(m,e,k,j,i));
    }
  }
  if (x.cub) {   // (cubed sphere: lateral sums over V_f; none on a seam face, whose
    // mirror chord is normal to it, so the two-point term is all of n.div(D E))
    if (!CleanBoundaryFace(x, dir, m, k, j, i) && !CleanCsSeamFace(x, dir, m, k, j, i)) {
      div += CleanStagCsLat(x, dir, m, k, j, i, kl, jl, il) +
             CleanStagCsLat(x, dir, m, k, j, i, k, j, i);
    }
    return div/fm.v;
  }
  if (!x.multi_d || CleanBoundaryFace(x, dir, m, k, j, i)) return div;
  // (clean_stag_lat_sums: lateral sums Y over V_f; sph phi faces keep the torque form)
  if (x.ys && !(x.sph && dir == 2)) {
    const Real vf = x.sph ? 0.5*(x.vol(m,k,j,i) + x.vol(m,kl,jl,il)) : 1.0;
    return div + (CleanStagCsLat(x, dir, m, k, j, i, kl, jl, il) +
                  CleanStagCsLat(x, dir, m, k, j, i, k, j, i))/vf;
  }
  for (int b = 0; b < 3; ++b) {
    if (b == dir || (b == 1 && !x.multi_d) || (b == 2 && !x.three_d)) continue;
    const int ab = (dir + b == 1) ? 3 : ((dir + b == 2) ? 4 : 5);   // D12, D13, D23
    const int dk = (b == 2), dj = (b == 1), di = (b == 0);
    auto de = [&](const int kk, const int jj, const int ii) {
      return x.cs(m,ab,kk,jj,ii)*x.w(m,e,kk,jj,ii);
    };
    // (central; one sided where CleanTransSide says the ghosts do not carry the field)
    const int up = CleanTransSide(x, b, m, k, j, i, 1);
    const int dn = CleanTransSide(x, b, m, k, j, i, 0);
    // (torque form: dual r faces i, i+1 and theta faces j, j+1)
    if (x.sph && dir == 2) {
      const Real vr = x.vol(m,k,j,i)*x.x1v(m,i)*x.tg(m,2,j);
      auto pe = [&](const int jj, const int ii) {   // (D E)_{b,ph} at a dual face
        const int jh = (jj > j && !up) ? j : jj, ih = (ii > i && !up) ? i : ii;
        const int jo = jj - (b == 1), io = ii - (b == 0);
        const int jw = (jo < j && !dn) ? j : jo, iw = (io < i && !dn) ? i : io;
        return 0.25*(de(k,jh,ih) + de(k-1,jh,ih) + de(k,jw,iw) + de(k-1,jw,iw));
      };
      if (b == 0) {
        div += (x.area[0](m,k,j,i+1)*x.x1f(m,i+1)*pe(j,i+1) -
                x.area[0](m,k,j,i)*x.x1f(m,i)*pe(j,i))*x.tg(m,2,j)/vr;
      } else {
        div += (x.area[1](m,k,j+1,i)*x.tg(m,3,j+1)*pe(j+1,i) -
                x.area[1](m,k,j,i)*x.tg(m,3,j)*pe(j,i))*x.x1v(m,i)/vr;
      }
      continue;
    }
    if (up + dn == 0) continue;
    const Real rb = !x.sph ? 1.0 : ((dir == 0) ? x.x1f(m,i) : x.x1v(m,i));
    const Real w = (up + dn == 2) ? 1.0 : 2.0;   // (the span: two cells or one)
    const Real ib = !x.sph ? w*x.idx[b] :
                    ((b == 0) ? 2.0/(x.x1v(m,i+up) - x.x1v(m,i-dn)) :
                    ((b == 1) ? 2.0/(rb*(x.x2v(m,j+up) - x.x2v(m,j-dn))) :
                    w*x.idx[2]/(rb*((dir == 1) ? x.tg(m,3,j) : x.tg(m,2,j)))));
    div += 0.25*ib*(de(k+up*dk,j+up*dj,i+up*di) + de(kl+up*dk,jl+up*dj,il+up*di)
                    - de(k-dn*dk,j-dn*dj,i-dn*di) - de(kl-dn*dk,jl-dn*dj,il-dn*di));
  }
  return div;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFaceCross
//! \brief flux of the off-diagonal Eddington terms, -(1-w) theta h c^2 d_b(D_{dir,b} E)
//! for b != dir.  E is work vector e.

KOKKOS_INLINE_FUNCTION
Real CleanFaceCross(const CleanArgs &x, const int e, const int dir, const int m,
                    const int k, const int j, const int i) {
  if (!x.multi_d || (x.xq >= 0 && !x.hi) || CleanBoundaryFace(x, dir, m, k, j, i)) {
    return 0.0;
  }
  const Real w = x.hi ? 0.0 : x.fc[dir](m,CLF_W,k,j,i);
  const Real t0 = x.fc[dir](m,CLF_TH,k,j,i);   // (clean_stag_vimp: theta_d, as
  const Real qd = x.drag ? x.fc[dir](m,CLF_DB,k,j,i) : 0.0;   //  CleanFaceWeights)
  const Real th = t0 + qd*(1.0 - t0);
  return (x.fwc ? x.fc[dir](m,CLF_QL+2,k,j,i) : -(1.0 - w)*th*x.h*x.c*x.c)*
         CleanFaceStress(x, e, dir, m, k, j, i, false);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanStagAdvectCorr
//! \brief clean_stag_upwind: the lagged rest of the upwind share phi of a face,
//!   A_f = a_up E_up + [a_f E_f - a_up E_up]^k,  a_f = (a_L + a_R)/2,
//! donor cell in the matrix (M-matrix at any beta tau), rest at E^k (CLW_EK); E_f the
//! plain central mean

KOKKOS_INLINE_FUNCTION
Real CleanStagAdvectCorr(const CleanArgs &x, const int dir, const int m, const int k,
                         const int j, const int i) {
  const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
  auto e = [&](const int o) {return x.w(m,CLW_EK,k+o*dk,j+o*dj,i+o*di);};
  const Real ul = CleanStagN(x, CLW_A1, dir, m, k, j, i, k-dk, j-dj, i-di);
  const Real ur = CleanStagN(x, CLW_A1, dir, m, k, j, i, k, j, i);
  return 0.5*((ul + ur > 0.0) ? ur*e(0) - ul*e(-1) : ul*e(-1) - ur*e(0));
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanThinFaceDG
//! \brief the lagged thin-face flux w DG, DG = q c (F_y - h+ J_L - h- J_R) of
//! CleanStagBracket (0 without it)

KOKKOS_INLINE_FUNCTION
Real CleanThinFaceDG(const CleanArgs &x, const int dir, const int m, const int k,
                     const int j, const int i, const Real w) {
  return (x.fdg >= 0) ? w*x.fc[dir](m,x.fdg,k,j,i) : 0.0;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFaceSource
//! \brief the E-independent face flux: (1-w) theta G (+ lagged terms) inside; the
//! incoming (beam, bath) or imposed flux on physical boundaries

CLEAN_HIP_INLINE KOKKOS_INLINE_FUNCTION
Real CleanFaceSource(const CleanArgs &x, const int dir, const int m, const int k,
                     const int j, const int i, const Real *phw = nullptr) {
  const int pos = (dir == 0) ? i : ((dir == 1) ? j : k);
  const int lo = (dir == 0) ? x.is : ((dir == 1) ? x.js : x.ks);
  const int hi = (dir == 0) ? x.ie + 1 : ((dir == 1) ? x.je + 1 : x.ke + 1);
  if (pos == lo && CleanBc(x, m, 2*dir) != VET_IBC_PERIODIC) {
    if (CleanBc(x, m, 2*dir) == VET_IBC_REFLECT) return 0.0;
    if (CleanBc(x, m, 2*dir) == VET_IBC_FLUX) return dir ? x.ifx2[0](m,k,i) : x.iflux[0];
    if (x.use_beam && dir == 0 && x.beam(m,k-x.ks,j-x.js) > 0.0) {
      return x.beam(m,k-x.ks,j-x.js);
    }
    const bool up = x.upw && x.feedback && !x.hi;
    const Real ab = up ? fmax(CleanStagN(x, CLW_A1, dir, m, k, j, i, k, j, i), 0.0) : 0.0;
    return 0.25*x.c*x.ebath[2*dir] + ab*x.ebath[2*dir];
  }
  if (pos == hi && CleanBc(x, m, 2*dir+1) != VET_IBC_PERIODIC) {
    if (CleanBc(x, m, 2*dir+1) == VET_IBC_REFLECT) return 0.0;
    if (CleanBc(x, m, 2*dir+1) == VET_IBC_FLUX) return dir ? x.ifx2[1](m,k,i) :
        x.iflux[1];
    const int kl = k - (dir == 2), jl = j - (dir == 1), il = i - (dir == 0);
    const bool up = x.upw && x.feedback && !x.hi;
    const Real ab = up ? fmin(CleanStagN(x, CLW_A1, dir, m, k, j, i, kl, jl, il), 0.0)
                       : 0.0;
    return -0.25*x.c*x.ebath[2*dir+1] + ab*x.ebath[2*dir+1];
  }
  // (clean_stag_xdonor: explicit share of the lagged fluxes, CleanStagXSplit)
  if (x.xq >= 0 && !x.hi) return x.fc[dir](m,x.xq+2,k,j,i);
  const Real w = x.hi ? 0.0 : x.fc[dir](m,CLF_W,k,j,i);
  const int kl = k - (dir == 2), jl = j - (dir == 1), il = i - (dir == 0);
  const Real th = (1.0 - w)*x.fc[dir](m,CLF_TH,k,j,i);
  // (clean_stag_vimp: G and the beta term take theta (1 - b), and b Q joins them)
  const Real qd = x.drag ? x.fc[dir](m,CLF_DB,k,j,i) : 0.0, tg = th*(1.0 - qd);
  Real s = tg*x.fc[dir](m,CLF_G,k,j,i);
  if (x.fd && !x.hi) s += th*x.fc[dir](m,CLF_N+3,k,j,i);   // (fkin: of F^H, the last F^H)
  if (CleanStagMemUp(x, dir, m, k, j, i) != 0.0) s = 0.0;   // (gflux: upwind, implicit)
  if (x.feedback) {   // (the emission part of the beta term, CleanFaceWeights)
    s += 0.25*tg*x.c*(CleanStagN(x, CLW_BT, dir, m, k, j, i, kl, jl, il) +
                      CleanStagN(x, CLW_BT, dir, m, k, j, i, k, j, i))*
         (x.w(m,CLW_RE,kl,jl,il) + x.w(m,CLW_RE,k,j,i));
    // (less the emission's linear response at E_k, CleanFaceWeights)
    {
      auto rk = [&](const int kk, const int jj, const int ii) {
        return x.w(m,CLW_GX,kk,jj,ii)*x.op(m,VET_OP_E,kk,jj,ii)*x.w(m,CLW_EK,kk,jj,ii);
      };
      s -= 0.25*tg*x.c*(CleanStagN(x, CLW_BT, dir, m, k, j, i, kl, jl, il) +
                        CleanStagN(x, CLW_BT, dir, m, k, j, i, k, j, i))*x.h*x.c*
           (rk(kl, jl, il) + rk(k, j, i));
    }
  }
  if (x.drag) s += (1.0 - w)*qd*x.fc[dir](m,CLF_DQ,k,j,i);
  if (!x.hi) s += CleanThinFaceDG(x, dir, m, k, j, i, w);   // (CleanStagBracket's DG)
  // (the advective share 1 - theta_d of an upwind face, CleanFaceWeights)
  if (x.upw && x.feedback) {
    Real al, ar, phi = 0.0;   // (phw: the caller's phi of these args, CleanFaceFlux)
    if (phw)
      phi = *phw;
    else
      CleanFaceWeights(x, dir, m, k, j, i, al, ar, &phi);
    if (phi > 0.0) {
      s += phi*((1.0 - w) - th)*(1.0 - qd)*CleanStagAdvectCorr(x, dir, m, k, j, i);
    }
  }
  return s;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStagXSplit
//! \brief clean_stag_xdonor: split the lagged cross flux C(E_k) (slot CLF_XC) and face
//! source S (from xr) into a donor part q E_donor (q <= c; slots xq, xq + 1) and an
//! explicit part (xq + 2).  The explicit share s keeps the outflow h (A/V) |s L| <=
//! CLN_MEMUP E_old,donor, so the RHS stays >= 0 and the operator an M-matrix

KOKKOS_INLINE_FUNCTION
void CleanStagXSplit(const CleanArgs &x, const CleanArgs &xr, const int dir, const int m,
                     const int k, const int j, const int i) {
  const auto &f = x.fc[dir];
  const Real sr = CleanFaceSource(xr, dir, m, k, j, i);
  Real ql = 0.0, qr = 0.0, ex = 0.0;
  if (CleanBoundaryFace(x, dir, m, k, j, i)) {
    ex = sr;
  } else {
    const int kl = k - (dir == 2), jl = j - (dir == 1), il = i - (dir == 0);
    // (explicit share s of a lagged flux l: outflow at most CLN_MEMUP E_old,donor)
    auto share = [&](const Real l) {
      if (l == 0.0) return 1.0;
      const bool up = (l > 0.0);
      const Real eo = fmax(up ? x.w(m,CLW_EO,kl,jl,il) : x.w(m,CLW_EO,k,j,i), 0.0);
      // (cubed sphere: A_f over the owner's V, so a seam face takes one share)
      const Real out = x.h*fabs(l)*(!x.cub ? (up ? CleanDivW(x, dir, m, kl, jl, il, 1) :
                                                   CleanDivW(x, dir, m, k, j, i, 0)) :
                     ((dir == 0) ? x.ae0 : x.area[dir])(m,k,j,i)/
                     (up ? x.g(m,CLG_F+3,kl,jl,il) : x.g(m,CLG_F+3,k,j,i)));
      return (out > CLN_MEMUP*eo) ? CLN_MEMUP*eo/out : 1.0;
    };
    // (donor part li: q = |li|/E_donor <= c, E_k for the cross flux (exact at the fixed
    // point), E_old for a face source (one linear map for the step))
    auto put = [&](const Real li, const bool src) {
      if (li == 0.0) return;
      const bool up = (li > 0.0);
      const int e = src ? CLW_EO : CLW_EK;
      const Real ed = up ? x.w(m,e,kl,jl,il) : x.w(m,e,k,j,i);
      const Real q = (fabs(li) < x.c*ed) ? fabs(li)/ed : x.c;
      if (up) ql += q;
      if (!up) qr -= q;
    };
    const Real cf = f(m,CLF_XC,k,j,i);
    const Real s = share(cf + sr);
    ex += s*(cf + sr);
    put((1.0 - s)*cf, false);
    put((1.0 - s)*sr, true);
  }
  f(m,x.xq,k,j,i) = ql;
  f(m,x.xq+1,k,j,i) = qr;
  f(m,x.xq+2,k,j,i) = ex;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanDirCase
//! \brief f(D) with dir as a compile-time constant (fewer index ops, bitwise)

template <typename F>
KOKKOS_FORCEINLINE_FUNCTION
void CleanDirCase(const int dir, const F &f) {
  if (dir == 0) {
    f(std::integral_constant<int, 0>());
  } else if (dir == 1) {
    f(std::integral_constant<int, 1>());
  } else {
    f(std::integral_constant<int, 2>());
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFaceFlux
//! \brief flux through face (k,j,i) normal to dir for E in work vector e; with_source
//! adds the E-independent part; D >= 0: dir = D at compile time (CleanDirCase)

template <int D = -1>
CLEAN_DEVICE_NOINLINE CLEAN_HIP_INLINE KOKKOS_INLINE_FUNCTION
Real CleanFaceFlux(const CleanArgs &x, const int e, const int dn, const int m,
                   const int k, const int j, const int i, const bool with_source) {
  const int dir = (D < 0) ? dn : D;
  const int kl = k - (dir == 2), jl = j - (dir == 1), il = i - (dir == 0);
  Real al, ar, phi = 0.0;   // (phi: the source's upwind share, weights once)
  CleanFaceWeights(x, dir, m, k, j, i, al, ar, &phi);
  Real f = CleanFaceCross(x, e, dir, m, k, j, i);
  if (al != 0.0) f += al*x.w(m,e,kl,jl,il);
  if (ar != 0.0) f += ar*x.w(m,e,k,j,i);
  if (with_source) f += CleanFaceSource(x, dir, m, k, j, i, &phi);
  return f;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFaceFluxCall, CleanFaceStressCall
//! \brief the face flux and stress as HIP calls, for paths a pass does not take

CLEAN_HIP_NOINLINE KOKKOS_INLINE_FUNCTION
Real CleanFaceFluxCall(const CleanArgs &x, const int e, const int dir, const int m,
                       const int k, const int j, const int i, const bool with_source) {
  return CleanFaceFlux(x, e, dir, m, k, j, i, with_source);
}

CLEAN_HIP_NOINLINE KOKKOS_INLINE_FUNCTION
Real CleanFaceStressCall(const CleanArgs &x, const int e, const int dir, const int m,
                         const int k, const int j, const int i, const bool hi) {
  return CleanFaceStress(x, e, dir, m, k, j, i, hi);
}

}  // namespace radvet

#endif  // RAD_VET_RAD_VET_FACES_HPP_
