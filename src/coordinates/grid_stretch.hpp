#ifndef COORDINATES_GRID_STRETCH_HPP_
#define COORDINATES_GRID_STRETCH_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file grid_stretch.hpp
//! \brief Coordinate stretching maps for the radial (x1) and polar (x2) directions of the
//! spherical grids, shared by Coordinates and the problem generators so both build the
//! same grid.

#include <cmath>

#include "athena.hpp"

// Layout of Mesh::fStretchRPoly (see StretchRPoly): NSTRETCH_R_PCOEF polynomial
// coefficients, NSTRETCH_R_BUMP bumps of (amplitude, centre, width), NSTRETCH_R_NPLAT
// plateaus of (amplitude, xa, xb, width), the log tail (s, xi_j, r1_old) and the
// geometric head (xi_h, s, xi_j, q, nh).
#define NSTRETCH_R_PCOEF 8
#define NSTRETCH_R_BUMP 2
#define NSTRETCH_R_NPLAT 2
#define NSTRETCH_R_PLAT (4*NSTRETCH_R_NPLAT)
#define NSTRETCH_R_TAIL 3
#define NSTRETCH_R_HEAD 5
#define NSTRETCH_R_POLY (NSTRETCH_R_PCOEF + 3*NSTRETCH_R_BUMP + NSTRETCH_R_PLAT + \
                         NSTRETCH_R_TAIL + NSTRETCH_R_HEAD)
#define NSTRETCH_R_TAIL_OFF (NSTRETCH_R_PCOEF + 3*NSTRETCH_R_BUMP + NSTRETCH_R_PLAT)
#define NSTRETCH_R_HEAD_OFF (NSTRETCH_R_TAIL_OFF + NSTRETCH_R_TAIL)

//----------------------------------------------------------------------------------------
//! \fn Real LogCoshStable
//! \brief ln cosh(x) without overflow: |x| + ln(1 + exp(-2|x|)) - ln 2
KOKKOS_INLINE_FUNCTION
Real LogCoshStable(const Real x) {
  const Real ax = fabs(x);
  return ax + log1p(exp(-2.0*ax)) - 0.69314718055994530942;
}

//----------------------------------------------------------------------------------------
//! \fn Real StretchRPlateauG
//! \brief Integral G(xi) of the plateau profile g = (tanh((xi-xa)/w) - tanh((xi-xb)/w))/2
//! used by StretchRPoly: du/dxi += a [g(xi) - G(1)], u += a [G(xi) - xi G(1)], so that
//! u(0) = 0 and u(1) = 1. a < 0 gives finer cells on [xa, xb].
KOKKOS_INLINE_FUNCTION
Real StretchRPlateauG(const Real xa, const Real xb, const Real w, const Real xi) {
  return 0.5*w*(LogCoshStable((xi - xa)/w) - LogCoshStable(-xa/w)
                - LogCoshStable((xi - xb)/w) + LogCoshStable(-xb/w));
}

//----------------------------------------------------------------------------------------
//! \fn void StretchR
//! \brief Exponential radial stretch r = r0 + (r1-r0) (1 - exp(-a xi))/(1 - exp(-a)).

KOKKOS_INLINE_FUNCTION
void StretchR(const Real a, const Real r0, const Real r1, Real &r) {
  Real xi = (r-r0)/(r1-r0);
  Real denom = 1.0 - exp(-a);
  r = r0 + (r1 - r0)*(1.0 - exp(-a*xi))/denom;
}
//--------------------------------------------------------------------------------------
//! \fn Real StretchRPolyU
//! \brief Map u(xi) of the polynomial radial stretch r = r0 + (r1 - r0) u(xi),
//! xi = (r - r0)/(r1 - r0),
//! u(xi) = xi + sum_k c_k xi^k (1 - xi), plus optional bumps and plateaus. u(0) = 0 and
//! u(1) = 1 for any coefficients. Allows a non-monotonic cell width (e.g. constant cells
//! per pressure scale height). Mesh checks that the map is strictly increasing.
KOKKOS_INLINE_FUNCTION
Real StretchRPolyU(const Real *c, const Real xi) {
  Real u = xi;
  Real xik = xi;                      // xi^k, built up as k increases
  for (int k=1; k<=NSTRETCH_R_PCOEF; ++k) {
    u += c[k-1]*xik*(1.0-xi);
    xik *= xi;
  }
  // optional bumps: du/dxi changes by a sech^2((xi - xb)/w), corrected so u(0) = 0 and
  // u(1) = 1; a < 0 gives finer cells near xb
  for (int b=0; b<NSTRETCH_R_BUMP; ++b) {
    const Real a = c[NSTRETCH_R_PCOEF + 3*b];
    if (a == 0.0) continue;
    const Real xb = c[NSTRETCH_R_PCOEF + 3*b + 1], w = c[NSTRETCH_R_PCOEF + 3*b + 2];
    u += a*w*(tanh((xi - xb)/w) - (1.0 - xi)*tanh(-xb/w) - xi*tanh((1.0 - xb)/w));
  }
  // optional plateaus (see StretchRPlateauG)
  for (int q=0; q<NSTRETCH_R_NPLAT; ++q) {
    const Real *cp = c + NSTRETCH_R_PCOEF + 3*NSTRETCH_R_BUMP + 4*q;
    if (cp[0] != 0.0) {
      u += cp[0]*(StretchRPlateauG(cp[1], cp[2], cp[3], xi)
                  - xi*StretchRPlateauG(cp[1], cp[2], cp[3], 1.0));
    }
  }
  return u;
}

//----------------------------------------------------------------------------------------
//! \fn void StretchRPoly
//! \brief Polynomial radial stretch with an optional geometric head or log tail, which
//! keep the faces of an older run's grid (nin cells) and extend the domain:
//! - head (xi_h > 0): faces of the old grid kept from face ij up; below, nh cells of
//!   constant width ratio q, r = r0 + (r_j - r0) (1 - q^k)/(1 - q^nh), k = (xi/xi_h) nh.
//! - tail (s > 0): faces of the old grid kept up to face ij; above, constant face ratio,
//!   r = r_j (r1/r_j)^((xi - xi_j)/(1 - xi_j)).
//! The slot values are set by Mesh.
KOKKOS_INLINE_FUNCTION
void StretchRPoly(const Real *c, const Real r0, const Real r1, Real &r) {
  const Real xi = (r-r0)/(r1-r0);
  const Real *chd = c + NSTRETCH_R_HEAD_OFF;
  if (chd[0] > 0.0) {
    const Real xh = chd[0], s = chd[1], xj = chd[2], q = chd[3], nh = chd[4];
    if (xi <= xh) {
      const Real rj = r0 + (r1 - r0)*StretchRPolyU(c, xj);
      const Real k = (xi/xh)*nh;
      const Real lq = log(q);
      const Real fr = (lq == 0.0) ? k/nh : expm1(k*lq)/expm1(nh*lq);
      r = r0 + (rj - r0)*fr;
    } else {
      r = r0 + (r1 - r0)*StretchRPolyU(c, xj + (xi - xh)*s);
    }
    return;
  }
  const Real *ct = c + NSTRETCH_R_TAIL_OFF;
  if (ct[0] > 0.0) {
    const Real s = ct[0], xj = ct[1], r1o = ct[2];
    if (xi <= xj) {
      r = r0 + (r1o - r0)*StretchRPolyU(c, xi*s);
    } else {
      const Real rj = r0 + (r1o - r0)*StretchRPolyU(c, xj*s);
      r = rj*exp(log(r1/rj)*(xi - xj)/(1.0 - xj));
    }
    return;
  }
  r = r0 + (r1-r0)*StretchRPolyU(c, xi);
}
//----------------------------------------------------------------------------------------
//! \fn void ApplyRStretch
//! \brief Applies the active radial stretch to an unstretched position from
//! CellCenterX/LeftEdgeX.

KOKKOS_INLINE_FUNCTION
void ApplyRStretch(const bool str_r, const Real fstr_r, const bool str_rp, const Real *c,
                   const Real r0, const Real r1, Real &r) {
  if (str_r) {
    StretchR(fstr_r, r0, r1, r);
  } else if (str_rp) {
    StretchRPoly(c, r0, r1, r);
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real RadialCentroid
//! \brief Volume centroid of the shell [r_l, r_r], stored in Coordinates::x1v on both
//! spherical grids.

KOKKOS_INLINE_FUNCTION
Real RadialCentroid(const Real r_l, const Real r_r) {
  const Real q = r_l/r_r;
  return 0.25*(q*q + 1.0)/((1.0/3.0)*(q*q + q + 1.0))*(r_r + r_l);
}

//----------------------------------------------------------------------------------------
//! \fn void StretchTheta
//! \brief Sinh polar-angle stretch t = pi/2 (1 + sinh(a (2 t/pi - 1))/sinh(a)).

KOKKOS_INLINE_FUNCTION
void StretchTheta(const Real a, Real &t) {
  Real xi = t/M_PI;
  t = M_PI/2.0*(1.0+sinh(a*(2.0*xi-1.0))/sinh(a));
}

//----------------------------------------------------------------------------------------
//! \fn Real ThetaPlateauHalf
//! \brief Polar-angle offset from the mesh centre for e >= 0 cells of the theta-plateau
//! map: a uniform core of ec cells of width d, then Np geometric cap cells (ratio q)
//! ending exactly at the mesh edge, mirrored about the edge in the ghost zones.
//! p[] = Mesh::fStretchThetaPlat: (x2min, x2max, nx2/2, ec, d, q, Np).
KOKKOS_INLINE_FUNCTION
Real ThetaPlateauHalf(const Real *p, Real e) {   // e >= 0, e <= 2*nh
  const Real nh = p[2], ec = p[3], d = p[4], q = p[5], np = p[6];
  const Real half = 0.5*(p[1] - p[0]);
  const Real hc = d*ec;
  const bool mirror = (e > nh);           // ghost beyond the mesh edge: mirror image
  if (mirror) e = 2.0*nh - e;
  Real h;
  if (e <= ec) {
    h = d*e;
  } else {
    const Real lq = log(q);
    h = hc + (half - hc)*expm1((e - ec)*lq)/expm1(np*lq);
  }
  return mirror ? 2.0*half - h : h;
}

//----------------------------------------------------------------------------------------
//! \fn void StretchThetaPlateau
//! \brief Theta map with a uniform equatorial core and geometric caps
//! (mesh/theta_plateau), symmetric about the mesh centre.

KOKKOS_INLINE_FUNCTION
void StretchThetaPlateau(const Real *p, Real &t) {
  const Real mid = 0.5*(p[0] + p[1]);
  const Real e = (t - p[0])/(p[1] - p[0])*(2.0*p[2]) - p[2];
  const Real h = ThetaPlateauHalf(p, fabs(e));
  t = (e < 0.0) ? mid - h : mid + h;
}

//----------------------------------------------------------------------------------------
//! \fn void ApplyThetaStretch
//! \brief Applies the theta-plateau map if p is set (p[2] > 0), else the sinh map.

KOKKOS_INLINE_FUNCTION
void ApplyThetaStretch(const Real a, const Real *p, Real &t) {
  if (p != nullptr && p[2] > 0.0) {
    StretchThetaPlateau(p, t);
  } else {
    StretchTheta(a, t);
  }
}

#endif // COORDINATES_GRID_STRETCH_HPP_
