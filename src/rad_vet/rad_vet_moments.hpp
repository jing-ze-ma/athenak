#ifndef RAD_VET_RAD_VET_MOMENTS_HPP_
#define RAD_VET_RAD_VET_MOMENTS_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_moments.hpp
//! \brief kernels shared by the closure files: moments into the closure slots,
//! realizability, face and local vectors, direction rows.

#include <cstdint>

#include "athena.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_momentum.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void CleanStoreMoments
//! \brief a cell's closure from its ray moments q (CLN_*, not normalised): D = K/J,
//! h+- = H+-/J, k+ = K+/J, C1 = <|I - S|>/J; J = 0: isotropic, marked for the fill

KOKKOS_INLINE_FUNCTION
void CleanStoreMoments(const DvceArray5D<Real> &cs, const DvceArray5D<Real> &valid,
                       const Real *q, const int m, const int k, const int j,
                       const int i) {
  const int ns = cs.extent_int(1);   // (the staggered layout: no k+, CleanInit)
  if (q[CLN_J] > 0.0) {
    const Real ij = 1.0/q[CLN_J];
    for (int a = 0; a < 6; ++a) cs(m,CLN_D11+a,k,j,i) = q[CLN_D11+a]*ij;
    for (int a = 0; a < 3; ++a) {
      cs(m,CLN_HP+a,k,j,i) = fmax(q[CLN_HP+a]*ij, 0.0);
      cs(m,CLN_HM+a,k,j,i) = fmin(q[CLN_HM+a]*ij, 0.0);
    }
    for (int a = 0; a < 9 && CLN_KP + a < ns; ++a) cs(m,CLN_KP+a,k,j,i) = q[CLN_KP+a]*ij;
    cs(m,CLN_C1,k,j,i) = fmax(q[CLN_C1]*ij, 0.0);
    cs(m,CLN_J,k,j,i) = q[CLN_J];
    valid(m,0,k,j,i) = 1.0;
  } else {
    for (int a = 0; a < 3; ++a) {
      cs(m,CLN_D11+a,k,j,i) = 1.0/3.0;
      cs(m,CLN_D11+3+a,k,j,i) = 0.0;
      cs(m,CLN_HP+a,k,j,i) = 0.25;
      cs(m,CLN_HM+a,k,j,i) = -0.25;
    }
    for (int a = 0; a < 9 && CLN_KP + a < ns; ++a) {
      cs(m,CLN_KP+a,k,j,i) = (ns == CLN_NHO) ? ((a < 3) ? 0.25 : -0.25) :
                             ((a%4 == 0) ? 1.0/6.0 : 0.0);
    }
    cs(m,CLN_C1,k,j,i) = 1.0;
    cs(m,CLN_J,k,j,i) = 0.0;
    valid(m,0,k,j,i) = 0.0;
  }
  // (a fixed trip count: q's indices known at compile time, q in registers)
  for (int n = CLN_F2; n < CLN_KP; ++n) {
    if (n < ns) cs(m,n,k,j,i) = q[n];
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanFaceVector
//! \brief exact vector area a of the low face of cell (k,j,i) normal to dir (spherical
//! polar exact, cubed sphere CleanCsGeometry's, Cartesian along e_dir)

KOKKOS_INLINE_FUNCTION
void CleanFaceVector(const CleanArgs &x, const DvceArray2D<Real> &x3f, const int dir,
                     const int m, const int k, const int j, const int i, Real a[3]) {
  if (x.cub) {   // the exact closed face vectors of CleanCsGeometry
    for (int c = 0; c < 3; ++c)
      a[c] = x.area[dir](m,k,j,i)*x.g.At<CLG_N>(m,3*dir+c,k,j,i);
    return;
  }
  if (!x.sph) {
    for (int c = 0; c < 3; ++c) a[c] = (c == dir) ? x.area[dir](m,k,j,i) : 0.0;
    return;
  }
  // (a theta or phi face reads only its own angle: no j+1 (k+1) edge in 1-D (2-D))
  const Real r1 = x.x1f(m,i), r2 = x.x1f(m,i+1), t1 = x.x2f(m,j);
  const Real t2 = (dir == 1) ? t1 : x.x2f(m,j+1);
  const Real p1 = x3f(m,k), p2 = (dir == 2) ? p1 : x3f(m,k+1);
  const Real hr = 0.5*(r2*r2 - r1*r1);
  if (dir == 0) {
    const Real s2 = 0.5*(t2 - t1) - 0.25*(sin(2.0*t2) - sin(2.0*t1));
    a[0] = r1*r1*s2*(sin(p2) - sin(p1));
    a[1] = r1*r1*s2*(cos(p1) - cos(p2));
    a[2] = r1*r1*0.5*(SQR(sin(t2)) - SQR(sin(t1)))*(p2 - p1);
  } else if (dir == 1) {
    a[0] = sin(t1)*hr*cos(t1)*(sin(p2) - sin(p1));
    a[1] = sin(t1)*hr*cos(t1)*(cos(p1) - cos(p2));
    a[2] = -sin(t1)*hr*sin(t1)*(p2 - p1);
  } else {
    a[0] = -hr*(t2 - t1)*sin(p1);
    a[1] = hr*(t2 - t1)*cos(p1);
    a[2] = 0.0;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanLocalBasis
//! \brief the local orthonormal basis (r_hat, theta_hat, phi_hat) at (th, ph)

KOKKOS_INLINE_FUNCTION
void CleanLocalBasis(const Real th, const Real ph, Real e[3][3]) {
  const Real st = sin(th), ct = cos(th), sp = sin(ph), cp = cos(ph);
  e[0][0] = st*cp; e[0][1] = st*sp; e[0][2] = ct;
  e[1][0] = ct*cp; e[1][1] = ct*sp; e[1][2] = -st;
  e[2][0] = -sp;   e[2][1] = cp;    e[2][2] = 0.0;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanThickI
//! \brief clean_ray_cut: intensity along n in a cell below the cut, the diffusion limit
//! I = J + 3 n . H of its E en and F f in the frame e (CleanThickFrame)

KOKKOS_INLINE_FUNCTION
void CleanThickFrame(const CleanArgs &x, const DvceArray2D<Real> &x3v, const int m,
                     const int k, const int j, const int i, Real e[3][3]) {
  if (x.cub) {
    CleanCsFrame(x, m, k, j, i, e);
  } else {
    CleanLocalBasis(x.x2v(m,j), x3v(m,k), e);
  }
}

KOKKOS_INLINE_FUNCTION
Real CleanThickI(const Real e[3][3], const Real f[3], const Real en, const Real ic,
                 const Real n[3]) {
  Real nf = 0.0;
  for (int b = 0; b < 3; ++b) nf += (n[0]*e[b][0] + n[1]*e[b][1] + n[2]*e[b][2])*f[b];
  return fmax(en + 3.0*ic*nf, 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn bool CleanRealizable
//! \brief true if moments q (J > 0) are realisable: |H| <= J, K - H H^T/J PSD (to 1e-12)

KOKKOS_INLINE_FUNCTION
bool CleanRealizable(const CleanArgs &x, const Real *q, const int m, const int k,
                     const int j, const int i) {
  const Real ij = 1.0/q[CLN_J];
  Real s[3], h[3], a[3][3];
  for (int b = 0; b < 3; ++b) s[b] = (q[CLN_HP+b] + q[CLN_HM+b])*ij;
  CleanRayFlux(x, m, k, j, i, s, h);
  a[0][0] = q[CLN_D11]*ij;
  a[1][1] = q[CLN_D11+1]*ij;
  a[2][2] = q[CLN_D11+2]*ij;
  a[0][1] = a[1][0] = q[CLN_D11+3]*ij;
  a[0][2] = a[2][0] = q[CLN_D11+4]*ij;
  a[1][2] = a[2][1] = q[CLN_D11+5]*ij;
  for (int b = 0; b < 3; ++b) {
    for (int c = 0; c < 3; ++c) a[b][c] -= h[b]*h[c];
  }
  // all principal minors of D - h h^T non-negative
  const Real eps = 1.0e-12;
  const Real m01 = a[0][0]*a[1][1] - a[0][1]*a[0][1];
  const Real m02 = a[0][0]*a[2][2] - a[0][2]*a[0][2];
  const Real m12 = a[1][1]*a[2][2] - a[1][2]*a[1][2];
  const Real det = a[0][0]*m12 - a[0][1]*(a[0][1]*a[2][2] - a[1][2]*a[0][2])
                   + a[0][2]*(a[0][1]*a[1][2] - a[1][1]*a[0][2]);
  return (h[0]*h[0] + h[1]*h[1] + h[2]*h[2] <= 1.0 + eps) &&
         (a[0][0] >= -eps) && (a[1][1] >= -eps) && (a[2][2] >= -eps) &&
         (m01 >= -eps) && (m02 >= -eps) && (m12 >= -eps) && (det >= -eps);
}

// one cell's directions as a strided row: the address once per cell (dummy: stride 0)
struct CleanRow {
  Real *p;
  int64_t s;
  KOKKOS_INLINE_FUNCTION Real operator[](const int d) const {return p[d*s];}
};

KOKKOS_INLINE_FUNCTION
CleanRow VgdRow(const RadiationVet::VgdRag &vi, const int m, const int kb, const int jb,
                const int i) {
  Real *p = &vi(m,0,kb,jb,i);
  return CleanRow{p, (p == vi.d.data() + vi.dummy) ? 0 : 1};
}

}  // namespace radvet

#endif  // RAD_VET_RAD_VET_MOMENTS_HPP_
