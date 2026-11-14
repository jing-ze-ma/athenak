#ifndef BVALS_CS_SEAM_FC3_HPP_
#define BVALS_CS_SEAM_FC3_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cs_seam_fc3.hpp
//! \brief Weights and gnomonic measures for the third-order cubed-sphere seam ghost face
//! fields (<mesh>/cs_seam_fc_order = 3, FCSeamPack in bvals_fc.cpp). Each ghost face
//! value is the average of B.n_dst over the receiver's own face (2-point Gauss weighted
//! by the arc or solid-angle density) of third-order point values reconstructed from the
//! source face/cell integrals. The radial extent factors out, so all measures are per
//! unit shell.

#include <cmath>
#include "athena.hpp"

namespace cs_fc3 {

//! 2-point Gauss abscissa on [-1,1] (equal weights)
constexpr Real kGauss = 0.57735026918962576451;

//----------------------------------------------------------------------------------------
//! \fn void cs_fc3::QuadDensW
//! \brief Weights of the quadratic density q(t) = w0 I0 + w1 I1 + w2 I2 whose integrals
//! over the cells [-1.5,-0.5], [-0.5,0.5], [0.5,1.5] are I0, I1, I2.

KOKKOS_FORCEINLINE_FUNCTION
void QuadDensW(const Real t, Real w[3]) {
  const Real c = 0.5*(t*t - 1.0/12.0);
  w[0] = c - 0.5*t;
  w[1] = 1.0 - 2.0*c;
  w[2] = c + 0.5*t;
}

//----------------------------------------------------------------------------------------
//! \fn void cs_fc3::Lag4W
//! \brief Cubic Lagrange weights through the nodes 0,1,2,3 at t.

KOKKOS_FORCEINLINE_FUNCTION
void Lag4W(const Real t, Real w[4]) {
  w[0] = -(t-1.0)*(t-2.0)*(t-3.0)/6.0;
  w[1] =  t*(t-2.0)*(t-3.0)/2.0;
  w[2] = -t*(t-1.0)*(t-3.0)/2.0;
  w[3] =  t*(t-1.0)*(t-2.0)/6.0;
}

//----------------------------------------------------------------------------------------
//! \fn Real cs_fc3::Arc
//! \brief Great-circle arc of the gnomonic line at fixed coordinate c (C = tan c) between
//! the tangents A0 < A1 of the other coordinate.

KOKKOS_FORCEINLINE_FUNCTION
Real Arc(const Real C, const Real A0, const Real A1) {
  const Real s = sqrt(1.0 + C*C);
  return atan(A1/s) - atan(A0/s);
}

//----------------------------------------------------------------------------------------
//! \fn Real cs_fc3::ArcDens
//! \brief d(arc)/d(a) along the line at fixed c (C = tan c, A = tan a).

KOKKOS_FORCEINLINE_FUNCTION
Real ArcDens(const Real C, const Real A) {
  return sqrt(1.0 + C*C)*(1.0 + A*A)/(1.0 + C*C + A*A);
}

//----------------------------------------------------------------------------------------
//! \fn Real cs_fc3::Solid
//! \brief Exact solid angle of the gnomonic cell [X0,X1] x [Y0,Y1] (tangents).

KOKKOS_FORCEINLINE_FUNCTION
Real Solid(const Real X0, const Real X1, const Real Y0, const Real Y1) {
  auto F = [](const Real x, const Real y) { return atan(x*y/sqrt(1.0 + x*x + y*y)); };
  return fabs(F(X1,Y1) - F(X0,Y1) - F(X1,Y0) + F(X0,Y0));
}

//----------------------------------------------------------------------------------------
//! \fn Real cs_fc3::SolidDens
//! \brief Solid-angle density d(Omega)/(dx dy) at (X, Y) = (tan x, tan y).

KOKKOS_FORCEINLINE_FUNCTION
Real SolidDens(const Real X, const Real Y) {
  const Real d = 1.0 + X*X + Y*Y;
  return (1.0 + X*X)*(1.0 + Y*Y)/(d*sqrt(d));
}


//----------------------------------------------------------------------------------------
//! \fn Real cs_fc3::SmoothW
//! \brief Smoothness weight w = max(0, 1 - r) of a higher-order correction, with
//! r = |fa - 2 fb + fc|/(|fb - fa| + |fc - fb|) from three equally spaced samples: ~1 on
//! smooth data, 0 across a jump.

KOKKOS_FORCEINLINE_FUNCTION
Real SmoothW(const Real fa, const Real fb, const Real fc) {
  const Real d1 = fb - fa, d2 = fc - fb;
  const Real s = fabs(d1) + fabs(d2);
  if (!(s > 0.0)) return 1.0;               // constant data
  const Real r = fabs(d2 - d1)/s;
  return fmax(0.0, 1.0 - r);
}

//----------------------------------------------------------------------------------------
//! \fn Real cs_fc3::CurvCorrect
//! \brief Raises a linear-exact value to quadratic exactness by subtracting coef times
//! the second difference of (fa, fb, fc), gated by SmoothW and bounded to [lo, hi].

KOKKOS_FORCEINLINE_FUNCTION
Real CurvCorrect(const Real base, const Real fa, const Real fb, const Real fc,
                 const Real coef, const Real lo, const Real hi) {
  const Real v = base - coef*SmoothW(fa, fb, fc)*(fa - 2.0*fb + fc);
  return fmin(hi, fmax(lo, v));
}

}  // namespace cs_fc3
#endif  // BVALS_CS_SEAM_FC3_HPP_
