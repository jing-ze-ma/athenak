#ifndef BVALS_CS_SEAM_FC3_HPP_
#define BVALS_CS_SEAM_FC3_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
//========================================================================================
//! \file cs_seam_fc3.hpp
//! \brief fcorder3-1103 <mesh>/cs_seam_fc_order = 3: the shared weights and gnomonic
//! measures of the third-order cubed-sphere seam ghost face fields (FCSeamPack::srcval3
//! in bvals_fc.cpp).  One set serves the along-seam AND the seam-normal direction.
//!
//! WHY.  b.x1f/x2f/x3f are face AVERAGES (flux / exact area, coordinates.cpp).  A ghost
//! face of the receiving chart is not a face of the sending chart: beyond a seam the
//! receiver's x2-type face (constant seam-normal angle) lies on a source face surface but
//! spans a different along-seam interval, its other angular face is an oblique curve in
//! the source chart, and its radial face a skewed quadrilateral.  The order-2 fill treats
//! the source averages as point values and resamples them, an O(h^2) error per ghost face
//! that the per-panel curl differences over one cell: a FIRST-order seam current.  Here
//! every ghost face is the average over ITS OWN face of third-order point values:
//!   point values  -- from the source face/cell INTEGRALS (average x exact measure): the
//!                    quadratic density with the three cell integrals (QuadDensW), one or
//!                    two directions, divided by the gnomonic measure density there; across
//!                    the staggered direction a cubic through 4 face lines (Lag4W);
//!   face average  -- 2-point Gauss (per direction) of B.n_dst over the receiver's face,
//!                    weighted by its own arc / solid-angle density.
//! The radial extent is common to both charts and factors out exactly (angular faces
//! r-moment x arc, radial faces r^2 x solid angle), so all of it is per unit shell.

#include <cmath>
#include "athena.hpp"

namespace cs_fc3 {

//! 2-point Gauss abscissa on [-1,1] (equal weights)
constexpr Real kGauss = 0.57735026918962576451;

//! Density weights: the quadratic q whose integrals over the cells [-1.5,-0.5],
//! [-0.5,0.5], [0.5,1.5] (cell units, middle cell centred at 0) are I0, I1, I2, at t:
//! q(t) = w0 I0 + w1 I1 + w2 I2 (exact for quadratic densities; third order).
KOKKOS_FORCEINLINE_FUNCTION
void QuadDensW(const Real t, Real w[3]) {
  const Real c = 0.5*(t*t - 1.0/12.0);
  w[0] = c - 0.5*t;
  w[1] = 1.0 - 2.0*c;
  w[2] = c + 0.5*t;
}

//! Cubic Lagrange through nodes 0,1,2,3 at t
KOKKOS_FORCEINLINE_FUNCTION
void Lag4W(const Real t, Real w[4]) {
  w[0] = -(t-1.0)*(t-2.0)*(t-3.0)/6.0;
  w[1] =  t*(t-2.0)*(t-3.0)/2.0;
  w[2] = -t*(t-1.0)*(t-3.0)/2.0;
  w[3] =  t*(t-1.0)*(t-2.0)/6.0;
}

//! Great-circle arc (radians) of the gnomonic line at fixed coordinate c (C = tan c)
//! between the other coordinate's tangents A0 < A1: the code's dth_* (coordinates.cpp)
//! in a form without the acos cancellation.
KOKKOS_FORCEINLINE_FUNCTION
Real Arc(const Real C, const Real A0, const Real A1) {
  const Real s = sqrt(1.0 + C*C);
  return atan(A1/s) - atan(A0/s);
}

//! d(arc)/d(a) along the line at fixed c (C = tan c, A = tan a)
KOKKOS_FORCEINLINE_FUNCTION
Real ArcDens(const Real C, const Real A) {
  return sqrt(1.0 + C*C)*(1.0 + A*A)/(1.0 + C*C + A*A);
}

//! exact solid angle of the gnomonic cell [X0,X1] x [Y0,Y1] (tangents), any order
KOKKOS_FORCEINLINE_FUNCTION
Real Solid(const Real X0, const Real X1, const Real Y0, const Real Y1) {
  auto F = [](const Real x, const Real y) { return atan(x*y/sqrt(1.0 + x*x + y*y)); };
  return fabs(F(X1,Y1) - F(X0,Y1) - F(X1,Y0) + F(X0,Y0));
}

//! d(Omega)/(dx dy) at (X, Y) = (tan x, tan y)
KOKKOS_FORCEINLINE_FUNCTION
Real SolidDens(const Real X, const Real Y) {
  const Real d = 1.0 + X*X + Y*Y;
  return (1.0 + X*X)*(1.0 + Y*Y)/(d*sqrt(d));
}


//! cornerfix3-1103: SHARED smoothness gate for a higher-order correction (one function
//! for the seam samples and the cube-vertex rule).  fa, fb, fc are three equally spaced
//! samples of the field along the correction's direction; their second difference
//! dd = fa - 2 fb + fc carries the curvature.  The smoothness ratio
//!   r = |dd| / (|fb - fa| + |fc - fb|)
//! is O(h) on smooth data and ~1 across a jump (or at an extremum the samples resolve
//! poorly), so the weight w = max(0, 1 - 2 r) keeps the full correction where the data
//! are smooth and switches it off continuously at a jump (no threshold flip-flop).
KOKKOS_FORCEINLINE_FUNCTION
Real SmoothW(const Real fa, const Real fb, const Real fc) {
  const Real d1 = fb - fa, d2 = fc - fb;
  const Real s = fabs(d1) + fabs(d2);
  if (!(s > 0.0)) return 1.0;               // constant data: the correction is 0 anyway
  const Real r = fabs(d2 - d1)/s;
  return fmax(0.0, 1.0 - 2.0*r);
}

//! cornerfix3-1103: a linear-exact base value `base` raised to quadratic exactness by
//! subtracting coef x (second difference of fa, fb, fc), gated by SmoothW and finally
//! bounded to [lo, hi] (the range of every value the estimate used): third-order point
//! value on smooth data, the base value at a jump, never an overshoot.
KOKKOS_FORCEINLINE_FUNCTION
Real CurvCorrect(const Real base, const Real fa, const Real fb, const Real fc,
                 const Real coef, const Real lo, const Real hi) {
  const Real v = base - coef*SmoothW(fa, fb, fc)*(fa - 2.0*fb + fc);
  return fmin(hi, fmax(lo, v));
}

}  // namespace cs_fc3
#endif  // BVALS_CS_SEAM_FC3_HPP_
