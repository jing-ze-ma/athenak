#ifndef RECONSTRUCT_WB_PLM_BLEND_HPP_
#define RECONSTRUCT_WB_PLM_BLEND_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file wb_plm_blend.hpp
//! \brief Continuous fallback and range bound for the well-balanced x1 PLM
//! (Hydro/MHD::GridPiecewiseLinearX1 and GridPiecewiseLinearDerX1).
//!
//! The default rule takes the plain-PLM faces whenever a well-balanced face is negative
//! or the data have an extremum (sL*sR <= 0).  Both tests are threshold switches: near a
//! flat top 1 ulp of input can flip a face between the two reconstructions.  With
//! <hydro|mhd>/wb_plm_smooth = true the faces are instead blended,
//!   q = w q_wb + (1 - w) q_plm,
//!   w = S(th / kWbThExt) * S(min(q_wb,L, q_wb,R) / (kWbPosFrac |q_i|)),
//!   th = 2 sL sR / (sL^2 + sR^2) in [-1, 1],   S(x) = smoothstep on [0, 1].
//! w is exactly 1 where the default rule keeps the well-balanced faces with margin (so a
//! hydrostatic column stays balanced to round-off), exactly 0 where it falls back, and
//! C1 in the data in between.  The blend is positive because w > 0 needs both
//! well-balanced faces > 0.

#include <cmath>
#include "athena.hpp"

namespace wbplm {
constexpr Real kWbThExt = 0.2;     // th at which the extremum factor reaches 1
constexpr Real kWbPosFrac = 0.1;   // face/cell value at which the positivity factor is 1

//----------------------------------------------------------------------------------------
//! \fn Real wbplm::SmoothStep01
//! \brief smoothstep t^2 (3 - 2t) of x clamped to [0, 1]

KOKKOS_INLINE_FUNCTION
Real SmoothStep01(const Real x) {
  const Real t = fmin(fmax(x, static_cast<Real>(0.0)), static_cast<Real>(1.0));
  return t*t*(3.0 - 2.0*t);
}

//----------------------------------------------------------------------------------------
//! \fn Real wbplm::WbPlmWeight
//! \brief weight of the well-balanced faces (qwl = q_{i+1/2,L}, qwr = q_{i-1/2,R})

KOKKOS_INLINE_FUNCTION
Real WbPlmWeight(const Real q_i, const Real sL, const Real sR, const Real qwl,
                 const Real qwr) {
  // extremum factor
  const Real den = sL*sL + sR*sR;
  const Real th = (den > 0.0) ? 2.0*sL*sR/den : 0.0;
  Real w = SmoothStep01(th/kWbThExt);
  // positivity factor
  const Real qmin = fmin(qwl, qwr);
  const Real sc = kWbPosFrac*fabs(q_i);
  w *= (sc > 0.0) ? SmoothStep01(qmin/sc) : ((qmin > 0.0) ? 1.0 : 0.0);
  return w;   // NaN inputs give 0: fmax(NaN, 0) = 0 in SmoothStep01
}

//----------------------------------------------------------------------------------------
//! \fn void wbplm::WbPlmBlend
//! \brief blends the well-balanced faces (ql, qr) in place with the plain-PLM faces
//! (pl, pr)

KOKKOS_INLINE_FUNCTION
void WbPlmBlend(const Real q_i, const Real sL, const Real sR, const Real pl,
                const Real pr, Real &ql, Real &qr) {
  const Real w = WbPlmWeight(q_i, sL, sR, ql, qr);
  if (w >= 1.0) return;
  if (w <= 0.0) {
    ql = pl;
    qr = pr;
    return;
  }
  ql = w*ql + (1.0 - w)*pl;
  qr = w*qr + (1.0 - w)*pr;
}

//----------------------------------------------------------------------------------------
//! \fn void wbplm::WbPlmRange
//! \brief <hydro|mhd>/wb_plm_range = true: bounds both faces by the data range of the
//! 3-cell stencil.  The dynamic background is a hydrostatic profile through cell i, so
//! where the column is far from hydrostatic (floor-cold or outflowing top cells) its
//! faces can leave the data range.  In a hydrostatic column the background faces lie
//! between the neighbours and the bound never acts.

KOKKOS_INLINE_FUNCTION
void WbPlmRange(const Real qim1, const Real qi, const Real qip1, Real &ql, Real &qr) {
  const Real lo = fmin(qim1, fmin(qi, qip1)), hi = fmax(qim1, fmax(qi, qip1));
  ql = fmin(fmax(ql, lo), hi);
  qr = fmin(fmax(qr, lo), hi);
}
}  // namespace wbplm
#endif  // RECONSTRUCT_WB_PLM_BLEND_HPP_
