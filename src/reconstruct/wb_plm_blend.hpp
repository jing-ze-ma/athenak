#ifndef RECONSTRUCT_WB_PLM_BLEND_HPP_
#define RECONSTRUCT_WB_PLM_BLEND_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file wb_plm_blend.hpp
//! \brief wbplm-1103 <hydro|mhd>/wb_plm_smooth (default false = bitwise legacy switch):
//! a CONTINUOUS replacement for the fallback of the well-balanced x1 PLM
//! (Hydro/MHD::GridPiecewiseLinearX1 and GridPiecewiseLinearDerX1).
//!
//! The legacy rule takes the plain-PLM faces whenever a well-balanced face is negative or
//! the raw data has an extremum (sL*sR <= 0, unless plm_c1).  Both tests are threshold
//! flips: at an outflow1 top the ghost repeats the cell to a few ulp, sR is 0 +- roundoff,
//! and 1 ulp of input switches the top face between the two reconstructions (measured on
//! dhj_r1m v8: face pressure 1.7e-2 vs 7.3e-4, B 1e-16 -> 2e-6 in one cycle).
//!
//! The blend uses the faces q = w q_wb + (1 - w) q_plm with
//!   w = S(th / kWbThExt) * S(min(q_wb,L, q_wb,R) / (kWbPosFrac |q_i|)),
//!   th = 2 sL sR / (sL^2 + sR^2) in [-1, 1]   (the same smoothness indicator as plm_c1),
//!   S(x) = smoothstep, 0 for x <= 0 and 1 for x >= 1.
//! w is exactly 1 where the legacy rule kept the well-balanced faces with margin (the
//! raw slopes within a factor ~10 of each other, both faces above a tenth of the cell
//! value), so hydrostatic columns stay well balanced bitwise; it is exactly 0 wherever the
//! legacy rule fell back; in between it is C1 in the data.  The extremum factor is 1
//! under plm_c1, as in the legacy rule.  Both face sets are positive where they are used
//! (w > 0 needs both well-balanced faces > 0), so the blend is too.
#include <cmath>
#include "athena.hpp"

namespace wbplm {
constexpr Real kWbThExt = 0.2;     // th at which the extremum factor reaches 1
constexpr Real kWbPosFrac = 0.1;   // face / cell value at which the positivity factor is 1

KOKKOS_INLINE_FUNCTION
Real SmoothStep01(const Real x) {
  const Real t = fmin(fmax(x, static_cast<Real>(0.0)), static_cast<Real>(1.0));
  return t*t*(3.0 - 2.0*t);
}

//! weight of the well-balanced faces (qwl = q_{i+1/2,L}, qwr = q_{i-1/2,R})
KOKKOS_INLINE_FUNCTION
Real WbPlmWeight(const Real q_i, const Real sL, const Real sR, const Real qwl,
                 const Real qwr, const bool c1) {
  Real w = 1.0;
  if (!c1) {
    const Real den = sL*sL + sR*sR;
    const Real th = (den > 0.0) ? 2.0*sL*sR/den : 0.0;
    w = SmoothStep01(th/kWbThExt);
  }
  const Real qmin = fmin(qwl, qwr);
  const Real sc = kWbPosFrac*fabs(q_i);
  w *= (sc > 0.0) ? SmoothStep01(qmin/sc) : ((qmin > 0.0) ? 1.0 : 0.0);
  return (isfinite(w) ? w : static_cast<Real>(0.0));
}

//! in place: (ql, qr) well-balanced faces -> blended with the plain-PLM faces (pl, pr)
KOKKOS_INLINE_FUNCTION
void WbPlmBlend(const Real q_i, const Real sL, const Real sR, const Real pl,
                const Real pr, const bool c1, Real &ql, Real &qr) {
  const Real w = WbPlmWeight(q_i, sL, sR, ql, qr, c1);
  if (w >= 1.0) return;
  if (w <= 0.0) { ql = pl; qr = pr; return; }
  ql = w*ql + (1.0 - w)*pl;
  qr = w*qr + (1.0 - w)*pr;
}
}  // namespace wbplm
#endif  // RECONSTRUCT_WB_PLM_BLEND_HPP_
