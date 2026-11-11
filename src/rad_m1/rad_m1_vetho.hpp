#ifndef RAD_M1_RAD_M1_VETHO_HPP_
#define RAD_M1_RAD_M1_VETHO_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_vetho.hpp
//! \brief scho-1111 (EXPLORATION, opt-in): higher-order monotone upwind-plane
//! interpolation of the vet_sc intensity and a BESSER source along the segment.
//!
//! <rad_m1>/vet_sc_interp = bezier2 | mono3 (read only when named; absent = the
//! bilinear stencil, bitwise).  The upwind intensity on the plane of cell centres is
//! interpolated along x2 (tensor product with x3 in 3-D) from the 4-point stencil
//! j-1, j, j+1, j+2 around the foot:
//!   bezier2: quadratic Bezier on [f0, f1] with the control point
//!            c = (f0 + d0/2 + f1 - d1/2)/2, d = central differences, clipped to
//!            [min(f0, f1), max(f0, f1)] (de la Cruz Rodriguez & Piskunov 2013;
//!            Hennicker et al. 2020): exact for quadratics, bounded by f0, f1.
//!   mono3:   cubic Hermite with the Fritsch-Butland (harmonic) slopes, zero at an
//!            extremum (PCHIP-type): monotone, bounded by f0, f1.
//! Both: I >= 0 and no new extrema inside the interval; a linear profile is exact.
//! chi and S on the upwind plane keep the vet_interp rule (lin / log); at a material
//! edge (vet_interp corner / auto) the per-point segment OUTCOMES are interpolated with
//! the same 1-D scheme.
//! <rad_m1>/vet_sc_src = besser (read only when named): the source along the segment is
//! the quadratic Bezier of BESSER (Stepan & Trujillo Bueno 2013) through S_u, S_0 and the
//! downwind S_d (log / lin stencil on the downwind plane), control point clipped to
//! [S_u, S_0] (and to [S_0, S_d] on the downwind side); used on the smooth (non-corner)
//! segments; corner segments and the last layer keep the linear-in-tau segment.

#include <cmath>

#include "athena.hpp"

// 1-D interpolation at fraction a in [0, 1) of [f0, f1], neighbours fm (-1) and f2 (+2)
KOKKOS_INLINE_FUNCTION
Real VetHo1(const int ho, const Real fm, const Real f0, const Real f1, const Real f2,
            const Real a) {
  const Real lo = fmin(f0, f1), hi = fmax(f0, f1);
  if (ho == 1) {
    const Real d0 = 0.5*(f1 - fm), d1 = 0.5*(f2 - f0);
    Real c = 0.5*((f0 + 0.5*d0) + (f1 - 0.5*d1));
    c = fmin(fmax(c, lo), hi);
    const Real b = 1.0 - a;
    return b*b*f0 + 2.0*a*b*c + a*a*f1;
  }
  const Real dm = f0 - fm, dc = f1 - f0, dp = f2 - f1;
  const Real d0 = (dm*dc > 0.0) ? 2.0*dm*dc/(dm + dc) : 0.0;
  const Real d1 = (dc*dp > 0.0) ? 2.0*dc*dp/(dc + dp) : 0.0;
  const Real a2 = a*a, a3 = a2*a;
  const Real v = (2.0*a3 - 3.0*a2 + 1.0)*f0 + (a3 - 2.0*a2 + a)*d0
                 + (3.0*a2 - 2.0*a3)*f1 + (a3 - a2)*d1;
  return fmin(fmax(v, lo), hi);   // round-off guard (the scheme is bounded)
}

// the quadratic-Bezier SC weights of one segment of optical depth t:
// I_0 = I_u e^-t + wu S_u + w0 S_0 + wc C
KOKKOS_INLINE_FUNCTION
void VetBesserW(const Real t, Real &ex, Real &wu, Real &w0, Real &wc) {
  ex = exp(-t);
  if (t < 1.0e-2) {
    wu = t*(1.0/3.0 - t*(0.25 - t*(0.1 - t/36.0)));
    w0 = t*(1.0/3.0 - t*(1.0/12.0 - t*(1.0/60.0 - t/360.0)));
    wc = t*(1.0/3.0 - t*(1.0/6.0 - t*(0.05 - t/90.0)));
  } else {
    const Real t2 = t*t;
    wu = (2.0 - ex*(t2 + 2.0*t + 2.0))/t2;
    w0 = 1.0 - 2.0*(ex + t - 1.0)/t2;
    wc = 2.0*(t - 2.0 + ex*(t + 2.0))/t2;
  }
}

// BESSER segment: upwind (iu, su, hu = tau of the segment), centre s0, downwind sd at
// optical depth hd beyond the centre (hd <= 0: no downwind point -> C on the upwind
// linear ramp, i.e. the linear-in-tau segment)
KOKKOS_INLINE_FUNCTION
Real VetSegBesser(const Real iu, const Real su, const Real s0, const Real sd,
                  const Real hu, const Real hd) {
  Real c = 0.5*(su + s0);
  if (hd > 0.0 && hu > 0.0) {
    const Real dm = (s0 - su)/hu, dp = (sd - s0)/hd;
    if (dm*dp > 0.0) {
      Real yp = (hu*dp + hd*dm)/(hu + hd);
      c = s0 - 0.5*hu*yp;
      const Real cd = s0 + 0.5*hd*yp;
      if (c < fmin(su, s0) || c > fmax(su, s0)) {
        c = su;
      } else if (cd < fmin(s0, sd) || cd > fmax(s0, sd)) {
        yp = (sd - s0)/(0.5*hd);
        c = s0 - 0.5*hu*yp;
        c = fmin(fmax(c, fmin(su, s0)), fmax(su, s0));
      }
    } else {
      c = s0;
    }
  }
  Real ex, wu, w0, wc;
  VetBesserW(hu, ex, wu, w0, wc);
  return iu*ex + wu*su + w0*s0 + wc*c;
}

#endif // RAD_M1_RAD_M1_VETHO_HPP_
