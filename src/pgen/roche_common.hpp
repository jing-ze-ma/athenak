#ifndef PGEN_ROCHE_COMMON_HPP_
#define PGEN_ROCHE_COMMON_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file roche_common.hpp
//! \brief Roche (restricted three-body) potential helpers in a frame CENTRED ON BODY 1
//! (the body the grid is centred on) and corotating with the circular orbit, shared by
//! problem generators (roche-dhj-1031: deep_hot_jupiter_rt.cpp, problem/roche).
//!
//! FRAME.  Origin at body 1 (mass M1, G M1 = gm1), companion (G M2 = gm2) at distance a
//! along the unit vector e_s, z = orbital (= spin, synchronous) axis, frame angular
//! velocity Omega about the axis through the centre of mass x_cm = a M2/(M1 + M2) e_s.
//! A pure translation of the CM frame to body 1, so
//!     Phi = -G M1/r - G M2/|x - a e_s| - Omega^2 [(x - x_cm)_perp^2]/2.
//! Expanding the centrifugal term about body 1:
//!     -Omega^2 R^2/2 + Omega^2 x_cm (x.e_s) - Omega^2 x_cm^2/2,
//! and Omega^2 x_cm = G M2/a^2 (Kepler, Omega^2 a^3 = G (M1 + M2)) is the INDIRECT term
//! (minus body 1's acceleration toward the companion).  Written ADDITIVELY here:
//!     Phi = Phi_1(r) - Omega^2 R^2/2 + TidePot,
//!     TidePot = -G M2/d + G M2/a + (G M2/a^2) s,   s = x.e_s, d = |x - a e_s|,
//! with the constant chosen so TidePot(0) = 0.  TidePot = -(G M2/a^3) r^2 (3 mu^2 - 1)/2
//! + O(r^3/a^4) (mu = s/r): the tidal tensor; gm2 = 0 returns exactly 0, so the caller's
//! Phi_1 - Omega^2 R^2/2 (deep_hot_jupiter_rt TotPotAt, problem/rot_potential) is
//! reproduced BITWISE.  The indirect term uses G M2/a^2 itself, so the frame is
//! consistent for any gm2; the CM rotation is consistent only for the Keplerian Omega
//! (the caller derives gm2 from Omega and a by default).
//!
//! The accretor (ry_per_accretor.cpp RocheParams/RochePotXYZ/RocheGeffXYZ, centrifugal
//! about x_cm written directly) is the same potential up to the constant; it can move
//! to these helpers once csacc-r-1031 has landed (kept unchanged there for now).

#include <math.h>

#include "athena.hpp"
#include "coordinates/cubed_sphere.hpp"

namespace roche {

//! the tidal (companion + indirect) potential at distance r from body 1, mu = cosine of
//! the angle to the companion direction; 0 at the origin, exactly 0 for gm2 = 0
KOKKOS_INLINE_FUNCTION
Real TidePot(const Real gm2, const Real a, const Real r, const Real mu) {
  if (gm2 == 0.0) return 0.0;
  const Real d = sqrt(fmax(r*r + a*a - 2.0*a*r*mu, 1.0e-30*a*a));
  return -gm2/d + gm2/a + gm2*r*mu/(a*a);
}

//! grad Phi (Cartesian, s along the companion, z along the spin axis, y completing a
//! right-handed or left-handed triple: only |y| enters) and |grad Phi| = the effective
//! gravity, for Phi = -gm1/r - om2 (s^2 + y^2)/2 + TidePot.
KOKKOS_INLINE_FUNCTION
Real GeffPlanet(const Real gm1, const Real om2, const Real gm2, const Real a,
                const Real s, const Real y, const Real z,
                Real *gs = nullptr, Real *gy = nullptr, Real *gz = nullptr) {
  const Real r = sqrt(fmax(s*s + y*y + z*z, 1.0e-300));
  const Real ir3 = gm1/(r*r*r);
  Real ds = ir3*s - om2*s, dy = ir3*y - om2*y, dz = ir3*z;
  if (gm2 != 0.0) {
    const Real d = sqrt(fmax(SQR(s - a) + y*y + z*z, 1.0e-30*a*a));
    const Real id3 = gm2/(d*d*d);
    ds += id3*(s - a) + gm2/(a*a);
    dy += id3*y;
    dz += id3*z;
  }
  if (gs != nullptr) {*gs = ds; *gy = dy; *gz = dz;}
  return sqrt(ds*ds + dy*dy + dz*dz);
}

//! grad TidePot in the same (s, y, z) frame (s toward the companion); all zero for
//! gm2 = 0.  The tidal ACCELERATION is minus this.
KOKKOS_INLINE_FUNCTION
void TideGrad(const Real gm2, const Real a, const Real s, const Real y, const Real z,
              Real &gs, Real &gy, Real &gz) {
  if (gm2 == 0.0) {gs = 0.0; gy = 0.0; gz = 0.0; return;}
  const Real d = sqrt(fmax(SQR(s - a) + y*y + z*z, 1.0e-30*a*a));
  const Real id3 = gm2/(d*d*d);
  gs = id3*(s - a) + gm2/(a*a);
  gy = id3*y;
  gz = id3*z;
}

//! the arc (radians, on the unit sphere) between the cubed-sphere panel points (a2, a3)
//! and (b2, b3), code coordinates x2, x3 in [-1, 1] (angle pi/4 x).  Same as the
//! accretor's CsArc (csacc-1031).
KOKKOS_INLINE_FUNCTION
Real CsArc(const int pnl, const Real a2, const Real a3, const Real b2, const Real b3) {
  Real qa[3], qb[3];
  cubed_sphere::PanelToCart(pnl, 0.25*M_PI*a2, 0.25*M_PI*a3, qa);
  cubed_sphere::PanelToCart(pnl, 0.25*M_PI*b2, 0.25*M_PI*b3, qb);
  const Real cx = qa[1]*qb[2] - qa[2]*qb[1], cy = qa[2]*qb[0] - qa[0]*qb[2];
  const Real cz = qa[0]*qb[1] - qa[1]*qb[0];
  return atan2(sqrt(cx*cx + cy*cy + cz*cz), qa[0]*qb[0] + qa[1]*qb[1] + qa[2]*qb[2]);
}

//! the L1 point between body 1 and the companion on the line of centres: distance from
//! body 1 (Newton on dPhi/ds = 0 along s, host only).  gm1, gm2, a, om2 as above.
inline Real L1Distance(const Real gm1, const Real gm2, const Real a, const Real om2) {
  // Hill guess
  Real s = a*std::cbrt(gm1/(3.0*(gm1 + gm2)));
  for (int it=0; it<100; ++it) {
    // dPhi/ds on the axis, 0 < s < a
    const Real f = gm1/(s*s) - gm2/SQR(a - s) - om2*s + gm2/(a*a);
    const Real df = -2.0*gm1/(s*s*s) - 2.0*gm2/((a - s)*(a - s)*(a - s)) - om2;
    const Real ds = f/df;
    s -= ds;
    if (std::fabs(ds) < 1.0e-14*s) break;
  }
  return s;
}

}  // namespace roche

#endif  // PGEN_ROCHE_COMMON_HPP_
