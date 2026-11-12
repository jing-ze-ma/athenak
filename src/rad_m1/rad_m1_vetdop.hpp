#ifndef RAD_M1_RAD_M1_VETDOP_HPP_
#define RAD_M1_RAD_M1_VETDOP_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_vetdop.hpp
//! \brief vetdop-1112 (<rad_m1>/vet_doppler = boost, read only when named, default off =
//! bitwise): the Eddington tensor of a STATIC formal solution (vet_gd / vet_col /
//! vet_sc: no Doppler shift or aberration) taken as the COMOVING-frame Eddington tensor
//! D0 = P0/E0 and carried to the lab frame by the exact Lorentz transformation of the
//! radiation stress tensor (design /viper/ptmp2/jinma/rc18_1110/VET_DOPPLER_1112.md).
//!
//! T^{mu nu} = [[E, F/c], [F/c, P]].  The comoving T0 follows from the lab one by the
//! boost L(-beta), the lab one from T0 by L(+beta), L^0_0 = gamma, L^0_i = L^i_0 =
//! gamma beta_i, L^i_j = delta_ij + gamma^2/(gamma + 1) beta_i beta_j.  The lab E and F
//! are the EVOLVED moments; the lab P is the unknown.  Fixed point: P_lab such that the
//! comoving P0 of (E, F, P_lab) equals D0 E0.  Iteration: (E, F, D E) -> L(-b) ->
//! (E0, F0) -> (E0, F0, D0 E0) -> L(+b) -> D = P'/E'.  The map contracts by O(beta) per
//! pass and the pass from D = D0 is already exact to O(beta), so nit passes leave
//! O(beta^(nit+1)).  D has trace 1 exactly (T is trace-free under the boost).
//! At beta = 0 D0 is returned unchanged (bitwise).

#include <cmath>

#include "athena.hpp"

namespace radm1 {

// T' = L T L^T for symmetric 4x4 T and L = L(b) (b = v/c, |b| < 1)
KOKKOS_INLINE_FUNCTION
void M1VdopBoost4(const Real b[3], const Real t[4][4], Real o[4][4]) {
  const Real b2 = b[0]*b[0] + b[1]*b[1] + b[2]*b[2];
  const Real g = 1.0/sqrt(1.0 - b2);
  const Real gq = g*g/(g + 1.0);   // (gamma - 1)/beta^2 without the cancellation
  Real l[4][4];
  l[0][0] = g;
  for (int a = 0; a < 3; ++a) {
    l[0][a+1] = g*b[a];
    l[a+1][0] = g*b[a];
    for (int c = 0; c < 3; ++c) {
      l[a+1][c+1] = ((a == c) ? 1.0 : 0.0) + gq*b[a]*b[c];
    }
  }
  Real lt[4][4];
  for (int p = 0; p < 4; ++p) {
    for (int q = 0; q < 4; ++q) {
      Real s = 0.0;
      for (int r = 0; r < 4; ++r) {s += l[p][r]*t[r][q];}
      lt[p][q] = s;
    }
  }
  for (int p = 0; p < 4; ++p) {
    for (int q = p; q < 4; ++q) {
      Real s = 0.0;
      for (int r = 0; r < 4; ++r) {s += lt[p][r]*l[q][r];}
      o[p][q] = s;
      o[q][p] = s;
    }
  }
}

//! d0: comoving D0 (11 22 33 12 13 23); b: v/c; e: lab E; h: lab F/c (E units);
//! nit: passes (>= 1); d: the lab D (same order).  Returns 0 = boosted, 1 = beta 0
//! (d = d0), 2 = fallback (E <= 0, |b| >= 1 or a non-positive E0 / E': d = d0).
//! fcap: the comoving F0 of each pass capped at |F0| <= E0 (the realisable branch).
KOKKOS_INLINE_FUNCTION
int M1VetDopplerBoost(const Real d0[6], const Real b[3], const Real e, const Real h[3],
                      const int nit, const bool fcap, Real d[6]) {
  for (int n = 0; n < 6; ++n) {d[n] = d0[n];}
  const Real b2 = b[0]*b[0] + b[1]*b[1] + b[2]*b[2];
  if (b2 == 0.0) {return 1;}
  if (!(e > 0.0) || !(b2 < 1.0)) {return 2;}
  const Real bm[3] = {-b[0], -b[1], -b[2]};
  constexpr int ia[6] = {0, 1, 2, 0, 0, 1};
  constexpr int ib[6] = {0, 1, 2, 1, 2, 2};
  Real dc[6];
  for (int n = 0; n < 6; ++n) {dc[n] = d0[n];}
  for (int it = 0; it < nit; ++it) {
    Real t[4][4], t0[4][4];
    t[0][0] = e;
    for (int a = 0; a < 3; ++a) {
      t[0][a+1] = h[a];
      t[a+1][0] = h[a];
    }
    for (int n = 0; n < 6; ++n) {
      t[ia[n]+1][ib[n]+1] = dc[n]*e;
      t[ib[n]+1][ia[n]+1] = dc[n]*e;
    }
    M1VdopBoost4(bm, t, t0);
    const Real e0 = t0[0][0];
    if (!(e0 > 0.0)) {return 2;}
    Real f0[3] = {t0[0][1], t0[0][2], t0[0][3]};
    if (fcap) {
      const Real fm = sqrt(f0[0]*f0[0] + f0[1]*f0[1] + f0[2]*f0[2]);
      if (fm > e0) {
        const Real s = e0/fm;
        for (int a = 0; a < 3; ++a) {f0[a] *= s;}
      }
    }
    Real s[4][4], o[4][4];
    s[0][0] = e0;
    for (int a = 0; a < 3; ++a) {
      s[0][a+1] = f0[a];
      s[a+1][0] = f0[a];
    }
    for (int n = 0; n < 6; ++n) {
      s[ia[n]+1][ib[n]+1] = d0[n]*e0;
      s[ib[n]+1][ia[n]+1] = d0[n]*e0;
    }
    M1VdopBoost4(b, s, o);
    if (!(o[0][0] > 0.0)) {return 2;}
    const Real ie = 1.0/o[0][0];
    for (int n = 0; n < 6; ++n) {dc[n] = o[ia[n]+1][ib[n]+1]*ie;}
  }
  for (int n = 0; n < 6; ++n) {d[n] = dc[n];}
  return 0;
}

//! vet_doppler = local: the exponent 3 + y of the lab Doppler factor
//! (gamma (1 - n.beta))^-(3 + y) applied to a static intensity I of the cell (S the
//! cell's source function, same units).  Exact for a ray through a uniform medium with
//! the cell's n.beta: I = S D^4 (1 - exp(-tau/D)), so dI/I = (3 + y) n.beta with
//! y = 1 - tau e^-tau/(1 - e^-tau) (0 thin, 1 thick); tau from I = S (1 - e^-tau),
//! i.e. q = I/S, y = 1 + (1 - q) ln(1 - q)/q (q >= 1 or S <= 0: y = 1, the ray sees
//! a brighter, thick emitter).  ym: 0 that estimate, 1 y = 1, 2 y = 0.
KOKKOS_INLINE_FUNCTION
Real VetDopplerLocalY(const Real iv, const Real s, const int ym) {
  if (ym == 1) {return 1.0;}
  if (ym == 2) {return 0.0;}
  if (!(s > 0.0)) {return 1.0;}
  const Real q = iv/s;
  if (!(q < 1.0)) {return 1.0;}
  if (q < 1.0e-6) {return 0.5*fmax(q, 0.0);}
  return 1.0 + (1.0 - q)*log1p(-q)/q;
}

}  // namespace radm1

#endif  // RAD_M1_RAD_M1_VETDOP_HPP_
