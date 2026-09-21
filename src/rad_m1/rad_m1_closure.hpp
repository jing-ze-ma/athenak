#ifndef RAD_M1_RAD_M1_CLOSURE_HPP_
#define RAD_M1_RAD_M1_CLOSURE_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_closure.hpp
//! \brief Levermore (1984) M1 closure, the flux limiter and the Skinner & Ostriker
//! (2013, eq. 41a) closed-form wave speeds, as device-inline functions.
//!
//! Conventions (docs/dev/rad_m1_design.md sect. 1 and 3): E is the lab-frame radiation
//! energy density, F_i the lab-frame flux, f_i = F_i/(c E) the reduced flux with
//! |f| <= 1, n_i = f_i/|f| and
//!
//!   chi  = (3 + 4 f^2)/(5 + 2 sqrt(4 - 3 f^2))
//!   P_ij = E [ (1-chi)/2 delta_ij + (3 chi - 1)/2 n_i n_j ]
//!
//! All of this is algebraic in flat space: there is no root find anywhere in here.

#include <math.h>
#include "athena.hpp"

namespace radm1 {

// <rad_m1>/thick_flux: which asymptotic-preserving correction the face flux carries
// (design sect. 3).  Device code branches on these, so they are plain ints.
constexpr int M1_THICK_NONE   = 0;   // plain HLL
constexpr int M1_THICK_APHLL  = 1;   // Bloch et al. (2021) alpha on the E-flux only
constexpr int M1_THICK_SCALED = 2;   // Jiang (2021) App. A: both wave speeds x eps(tau)

//----------------------------------------------------------------------------------------
//! \fn M1Chi
//! \brief Levermore closure factor chi as a function of the reduced flux magnitude.
//! chi(0) = 1/3 (isotropic), chi(1) = 1 (free streaming).

KOKKOS_INLINE_FUNCTION
Real M1Chi(const Real fnorm) {
  Real f2 = fnorm*fnorm;
  if (f2 > 1.0) {f2 = 1.0;}
  Real s = sqrt(fmax(4.0 - 3.0*f2, 0.0));
  return (3.0 + 4.0*f2)/(5.0 + 2.0*s);
}

//----------------------------------------------------------------------------------------
//! \fn M1ApplyLimits
//! \brief Apply the energy floor and the flux limit |F| <= c E to one cell state.  The
//! flux is scaled by (c E/|F|), NOT by the square of that ratio.  Returns true if either
//! limit was active.

KOKKOS_INLINE_FUNCTION
bool M1ApplyLimits(const Real cl, const Real e_floor,
                   Real &e, Real &f1, Real &f2, Real &f3) {
  bool limited = false;
  if (!(e > e_floor)) {
    e = e_floor;
    limited = true;
  }
  Real fmag = sqrt(f1*f1 + f2*f2 + f3*f3);
  Real fmax = cl*e;
  if (fmag > fmax) {
    Real scale = (fmag > 0.0) ? (fmax/fmag) : 0.0;
    f1 *= scale;
    f2 *= scale;
    f3 *= scale;
    limited = true;
  }
  return limited;
}

//----------------------------------------------------------------------------------------
//! \fn M1ReducedFlux
//! \brief Given (E, F_i) return the reduced flux f_i = F_i/(c E) clipped to |f| <= 1,
//! and its magnitude.  E is assumed already floored.

KOKKOS_INLINE_FUNCTION
void M1ReducedFlux(const Real cl, const Real e, const Real ff1, const Real ff2,
                   const Real ff3, Real &r1, Real &r2, Real &r3, Real &rnorm) {
  Real inv = 1.0/(cl*e);
  r1 = ff1*inv;
  r2 = ff2*inv;
  r3 = ff3*inv;
  rnorm = sqrt(r1*r1 + r2*r2 + r3*r3);
  if (rnorm > 1.0) {
    Real scale = 1.0/rnorm;
    r1 *= scale;
    r2 *= scale;
    r3 *= scale;
    rnorm = 1.0;
  }
}

//----------------------------------------------------------------------------------------
//! \fn M1PressureCol
//! \brief One column of the Eddington tensor: P_{i,d} for i = 1,2,3 at fixed direction
//! index d (1,2,3), given E, the reduced flux (r1,r2,r3), its magnitude and chi.

KOKKOS_INLINE_FUNCTION
void M1PressureCol(const int d, const Real e, const Real r1, const Real r2,
                   const Real r3, const Real rnorm, const Real chi,
                   Real &p1d, Real &p2d, Real &p3d) {
  Real diag = 0.5*(1.0 - chi)*e;
  Real aniso = 0.5*(3.0*chi - 1.0)*e;
  // n_i n_d with n = f/|f|.  NOTE the unit vector is formed FIRST: writing this as
  // r_i r_d/|f|^2 makes a denormal |f| underflow the numerator while 1/|f|^2 overflows,
  // and 0*inf is a NaN (seen in the vacuum cells of the beam test).
  Real inv = (rnorm > 1.0e-100) ? (1.0/rnorm) : 0.0;
  Real n1 = r1*inv, n2 = r2*inv, n3 = r3*inv;
  Real nd = (d == 1) ? n1 : ((d == 2) ? n2 : n3);
  p1d = aniso*n1*nd;
  p2d = aniso*n2*nd;
  p3d = aniso*n3*nd;
  if (d == 1) {p1d += diag;} else if (d == 2) {p2d += diag;} else {p3d += diag;}
}

//----------------------------------------------------------------------------------------
//! \fn M1WaveSpeeds
//! \brief Skinner & Ostriker (2013) eq. 41a closed-form eigenvalues of the M1 system,
//! in units of the true speed of light:
//!
//!   lam_pm/c = { mu f +- sqrt[ (2/3)(4 - 3 f^2 - s) + 2 mu^2 (2 - f^2 - s) ] }/s,
//!   s = sqrt(4 - 3 f^2),   mu = cos(angle between F and the face normal).
//!
//! Checks: f = 0 gives +-1/sqrt(3); f = 1, mu = 1 gives +1 (degenerate).

KOKKOS_INLINE_FUNCTION
void M1WaveSpeeds(const Real fnorm, const Real mu, Real &lam_m, Real &lam_p) {
  Real f2 = fnorm*fnorm;
  if (f2 > 1.0) {f2 = 1.0;}
  Real s = sqrt(fmax(4.0 - 3.0*f2, 0.0));
  if (!(s > 1.0e-12)) {
    // f -> 1 in the face-normal direction: the two waves merge at mu c
    lam_m = mu;
    lam_p = mu;
    return;
  }
  Real rad = (2.0/3.0)*(4.0 - 3.0*f2 - s) + 2.0*mu*mu*(2.0 - f2 - s);
  rad = sqrt(fmax(rad, 0.0));
  lam_m = (mu*fnorm - rad)/s;
  lam_p = (mu*fnorm + rad)/s;
  // the M1 system is hyperbolic with |lam| <= c; guard against roundoff
  lam_m = fmax(lam_m, -1.0);
  lam_p = fmin(lam_p,  1.0);
}

//----------------------------------------------------------------------------------------
//! \fn M1HLLFlux
//! \brief HLL flux of (E, F_1, F_2, F_3) through a face whose normal is the coordinate
//! direction ivx (1, 2 or 3), from the reconstructed left/right states.  The evolved
//! system (design sect. 1) is
//!
//!   dE/dt   + (chat/c) d_d F_d     = 0
//!   dF_i/dt + chat c   d_d P_id    = 0
//!
//! which is the true M1 system with every characteristic speed multiplied by chat/c;
//! hence the eigenvalues above are used scaled by chat.  thick_flux = none: this is
//! plain HLL, with b_L = min(0, lam_-) and b_R = max(0, lam_+) over BOTH states.
//!
//! THICK LIMIT (design sect. 3), selected by `thick` and supplied with the face optical
//! depth tau_face = (1/2)[(rho kappa_F+rho kappa_s)_L + (..)_R] dx of THIS direction:
//!
//!  * M1_THICK_APHLL: alpha = 1/[1 - 3 tau_face (1-f^2) lam_+ lam_-/(c(lam_+ - lam_-))]
//!    (Bloch et al. 2021 eq. 25, f the face mean of |f|) blends the HLL E-flux with the
//!    compact two-point diffusion flux,
//!
//!      F_E = alpha F_E^HLL(reconstructed states) + (1 - alpha) F_diff,
//!      F_diff = -(chat/c) c (E_{i+1} - E_i)/(3 tau_face)    [CELL-CENTRE E, not the
//!                                                            reconstructed face states]
//!
//!    which is what survives a SECOND-ORDER reconstruction: with PLM the face jump
//!    E_R - E_L is O(dx^2) in smooth regions, so the dissipation term of Berthon's
//!    original AP-HLL flux -- which is what carries the physical diffusion there --
//!    vanishes and that scheme grossly UNDER-diffuses.  The
//!    blend is NOT the same scheme as Berthon's: with piecewise-constant states
//!    alpha F^HLL is ALREADY exactly the physical flux (the identity of design sect. 3),
//!    so adding (1-alpha) F_diff on top double counts and the diffusion comes out a
//!    factor 1 + (sqrt3/2)tau/(1 + (sqrt3/2)tau) -> 2 too large.  The two forms are
//!    therefore selected by the reconstruction: `berthon` (set by the caller for
//!    reconstruct = dc) is alpha F^HLL alone, the blend is used with PLM.  The
//!    diffusion limit of the blend is dE/dt = div[(chat/(3 rho kappa)) grad E], i.e. the
//!    physical diffusivity c/(3 rho kappa) carrying the (chat/c) of the evolved E
//!    equation.  F_diff uses the Eddington value 1/3 deliberately: the (1-f^2) guard and
//!    a small tau_face both send alpha -> 1 wherever the closure departs from 1/3.
//!    The F-flux (the pressure flux) stays plain HLL with the reconstructed states,
//!    which is what keeps -grad P_rad on the gas (sect. 3).
//!  * M1_THICK_SCALED: both HLL wave speeds are multiplied by
//!    eps = sqrt[(1 - exp(-tau_c^2))/tau_c^2], tau_c = scaled_pref*tau_face, in BOTH
//!    fluxes (Jiang 2021 App. A in moment form).
//!
//! MILESTONE 1c HOOK: the advective enthalpy-flux split for moving media replaces the
//! per-side normal flux entering the E-flux by the comoving F0_n = F_n - (v E + v.P)_n
//! and adds the upwinded full enthalpy flux A afterwards; the two marked locals f0ln /
//! f0rn and the separate assembly of flx[0] are there so that lands as a local change.

KOKKOS_INLINE_FUNCTION
void M1HLLFlux(const int ivx, const Real cl, const Real chat, const bool eddington,
               const Real el, const Real fl1, const Real fl2, const Real fl3,
               const Real er, const Real fr1, const Real fr2, const Real fr3,
               const int thick, const Real tau_face, const Real scaled_pref,
               const bool berthon, const Real ecl, const Real ecr, Real *flx) {
  // reduced fluxes and closure factors
  Real rl1, rl2, rl3, rlnorm, rr1, rr2, rr3, rrnorm;
  M1ReducedFlux(cl, el, fl1, fl2, fl3, rl1, rl2, rl3, rlnorm);
  M1ReducedFlux(cl, er, fr1, fr2, fr3, rr1, rr2, rr3, rrnorm);
  Real chil = eddington ? (1.0/3.0) : M1Chi(rlnorm);
  Real chir = eddington ? (1.0/3.0) : M1Chi(rrnorm);

  // Eddington-tensor column in the face-normal direction
  Real pl1, pl2, pl3, pr1, pr2, pr3;
  M1PressureCol(ivx, el, rl1, rl2, rl3, rlnorm, chil, pl1, pl2, pl3);
  M1PressureCol(ivx, er, rr1, rr2, rr3, rrnorm, chir, pr1, pr2, pr3);

  // physical fluxes
  Real fln = (ivx == 1) ? fl1 : ((ivx == 2) ? fl2 : fl3);
  Real frn = (ivx == 1) ? fr1 : ((ivx == 2) ? fr2 : fr3);
  Real fxl[4], fxr[4];
  fxl[0] = (chat/cl)*fln;
  fxr[0] = (chat/cl)*frn;
  fxl[1] = chat*cl*pl1;
  fxl[2] = chat*cl*pl2;
  fxl[3] = chat*cl*pl3;
  fxr[1] = chat*cl*pr1;
  fxr[2] = chat*cl*pr2;
  fxr[3] = chat*cl*pr3;

  // signal speeds
  Real bl, br;
  if (eddington) {
    br = chat/sqrt(3.0);
    bl = -br;
  } else {
    Real mul = (rlnorm > 0.0) ?
               (((ivx == 1) ? rl1 : ((ivx == 2) ? rl2 : rl3))/rlnorm) : 0.0;
    Real mur = (rrnorm > 0.0) ?
               (((ivx == 1) ? rr1 : ((ivx == 2) ? rr2 : rr3))/rrnorm) : 0.0;
    Real lml, lpl, lmr, lpr;
    M1WaveSpeeds(rlnorm, mul, lml, lpl);
    M1WaveSpeeds(rrnorm, mur, lmr, lpr);
    bl = chat*fmin(fmin(lml, lmr), 0.0);
    br = chat*fmax(fmax(lpl, lpr), 0.0);
  }

  // thick-limit corrections.  alpha = 1 and eps = 1 leave the plain HLL flux, and the
  // arithmetic below is then bit-identical to milestone 1a.
  Real alpha = 1.0;
  if (thick == M1_THICK_SCALED && tau_face > 0.0) {
    Real tc = scaled_pref*tau_face;
    // eps = sqrt[(1 - exp(-tc^2))/tc^2]; the series is used where the ratio is 0/0
    Real eps = 1.0;
    if (tc > 1.0e-4) {
      eps = sqrt(-expm1(-tc*tc))/tc;
    } else {
      eps = sqrt(fmax(1.0 - 0.5*tc*tc, 0.0));
    }
    bl *= eps;
    br *= eps;
  } else if (thick == M1_THICK_APHLL && tau_face > 0.0) {
    Real fbar = 0.5*(rlnorm + rrnorm);
    Real guard = fmax(1.0 - fbar*fbar, 0.0);
    // lam_+ lam_-/(lam_+ - lam_-) in units of the true light speed
    Real lp = br/chat;
    Real lm = bl/chat;
    Real den = 1.0 - 3.0*tau_face*guard*lp*lm/(lp - lm + 1.0e-300);
    alpha = 1.0/fmax(den, 1.0);
  }

  // HLL
  Real ul[4], ur[4];
  ul[0] = el;
  ul[1] = fl1;
  ul[2] = fl2;
  ul[3] = fl3;
  ur[0] = er;
  ur[1] = fr1;
  ur[2] = fr2;
  ur[3] = fr3;
  Real invb = 1.0/(br - bl + 1.0e-300);
  for (int n=1; n<4; ++n) {
    flx[n] = (br*fxl[n] - bl*fxr[n] + br*bl*(ur[n] - ul[n]))*invb;
  }
  // the E-flux, assembled separately: it is the one the thick-limit blend touches, and
  // the one milestone 1c splits into a comoving part plus an upwinded enthalpy flux A.
  Real f0ln = fxl[0];   // 1c: (chat/c)(F_n - (v E + v.P)_n) of the LEFT state
  Real f0rn = fxr[0];   // 1c: the same of the RIGHT state
  Real fe = (br*f0ln - bl*f0rn + br*bl*(ur[0] - ul[0]))*invb;
  if (thick == M1_THICK_APHLL && tau_face > 0.0 && alpha < 1.0) {
    if (berthon) {
      fe *= alpha;
    } else {
      // (chat/c)*[-c (E_R - E_L)/(3 tau_face)] = -chat (E_R - E_L)/(3 tau_face)
      Real fdiff = -chat*(ecr - ecl)/(3.0*tau_face);
      Real fadv = (br*f0ln - bl*f0rn)*invb;
      Real fdis = br*bl*(ur[0] - ul[0])*invb;
      // The HLL DISSIPATION carries weight alpha^2, not alpha.  It is spurious once
      // F_diff supplies the physical diffusion, and it has to vanish FASTER than that
      // flux does: D_HLL ~ c dE while F_diff ~ c dE/tau, so a weight alpha ~ 1/tau
      // leaves a residual of the same order as the physical flux.  That is invisible on
      // a smooth profile, where the reconstructed jump dE is O(dx^2), but at a KINK --
      // an opacity jump (T3b), a limiter-clipped extremum (T3c) -- dE is O(dx) and the
      // residual is O(F): 65 % of the flux at the 1e3 jump of T3b, measured, before
      // this was squared.  alpha^2 ~ 1/tau^2 removes it and leaves the thin limit
      // (alpha -> 1) untouched.
      fe = alpha*fadv + alpha*alpha*fdis + (1.0 - alpha)*fdiff;
    }
  }
  flx[0] = fe;          // 1c: + A_upwind
}

} // namespace radm1
#endif // RAD_M1_RAD_M1_CLOSURE_HPP_
