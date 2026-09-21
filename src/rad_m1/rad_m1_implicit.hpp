#ifndef RAD_M1_RAD_M1_IMPLICIT_HPP_
#define RAD_M1_RAD_M1_IMPLICIT_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_implicit.hpp
//! \brief constants and device-inline helpers of the IMPLICIT x1 column transport,
//! <rad_m1>/transport = implicit_x1 (docs/dev/rad_m1_implicit_design.md, milestone 3a).
//!
//! Nothing in here is reachable from the explicit path: the file is included only by
//! rad_m1_implicit.cpp and by the small delimited hooks of rad_m1.cpp / rad_m1_tasks.cpp
//! / rad_m1_newdt.cpp, so `transport = explicit` stays bitwise what it was.

#include <math.h>

#include "athena.hpp"

namespace radm1 {

// <rad_m1>/transport
constexpr int M1_TRANSPORT_EXPLICIT   = 0;   // stages 1-2: PD-ARS, sub-cycled
constexpr int M1_TRANSPORT_IMPLICIT_X1 = 1;  // stage 3a: backward Euler on x1 columns

// x1 boundary conditions of the implicit solve, in FACE-FLUX form (design sect. 3).
// The face value used is the TOTAL normal flux F at the boundary face.
constexpr int M1_IBC_MARSHAK  = 0;   // free surface, F_f = +- c*marshak_q*E_boundary
constexpr int M1_IBC_FLUX     = 1;   // imposed flux, F_f = implicit_flux_x1min/max
constexpr int M1_IBC_REFLECT  = 2;   // F_f = 0
constexpr int M1_IBC_PERIODIC = 3;   // cyclic tridiagonal (Sherman-Morrison)
constexpr int M1_IBC_EFIX     = 4;   // Dirichlet: the boundary CELL keeps the E it has
                                     // at the start of every step, i.e. it is frozen at
                                     // its initial value.  This is the face-flux form of
                                     // the fixed-E ghost the explicit tests T3b and T6
                                     // use, and the only boundary that anchors the LEVEL
                                     // of E in a pure-scattering column (imposed flux on
                                     // both ends leaves the operator singular).

// components of the implicit work array RadiationM1::iw (m,n,k,j,i).  All are per-cell
// and live only for the duration of one solve, except EN/EGN which carry the
// start-of-step state through the Picard loop.
constexpr int M1_IW_EN   = 0;   // E^n, the radiation energy at the start of the step
constexpr int M1_IW_EP   = 1;   // the current Picard iterate of E'
constexpr int M1_IW_TP   = 2;   // the current Picard iterate of T'
constexpr int M1_IW_EGN  = 3;   // rho e^n, the gas INTERNAL energy density at the start
constexpr int M1_IW_WCHI = 4;   // w = P_11/E of the lagged closure (= chi in 1-D)
constexpr int M1_IW_ADV  = 5;   // a = v1 (1 + w), so that A = v E + v.P = a E
constexpr int M1_IW_DE0  = 6;   // E0 - E, the O(beta) comoving correction (explicit)
constexpr int M1_IW_G0   = 7;   // g0 = rho(kappa_E E0 - kappa_P a T^4), lagged
constexpr int M1_IW_V1   = 8;   // the x1 gas velocity (0 without <hydro>)
constexpr int M1_IW_SRCB = 9;   // the emission/absorption addition to the DIAGONAL
constexpr int M1_IW_SRCR = 10;  // the emission/absorption addition to the RHS
constexpr int M1_IW_TA   = 11;  // tridiagonal lower diagonal
constexpr int M1_IW_TB   = 12;  // tridiagonal diagonal
constexpr int M1_IW_TC   = 13;  // tridiagonal upper diagonal
constexpr int M1_IW_TR   = 14;  // tridiagonal right-hand side
constexpr int M1_IW_S1   = 15;  // Thomas scratch (c')
constexpr int M1_IW_S2   = 16;  // Thomas solution
constexpr int M1_IW_S3   = 17;  // second Thomas solution (cyclic / Sherman-Morrison)
constexpr int M1_IW_RES  = 18;  // per-cell Picard residual
constexpr int M1_IW_F1   = 19;  // the derived cell-centred x1 flux of the iterate
constexpr int M1_NIW = 20;

// safeguarded root find for T' inside the Picard loop
constexpr int  M1_IMPL_TMAXIT = 100;
constexpr Real M1_IMPL_TRTOL  = 1.0e-12;

//----------------------------------------------------------------------------------------
//! \fn M1ImplTemperature
//! \brief solve the LOCAL gas energy equation of the implicit transport scheme,
//!
//!   y(T) = rho e(rho,T) + c dt rho kappa_P a T^4 - [ rho e^n + c dt rho kappa_E E0' ]
//!        = 0,
//!
//! for T at fixed E' (E0' = E' + de0), by the same safeguarded Newton (`rtsafe`) that
//! the explicit coupling uses -- bracket by walking, then a Newton step accepted only
//! when it stays inside the bracket AND shrinks it at least as fast as bisection.  y is
//! strictly increasing in T, so the bracket is unique.
//!
//! This is NOT the explicit path's function (there E'(T) is eliminated locally and the
//! invariant is the bracket variable), so it cannot call rad_m1_coupling.cpp's copy;
//! the ALGORITHM is the same and the explicit arithmetic is left untouched.
//!
//! Returns the number of function evaluations; `ok` is false if no bracket was found.

template <class EosT>
KOKKOS_INLINE_FUNCTION
int M1ImplTemperature(const EosT &eos, const Real dd, const Real tguess,
                      const Real egn, const Real cdtkp_a, const Real rhs_abs,
                      Real &tout, bool &ok) {
  // y(T) = ee(T) + cdtkp_a*T^4 - (egn + rhs_abs)
  Real target = egn + rhs_abs;
  Real tlo = tguess, thi = tguess;
  Real ylo, yhi;
  int nev = 0;
  {
    Real ee, pp, cr, ct, cv;
    eos.ThermoAt(dd, tguess, ee, pp, cr, ct, cv);
    Real t2 = tguess*tguess;
    Real y0 = ee + cdtkp_a*t2*t2 - target;
    ylo = y0;
    yhi = y0;
    ++nev;
  }
  ok = true;
  if (ylo > 0.0) {
    for (int it=0; it<80 && ylo > 0.0; ++it) {
      tlo *= 0.5;
      Real ee, pp, cr, ct, cv;
      eos.ThermoAt(dd, tlo, ee, pp, cr, ct, cv);
      Real t2 = tlo*tlo;
      ylo = ee + cdtkp_a*t2*t2 - target;
      ++nev;
    }
    ok = (ylo <= 0.0);
  } else if (yhi < 0.0) {
    for (int it=0; it<80 && yhi < 0.0; ++it) {
      thi *= 2.0;
      Real ee, pp, cr, ct, cv;
      eos.ThermoAt(dd, thi, ee, pp, cr, ct, cv);
      Real t2 = thi*thi;
      yhi = ee + cdtkp_a*t2*t2 - target;
      ++nev;
    }
    ok = (yhi >= 0.0);
  }
  if (!ok) {
    tout = tguess;
    return nev;
  }
  Real tp = 0.5*(tlo + thi);
  if (tlo == thi) {tp = tlo;}
  Real dxold = fabs(thi - tlo);
  Real dx = dxold;
  bool conv = (tlo == thi);
  for (int it=0; it<M1_IMPL_TMAXIT && !conv; ++it) {
    Real ee, pp, cr, ct, cv;
    eos.ThermoAt(dd, tp, ee, pp, cr, ct, cv);
    ++nev;
    Real t3 = tp*tp*tp;
    Real y = ee + cdtkp_a*tp*t3 - target;
    if (y > 0.0) {thi = tp;} else {tlo = tp;}
    Real dy = dd*cv + 4.0*cdtkp_a*t3;
    dxold = dx;
    Real tn = (dy > 0.0) ? (tp - y/dy) : tp;
    if (!(tn > tlo && tn < thi) || (fabs(2.0*y) > fabs(dxold*dy))) {
      tn = 0.5*(tlo + thi);
    }
    dx = fabs(tn - tp);
    conv = (dx <= M1_IMPL_TRTOL*fabs(tn)) || (fabs(thi - tlo) <= M1_IMPL_TRTOL*fabs(tp));
    tp = tn;
  }
  tout = tp;
  return nev;
}

} // namespace radm1
#endif // RAD_M1_RAD_M1_IMPLICIT_HPP_
