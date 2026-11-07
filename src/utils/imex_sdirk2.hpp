#ifndef UTILS_IMEX_SDIRK2_HPP_
#define UTILS_IMEX_SDIRK2_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file imex_sdirk2.hpp
//! \brief ckhesdirk-1102: the tableau of the IMEX second-order step shared by
//! <rad_m1>/time_scheme = hesdirk2 (rad_m1_time2.cpp) and <problem>/ck_time_scheme =
//! hesdirk2 (deep_hot_jupiter_rt_src.cpp, DhjCkT2Stage).  The explicit part is AthenaK's
//! rk2 (Heun); the implicit part the stiffly accurate, L-stable ESDIRK
//!   A = [[0,0,0],[1-g,g,0],[1/2,b2,g]],  g = 1 - 1/sqrt 2,  b2 = g/(2(1-g)),
//! FSAL.  Each implicit stage is a backward-Euler solve with dt -> g dt on the old vector
//!   (stage start) + c1 K1 + c2 K2,   stage 1: c1 = a21 dt,   stage 2: c1 = c31 dt,
//! c2 = c32 dt (the Heun average w0 U^n + w1 Y2 of the explicit part supplies the rest of
//! the row), and its slope is K = (Y - old vector)/(g dt).  trbdf2: the ESDIRK form of
//! TR-BDF2 with the RK2 of c2 = 2 - sqrt 2 (M1 only).  The expressions are those of
//! rad_m1_time2.cpp before the factorisation (bitwise).

#include <cmath>

#include "athena.hpp"

namespace imex_sdirk2 {

struct Tableau {
  Real g;
  Real a21;          // stage-1 increment coefficient of K1
  Real w0, w1;       // explicit stage-2 weights of U^n and Y2
  Real c31, c32;     // stage-2 increment coefficients of K1, K2
};

inline Real Gamma() {return 1.0 - 1.0/std::sqrt(2.0);}

// tab = 0: hesdirk2, 1: trbdf2, 2: sdirk2
inline Tableau Make(const int tab) {
  Tableau t;
  const Real g = Gamma();
  t.g = g;
  if (tab == 0) {
    const Real b2 = g/(2.0*(1.0 - g));
    t.a21 = 1.0 - g;
    t.w0 = 0.5;
    t.w1 = 0.5;
    t.c31 = 0.5*g;
    t.c32 = b2 - 0.5*g;
  } else if (tab == 2) {
    // seamE-1106 sdirk2: Alexander's two-stage L-stable SDIRK (g = 1 - 1/sqrt 2): no K1 in
    // either stage.  With the Heun average of the start state, stage 2's start is
    // U^n + (1 - g) dt K2 + (the implicit g dt K3): c32 = (1 - g) - g/2.  The hydro keeps
    // Heun (c2 = 1), so the implicit stage 1 at c = g is not the hydro's stage time.
    t.a21 = 0.0;
    t.w0 = 0.5;
    t.w1 = 0.5;
    t.c31 = 0.0;
    t.c32 = (1.0 - g) - 0.5*g;
} else {
    const Real w = 0.25*std::sqrt(2.0);
    t.a21 = g;
    t.w0 = 0.25;
    t.w1 = 0.75;
    t.c31 = w - 0.25*g;
    t.c32 = w - 0.25*g;
  }
  return t;
}

// the increment coefficients of stage `stage` (1 or 2) for the step dt
inline void IncCoeffs(const Tableau &t, const int stage, const Real dt, Real &c1,
                      Real &c2) {
  if (stage == 1) {
    c1 = t.a21*dt;
    c2 = 0.0;
  } else {
    c1 = t.c31*dt;
    c2 = t.c32*dt;
  }
}

}  // namespace imex_sdirk2

#endif  // UTILS_IMEX_SDIRK2_HPP_
