#ifndef UTILS_IMEX_SDIRK2_HPP_
#define UTILS_IMEX_SDIRK2_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file imex_sdirk2.hpp
//! \brief The tableau of the IMEX second-order step (rad_vet_imex.cpp): the explicit part
//! is AthenaK's rk2 (Heun), the implicit part Alexander's two-stage L-stable SDIRK2,
//! g = 1 - 1/sqrt 2.  Each implicit stage is a backward-Euler solve with dt -> g dt on
//! the old vector (stage start) + c1 K1 + c2 K2: stage 1 c1 = a21 dt = 0 (no K1);
//! stage 2 starts from the Heun average w0 U^n + w1 Y2 with c2 = c32 dt, c32 =
//! (1 - g) - g/2.  Its slope is K = (Y - old vector)/(g dt).  The gas of the stage-1
//! solve is at t^n + g dt (Time2GasTime).

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

inline Tableau Make() {
  Tableau t;
  const Real g = Gamma();
  t.g = g;
  t.a21 = 0.0;
  t.w0 = 0.5;
  t.w1 = 0.5;
  t.c31 = 0.0;
  t.c32 = (1.0 - g) - 0.5*g;
  return t;
}

}  // namespace imex_sdirk2

#endif  // UTILS_IMEX_SDIRK2_HPP_
