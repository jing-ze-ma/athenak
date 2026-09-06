#ifndef UTILS_ROSSELAND_HPP_
#define UTILS_ROSSELAND_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rosseland.hpp
//! \brief the Freedman et al. (2014, ApJS 214, 25) analytic fit to the Rosseland mean
//! gas opacity, kappa_R(T, p, [M/H]) in cm^2/g, T in K, p in dyn/cm^2.  Their eqs. 4-5
//! with the Table 2 coefficients; clamped to the fit's range (75-4000 K, 1-3e8 dyn/cm^2)
//! outside it.  Shared by the deep_hot_jupiter_rt problem generator (its grey picket
//! fence) and the radiative-conductivity branch of the conduction module.

#include <math.h>

#include "athena.hpp"

KOKKOS_INLINE_FUNCTION
Real RosselandFreedman2014(const Real T, const Real p, const Real met) {
  Real T1 = fmin(fmax(T, 75.0), 4000.0);
  Real p1 = fmin(fmax(p, 1.0), 3.0e8);
  const Real lgT = log10(T1);
  const Real lgp = log10(p1);
  const Real c1 = 10.602, c2 = 2.882, c3 = 6.09e-15, c4 = 2.954, c5 = -2.526;
  const Real c6 = 0.843, c7 = -5.490;
  Real c8, c9, c10, c11, c12, c13;
  if (T1 < 800.0) {
    c8 = -14.051; c9 = 3.055; c10 = 0.024; c11 = 1.877; c12 = -0.445; c13 = 0.8321;
  } else {
    c8 = 82.241; c9 = -55.456; c10 = 8.754; c11 = 0.7048; c12 = -0.0414; c13 = 0.8321;
  }
  const Real lgkl = c1*atan(lgT-c2) - c3/(lgp+c4)*exp(SQR(lgT-c5)) + c6*met + c7;
  const Real lgkh = c8 + c9*lgT + c10*SQR(lgT) + lgp*(c11+c12*lgT)
                  + c13*met*(0.5+1.0/M_PI*atan((lgT-2.5)/0.2));
  return pow(10.0,lgkl) + pow(10.0,lgkh);
}

#endif // UTILS_ROSSELAND_HPP_
