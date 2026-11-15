#ifndef EOS_C2P_FLOORS_IDEAL_HPP_
#define EOS_C2P_FLOORS_IDEAL_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file c2p_floors_ideal.hpp
//! \brief ideal-gas energy floors of the hydro/MHD inversion on e = E - KE (- ME) with a
//! final kinetic energy: pressure floor e >= pfloor/(gamma - 1), temperature floor (tflr)
//! and entropy floor (sflr).  Used by SingleC2P_IdealHyd/MHD (own = !defer_cons_floors,
//! tflr = sflr = true) and on the cubed sphere by GnomonicRaiseVelFloors /
//! GnomonicRaiseVelMHDFloors (own = true, tflr = sflr = false).  With own, etot is
//! rebuilt as e + e_k (+ e_m with kMag).

#include <math.h>

#include "athena.hpp"
#include "eos/eos.hpp"

template <bool kMag>
KOKKOS_INLINE_FUNCTION
void IdealFloorsPost(const EOS_Data &eos, const Real d, const Real di, const bool own,
                     const bool tflr, const bool sflr, const Real e_k, const Real e_m,
                     Real &eint, Real &etot, bool &efloor_used, bool &tfloor_used) {
  Real efloor = eos.pfloor/(eos.gamma - 1.0);
  Real gm1 = eos.gamma - 1.0;
  if (eint < efloor) {
    eint = efloor;
    if (own) {
      if constexpr (kMag) {
        etot = efloor + e_k + e_m;
      } else {
        etot = efloor + e_k;
      }
    }
    efloor_used = true;
  }
  // apply temperature floor
  if (tflr) {
    const Real tfloor = eos.tfloor;
    if (gm1*eint*di < tfloor) {
      eint = d*tfloor/gm1;
      if (own) {
        if constexpr (kMag) {
          etot = eint + e_k + e_m;
        } else {
          etot = eint + e_k;
        }
      }
      tfloor_used = true;
    }
  }
  // apply entropy floor
  if (sflr) {
    const Real sfloor = eos.sfloor;
    Real spe_over_eps = gm1/pow(d, gm1);
    Real spe = spe_over_eps*eint*di;
    if (spe <= sfloor) {
      eint = d*sfloor/spe_over_eps;
      efloor_used = true;
    }
  }
}

#endif  // EOS_C2P_FLOORS_IDEAL_HPP_
