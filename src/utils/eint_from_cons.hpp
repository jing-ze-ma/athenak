#ifndef UTILS_EINT_FROM_CONS_HPP_
#define UTILS_EINT_FROM_CONS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file eint_from_cons.hpp
//! \brief the internal energy density of a cell, taken from the CONSERVED state.
//!
//! WHY THIS EXISTS AS A SHARED FUNCTION.  w0(IEN) is whatever the last ConToPrim wrote,
//! i.e. the state at the END of the previous stage.  Every operator that runs inside a
//! stage -- RKUpdate, the source terms, the implicit radial conduction, the two-stream
//! apply -- sees a u0 that has already moved, and in a smooth cell u0 and w0 agree to
//! O(dt) so the difference never mattered.  In a runaway cell it does: the hydro flux
//! divergence can change e by ~100 % in one stage, and an operator that evaluates its
//! own limiter, positivity guard or equilibrium on the stale w0 is working from a number
//! that no longer exists.  That is how the red-giant photospheric loop opens -- the RT
//! removes more energy than the cell has, the floor creates energy repairing it, and the
//! next Riemann solve sees a 0.1 K cell against a 1e4 K neighbour.
//!
//! The extraction is exactly the one ConToPrim performs, in the same order, so that a
//! caller which switches from w0 to this expression changes nothing in a smooth cell:
//!   * the gravitational term rho*Phi, when <hydro>/etotgrav is on (see AddGravEtot);
//!   * the kinetic energy formed with the NON-ORTHOGONAL cubed-sphere metric, i.e. with
//!     the covariant momentum properly raised (see GnomonicEquiangleRaiseVel) -- summing
//!     the squares of the stored components instead is wrong by the cross term and the
//!     error is largest exactly at the panel corners, where these cells live.
//! Callers: Conduction::ImplicitRadialUpdate and, under <problem>/rt_use_cons, the
//! two-stream apply and its temperature/opacity precompute.

#include "athena.hpp"

//----------------------------------------------------------------------------------------
//! \fn Real EintFromCons
//! \brief e = u0(IEN) - rho*Phi - KE, with KE on the gnomonic metric when cs is set.
//! `cosc` is cos_cell(m,k,j) (ignored unless cs); `phi` is phicc0(m,k,j,i) (ignored
//! unless etg).  Returns whatever the arithmetic gives, including a non-positive value:
//! deciding what to do about that belongs to the caller.

KOKKOS_INLINE_FUNCTION
Real EintFromCons(const DvceArray5D<Real> &u0, const int m, const int k, const int j,
                  const int i, const Real cosc, const bool cs, const bool etg,
                  const Real phi) {
  const Real d = u0(m,IDN,k,j,i);
  Real ekin;
  if (cs) {
    const Real c = cosc;
    const Real det = 1.0 - c*c;
    const Real q1 = u0(m,IM1,k,j,i);
    const Real q2 = u0(m,IM2,k,j,i);
    const Real q3 = u0(m,IM3,k,j,i);
    ekin = 0.5*(q1*(q1/d) + q2*((q2 - c*q3)/(d*det)) + q3*((q3 - c*q2)/(d*det)));
  } else {
    ekin = 0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i))
                + SQR(u0(m,IM3,k,j,i)))/d;
  }
  Real ei = u0(m,IEN,k,j,i) - ekin;
  if (etg) ei -= d*phi;
  return ei;
}

#endif  // UTILS_EINT_FROM_CONS_HPP_
