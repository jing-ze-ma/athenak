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
//!   * in MHD, the MAGNETIC energy: u0(IEN) holds E = e + KE + ME, so an extraction
//!     that stops at KE returns e + ME and is wrong by the whole field energy -- at
//!     plasma beta ~ 1 that is a factor of two in the temperature the RT and the
//!     conduction then work from.  See MagEnergyCC below for which field to use.
//! Callers: Conduction::ImplicitRadialUpdate and, under <problem>/rt_use_cons, the
//! two-stream apply and its temperature/opacity precompute.

#include "athena.hpp"

//----------------------------------------------------------------------------------------
//! \fn Real EintFromCons
//! \brief e = u0(IEN) - rho*Phi - KE, with KE on the gnomonic metric when cs is set.
//! `cosc` is cos_cell(m,k,j) (ignored unless cs); `phi` is phicc0(m,k,j,i) (ignored
//! unless etg).  Returns whatever the arithmetic gives, including a non-positive value:
//! deciding what to do about that belongs to the caller.

//----------------------------------------------------------------------------------------
//! \fn Real MagEnergyCC
//! \brief 0.5|B|^2 from the CELL-CENTRED field bcc0, on every grid.
//!
//! WHY THE PLAIN SUM OF SQUARES IS RIGHT HERE, INCLUDING ON THE CUBED SPHERE.  b0.x*f
//! holds the flux density through its own face, so the triple ConsToPrim averages is
//! (B.rhat, B.nhat_xi, B.nhat_eta), and those last two directions are NOT orthogonal
//! (nhat_xi.nhat_eta = -c): summing their squares would be wrong by the cross term, the
//! error largest at a panel corner.  But MHD::ConToPrim calls
//! Coordinates::GnomonicEquiangleRaiseVelMHD immediately afterwards on the cubed sphere,
//! and that REBUILDS bcc0 in the orthonormal frame {rhat, e_xi, (e_eta - c e_xi)/s} --
//! precisely so that every consumer which squares and sums bcc0 is correct with no
//! further change.  So bcc0 needs no metric here, and this is bitwise the magnetic
//! energy that GnomonicEquiangleRaiseVelMHD (or, off the cubed sphere, ConsToPrim)
//! subtracted, which is what makes an e extracted with it consistent with w0(IEN).
//!
//! WHICH FIELD.  bcc0 is the cell-centred form of the CURRENT b0, exactly: the last
//! ConToPrim built it from those faces and nothing writes b0 again until MHD::CT.  Any
//! operator that runs inside a stage BEFORE CT -- the source terms, the two-stream
//! apply, the implicit radial conduction -- therefore has in bcc0 the field that
//! belongs to the b0 it can see, and paying for a re-average from the faces would
//! return the same numbers.  An operator placed AFTER CT would instead have to rebuild
//! from the post-CT faces with the ConsToPrim weights; no caller does that today.

KOKKOS_INLINE_FUNCTION
Real MagEnergyCC(const DvceArray5D<Real> &bcc, const int m, const int k, const int j,
                 const int i) {
  return 0.5*(SQR(bcc(m,IBX,k,j,i)) + SQR(bcc(m,IBY,k,j,i)) + SQR(bcc(m,IBZ,k,j,i)));
}

//----------------------------------------------------------------------------------------
//! \fn Real EintFromCons
//! \brief the MHD form: e = u0(IEN) - KE - emag - rho*Phi.  `emag` is MagEnergyCC of
//! the same cell (0 in hydro, where the hydro overload below is the one to call).  The
//! subtractions are kept separate and in this order because that is the order
//! GnomonicEquiangleRaiseVelMHD and ConsToPrim use, so the result agrees with theirs in
//! the last bits rather than merely to round-off.

KOKKOS_INLINE_FUNCTION
Real EintFromCons(const DvceArray5D<Real> &u0, const int m, const int k, const int j,
                  const int i, const Real cosc, const bool cs, const bool etg,
                  const Real phi, const Real emag) {
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
  Real ei = u0(m,IEN,k,j,i) - ekin - emag;
  if (etg) ei -= d*phi;
  return ei;
}

//----------------------------------------------------------------------------------------
//! \fn Real EintFromCons
//! \brief the HYDRO form, unchanged: e = u0(IEN) - KE - rho*Phi.  Kept as its own
//! overload rather than defaulting emag to zero so that the hydro arithmetic really is
//! the arithmetic that stood here, with no subtraction of a zero to flip a signed zero.

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
