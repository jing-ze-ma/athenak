#ifndef EOS_GENERAL_C2P_HYD_HPP_
#define EOS_GENERAL_C2P_HYD_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file general_c2p_hyd.hpp
//! \brief Inline conserved->primitive conversion for non-relativistic hydrodynamics with
//! a general EOS.
//!
//! The conserved inversion itself is identical to the ideal-gas case and stays purely
//! algebraic: e = E - E_kin involves no EOS at all. What differs is
//!   (1) the floors, which are specified in p and T but applied to e, and therefore need
//!       the (density-dependent) EOS inverses rather than a constant 1/(gamma-1); and
//!   (2) the derived quantities p and Gamma_1, which are evaluated here, ONCE per cell,
//!       and stored so that the Riemann solvers never have to call the EOS.
//!
//! This routine performs EXACTLY ONE temperature solve per cell in the common case, which
//! for a general EOS is the only expensive operation in the whole time step. Everything
//! else here -- p, Gamma_1, and both floor tests -- is cheap once T is known. Keep it
//! that way: an EOS accessor called without a temperature argument re-solves the root
//! find.

//----------------------------------------------------------------------------------------
//! \fn void SingleC2P_GeneralHyd()
//! \brief Converts a single state of conserved variables into primitive variables for
//! non-relativistic hydrodynamics with a general EOS, and evaluates the derived
//! thermodynamic quantities that are subsequently reconstructed to interfaces.
//! Conserved = (d,M1,M2,M3,E), Primitive = (d,vx,vy,vz,e)
//! where E=total energy density and e=internal energy density.
//! `tguess` warm starts the temperature solve (the cached T in this cell from the
//! previous stage); the solved temperature is returned in `temp`.

KOKKOS_INLINE_FUNCTION
void SingleC2P_GeneralHyd(HydCons1D &u, const EOS_Data &eos, HydPrim1D &w,
                          const Real tguess, Real &temp, Real &pgas, Real &g1,
                          bool &dfloor_used, bool &efloor_used, bool &tfloor_used) {
  // apply density floor, without changing momentum or energy
  if (u.d < eos.dfloor) {
    u.d = eos.dfloor;
    dfloor_used = true;
  }
  w.d = u.d;

  // compute velocities
  Real di = 1.0/u.d;
  w.vx = di*u.mx;
  w.vy = di*u.my;
  w.vz = di*u.mz;

  // set internal energy, apply floor, correct total energy (if needed)
  Real e_k = 0.5*di*(SQR(u.mx) + SQR(u.my) + SQR(u.mz));
  w.e = (u.e - e_k);

  // Solve for the temperature. For a general EOS this is the one expensive EOS call in
  // the time step, so it happens here, once, and everything below reuses the result. A
  // badly under-resolved cell can leave e non-positive, where the root find has nothing
  // to converge to; that case skips the solve and goes straight to the pressure floor.
  bool e_positive = (w.e > 0.0);
  temp = e_positive ? eos.Temperature(w.d, w.e, tguess) : -1.0;

  // Apply the pressure floor. The test is on p rather than on e so that the inverse
  // e(d,pfloor) -- density dependent, and a root find for a general EOS -- is evaluated
  // only in the rare cells where the floor actually trips, instead of in every cell as a
  // constant would be. See EOS_Data::BelowPressureFloor().
  if (!e_positive || eos.BelowPressureFloor(w.d, w.e, temp)) {
    // the three-argument form hands back the temperature the inversion solved for,
    // so the floored cell does not pay for a second root find
    w.e = eos.EnergyFromPressure(w.d, eos.pfloor, temp);
    u.e = w.e + e_k;
    efloor_used = true;
  }

  // Apply the temperature floor. Free to test now that T is known, and e(d,tfloor) is a
  // direct evaluation rather than an inversion, so this costs no root find either.
  if (temp < eos.tfloor) {
    w.e = eos.EnergyFromTemperature(w.d, eos.tfloor);
    u.e = w.e + e_k;
    temp = eos.tfloor;
    tfloor_used = true;
  }

  // Apply the entropy floor. Only meaningful while the general EOS evaluates a gamma law,
  // where it is the ideal-gas floor verbatim; under a tabulated EOS this is a no-op and a
  // run that sets sfloor is refused at startup. See EOS_Data::ApplyEntropyFloor().
  if (eos.ApplyEntropyFloor(w.d, di, w.e)) {
    temp = eos.Temperature(w.d, w.e, temp);
    efloor_used = true;
  }

  // Evaluate the derived thermodynamic quantities that the Riemann solvers will need.
  // Both are cheap: the temperature they depend on has already been solved for above.
  pgas = eos.Pressure(w.d, w.e, temp);
  g1 = eos.Gamma1(w.d, w.e, temp);
  return;
}

#endif // EOS_GENERAL_C2P_HYD_HPP_
