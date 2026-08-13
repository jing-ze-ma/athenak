#ifndef EOS_GENERAL_C2P_MHD_HPP_
#define EOS_GENERAL_C2P_MHD_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file general_c2p_mhd.hpp
//! \brief Inline conserved->primitive conversion for non-relativistic MHD with a general
//! EOS. See general_c2p_hyd.hpp for the rationale; the only difference here is that the
//! magnetic energy is also subtracted off before the internal energy is obtained.

//----------------------------------------------------------------------------------------
//! \fn void SingleC2P_GeneralMHD()
//! \brief Converts a single state of conserved variables into primitive variables for
//! non-relativistic MHD with a general EOS, and evaluates the derived thermodynamic
//! quantities that are subsequently reconstructed to interfaces. Note the input CONSERVED
//! state contains cell-centered magnetic fields, but the PRIMITIVE state returned through
//! the arguments does not. As in the hydro version, `tguess` warm starts the temperature
//! solve and the solved temperature is returned in `temp`; there is exactly one such
//! solve per cell in the common case.

KOKKOS_INLINE_FUNCTION
void SingleC2P_GeneralMHD(MHDCons1D &u, const EOS_Data &eos, HydPrim1D &w,
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

  // set internal energy, apply floor, correcting total energy
  Real e_k = 0.5*di*(SQR(u.mx) + SQR(u.my) + SQR(u.mz));
  Real e_m = 0.5*(SQR(u.bx) + SQR(u.by) + SQR(u.bz));
  w.e = (u.e - e_k - e_m);

  // Solve for the temperature, once per cell; see general_c2p_hyd.hpp for why this is the
  // only expensive EOS call here and why everything below reuses its result.
  bool e_positive = (w.e > 0.0);
  temp = e_positive ? eos.Temperature(w.d, w.e, tguess) : -1.0;

  // Apply the pressure floor, testing on p so that the inversion e(d,pfloor) is performed
  // only in the rare cells where the floor trips.
  if (!e_positive || eos.BelowPressureFloor(w.d, w.e, temp)) {
    // three-argument form: reuses the temperature the inversion solved for (see
    // general_c2p_hyd.hpp) instead of re-solving the root find
    w.e = eos.EnergyFromPressure(w.d, eos.pfloor, temp);
    u.e = w.e + e_k + e_m;
    efloor_used = true;
  }

  // Apply the temperature floor, now free to test and free to apply.
  if (temp < eos.tfloor) {
    w.e = eos.EnergyFromTemperature(w.d, eos.tfloor);
    u.e = w.e + e_k + e_m;
    temp = eos.tfloor;
    tfloor_used = true;
  }

  // Apply the entropy floor; a no-op under a tabulated EOS, which refuses sfloor at
  // startup instead. See general_c2p_hyd.hpp and EOS_Data::ApplyEntropyFloor().
  if (eos.ApplyEntropyFloor(w.d, di, w.e)) {
    temp = eos.Temperature(w.d, w.e, temp);
    efloor_used = true;
  }

  // Derived thermodynamic quantities, cheap now that T is known (see general_c2p_hyd.hpp)
  pgas = eos.Pressure(w.d, w.e, temp);
  g1 = eos.Gamma1(w.d, w.e, temp);
  return;
}

#endif // EOS_GENERAL_C2P_MHD_HPP_
