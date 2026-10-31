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
//! \fn void GeneralMHDFloorsPost()
//! \brief share3-1031: THE energy floors of the general-EOS MHD inversion on an internal
//! energy e = E - KE - ME with final kinetic and magnetic energies: the one temperature
//! solve, the pressure floor, the temperature floor, the entropy floor (`sflr`) and the
//! p, Gamma_1 re-evaluation.  Shared by SingleC2P_GeneralMHD (own = !defer_cons_floors,
//! sflr = true) and, on the cubed sphere, GnomonicRaiseVelMHDFloors (own = true, sflr =
//! false: the entropy floor that path has never applied).  `own`: the total `etot` is
//! rebuilt as e + e_k + e_m at each floor.  `stale`: e moved (p, Gamma_1 re-evaluated).
//! The arithmetic is the former body of SingleC2P_GeneralMHD, in its order.

KOKKOS_INLINE_FUNCTION
void GeneralMHDFloorsPost(const EOS_Data &eos, const Real d, const Real di,
                          const Real tguess, const bool own, const bool sflr,
                          const Real e_k, const Real e_m, Real &eint, Real &etot,
                          Real &temp, Real &pgas, Real &g1, bool &efloor_used,
                          bool &tfloor_used, bool &tclamp_used, bool &stale) {
  bool e_positive = (eint > 0.0);
  temp = -1.0;

  stale = !e_positive;
  if (e_positive) {
    eos.TemperaturePressureGamma1(d, eint, tguess, temp, pgas, g1, tclamp_used);
  }

  if (!e_positive || pgas < eos.pfloor) {
    eint = eos.EnergyFromPressure(d, eos.pfloor, temp);
    if (own) etot = eint + e_k + e_m;
    efloor_used = true;
    stale = true;
  }

  if (temp < eos.tfloor) {
    eint = eos.EnergyFromTemperature(d, eos.tfloor);
    if (own) etot = eint + e_k + e_m;
    temp = eos.tfloor;
    tfloor_used = true;
    stale = true;
  }

  if (sflr && eos.ApplyEntropyFloor(d, di, eint)) {
    temp = eos.Temperature(d, eint, temp);
    efloor_used = true;
    stale = true;
  }

  if (stale) {
    eos.PressureAndGamma1(d, eint, temp, pgas, g1);
  }
}

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
                          bool &dfloor_used, bool &efloor_used, bool &tfloor_used,
                          bool &vceil_used, bool &vceil_test, bool &tclamp_used) {
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
  // THE VELOCITY CEILING: identical to the block in ideal_c2p_mhd.hpp; see eos.hpp.
  // <mhd>/vceil_thermalise leaves the TOTAL energy alone, so that w.e below picks the
  // clipped kinetic energy up as internal energy; the magnetic energy is not involved
  // either way.  The converted energy is NOT accumulated here: this inversion has no
  // energy accumulator (see EOS_Data::vceil_thermalise).
  if (eos.vceil > 0.0) {
    const Real vsq = SQR(w.vx) + SQR(w.vy) + SQR(w.vz);
    if (vsq > SQR(eos.vceil)) {
      vceil_test = true;
      if (!eos.defer_cons_floors) {
        const Real fs = eos.vceil/sqrt(vsq);
        u.mx *= fs; u.my *= fs; u.mz *= fs;
        w.vx *= fs; w.vy *= fs; w.vz *= fs;
        if (!eos.vceil_thermalise) u.e -= (1.0 - fs*fs)*e_k;
        e_k *= fs*fs;
        vceil_used = true;
      }
    }
  }
  w.e = (u.e - e_k - e_m);

  // the floors on the final kinetic and magnetic energies (shared with the cubed
  // sphere's GnomonicRaiseVelMHDFloors)
  bool stale;
  GeneralMHDFloorsPost(eos, w.d, di, tguess, !eos.defer_cons_floors, true, e_k, e_m, w.e,
                       u.e, temp, pgas, g1, efloor_used, tfloor_used, tclamp_used, stale);
  return;
}

#endif // EOS_GENERAL_C2P_MHD_HPP_
