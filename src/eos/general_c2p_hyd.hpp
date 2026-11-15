#ifndef EOS_GENERAL_C2P_HYD_HPP_
#define EOS_GENERAL_C2P_HYD_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file general_c2p_hyd.hpp
//! \brief Inline conserved->primitive conversion for non-relativistic hydrodynamics with
//! a general EOS.  The inversion e = E - E_kin is algebraic as for an ideal gas; the
//! floors need the EOS inverses, and p and Gamma_1 are evaluated here once per cell for
//! the Riemann solvers.  The common case costs exactly one temperature solve per cell.

//----------------------------------------------------------------------------------------
//! \fn void SingleC2P_GeneralHydLegacy()
//! \brief the inversion without the floor switches (EOS_Data::floors_legacy), a separate
//! function so its kernel keeps the original arithmetic bitwise.  Changes meant for the
//! default path must be made in both functions.

KOKKOS_INLINE_FUNCTION
void SingleC2P_GeneralHydLegacy(HydCons1D &u, const EOS_Data &eos, HydPrim1D &w,
                          const Real tguess, Real &temp, Real &pgas, Real &g1,
                          bool &dfloor_used, bool &efloor_used, bool &tfloor_used,
                          bool &tclamp_used) {
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

  // solve for the temperature once; a non-positive e skips the solve and goes to the
  // pressure floor
  bool e_positive = (w.e > 0.0);
  temp = -1.0;

  // pressure floor, tested on p so that the inverse e(d, pfloor) is evaluated only where
  // the floor trips (EOS_Data::BelowPressureFloor).  One table evaluation gives the floor
  // test and the p, Gamma_1 below; stale marks paths where a floor moved the state
  bool stale = !e_positive;
  if (e_positive) {
    // the inversion and the evaluation share their logarithms
    eos.TemperaturePressureGamma1(w.d, w.e, tguess, temp, pgas, g1, tclamp_used);
  }

  if (!e_positive || pgas < eos.pfloor) {
    // also returns the temperature, so the floored cell needs no second root find
    w.e = eos.EnergyFromPressure(w.d, eos.pfloor, temp);
    if (!eos.defer_cons_floors) u.e = w.e + e_k;
    efloor_used = true;
    stale = true;
  }

  // temperature floor (e(d, tfloor) is a direct evaluation)
  if (temp < eos.tfloor) {
    w.e = eos.EnergyFromTemperature(w.d, eos.tfloor);
    if (!eos.defer_cons_floors) u.e = w.e + e_k;
    temp = eos.tfloor;
    tfloor_used = true;
    stale = true;
  }

  // entropy floor, gamma law only (EOS_Data::ApplyEntropyFloor)
  if (eos.ApplyEntropyFloor(w.d, di, w.e)) {
    temp = eos.Temperature(w.d, w.e, temp);
    efloor_used = true;
    stale = true;
  }

  // derived quantities for the Riemann solvers, at the solved temperature
  if (stale) {
    eos.PressureAndGamma1(w.d, w.e, temp, pgas, g1);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void GeneralHydFloorsPost()
//! \brief energy floors of the general-EOS hydro inversion on e = E - KE with a final
//! kinetic energy: the temperature solve, the pressure floor (paid from the kinetic
//! energy under efloor_from_ekin), the temperature floor, the entropy floor (sflr), the
//! efloor_as_tfloor rebuild and the p, Gamma_1 re-evaluation.  Used by
//! SingleC2P_GeneralHyd (own = !defer_cons_floors, sflr = true) and on the cubed sphere
//! by GnomonicRaiseVelFloors (own = true, sflr = false).  With own the floors account the
//! energy they create in de, rebuild etot = e + ekin and may pay from the kinetic energy
//! (fvk is the momentum factor, mom_scaled true).  floor_fired: a floor (not the clamp
//! repair) moved e; stale: e moved and p, Gamma_1 were re-evaluated.

KOKKOS_INLINE_FUNCTION
void GeneralHydFloorsPost(const EOS_Data &eos, const Real d, const Real di,
                          const Real tguess, const bool own, const bool sflr,
                          Real &eint, Real &ekin, Real &etot, Real &fvk,
                          Real &temp, Real &pgas, Real &g1,
                          bool &efloor_used, bool &tfloor_used, bool &tclamp_used,
                          bool &tset_used, bool &mom_scaled, bool &floor_fired,
                          bool &stale, Real &de) {
  // solve for the temperature once; a non-positive e skips the solve and goes to the
  // pressure floor
  bool e_positive = (eint > 0.0);
  temp = -1.0;

  // pressure floor, tested on p so that the inverse e(d, pfloor) is evaluated only where
  // the floor trips (EOS_Data::BelowPressureFloor).  One table evaluation gives the floor
  // test and the p, Gamma_1 below; stale marks paths where a floor moved the state
  stale = !e_positive;
  if (e_positive) {
    // the inversion and the evaluation share their logarithms
    eos.TemperaturePressureGamma1(d, eint, tguess, temp, pgas, g1, tclamp_used);
  }

  if (!e_positive || pgas < eos.pfloor) {
    // also returns the temperature, so the floored cell needs no second root find
    const Real efl = eos.EnergyFromPressure(d, eos.pfloor, temp);
    // efloor_from_ekin: hold E fixed and rescale the momentum, so energy is created only
    // when E itself is below the floor.  Not with own = false (cubed sphere ConsToPrim),
    // where KE is not the metric kinetic energy; GnomonicRaiseVelFloors applies it there.
    if (eos.efloor_from_ekin && own && ekin > 0.0 && (etot - efl) < ekin) {
      Real ek_new = etot - efl;
      if (!(ek_new > 0.0)) ek_new = 0.0;
      fvk = sqrt(ek_new/ekin);
      de += (efl + ek_new) - etot;
      ekin = ek_new;
      mom_scaled = true;
    } else if (own) {
      de += efl - eint;
    }
    eint = efl;
    if (own) etot = eint + ekin;
    efloor_used = true;
    stale = true;
  }

  // temperature floor (e(d, tfloor) is a direct evaluation)
  if (temp < eos.tfloor) {
    const Real etf = eos.EnergyFromTemperature(d, eos.tfloor);
    if (own) de += etf - eint;
    eint = etf;
    if (own) etot = eint + ekin;
    temp = eos.tfloor;
    tfloor_used = true;
    stale = true;
  }

  // entropy floor, gamma law only (EOS_Data::ApplyEntropyFloor)
  if (sflr && eos.ApplyEntropyFloor(d, di, eint)) {
    temp = eos.Temperature(d, eint, temp);
    efloor_used = true;
    stale = true;
  }
  // whether a floor fired (FOFC flags cells on it; the table clamp is not a floor)
  floor_fired = stale;

  // efloor_as_tfloor: after a floor or a table clamp, set e = e(rho, T) at the
  // temperature settled on, so re-inverting the cell returns this T (raises e at the low
  // table edge, lowers it at the high edge)
  if (eos.efloor_as_tfloor && (stale || tclamp_used)) {
    const Real tset = (temp > eos.tfloor) ? temp : eos.tfloor;
    const Real ets = eos.EnergyFromTemperature(d, tset);
    if (own) de += ets - eint;
    eint = ets;
    if (own) etot = eint + ekin;
    temp = tset;
    tset_used = true;
    stale = true;
  }

  // derived quantities for the Riemann solvers, at the solved temperature
  if (stale) {
    eos.PressureAndGamma1(d, eint, temp, pgas, g1);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void SingleC2P_GeneralHyd()
//! \brief Converts a single state of conserved variables into primitive variables for
//! non-relativistic hydrodynamics with a general EOS, and evaluates the derived
//! quantities. Conserved = (d,M1,M2,M3,E), Primitive = (d,vx,vy,vz,e).  tguess
//! warm-starts the temperature solve and temp returns it.  efloor_de accumulates the
//! internal energy created by the energy floor and mom_scaled flags a momentum rescale
//! (efloor_from_ekin; the caller writes the momentum back); vceil_de accumulates the
//! kinetic energy removed by the velocity ceiling.

//! `efloor_de` accumulates the internal energy density this call CREATED at the energy
//! floor, and `mom_scaled` says whether the momentum was rescaled to pay for it
//! (<block>/efloor_from_ekin), in which case the caller must write u.mx/my/mz back.

//! `vceil_de` accumulates the kinetic energy density the VELOCITY CEILING removed from
//! the conserved total.

KOKKOS_INLINE_FUNCTION
void SingleC2P_GeneralHyd(HydCons1D &u, const EOS_Data &eos, HydPrim1D &w,
                          const Real tguess, Real &temp, Real &pgas, Real &g1,
                          bool &dfloor_used, bool &efloor_used, bool &tfloor_used,
                          Real &efloor_de, bool &mom_scaled, Real &dfloor_fv,
                          bool &vceil_used, bool &vceil_test, bool &tclamp_used,
                          bool &tset_used, Real &vceil_de) {
  // density floor: by default raise u.d only; dfloor_keep_velocity scales the momentum by
  // fv = d_old/dfloor and removes the lost kinetic energy (on the cubed sphere the energy
  // correction is deferred to GnomonicEquiangleRaiseVel through dfloor_fv)
  dfloor_fv = 1.0;
  if (u.d < eos.dfloor) {
    if (eos.dfloor_keep_velocity) {
      const Real fv = (u.d > 0.0) ? u.d/eos.dfloor : 0.0;
      if (!eos.defer_cons_floors) {
        const Real ke_old = (u.d > 0.0)
            ? 0.5*(SQR(u.mx) + SQR(u.my) + SQR(u.mz))/u.d : 0.0;
        if (eos.dfloor_keep_velocity) u.e -= (1.0 - fv*fv*fv)*ke_old;
      }
      if (eos.dfloor_keep_velocity) {
        u.mx *= fv;
        u.my *= fv;
        u.mz *= fv;
        mom_scaled = true;
      }
      dfloor_fv = fv;
    }
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
  // velocity ceiling (EOS_Data::vceil): scale the momentum by vceil/|v| and remove
  // (1 - fs^2) KE from the total; deferred to GnomonicEquiangleRaiseVel on the cubed
  // sphere. vceil_test records the decision alone, for the FOFC floor test
  if (eos.vceil > 0.0) {
    const Real vsq = SQR(w.vx) + SQR(w.vy) + SQR(w.vz);
    if (vsq > SQR(eos.vceil)) {
      vceil_test = true;
      if (!eos.defer_cons_floors) {
        const Real fs = eos.vceil/sqrt(vsq);
        u.mx *= fs; u.my *= fs; u.mz *= fs;
        w.vx *= fs; w.vy *= fs; w.vz *= fs;
        u.e -= (1.0 - fs*fs)*e_k;
        vceil_de += (1.0 - fs*fs)*e_k;
        e_k *= fs*fs;
        vceil_used = true;
        mom_scaled = true;
      }
    }
  }
  w.e = (u.e - e_k);

  // the floors on the final kinetic energy (shared with the cubed sphere's
  // GnomonicRaiseVelFloors)
  Real fvk = 1.0;
  bool msc = false, ffired, stale;
  GeneralHydFloorsPost(eos, w.d, di, tguess, !eos.defer_cons_floors, true, w.e, e_k, u.e,
                       fvk, temp, pgas, g1, efloor_used, tfloor_used, tclamp_used,
                       tset_used, msc, ffired, stale, efloor_de);
  if (msc) {
    u.mx *= fvk; u.my *= fvk; u.mz *= fvk;
    w.vx *= fvk; w.vy *= fvk; w.vz *= fvk;
    mom_scaled = true;
  }
  return;
}

#endif // EOS_GENERAL_C2P_HYD_HPP_
