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
//! \fn void SingleC2P_GeneralHydLegacy()
//! \brief THE PRE-FLOOR-SWITCH inversion, kept verbatim.  See the note on
//! EOS_Data::floors_legacy: with none of the floor switches set SingleC2P_GeneralHyd()
//! below is algebraically this routine, but its extra tests and temporaries change how
//! the device compiler contracts the arithmetic, so a run built with the switches
//! available but unused would not reproduce one built without them.  It is a separate
//! function rather than a branch inside the other one precisely so that the kernel that
//! calls it inlines this code and nothing else.  ConsToPrim() selects between the two on
//! eos.floors_legacy; any change meant to apply by DEFAULT has to be made in both.

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

  // Solve for the temperature. For a general EOS this is the one expensive EOS call in
  // the time step, so it happens here, once, and everything below reuses the result. A
  // badly under-resolved cell can leave e non-positive, where the root find has nothing
  // to converge to; that case skips the solve and goes straight to the pressure floor.
  bool e_positive = (w.e > 0.0);
  temp = -1.0;

  // Apply the pressure floor. The test is on p rather than on e so that the inverse
  // e(d,pfloor) -- density dependent, and a root find for a general EOS -- is evaluated
  // only in the rare cells where the floor actually trips, instead of in every cell as a
  // constant would be. See EOS_Data::BelowPressureFloor().
  // ONE table evaluation now serves three consumers -- the pressure-floor test below, and
  // the p and Gamma_1 the Riemann solvers need at the bottom -- because all three are at
  // the same (d,T). `stale` marks the rare paths where a floor moves the state afterwards
  // and they have to be redone.
  bool stale = !e_positive;
  if (e_positive) {
    // fused: the inversion and the evaluation that follows it share their logarithms
    eos.TemperaturePressureGamma1(w.d, w.e, tguess, temp, pgas, g1, tclamp_used);
  }

  if (!e_positive || pgas < eos.pfloor) {
    // the three-argument form hands back the temperature the inversion solved for,
    // so the floored cell does not pay for a second root find
    w.e = eos.EnergyFromPressure(w.d, eos.pfloor, temp);
    if (!eos.defer_cons_floors) u.e = w.e + e_k;
    efloor_used = true;
    stale = true;
  }

  // Apply the temperature floor. Free to test now that T is known, and e(d,tfloor) is a
  // direct evaluation rather than an inversion, so this costs no root find either.
  if (temp < eos.tfloor) {
    w.e = eos.EnergyFromTemperature(w.d, eos.tfloor);
    if (!eos.defer_cons_floors) u.e = w.e + e_k;
    temp = eos.tfloor;
    tfloor_used = true;
    stale = true;
  }

  // Apply the entropy floor. Only meaningful while the general EOS evaluates a gamma law,
  // where it is the ideal-gas floor verbatim; under a tabulated EOS this is a no-op and a
  // run that sets sfloor is refused at startup. See EOS_Data::ApplyEntropyFloor().
  if (eos.ApplyEntropyFloor(w.d, di, w.e)) {
    temp = eos.Temperature(w.d, w.e, temp);
    efloor_used = true;
    stale = true;
  }

  // Evaluate the derived thermodynamic quantities that the Riemann solvers will need.
  // Both are cheap: the temperature they depend on has already been solved for above.
  if (stale) {
    eos.PressureAndGamma1(w.d, w.e, temp, pgas, g1);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void SingleC2P_GeneralHyd()
//! \brief Converts a single state of conserved variables into primitive variables for
//! non-relativistic hydrodynamics with a general EOS, and evaluates the derived
//! thermodynamic quantities that are subsequently reconstructed to interfaces.
//! Conserved = (d,M1,M2,M3,E), Primitive = (d,vx,vy,vz,e)
//! where E=total energy density and e=internal energy density.
//! `tguess` warm starts the temperature solve (the cached T in this cell from the
//! previous stage); the solved temperature is returned in `temp`.

//! `efloor_de` accumulates the internal energy density this call CREATED at the energy
//! floor, and `mom_scaled` says whether the momentum was rescaled to pay for it
//! (<block>/efloor_from_ekin), in which case the caller must write u.mx/my/mz back.

KOKKOS_INLINE_FUNCTION
void SingleC2P_GeneralHyd(HydCons1D &u, const EOS_Data &eos, HydPrim1D &w,
                          const Real tguess, Real &temp, Real &pgas, Real &g1,
                          bool &dfloor_used, bool &efloor_used, bool &tfloor_used,
                          Real &efloor_de, bool &mom_scaled, Real &dfloor_fv,
                          bool &vceil_used, bool &vceil_test, bool &tclamp_used) {
  // THE DENSITY FLOOR.  Default: raise u.d and leave m and E alone -- which changes the
  // velocity and creates internal energy, because the kinetic share of the fixed total
  // drops when rho goes up.  <block>/dfloor_keep_velocity instead scales the momentum by
  // fv = d_old/dfloor so KE -> fv^3 KE, and removes exactly that from the total energy,
  // leaving the internal energy untouched.  On the cubed sphere (defer_cons_floors) the
  // e_k below is the ORTHONORMAL kinetic energy and is NOT the real one, so the energy
  // correction is deferred to GnomonicEquiangleRaiseVel, which has the metric; it is
  // handed fv through dfloor_fv and applies (1/fv^3 - 1) times the metric KE there.
  dfloor_fv = 1.0;
  if (u.d < eos.dfloor) {
    if (eos.dfloor_keep_velocity || eos.dfloor_keep_temperature) {
      const Real fv = (u.d > 0.0) ? u.d/eos.dfloor : 0.0;
      if (!eos.defer_cons_floors) {
        const Real ke_old = (u.d > 0.0)
            ? 0.5*(SQR(u.mx) + SQR(u.my) + SQR(u.mz))/u.d : 0.0;
        const Real eint_old = u.e - ke_old;
        if (eos.dfloor_keep_velocity) u.e -= (1.0 - fv*fv*fv)*ke_old;
        // <block>/dfloor_keep_temperature: e -> fv e, which holds e/rho at the value the
        // plain floor would have inflated.  EXACT only for an ideal gas; the tabulated
        // EOS is not called here (see the note in eos.hpp).
        if (eos.dfloor_keep_temperature && fv > 0.0) u.e -= (1.0 - fv)*eint_old;
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
  // THE VELOCITY CEILING: see the note in eos.hpp.  Scale the momentum by
  // fs = vceil/|v| and take (1 - fs^2) KE out of the total, which leaves the internal
  // energy exactly unchanged.  On the cubed sphere (defer_cons_floors) |v| and e_k here
  // are the ORTHONORMAL ones and are not the real ones, so the ceiling is applied by
  // GnomonicEquiangleRaiseVel instead, which owns the metric.
  // `vceil_test` is the DECISION alone, recorded even where the ceiling itself is
  // deferred, so that the FOFC floor-test pass can flag the cell.  It performs no write.
  if (eos.vceil > 0.0) {
    const Real vsq = SQR(w.vx) + SQR(w.vy) + SQR(w.vz);
    if (vsq > SQR(eos.vceil)) {
      vceil_test = true;
      if (!eos.defer_cons_floors) {
        const Real fs = eos.vceil/sqrt(vsq);
        u.mx *= fs; u.my *= fs; u.mz *= fs;
        w.vx *= fs; w.vy *= fs; w.vz *= fs;
        u.e -= (1.0 - fs*fs)*e_k;
        e_k *= fs*fs;
        vceil_used = true;
        mom_scaled = true;
      }
    }
  }
  w.e = (u.e - e_k);

  // Solve for the temperature. For a general EOS this is the one expensive EOS call in
  // the time step, so it happens here, once, and everything below reuses the result. A
  // badly under-resolved cell can leave e non-positive, where the root find has nothing
  // to converge to; that case skips the solve and goes straight to the pressure floor.
  bool e_positive = (w.e > 0.0);
  temp = -1.0;

  // Apply the pressure floor. The test is on p rather than on e so that the inverse
  // e(d,pfloor) -- density dependent, and a root find for a general EOS -- is evaluated
  // only in the rare cells where the floor actually trips, instead of in every cell as a
  // constant would be. See EOS_Data::BelowPressureFloor().
  // ONE table evaluation now serves three consumers -- the pressure-floor test below, and
  // the p and Gamma_1 the Riemann solvers need at the bottom -- because all three are at
  // the same (d,T). `stale` marks the rare paths where a floor moves the state afterwards
  // and they have to be redone.
  bool stale = !e_positive;
  if (e_positive) {
    // fused: the inversion and the evaluation that follows it share their logarithms
    eos.TemperaturePressureGamma1(w.d, w.e, tguess, temp, pgas, g1, tclamp_used);
  }

  if (!e_positive || pgas < eos.pfloor) {
    // the three-argument form hands back the temperature the inversion solved for,
    // so the floored cell does not pay for a second root find
    const Real efl = eos.EnergyFromPressure(w.d, eos.pfloor, temp);
    // PAY FROM THE KINETIC ENERGY FIRST.  Rebuilding u.e as efl + e_k with e_k untouched
    // donates efl - w.e, which is the whole kinetic energy whenever the update left
    // w.e = u.e - e_k negative.  Holding u.e fixed and rescaling the momentum instead
    // creates nothing unless the TOTAL energy is itself below the floor.
    //
    // NOT under defer_cons_floors (the cubed sphere): there e_k above is the ORTHONORMAL
    // kinetic energy, which is not the kinetic energy at all on a non-orthogonal tangent
    // basis, and this floor is provisional -- Coordinates::GnomonicEquiangleRaiseVel
    // re-applies it to the metric-correct energy and does the payment and the accounting
    // there.  Touching the conserved momentum here would freeze the wrong KE into it.
    if (eos.efloor_from_ekin && !eos.defer_cons_floors &&
        e_k > 0.0 && (u.e - efl) < e_k) {
      Real ek_new = u.e - efl;
      if (!(ek_new > 0.0)) ek_new = 0.0;
      const Real fv = sqrt(ek_new/e_k);
      u.mx *= fv; u.my *= fv; u.mz *= fv;
      w.vx *= fv; w.vy *= fv; w.vz *= fv;
      efloor_de += (efl + ek_new) - u.e;
      e_k = ek_new;
      mom_scaled = true;
    } else if (!eos.defer_cons_floors) {
      efloor_de += efl - w.e;
    }
    w.e = efl;
    if (!eos.defer_cons_floors) u.e = w.e + e_k;
    efloor_used = true;
    stale = true;
  }

  // Apply the temperature floor. Free to test now that T is known, and e(d,tfloor) is a
  // direct evaluation rather than an inversion, so this costs no root find either.
  if (temp < eos.tfloor) {
    const Real etf = eos.EnergyFromTemperature(w.d, eos.tfloor);
    if (!eos.defer_cons_floors) efloor_de += etf - w.e;
    w.e = etf;
    if (!eos.defer_cons_floors) u.e = w.e + e_k;
    temp = eos.tfloor;
    tfloor_used = true;
    stale = true;
  }

  // Apply the entropy floor. Only meaningful while the general EOS evaluates a gamma law,
  // where it is the ideal-gas floor verbatim; under a tabulated EOS this is a no-op and a
  // run that sets sfloor is refused at startup. See EOS_Data::ApplyEntropyFloor().
  if (eos.ApplyEntropyFloor(w.d, di, w.e)) {
    temp = eos.Temperature(w.d, w.e, temp);
    efloor_used = true;
    stale = true;
  }

  // Evaluate the derived thermodynamic quantities that the Riemann solvers will need.
  // Both are cheap: the temperature they depend on has already been solved for above.
  if (stale) {
    eos.PressureAndGamma1(w.d, w.e, temp, pgas, g1);
  }
  return;
}

#endif // EOS_GENERAL_C2P_HYD_HPP_
