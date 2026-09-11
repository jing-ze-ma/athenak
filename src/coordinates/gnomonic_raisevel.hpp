#ifndef COORDINATES_GNOMONIC_RAISEVEL_HPP_
#define COORDINATES_GNOMONIC_RAISEVEL_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file gnomonic_raisevel.hpp
//! \brief the per-cell body of Coordinates::GnomonicEquiangleRaiseVel, as a free
//! device function.
//!
//! WHY IT LIVES HERE.  On the cubed sphere the conserved-to-primitive inversion is split
//! in two: EquationOfState::ConsToPrim does the part that needs no metric, and
//! GnomonicEquiangleRaiseVel does the part that does -- raising the covariant momentum,
//! forming the metric kinetic energy, applying the velocity CEILING and RE-APPLYING the
//! pressure/temperature/energy floors against the corrected internal energy
//! (EOS_Data::defer_cons_floors).  A cell that needs one of those deferred floors is
//! therefore invisible to ConsToPrim, and the FOFC floor-TEST pass, which is nothing but
//! a ConsToPrim over a trial state, never flagged it.  Factoring the body out lets the
//! test pass run exactly the same arithmetic on the trial state (Hydro::FOFC), with
//! `testonly` suppressing nothing: the caller simply discards the state instead of
//! writing it back.
//!
//! EVERYTHING IS LOCAL.  The caller loads d,m1,m2,m3,etot out of its conserved array and
//! writes the (possibly modified) m1,m2,m3,etot back.  Writing them back unconditionally
//! is bit-for-bit what the old kernel did, which wrote them only when something fired:
//! when nothing fires the locals still hold the bits that were loaded.

#include "athena.hpp"
#include "eos/eos.hpp"

//----------------------------------------------------------------------------------------
//! \fn void GnomonicRaiseVelFloors
//! \brief metric velocity raise + deferred ceiling/floors for one cubed-sphere cell.
//!
//! \param[in]     c        cos_cell, the angle between the two tangent basis vectors
//! \param[in]     fv       dfl_fv for this cell (the density floor's d_old/dfloor)
//! \param[in]     twarm    the temperature ConsToPrim found, used to warm-start the
//!                         general-EOS root find
//! \param[in]     act      true in ACTIVE cells only; gates the velocity ceiling
//! \param[in,out] m1,m2,m3 covariant momentum in, (possibly rescaled) momentum out
//! \param[in,out] etot     conserved total energy density
//! \param[out]    v1,v2,v3 CONTRAVARIANT primitive velocity
//! \param[out]    eint     internal energy density, after the deferred floors
//! \param[out]    ceil_used, floored  whether the ceiling / any floor fired
//! \param[out]    de       energy CREATED by the floors (negative if removed)

KOKKOS_INLINE_FUNCTION
void GnomonicRaiseVelFloors(const Real c, const EOS_Data &eos_, const bool gen_,
    const bool keepv_, const bool keept_, const Real vceil_, const bool act_,
    const Real fv, const Real twarm, const Real d,
    Real &m1, Real &m2, Real &m3, Real &etot,
    Real &v1, Real &v2, Real &v3, Real &eint,
    Real &pnew, Real &g1new, Real &temp,
    bool &ceil_used, bool &floored, Real &de) {
  ceil_used = false;
  floored = false;
  de = 0.0;
  pnew = 0.0;
  g1new = 0.0;
  temp = -1.0;
  const Real det = 1.0 - c*c;
  // v^i = g^{ij} m_j / rho, with the metric acting on the ANGULAR pair only
  v1 = m1/d;
  v2 = (m2 - c*m3)/(d*det);
  v3 = (m3 - c*m2)/(d*det);
  // KE = 0.5 rho g_ij v^i v^j = 0.5 m_i v^i, which is the cross-term-correct form.
  Real ekin = 0.5*(m1*v1 + m2*v2 + m3*v3);
  // <hydro>/dfloor_keep_velocity: ConsToPrim scaled this cell's momentum by
  // fv = d_old/dfloor and raised rho to dfloor, so the kinetic energy it now carries is
  // fv^3 of what it had.  The matching reduction of the TOTAL energy could not be made
  // there -- the orthonormal e_k in the c2p is not the kinetic energy on this grid --
  // so it is made here from the metric-correct ekin: KE_old = ekin/fv^3, and removing
  // KE_old - ekin leaves the internal energy exactly where it was.
  if (keepv_) {
    if (fv > 0.0 && fv < 1.0) {
      etot -= ekin*(1.0/(fv*fv*fv) - 1.0);
    } else if (!(fv > 0.0)) {
      etot -= ekin;     // the momentum was zeroed: remove all of it
    }
  }
  // <hydro>/vceil, the deferred half.  |v|^2 = g_ij v^i v^j = 2 KE/rho with the
  // metric-correct kinetic energy formed above, which is the quantity ConsToPrim
  // cannot build on this grid.  Scale the momentum by fs = vceil/|v| and remove
  // (1 - fs^2) KE from the conserved total, so the INTERNAL energy below is exactly
  // what it would have been.
  if (vceil_ > 0.0 && act_ && ekin > 0.0 && d > 0.0) {
    const Real vsq = 2.0*ekin/d;
    if (vsq > vceil_*vceil_) {
      const Real fs = vceil_/sqrt(vsq);
      m1 *= fs; m2 *= fs; m3 *= fs;
      v1 *= fs; v2 *= fs; v3 *= fs;
      etot -= (1.0 - fs*fs)*ekin;
      ekin *= fs*fs;
      ceil_used = true;
    }
  }
  eint = etot - ekin;
  // <hydro>/dfloor_keep_temperature, the deferred half.  ConsToPrim raised rho to
  // dfloor from d_old = fv*dfloor and left the internal energy alone; scale it by fv
  // here, once the metric kinetic energy (and the ceiling above) have fixed what the
  // internal energy actually is.  fv = 1 in every cell the floor did not touch, and
  // fv = 0 marks d_old <= 0, which keeps the old behaviour.  The energy floors below
  // then run on the corrected value.
  if (keept_) {
    if (fv > 0.0 && fv < 1.0) {
      eint *= fv;
      etot = eint + ekin;
    }
  }
  bool mom_scaled = false;
  // -----------------------------------------------------------------------------------
  // RE-APPLY THE FLOORS.  ConsToPrim floored the state it inverted, but that state
  // carried an ORTHONORMAL kinetic energy; the metric cross term above MOVES the
  // internal energy, so a cell ConsToPrim left comfortably above the floor can land
  // below it -- or below zero -- here.  Leaving that unfloored hands a non-positive
  // internal energy straight to the tabulated inversion below, which is undefined
  // there and returns NaN in T, p and Gamma_1; the NaN then leaves the cell through
  // the reconstruction stencil and takes the whole grid down within ~100 cycles.
  // The sequence deliberately mirrors SingleC2P_GeneralHyd, including its guard
  // against inverting a non-positive energy.  NOTE the pressure floor alone is NOT
  // sufficient under a tabulated EOS: at upper-atmosphere densities e(d,pfloor) lies
  // far below the table's lowest temperature, so it is the TEMPERATURE floor that
  // actually keeps the lookup in range.  The conserved state is updated in step,
  // exactly as ConsToPrim updates cons when a floor fires, so the conserved energy
  // cannot keep sinking and re-trip the floor on every cycle.
  if (gen_) {
    bool stale = false;
    const bool e_positive = (eint > 0.0);
    stale = !e_positive;
    if (e_positive) {
      eos_.TemperaturePressureGamma1(d, eint, twarm, temp, pnew, g1new);
    }
    if (!e_positive || pnew < eos_.pfloor) {
      const Real efl = eos_.EnergyFromPressure(d, eos_.pfloor, temp);
      // See EOS_Data::efloor_from_ekin: rebuilding the conserved energy as efl + ekin
      // with ekin untouched donates efl - eint, which is the cell's whole kinetic
      // energy whenever the update left eint negative.  Paying out of the kinetic
      // energy instead holds the conserved total fixed.
      const Real etot0 = etot;
      if (eos_.efloor_from_ekin && ekin > 0.0 && (etot0 - efl) < ekin) {
        Real ek_new = etot0 - efl;
        if (!(ek_new > 0.0)) ek_new = 0.0;
        const Real fvk = sqrt(ek_new/ekin);
        m1 *= fvk; m2 *= fvk; m3 *= fvk;
        v1 *= fvk; v2 *= fvk; v3 *= fvk;
        de += (efl + ek_new) - etot0;
        ekin = ek_new;
        mom_scaled = true;
      } else {
        de += efl - eint;
      }
      eint = efl;
      stale = true;
    }
    if (temp < eos_.tfloor) {
      const Real etf = eos_.EnergyFromTemperature(d, eos_.tfloor);
      de += etf - eint;
      eint = etf;
      temp = eos_.tfloor;
      stale = true;
    }
    if (stale) {
      eos_.PressureAndGamma1(d, eint, temp, pnew, g1new);
      etot = eint + ekin;
      floored = true;
    }
    if (mom_scaled) { floored = true; }
  } else {
    const Real eold = eint;
    eos_.ApplyEnergyFloor(d, eint);
    if (eint != eold) {
      de += eint - eold;
      etot = eint + ekin;
      floored = true;
    }
  }
  return;
}
#endif // COORDINATES_GNOMONIC_RAISEVEL_HPP_
