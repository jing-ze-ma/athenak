#ifndef COORDINATES_GNOMONIC_RAISEVEL_HPP_
#define COORDINATES_GNOMONIC_RAISEVEL_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file gnomonic_raisevel.hpp
//! \brief Per-cell metric part of the cons-to-prim inversion on the cubed sphere, used by
//! Coordinates::GnomonicEquiangleRaiseVel and by the FOFC floor test. It raises the
//! covariant momentum, forms the metric kinetic energy, applies the velocity ceiling and
//! re-applies the energy/temperature floors to the corrected internal energy.

#include "athena.hpp"
#include "eos/eos.hpp"
#include "eos/general_c2p_hyd.hpp"
#include "eos/c2p_floors_ideal.hpp"

//----------------------------------------------------------------------------------------
//! \fn void GnomonicRaiseVelFloors
//! \brief Metric velocity raise plus the deferred ceiling and floors for one cell. The
//! caller loads d,m1,m2,m3,etot from the conserved array and writes m1,m2,m3,etot back.
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
//! \param[out]    tset_used  whether <hydro>/efloor_as_tfloor rebuilt e from T
//! \param[out]    de       energy CREATED by the floors (negative if removed)
//! \param[out]    vde      kinetic energy the VELOCITY CEILING removed from the
//!                         conserved total

KOKKOS_INLINE_FUNCTION
void GnomonicRaiseVelFloors(const Real c, const EOS_Data &eos_, const bool gen_,
    const bool keepv_, const Real vceil_, const bool act_,
    const Real fv, const Real twarm, const Real d,
    Real &m1, Real &m2, Real &m3, Real &etot,
    Real &v1, Real &v2, Real &v3, Real &eint,
    Real &pnew, Real &g1new, Real &temp,
    bool &ceil_used, bool &floored, bool &tset_used, Real &de, Real &vde) {
  ceil_used = false;
  floored = false;
  tset_used = false;
  de = 0.0;
  vde = 0.0;
  pnew = 0.0;
  g1new = 0.0;
  temp = -1.0;
  const Real det = 1.0 - c*c;
  // contravariant velocity v^i = g^{ij} m_j / rho (metric acts on the angular pair)
  v1 = m1/d;
  v2 = (m2 - c*m3)/(d*det);
  v3 = (m3 - c*m2)/(d*det);
  // kinetic energy 0.5 rho g_ij v^i v^j = 0.5 m_i v^i
  Real ekin = 0.5*(m1*v1 + m2*v2 + m3*v3);
  // dfloor_keep_velocity: ConsToPrim scaled the momentum by fv = d_old/dfloor, so remove
  // the kinetic energy lost, KE_old - ekin = ekin*(1/fv^3 - 1), from the total energy
  if (keepv_) {
    if (fv > 0.0 && fv < 1.0) {
      etot -= ekin*(1.0/(fv*fv*fv) - 1.0);
    } else if (!(fv > 0.0)) {
      etot -= ekin;     // momentum was zeroed
    }
  }
  // velocity ceiling: scale the momentum by fs = vceil/|v| with |v|^2 = 2 KE/rho, and
  // remove (1 - fs^2) KE from the total energy
  if (vceil_ > 0.0 && act_ && ekin > 0.0 && d > 0.0) {
    const Real vsq = 2.0*ekin/d;
    if (vsq > vceil_*vceil_) {
      const Real fs = vceil_/sqrt(vsq);
      m1 *= fs; m2 *= fs; m3 *= fs;
      v1 *= fs; v2 *= fs; v3 *= fs;
      etot -= (1.0 - fs*fs)*ekin;
      vde += (1.0 - fs*fs)*ekin;
      ekin *= fs*fs;
      ceil_used = true;
    }
  }
  eint = etot - ekin;
  bool mom_scaled = false;
  // Re-apply the floors to the internal energy corrected by the metric kinetic energy,
  // updating the conserved state consistently (as in SingleC2P_GeneralHyd/IdealHyd).
  if (gen_) {
    // general EOS: energy/temperature floors of GeneralHydFloorsPost (no entropy floor)
    bool tclamped = false, efl_used = false, tfl_used = false, stale = false;
    bool floor_fired = false;
    Real fvk = 1.0;
    GeneralHydFloorsPost(eos_, d, 1.0/d, twarm, true, false, eint, ekin, etot, fvk, temp,
                         pnew, g1new, efl_used, tfl_used, tclamped, tset_used,
                         mom_scaled, floor_fired, stale, de);
    if (mom_scaled) {
      m1 *= fvk; m2 *= fvk; m3 *= fvk;
      v1 *= fvk; v2 *= fvk; v3 *= fvk;
    }
    if (stale) {
      etot = eint + ekin;
      floored = floor_fired;
    }
    if (mom_scaled) { floored = true; }
  } else {
    // ideal EOS: pressure floor of IdealFloorsPost
    const Real eold = eint;
    bool efl_used = false, tfl_used = false;
    IdealFloorsPost<false>(eos_, d, 1.0/d, true, false, false, ekin, 0.0, eint, etot,
                           efl_used, tfl_used);
    if (efl_used) {
      de += eint - eold;
      floored = true;
    }
  }
  return;
}
#endif // COORDINATES_GNOMONIC_RAISEVEL_HPP_
