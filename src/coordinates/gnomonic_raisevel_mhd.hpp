#ifndef COORDINATES_GNOMONIC_RAISEVEL_MHD_HPP_
#define COORDINATES_GNOMONIC_RAISEVEL_MHD_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file gnomonic_raisevel_mhd.hpp
//! \brief the per-cell body of Coordinates::GnomonicEquiangleRaiseVelMHD, as a free
//! device function.
//!
//! WHY IT LIVES HERE.  Exactly the reason gnomonic_raisevel.hpp exists for hydro: on the
//! cubed sphere the conserved-to-primitive inversion is split in two, and the half that
//! needs the metric -- raising the covariant momentum, forming the metric kinetic energy
//! and the ORTHONORMAL magnetic energy, applying the velocity CEILING and re-applying the
//! pressure/temperature/energy floors to the corrected internal energy -- lives in
//! Coordinates, not in ConsToPrim (EOS_Data::defer_cons_floors).  A cell that needs one
//! of those deferred floors is therefore invisible to ConsToPrim, and the FOFC
//! floor-TEST pass, which is nothing but a ConsToPrim over a trial state, never flagged
//! it.  Factoring the body out lets MHD::FOFC run exactly the same arithmetic on its
//! trial state and discard it.
//!
//! EVERYTHING IS LOCAL.  The caller forms the cell-centred field itself (it is rebuilt
//! from the faces, which this function has no business reading), loads d,m1,m2,m3,etot
//! out of its conserved array, and writes m1,m2,m3,etot back only when `ceil_used` or
//! `floored` says something changed -- which is bit-for-bit what the in-line kernel did.

#include "athena.hpp"
#include "eos/eos.hpp"

//----------------------------------------------------------------------------------------
//! \fn void GnomonicRaiseVelMHDFloors
//! \brief metric velocity raise + deferred ceiling/floors for one cubed-sphere MHD cell.
//!
//! \param[in]     c        cos_cell, the angle between the two tangent basis vectors
//! \param[in]     twarm    the temperature ConsToPrim found, used to warm-start the
//!                         general-EOS root find
//! \param[in]     act_     true in ACTIVE cells only; gates the velocity ceiling
//! \param[in]     bx,by,bz the cell-centred field in the ORTHONORMAL frame
//! \param[in,out] m1,m2,m3 covariant momentum in, (possibly rescaled) momentum out
//! \param[in,out] etot     conserved total energy density
//! \param[out]    v1,v2,v3 CONTRAVARIANT primitive velocity
//! \param[out]    eint     internal energy density, after the deferred floors
//! \param[out]    ceil_used, floored  whether the ceiling / any floor fired

KOKKOS_INLINE_FUNCTION
void GnomonicRaiseVelMHDFloors(const Real c, const EOS_Data &eos_, const bool gen_,
    const Real vceil_, const bool act_, const Real twarm, const Real d,
    const Real bx, const Real by, const Real bz,
    Real &m1, Real &m2, Real &m3, Real &etot,
    Real &v1, Real &v2, Real &v3, Real &eint,
    Real &pnew, Real &g1new, Real &temp,
    bool &ceil_used, bool &floored) {
  ceil_used = false;
  floored = false;
  pnew = 0.0;
  g1new = 0.0;
  temp = -1.0;
  const Real det = 1.0 - c*c;
  // v^i = g^{ij} m_j / rho, with the metric acting on the ANGULAR pair only
  v1 = m1/d;
  v2 = (m2 - c*m3)/(d*det);
  v3 = (m3 - c*m2)/(d*det);
  // Keep the two energies separate and subtract them in this order: folding them into a
  // single (kinetic + magnetic) sum re-associates the rounding and perturbs every
  // ideal-EOS cubed-sphere answer in the last bits, for nothing.
  Real ekin = 0.5*(m1*v1 + m2*v2 + m3*v3);
  // the frame is orthonormal, so the magnetic energy IS the sum of squares
  const Real emag = 0.5*(bx*bx + by*by + bz*bz);
  // <mhd>/vceil, the deferred half.  |v|^2 = g_ij v^i v^j = 2 KE/rho with the
  // metric-correct kinetic energy formed above, which is the quantity ConsToPrim cannot
  // build on this grid.  Scale the momentum by fs = vceil/|v| and remove (1 - fs^2) KE
  // from the conserved total; the MAGNETIC energy is not involved, so the internal
  // energy below is exactly what it would have been -- the same bookkeeping the
  // Cartesian MHD c2p does in ideal_c2p_mhd.hpp.
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
  eint = etot - ekin - emag;
  // Re-apply the floors to the CORRECTED internal energy; see the extended note in
  // gnomonic_raisevel.hpp.  Without this the tabulated inversion below is handed a
  // non-positive energy and returns NaN.
  if (gen_) {
    const bool e_positive = (eint > 0.0);
    bool stale = !e_positive;
    if (e_positive) {
      eos_.TemperaturePressureGamma1(d, eint, twarm, temp, pnew, g1new);
    }
    if (!e_positive || pnew < eos_.pfloor) {
      eint = eos_.EnergyFromPressure(d, eos_.pfloor, temp);
      stale = true;
    }
    if (temp < eos_.tfloor) {
      eint = eos_.EnergyFromTemperature(d, eos_.tfloor);
      temp = eos_.tfloor;
      stale = true;
    }
    if (stale) {
      eos_.PressureAndGamma1(d, eint, temp, pnew, g1new);
      etot = eint + ekin + emag;
      floored = true;
    }
  } else {
    const Real eold = eint;
    eos_.ApplyEnergyFloor(d, eint);
    if (eint != eold) {
      etot = eint + ekin + emag;
      floored = true;
    }
  }
  return;
}
#endif // COORDINATES_GNOMONIC_RAISEVEL_MHD_HPP_
