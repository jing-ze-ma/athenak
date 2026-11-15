#ifndef COORDINATES_GNOMONIC_RAISEVEL_MHD_HPP_
#define COORDINATES_GNOMONIC_RAISEVEL_MHD_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file gnomonic_raisevel_mhd.hpp
//! \brief Per-cell metric part of the MHD cons-to-prim inversion on the cubed sphere,
//! used by Coordinates::GnomonicEquiangleRaiseVelMHD and by the MHD FOFC floor test.
//! The caller forms the cell-centred field, loads d,m1,m2,m3,etot and writes back
//! m1,m2,m3,etot when the ceiling or a floor fired.

#include "athena.hpp"
#include "eos/eos.hpp"
#include "eos/general_c2p_mhd.hpp"
#include "eos/c2p_floors_ideal.hpp"

//----------------------------------------------------------------------------------------
//! \fn void GnomonicRaiseVelMHDFloors
//! \brief Metric velocity raise plus the deferred ceiling and floors for one MHD cell.
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
  // contravariant velocity v^i = g^{ij} m_j / rho (metric acts on the angular pair)
  v1 = m1/d;
  v2 = (m2 - c*m3)/(d*det);
  v3 = (m3 - c*m2)/(d*det);
  // kinetic and magnetic energies are kept separate and subtracted in this order
  Real ekin = 0.5*(m1*v1 + m2*v2 + m3*v3);
  // magnetic energy in the orthonormal frame
  const Real emag = 0.5*(bx*bx + by*by + bz*bz);
  // velocity ceiling: scale the momentum by fs = vceil/|v| with |v|^2 = 2 KE/rho and
  // remove (1 - fs^2) KE from the total energy
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
  // re-apply the floors to the internal energy corrected by the metric kinetic energy
  if (gen_) {
    // general EOS: energy/temperature floors of GeneralMHDFloorsPost (no entropy floor)
    // share3-1031: the floors of SingleC2P_GeneralMHD (GeneralMHDFloorsPost) on the
    // metric-correct energy; own = true; the entropy floor is not applied on this path
    // (sflr = false, as before); the table clamp is not counted here (as before)
    bool efl_used = false, tfl_used = false, tclamped = false, stale = false;
    GeneralMHDFloorsPost(eos_, d, 1.0/d, twarm, true, false, ekin, emag, eint, etot, temp,
                         pnew, g1new, efl_used, tfl_used, tclamped, stale);
    if (stale) {
      etot = eint + ekin + emag;
      floored = true;
    }
  } else {
    // ideal EOS: pressure floor of IdealFloorsPost
    // share3-1031: the ideal floors of SingleC2P_IdealMHD (IdealFloorsPost); on this path
    // only the pressure floor (tflr = sflr = false, as before)
    bool efl_used = false, tfl_used = false;
    IdealFloorsPost<true>(eos_, d, 1.0/d, true, false, false, ekin, emag, eint, etot,
                          efl_used, tfl_used);
    if (efl_used) {floored = true;}
  }
  return;
}
#endif // COORDINATES_GNOMONIC_RAISEVEL_MHD_HPP_
