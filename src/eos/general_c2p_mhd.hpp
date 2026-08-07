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
//! the arguments does not.

KOKKOS_INLINE_FUNCTION
void SingleC2P_GeneralMHD(MHDCons1D &u, const EOS_Data &eos, HydPrim1D &w,
                          Real &pgas, Real &g1,
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

  // apply pressure floor. Density dependent for a general EOS, so unlike the ideal-gas
  // path it cannot be precomputed as a constant.
  Real efloor = eos.EnergyFromPressure(w.d, eos.pfloor);
  if (w.e < efloor) {
    w.e = efloor;
    u.e = efloor + e_k + e_m;
    efloor_used = true;
  }

  // apply temperature floor, likewise density dependent
  Real e_tfloor = eos.EnergyFromTemperature(w.d, eos.tfloor);
  if (w.e < e_tfloor) {
    w.e = e_tfloor;
    u.e = w.e + e_k + e_m;
    tfloor_used = true;
  }

  // TODO(stage3): the ideal-gas path additionally applies an entropy floor, which needs
  // a general-EOS entropy function. Not applied here; with the default sfloor of FLT_MIN
  // it never triggers.

  // Evaluate the derived thermodynamic quantities once per cell (see general_c2p_hyd.hpp)
  pgas = eos.Pressure(w.d, w.e);
  g1 = eos.Gamma1(w.d, w.e);
  return;
}

#endif // EOS_GENERAL_C2P_MHD_HPP_
