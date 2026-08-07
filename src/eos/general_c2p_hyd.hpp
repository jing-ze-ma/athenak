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

//----------------------------------------------------------------------------------------
//! \fn void SingleC2P_GeneralHyd()
//! \brief Converts a single state of conserved variables into primitive variables for
//! non-relativistic hydrodynamics with a general EOS, and evaluates the derived
//! thermodynamic quantities that are subsequently reconstructed to interfaces.
//! Conserved = (d,M1,M2,M3,E), Primitive = (d,vx,vy,vz,e)
//! where E=total energy density and e=internal energy density.

KOKKOS_INLINE_FUNCTION
void SingleC2P_GeneralHyd(HydCons1D &u, const EOS_Data &eos, HydPrim1D &w,
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

  // set internal energy, apply floor, correct total energy (if needed)
  Real e_k = 0.5*di*(SQR(u.mx) + SQR(u.my) + SQR(u.mz));
  w.e = (u.e - e_k);

  // apply pressure floor. Unlike the ideal-gas case this cannot be hoisted out of the
  // loop as a constant: e(p) depends on density for a general EOS.
  Real efloor = eos.EnergyFromPressure(w.d, eos.pfloor);
  if (w.e < efloor) {
    w.e = efloor;
    u.e = efloor + e_k;
    efloor_used = true;
  }

  // apply temperature floor, likewise density dependent
  Real e_tfloor = eos.EnergyFromTemperature(w.d, eos.tfloor);
  if (w.e < e_tfloor) {
    w.e = e_tfloor;
    u.e = w.e + e_k;
    tfloor_used = true;
  }

  // TODO(stage3): the ideal-gas path additionally applies an entropy floor using the
  // ideal-gas specific entropy s = (gamma-1)*e/(d*d^(gamma-1)). A general EOS needs its
  // own entropy function, so that floor is not applied here. With the default sfloor of
  // FLT_MIN it never triggers; guard against a user setting it explicitly.

  // Evaluate the derived thermodynamic quantities. For a general EOS this is the single
  // expensive step (Temperature() is a root find), which is precisely why it is done once
  // per cell here rather than at every interface inside the Riemann solvers.
  pgas = eos.Pressure(w.d, w.e);
  g1 = eos.Gamma1(w.d, w.e);
  return;
}

#endif // EOS_GENERAL_C2P_HYD_HPP_
