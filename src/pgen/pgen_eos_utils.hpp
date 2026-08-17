#ifndef PGEN_PGEN_EOS_UTILS_HPP_
#define PGEN_PGEN_EOS_UTILS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file pgen_eos_utils.hpp
//! \brief small EOS-aware conversions shared by the gravitationally stratified problem
//! generators (hydrostatic atmospheres, convection, hot Jupiters, hot bubbles).
//!
//! These problem generators were written against an ideal gas and convert between
//! pressure, internal energy and temperature inline: `p*igm1`, `e*gm1`, `p/(Rgas d)`.
//! Every one of those is an ideal-gas identity. Under a general EOS the primitive stored
//! in w0(IEN) is still the internal energy, but pressure is no longer proportional to it
//! and the temperature no longer follows from a fixed mean molecular weight.
//!
//! Each helper keeps the ORIGINAL ideal-gas arithmetic on the ideal branch -- not the
//! EOS's algebraically equal form -- so that existing ideal-gas runs are reproduced bit
//! for bit. Multiplying by a reciprocal and dividing differ in the last bit, and these
//! problems are well-balanced ones where the velocity is the residual of a near
//! cancellation, so last bits are visible.

#include "athena.hpp"
#include "eos/eos.hpp"

namespace pgen_eos {

//----------------------------------------------------------------------------------------
//! \fn Real EintFromP
//! \brief internal energy density from pressure.
KOKKOS_INLINE_FUNCTION
Real EintFromP(const EOS_Data &eos, const Real igm1, const Real d, const Real p) {
  return (eos.IsGeneral()) ? eos.EnergyFromPressure(d, p) : p*igm1;
}

//----------------------------------------------------------------------------------------
//! \fn Real PresFromEint
//! \brief pressure from the primitive pair (d,e).
//!
//! NOTE the general branch is `Pressure(d,e)`, which SOLVES for the temperature. If the
//! caller also wants T -- radiative transfer and cooling always do -- use
//! PresTempFromEint below instead, which gets both from one solve, warm started.
KOKKOS_INLINE_FUNCTION
Real PresFromEint(const EOS_Data &eos, const Real gm1, const Real d, const Real e) {
  return (eos.IsGeneral()) ? eos.Pressure(d, e) : e*gm1;
}

//----------------------------------------------------------------------------------------
//! \fn void PresTempFromEint
//! \brief pressure AND kelvin temperature from (d,e), with ONE temperature solve.
//!
//! The obvious spelling of this -- `PresFromEint` followed by `TempKelvin` -- costs the
//! general EOS TWO cold-start root finds per call, because each of `Pressure(d,e)` and
//! `Temperature(d,e)` inverts e(d,T) independently and neither is given a starting
//! point. Solving once and passing the temperature to `Pressure(d,e,t)` removes one
//! inversion outright, and `tguess` (normally Hydro/MHD::wtemp for the cell, which
//! ConsToPrim has already solved for this very state) cuts the survivor to a couple of
//! Newton steps. A non-positive or stale guess is safe: SolveTemperature keeps its full
//! bracket and converges to the same root regardless.
KOKKOS_INLINE_FUNCTION
void PresTempFromEint(const EOS_Data &eos, const Real gm1, const Real Rgas, const Real d,
                      const Real e, const Real tguess, Real &p, Real &tk) {
  if (eos.IsGeneral()) {
    Real tc = eos.Temperature(d, e, tguess);
    p  = eos.Pressure(d, e, tc);
    tk = tc*eos.temp_cgs;
  } else {
    p  = e*gm1;
    tk = p/Rgas/d;
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real TempKelvin
//! \brief temperature in KELVIN, which is what cooling functions and radiative transfer
//! need. The ideal branch is the problem's own p/(Rgas d), with Rgas carrying the fixed
//! mean molecular weight; the general branch asks the EOS, where composition lives
//! instead, and converts with eos.temp_cgs. Pressure is passed in because call sites have
//! normally evaluated it already.
//!
//! `tguess` is a CODE-temperature warm start for the general branch, normally
//! Hydro/MHD::wtemp for the same cell -- the temperature ConsToPrim already solved for
//! this very state. Without it every call brackets from scratch, which inside a per-cell
//! radiative-transfer sweep is the difference between two Newton steps and dozens. A
//! non-positive value means "no guess"; a stale one costs iterations but never accuracy,
//! since the solver converges to the same root either way. The ideal branch ignores it.
KOKKOS_INLINE_FUNCTION
Real TempKelvin(const EOS_Data &eos, const Real Rgas, const Real d, const Real e,
                const Real p, const Real tguess = -1.0) {
  return (eos.IsGeneral()) ? eos.Temperature(d, e, tguess)*eos.temp_cgs : p/Rgas/d;
}

//----------------------------------------------------------------------------------------
//! \fn Real DensFromPT
//! \brief density from pressure and a KELVIN temperature. Hydrostatic backgrounds are
//! built by marching in (p,T), which is the one direction the general EOS interface
//! cannot evaluate from the primitive pair -- hence DensityFromPressureTemperature, which
//! root finds on density. Setup only; never called inside a time step.
KOKKOS_INLINE_FUNCTION
Real DensFromPT(const EOS_Data &eos, const Real Rgas, const Real p, const Real tk) {
  return (eos.IsGeneral()) ? eos.DensityFromPressureTemperature(p, tk/eos.temp_cgs)
                           : p/Rgas/tk;
}

//----------------------------------------------------------------------------------------
//! \fn Real GradAd
//! \brief adiabatic temperature gradient (dln T/dln p)_s at a point specified by pressure
//! and a KELVIN temperature. grad_ad = (Gamma_3-1)/Gamma_1 with Gamma_3-1 =
//! p chi_T/(d T c_v), which collapses to (gamma-1)/gamma for an ideal gas. This is the
//! one place in a stratified initial condition where a real EOS changes the answer
//! qualitatively: grad_ad drops sharply through the H2 dissociation and H ionization
//! zones, which is what sets where the atmosphere is convectively unstable.
KOKKOS_INLINE_FUNCTION
Real GradAd(const EOS_Data &eos, const Real gamma, const Real Rgas, const Real p,
            const Real tk) {
  if (!eos.IsGeneral()) return ((gamma-1.0)/gamma);
  Real tc = tk/eos.temp_cgs;
  Real d = eos.DensityFromPressureTemperature(p, tc);
  Real e = eos.EnergyFromTemperature(d, tc);
  Real g3m1 = p*eos.ChiT(d, e)/(d*tc*eos.SpecificHeatCv(d, e));
  return (g3m1/eos.Gamma1(d, e));
}

//----------------------------------------------------------------------------------------
//! \fn Real EintFromDensT
//! \brief internal energy density from density and a KELVIN temperature, for atmospheres
//! specified by a temperature profile rather than a pressure one.
KOKKOS_INLINE_FUNCTION
Real EintFromDensT(const EOS_Data &eos, const Real Rgas, const Real igm1, const Real d,
                   const Real tk) {
  return (eos.IsGeneral()) ? eos.EnergyFromTemperature(d, tk/eos.temp_cgs)
                           : d*Rgas*tk*igm1;
}

//----------------------------------------------------------------------------------------
//! \fn Real HostGamma1FromP
//! \brief Gamma_1 at the state (d,p), evaluated where the EOS lives.
//!
//! The helpers above are all KOKKOS_INLINE_FUNCTIONs meant to be called from inside a
//! kernel. That is not a style preference: the tabulated EOS keeps its table in a
//! DvceArray, so on an accelerator build it CANNOT be evaluated from host code. A problem
//! generator that needs one EOS number on the host -- the sound speed of a background
//! state, to build a wave eigenvector -- has to fetch it through the device, which is
//! what this does. It launches a one-element kernel, so it is for setup only; never call
//! it in a loop.
inline Real HostGamma1FromP(const EOS_Data &eos, const Real d, const Real p) {
  if (!eos.IsGeneral()) return eos.gamma;
  DvceArray1D<Real> dout("pgen_eos_scalar", 1);
  par_for("pgen_eos_eval", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(int) {
    dout(0) = eos.Gamma1(d, eos.EnergyFromPressure(d, p));
  });
  auto hout = Kokkos::create_mirror_view(dout);
  Kokkos::deep_copy(hout, dout);
  return hout(0);
}

}  // namespace pgen_eos

#endif  // PGEN_PGEN_EOS_UTILS_HPP_
