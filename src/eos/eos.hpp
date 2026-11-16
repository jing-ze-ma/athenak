#ifndef EOS_EOS_HPP_
#define EOS_EOS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file eos.hpp
//! \brief Contains data and functions that implement conserved->primitive variable
//! conversion for various EOS (e.g. ideal gas, isothermal, etc.), for various fluids
//! (Hydro, MHD, etc.), and for non-relativistic and relativistic flows.

//#include <cmath>
#include <math.h>
#include <string>

#include "athena.hpp"
#include "mesh/meshblock.hpp"
#include "parameter_input.hpp"
#include "eos/eos_table.hpp"

//----------------------------------------------------------------------------------------
//! \enum EOSType
//! \brief equation of state evaluated by the general EOS interface of EOS_Data; ideal
//! reproduces the ideal-gas expressions exactly

enum class EOSType {ideal, general};

//----------------------------------------------------------------------------------------
//! \struct EOSData
//! \brief container for EOS variables and functions needed inside kernels. Storing
//! everything in a container makes them easier to capture, and pass to inline functions,
//! inside kernels.

struct EOS_Data {
  Real gamma;        // ratio of specific heats for ideal gas
  Real iso_cs;       // isothermal sound speed
  bool is_ideal;     // flag to denote ideal gas EOS
  Real dfloor, pfloor, tfloor, sfloor;  // density, pressure, temperature, entropy floors
  // cubed sphere: ConsToPrim clamps the primitives but leaves the conserved floor writes
  // to Coordinates::GnomonicEquiangleRaiseVel{,MHD}.  The covariant momentum and
  // face-normal B triple do not give the true kinetic and magnetic energies there, so a
  // floor written back from ConsToPrim would freeze that error into u.e
  bool defer_cons_floors = false;
  // <block>/dfloor_keep_velocity: the density floor scales the momentum by d_old/dfloor
  // and removes the kinetic energy lost, so the internal energy is unchanged (default
  // false: only u.d is raised)
  bool dfloor_keep_velocity = false;

  // <hydro|mhd>/vceil: velocity ceiling of the Newtonian inversions (0 = off).  Where
  // |v| > vceil the momentum is scaled by vceil/|v| and the kinetic energy removed from
  // the total, so the internal energy is unchanged.  On the cubed sphere it is applied in
  // Coordinates::GnomonicEquiangleRaiseVel{,MHD}.  Counted in EventCounters::neos_vceil
  Real vceil = 0.0;

  // <mhd>/valfven_max_cgs: Alfven-speed density floor [cm/s] (0 = off).  Where
  // B^2/rho > valfven_max^2 the density is raised to B^2/valfven_max^2 like dfloor
  // (momentum and energy unchanged).  Stored as 1/valfven_max^2 in code units
  Real va2max_inv = 0.0;

  // <mhd>/valfven_noinfall_fac (0 = off): where rho <= fac B^2/valfven_max^2 the inward
  // radial momentum and its kinetic energy are removed, so gas held up by the Alfven
  // floor does not rain into the atmosphere.  Needs a radial x1 (spherical polar or cubed
  // sphere)
  Real va_noinfall_fac = 0.0;

  // <mhd>/floor_tracer: scalar index of the floor-mass tracer (-1 = off); every
  // density-floor increment is added to it
  int ifltr = -1;
  // <mhd>/xe_tracer: scalar index of the electron-fraction tracer (-1 = off); scaled with
  // the density when a floor raises it
  int ixe = -1;

  // <mhd>/hlld_bx_zero_tol: the HLLD family drops the rotational (**) states where
  // 0.5 Bx^2 < tol pt*, i.e. for plasma beta above ~2/tol.  The default 1.0e-4 is the
  // AthenaK/Athena++ value; ~1.0e-8 keeps the ** states at very high beta
  Real hlld_bx_zero_tol = 1.0e-4;

  // <hydro>/lhllc_x1_phi_min, <mhd>/lhlld_x1_phi_min: floor on the low-Mach factor
  // phi = chi (2 - chi) of LHLLC/LHLLD on x1 faces only (radial on curvilinear grids);
  // >= 1 gives full HLLC/HLLD radially, <= 0 = off
  Real x1_phi_min = 0.0;

  // <mhd>/lhlld_emf_eta_pos: non-negative numerical resistivity of the LHLLD
  // transverse-field flux.  The upwind part of the HLLD B_t flux is a resistivity
  // eta_face = (|S_M| + du/2) dx/2, negative on a strongly converging face.  With the
  // switch on, dF_Bt = -(1/2) max(0, -(|S_M| + du/2)) dB_t is added there, with the
  // matching Poynting term in the energy flux (the removed magnetic energy becomes heat)
  bool lhlld_emf_eta_pos = false;
  // <hydro>/lhllc_cb_sensor, <mhd>/lhlld_cb_sensor: odd-even (2-cell) sensor of
  // LHLLC/LHLLD (<= 0 = off).  From the normal mass flux q = rho u_n at cells i-2..i+1
  // the odd part of the face jump is J = (2 d2 - d1 - d3)/4; s = clamp(J/d2, 0, 1)^p,
  // p = cb_sensor, and phi -> phi + (1 - phi) s, so the face dissipation lies between
  // LHLLC/LHLLD (s = 0) and HLLC/HLLD (s = 1) (utils/cb_sensor.hpp)
  Real cb_sensor = 0.0;
  // <hydro>/lhllc_cb_cont, <mhd>/lhlld_cb_cont: continuous sensor fraction
  // s = clamp(J d2/max(d2^2, J^2), 0, 1); equal to clamp(J/d2, 0, 1) where |J| <= |d2|
  bool cb_cont = false;
  // <hydro|mhd>/wb_plm_smooth: continuous blend of the well-balanced and plain x1 PLM
  // faces (reconstruct/wb_plm_blend.hpp)
  bool wb_plm_smooth = false;
  // <hydro|mhd>/wb_plm_range: bound the well-balanced x1 faces by the 3-cell data range
  // (reconstruct/wb_plm_blend.hpp)
  bool wb_plm_range = false;
  // <mhd>/lhlld_chi_normal: c_u of LHLLD from the face-normal velocity only
  // (|u|^2 -> u_n^2), as in LHLLC; also used by the UCT face store
  bool chi_normal = false;
  // <mhd>/lhlld_alfmach: Alfven-Mach scaled floor of the LHLLD c_u.  The Alfven speeds
  // are scaled by min(1, |u|/c_a), so chi = O(|u|/c_f) for M_A < 1; unchanged for
  // M_A >= 1
  bool alfmach = false;

  // <block>/efloor_from_ekin: the energy floor is paid from the kinetic energy (total
  // energy fixed, momentum rescaled); energy is created only when E < e_floor
  bool efloor_from_ekin = false;

  // <block>/efloor_as_tfloor: whenever a floor fires or the table inversion is clamped,
  // the internal energy is reset to e(rho, T) at the temperature settled on (at least
  // tfloor) and the conserved energy corrected to match, so the floored cell is a
  // consistent table state.  ReadEOS_Params raises tfloor to the table's lowest
  // temperature.  Counted in EventCounters::neos_tset
  bool efloor_as_tfloor = false;

  // true when none of dfloor_keep_velocity, vceil, efloor_from_ekin, efloor_as_tfloor is
  // on; kernels whose floor code was regrouped for those switches then run their original
  // form, which keeps default runs bitwise identical (set in ReadEOS_Params)
  bool floors_legacy = true;
  Real gamma_max;    // ceiling on Lorentz factor in SR/GR
  // <hydro>/ausm_mcut, ausm_mcut_p: AUSM+-up cut-off Mach numbers of the velocity- and
  // pressure-diffusion terms (Liou 2006; mcut_p = 1 removes 1/f_a from M_p,
  // Edelmann+2021)
  Real ausm_mcut = 1.0e-13, ausm_mcut_p = 1.0;
  // ideal-gas run with the composition table built only for the electron fraction
  // (ohmic_resistivity = eos): tbl holds valid surfaces but tbl.active is false, so all
  // thermodynamics stays ideal (EquationOfState::BuildElectronFractionTable)
  bool xe_from_table = false;

  EOSType eos_type = EOSType::ideal;  // EOS evaluated by the general interface below

  // code -> cgs unit factors, used by the general EOS (copied from units::Units, which
  // is host-only)
  Real dens_cgs = 1.0;   // g/cm^3       per unit code density
  Real pres_cgs = 1.0;   // erg/cm^3     per unit code pressure or energy density
  Real temp_cgs = 1.0;   // K            per unit code temperature

  // tabulated EOS; tbl.active = false makes the general interface evaluate the ideal gas
  EOSTable tbl;

  //! \fn bool IsGeneral
  //! \brief true when p and Gamma_1 are precomputed into the derived variables and
  //! reconstructed instead of evaluated inside the Riemann solvers
  KOKKOS_INLINE_FUNCTION
  bool IsGeneral() const {
    return (eos_type == EOSType::general);
  }

  //--------------------------------------------------------------------------------------
  // General EOS interface of non-relativistic hydro and MHD, in terms of the primitive
  // pair (d, e) = (density, internal energy density).  Each function reduces exactly to
  // the ideal-gas expression for EOSType::ideal.  Recovering T from (d, e) is the
  // expensive root find, so most accessors come in a (d, e, T) form for hot kernels (T
  // solved once per cell by ConsToPrim and cached in wtemp) and a (d, e) form for
  // setup-time callers.

  //--------------------------------------------------------------------------------------
  //! \fn Real Temperature
  //! \brief temperature T(d, e) in code units (p/d for an ideal gas); tguess > 0
  //! warm-starts the table inversion
  KOKKOS_INLINE_FUNCTION
  Real Temperature(const Real d, const Real e, const Real tguess = -1.0) const {
    if (tbl.active) {
      Real tg = (tguess > 0.0) ? tguess*temp_cgs : -1.0;
      return tbl.SolveTemperature(d*dens_cgs, e*pres_cgs, tg)/temp_cgs;
    }
    return ((gamma-1.0)*e/d);
  }

  //! \fn Real Pressure
  //! \brief gas pressure p(d,e), from an already-solved temperature
  KOKKOS_INLINE_FUNCTION
  Real Pressure(const Real d, const Real e, const Real t) const {
    if (tbl.active) {
      EOSThermoState s;
      tbl.Eval(d*dens_cgs, t*temp_cgs, s);
      return (s.p/pres_cgs);
    }
    return ((gamma-1.0)*e);
  }

  //! \fn Real Pressure
  //! \brief gas pressure p(d, e); solves for the temperature (setup-time use)
  KOKKOS_INLINE_FUNCTION
  Real Pressure(const Real d, const Real e) const {
    return Pressure(d, e, Temperature(d,e));
  }

  //! \fn Real Gamma1
  //! \brief first adiabatic exponent Gamma_1 from an already-solved temperature (gamma
  //! for an ideal gas); replaces gamma in the wave speeds of the Riemann solvers
  KOKKOS_INLINE_FUNCTION
  Real Gamma1(const Real d, const Real e, const Real t) const {
    if (tbl.active) {
      EOSThermoState s;
      Real rho = d*dens_cgs;
      Real tk = t*temp_cgs;
      tbl.Eval(rho, tk, s);
      // Gamma_1 = chi_rho + p chi_T^2/(rho T c_v), all from one table evaluation
      return (s.chi_rho + s.p*s.chi_t*s.chi_t/(rho*tk*s.cv));
    }
    return gamma;
  }

  //! \fn Real Gamma1
  //! \brief Gamma_1(d, e); solves for the temperature (setup-time use)
  KOKKOS_INLINE_FUNCTION
  Real Gamma1(const Real d, const Real e) const {
    return Gamma1(d, e, Temperature(d,e));
  }

  //! \fn Real EnergyFromPressure
  //! \brief internal energy density e(d, p), the inverse of Pressure (a root find for a
  //! general EOS: call it only where a floor has been shown to trip)
  KOKKOS_INLINE_FUNCTION
  Real EnergyFromPressure(const Real d, const Real p) const {
    if (tbl.active) {
      Real rho = d*dens_cgs;
      Real tk = tbl.SolveTemperatureFromP(rho, p*pres_cgs, -1.0);
      EOSThermoState s;
      tbl.Eval(rho, tk, s);
      return (s.e/pres_cgs);
    }
    return (p/(gamma-1.0));
  }

  //! \fn Real EnergyFromPressure
  //! \brief as above, also returning the temperature t of the inversion.  The ideal
  //! branch forms T from e so that it matches a following Temperature(d, e) bitwise
  KOKKOS_INLINE_FUNCTION
  Real EnergyFromPressure(const Real d, const Real p, Real &t) const {
    if (tbl.active) {
      Real rho = d*dens_cgs;
      Real tk = tbl.SolveTemperatureFromP(rho, p*pres_cgs, -1.0);
      EOSThermoState s;
      tbl.Eval(rho, tk, s);
      t = tk/temp_cgs;
      return (s.e/pres_cgs);
    }
    Real e = p/(gamma-1.0);
    t = (gamma-1.0)*e/d;
    return e;
  }

  //! \fn void PressureAndGamma1
  //! \brief pressure and Gamma_1 at an already-solved temperature from one table
  //! evaluation (bitwise equal to the separate calls)
  KOKKOS_INLINE_FUNCTION
  void PressureAndGamma1(const Real d, const Real e, const Real t,
                         Real &p, Real &g1) const {
    if (tbl.active) {
      EOSThermoState s;
      Real rho = d*dens_cgs;
      Real tk = t*temp_cgs;
      tbl.EvalNoMu(rho, tk, s);
      p = s.p/pres_cgs;
      g1 = (s.chi_rho + s.p*s.chi_t*s.chi_t/(rho*tk*s.cv));
      return;
    }
    p = (gamma-1.0)*e;
    g1 = gamma;
    return;
  }

  //! \fn void TemperaturePressureGamma1
  //! \brief the EOS step of ConsToPrim: solve e(d, T) for T, then p and Gamma_1 at that
  //! T, chained in the log domain.  tclamped is true when the inversion was held at the
  //! table's lowest or highest temperature (EOSTable::ClampLogT)
  KOKKOS_INLINE_FUNCTION
  void TemperaturePressureGamma1(const Real d, const Real e, const Real tguess,
                                 Real &t, Real &p, Real &g1, bool &tclamped) const {
    if (tbl.active) {
      const Real rho = d*dens_cgs;
      const Real lrho = log10(rho);
      const Real zg = (tguess > 0.0) ? log10(tguess*temp_cgs) : -1.0e30;
      const Real ltk = tbl.SolveLogTemperature(lrho, log10(e*pres_cgs), zg, tclamped);
      const Real tk = EOSTable::Pow10(ltk);
      EOSThermoState s;
      tbl.EvalFromLogs<false>(lrho, ltk, rho, tk, s);
      t = tk/temp_cgs;
      p = s.p/pres_cgs;
      g1 = (s.chi_rho + s.p*s.chi_t*s.chi_t/(rho*tk*s.cv));
      return;
    }
    t = (gamma-1.0)*e/d;
    p = (gamma-1.0)*e;
    g1 = gamma;
    return;
  }

  //! \fn void TemperaturePressureGamma1
  //! \brief as above, without the table-range clamp flag
  KOKKOS_INLINE_FUNCTION
  void TemperaturePressureGamma1(const Real d, const Real e, const Real tguess,
                                 Real &t, Real &p, Real &g1) const {
    bool tclamped = false;
    TemperaturePressureGamma1(d, e, tguess, t, p, g1, tclamped);
    return;
  }

  //! \fn bool BelowPressureFloor
  //! \brief true when p(d, e) < pfloor at an already-solved temperature; equivalent to
  //! e < EnergyFromPressure(d, pfloor) without its root find
  KOKKOS_INLINE_FUNCTION
  bool BelowPressureFloor(const Real d, const Real e, const Real t) const {
    return (Pressure(d, e, t) < pfloor);
  }

  //! \fn Real TableTempMin
  //! \brief lowest tabulated temperature in code units (-1 without a table)
  KOKKOS_INLINE_FUNCTION
  Real TableTempMin() const {
    return tbl.active ? (EOSTable::Pow10(tbl.ymin)/temp_cgs) : -1.0;
  }
  //! \fn Real TableTempMax
  //! \brief highest tabulated temperature in code units (-1 without a table)
  KOKKOS_INLINE_FUNCTION
  Real TableTempMax() const {
    return tbl.active ? (EOSTable::Pow10(tbl.ymax)/temp_cgs) : -1.0;
  }

  //! \fn Real EnergyFloorBound
  //! \brief cheap upper bound on EnergyFromPressure(d, pfloor) (exact for an ideal gas):
  //! cells with e above it skip the floor's root find
  KOKKOS_INLINE_FUNCTION
  Real EnergyFloorBound(const Real d) const {
    if (tbl.active) {
      return tbl.EnergyFloorBound(d);
    }
    return (pfloor/(gamma-1.0));
  }

  //! \fn void ApplyEnergyFloor
  //! \brief raises e to e(d, pfloor) if it lies below it (reconstructed face states)
  KOKKOS_INLINE_FUNCTION
  void ApplyEnergyFloor(const Real d, Real &e) const {
    if (e < EnergyFloorBound(d)) {
      e = fmax(e, EnergyFromPressure(d, pfloor));
    }
  }

  //! \fn bool ApplyEntropyFloor
  //! \brief raises e so that K = p/d^gamma reaches sfloor; returns whether it did.  di =
  //! 1/d. Gamma-law only: a tabulated EOS returns false (BuildGeneralEOS refuses sfloor
  //! with general_eos = table).  As in the ideal-gas floor, the total energy is not
  //! corrected
  KOKKOS_INLINE_FUNCTION
  bool ApplyEntropyFloor(const Real d, const Real di, Real &e) const {
    if (tbl.active) return false;
    Real gm1 = gamma - 1.0;
    Real spe_over_eps = gm1/pow(d, gm1);
    Real spe = spe_over_eps*e*di;
    if (spe <= sfloor) {
      e = d*sfloor/spe_over_eps;
      return true;
    }
    return false;
  }

  //! \fn void ThermoAt
  //! \brief e, p, chi_rho, chi_T and c_v at (d, T) in code units from one table
  //! evaluation (identical to the single accessors)
  KOKKOS_INLINE_FUNCTION
  void ThermoAt(const Real d, const Real t, Real &e, Real &p, Real &chi_rho,
                Real &chi_t, Real &cv) const {
    if (tbl.active) {
      const Real rho = d*dens_cgs;
      const Real tk = t*temp_cgs;
      EOSThermoState s;
      tbl.EvalNoMu(rho, tk, s);   // mu is not wanted: one surface fewer
      e = s.e/pres_cgs;
      p = s.p/pres_cgs;
      chi_rho = s.chi_rho;
      chi_t = s.chi_t;
      cv = s.cv*dens_cgs*temp_cgs/pres_cgs;
      return;
    }
    e = d*t/(gamma-1.0);
    p = (gamma-1.0)*e;
    chi_rho = 1.0;
    chi_t = 1.0;
    cv = 1.0/(gamma-1.0);
    return;
  }

  //! \fn Real EnergyFromTemperature
  //! \brief internal energy density e(d, T), the inverse of Temperature
  KOKKOS_INLINE_FUNCTION
  Real EnergyFromTemperature(const Real d, const Real t) const {
    if (tbl.active) {
      const Real rho = d*dens_cgs;
      const Real tk = t*temp_cgs;
      EOSThermoState s;
      tbl.Eval(rho, tk, s);
      return (s.e/pres_cgs);
    }
    return (d*t/(gamma-1.0));
  }

  //! \fn Real MeanMolecularWeight
  //! \brief mean molecular weight mu (mass per particle in m_u); 1 for the ideal gas
  //! (mu_ref = 1 convention)
  KOKKOS_INLINE_FUNCTION
  Real MeanMolecularWeight(const Real d, const Real e, const Real t) const {
    if (tbl.active) {
      EOSThermoState s;
      tbl.Eval(d*dens_cgs, t*temp_cgs, s);
      return s.mu;
    }
    return 1.0;
  }

  //! \fn Real MeanMolecularWeight
  //! \brief mu(d, e); solves for the temperature (setup-time use)
  KOKKOS_INLINE_FUNCTION
  Real MeanMolecularWeight(const Real d, const Real e) const {
    return MeanMolecularWeight(d, e, Temperature(d,e));
  }

  //! \fn Real ElectronFraction
  //! \brief free electrons per particle n_e/n_tot from the table composition (Saha); 0
  //! for an ideal gas.  Below ~4000 K the electrons come from the metals
  //! (eos_metal_ionization)
  KOKKOS_INLINE_FUNCTION
  Real ElectronFraction(const Real d, const Real e, const Real t) const {
    if (tbl.active) {
      return tbl.ElectronFractionCgs(d*dens_cgs, t*temp_cgs);
    }
    return 0.0;
  }

  //! \fn bool MetalIonization
  //! \brief were the metal electron donors included when the EOS table was built?
  KOKKOS_INLINE_FUNCTION
  bool MetalIonization() const { return (tbl.active && tbl.metal_ion); }

  //! \fn Real MetalMetallicity
  //! \brief [M/H] in dex that those donors were built with
  KOKKOS_INLINE_FUNCTION
  Real MetalMetallicity() const { return tbl.metal_mh; }

  //! \fn Real ElectronFraction
  //! \brief n_e/n_tot at (d, e); solves for the temperature (setup-time use)
  KOKKOS_INLINE_FUNCTION
  Real ElectronFraction(const Real d, const Real e) const {
    return ElectronFraction(d, e, Temperature(d,e));
  }

  //! \fn bool HasElectronFraction
  //! \brief true when an electron fraction is available (general EOS table, or a table
  //! built only for it on top of an ideal gas)
  KOKKOS_INLINE_FUNCTION
  bool HasElectronFraction() const { return (tbl.active || xe_from_table); }

  //! \fn Real ElectronFractionKelvin
  //! \brief n_e/n_tot from a code density and a temperature in kelvin; the density scale
  //! is the table's, so it also serves ideal-gas runs with xe_from_table
  KOKKOS_INLINE_FUNCTION
  Real ElectronFractionKelvin(const Real d, const Real t_kelvin) const {
    if (tbl.active || xe_from_table) {
      return tbl.ElectronFractionCgs(d*tbl.dens_cgs, t_kelvin);
    }
    return 0.0;
  }

  //! \fn Real DensityFromPressureTemperature
  //! \brief density d(p, T) (a root find for a table; setup-time use, e.g. hydrostatic
  //! backgrounds built in (p, T))
  KOKKOS_INLINE_FUNCTION
  Real DensityFromPressureTemperature(const Real p, const Real t) const {
    if (tbl.active) {
      return tbl.SolveDensity(p*pres_cgs, t*temp_cgs, -1.0)/dens_cgs;
    }
    return (p/t);
  }

  //! \fn Real SpecificHeatCv
  //! \brief specific heat at constant volume c_v = (de/dT)/d
  KOKKOS_INLINE_FUNCTION
  Real SpecificHeatCv(const Real d, const Real e, const Real t) const {
    if (tbl.active) {
      EOSThermoState s;
      tbl.Eval(d*dens_cgs, t*temp_cgs, s);
      return (s.cv*dens_cgs*temp_cgs/pres_cgs);
    }
    return (1.0/(gamma-1.0));
  }

  //! \fn Real SpecificHeatCv
  //! \brief c_v(d, e); solves for the temperature (setup-time use)
  KOKKOS_INLINE_FUNCTION
  Real SpecificHeatCv(const Real d, const Real e) const {
    return SpecificHeatCv(d, e, Temperature(d,e));
  }

  //! \fn Real ChiRho
  //! \brief chi_rho = (dln p/dln d) at constant temperature (not Gamma_1)
  KOKKOS_INLINE_FUNCTION
  Real ChiRho(const Real d, const Real e, const Real t) const {
    if (tbl.active) {
      EOSThermoState s;
      tbl.Eval(d*dens_cgs, t*temp_cgs, s);
      return s.chi_rho;
    }
    return 1.0;
  }

  //! \fn Real ChiRho
  //! \brief chi_rho(d, e); solves for the temperature (setup-time use)
  KOKKOS_INLINE_FUNCTION
  Real ChiRho(const Real d, const Real e) const {
    return ChiRho(d, e, Temperature(d,e));
  }

  //! \fn Real ChiT
  //! \brief chi_T = (dln p/dln T) at constant density
  KOKKOS_INLINE_FUNCTION
  Real ChiT(const Real d, const Real e, const Real t) const {
    if (tbl.active) {
      EOSThermoState s;
      tbl.Eval(d*dens_cgs, t*temp_cgs, s);
      return s.chi_t;
    }
    return 1.0;
  }

  //! \fn Real ChiT
  //! \brief chi_T(d, e); solves for the temperature (setup-time use)
  KOKKOS_INLINE_FUNCTION
  Real ChiT(const Real d, const Real e) const {
    return ChiT(d, e, Temperature(d,e));
  }

  //! \fn Real Enthalpy
  //! \brief specific enthalpy h = (e + p)/d from an already-solved temperature
  KOKKOS_INLINE_FUNCTION
  Real Enthalpy(const Real d, const Real e, const Real t) const {
    return ((e + Pressure(d,e,t))/d);
  }

  //! \fn Real Enthalpy
  //! \brief h(d, e); solves for the temperature (setup-time use)
  KOKKOS_INLINE_FUNCTION
  Real Enthalpy(const Real d, const Real e) const {
    return Enthalpy(d, e, Temperature(d,e));
  }

  //! \fn Real SoundSpeedFromP
  //! \brief adiabatic sound speed from precomputed pressure and Gamma_1
  KOKKOS_INLINE_FUNCTION
  Real SoundSpeedFromP(const Real d, const Real p, const Real g1) const {
    return sqrt(g1*p/d);
  }

  //! \fn Real FastSpeedFromP
  //! \brief fast magnetosonic speed from precomputed pressure and Gamma_1
  KOKKOS_INLINE_FUNCTION
  Real FastSpeedFromP(const Real d, const Real p, const Real g1,
                      const Real bx, const Real by, const Real bz) const {
    Real asq = g1*p;
    Real ct2 = by*by + bz*bz;
    Real qsq = bx*bx + ct2 + asq;
    Real tmp = bx*bx + ct2 - asq;
    return sqrt(0.5*(qsq + sqrt(tmp*tmp + 4.0*asq*ct2))/d);
  }

  //! \fn Real SoundSpeed
  //! \brief adiabatic sound speed sqrt(Gamma_1 p/d) at (d, e) (setup-time use; kernels
  //! use SoundSpeedFromP)
  KOKKOS_INLINE_FUNCTION
  Real SoundSpeed(const Real d, const Real e) const {
    Real t = Temperature(d,e);
    return sqrt(Gamma1(d,e,t)*Pressure(d,e,t)/d);
  }

  //! \fn Real FastSpeed
  //! \brief fast magnetosonic speed at (d, e) (setup-time use; kernels use
  //! FastSpeedFromP)
  KOKKOS_INLINE_FUNCTION
  Real FastSpeed(const Real d, const Real e,
                 const Real bx, const Real by, const Real bz) const {
    Real t = Temperature(d,e);
    Real asq = Gamma1(d,e,t)*Pressure(d,e,t);
    Real ct2 = by*by + bz*bz;
    Real qsq = bx*bx + ct2 + asq;
    Real tmp = bx*bx + ct2 - asq;
    return sqrt(0.5*(qsq + sqrt(tmp*tmp + 4.0*asq*ct2))/d);
  }

  //--------------------------------------------------------------------------------------
  // Ideal-gas interface, used by the SR/GR solvers

  // IDEAL GAS PRESSURE: converts primitive variable (either internal energy density e
  // or temperature e/d) into pressure.
  KOKKOS_INLINE_FUNCTION
  Real IdealGasPressure(const Real eint) const {
    return ((gamma-1.0)*eint);
  }

  // NON-RELATIVISTIC IDEAL GAS HYDRO: inlined sound speed function
  KOKKOS_INLINE_FUNCTION
  Real IdealHydroSoundSpeed(const Real d, const Real p) const {
    return sqrt(gamma*p/d);
  }

  // NON-RELATIVISTIC IDEAL GAS MHD: inlined fast magnetosonic speed function
  KOKKOS_INLINE_FUNCTION
  Real IdealMHDFastSpeed(const Real d, const Real p,
                         const Real bx, const Real by, const Real bz) const {
    Real asq = gamma*p;
    Real ct2 = by*by + bz*bz;
    Real qsq = bx*bx + ct2 + asq;
    Real tmp = bx*bx + ct2 - asq;
    return sqrt(0.5*(qsq + sqrt(tmp*tmp + 4.0*asq*ct2))/d);
  }

  // NON-RELATIVISTIC ISOTHERMAL MHD: inlined fast magnetosonic speed function
  KOKKOS_INLINE_FUNCTION
  Real IdealMHDFastSpeed(const Real d,
                         const Real bx, const Real by, const Real bz) const {
    Real asq = (iso_cs*iso_cs)*d;
    Real ct2 = by*by + bz*bz;
    Real qsq = bx*bx + ct2 + asq;
    Real tmp = bx*bx + ct2 - asq;
    return sqrt(0.5*(qsq + sqrt(tmp*tmp + 4.0*asq*ct2))/d);
  }

  // SPECIAL RELATIVISTIC IDEAL GAS HYDRO: inlined maximal sound wave speeds function
  // Inputs:
  //   d: density in comoving frame
  //   p: gas pressure
  //   ux: x-component of 4-velocity u^x
  //   lor: Lorentz factor \gamma
  // Outputs:
  //   l_p/m: most positive/negative wavespeed
  // Reference:
  //   Del Zanna et al, A&A 473, 11 (2007) (eq. 76)
  KOKKOS_INLINE_FUNCTION
  void IdealSRHydroSoundSpeeds(const Real d, const Real p, const Real ux, const Real lor,
                               Real& l_p, Real& l_m) const {
    Real cs2 = gamma*p / (d + gamma*p/(gamma - 1.0));  // (DZB 73)
    Real v2 = 1.0 - 1.0/(lor*lor);
    auto const p1 = (ux/lor) * (1.0 - cs2);
    auto const tmp = sqrt(cs2 * ((1.0-v2*cs2) - p1*(ux/lor))) / lor;
    auto const invden = 1.0/(1.0 - v2*cs2);

    l_p = (p1 + tmp) * invden;
    l_m = (p1 - tmp) * invden;
  }

  // SPECIAL RELATIVISTIC IDEAL GAS MHD: inlined maximal fast magnetosonic wave speeds fn
  // arguments same or SR hydro version, with the addition of b_sq = b_\mu b_\mu
  // Reference:
  //   Del Zanna et al, A&A 473, 11 (2007) (eq. 76)
  KOKKOS_INLINE_FUNCTION
  void IdealSRMHDFastSpeeds(const Real d, const Real p, const Real ux, const Real lor,
                            const Real b_sq, Real& l_p, Real& l_m) const {
    // Calculate comoving fast magnetosonic speed
    Real w = d + gamma*p/(gamma - 1.0);
    Real cs_sq = gamma*p/w;                            // (DZB 73)
    Real va_sq = b_sq / (b_sq + w);                    // (DZB 73)
    Real cms_sq = cs_sq + va_sq - cs_sq * va_sq;       // (DZB 72)

    Real v2 = 1.0 - 1.0/(lor*lor);
    auto const p1 = (ux/lor) * (1.0 - cms_sq);
    auto const tmp = sqrt(cms_sq * ((1.0-v2*cms_sq) - p1*(ux/lor))) / lor;
    auto const invden = 1.0/(1.0 - v2*cms_sq);

    l_p = (p1 + tmp) * invden;
    l_m = (p1 - tmp) * invden;
  }

  // GENERAL RELATIVISTIC IDEAL GAS HYDRO: inlined maximal sound wave speeds fn
  // Inputs:
  //  - d: density in comoving frame
  //  - p: gas pressure
  //  - u0,u1: 4-velocity components u^0, u^1
  //  - g00,g01,g11: metric components g^00, g^01, g^11
  // Outputs:
  //  - l_p/l_m: most positive/negative wavespeed
  // Notes:
  //  - Follows same general procedure as vchar() in phys.c in Harm.
  //  - Variables are named as though 1 is normal direction.
  KOKKOS_INLINE_FUNCTION
  void IdealGRHydroSoundSpeeds(const Real d, const Real p, const Real u0, const Real u1,
                               const Real g00, const Real g01, const Real g11,
                               Real& l_p, Real& l_m) const {
    // Parameters and constants
    const Real discriminant_tol = -1.0e-10;  // values between this and 0 are considered 0

    // Calculate comoving sound speed
    Real cs_sq = gamma * p / (d + gamma*p/(gamma - 1.0));

    // Set sound speeds in appropriate coordinates
    Real a = SQR(u0) - (g00 + SQR(u0)) * cs_sq;
    Real b = -2.0 * (u0*u1 - (g01 + u0*u1) * cs_sq);
    Real c = SQR(u1) - (g11 + SQR(u1)) * cs_sq;
    Real dis = SQR(b) - 4.0*a*c;
    if (dis < 0.0 && dis > discriminant_tol) {
      dis = 0.0;
    }
    // TODO(@pdmullen): fmax(dis, 0.0) prevents NaNs (see Issue #7), but this should be
    // eliminated after enforcing positivity on recon L/R densities and pressures
    Real dis_sqrt = sqrt(fmax(dis, 0.0));
    Real root_1 = (-b + dis_sqrt) / (2.0*a);
    Real root_2 = (-b - dis_sqrt) / (2.0*a);
    if (root_1 > root_2) {
      l_p = root_1;
      l_m = root_2;
    } else {
      l_p = root_2;
      l_m = root_1;
    }
  }

  // GENERAL RELATIVISTIC IDEAL GAS MHD: inlined maximal fast magnetosonic wave speeds fn
  // Inputs:
  //  - d: density in comoving frame
  //  - h: enthalpy per unit volume
  //  - p: gas pressure
  //  - u0, u1: contravariant components of 4-velocity
  //  - b_sq: b_\mu b^\mu
  //  - g00, g01, g11: contravariant components of metric (-1, 0, 1 in SR)
  // Outputs:
  //  - l_p/l_m: most positive/negative wavespeed
  // Notes:
  //  - Follows same general procedure as vchar() in phys.c in Harm.
  //  - Variables are named as though 1 is normal direction.
  KOKKOS_INLINE_FUNCTION
  void IdealGRMHDFastSpeeds(const Real d, const Real p, const Real u0, const Real u1,
                            const Real b_sq, const Real g00, const Real g01,
                            const Real g11, Real& l_p, Real& l_m) const {
    // Calculate comoving fast magnetosonic speed
    Real w = d + gamma*p/(gamma - 1.0);
    Real cs_sq = gamma * p / w;
    Real va_sq = b_sq / (b_sq + w);
    Real cms_sq = cs_sq + va_sq - cs_sq * va_sq;

    // Set fast magnetosonic speeds in appropriate coordinates
    Real a = SQR(u0) - (g00 + SQR(u0)) * cms_sq;
    Real b = -2.0 * (u0 * u1 - (g01 + u0 * u1) * cms_sq);
    Real c = SQR(u1) - (g11 + SQR(u1)) * cms_sq;
    Real a1 = b / a;
    Real a0 = c / a;
    Real s = fmax(SQR(a1) - 4.0 * a0, 0.0);
    s = sqrt(s);
    l_p = (a1 >= 0.0) ? -2.0 * a0 / (a1 + s) : (-a1 + s) / 2.0;
    l_m = (a1 >= 0.0) ? (-a1 - s) / 2.0 : -2.0 * a0 / (a1 - s);
  }
};

//----------------------------------------------------------------------------------------
//! \class EquationOfState
//! \brief Abstract base class for EOS.

class EquationOfState {
 public:
  EquationOfState(std::string block, MeshBlockPack *pp, ParameterInput *pin);
  virtual ~EquationOfState() = default;

  MeshBlockPack* pmy_pack;
  EOS_Data eos_data;

  // builds what the general EOS interface evaluates (gamma law or table); called from the
  // GeneralHydro/GeneralMHD constructors after the unit scales are set
  void BuildGeneralEOS(std::string block, ParameterInput *pin);

  // builds the composition table of an ideal-gas run for ohmic_resistivity = eos only
  // (no-op otherwise)
  void BuildElectronFractionTable(std::string block, ParameterInput *pin);

  // virtual functions to convert cons to prim in either Hydro or MHD (depending on
  // arguments), overwritten in derived eos classes
  virtual void ConsToPrim(DvceArray5D<Real> &cons, DvceArray5D<Real> &prim,
                          const bool only_testfloors,
                          const int il, const int iu, const int jl, const int ju,
                          const int kl, const int ku);
  virtual void ConsToPrim(DvceArray5D<Real> &cons, const DvceFaceFld4D<Real> &b,
                          DvceArray5D<Real> &prim, DvceArray5D<Real> &bcc,
                          const bool only_testfloors,
                          const int il, const int iu, const int jl, const int ju,
                          const int kl, const int ku);

  // virtual functions to convert prim to cons in either Hydro or MHD (depending on
  // arguments), overwritten in derived eos classes.
  virtual void PrimToCons(const DvceArray5D<Real> &prim, DvceArray5D<Real> &cons,
                          const int il, const int iu, const int jl, const int ju,
                          const int kl, const int ku);
  virtual void PrimToCons(const DvceArray5D<Real> &prim, const DvceArray5D<Real> &bcc,
                          DvceArray5D<Real> &cons, const int il, const int iu,
                          const int jl, const int ju, const int kl, const int ku);
};

//----------------------------------------------------------------------------------------
//! \class IsothermalHydro
//! \brief Derived class for isothermal EOS in nonrelativistic Hydro

class IsothermalHydro : public EquationOfState {
 public:
  // Following suppress warnings that MHD versions are not over-ridden
  using EquationOfState::ConsToPrim;
  using EquationOfState::PrimToCons;

  IsothermalHydro(MeshBlockPack *pp, ParameterInput *pin);
  void ConsToPrim(DvceArray5D<Real> &cons, DvceArray5D<Real> &prim,
                  const bool only_testfloors,
                  const int il, const int iu, const int jl, const int ju,
                  const int kl, const int ku) override;
  void PrimToCons(const DvceArray5D<Real> &prim, DvceArray5D<Real> &cons,
                  const int il, const int iu, const int jl, const int ju,
                  const int kl, const int ku) override;
};

//----------------------------------------------------------------------------------------
//! \class IdealHydro
//! \brief Derived class for ideal gas EOS in nonrelativistic hydro

class IdealHydro : public EquationOfState {
 public:
  // Following suppress warnings that MHD versions are not over-ridden
  using EquationOfState::ConsToPrim;
  using EquationOfState::PrimToCons;

  IdealHydro(MeshBlockPack *pp, ParameterInput *pin);
  void ConsToPrim(DvceArray5D<Real> &cons, DvceArray5D<Real> &prim,
                  const bool only_testfloors,
                  const int il, const int iu, const int jl, const int ju,
                  const int kl, const int ku) override;
  // inversion with the EOS_Data floor switches on (ideal_hyd_floors.cpp)
  void ConsToPrimFloors(DvceArray5D<Real> &cons, DvceArray5D<Real> &prim,
                        const bool only_testfloors,
                        const int il, const int iu, const int jl, const int ju,
                        const int kl, const int ku);
  void PrimToCons(const DvceArray5D<Real> &prim, DvceArray5D<Real> &cons,
                  const int il, const int iu, const int jl, const int ju,
                  const int kl, const int ku) override;
};

//----------------------------------------------------------------------------------------
//! \class GeneralHydro
//! \brief Derived class for a general EOS in nonrelativistic hydro; p and Gamma_1 are
//! evaluated once per cell in ConsToPrim and stored in Hydro::wder

class GeneralHydro : public EquationOfState {
 public:
  // Following suppress warnings that MHD versions are not over-ridden
  using EquationOfState::ConsToPrim;
  using EquationOfState::PrimToCons;

  GeneralHydro(MeshBlockPack *pp, ParameterInput *pin);
  void ConsToPrim(DvceArray5D<Real> &cons, DvceArray5D<Real> &prim,
                  const bool only_testfloors,
                  const int il, const int iu, const int jl, const int ju,
                  const int kl, const int ku) override;
  // inversion with the EOS_Data floor switches on (general_hyd_floors.cpp)
  void ConsToPrimFloors(DvceArray5D<Real> &cons, DvceArray5D<Real> &prim,
                        const bool only_testfloors,
                        const int il, const int iu, const int jl, const int ju,
                        const int kl, const int ku);
  // inversion of the first ConToPrim after a restart: the temperature and (p, Gamma_1)
  // are taken from the restored wtemp / wder instead of being solved for
  // (Hydro::c2p_freeze_derived; general_hyd_frozen.cpp)
  void ConsToPrimFrozen(DvceArray5D<Real> &cons, DvceArray5D<Real> &prim,
                        const int il, const int iu, const int jl, const int ju,
                        const int kl, const int ku);
  void PrimToCons(const DvceArray5D<Real> &prim, DvceArray5D<Real> &cons,
                  const int il, const int iu, const int jl, const int ju,
                  const int kl, const int ku) override;
};

//----------------------------------------------------------------------------------------
//! \class IdealSRHydro
//! \brief Derived class for ideal gas EOS in special relativistic Hydro

class IdealSRHydro : public EquationOfState {
 public:
  // Following suppress warnings that MHD versions are not over-ridden
  using EquationOfState::ConsToPrim;
  using EquationOfState::PrimToCons;

  IdealSRHydro(MeshBlockPack *pp, ParameterInput *pin);
  void ConsToPrim(DvceArray5D<Real> &cons, DvceArray5D<Real> &prim,
                  const bool only_testfloors,
                  const int il, const int iu, const int jl, const int ju,
                  const int kl, const int ku) override;
  void PrimToCons(const DvceArray5D<Real> &prim, DvceArray5D<Real> &cons,
                  const int il, const int iu, const int jl, const int ju,
                  const int kl, const int ku) override;
};

//----------------------------------------------------------------------------------------
//! \class IdealGRHydro
//! \brief Derived class for ideal gas EOS in general relativistic Hydro

class IdealGRHydro : public EquationOfState {
 public:
  // Following suppress warnings that MHD versions are not over-ridden
  using EquationOfState::ConsToPrim;
  using EquationOfState::PrimToCons;

  IdealGRHydro(MeshBlockPack *pp, ParameterInput *pin);
  void ConsToPrim(DvceArray5D<Real> &cons, DvceArray5D<Real> &prim,
                  const bool only_testfloors,
                  const int il, const int iu, const int jl, const int ju,
                  const int kl, const int ku) override;
  void PrimToCons(const DvceArray5D<Real> &prim, DvceArray5D<Real> &cons,
                  const int il, const int iu, const int jl, const int ju,
                  const int kl, const int ku) override;
};

//----------------------------------------------------------------------------------------
//! \class IsothermalMHD
//! \brief Derived class for isothermal EOS in nonrelativistic MHD

class IsothermalMHD : public EquationOfState {
 public:
  // Following suppress warnings that Hydro versions are not over-ridden
  using EquationOfState::ConsToPrim;
  using EquationOfState::PrimToCons;

  IsothermalMHD(MeshBlockPack *pp, ParameterInput *pin);
  void ConsToPrim(DvceArray5D<Real> &cons, const DvceFaceFld4D<Real> &b,
                  DvceArray5D<Real> &prim, DvceArray5D<Real> &bcc,
                  const bool only_testfloors,
                  const int il, const int iu, const int jl, const int ju,
                  const int kl, const int ku) override;
  void PrimToCons(const DvceArray5D<Real> &prim, const DvceArray5D<Real> &bcc,
                  DvceArray5D<Real> &cons, const int il, const int iu,
                  const int jl, const int ju, const int kl, const int ku) override;
};

//----------------------------------------------------------------------------------------
//! \class IdealMHD
//! \brief Derived class for ideal gas EOS in nonrelativistic MHD

class IdealMHD : public EquationOfState {
 public:
  // Following suppress warnings that Hydro versions are not over-ridden
  using EquationOfState::ConsToPrim;
  using EquationOfState::PrimToCons;

  IdealMHD(MeshBlockPack *pp, ParameterInput *pin);
  void ConsToPrim(DvceArray5D<Real> &cons, const DvceFaceFld4D<Real> &b,
                  DvceArray5D<Real> &prim, DvceArray5D<Real> &bcc,
                  const bool only_testfloors,
                  const int il, const int iu, const int jl, const int ju,
                  const int kl, const int ku) override;
  void PrimToCons(const DvceArray5D<Real> &prim, const DvceArray5D<Real> &bcc,
                  DvceArray5D<Real> &cons, const int il, const int iu,
                  const int jl, const int ju, const int kl, const int ku) override;
};

//----------------------------------------------------------------------------------------
//! \class GeneralMHD
//! \brief Derived class for a general EOS in nonrelativistic MHD; p and Gamma_1 are
//! evaluated once per cell in ConsToPrim and stored in MHD::wder

class GeneralMHD : public EquationOfState {
 public:
  // Following suppress warnings that Hydro versions are not over-ridden
  using EquationOfState::ConsToPrim;
  using EquationOfState::PrimToCons;

  GeneralMHD(MeshBlockPack *pp, ParameterInput *pin);
  void ConsToPrim(DvceArray5D<Real> &cons, const DvceFaceFld4D<Real> &b,
                  DvceArray5D<Real> &prim, DvceArray5D<Real> &bcc,
                  const bool only_testfloors,
                  const int il, const int iu, const int jl, const int ju,
                  const int kl, const int ku) override;
  // MHD counterpart of GeneralHydro::ConsToPrimFrozen (MHD::c2p_freeze_derived;
  // general_mhd_frozen.cpp)
  void ConsToPrimFrozen(DvceArray5D<Real> &cons, const DvceFaceFld4D<Real> &b,
                        DvceArray5D<Real> &prim, DvceArray5D<Real> &bcc,
                        const int il, const int iu, const int jl, const int ju,
                        const int kl, const int ku);
  void PrimToCons(const DvceArray5D<Real> &prim, const DvceArray5D<Real> &bcc,
                  DvceArray5D<Real> &cons, const int il, const int iu,
                  const int jl, const int ju, const int kl, const int ku) override;
};

//----------------------------------------------------------------------------------------
//! \class IdealSRMHD
//! \brief Derived class for ideal gas EOS in special relativistic MHD

class IdealSRMHD : public EquationOfState {
 public:
  // Following suppress warnings that hydro versions are not over-ridden
  using EquationOfState::ConsToPrim;
  using EquationOfState::PrimToCons;

  IdealSRMHD(MeshBlockPack *pp, ParameterInput *pin);
  void ConsToPrim(DvceArray5D<Real> &cons, const DvceFaceFld4D<Real> &b,
                  DvceArray5D<Real> &prim, DvceArray5D<Real> &bcc,
                  const bool only_testfloors,
                  const int il, const int iu, const int jl, const int ju,
                  const int kl, const int ku) override;
  void PrimToCons(const DvceArray5D<Real> &prim, const DvceArray5D<Real> &bcc,
                  DvceArray5D<Real> &cons, const int il, const int iu,
                  const int jl, const int ju, const int kl, const int ku) override;
};

//----------------------------------------------------------------------------------------
//! \class IdealGRMHD
//! \brief Derived class for ideal gas EOS in general relativistic MHD

class IdealGRMHD : public EquationOfState {
 public:
  // Following suppress warnings that MHD versions are not over-ridden
  using EquationOfState::ConsToPrim;
  using EquationOfState::PrimToCons;

  IdealGRMHD(MeshBlockPack *pp, ParameterInput *pin);
  void ConsToPrim(DvceArray5D<Real> &cons, const DvceFaceFld4D<Real> &b,
                  DvceArray5D<Real> &prim, DvceArray5D<Real> &bcc,
                  const bool only_testfloors,
                  const int il, const int iu, const int jl, const int ju,
                  const int kl, const int ku) override;
  void PrimToCons(const DvceArray5D<Real> &prim, const DvceArray5D<Real> &bcc,
                  DvceArray5D<Real> &cons, const int il, const int iu,
                  const int jl, const int ju, const int kl, const int ku) override;
};

#endif // EOS_EOS_HPP_
