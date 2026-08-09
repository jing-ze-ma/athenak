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

//----------------------------------------------------------------------------------------
//! \enum EOSType
//! \brief Selects which equation of state the general EOS interface in EOS_Data
//! evaluates. EOSType::ideal reproduces the ideal-gas expressions identically, so runs
//! that select it are numerically unchanged by the general interface.

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
  Real gamma_max;    // ceiling on Lorentz factor in SR/GR
  EOSType eos_type = EOSType::ideal;  // EOS evaluated by the general interface below

  // Unit conversion factors, code units -> cgs. An ideal gas is scale free and ignores
  // these, but a general EOS (ionization, dissociation, radiation pressure) depends on
  // ABSOLUTE density and temperature, so it must know the physical scale. All three are
  // 1.0 when the run itself is in cgs units, which is the recommended setup. They are
  // copied out of the units::Units class at construction because that class is host-only
  // and cannot be captured in a device lambda.
  Real dens_cgs = 1.0;   // g/cm^3       per unit code density
  Real pres_cgs = 1.0;   // erg/cm^3     per unit code pressure or energy density
  Real temp_cgs = 1.0;   // K            per unit code temperature

  //! \fn bool IsGeneral
  //! \brief true when the EOS is expensive enough that p and Gamma_1 must be precomputed
  //! into the derived-variable arrays and reconstructed, rather than evaluated inside the
  //! Riemann solvers. False for an ideal gas, which keeps the original code path.
  KOKKOS_INLINE_FUNCTION
  bool IsGeneral() const {
    return (eos_type == EOSType::general);
  }

  //--------------------------------------------------------------------------------------
  // GENERAL EOS INTERFACE for non-relativistic Hydro and MHD.
  //
  // These are the only EOS entry points non-relativistic kernels should call. They are
  // written in terms of the primitive pair (d, e) = (density, internal energy density),
  // which is what is actually stored in prim(IDN,...) and prim(IEN,...). Pressure is a
  // derived quantity: for a general EOS it depends on BOTH d and e, which is why these
  // take a density argument where the ideal-gas versions below do not.
  //
  // Every function reduces exactly to the ideal-gas expression when
  // eos_type == EOSType::ideal, so existing ideal-gas runs are bitwise unchanged.

  //! \fn Real Pressure
  //! \brief gas pressure p(d,e)
  KOKKOS_INLINE_FUNCTION
  Real Pressure(const Real d, const Real e) const {
    // TODO(stage3): add EOSType::general branch (analytic partial-ionization EOS)
    return ((gamma-1.0)*e);
  }

  //! \fn Real Gamma1
  //! \brief first adiabatic exponent Gamma_1 = (dln p/dln d) at constant entropy. This is
  //! the exponent that sets the sound speed, and is what replaces the constant `gamma` in
  //! the wave-structure algebra of the HLLC/HLLD/Roe solvers. For an ideal gas it is the
  //! constant ratio of specific heats; for a general EOS it varies with (d,e) and dips
  //! well below 5/3 inside partial-ionization zones.
  KOKKOS_INLINE_FUNCTION
  Real Gamma1(const Real d, const Real e) const {
    // TODO(stage3): add EOSType::general branch
    return gamma;
  }

  //! \fn Real EnergyFromPressure
  //! \brief internal energy density e(d,p); inverse of Pressure(). Needed for pressure
  //! floors, which are specified in terms of p but applied to e.
  KOKKOS_INLINE_FUNCTION
  Real EnergyFromPressure(const Real d, const Real p) const {
    // TODO(stage3): add EOSType::general branch
    return (p/(gamma-1.0));
  }

  //! \fn Real Temperature
  //! \brief temperature T(d,e) in code units (p/d for an ideal gas)
  KOKKOS_INLINE_FUNCTION
  Real Temperature(const Real d, const Real e) const {
    // TODO(stage3): add EOSType::general branch
    return ((gamma-1.0)*e/d);
  }

  //! \fn Real EnergyFromTemperature
  //! \brief internal energy density e(d,T); inverse of Temperature(). Needed for the
  //! temperature floor, which is specified in T but must be applied to e.
  KOKKOS_INLINE_FUNCTION
  Real EnergyFromTemperature(const Real d, const Real t) const {
    // TODO(stage3): add EOSType::general branch
    return (d*t/(gamma-1.0));
  }

  //! \fn Real MeanMolecularWeight
  //! \brief mean molecular weight mu(d,e), i.e. mass per particle in units of m_u. This
  //! is what converts a mass density into a number density, which modules that do
  //! per-particle microphysics (resistivity, cooling) need. Under the mu_ref = 1
  //! convention the ideal-gas interface carries no composition at all, so it returns 1;
  //! a general EOS returns the value its own composition implies, dropping from ~2.3 in
  //! molecular H2/He to ~0.6 once hydrogen is ionized.
  KOKKOS_INLINE_FUNCTION
  Real MeanMolecularWeight(const Real d, const Real e) const {
    // TODO(stage3): add EOSType::general branch
    return 1.0;
  }

  //! \fn Real DensityFromPressureTemperature
  //! \brief density d(p,T); the (p,T) -> d inversion. Unlike everything else here this is
  //! NOT a function of the primitive pair, and no kernel needs it: it exists for problem
  //! generators that build a hydrostatic background by integrating dln p/dz in (p,T)
  //! space and then have to close it with a density. For an ideal gas p = d T under the
  //! mu_ref = 1 convention, so the inverse is exact; a general EOS has to root find on d,
  //! which is why this is used at setup only and never inside a time step.
  KOKKOS_INLINE_FUNCTION
  Real DensityFromPressureTemperature(const Real p, const Real t) const {
    // TODO(stage3): add EOSType::general branch (root find on d at fixed p and T)
    return (p/t);
  }

  //! \fn Real SpecificHeatCv
  //! \brief specific heat at constant volume c_v(d,e) = (de/dT)/d. Used by thermal
  //! conduction, which currently hardwires the ideal-gas value 1/(gamma-1).
  KOKKOS_INLINE_FUNCTION
  Real SpecificHeatCv(const Real d, const Real e) const {
    // TODO(stage3): add EOSType::general branch
    return (1.0/(gamma-1.0));
  }

  //! \fn Real ChiRho
  //! \brief chi_rho = (dln p/dln d) at constant TEMPERATURE. Note this is NOT Gamma_1,
  //! which is the same derivative at constant entropy; the two coincide only for an
  //! isothermal gas. Needed by the well-balanced scheme, which integrates a hydrostatic
  //! background at fixed T, and by Gamma_1 itself through
  //! Gamma_1 = chi_rho + p chi_T^2/(d T c_v).
  KOKKOS_INLINE_FUNCTION
  Real ChiRho(const Real d, const Real e) const {
    // TODO(stage3): add EOSType::general branch
    return 1.0;
  }

  //! \fn Real ChiT
  //! \brief chi_T = (dln p/dln T) at constant density. Unity for an ideal gas; it departs
  //! from unity wherever the number of particles depends on temperature (dissociation,
  //! ionization) or radiation pressure contributes.
  KOKKOS_INLINE_FUNCTION
  Real ChiT(const Real d, const Real e) const {
    // TODO(stage3): add EOSType::general branch
    return 1.0;
  }

  //! \fn Real Enthalpy
  //! \brief specific enthalpy h = (e + p)/d. Along a hydrostatic isentrope h + Phi is
  //! constant, which is what makes the isentropic well-balanced background integrable.
  KOKKOS_INLINE_FUNCTION
  Real Enthalpy(const Real d, const Real e) const {
    return ((e + Pressure(d,e))/d);
  }

  //! \fn Real SoundSpeedFromP
  //! \brief adiabatic sound speed from an ALREADY EVALUATED pressure and Gamma_1. Used
  //! wherever those were precomputed in ConsToPrim (timestep, FOFC, Riemann solvers), so
  //! that the expensive EOS evaluation is not repeated.
  KOKKOS_INLINE_FUNCTION
  Real SoundSpeedFromP(const Real d, const Real p, const Real g1) const {
    return sqrt(g1*p/d);
  }

  //! \fn Real FastSpeedFromP
  //! \brief fast magnetosonic speed from an ALREADY EVALUATED pressure and Gamma_1
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
  //! \brief adiabatic sound speed c_s(d,e) = sqrt(Gamma_1 p/d)
  KOKKOS_INLINE_FUNCTION
  Real SoundSpeed(const Real d, const Real e) const {
    return sqrt(Gamma1(d,e)*Pressure(d,e)/d);
  }

  //! \fn Real FastSpeed
  //! \brief fast magnetosonic speed for a general EOS, from (d,e) and the field
  KOKKOS_INLINE_FUNCTION
  Real FastSpeed(const Real d, const Real e,
                 const Real bx, const Real by, const Real bz) const {
    Real asq = Gamma1(d,e)*Pressure(d,e);
    Real ct2 = by*by + bz*bz;
    Real qsq = bx*bx + ct2 + asq;
    Real tmp = bx*bx + ct2 - asq;
    return sqrt(0.5*(qsq + sqrt(tmp*tmp + 4.0*asq*ct2))/d);
  }

  //--------------------------------------------------------------------------------------
  // IDEAL GAS interface. Retained for the SR/GR solvers, which are out of scope for the
  // general EOS; non-relativistic kernels should use the general interface above.

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
  void PrimToCons(const DvceArray5D<Real> &prim, DvceArray5D<Real> &cons,
                  const int il, const int iu, const int jl, const int ju,
                  const int kl, const int ku) override;
};

//----------------------------------------------------------------------------------------
//! \class GeneralHydro
//! \brief Derived class for a general EOS in nonrelativistic hydro. Evaluates pressure
//! and Gamma_1 once per cell in ConsToPrim and stores them in Hydro::wder, from where
//! they are reconstructed to interfaces.

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
//! \brief Derived class for a general EOS in nonrelativistic MHD. Evaluates pressure and
//! Gamma_1 once per cell in ConsToPrim and stores them in MHD::wder, from where they are
//! reconstructed to interfaces.

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
