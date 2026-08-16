//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file eos.cpp
//! \brief implements constructor and some fns for EquationOfState abstract base class

#include <float.h>

#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "parameter_input.hpp"
#include "units/units.hpp"
#include "eos/eos.hpp"

//----------------------------------------------------------------------------------------
// EquationOfState constructor

EquationOfState::EquationOfState(std::string bk, MeshBlockPack* pp, ParameterInput *pin) :
    pmy_pack(pp) {
  // Whether the input actually SET tfloor has to be asked before GetOrAddReal below,
  // which adds the parameter with its default and would make the question answer itself.
  const bool tfloor_set = pin->DoesParameterExist(bk, "tfloor");
  eos_data.dfloor = pin->GetOrAddReal(bk,"dfloor",(FLT_MIN));
  eos_data.pfloor = pin->GetOrAddReal(bk,"pfloor",(FLT_MIN));
  eos_data.tfloor = pin->GetOrAddReal(bk,"tfloor",(FLT_MIN));
  eos_data.sfloor = pin->GetOrAddReal(bk,"sfloor",(FLT_MIN));

  // <block>/tfloor_kelvin -- the SAME floor, stated in kelvin instead of code units.
  //
  // `tfloor` is a floor on the CODE temperature, and what that means changes with the
  // EOS: for an ideal gas it is p/d, whose kelvin scale carries the fixed <units>/mu, and
  // a general EOS works at mu_ref = 1 and keeps composition in the EOS, so its scale is
  // that same factor with the mu divided back out. The two differ by mu, which is how a
  // value carried across from an ideal-gas input silently becomes a floor thousands of
  // times too low or too high. Stating the floor in kelvin removes the trap and lets an
  // ideal and a general input hold the SAME number.
  //
  // Additive on purpose: `tfloor` keeps its meaning exactly, so no existing input changes
  // behaviour. Setting both is an error rather than a precedence rule.
  if (pin->DoesParameterExist(bk, "tfloor_kelvin")) {
    Real tk = pin->GetReal(bk, "tfloor_kelvin");
    if (tfloor_set) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<" << bk << "> sets both tfloor and tfloor_kelvin; they "
                << "are the same floor in different units. Set one." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // The conversion is a unit system, so there has to be one.
    if (pp->punit == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<" << bk << ">/tfloor_kelvin needs a <units> block to "
                << "convert kelvin into code temperature" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // Relativistic and isothermal EOSs are excluded: the isothermal C2P has no
    // temperature floor to set, and under SR/GR `tfloor` is not this quantity. Refuse
    // rather than convert a number that will not be used the way the name suggests.
    std::string eqn = pin->GetString(bk, "eos");
    bool general = (eqn.compare("general") == 0);
    bool ideal = (eqn.compare("ideal") == 0);
    bool rel = pin->GetOrAddBoolean("coord","special_rel",false) ||
               pin->GetOrAddBoolean("coord","general_rel",false);
    if (rel || !(general || ideal)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<" << bk << ">/tfloor_kelvin is supported for the "
                << "non-relativistic ideal and general EOSs only; this run has eos = "
                << eqn << (rel ? " and is relativistic" : "") << ". Use tfloor."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // kelvin per unit code temperature. Units::temperature_cgs() already folds in
    // <units>/mu, which is what the ideal gas wants because its code temperature is p/d;
    // the general EOS divides it back out for mu_ref = 1, exactly as EOS_Data::temp_cgs
    // is built in the GeneralHydro/GeneralMHD constructors.
    Real k_per_code = pp->punit->temperature_cgs();
    if (general) k_per_code /= pp->punit->mu();
    eos_data.tfloor = tk/k_per_code;
    if (global_variable::my_rank == 0) {
      std::cout << "<" << bk << ">/tfloor_kelvin = " << tk << " K -> tfloor = "
                << eos_data.tfloor << " in code temperature (" << k_per_code
                << " K per unit, " << (general ? "general" : "ideal") << " EOS"
                << (general ? "" : ", via <units>/mu") << ")" << std::endl;
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::BuildGeneralEOS()
//! \brief select what the general EOS interface evaluates, and build the table if asked.
//!
//! Two modes, chosen by <block>/general_eos:
//!   gamma  (default) -- the interface evaluates the same ideal gamma law as the ideal
//!                       EOS. This is not a physics option so much as a control: it makes
//!                       the general code path -- reconstruction of p and Gamma_1, the
//!                       general Riemann solvers, the general well-balanced background --
//!                       exercisable against an answer that is known exactly, which is
//!                       what the general-vs-ideal regression tests assert.
//!   table            -- the tabulated analytic EOS: H2 dissociation, H and He
//!                       ionization, and optionally radiation pressure.
//!
//! Called only from the GeneralHydro/GeneralMHD constructors, and only after dens_cgs,
//! pres_cgs and temp_cgs have been set, because the table is built in physical units.

void EquationOfState::BuildGeneralEOS(std::string block, ParameterInput *pin) {
  std::string mode = pin->GetOrAddString(block, "general_eos", "gamma");
  if (mode.compare("gamma") == 0) {
    eos_data.tbl.active = false;
  } else if (mode.compare("table") == 0) {
    // The entropy floor is a floor on the ideal-gas entropy variable p/d^gamma, which is
    // not an invariant of a tabulated EOS and cannot be given a meaning here. Refuse it
    // rather than apply a floor on a quantity that means nothing, or silently drop one
    // the user asked for. (The default FLT_MIN can never trigger, so it is allowed
    // through.) See EOS_Data::ApplyEntropyFloor().
    if (eos_data.sfloor > static_cast<Real>(FLT_MIN)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<" << block << ">/sfloor is set, but the entropy floor "
                << "is defined in terms of the ideal-gas entropy variable p/d^gamma, "
                << "which a tabulated EOS has no analogue for. Use <" << block
                << ">/pfloor or tfloor instead." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    BuildEOSTable(eos_data.tbl, pin, block, eos_data.dens_cgs, eos_data.pres_cgs,
                  eos_data.temp_cgs, eos_data.pfloor);
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<" << block << ">/general_eos = '" << mode << "' not recognized; "
              << "must be 'gamma' or 'table'" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ConsToPrim()
//! \brief No-Op versions of hydro and MHD conservative to primitive functions.
//! Required because each derived class overrides only one.

void EquationOfState::ConsToPrim(DvceArray5D<Real> &cons, DvceArray5D<Real> &prim,
                                 const bool only_testfloors,
                                 const int il, const int iu, const int jl, const int ju,
                                 const int kl, const int ku) {
}

void EquationOfState::ConsToPrim(DvceArray5D<Real> &cons, const DvceFaceFld4D<Real> &b,
                                 DvceArray5D<Real> &prim, DvceArray5D<Real> &bcc,
                                 const bool only_testfloors,
                                 const int il, const int iu, const int jl, const int ju,
                                 const int kl, const int ku) {
}

//----------------------------------------------------------------------------------------
//! \fn void PrimToCon()
//! \brief No-Op versions of hydro and MHD primitive to conservative functions.
//! Required because each derived class overrides only one.

void EquationOfState::PrimToCons(const DvceArray5D<Real> &prim, DvceArray5D<Real> &cons,
                                 const int il, const int iu, const int jl, const int ju,
                                 const int kl, const int ku) {
}
void EquationOfState::PrimToCons(const DvceArray5D<Real> &prim,
                                 const DvceArray5D<Real> &bcc, DvceArray5D<Real> &cons,
                                 const int il, const int iu, const int jl, const int ju,
                                 const int kl, const int ku) {
}
