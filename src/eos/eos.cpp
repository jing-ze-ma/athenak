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
#include "coordinates/coordinates.hpp"
#include "parameter_input.hpp"
#include "units/units.hpp"
#include "eos/eos.hpp"
#include "utils/recipe.hpp"

//----------------------------------------------------------------------------------------
// EquationOfState constructor

EquationOfState::EquationOfState(std::string bk, MeshBlockPack* pp, ParameterInput *pin) :
    pmy_pack(pp) {
  // whether the input set tfloor, asked before GetOrAddReal adds it; a default recorded
  // in a restart file does not count as set
  const bool tfloor_set = pin->DoesParameterExist(bk, "tfloor") &&
                          !pin->IsParameterDefaulted(bk, "tfloor");
  eos_data.dfloor = pin->GetOrAddReal(bk,"dfloor",(FLT_MIN));
  eos_data.dfloor_keep_velocity =
      pin->GetOrAddBoolean(bk,"dfloor_keep_velocity",false);
  // dfloor_keep_velocity is implemented in the hydro inversions only
  if (eos_data.dfloor_keep_velocity && bk.compare("hydro") != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<" << bk << ">/dfloor_keep_velocity is implemented only "
              << "for hydrodynamics" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // velocity ceiling (EOS_Data::vceil): Newtonian hydro and MHD only
  eos_data.vceil = pin->GetOrAddReal(bk,"vceil",0.0);
  if (eos_data.vceil > 0.0 &&
      (bk.compare("hydro") != 0 && bk.compare("mhd") != 0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<" << bk << ">/vceil is implemented only "
              << "for hydrodynamics and MHD" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (eos_data.vceil > 0.0 && bk.compare("mhd") == 0) {
    if (pp->pcoord->is_special_relativistic || pp->pcoord->is_general_relativistic) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<mhd>/vceil is implemented only for NEWTONIAN MHD; "
                << "the relativistic inversions use gamma_max instead" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // Alfven-speed density floor (EOS_Data::va2max_inv): Newtonian MHD only
  if (pin->DoesParameterExist(bk,"valfven_max_cgs")) {
    if (bk.compare("mhd") != 0 || pp->pcoord->is_special_relativistic ||
        pp->pcoord->is_general_relativistic) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<" << bk << ">/valfven_max_cgs is implemented only for "
                << "Newtonian <mhd>" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    const Real va_cgs = pin->GetReal(bk,"valfven_max_cgs");
    const Real vunit = (pp->punit != nullptr) ? pp->punit->velocity_cgs() : 1.0;
    if (va_cgs > 0.0) {
      const Real va = va_cgs/vunit;
      eos_data.va2max_inv = 1.0/(va*va);
    }
    if (global_variable::my_rank == 0) {
      std::cout << "<mhd>/valfven_max_cgs = " << va_cgs << " [cm/s] = " << va_cgs/vunit
                << " [code]: density floor rho >= B^2/v_A,max^2" << std::endl;
    }
  }
  // no-infall companion of the Alfven floor (EOS_Data::va_noinfall_fac); needs
  // valfven_max_cgs and a radial x1
  if (pin->DoesParameterExist(bk,"valfven_noinfall_fac")) {
    const Real fac = pin->GetReal(bk,"valfven_noinfall_fac");
    if (fac > 0.0) {
      if (eos_data.va2max_inv <= 0.0 || fac < 1.0 ||
          !(pp->pmesh->use_cubed_sphere || pp->pmesh->use_spherical_polar)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<" << bk << ">/valfven_noinfall_fac needs "
                  << "valfven_max_cgs > 0, a value >= 1 (or 0 = off) and a radial x1 "
                  << "(cubed sphere or spherical polar)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      eos_data.va_noinfall_fac = fac;
      if (global_variable::my_rank == 0) {
        std::cout << "<mhd>/valfven_noinfall_fac = " << fac << ": inward radial momentum "
                  << "removed where rho <= " << fac << " B^2/v_A,max^2" << std::endl;
      }
    }
  }
  // HLLD ** state tolerance (EOS_Data::hlld_bx_zero_tol), MHD only; recorded with its
  // default in <mhd> only, so a <hydro> restart is not refused
  const bool hlld_tol_set = pin->DoesParameterExist(bk,"hlld_bx_zero_tol");
  if (bk.compare("mhd") == 0) {
    eos_data.hlld_bx_zero_tol = pin->GetOrAddReal(bk,"hlld_bx_zero_tol",1.0e-4);
  }
  if (hlld_tol_set && bk.compare("mhd") != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<" << bk << ">/hlld_bx_zero_tol is implemented only "
              << "for MHD" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // radial low-Mach floor (EOS_Data::x1_phi_min)
  {
    const bool is_hyd = (bk.compare("hydro") == 0);
    const std::string key = is_hyd ? "lhllc_x1_phi_min" : "lhlld_x1_phi_min";
    const std::string bad = is_hyd ? "lhlld_x1_phi_min" : "lhllc_x1_phi_min";
    if (pin->DoesParameterExist(bk, bad)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<" << bk << ">/" << bad << " is not valid in this "
                << "block; use <" << bk << ">/" << key << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (pin->DoesParameterExist(bk, key)) {
      eos_data.x1_phi_min = pin->GetReal(bk, key);
    }
  }
  // non-negative EMF resistivity of LHLLD (EOS_Data::lhlld_emf_eta_pos), default true in
  // <mhd>, MHD only
  if (bk.compare("mhd") == 0) eos_data.lhlld_emf_eta_pos = true;
  if (pin->DoesParameterExist(bk, "lhlld_emf_eta_pos")) {
    if (bk.compare("mhd") != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<" << bk << ">/lhlld_emf_eta_pos is implemented only "
                << "for MHD (LHLLD)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    eos_data.lhlld_emf_eta_pos = pin->GetBoolean(bk, "lhlld_emf_eta_pos");
  }
  // odd-even sensor of LHLLC/LHLLD (EOS_Data::cb_sensor), default 1; 0 = off
  eos_data.cb_sensor = 1.0;
  {
    const bool is_hyd = (bk.compare("hydro") == 0);
    const std::string key = is_hyd ? "lhllc_cb_sensor" : "lhlld_cb_sensor";
    const std::string bad = is_hyd ? "lhlld_cb_sensor" : "lhllc_cb_sensor";
    if (pin->DoesParameterExist(bk, bad)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<" << bk << ">/" << bad << " is not valid in this "
                << "block; use <" << bk << ">/" << key << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (pin->DoesParameterExist(bk, key)) {
      eos_data.cb_sensor = pin->GetReal(bk, key);
    }
    // continuous sensor fraction (EOS_Data::cb_cont)
    const std::string kc = is_hyd ? "lhllc_cb_cont" : "lhlld_cb_cont";
    const std::string bc = is_hyd ? "lhlld_cb_cont" : "lhllc_cb_cont";
    if (pin->DoesParameterExist(bk, bc)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<" << bk << ">/" << bc << " is not valid in this "
                << "block; use <" << bk << ">/" << kc << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // fresh-start default (recipe::NewDefaults): cb_cont = true with rsolver =
    // lhllc/lhlld, recorded for the restarts; a restart input without the key keeps false
    {
      const std::string rs = is_hyd ? "lhllc" : "lhlld";
      if (!pin->DoesParameterExist(bk, kc) && pin->DoesParameterExist(bk, "rsolver") &&
          pin->GetString(bk, "rsolver") == rs) {
        if (recipe::NewDefaults(pin)) {
          pin->GetOrAddBoolean(bk, kc, true);
          recipe::Log("<" + bk + ">/" + kc + " = true (default with rsolver = " + rs +
                      ")");
        } else {
          recipe::Log("<" + bk + ">/" + kc + " = false (" + std::string(
              recipe::GetMode(pin) == recipe::Mode::none ? "recipe none" :
              "restart input has no key") + ": the old default)");
        }
      }
    }
    if (pin->DoesParameterExist(bk, kc)) {
      eos_data.cb_cont = pin->GetBoolean(bk, kc);
    }
    // well-balanced face range bound (EOS_Data::wb_plm_range)
    if (pin->DoesParameterExist(bk, "wb_plm_range")) {
      eos_data.wb_plm_range = pin->GetBoolean(bk, "wb_plm_range");
      if (eos_data.wb_plm_range &&
          !(pp->pmesh->use_cubed_sphere || pp->pmesh->use_spherical_polar)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<" << bk << ">/wb_plm_range is implemented in the "
                  << "non-uniform x1 well-balanced PLM only" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    // well-balanced face blend (EOS_Data::wb_plm_smooth)
    if (pin->DoesParameterExist(bk, "wb_plm_smooth")) {
      eos_data.wb_plm_smooth = pin->GetBoolean(bk, "wb_plm_smooth");
      if (eos_data.wb_plm_smooth &&
          !(pp->pmesh->use_cubed_sphere || pp->pmesh->use_spherical_polar)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<" << bk << ">/wb_plm_smooth is implemented in the "
                  << "non-uniform x1 well-balanced PLM only (cubed sphere, spherical "
                  << "polar)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    // fresh-start default: lhlld_chi_normal = true with rsolver = lhlld (as cb_cont
    // above); LHLLC always uses the face-normal chi
    if (!is_hyd && !pin->DoesParameterExist(bk, "lhlld_chi_normal") &&
        pin->DoesParameterExist(bk, "rsolver") &&
        pin->GetString(bk, "rsolver") == "lhlld") {
      if (recipe::NewDefaults(pin)) {
        pin->GetOrAddBoolean(bk, "lhlld_chi_normal", true);
        recipe::Log("<mhd>/lhlld_chi_normal = true (default with rsolver = lhlld)");
      } else {
        recipe::Log("<mhd>/lhlld_chi_normal = false (" + std::string(
            recipe::GetMode(pin) == recipe::Mode::none ? "recipe none" :
            "restart input has no key") + ": the old default)");
      }
    }
    if (pin->DoesParameterExist(bk, "lhlld_chi_normal")) {
      if (is_hyd) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<hydro>/lhlld_chi_normal is MHD (LHLLD) only"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      eos_data.chi_normal = pin->GetBoolean(bk, "lhlld_chi_normal");
    }
    // fresh-start default: lhlld_alfmach = true with rsolver = lhlld (as cb_cont above)
    if (!is_hyd && !pin->DoesParameterExist(bk, "lhlld_alfmach") &&
        pin->DoesParameterExist(bk, "rsolver") &&
        pin->GetString(bk, "rsolver") == "lhlld") {
      if (recipe::NewDefaults(pin)) {
        pin->GetOrAddBoolean(bk, "lhlld_alfmach", true);
        recipe::Log("<mhd>/lhlld_alfmach = true (default with rsolver = lhlld)");
      } else {
        recipe::Log("<mhd>/lhlld_alfmach = false (" + std::string(
            recipe::GetMode(pin) == recipe::Mode::none ? "recipe none" :
            "restart input has no key") + ": the old default)");
      }
    }
    // Alfven-Mach scaled floor of LHLLD (EOS_Data::alfmach)
    if (pin->DoesParameterExist(bk, "lhlld_alfmach")) {
      if (is_hyd) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<hydro>/lhlld_alfmach is MHD (LHLLD) only"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (pin->GetOrAddString(bk, "rsolver", "none").compare("lhlld") != 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<mhd>/lhlld_alfmach needs <mhd>/rsolver = lhlld"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      eos_data.alfmach = pin->GetBoolean(bk, "lhlld_alfmach");
    }
  }
  eos_data.pfloor = pin->GetOrAddReal(bk,"pfloor",(FLT_MIN));
  eos_data.tfloor = pin->GetOrAddReal(bk,"tfloor",(FLT_MIN));

  // energy floor paid from the kinetic energy (EOS_Data::efloor_from_ekin)
  eos_data.efloor_from_ekin = pin->GetOrAddBoolean(bk,"efloor_from_ekin",false);
  // floored states rebuilt from their temperature (EOS_Data::efloor_as_tfloor)
  eos_data.efloor_as_tfloor = pin->GetOrAddBoolean(bk,"efloor_as_tfloor",false);
  // both are implemented in the hydro general-EOS inversion only (general_c2p_hyd.hpp)
  if ((eos_data.efloor_from_ekin || eos_data.efloor_as_tfloor) &&
      bk.compare("hydro") != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<" << bk << ">/"
              << (eos_data.efloor_from_ekin ? "efloor_from_ekin" : "efloor_as_tfloor")
              << " is implemented only for hydrodynamics (general_c2p_hyd.hpp)"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  eos_data.sfloor = pin->GetOrAddReal(bk,"sfloor",(FLT_MIN));
  // cubed sphere: conserved floor writes deferred to GnomonicEquiangleRaiseVel{,MHD}
  eos_data.defer_cons_floors = pp->pmesh->use_cubed_sphere;

  // original floor kernels when no floor switch is on (EOS_Data::floors_legacy)
  eos_data.floors_legacy = !(eos_data.dfloor_keep_velocity ||
                             (eos_data.vceil > 0.0) ||
                             eos_data.efloor_from_ekin ||
                             eos_data.efloor_as_tfloor);

  // <block>/tfloor_kelvin: tfloor in kelvin.  The code temperature of an ideal gas
  // carries <units>/mu while the general EOS works at mu_ref = 1, so a kelvin value means
  // the same in both; setting tfloor as well is an error
  if (pin->DoesParameterExist(bk, "tfloor_kelvin")) {
    Real tk = pin->GetReal(bk, "tfloor_kelvin");
    if (tfloor_set) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<" << bk << "> sets both tfloor and tfloor_kelvin; they "
                << "are the same floor in different units. Set one." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (pp->punit == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<" << bk << ">/tfloor_kelvin needs a <units> block to "
                << "convert kelvin into code temperature" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // non-relativistic ideal and general EOS only
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
    // kelvin per unit code temperature; the general EOS divides <units>/mu back out, as
    // for EOS_Data::temp_cgs
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
//! \brief selects what the general EOS interface evaluates, <block>/general_eos:
//!   gamma (default): the ideal gamma law through the general code path (reference for
//!                    the general-vs-ideal tests)
//!   table:           the tabulated analytic EOS (H2 dissociation, H/He ionization,
//!                    optional radiation pressure)
//! Called from the GeneralHydro/GeneralMHD constructors after the unit scales are set.

void EquationOfState::BuildGeneralEOS(std::string block, ParameterInput *pin) {
  std::string mode = pin->GetOrAddString(block, "general_eos", "gamma");
  if (mode.compare("gamma") == 0) {
    eos_data.tbl.active = false;
  } else if (mode.compare("table") == 0) {
    // the entropy floor p/d^gamma has no meaning for a table; refuse a set sfloor
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
    // efloor_as_tfloor: ClampLogT holds every inversion at the lowest tabulated row, so a
    // tfloor below the table can never fire; raise it to the table edge
    if (eos_data.efloor_as_tfloor) {
      const Real tmin = eos_data.TableTempMin();
      if (tmin > 0.0 && eos_data.tfloor < tmin) {
        if (global_variable::my_rank == 0) {
          std::cout << "General EOS: <" << block << ">/efloor_as_tfloor raised tfloor "
                    << eos_data.tfloor*eos_data.temp_cgs << " K -> "
                    << tmin*eos_data.temp_cgs << " K, the table's lowest tabulated "
                    << "temperature; a floor below it is unreachable (the inversion is "
                    << "clamped to that row first)." << std::endl;
        }
        eos_data.tfloor = tmin;
      }
    }
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<" << block << ">/general_eos = '" << mode << "' not recognized; "
              << "must be 'gamma' or 'table'" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::BuildElectronFractionTable()
//! \brief builds the composition table of an ideal-gas run for the electron fraction of
//! ohmic_resistivity = eos (Saha over H and the metal donors) without making the
//! thermodynamics tabulated.  The temperature remains the ideal-gas one (fixed mu), so
//! x_e is only as good as that temperature.

void EquationOfState::BuildElectronFractionTable(std::string block, ParameterInput *pin) {
  // only for ohmic_resistivity = eos; probed without GetOrAdd, which would add the key
  if (!pin->DoesParameterExist(block, "ohmic_resistivity")) return;
  if (pin->GetString(block, "ohmic_resistivity").compare("eos") != 0) return;

  if (pmy_pack->punit == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<" << block << ">/ohmic_resistivity = eos with an ideal gas builds a "
              << "composition table, which needs a <units> block to know the physical "
              << "density and temperature scale" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  Real dens_cgs = pmy_pack->punit->density_cgs();
  Real pres_cgs = pmy_pack->punit->pressure_cgs();
  Real temp_cgs = pmy_pack->punit->temperature_cgs()/pmy_pack->punit->mu();

  BuildEOSTable(eos_data.tbl, pin, block, dens_cgs, pres_cgs, temp_cgs, eos_data.pfloor);

  // thermodynamics stay ideal: only ElectronFractionKelvin reads the table
  // (xe_from_table)
  eos_data.tbl.active = false;
  eos_data.xe_from_table = true;

  // the resistivity's ideal-gas T = p/(Rgas rho) is in kelvin only for a cgs run
  if (global_variable::my_rank == 0) {
    std::cout << "<" << block << ">: ideal gas + tabulated electron fraction for "
              << "ohmic_resistivity = eos; the thermodynamics stay ideal. Table indexed "
              << "at " << dens_cgs << " g/cm^3 per code density, and the resistivity's "
              << "own T = p/(Rgas rho) must already be in kelvin." << std::endl;
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
