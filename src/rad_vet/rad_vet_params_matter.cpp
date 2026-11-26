//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_params_matter.cpp
//! \brief keys of the matter side: opacities, scattering, coupling, floor caps, force
//! reference.
//!
//!  called from: rad_vet.cpp

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "units/units.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "utils/recipe.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ReadMatterKeys
//! \brief opacities, scattering, coupling to the gas, floor caps, force reference

void RadiationVet::ReadMatterKeys(ParameterInput *pin) {
  const bool nw = !global_variable::restart_run;
  // the fresh default of a scheme key, on a restart without it too (fix)
  auto fix = [&](const char *k, const char *v, const char *o) {
    return nw || recipe::Fix(pin, "rad_vet", k, v, o);
  };
  // opacities per unit mass (code units); kappa_e and kappa_f default to kappa_p
  {std::string op = pin->GetOrAddString("rad_vet","opacity",
      (nw && pin->DoesParameterExist("problem","star_opac_table")) ? "table" : "const");
  if (op.compare("const") == 0) {
    opacity_type = VET_OPAC_CONST;
  } else if (op.compare("powerlaw") == 0) {
    opacity_type = VET_OPAC_POWERLAW;
  } else if (op.compare("user") == 0) {
    opacity_type = VET_OPAC_USER;
  } else if (op.compare("table") == 0) {
    opacity_type = VET_OPAC_TABLE;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/opacity = '" << op << "' is not a valid choice "
      << "(const | powerlaw | user | table)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  }
  // table lookup units: T[K] = T_code*temp_unit_kelvin, rho_unit_cgs, kappa_unit
  otab.tunit = pin->GetOrAddReal("rad_vet","temp_unit_kelvin",1.0);
  otab.dunit = pin->GetOrAddReal("rad_vet","rho_unit_cgs",1.0);
  otab.kunit = pin->GetOrAddReal("rad_vet","kappa_unit",1.0);
  if (!(otab.tunit > 0.0) || !(otab.dunit > 0.0) || !(otab.kunit > 0.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/temp_unit_kelvin, rho_unit_cgs and kappa_unit must all "
      << "be positive" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // opac_abs_rho_max (read only when named): absorption mask of the table opacity; cells
  // with code density below it get kappa_P = kappa_E = 0 and the transport opacity
  // opac_abs_kappa_s [cm^2/g]; applied at every table lookup (VetTableOpacities)
  if (pin->DoesParameterExist("rad_vet","opac_abs_rho_max")) {
    otab.amask_rho = pin->GetReal("rad_vet","opac_abs_rho_max");
    const Real kes = pin->DoesParameterExist("rad_vet","opac_abs_kappa_s") ?
                     pin->GetReal("rad_vet","opac_abs_kappa_s") : 0.34;
    otab.amask_kes = kes*otab.kunit;
    if (otab.amask_rho < 0.0 || kes < 0.0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_vet>/opac_abs_rho_max and opac_abs_kappa_s must be >= 0"
        << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  kappa_p = pin->GetOrAddReal("rad_vet","kappa_p",0.0);
  kappa_e = pin->GetOrAddReal("rad_vet","kappa_e",kappa_p);
  kappa_f = pin->GetOrAddReal("rad_vet","kappa_f",kappa_p);
  kappa_s = pin->GetOrAddReal("rad_vet","kappa_s",0.0);
  opac_rho_ref = pin->GetOrAddReal("rad_vet","rho_ref",1.0);
  opac_t_ref = pin->GetOrAddReal("rad_vet","t_ref",1.0);
  opac_a = pin->GetOrAddReal("rad_vet","opac_a",0.0);
  opac_b = pin->GetOrAddReal("rad_vet","opac_b",0.0);
  if (!(opac_rho_ref > 0.0) || !(opac_t_ref > 0.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/rho_ref and t_ref must be positive" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (kappa_p < 0.0 || kappa_e < 0.0 || kappa_f < 0.0 || kappa_s < 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet> opacities must be non-negative" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  opac_zero = (kappa_p == 0.0 && kappa_e == 0.0 && kappa_f == 0.0 && kappa_s == 0.0 &&
               opacity_type != VET_OPAC_USER && opacity_type != VET_OPAC_TABLE);
  // vet_scatter: electron scattering.
  //  (i) opacity split: table os = kappa_e, of = kappa_a = kappa_R - kappa_e (guarded,
  //      kappa_a >= vet_kappa_floor), kappa_e = vet_kappa_es (const) or
  //      sigma_T x_e/(mu m_u) (eos); other laws: kappa_f, kappa_s. rho kappa_e goes to
  //      opac(VET_OP_S); the extinction stays the total VET_OP_T; scattering exchanges no
  //      energy.
  //  (ii) closure source VetScatterMix: the ma form,
  //      vet_scatter_j = relaxed | current.
  // Default on for a fresh run with rays, opacity =
  // table and a tabulated general EOS (kappa_e = eos); else off.
  const bool fresh = !global_variable::restart_run;
  bool gen_eos_tab = false;   // <hydro|mhd>/eos = general with general_eos = table
  for (const char *bk : {"hydro", "mhd"}) {
    if (pin->DoesBlockExist(bk) && pin->DoesParameterExist(bk,"eos") &&
        pin->GetString(bk,"eos").compare("general") == 0 &&
        pin->DoesParameterExist(bk,"general_eos") &&
        pin->GetString(bk,"general_eos").compare("table") == 0) {
      gen_eos_tab = true;
    }
  }
  const bool vet_rays = !sph_geom || VetCurvedRays(pin);
  recipe::Old(pin, "rad_vet", "vet_scatter", "false");
  vscat = pin->GetOrAddBoolean("rad_vet","vet_scatter", fresh && vet_rays &&
                               opacity_type == VET_OPAC_TABLE);
  if (vscat) {
    std::string km = "const";
    if (opacity_type == VET_OPAC_TABLE) {
      // eos where a tabulated general EOS supplies x_e
      km = pin->GetOrAddString("rad_vet","vet_scatter_kappa_e",
          (gen_eos_tab && fix("vet_scatter_kappa_e", "eos", "const")) ? "eos" : "const");
      if (km.compare("const") == 0) {
        otab.scat = 1;
      } else if (km.compare("eos") == 0) {
        otab.scat = 2;
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "<rad_vet>/vet_scatter_kappa_e = '" << km << "' not "
          << "implemented (const | eos)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      // cm^2/g -> code opacity per unit mass
      otab.kes = otab.kunit*pin->GetOrAddReal("rad_vet","vet_kappa_es",0.34);
      otab.kfl = otab.kunit*pin->GetOrAddReal("rad_vet","vet_kappa_floor",1.0e-5);
      // sigma_T/m_u [cm^2/g] (CODATA 2018)
      otab.kxe = otab.kunit*(6.6524587321e-25/1.66053906660e-24);
    }
    if (global_variable::my_rank == 0) {
      std::cout << "rad_vet: vet_scatter on, the ma form, J of the relaxed source";
      if (opacity_type == VET_OPAC_TABLE) {
        std::cout << ", table split kappa_e = " << km << " (kappa_e " << otab.kes
                  << ", floor " << otab.kfl << ", sigma_T/m_u " << otab.kxe
                  << ", code units)";
      } else {
        std::cout << ", kappa_e = kappa_s " << kappa_s;
      }
      std::cout << std::endl;
    }
  }

  // matter coupling: reads rho, v, gas energy from the fluid's conserved u0 and writes
  // u0(IEN), u0(IM1..3) back; needs <hydro> or <mhd> unless every opacity vanishes
  const bool have_mhd = pin->DoesBlockExist("mhd");
  bool have_hydro = pin->DoesBlockExist("hydro") || have_mhd;
  if (have_mhd && pin->DoesBlockExist("hydro")) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet> with both <hydro> and <mhd> (ion-neutral) is not "
      << "supported: the module couples to ONE fluid" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // <mhd> on the cubed sphere too (EmagBuild takes the C2P's orthonormal |B|^2)
  const std::string fblk = have_mhd ? "mhd" : "hydro";
  fl_on = have_hydro;
  fl_mhd = have_mhd;
  coupling = pin->GetOrAddBoolean("rad_vet","coupling",!opac_zero);
  gas_feedback = pin->GetOrAddBoolean("rad_vet","gas_feedback",true);
  // debug switch (see rad_vet.hpp)
  dbg_gas_force = pin->GetOrAddBoolean("rad_vet","dbg_gas_force",true);
  if (!dbg_gas_force && global_variable::my_rank == 0) {
    std::cout << "### WARNING: <rad_vet> DEBUG switch dbg_gas_force = false is active"
              << " -- this is not a physical configuration" << std::endl;
  }
  // removed debug switches: only their production values
  {const std::pair<const char *, bool> rm[] = {{"opac_freeze", false},
       {"dbg_gas_heat", true}, {"dbg_gas_force_trans", true}};
  for (const auto &r : rm) {
    if (pin->DoesParameterExist("rad_vet", r.first) &&
        pin->GetBoolean("rad_vet", r.first) != r.second) {
      std::cout << "### FATAL ERROR: <rad_vet>/" << r.first << " = "
                << pin->GetString("rad_vet", r.first) << " was removed (only "
                << (r.second ? "true" : "false") << ")" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }}
  if (coupling && !have_hydro) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet> matter coupling requires a <hydro> or <mhd> block"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (!opac_zero && !coupling) {
    std::cout << "### WARNING: <rad_vet> has non-zero opacities but coupling = false: "
      << "the opacity enters the transport only" << std::endl;
  }
  // an isothermal gas has no energy equation to couple to
  if (have_hydro && (coupling || !opac_zero)) {
    std::string heos = pin->GetOrAddString(fblk,"eos","ideal");
    if (heos.compare("isothermal") == 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_vet> needs a hydro EOS with an energy equation, not "
        << "'isothermal'" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // gas-only EOS
  if (have_hydro && pin->GetOrAddBoolean(fblk,"eos_radiation",false)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet> requires a GAS-ONLY EOS: <hydro>/eos_radiation must be "
      << "false (the radiation energy and pressure are evolved, not in the EOS)"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // radiation constant in code units (equilibrium E = arad*T^4, T the EOS's code
  // temperature), required with emission/absorption
  if (coupling && (kappa_p > 0.0 || kappa_e > 0.0 ||
                   opacity_type == VET_OPAC_USER ||
                   opacity_type == VET_OPAC_TABLE)) {
    arad = pin->GetReal("rad_vet","arad");
    if (!(arad > 0.0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_vet>/arad must be positive" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else {
    arad = pin->GetOrAddReal("rad_vet","arad",0.0);
  }
  // floor-gas cap (read only when named)
  if (pin->DoesParameterExist("rad_vet","floorcap_rho")) {
    fcap_rho = pin->GetReal("rad_vet","floorcap_rho");
    if (pin->DoesParameterExist("rad_vet","floorcap_v")) {
      fcap_v = pin->GetReal("rad_vet","floorcap_v");
    }
    // floorcap_lz (default false): the velocity cap scales the momentum without
    // L_z only (phi_hat: spherical polar e3, Cartesian (-y, x)/R, the cubed sphere
    // GnomonicPhiHat; the part along it is kept, also where it alone exceeds
    // floorcap_v), so that L_z is exact
    fcap_lz = pin->GetOrAddBoolean("rad_vet","floorcap_lz",fcap_lz);
    // floorcap_lz_spin (default false; needs floorcap_lz, floorcap_mesh_vol): the
    // cap scales the whole velocity, so |v| <= floorcap_v, and gives the L_z it removes
    // to the gas at or above floorcap_rho as one rigid rotation (FloorGasCap)
    fcap_spin = pin->GetOrAddBoolean("rad_vet","floorcap_lz_spin",fcap_spin);
    if (pin->DoesParameterExist("rad_vet","floorcap_report")) {
      fcap_report = pin->GetInteger("rad_vet","floorcap_report");
    }
    // floorcap_mesh_vol: counters weighted by the mesh cell volume, |dp| and L_z in the
    // frame of y (CleanVol, CleanLzArm)
    {
      fcap_mvol = pin->GetOrAddBoolean("rad_vet","floorcap_mesh_vol",
                                       fix("floorcap_mesh_vol", "true", "false"));
    }
    if (fcap_spin && !(fcap_lz && fcap_mvol)) {
      std::cout << "### FATAL ERROR: <rad_vet>/floorcap_lz_spin needs floorcap_lz and "
                << "floorcap_mesh_vol" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (fcap_rho > 0.0 && (!coupling || !(arad > 0.0))) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/floorcap_rho needs coupling and arad > 0"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (fcap_rho > 0.0 && fl_mhd) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/floorcap_rho is not implemented for <mhd> "
                << "(hydro only)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    Kokkos::realloc(fcap_cnt, 6);
    if (global_variable::my_rank == 0 && fcap_rho > 0.0) {
      std::cout << "<rad_vet> floorcap: rho < " << fcap_rho << ": |v| <= " << fcap_v
                << " (0 = no cap), T <= T_rad, clipped KE -> E_rad" << std::endl;
    }
  }
  // unit audit: print each code-unit constant next to the value <units> implies and warn
  // on a mismatch (diagnostic; the inputs are used unchanged)
  if (pmy_pack->punit != nullptr && global_variable::my_rank == 0) {
    const auto *pu = pmy_pack->punit;
    std::string eblk = pin->DoesBlockExist("mhd") ? "mhd" : "hydro";
    const bool gen = pin->DoesParameterExist(eblk, "eos") &&
                     (pin->GetString(eblk, "eos").compare("general") == 0);
    const Real tu = pu->temperature_cgs()/(gen ? pu->mu() : 1.0);
    const Real cexp = units::Units::speed_of_light_cgs/pu->velocity_cgs();
    const Real aexp = units::Units::rad_constant_cgs*SQR(SQR(tu))/pu->pressure_cgs();
    const Real dexp = pu->density_cgs();
    const Real kexp = pu->density_cgs()*pu->length_cgs();
    auto line = [](const char *key, Real val, Real ex, const char *unit) {
      const Real rel = std::abs(val/ex - 1.0);
      std::cout << "<rad_vet>/" << key << " = " << val << " [code] ; <units> implies "
                << ex << " (" << unit << "), rel. diff " << rel
                << (rel > 1.0e-5 ? "  <-- WARNING: inconsistent with <units>" : "")
                << std::endl;
    };
    line("c_light", c_light, cexp, "c in code velocity");
    if (arad > 0.0) line("arad", arad, aexp, "a_rad T_unit^4 / p_unit");
    line("temp_unit_kelvin", otab.tunit, tu, "K per code temperature");
    line("rho_unit_cgs", otab.dunit, dexp, "g/cm^3 per code density");
    line("kappa_unit", otab.kunit, kexp, "code opacity per cm^2/g");
  }

  // force_reference: wb_arad = an external well-balanced scheme delivers rho*arad_ref to
  // the gas (SetForceReference, required); clean = the clean solve's own force of its
  // first coupled solve, the hydro's static well-balanced background
  {std::string fr = pin->GetOrAddString("rad_vet","force_reference","none");
  if (fr.compare("none") == 0) {
    force_ref = VET_FREF_NONE;
  } else if (fr.compare("wb_arad") == 0) {
    force_ref = VET_FREF_WB_ARAD;
  } else if (fr.compare("clean") == 0) {
    force_ref = VET_FREF_CLEAN;
  } else if (fr.compare("wb_arad_md") == 0) {
    force_ref = VET_FREF_WB_MD;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/force_reference = '" << fr << "' is not a valid choice "
      << "(none | wb_arad | clean or wb_arad_md with scheme = clean)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  }
  // force_reference_work (wb_arad): auto = split when the pgen gives the gas the
  // reference work (fref_wsplit_ok before SetForceReference), else full
  fref_wsplit = false;
  fref_wsplit_auto = false;
  const bool fwb = (force_ref == VET_FREF_WB_ARAD) || (force_ref == VET_FREF_WB_MD);
  if (fwb || pin->DoesParameterExist("rad_vet", "force_reference_work")) {
    const std::string fw = pin->GetOrAddString("rad_vet", "force_reference_work",
                                               "auto");
    if (fw.compare("split") == 0) {
      fref_wsplit = fwb;
    } else if (fw.compare("auto") == 0) {
      fref_wsplit = fwb;
      fref_wsplit_auto = fref_wsplit;
    } else if (fw.compare("full") != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_vet>/force_reference_work = '" << fw << "' is not a valid "
        << "choice (auto | full | split)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // the hydro integrator's stage weights (driver.cpp), for fref_wacc
  {
    const std::string hi = pin->DoesParameterExist("time", "integrator") ?
                           pin->GetString("time", "integrator") : "rk2";
    fref_hnst = 0;
    if (hi.compare("rk1") == 0) {
      fref_hnst = 1;
      fref_hgam0[0] = 0.0;  fref_hbeta[0] = 1.0;
    } else if (hi.compare("rk2") == 0) {
      fref_hnst = 2;
      fref_hgam0[0] = 0.0;  fref_hbeta[0] = 1.0;
      fref_hgam0[1] = 0.5;  fref_hbeta[1] = 0.5;
      if (pin->DoesParameterExist("rad_vet", "time2_tableau") &&
          pin->GetString("rad_vet", "time2_tableau").compare("trbdf2") == 0) {
        const Real c2 = 2.0 - std::sqrt(2.0);
        fref_hbeta[0] = c2;
        fref_hgam0[1] = 0.25;  fref_hbeta[1] = 0.5/c2;
      }
    } else if (hi.compare("rk3") == 0) {
      fref_hnst = 3;
      fref_hgam0[0] = 0.0;      fref_hbeta[0] = 1.0;
      fref_hgam0[1] = 0.25;     fref_hbeta[1] = 0.25;
      fref_hgam0[2] = 2.0/3.0;  fref_hbeta[2] = 2.0/3.0;
    }
  }
  if (force_ref != VET_FREF_NONE && !(coupling && gas_feedback)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/force_reference needs coupling = true and "
      << "gas_feedback = true: it modifies the momentum handed to the gas" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // RSLA start-up check
  rsla_warn = pin->GetOrAddReal("rad_vet","rsla_warn",0.1);
  rsla_vmax = pin->GetOrAddReal("rad_vet","rsla_vmax",-1.0);
  rsla_force = pin->GetOrAddBoolean("rad_vet","rsla_force",false);
  rsla_done = false;
}

}  // namespace radvet
