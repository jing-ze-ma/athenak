//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd.cpp
//! \brief implementation of MHD class constructor and assorted functions

#include <cstdlib>
#include <iostream>
#include <string>
#include <algorithm>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "diffusion/viscosity.hpp"
#include "diffusion/resistivity.hpp"
#include "diffusion/conduction.hpp"
#include "srcterms/srcterms.hpp"
#include "shearing_box/shearing_box.hpp"
#include "shearing_box/orbital_advection.hpp"
#include "bvals/bvals.hpp"
#include "mhd/mhd.hpp"
#include "hydro/hydro_am.hpp"
#include "utils/recipe.hpp"
#include "units/units.hpp"

namespace {
//----------------------------------------------------------------------------------------
//! tracer-1101 defaults (user 11-01).  Each applies only when its key is NOT named and
//! recipe::NewDefaults(pin) holds, i.e. on a FRESH start (or <job>/recipe = prod1019);
//! a restart whose embedded input lacks the key keeps the old behaviour (off).  The
//! resolved value is written into the ParameterInput (GetOrAdd), so it is embedded in
//! every restart the run writes and the restart continues with it.  An explicit key in
//! the deck or on the command line always wins.  Every decision goes to the recipe log.
//!   <mhd>/valfven_noinfall_fac = 4    when valfven_max_cgs > 0 and x1 is radial
//!   <mhd>/floor_tracer = true         when valfven_max_cgs > 0
//!   <mhd>/xe_tracer = true            (DefaultXeTracer, after the EOS exists) when
//!                                     ohmic_resistivity = eos or ambipolar = eos and the
//!                                     EOS has the electron-fraction table (else false)
//!   <mhd>/ambipolar = eos             (DefaultAmbipolar, after DefaultXeTracer) when
//!                                     ohmic_resistivity = eos and xe_tracer is on

void LogDefault(ParameterInput *pin, const std::string &key, const std::string &on) {
  if (recipe::NewDefaults(pin)) {
    recipe::Log(key + " = " + on);
  } else {
    recipe::Log(key + " = off (" + std::string(
        recipe::GetMode(pin) == recipe::Mode::none ? "recipe none" :
        "restart input has no key") + ": the old default)");
  }
}

void DefaultAlfvenFloorKeys(MeshBlockPack *ppack, ParameterInput *pin) {
  if (!pin->DoesParameterExist("mhd","valfven_max_cgs")) return;
  if (!(pin->GetReal("mhd","valfven_max_cgs") > 0.0)) return;
  const bool radial = ppack->pmesh->use_cubed_sphere || ppack->pmesh->use_spherical_polar;
  if (radial && !pin->DoesParameterExist("mhd","valfven_noinfall_fac")) {
    LogDefault(pin, "<mhd>/valfven_noinfall_fac", "4 (default with valfven_max_cgs > 0)");
    if (recipe::NewDefaults(pin)) pin->GetOrAddReal("mhd","valfven_noinfall_fac",4.0);
  }
  if (!pin->DoesParameterExist("mhd","floor_tracer")) {
    LogDefault(pin, "<mhd>/floor_tracer", "true (default with valfven_max_cgs > 0)");
    if (recipe::NewDefaults(pin)) pin->GetOrAddBoolean("mhd","floor_tracer",true);
  }
}

void DefaultXeTracer(ParameterInput *pin, const EOS_Data &eos) {
  if (pin->DoesParameterExist("mhd","xe_tracer")) return;
  const bool ohm_eos = pin->DoesParameterExist("mhd","ohmic_resistivity") &&
                       pin->GetString("mhd","ohmic_resistivity").compare("eos") == 0;
  const bool ad_eos = pin->DoesParameterExist("mhd","ambipolar") &&
                      pin->GetString("mhd","ambipolar").compare("eos") == 0;
  if (!(ohm_eos || ad_eos)) return;
  if (!eos.HasElectronFraction()) {
    recipe::Log("<mhd>/xe_tracer = off (no EOS electron-fraction table)");
    return;
  }
  LogDefault(pin, "<mhd>/xe_tracer", "true (default with ohmic_resistivity or ambipolar "
             "= eos)");
  if (recipe::NewDefaults(pin)) pin->GetOrAddBoolean("mhd","xe_tracer",true);
}

//! adeos-1101 (user 11-01): <mhd>/ambipolar = eos on a fresh start (same rules as
//! above) when ohmic_resistivity = eos AND the x_e tracer is on (after DefaultXeTracer).
//! AD with the local Saha x_e costs +31 % RKG stages: not defaulted without the tracer.
void DefaultAmbipolar(ParameterInput *pin) {
  if (pin->DoesParameterExist("mhd","ambipolar")) return;
  const bool ohm_eos = pin->DoesParameterExist("mhd","ohmic_resistivity") &&
                       pin->GetString("mhd","ohmic_resistivity").compare("eos") == 0;
  const bool xe_on = pin->DoesParameterExist("mhd","xe_tracer") &&
                     pin->GetBoolean("mhd","xe_tracer");
  if (!ohm_eos || !xe_on) return;
  LogDefault(pin, "<mhd>/ambipolar", "eos (default with ohmic_resistivity = eos and "
             "xe_tracer)");
  if (recipe::NewDefaults(pin)) pin->GetOrAddString("mhd","ambipolar","eos");
}
}  // namespace

namespace mhd {

void DtReportInit(ParameterInput *pin, MeshBlockPack *pp, const bool rss);   // dtdiag-1103, mhd_dtreport.cpp
//----------------------------------------------------------------------------------------
// constructor, initializes data structures and parameters

MHD::MHD(MeshBlockPack *ppack, ParameterInput *pin) :
    pmy_pack(ppack),
    u0("cons",1,1,1,1,1),
    w0("prim",1,1,1,1,1),
    b0("B_fc",1,1,1,1),
    bcc0("B_cc",1,1,1,1,1),
    coarse_u0("ccons",1,1,1,1,1),
    coarse_w0("cprim",1,1,1,1,1),
    coarse_b0("cB_fc",1,1,1,1),
    u1("cons1",1,1,1,1,1),
    b1("B_fc1",1,1,1,1),
    uflx("uflx",1,1,1,1,1),
    efld("efld",1,1,1,1),
    wsaved("wsaved",1,1,1,1,1),
    bccsaved("bccsaved",1,1,1,1,1),
    e3x1("e3x1",1,1,1,1),
    e2x1("e2x1",1,1,1,1),
    e1x2("e1x2",1,1,1,1),
    e3x2("e3x2",1,1,1,1),
    e2x3("e2x3",1,1,1,1),
    e1x3("e1x3",1,1,1,1),
    e1_cc("e1_cc",1,1,1,1),
    e2_cc("e2_cc",1,1,1,1),
    e3_cc("e3_cc",1,1,1,1),
    inner_local("inner",1),
    outer_local("outer",1),
    phi0("phi_fc",1,1,1,1),
    phicc0("phi_cc",1,1,1,1),
    phicc_wb("phi_cc",1,1,1,1),
    phi_wb_x1f("phi_fc",1,1,1,1),
    u0wb("conswb",1,1,1,1,1),
    w0wb("primwb",1,1,1,1,1),
    w0facewb("primfwb",1,1,1,1,1),
    pwb("preswb",1,1,1,1),
    pfacewb("presfwb",1,1,1,1),
    utest("utest",1,1,1,1,1),
    bcctest("bcctest",1,1,1,1,1),
    fofc("fofc",1,1,1,1),
    fofc_cnt("fofc_cnt",1,1,1,1) {
  // cs-open-1018: the angular-momentum module (read only when named)
  // amdef-1019: am_transport defaults on (cs, multi-D sp, fresh dynamic runs)
  {
    // (DefaultOn first: it is false at once when am_transport itself is named)
    const bool amd = hydro::AmModule::DefaultOn(ppack, pin, "mhd");
    if (amd || hydro::AmModule::Named(pin, "mhd")) {
      pam = new hydro::AmModule(ppack, pin, "mhd");
      pam->am_default = amd;
    }
  }
  // Total number of MeshBlocks on this rank to be used in array dimensioning
  int nmb = std::max((ppack->nmb_thispack), (ppack->pmesh->nmb_maxperrank));

  // tracer-1101: fresh-start defaults of the Alfven-floor keys (read by the EOS ctor)
  DefaultAlfvenFloorKeys(ppack, pin);

  // (1) construct EOS object (no default)
  std::string eqn_of_state = pin->GetString("mhd","eos");
  // ideal gas EOS
  if (eqn_of_state.compare("ideal") == 0) {
    if (pmy_pack->pcoord->is_special_relativistic) {
      peos = new IdealSRMHD(ppack, pin);
    } else if (pmy_pack->pcoord->is_general_relativistic) {
      peos = new IdealGRMHD(ppack, pin);
    } else {
      peos = new IdealMHD(ppack, pin);
    }
    nmhd = 5;

  // general EOS (non-relativistic only)
  } else if (eqn_of_state.compare("general") == 0) {
    if (pmy_pack->pcoord->is_special_relativistic ||
        pmy_pack->pcoord->is_general_relativistic) {
      std::cout <<"### FATAL ERROR in "<< __FILE__ <<" at line "<< __LINE__ << std::endl
                <<"<mhd> eos = general cannot be used with SR/GR"<< std::endl;
      std::exit(EXIT_FAILURE);
    }
    peos = new GeneralMHD(ppack, pin);
    nmhd = 5;

  // isothermal EOS
  } else if (eqn_of_state.compare("isothermal") == 0) {
    if (pmy_pack->pcoord->is_special_relativistic ||
        pmy_pack->pcoord->is_general_relativistic) {
      std::cout <<"### FATAL ERROR in "<< __FILE__ <<" at line "<< __LINE__ << std::endl
                <<"<mhd> eos = isothermal cannot be used with SR/GR"<< std::endl;
      std::exit(EXIT_FAILURE);
    } else {
      peos = new IsothermalMHD(ppack, pin);
      nmhd = 4;
    }

  // EOS string not recognized
  } else {
    std::cout <<"### FATAL ERROR in "<< __FILE__ <<" at line "<< __LINE__ << std::endl
              <<"<mhd> eos = '"<< eqn_of_state <<"' not implemented"<< std::endl;
    std::exit(EXIT_FAILURE);
  }

  // (2) Initialize scalars, diffusion, source terms
  nscalars = pin->GetOrAddInteger("mhd","nscalars",0);
  // tracer-1101: fresh-start default of <mhd>/xe_tracer (needs the EOS table)
  DefaultXeTracer(pin, peos->eos_data);
  // adeos-1101: fresh-start default <mhd>/ambipolar = eos (ohmic eos + x_e tracer)
  DefaultAmbipolar(pin);
  // tracer-1101: the passive-tracer slots after the user's scalars (mhd_tracers.hpp):
  // x_e (<mhd>/xe_tracer) first, the floor mass (<mhd>/floor_tracer) last.  Both keys
  // are read only when named; off, nscalars and every path are unchanged.
  tracers = ReadTracerLayout(pin);
  tracers.nuser = nscalars;
  nscalars = tracers.Total();
  xe_tracer = tracers.xe;
  ixe = tracers.xe ? (nmhd + tracers.XeOff()) : -1;
  floor_tracer = tracers.fl;
  ifltr = tracers.fl ? (nmhd + tracers.FlOff()) : -1;
  // duale-1112 (hydro/dual_energy.hpp): every key read only when named
  dual.Read(pin, "mhd", peos->eos_data, nmhd, tracers.DuOff());
  if (dual.retest) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<mhd>/fofc_retest is not implemented (hydro only)"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (floor_tracer) {
    peos->eos_data.ifltr = ifltr;
    if (pin->DoesParameterExist("mhd","floor_tracer_margin")) {
      fltr_margin = pin->GetReal("mhd","floor_tracer_margin");
    }
    if (pin->DoesParameterExist("mhd","floor_tracer_fmax")) {
      fltr_fmax = pin->GetReal("mhd","floor_tracer_fmax");
    }
    if (fltr_margin < 1.0 || fltr_fmax <= 0.0 || fltr_fmax > 1.0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<mhd>/floor_tracer needs floor_tracer_margin >= 1 and "
                << "0 < floor_tracer_fmax <= 1" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (global_variable::my_rank == 0) {
      std::cout << "<mhd>/floor_tracer: floor-added mass in scalar slot " << ifltr
                << " (nscalars " << nscalars << "); drained once per cycle where rho > "
                << fltr_margin << " x the floor density, at most " << fltr_fmax
                << " of the cell mass per cycle" << std::endl;
    }
  }

  // Viscosity (only constructed if needed)
  if (pin->DoesParameterExist("mhd","isotropic_viscosity")) {
    pvisc = new Viscosity("mhd", ppack, pin);
  } else {
    pvisc = nullptr;
  }

  // Resistivity (only constructed if needed)
  // ambipolar-1101: <mhd>/ambipolar (not none) also needs the resistive machinery
  if (pin->DoesParameterExist("mhd","ohmic_resistivity") ||
      (pin->DoesParameterExist("mhd","ambipolar") &&
       pin->GetString("mhd","ambipolar").compare("none") != 0)) {
    presist = new Resistivity(ppack, pin);
  } else {
    presist = nullptr;
  }
  // rkgfuse-1102 <mhd>/tracer_strang (read only when named).  prod1103b default (USER
  // 11-03): true on a FRESH start (recipe::NewDefaults) whenever rkg_strang_fuse is on
  // and a tracer (x_e or floor) exists and the key is not named; recorded (GetOrAdd)
  // for the run's restarts; a restart whose input lacks the key keeps false (bitwise).
  if (!pin->DoesParameterExist("mhd","tracer_strang") && (tracers.xe || tracers.fl) &&
      presist != nullptr && presist->use_rkg_sts && presist->rkg_strang &&
      presist->rkg_strang_fuse) {
    if (recipe::NewDefaults(pin)) {
      pin->GetOrAddBoolean("mhd","tracer_strang",true);
      recipe::Log("<mhd>/tracer_strang = true (default with rkg_strang_fuse and a tracer)");
    } else {
      recipe::Log("<mhd>/tracer_strang = false (" + std::string(
          recipe::GetMode(pin) == recipe::Mode::none ? "recipe none" :
          "restart input has no key") + ": the old default)");
    }
  }
  if (pin->DoesParameterExist("mhd","tracer_strang")) {   // rkgfuse-1102
    trc_strang = pin->GetBoolean("mhd","tracer_strang") && (tracers.xe || tracers.fl);
  }

  // xetracer-1101: advected electron fraction (mhd_tracers.cpp).  Read only when named;
  // off, nscalars and every path are unchanged.
  if (xe_tracer) {
    const bool eos_xe = (presist != nullptr) &&
        ((pin->DoesParameterExist("mhd", "ohmic_resistivity") &&
          pin->GetString("mhd", "ohmic_resistivity").compare("eos") == 0) ||
         presist->ad_type == 2);
    if (!eos_xe || !peos->eos_data.HasElectronFraction()) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<mhd>/xe_tracer = true needs ohmic_resistivity = eos "
                << "(or ambipolar = eos) and the EOS electron fraction table"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    presist->xe_idx = ixe;
    peos->eos_data.ixe = ixe;
    xe_init_pending = true;     // cleared by the restart reader if the file has it
    xe_par.alpha_fac = pin->GetOrAddReal("mhd", "xe_alpha_fac", 1.0);
    xe_par.time_cgs = (ppack->punit != nullptr) ? ppack->punit->time_cgs() : 1.0;
    xe_par.xh = pin->GetOrAddReal("mhd", "eos_xh", 0.7381);
    xe_par.metal_mh = peos->eos_data.tbl.metal_mh;
    xe_par.metal_tcond = pin->GetOrAddReal("mhd", "eos_metal_tcond", 0.0);
    xe_par.metal_cond = pin->GetOrAddBoolean("mhd", "eos_metal_condensation", false);
    xe_par.metal_ion = pin->GetOrAddBoolean("mhd", "eos_metal_ionization", false);
    xetracer::FillDonors(xe_par);
    // the tracer carries y = n_e/n_H (per H nucleus, a per-mass quantity); eta and the
    // ambipolar rho_i convert with n_H = xe_nhfac rho_cgs
    presist->xe_nhfac = xe_par.xh/(eos_cgs::a_hyd*eos_cgs::m_u);
    if (global_variable::my_rank == 0) {
      std::cout << "<mhd>/xe_tracer: advected y = n_e/n_H in scalar " << ixe
                << " (nscalars now "
                << nscalars << "), alpha = " << xe_par.alpha_fac << " x Verner Na+, "
                << "time unit " << xe_par.time_cgs << " s, [M/H] " << xe_par.metal_mh
                << ", condensation " << (xe_par.metal_cond ? "on" : "off") << std::endl;
    }
  }

  // <mhd>/wall_conducting = off | resistive | full (default) (condwall-1028): a REFLECT
  // face is a conducting wall for the field.  The resistive tangential edge EMF eta J_t
  // at a reflect wall grew the wall B.n by 1e3-3e4x, which then carried a Maxwell torque
  // (dLz 3e-5..7e-5, not converging) and, with the resistive Poynting flux, an energy
  // loss (dE/E -1.4 %) on the closed cs shell (csam_1028/RESULTS.md).  `resistive` zeroes
  // only the resistive part (ideal MHD bitwise unchanged); `full` also the ideal edge EMF
  // (the former cs-open-1018 <mhd>/wall_emf_zero, which this key replaces) and the
  // Riemann mass+energy fluxes through the wall.  Default `full` (user 10-29): it
  // changes ideal MHD with reflect walls (the wall B.n is frozen, E closes to round-off).
  if (pin->DoesParameterExist("mhd","wall_emf_zero")) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<mhd>/wall_emf_zero is replaced by <mhd>/wall_conducting "
              << "(off | resistive | full; full = the old wall_emf_zero = true)"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  {
    const std::string wc = pin->GetOrAddString("mhd","wall_conducting","full");
    if (wc == "off") {
      wall_conducting = 0;
    } else if (wc == "resistive") {
      wall_conducting = 1;
    } else if (wc == "full") {
      wall_conducting = 2;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<mhd>/wall_conducting = '" << wc
                << "'; expected off, resistive or full" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    for (int f=0; f<6; ++f) {
      if (ppack->pmesh->mesh_bcs[f] == BoundaryFlag::reflect) any_reflect_face = true;
    }
    if (presist != nullptr) {
      presist->wall_cond_reflect = any_reflect_face && (wall_conducting >= 1);
    }
  }

  // Thermal conduction (only constructed if needed)
  // isotropic_conduction = none: no conduction module at all, exactly as if the key
  // were absent (lets a restart drop it: ck-nq2, the ck two-stream down to the wall)
  if (pin->DoesParameterExist("mhd","isotropic_conduction") &&
      pin->GetString("mhd","isotropic_conduction").compare("none") != 0) {
    if (peos->eos_data.is_ideal) {
      pcond = new Conduction("mhd", ppack, pin);
    } else {
      std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__ << std::endl
                << "Thermal conduction in MHD requires ideal gas EOS" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else {
    pcond = nullptr;
  }

  // Source terms (if needed)
  if (pin->DoesBlockExist("mhd_srcterms")) {
    psrc = new SourceTerms("mhd_srcterms", ppack, pin);
  }

  // (3) read time-evolution option [already error checked in driver constructor]
  // Then initialize memory and algorithms for reconstruction and Riemann solvers
  std::string evolution_t = pin->GetString("time","evolution");

  // allocate memory for conserved and primitive variables
  // With AMR, maximum size of Views are limited by total device memory through an input
  // parameter, which in turn limits max number of MBs that can be created.
  {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    int ncells1 = indcs.nx1 + 2*(indcs.ng);
    int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
    int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
    Kokkos::realloc(u0,   nmb, (nmhd+nscalars), ncells3, ncells2, ncells1);
    Kokkos::realloc(w0,   nmb, (nmhd+nscalars), ncells3, ncells2, ncells1);
    // derived thermodynamic variables are only needed for an expensive (general) EOS
    if (peos->eos_data.IsGeneral()) {
      Kokkos::realloc(wder, nmb, NDERIVED, ncells3, ncells2, ncells1);
      Kokkos::realloc(wtemp, nmb, ncells3, ncells2, ncells1);
    }

    // allocate memory for face-centered and cell-centered magnetic fields
    Kokkos::realloc(bcc0,   nmb, 3, ncells3, ncells2, ncells1);
    Kokkos::realloc(b0.x1f, nmb, ncells3, ncells2, ncells1+1);
    Kokkos::realloc(b0.x2f, nmb, ncells3, ncells2+1, ncells1);
    Kokkos::realloc(b0.x3f, nmb, ncells3+1, ncells2, ncells1);
  }

  // allocate memory for conserved variables on coarse mesh
  if (ppack->pmesh->multilevel) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    int n_ccells1 = indcs.cnx1 + 2*(indcs.ng);
    int n_ccells2 = (indcs.cnx2 > 1)? (indcs.cnx2 + 2*(indcs.ng)) : 1;
    int n_ccells3 = (indcs.cnx3 > 1)? (indcs.cnx3 + 2*(indcs.ng)) : 1;
    Kokkos::realloc(coarse_u0, nmb, (nmhd+nscalars), n_ccells3, n_ccells2, n_ccells1);
    Kokkos::realloc(coarse_w0, nmb, (nmhd+nscalars), n_ccells3, n_ccells2, n_ccells1);
    Kokkos::realloc(coarse_b0.x1f, nmb, n_ccells3, n_ccells2, n_ccells1+1);
    Kokkos::realloc(coarse_b0.x2f, nmb, n_ccells3, n_ccells2+1, n_ccells1);
    Kokkos::realloc(coarse_b0.x3f, nmb, n_ccells3+1, n_ccells2, n_ccells1);
  }

  // allocate boundary buffers for conserved (cell-centered) and face-centered variables
  pbval_u = new MeshBoundaryValuesCC(ppack, pin, false);
  pbval_u->cs_rho_guard = ppack->pmesh->cs_seam_rho_guard;
  pbval_u->cs_lin_resample = ppack->pmesh->cs_seam_resample_linear;
  pbval_u->cs_blockend_ext_ = CsSeamBlockendExtend(ppack, pin);   // layoutbug-1103
  pbval_u->InitializeBuffers((nmhd+nscalars));
  pbval_b = new MeshBoundaryValuesFC(ppack, pin);
  pbval_b->InitializeBuffers(3);
  pbval_b->cs_blockend_ext_ = CsSeamBlockendExtend(ppack, pin);   // layoutbug-1103

  // Orbital advection and shearing box BCs (if requested in input file)
  if (pin->DoesBlockExist("shearing_box")) {
    porb_u = new OrbitalAdvectionCC(ppack, pin, (nmhd+nscalars));
    porb_b = new OrbitalAdvectionFC(ppack, pin);
    psbox_u = new ShearingBoxCC(ppack, pin, (nmhd+nscalars));
    psbox_b = new ShearingBoxFC(ppack, pin);
  } else {
    porb_u = nullptr;
    porb_b = nullptr;
    psbox_u = nullptr;
    psbox_b = nullptr;
  }
        
    // determine if etotgrav and local well-balanced scheme is enabled
    // amdef-1019: default true where supported (mhd::etotgrav_pgen_default)
    if (!pin->DoesParameterExist("mhd","etotgrav") && etotgrav_pgen_default &&
        pin->GetString("mhd","eos") != "isothermal" &&
        !pmy_pack->pmesh->multilevel) {
      if (recipe::NewDefaults(pin)) {
        pin->GetOrAddBoolean("mhd","etotgrav",true);
        recipe::Log("<mhd>/etotgrav = true (default where supported)");
      } else {
        recipe::Log("<mhd>/etotgrav = false (" + std::string(
            recipe::GetMode(pin) == recipe::Mode::none ? "recipe none" :
            "restart input has no key") + ": the old default)");
      }
    }
    use_etotgrav = pin->GetOrAddBoolean("mhd","etotgrav",false);
    if (pin->DoesParameterExist("mhd","rst_eint")) {   // t9rst-1113 (mhd.hpp)
      rst_eint_ideal = pin->GetBoolean("mhd","rst_eint");
    }
    // ckbval-1102 (read only when named; see MHD::SplitBvalFastOK in mhd_splitbval.cpp).
    // prod1103 default (USER 11-02/11-03): true on a FRESH start (recipe::NewDefaults), as
    // the hydro twin; ckbval_1102 + prod_1103 gates: r1m bitwise 40/300 and restart split
    // (on r1m it never takes the energy-only path: 0 fast refreshes, a no-op).  The value is
    // recorded (GetOrAdd); a restart whose input lacks the key keeps false; explicit wins.
    if (!pin->DoesParameterExist("mhd","split_bval_fast")) {
      if (recipe::NewDefaults(pin)) {
        pin->GetOrAddBoolean("mhd","split_bval_fast",true);
        recipe::Log("<mhd>/split_bval_fast = true (default on a fresh start)");
      } else {
        recipe::Log("<mhd>/split_bval_fast = false (" + std::string(
            recipe::GetMode(pin) == recipe::Mode::none ? "recipe none" :
            "restart input has no key") + ": the old default)");
      }
    }
    sbf_on = pin->DoesParameterExist("mhd","split_bval_fast") ?
             pin->GetBoolean("mhd","split_bval_fast") : false;
    if (sbf_on && (pmy_pack->pmesh->multilevel || pin->DoesBlockExist("shearing_box") ||
                   pin->DoesBlockExist("orbital_advection"))) {
      sbf_on = false;
    }
    // SMR/AMR is not supported with etotgrav: the coarse-buffer inversion used for
    // primitive prolongation (bvals/prolong_prims.cpp) reads e + rho*Phi as the internal
    // energy, the conserved prolongation of rho*Phi is removed again with the FINE
    // potential, and a regrid never rebuilds phicc0/phi0 on the new blocks.
    if (use_etotgrav && pmy_pack->pmesh->multilevel) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<mhd>/etotgrav = true is not supported with mesh "
                << "refinement (SMR/AMR): the potential energy is not handled at level "
                << "boundaries or on regrid." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    use_wellbalance_dynamic = pin->GetOrAddBoolean("mhd","wellbalance_dynamic",false);
    use_wb_x1 = pin->GetOrAddBoolean("mhd","wb_x1",false);
    // rebuild the well-balanced background cache at stage 1 of every N-th cycle instead
    // of every stage (0 = every stage).  A balanced state is unaffected; an evolving one
    // reconstructs against a background up to N cycles old, which is nothing against the
    // radiative timescale the profile changes on.  N = 1 halves the walk cost, N = 10
    // removes it.
    wb_cache_every = pin->GetOrAddInteger("mhd","wb_cache_every",0);
    // Ideal-gas x1 WB cache (BuildWBCache) safety keys, BSG arm 1 crash (10-02):
    //  <mhd>/wb_clamp_ideal: polytropic closed form limits the T DECREASE over a
    //    segment to t0/2, as the general-EOS walk does (getWBerho, clamp_dt);
    //  <mhd>/wb_guard_ideal: WBGuard flattens a non-finite, non-positive or absurd
    //    walked stencil to its anchor (backstop).
    // Both DEFAULT true on a fresh run where that path is active (wb_x1 with a
    // non-general EOS); a restart whose input lacks a key keeps the old default (off,
    // one line); an explicit value always wins.
    const bool wbi_path = use_wb_x1 && !peos->eos_data.IsGeneral();
    auto wbi_default = [&](const char *key) {
      if (!wbi_path || pin->DoesParameterExist("mhd", key)) return false;
      if (!global_variable::restart_run) return true;
      if (global_variable::my_rank == 0) {
        std::cout << "mhd: restart input has no mhd/" << key << ": keeping the old "
                  << "default (off); set it true to switch" << std::endl;
      }
      return false;
    };
    wb_guard_ideal = pin->GetOrAddBoolean("mhd","wb_guard_ideal",
                                          wbi_default("wb_guard_ideal"));
    wb_clamp_ideal = pin->GetOrAddBoolean("mhd","wb_clamp_ideal",
                                          wbi_default("wb_clamp_ideal"));
    use_wb_x2 = pin->GetOrAddBoolean("mhd","wb_x2",false);
    use_wb_rho = pin->GetOrAddBoolean("mhd","wb_rho",false);
    // wallwb-1105 (mhd.hpp): as <hydro>/wall_wb_onesided, fresh-start default with WB x1
    if (!pin->DoesParameterExist("mhd","wall_wb_onesided") &&
        use_wellbalance_dynamic && use_wb_x1) {
      if (recipe::NewDefaults(pin)) {
        pin->GetOrAddBoolean("mhd","wall_wb_onesided",true);
        recipe::Log("<mhd>/wall_wb_onesided = true (default with wellbalance_dynamic + "
                    "wb_x1; acts at closed / reflecting x1 walls)");
      } else {
        recipe::Log("<mhd>/wall_wb_onesided = false (" + std::string(
            recipe::GetMode(pin) == recipe::Mode::none ? "recipe none" :
            "restart input has no key") + ": the old wall cell)");
      }
    }
    if (pin->DoesParameterExist("mhd","wall_wb_onesided")) {
      wall_wb_onesided = pin->GetBoolean("mhd","wall_wb_onesided");
    }
    // switch the x1 well-balanced reconstruction off above this radius (0 = never); see
    // the declaration in mhd.hpp.  Only the position-aware x1 path (spherical polar and
    // cubed sphere, where x1v is the radius) honours it -- a Cartesian mesh has no radius
    // to compare against.  Mirror of <hydro>/wb_rmax.
    wb_rmax = pin->GetOrAddReal("mhd","wb_rmax",0.0);
    if (pin->DoesParameterExist("mhd","wb_poison_test")) wb_poison = pin->GetBoolean("mhd","wb_poison_test");
    // ...and off BELOW this radius (0 = never), same path, same caveat.
    wb_rmin = pin->GetOrAddReal("mhd","wb_rmin",0.0);
    // radiation-modified CFL signal speed (mhd_newdt.cpp), the MHD form of
    // <hydro>/rad_signal_speed with the same default: on with an IMPLICIT M1 transport
    // (<rad_m1>/transport = implicit | implicit_x1), off otherwise.  Read only when
    // named, so the parameter dump of a run without the key is unchanged.
    if (pin->DoesBlockExist("rad_m1") &&
        pin->DoesParameterExist("rad_m1","transport")) {
      std::string tr = pin->GetString("rad_m1","transport");
      rad_signal_speed = (tr.compare("implicit") == 0 || tr.compare("implicit_x1") == 0);
    }
    if (pin->DoesParameterExist("mhd","rad_signal_speed")) {
      rad_signal_speed = pin->GetBoolean("mhd","rad_signal_speed");
    }
    // dtdiag-1103: <mhd>/dt_report (mhd_dtreport.cpp); read only when named
    DtReportInit(pin, ppack, rad_signal_speed);
    // allocate array of flags used with etotgrav
    if (use_etotgrav || use_wellbalance_dynamic) {
      auto &indcs = pmy_pack->pmesh->mb_indcs;
      int ncells1 = indcs.nx1 + 2*(indcs.ng);
      int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
      int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
      Kokkos::realloc(phicc0, nmb, ncells3, ncells2, ncells1);
      Kokkos::realloc(phi0.x1f, nmb, ncells3, ncells2, ncells1+1);
      Kokkos::realloc(phi0.x2f, nmb, ncells3, ncells2+1, ncells1);
      Kokkos::realloc(phi0.x3f, nmb, ncells3+1, ncells2, ncells1);
    if (use_wellbalance_dynamic) {
      // 5 stencil states x 3 channels (d, e, p) of the local hydrostatic background
      Kokkos::realloc(wbq0, nmb, 15, ncells3, ncells2, ncells1);
    }
    }
    // the x1 well-balanced scheme sees the TRUE potential unless a problem generator
    // calls EnableWBEffectivePotential(): shallow copies, the very same allocations
    // (the hydro twin, hydro.cpp); the default path is bitwise unchanged
    if (!use_phi_wb) {
      phicc_wb = phicc0;
      phi_wb_x1f = phi0.x1f;
    }
    if (use_wellbalance_dynamic) {
      // select well-balanced scheme assumption (no default).  Test for compatibility of options
      std::string wb_opt = pin->GetString("mhd","wb_option");
      if (wb_opt.compare("isothermal") == 0) {
        wb_option = WBOption::isothermal;
      } else if (wb_opt.compare("isodensity") == 0) {
        wb_option = WBOption::isodensity;
      } else if (wb_opt.compare("isentropic") == 0) {
        wb_option = WBOption::isentropic;
      } else if (wb_opt.compare("adaptive") == 0) {
        wb_option = WBOption::adaptive;
      } else if (wb_opt.compare("polytropic") == 0) {
        wb_option = WBOption::polytropic;
      } else if (wb_opt.compare("isentropic_dt") == 0) {
        wb_option = WBOption::isentropic_dt;
      } else if (wb_opt.compare("adaptive_fast") == 0) {
        wb_option = WBOption::adaptive_fast;
      // Error for anything else
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<mhd> wb_option = '" << wb_opt
                  << "' not implemented for mhd" << std::endl;
        std::exit(EXIT_FAILURE);
      }
        // wball-1104 <mhd>/wb_dirs = x1 (default) | all (read only when named)
        if (pin->DoesParameterExist("mhd","wb_dirs")) {
          const std::string wd = pin->GetString("mhd","wb_dirs");
          if (wd.compare("all") == 0) {
            wb_all = true;
          } else if (wd.compare("x1") != 0) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "<mhd>/wb_dirs must be x1 or all" << std::endl;
            std::exit(EXIT_FAILURE);
          }
        }
        if (wb_all) {
          if (!use_wb_x1 || use_wb_x2 || pmy_pack->pmesh->one_d) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "<mhd>/wb_dirs = all needs wb_x1 = true, wb_x2 = "
                      << "false and a 2-D/3-D mesh" << std::endl;
            std::exit(EXIT_FAILURE);
          }
          wb_weight = pin->GetOrAddBoolean("mhd","wb_weight",true);
          wbw_r0 = pin->GetOrAddReal("mhd","wb_w_r0",0.1);
          wbw_r1 = pin->GetOrAddReal("mhd","wb_w_r1",0.3);
          wbw_m0 = pin->GetOrAddReal("mhd","wb_w_m0",0.2);
          wbw_m1 = pin->GetOrAddReal("mhd","wb_w_m1",0.5);
          wbw_eta = pin->GetOrAddReal("mhd","wb_w_eta",1.0e-4);
          auto &ic = pmy_pack->pmesh->mb_indcs;
          const int c1 = ic.nx1 + 2*ic.ng;
          const int c2 = (ic.nx2 > 1) ? (ic.nx2 + 2*ic.ng) : 1;
          const int c3 = (ic.nx3 > 1) ? (ic.nx3 + 2*ic.ng) : 1;
          const int nmbw = std::max((pmy_pack->nmb_thispack), (pmy_pack->pmesh->nmb_maxperrank));
          Kokkos::realloc(wbq2, nmbw, 15, c3, c2, c1);
          if (c3 > 1) Kokkos::realloc(wbq3, nmbw, 15, c3, c2, c1);
          Kokkos::realloc(wbw, nmbw, 7, c3, c2, c1);
          Kokkos::deep_copy(wbw, 1.0);
          if (pin->DoesParameterExist("mhd","wb_seam") && pin->GetBoolean("mhd","wb_seam")) {
            if (!pmy_pack->pmesh->use_cubed_sphere || pmy_pack->pmesh->multilevel) {
              std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                        << std::endl << "<mhd>/wb_seam needs a single-level cubed sphere"
                        << std::endl;
              std::exit(EXIT_FAILURE);
            }
            pbval_u->wb_seam = true;
            pbval_u->wbs_phi = phicc0;
            pbval_u->wbs_eos = peos->eos_data;
            pbval_u->wbs_opt = (wb_option == WBOption::polytropic) ? 3 : 1;
            pbval_u->wbs_etg = use_etotgrav;
            pbval_u->wbs_ien = IEN;
            pbval_u->wbs_bcc = bcc0;
          }
          if (global_variable::my_rank == 0) {
            std::cout << "mhd: wb_dirs = all (x2/x3 dynamic well balance), wb_weight = "
                      << (wb_weight ? "true" : "false") << std::endl;
          }
        }
        // select well-balanced scheme direction
        if (!use_wb_x1 && !use_wb_x2) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<mhd> wb_direction not set!" << std::endl;
          std::exit(EXIT_FAILURE);
        }
    }
          
    // determine if static wellbalance is enabled
    use_wellbalance_static = pin->GetOrAddBoolean("mhd","wellbalance_static",false);
    use_wellbalance_static_reconst_perturb = pin->GetOrAddBoolean("mhd","wellbalance_static_reconst",false);
  // the deviation reconstruction reads the background arrays, which only
  // wellbalance_static allocates -- the pair is one scheme, not two options
  if (use_wellbalance_static_reconst_perturb && !use_wellbalance_static) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<mhd>/wellbalance_static_reconst needs "
              << "<mhd>/wellbalance_static = true" << std::endl;
    std::exit(EXIT_FAILURE);
  }
    // allocate array of flags used with wellbalance
    if (use_wellbalance_static) {
      auto &indcs = pmy_pack->pmesh->mb_indcs;
      int ncells1 = indcs.nx1 + 2*(indcs.ng);
      int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
      int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
      Kokkos::realloc(u0wb, nmb, nmhd, ncells3, ncells2, ncells1);
      Kokkos::realloc(w0wb, nmb, nmhd, ncells3, ncells2, ncells1);
      Kokkos::realloc(w0facewb.x1f, nmb, nmhd, ncells3, ncells2, ncells1+1);
      Kokkos::realloc(w0facewb.x2f, nmb, nmhd, ncells3, ncells2+1, ncells1);
      Kokkos::realloc(w0facewb.x3f, nmb, nmhd, ncells3+1, ncells2, ncells1);
    // background pressure: one scalar per cell and per face, filled once (see
    // SetWbBackgroundPressure) instead of being re-derived from the background every
    // stage
    Kokkos::realloc(pwb, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(pfacewb.x1f, nmb, ncells3, ncells2, ncells1+1);
    Kokkos::realloc(pfacewb.x2f, nmb, ncells3, ncells2+1, ncells1);
    Kokkos::realloc(pfacewb.x3f, nmb, ncells3+1, ncells2, ncells1);
    }
        
    bool use_polar_boundary = pin->GetOrAddBoolean("mesh", "use_polar_boundary", false);
    if (use_polar_boundary) {
      auto &indcs = pmy_pack->pmesh->mb_indcs;
      int ncells1 = indcs.nx1 + 2*(indcs.ng);
      Kokkos::realloc(inner_local, ncells1);
      Kokkos::realloc(outer_local, ncells1);
    }
        
  // for time-evolving problems, continue to construct methods, allocate arrays
  if (evolution_t.compare("stationary") != 0) {
    // determine if FOFC is enabled
    // default false; true for the problem generators that set fofc_pgen_default
    // (mhd.hpp), except on a restart whose file lacks the key (old default kept)
    bool fofc_dflt = false;
    if (fofc_pgen_default && !pin->DoesParameterExist("mhd","fofc")) {
      if (global_variable::restart_run) {
        if (global_variable::my_rank == 0) {
          std::cout << "mhd: restart input has no mhd/fofc: keeping the old default "
                    << "(off); set it true to switch" << std::endl;
        }
      } else {
        fofc_dflt = true;
      }
    }
    use_fofc = pin->GetOrAddBoolean("mhd","fofc",fofc_dflt);
    if (use_fofc) {fofc_rep.Init(pin, "mhd");}
    // <mhd>/fofc_rsolver: see the note in mhd.hpp.  Default llf keeps every existing
    // run bit-for-bit unchanged.
    {
      std::string fofc_rs = pin->GetOrAddString("mhd","fofc_rsolver","llf");
      if (fofc_rs.compare("llf") == 0) {
        fofc_hlle = false;
      } else if (fofc_rs.compare("hlle") == 0) {
        fofc_hlle = true;
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<mhd>/fofc_rsolver = '" << fofc_rs
                  << "' not implemented; choose llf or hlle" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (fofc_hlle && (pmy_pack->pcoord->is_special_relativistic ||
                        pmy_pack->pcoord->is_general_relativistic)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<mhd>/fofc_rsolver=hlle is implemented only for "
                  << "NEWTONIAN MHD: there is no single-state HLLE for SR/GR MHD"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    // Diagnostic: the plain four-face average in place of the GS05 upwind corner EMF.
    // The GS05 derivative terms are fed by the transverse-face EMFs, which DO see the
    // (-1)^(i+k) CT null mode that the cell-centred field cannot; in a radial outflow
    // that biases the corner EMF downwind and grows the mode at ~v_r/dr in a polar row.
    use_bs_emf = pin->GetOrAddBoolean("mhd","bs_emf",false);
    // lhlld2-1102: smooth GS07 upwind weights in the corner EMF (EmfSmoothW); read only
    // when named, so key-absent runs keep the sign switch bitwise
    if (pin->DoesParameterExist("mhd","emf_smooth")) {
      emf_smooth = pin->GetReal("mhd","emf_smooth");
      if (emf_smooth < 0.0) emf_smooth = -1.0;
    }
    // <mhd>/emf_scheme = gs07 (default) | uct | uct_hlld, read only when named
    // (lhlld2-1102, ucthlld-1102).  Both UCT variants use the Mignone & Del Zanna 2021
    // composition formula (MD21 eq. 35) for the corner EMF, inside the GS07 corner kernel
    // (MHD::CornerE), with per-face coefficients stored by the flux sweeps (UctFaceStore,
    // mhd_fluxes.cpp): uct = UCT-HLL (MD21 eq. 30 coefficients, any of hlle/hlld/lhlld),
    // uct_hlld = UCT-HLLD (MD21 eqs. 40, 44-46; rsolver hlld or lhlld).  Cartesian 2-D/3-D,
    // the cubed sphere (gnomonic face frames, the existing seam/vertex EMF averaging) and
    // spherical polar (spdef-1103: orthonormal face frames = the Cartesian branch of the
    // shared UctE1/E2/E3; polar_hlle_rows faces take the HLL coefficients, the polar-row
    // polar_emf_diss term and the polar E_r azimuthal average act on the UCT EMF as on GS07).
    // prod1103 default (USER 11-03): emf_scheme = uct_hlld on a FRESH start
    // (recipe::NewDefaults) when rsolver = hlld or lhlld and the key is not named, on
    // every mesh UCT-HLLD supports (2-D/3-D Cartesian, cubed sphere; spherical polar:
    // UCT_SP_DEFAULT below, on since spdef-1103), Newtonian, ideal/general EOS, no fofc / emf_smooth / bs_emf,
    // a uniform mesh (UCT is not validated with refinement); elsewhere gs07, logged with
    // the reason.  The value is recorded (GetOrAdd) so the run's restarts carry it.  A
    // restart whose input lacks the key keeps gs07 (bitwise); an explicit key wins.
    if (!pin->DoesParameterExist("mhd","emf_scheme") &&
        pin->DoesParameterExist("mhd","rsolver") &&
        (pin->GetString("mhd","rsolver") == "hlld" ||
         pin->GetString("mhd","rsolver") == "lhlld")) {
      constexpr bool UCT_SP_DEFAULT = true;   // spdef-1103: sp port gated (prod_1103/sp)
      std::string no;
      if (ppack->pmesh->one_d) {
        no = "1-D (no corner EMF)";
      } else if (ppack->pmesh->use_spherical_polar && !UCT_SP_DEFAULT) {
        no = "spherical polar not supported";
      } else if (ppack->pcoord->is_special_relativistic ||
                 ppack->pcoord->is_general_relativistic ||
                 ppack->pcoord->is_dynamical_relativistic) {
        no = "relativistic MHD";
      } else if (!(peos->eos_data.is_ideal || peos->eos_data.IsGeneral())) {
        no = "isothermal EOS";
      } else if (use_fofc || emf_smooth >= 0.0 || use_bs_emf) {
        no = "fofc, emf_smooth or bs_emf set";
      } else if (ppack->pmesh->multilevel) {
        no = "mesh refinement";
      }
      if (!recipe::NewDefaults(pin)) {
        recipe::Log("<mhd>/emf_scheme = gs07 (" + std::string(
            recipe::GetMode(pin) == recipe::Mode::none ? "recipe none" :
            "restart input has no key") + ": the old default)");
      } else if (!no.empty()) {
        recipe::Log("<mhd>/emf_scheme = gs07 (uct_hlld default not applied: " + no + ")");
      } else {
        pin->GetOrAddString("mhd","emf_scheme","uct_hlld");
        recipe::Log("<mhd>/emf_scheme = uct_hlld (default with rsolver = " +
                    pin->GetString("mhd","rsolver") + ")");
      }
    }
    if (pin->DoesParameterExist("mhd","emf_scheme")) {
      const std::string es = pin->GetString("mhd","emf_scheme");
      if (es == "uct") {
        emf_uct = true;
      } else if (es == "uct_hlld") {
        emf_uct_hlld = true;
      } else if (es != "gs07") {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<mhd>/emf_scheme = '" << es << "' (gs07 | uct | "
                  << "uct_hlld)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if ((emf_uct || emf_uct_hlld) && (ppack->pmesh->one_d ||
                      ppack->pcoord->is_special_relativistic ||
                      ppack->pcoord->is_general_relativistic ||
                      ppack->pcoord->is_dynamical_relativistic ||
                      !(peos->eos_data.is_ideal || peos->eos_data.IsGeneral()) ||
                      use_fofc || emf_smooth >= 0.0 || use_bs_emf)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<mhd>/emf_scheme = " << es << " is implemented only "
                  << "for 2-D/3-D Cartesian, cubed-sphere or spherical-polar Newtonian MHD "
                  << "with an ideal or general EOS, without fofc, emf_smooth or bs_emf"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      // uctaudit-1103: UCT is not validated with mesh refinement (the fresh-start default
      // already skips it there); an explicit emf_scheme = uct | uct_hlld is refused too
      if ((emf_uct || emf_uct_hlld) && ppack->pmesh->multilevel) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<mhd>/emf_scheme = " << es << " is not supported with "
                  << "mesh refinement (SMR/AMR); use emf_scheme = gs07" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      // uctaudit-1103: emf_scheme = uct builds the corner EMF from the HLL coefficients
      // (MD21 eq. 30), which have no Alfven fan, so the lhlld_alfmach_rot scaling of the
      // fan dissipation (a TEST key) would act on the face fluxes only and be silently
      // absent from CT.  uct_hlld carries it (UctFaceStore); gs07 reads the face EMFs.
      if (emf_uct && peos->eos_data.alfmach_rot) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<mhd>/lhlld_alfmach_rot is not seen by the "
                  << "emf_scheme = uct corner EMF (HLL coefficients); use emf_scheme = "
                  << "uct_hlld or gs07" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    // layoutaudit-1103 (USER 11-03, from vertexacc-1103): the cell-centred cube-vertex
    // fill <mesh>/cs_vertex_fill_cc = true with the GS07 corner EMF (named, or the
    // fallback when uct / uct_hlld is not selected) or with MHD fofc collapses dt to
    // ~1e-28 at the cube vertex at low beta, or goes NaN by cycle 2.  Refused.
    // cornerfix-1103: with <mesh>/cs_vertex_3cell (fresh default) no EMF reads the corner
    // block any more (3-cell vertex rule, mirrored near-vertex GS07 edges); the
    // combination then runs (gates: cornerfix_1103/FATAL.txt) and is only noted.
    if (ppack->pmesh->cs_vertex_fill_cc && (!(emf_uct || emf_uct_hlld) || use_fofc) &&
        ppack->pmesh->cs_vertex_3cell) {
      if (global_variable::my_rank == 0) {
        std::cout << "### WARNING <mesh>/cs_vertex_fill_cc = true with "
                  << (use_fofc ? "<mhd>/fofc" : "the GS07 corner EMF")
                  << ": allowed under <mesh>/cs_vertex_3cell = true (the corner block is "
                  << "not read); refused without it" << std::endl;
      }
    } else if (ppack->pmesh->cs_vertex_fill_cc &&
               (!(emf_uct || emf_uct_hlld) || use_fofc)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<mesh>/cs_vertex_fill_cc = true is not supported with "
                << (use_fofc ? "<mhd>/fofc = true" : "the GS07 corner EMF (<mhd>/"
                    "emf_scheme = gs07, named or by fallback)")
                << ": the vertex dt collapses (~1e-28 at low beta) or the run goes NaN by "
                << "cycle 2 (vertexacc-1103).  Use emf_scheme = uct or uct_hlld without "
                << "fofc, or leave cs_vertex_fill_cc off." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // Fix for the above at the source: in the two POLAR cell rows (the only place the
    // (-1)^(i+k) mode closes on itself, because the polar ghost row is the k+N/2 image of
    // the row) add the UCT/Rusanov face-field dissipation to the corner e2, with the local
    // fast speed. Damps the face checkerboard at ~2 c_f/dx, which bcc-based upwinding
    // cannot do at all.
    // Default ON wherever there is a polar boundary: verified on CPU (bench/sp_cpu_ablate/
    // diss, 50 G -> 3e-4 G in 300 s) and on the GPU from scratch (bench/sp_diss, clean to
    // rot 1.13 where the control jumped 100x at rot 0.8).
  // wbplm-1103 (read only when named): see mhd.hpp top_ghost_prim_copy
  if (pin->DoesParameterExist("mhd","top_ghost_prim_copy")) {
    top_ghost_prim_copy = pin->GetBoolean("mhd","top_ghost_prim_copy");
  }
    use_polar_emf_diss = pin->GetOrAddBoolean("mhd","polar_emf_diss",
                                              pmy_pack->pmesh->use_polar_boundary);
    // On the cubed sphere the scheme is UNSTABLE below beta ~ 0.05 and the fallback is
    // what stops it, so this defaults ON there -- there is no correct answer below that
    // beta to preserve.  Off (0) by default on every other grid; share-1104 G8: a named
    // value > 0 acts there too (same face test and HLLE flux, mhd_fluxes.cpp).
    cs_lowbeta_fallback = pin->GetOrAddReal("mhd","cs_lowbeta_fallback",
                     (pmy_pack->pmesh->use_cubed_sphere ? 0.5 : 0.0));
    {
      // per-sweep mask; accept the earlier boolean spelling too
      std::string s = pin->GetOrAddString("mhd","polar_hlle_rows","6");
      if (s == "true") polar_hlle_rows = 7;
      else if (s == "false") polar_hlle_rows = 0;
      else polar_hlle_rows = std::stoi(s);
    }

    // select reconstruction method (default PLM)
    std::string xorder = pin->GetOrAddString("mhd","reconstruct","plm");
    if (xorder.compare("dc") == 0) {
      recon_method = ReconstructionMethod::dc;
    } else if (xorder.compare("plm") == 0) {
      recon_method = ReconstructionMethod::plm;
      // check that nghost > 2 with PLM+FOFC
      auto &indcs = pmy_pack->pmesh->mb_indcs;
      if (use_fofc && indcs.ng < 3) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "FOFC and " << xorder << " reconstruction requires at "
          << "least 3 ghost zones, but <mesh>/nghost=" << indcs.ng << std::endl;
        std::exit(EXIT_FAILURE);
      }
    } else if (xorder.compare("ppm4") == 0 ||
               xorder.compare("ppmx") == 0 ||
               xorder.compare("wenoz") == 0) {
      // check that nghost > 2
      auto &indcs = pmy_pack->pmesh->mb_indcs;
      if (indcs.ng < 3) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << xorder << " reconstruction requires at least 3 ghost zones, "
          << "but <mesh>/nghost=" << indcs.ng << std::endl;
        std::exit(EXIT_FAILURE);
      }
      // check that nghost > 3 with PPM4(or PPMX or WENOZ)+FOFC
      if (use_fofc && indcs.ng < 4) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "FOFC and " << xorder << " reconstruction requires at "
          << "least 4 ghost zones, but <mesh>/nghost=" << indcs.ng << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (xorder.compare("ppm4") == 0) {
        recon_method = ReconstructionMethod::ppm4;
      } else if (xorder.compare("ppmx") == 0) {
        recon_method = ReconstructionMethod::ppmx;
      } else if (xorder.compare("wenoz") == 0) {
        recon_method = ReconstructionMethod::wenoz;
      }
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<mhd>/recon = '" << xorder << "' not implemented"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // kriel-1103 <mhd>/uct_corner_recon = plm (default, bitwise) | match, read only when
    // named: match = the UCT corner reconstructions along the face (staggered B_q, the
    // half-cell change of v.n_p; mhd_corner_e.cpp UctE3M/E1M/E2M) use <mhd>/reconstruct
    // (ppm4 | ppmx | wenoz) instead of PLM.  emf_scheme = uct | uct_hlld, Cartesian only.
    if (pin->DoesParameterExist("mhd","uct_corner_recon")) {
      const std::string ucr = pin->GetString("mhd","uct_corner_recon");
      if (ucr == "match") {
        if (!(emf_uct || emf_uct_hlld) || pmy_pack->pmesh->use_cubed_sphere ||
            pmy_pack->pmesh->use_spherical_polar) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<mhd>/uct_corner_recon = match needs emf_scheme = uct"
                    << " | uct_hlld on a Cartesian mesh" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (recon_method == ReconstructionMethod::ppm4) {
          uct_corner_rec = 1;
        } else if (recon_method == ReconstructionMethod::ppmx) {
          uct_corner_rec = 2;
        } else if (recon_method == ReconstructionMethod::wenoz) {
          uct_corner_rec = 3;
        }   // plm / dc: the PLM corner (0)
      } else if (ucr != "plm") {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<mhd>/uct_corner_recon = '" << ucr << "' (plm | match)"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    // select RADIAL (x1) reconstruction on the cubed sphere / spherical polar (default
    // PLM, i.e. today's GridPiecewiseLinearX1, bitwise unchanged). `reconstruct` above
    // never reaches the x1 sweep on those two grids -- see the note on str_r1_ in
    // mhd_fluxes.cpp -- so this is a SEPARATE knob, not a fallback for it.
    {
      bool cs_or_sp = pmy_pack->pmesh->use_cubed_sphere ||
                      pmy_pack->pmesh->use_spherical_polar;
      std::string xorder1 = pin->GetOrAddString("mhd","reconstruct_x1","plm");
      auto &indcs = pmy_pack->pmesh->mb_indcs;
      if (xorder1.compare("plm") == 0) {
        recon_method_x1 = ReconstructionMethod::plm;
      } else if (xorder1.compare("ppmx") == 0) {
        // A. Mignone, JCP 270, 784 (2014) [arXiv:1404.0537] gives the published
        // curvilinear (volume-coordinate) generalization of PPM (his Eqs. 45-48,
        // implemented in reconstruct/mignone_curvilinear.hpp) but ONLY for the
        // CW84-style limiter, not for the Colella & Sekora extremum-preserving one
        // that this code's uniform-grid `ppmx` uses. Inventing that generalization
        // is out of scope, so refuse rather than silently fall back to something
        // unpublished; use reconstruct_x1=ppm4 (the curvilinear CW84-limiter PPM)
        // instead.
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "<mhd>/reconstruct_x1=ppmx is not implemented: Mignone (2014) "
          << "gives no curvilinear Colella-Sekora extremum-preserving PPM limiter; use "
          << "reconstruct_x1=ppm4 (curvilinear CW84 limiter) or plm" << std::endl;
        std::exit(EXIT_FAILURE);
      } else if (xorder1.compare("ppm4") == 0 ||
                 xorder1.compare("wenoz") == 0) {
        // The published, non-uniform-grid-EXACT reconstruction from A. Mignone, "High-
        // order conservative reconstruction schemes for finite volume methods in
        // cylindrical and spherical coordinates", JCP 270, 784 (2014), arXiv:1404.0537:
        // ppm4 = his curvilinear 4th-order PPM (Sec. 3.3, Eqs. 45-48); wenoz = his
        // curvilinear WENO3 (Sec. 3.2, third order -- there is no published curvilinear
        // WENO beyond third order, so this is a genuine order REDUCTION from the
        // uniform-grid 5th-order WENO-Z used elsewhere). See
        // reconstruct/mignone_curvilinear.hpp for the implementation and equation
        // references, and tests_cs_radial_recon/ for the measured order.
        if (indcs.ng < 3) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "<mhd>/reconstruct_x1=" << xorder1 << " requires at least 3 "
            << "ghost zones, but <mesh>/nghost=" << indcs.ng << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (use_fofc && indcs.ng < 4) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "FOFC and <mhd>/reconstruct_x1=" << xorder1 << " require at "
            << "least 4 ghost zones, but <mesh>/nghost=" << indcs.ng << std::endl;
          std::exit(EXIT_FAILURE);
        }
        // The Kappeli & Mishra dynamic radial well-balancing (GridPiecewiseLinearX1's
        // background-subtracted branch) is a PLM-specific scheme (it subtracts a
        // background then re-applies PLM_nonuniform to the deviation); there is no
        // ppm4/ppmx/wenoz counterpart, so refuse the combination rather than silently
        // dropping the well-balancing.
        if (use_wellbalance_dynamic && use_wb_x1) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "<mhd>/reconstruct_x1=" << xorder1 << " is not implemented "
            << "with <mhd>/wellbalance_dynamic + wb_x1 (the K&M radial well-balancing is "
            << "PLM-only); use <mhd>/reconstruct_x1=plm or turn wellbalance_dynamic off"
            << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (xorder1.compare("ppm4") == 0) {
          recon_method_x1 = ReconstructionMethod::ppm4;
        } else {
          recon_method_x1 = ReconstructionMethod::wenoz;
        }
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<mhd>/reconstruct_x1 = '" << xorder1
                  << "' not implemented" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (cs_or_sp && global_variable::my_rank == 0) {
        std::cout << "Radial (x1) reconstruction on the "
                  << (pmy_pack->pmesh->use_cubed_sphere ? "cubed sphere" : "spherical "
                      "polar") << " grid (<mhd>/reconstruct_x1) = " << xorder1
                  << std::endl;
      }
    }

    // select Riemann solver (no default).  Test for compatibility of options
    std::string rsolver = pin->GetString("mhd","rsolver");
    // Special relativistic solvers
    if (pmy_pack->pcoord->is_special_relativistic) {
      if (evolution_t.compare("dynamic") == 0) {
        if (rsolver.compare("llf") == 0) {
          rsolver_method = MHD_RSolver::llf_sr;
        } else if (rsolver.compare("hlle") == 0) {
          rsolver_method = MHD_RSolver::hlle_sr;
        // Error for anything else
        } else {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<mhd> rsolver = '" << rsolver << "' not implemented"
                    << " for SR dynamics" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "kinematic dynamics not implemented for SR" <<std::endl;
        std::exit(EXIT_FAILURE);
      }

    // General relativistic solvers
    } else if (pmy_pack->pcoord->is_general_relativistic) {
      if (evolution_t.compare("dynamic") == 0) {
        if (rsolver.compare("llf") == 0) {
          rsolver_method = MHD_RSolver::llf_gr;
        } else if (rsolver.compare("hlle") == 0) {
          rsolver_method = MHD_RSolver::hlle_gr;
        // Error for anything else
        } else {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<mhd> rsolver = '" << rsolver << "' not implemented"
                    << " for GR dynamics" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "kinematic dynamics not implemented for GR" <<std::endl;
        std::exit(EXIT_FAILURE);
      }

    // Non-relativistic dynamic solvers
    } else if (evolution_t.compare("dynamic") == 0) {
      // LLF solver
      if (rsolver.compare("llf") == 0) {
        rsolver_method = MHD_RSolver::llf;
      // HLLE solver
      } else if (rsolver.compare("hlle") == 0) {
        rsolver_method = MHD_RSolver::hlle;
      // HLLD solver
      } else if (rsolver.compare("hlld") == 0) {
        rsolver_method = MHD_RSolver::hlld;
      // LHLLD solver (HLLD with the low-Mach fix of Minoshima & Miyoshi 2021)
      } else if (rsolver.compare("lhlld") == 0) {
        rsolver_method = MHD_RSolver::lhlld;
      // Roe solver
      // } else if (rsolver.compare("roe") == 0) {
      //   rsolver_method = MHD_RSolver::roe;
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<mhd>/rsolver = '" << rsolver << "' not implemented"
                  << " for dynamic problems" << std::endl;
        std::exit(EXIT_FAILURE);
      }

    // Non-relativistic kinematic solver
    } else {
      // Advect solver
      if (rsolver.compare("advect") == 0) {
        rsolver_method = MHD_RSolver::advect;
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<mhd>/rsolver = '" << rsolver << "' not implemented"
                  << " for kinematic problems" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    // Final memory allocations
    {
      // allocate second registers
      auto &indcs = pmy_pack->pmesh->mb_indcs;
      int ncells1 = indcs.nx1 + 2*(indcs.ng);
      int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
      int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
      Kokkos::realloc(u1,     nmb, (nmhd+nscalars), ncells3, ncells2, ncells1);
      Kokkos::realloc(b1.x1f, nmb, ncells3, ncells2, ncells1+1);
      Kokkos::realloc(b1.x2f, nmb, ncells3, ncells2+1, ncells1);
      Kokkos::realloc(b1.x3f, nmb, ncells3+1, ncells2, ncells1);

      // allocate fluxes, electric fields
      Kokkos::realloc(uflx.x1f, nmb, (nmhd+nscalars), ncells3, ncells2, ncells1+1);
      Kokkos::realloc(uflx.x2f, nmb, (nmhd+nscalars), ncells3, ncells2+1, ncells1);
      Kokkos::realloc(uflx.x3f, nmb, (nmhd+nscalars), ncells3+1, ncells2, ncells1);
      Kokkos::realloc(efld.x1e, nmb, ncells3+1, ncells2+1, ncells1);
      Kokkos::realloc(efld.x2e, nmb, ncells3+1, ncells2, ncells1+1);
      Kokkos::realloc(efld.x3e, nmb, ncells3, ncells2+1, ncells1+1);

      // allocate scratch arrays for face- and cell-centered E used in CornerE
      Kokkos::realloc(e3x1, nmb, ncells3, ncells2, ncells1);
      Kokkos::realloc(e2x1, nmb, ncells3, ncells2, ncells1);
      Kokkos::realloc(e1x2, nmb, ncells3, ncells2, ncells1);
      Kokkos::realloc(e3x2, nmb, ncells3, ncells2, ncells1);
      Kokkos::realloc(e2x3, nmb, ncells3, ncells2, ncells1);
      Kokkos::realloc(e1x3, nmb, ncells3, ncells2, ncells1);
      Kokkos::realloc(e1_cc, nmb, ncells3, ncells2, ncells1);
      Kokkos::realloc(e2_cc, nmb, ncells3, ncells2, ncells1);
      Kokkos::realloc(e3_cc, nmb, ncells3, ncells2, ncells1);
      if (emf_uct || emf_uct_hlld) {   // UCT face store (UctFaceStore); rsolver now known
        if (!(rsolver_method == MHD_RSolver::hlld || rsolver_method == MHD_RSolver::lhlld ||
              (emf_uct && rsolver_method == MHD_RSolver::hlle))) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<mhd>/emf_scheme = uct_hlld needs rsolver = hlld or "
                    << "lhlld; uct needs hlle, hlld or lhlld" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        Kokkos::realloc(uct1, nmb, 6, ncells3, ncells2, ncells1);
        Kokkos::realloc(uct2, nmb, 6, ncells3, ncells2, ncells1);
        if (ppack->pmesh->three_d) Kokkos::realloc(uct3, nmb, 6, ncells3, ncells2, ncells1);
      }
      // prod1103 default (USER 11-03): cs_seam_mhd = shared_o on a FRESH start
      // (recipe::NewDefaults) of cubed-sphere MHD when the key is not named and the
      // shared seam modes are supported (CsSeamInit's conditions: 3-D, reconstruct = plm,
      // no x2 well-balance / static-perturbation reconstruction, a uniform mesh,
      // Newtonian, no fofc, and emf_scheme = uct | uct_hlld); else none, logged.  The
      // value is recorded (GetOrAdd); a restart whose input lacks the key keeps none
      // (bitwise); an explicit key wins.
      if (pmy_pack->pmesh->use_cubed_sphere &&
          !pin->DoesParameterExist("mhd","cs_seam_mhd") &&
          !pin->DoesParameterExist("mhd","cs_seam_emf_avg")) {
        std::string no;
        if (!pmy_pack->pmesh->three_d || recon_method != ReconstructionMethod::plm ||
            (use_wellbalance_dynamic && use_wb_x2) ||
            use_wellbalance_static_reconst_perturb || pmy_pack->pmesh->multilevel ||
            pmy_pack->pcoord->is_special_relativistic ||
            pmy_pack->pcoord->is_general_relativistic || use_fofc) {
          no = "needs 3-D, plm, no x2 well-balance, uniform mesh, Newtonian, no fofc";
        } else if (!(emf_uct || emf_uct_hlld)) {
          no = "needs emf_scheme = uct or uct_hlld";
        }
        if (!recipe::NewDefaults(pin)) {
          recipe::Log("<mhd>/cs_seam_mhd = none (" + std::string(
              recipe::GetMode(pin) == recipe::Mode::none ? "recipe none" :
              "restart input has no key") + ": the old default)");
        } else if (!no.empty()) {
          recipe::Log("<mhd>/cs_seam_mhd = none (shared_o default not applied: " + no +
                      ")");
        } else {
          pin->GetOrAddString("mhd","cs_seam_mhd","shared_o");
          recipe::Log("<mhd>/cs_seam_mhd = shared_o (default on the cubed sphere)");
        }
      }
      // ucthshare-1102 <mhd>/cs_seam_mhd (read only when named; mhd_cs_seam.cpp)
      if (pin->DoesParameterExist("mhd","cs_seam_mhd") ||
          pin->DoesParameterExist("mhd","cs_seam_emf_avg")) {
        CsSeamInit(pin);
      }

      // boris-1102 <mhd>/boris_c (mhd_boris.cpp)
      if (peos->eos_data.boris_c2inv > 0.0) BorisInit(pin);

      // allocate array of flags used with FOFC
      if (use_fofc) {
        int nvars = (pmy_pack->pcoord->is_dynamical_relativistic) ? nmhd+nscalars : nmhd;
        Kokkos::realloc(fofc,    nmb, ncells3, ncells2, ncells1);
        Kokkos::realloc(utest,   nmb, nvars, ncells3, ncells2, ncells1);
        Kokkos::realloc(bcctest, nmb, 3,    ncells3, ncells2, ncells1);
        Kokkos::deep_copy(fofc, false);
        // per-cell FOFC flag count for the `mhd_fofc` output variable (mhd.hpp)
        Kokkos::realloc(fofc_cnt, nmb, ncells3, ncells2, ncells1);
        Kokkos::deep_copy(fofc_cnt, 0.0);
      }
    }
  }
}

//----------------------------------------------------------------------------------------
// destructor

MHD::~MHD() {
  if (sbf_on && global_variable::my_rank == 0) {
    std::cout << "<mhd>/split_bval_fast: energy-only refreshes " << sbf_nfast
              << ", full refreshes " << sbf_nfull << " (rank 0)" << std::endl;
  }
  if (pam != nullptr) {delete pam;}
  delete peos;
  delete pbval_u;
  delete pbval_b;
  if (pbval_fs != nullptr) {delete pbval_fs;}
  if (pbval_us != nullptr) {delete pbval_us;}
  if (pseamfs != nullptr) {delete pseamfs;}
  if (pusown != nullptr) {delete pusown;}
  if (pusnbr != nullptr) {delete pusnbr;}
  CsSeamRelease();
  if (pvisc != nullptr) {delete pvisc;}
  if (presist!= nullptr) {delete presist;}
  if (pcond != nullptr) {delete pcond;}
  if (psrc!= nullptr) {delete psrc;}
}

//----------------------------------------------------------------------------------------
// SetSaveWBcc:  set flag to save primitives and cell-centered B field, e.g., for jcon

void MHD::SetSaveWBcc() {
  int nmb = std::max((pmy_pack->nmb_thispack), (pmy_pack->pmesh->nmb_maxperrank));
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;

  // allocated saved arrays for time derivatives
  Kokkos::realloc(wsaved,   nmb, (nmhd+nscalars), ncells3, ncells2, ncells1);
  Kokkos::realloc(bccsaved, nmb, 3,               ncells3, ncells2, ncells1);

  wbcc_saved = true;
}

} // namespace mhd
