//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hydro.cpp
//! \brief implementation of Hydro class constructor and assorted other functions

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "diffusion/viscosity.hpp"
#include "diffusion/conduction.hpp"
#include "srcterms/srcterms.hpp"
#include "shearing_box/shearing_box.hpp"
#include "shearing_box/orbital_advection.hpp"
#include "bvals/bvals.hpp"
#include "hydro/hydro.hpp"
#include "hydro/hydro_am.hpp"
#include "utils/recipe.hpp"
#include "hydro/hydro_mergelat.hpp"
#include "hydro/hydro_ringavg.hpp"

namespace hydro {
//----------------------------------------------------------------------------------------
// constructor, initializes data structures and parameters

Hydro::Hydro(MeshBlockPack *ppack, ParameterInput *pin) :
    pmy_pack(ppack),
    u0("cons",1,1,1,1,1),
    w0("prim",1,1,1,1,1),
    coarse_u0("ccons",1,1,1,1,1),
    coarse_w0("cprim",1,1,1,1,1),
    u1("cons1",1,1,1,1,1),
    uflx("uflx",1,1,1,1,1),
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
    fofc("fofc",1,1,1,1),
    fofc_cnt("fofc_cnt",1,1,1,1) {
  // Total number of MeshBlocks on this rank to be used in array dimensioning
  int nmb = std::max((ppack->nmb_thispack), (ppack->pmesh->nmb_maxperrank));

  // (1) construct EOS object (no default)
  std::string eqn_of_state = pin->GetString("hydro","eos");
  // ideal gas EOS
  if (eqn_of_state.compare("ideal") == 0) {
    if (pmy_pack->pcoord->is_special_relativistic) {
      peos = new IdealSRHydro(ppack, pin);
    } else if (pmy_pack->pcoord->is_general_relativistic) {
      peos = new IdealGRHydro(ppack, pin);
    } else {
      peos = new IdealHydro(ppack, pin);
    }
    nhydro = 5;
  // general EOS (non-relativistic only)
  } else if (eqn_of_state.compare("general") == 0) {
    if (pmy_pack->pcoord->is_special_relativistic ||
        pmy_pack->pcoord->is_general_relativistic) {
      std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__ << std::endl
                << "<hydro>/eos = general cannot be used with SR/GR" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    peos = new GeneralHydro(ppack, pin);
    nhydro = 5;
  // isothermal EOS
  } else if (eqn_of_state.compare("isothermal") == 0) {
    if (pmy_pack->pcoord->is_special_relativistic ||
        pmy_pack->pcoord->is_general_relativistic) {
      std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__ << std::endl
                << "<hydro>/eos = isothermal cannot be used with SR/GR" << std::endl;
      std::exit(EXIT_FAILURE);
    } else {
      peos = new IsothermalHydro(ppack, pin);
      nhydro = 4;
    }
  // EOS string not recognized
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<hydro>/eos = '" << eqn_of_state << "' not implemented" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // (2) Initialize scalars, diffusion, source terms
  nscalars = pin->GetOrAddInteger("hydro","nscalars",0);
  // share-1104 P6: the dt report keys (hydro_newdt.cpp; env ATHENA_DT_REPORT[_NPOL]
  // aliases), read only when named
  {
    extern int dt_report_key, dt_report_npol_key;
    if (pin->DoesParameterExist("hydro","dt_report")) {
      dt_report_key = pin->GetInteger("hydro","dt_report");
    }
    if (pin->DoesParameterExist("hydro","dt_report_npol")) {
      dt_report_npol_key = pin->GetInteger("hydro","dt_report_npol");
    }
  }
  // share-1104 P4: the passive-tracer layout shared with MHD (mhd/mhd_tracers.hpp): the
  // floor-mass tracer <hydro>/floor_tracer (read only when named; default off) is the
  // LAST scalar slot.  Off, nscalars and every path are unchanged.  Hydro has no x_e
  // tracer (its consumer, the resistivity, is MHD only).
  if (mhd::NamedTrue(pin, "hydro", "xe_tracer")) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<hydro>/xe_tracer is not implemented (MHD only: its "
              << "consumer is the resistivity); use <mhd>/xe_tracer" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  tracers = mhd::ReadTracerLayout(pin, "hydro");
  tracers.nuser = nscalars;
  nscalars = tracers.Total();
  floor_tracer = tracers.fl;
  ifltr = tracers.fl ? (nhydro + tracers.FlOff()) : -1;
  // duale-1112 (hydro/dual_energy.hpp): every key read only when named
  dual.Read(pin, "hydro", peos->eos_data, nhydro, tracers.DuOff());
  if (floor_tracer) {
    if (nhydro != 5 || pmy_pack->pcoord->is_special_relativistic ||
        pmy_pack->pcoord->is_general_relativistic) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<hydro>/floor_tracer needs a non-relativistic hydro "
                << "with an energy equation (eos ideal or general)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    peos->eos_data.ifltr = ifltr;
    if (pin->DoesParameterExist("hydro","floor_tracer_margin")) {
      fltr_margin = pin->GetReal("hydro","floor_tracer_margin");
    }
    if (pin->DoesParameterExist("hydro","floor_tracer_fmax")) {
      fltr_fmax = pin->GetReal("hydro","floor_tracer_fmax");
    }
    if (fltr_margin < 1.0 || fltr_fmax <= 0.0 || fltr_fmax > 1.0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<hydro>/floor_tracer needs floor_tracer_margin >= 1 and "
                << "0 < floor_tracer_fmax <= 1" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (global_variable::my_rank == 0) {
      std::cout << "<hydro>/floor_tracer: floor-added mass in scalar slot " << ifltr
                << " (nscalars " << nscalars << "); drained once per cycle where rho > "
                << fltr_margin << " x dfloor, at most " << fltr_fmax
                << " of the cell mass per cycle" << std::endl;
    }
  }

  // Viscosity (if requested in input file)
  if (pin->DoesParameterExist("hydro","isotropic_viscosity")) {
    pvisc = new Viscosity("hydro", ppack, pin);
  } else {
    pvisc = nullptr;
  }

  // Thermal conduction (if requested in input file)
  // isotropic_conduction = none: no conduction module at all, exactly as if the key
  // were absent (lets a restart drop it: ck-nq2, the ck two-stream down to the wall)
  if (pin->DoesParameterExist("hydro","isotropic_conduction") &&
      pin->GetString("hydro","isotropic_conduction").compare("none") != 0) {
    if (peos->eos_data.is_ideal) {
      pcond = new Conduction("hydro", ppack, pin);
    } else {
      std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__ << std::endl
                << "Thermal conduction in hydro requires ideal gas EOS" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else {
    pcond = nullptr;
  }

  // Source terms (if needed)
  if (pin->DoesBlockExist("hydro_srcterms")) {
    psrc = new SourceTerms("hydro_srcterms", ppack, pin);
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
    Kokkos::realloc(u0, nmb, (nhydro+nscalars), ncells3, ncells2, ncells1);
    Kokkos::realloc(w0, nmb, (nhydro+nscalars), ncells3, ncells2, ncells1);
    // derived thermodynamic variables are only needed for an expensive (general) EOS
    if (peos->eos_data.IsGeneral()) {
      Kokkos::realloc(wder, nmb, NDERIVED, ncells3, ncells2, ncells1);
      Kokkos::realloc(wtemp, nmb, ncells3, ncells2, ncells1);
    }
    // <hydro>/dfloor_keep_velocity with deferred conserved floors (the cubed sphere):
    // the metric-correct kinetic-energy correction is made by GnomonicEquiangleRaiseVel,
    // which reads fv from here.  Needed for both EOS flavours, so it sits outside the
    // general-EOS branch above.
    if ((peos->eos_data.dfloor_keep_velocity ||
         peos->eos_data.dfloor_keep_temperature) &&
        peos->eos_data.defer_cons_floors) {
      Kokkos::realloc(dfl_fv, nmb, ncells3, ncells2, ncells1);
      Kokkos::deep_copy(dfl_fv, 1.0);
    }
  }

  // allocate memory for conserved variables on coarse mesh
  if (ppack->pmesh->multilevel) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    int n_ccells1 = indcs.cnx1 + 2*(indcs.ng);
    int n_ccells2 = (indcs.cnx2 > 1)? (indcs.cnx2 + 2*(indcs.ng)) : 1;
    int n_ccells3 = (indcs.cnx3 > 1)? (indcs.cnx3 + 2*(indcs.ng)) : 1;
    Kokkos::realloc(coarse_u0, nmb, (nhydro+nscalars), n_ccells3, n_ccells2, n_ccells1);
    Kokkos::realloc(coarse_w0, nmb, (nhydro+nscalars), n_ccells3, n_ccells2, n_ccells1);
  }

  // allocate boundary buffers for conserved (cell-centered) variables
  pbval_u = new MeshBoundaryValuesCC(ppack, pin, false);
  pbval_u->cs_rho_guard = ppack->pmesh->cs_seam_rho_guard;
  pbval_u->cs_lin_resample = ppack->pmesh->cs_seam_resample_linear;
  pbval_u->cs_blockend_ext_ = CsSeamBlockendExtend(ppack, pin);   // layoutbug-1103
  pbval_u->InitializeBuffers((nhydro+nscalars));
  // csedge-1019 <mesh>/cs_seam_flux = shared: the seam face-state exchange (hydro.hpp)
  if (ppack->pmesh->use_cubed_sphere && ppack->pmesh->cs_seam_flux == 3) {
    cs_seam_shared = true;
    const int nfs = nhydro + nscalars + NDERIVED;
    pbval_fs = new MeshBoundaryValuesCC(ppack, pin, false);
    pbval_fs->seam_raw = true;
    pbval_fs->InitializeBuffers(nfs);
    auto &ind = ppack->pmesh->mb_indcs;
    const int n1 = ind.nx1 + 2*ind.ng, n2 = ind.nx2 + 2*ind.ng, n3 = ind.nx3 + 2*ind.ng;
    pseamfs = new DvceFaceFld5D<Real>("cs_seamfs", nmb, nfs, n3, n2, n1);
    Kokkos::realloc(pseamfs->x1f, nmb, nfs, 1, 1, 1);   // x1 faces never cross a seam
  }

  // Orbital advection and shearing box BCs (if requested in input file)
  if (pin->DoesBlockExist("shearing_box")) {
    porb_u = new OrbitalAdvectionCC(ppack, pin, (nhydro+nscalars));
    psbox_u = new ShearingBoxCC(ppack, pin, (nhydro+nscalars));
  } else {
    porb_u = nullptr;
    psbox_u = nullptr;
  }
        
  // determine if etotgrav and local well-balanced scheme is enabled
  // amdef-1019: default true where supported (hydro::etotgrav_pgen_default, an energy
  // equation, no refinement) under the new default set (utils/recipe.hpp)
  if (!pin->DoesParameterExist("hydro","etotgrav") && etotgrav_pgen_default &&
      pin->GetString("hydro","eos") != "isothermal" &&
      !pmy_pack->pmesh->multilevel) {
    if (recipe::NewDefaults(pin)) {
      pin->GetOrAddBoolean("hydro","etotgrav",true);
      recipe::Log("<hydro>/etotgrav = true (default where supported)");
    } else {
      recipe::Log("<hydro>/etotgrav = false (" + std::string(
          recipe::GetMode(pin) == recipe::Mode::none ? "recipe none" :
          "restart input has no key") + ": the old default)");
    }
  }
  // wbplm-1103 (read only when named): see hydro.hpp top_ghost_prim_copy
  if (pin->DoesParameterExist("hydro","top_ghost_prim_copy")) {
    top_ghost_prim_copy = pin->GetBoolean("hydro","top_ghost_prim_copy");
  }
  use_etotgrav = pin->GetOrAddBoolean("hydro","etotgrav",false);
  // ckbval-1102: the energy-only ghost refresh after the operator-split source (default
  // ON after the bitwise gates; <hydro>/split_bval_fast = false restores the full refresh; read
  // only when named; see Hydro::SplitBvalFastOK).  Not with refinement, shearing box /
  // orbital advection or the polar grid: their halos are not plain copies or seam
  // resamples, which is all the check below covers.
  sbf_on = pin->DoesParameterExist("hydro","split_bval_fast") ?
           pin->GetBoolean("hydro","split_bval_fast") : true;
  // prod1103 default (USER 11-02/11-03): split_c2p_once = true on a FRESH start
  // (recipe::NewDefaults; ckbval_1102: gates a/c/d bitwise where they must be, seed-noise
  // gate b: the change is at the round-off seed-spread level, worst bin 1.26x the 1e-14
  // and 1.61x the 1e-12 member spread; -11 % hydro ms/cycle).  The value is recorded
  // (GetOrAdd); a restart whose input lacks the key keeps false (bitwise); an explicit key
  // wins.  With it on, c2p_idempotent is a no-op (bitwise, ckbval_1102 gates_idem).
  // prod-1111 (USER 11-08): also true on a restart whose input lacks the key (recorded);
  // only <job>/recipe = none keeps the old false.  Name false for prod-1110.
  if (!pin->DoesParameterExist("hydro","split_c2p_once")) {
    if (recipe::GetMode(pin) != recipe::Mode::none) {
      pin->GetOrAddBoolean("hydro","split_c2p_once",true);
      recipe::Log("<hydro>/split_c2p_once = true (default, fresh start and restart)");
    } else {
      recipe::Log("<hydro>/split_c2p_once = false (recipe none: the old default)");
    }
  }
  scp_once = pin->DoesParameterExist("hydro","split_c2p_once") ?
             pin->GetBoolean("hydro","split_c2p_once") : false;
  // prod1104b default (USER 11-04) WITHDRAWN by defaults-1105 (USER 11-05: GPU timing showed
  // no gain, RSG -2 %, dhj +2 % = noise; it never acts in dhj hydro): c2p_idempotent is again
  // read only when named, default false on a fresh start.  A restart whose input names it
  // (runs started 11-04/05 recorded true) keeps that value; an explicit key wins.
  c2pi_on = pin->DoesParameterExist("hydro","c2p_idempotent") ?
            pin->GetBoolean("hydro","c2p_idempotent") : false;
  if (sbf_on && (pmy_pack->pmesh->multilevel ||
                 pin->DoesBlockExist("shearing_box") ||
                 pin->DoesBlockExist("orbital_advection"))) {
    sbf_on = false;
  }
  // SMR/AMR is not supported with etotgrav: the coarse-buffer inversion used for
  // primitive prolongation (bvals/prolong_prims.cpp) reads e + rho*Phi as the internal
  // energy, the conserved prolongation of rho*Phi is removed again with the FINE
  // potential, and a regrid never rebuilds phicc0/phi0 on the new blocks.
  if (use_etotgrav && pmy_pack->pmesh->multilevel) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<hydro>/etotgrav = true is not supported with mesh "
              << "refinement (SMR/AMR): the potential energy is not handled at level "
              << "boundaries or on regrid." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  use_wellbalance_dynamic = pin->GetOrAddBoolean("hydro","wellbalance_dynamic",false);
  use_wb_x1 = pin->GetOrAddBoolean("hydro","wb_x1",false);
  // csedge-1019 <hydro>/cs_seam_rsolver (cubed sphere only): "same" = the seam faces use
  // <hydro>/rsolver like every other face; "hlle" = they take the HLLE flux.  HLLC (and
  // its variants) resolve the contact, so a shear ACROSS a seam face is not damped at all
  // at S_M ~ 0; with the resampled seam ghosts that left a growing one-cell along-edge jet
  // in the edge row, seeded at the cube vertices (RSG g25 cs, csedge_1019/NOTE.md).
  {
    const std::string csr = pin->GetOrAddString("hydro","cs_seam_rsolver","same");
    if (csr != "same" && csr != "hlle") {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<hydro>/cs_seam_rsolver must be same or hlle" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    cs_seam_hlle = (csr == "hlle") && ppack->pmesh->use_cubed_sphere;
  }
  // rebuild the well-balanced background cache at stage 1 of every N-th cycle instead
  // of every stage (0 = every stage).  A balanced state is unaffected; an evolving one
  // reconstructs against a background up to N cycles old, which is nothing against the
  // radiative timescale the profile changes on.  N = 1 halves the walk cost, N = 10
  // removes it.
  wb_cache_every = pin->GetOrAddInteger("hydro","wb_cache_every",0);
  // Ideal-gas x1 WB cache (BuildWBCache) safety keys, BSG arm 1 crash (10-02):
  //  <hydro>/wb_clamp_ideal: polytropic closed form limits the T DECREASE over a
  //    segment to t0/2, as the general-EOS walk does (getWBerho, clamp_dt);
  //  <hydro>/wb_guard_ideal: WBGuard flattens a non-finite, non-positive or absurd
  //    walked stencil to its anchor (backstop).
  // Both DEFAULT true on a fresh run where that path is active (wb_x1 with a
  // non-general EOS); a restart whose input lacks a key keeps the old default (off,
  // one line); an explicit value always wins.
  const bool wbi_path = use_wb_x1 && !peos->eos_data.IsGeneral();
  auto wbi_default = [&](const char *key) {
    if (!wbi_path || pin->DoesParameterExist("hydro", key)) return false;
    if (!global_variable::restart_run) return true;
    if (global_variable::my_rank == 0) {
      std::cout << "hydro: restart input has no hydro/" << key << ": keeping the old "
                << "default (off); set it true to switch" << std::endl;
    }
    return false;
  };
  wb_guard_ideal = pin->GetOrAddBoolean("hydro","wb_guard_ideal",
                                        wbi_default("wb_guard_ideal"));
  wb_clamp_ideal = pin->GetOrAddBoolean("hydro","wb_clamp_ideal",
                                        wbi_default("wb_clamp_ideal"));
  scratch_level = pin->GetOrAddInteger("hydro","scratch_level",0);
  use_wb_x2 = pin->GetOrAddBoolean("hydro","wb_x2",false);
  use_wb_rho = pin->GetOrAddBoolean("hydro","wb_rho",false);
  // wallwb-1105 (hydro.hpp, reconstruct/wb_wall.hpp).  USER 11-05: default true on a
  // FRESH start (recipe::NewDefaults) whenever the x1 dynamic WB scheme is on and the
  // key is not named (it acts only at closed / reflecting x1 walls); recorded, so the
  // run's restarts carry it.  A restart whose input lacks the key keeps the old
  // behaviour (false, bitwise); an explicit key always wins.
  if (!pin->DoesParameterExist("hydro","wall_wb_onesided") &&
      use_wellbalance_dynamic && use_wb_x1) {
    // prod-1111 (USER 11-08): also on a restart whose input lacks the key (recorded);
    // only <job>/recipe = none keeps the old wall cell.  Name false for prod-1110.
    if (recipe::GetMode(pin) != recipe::Mode::none) {
      pin->GetOrAddBoolean("hydro","wall_wb_onesided",true);
      recipe::Log("<hydro>/wall_wb_onesided = true (default with wellbalance_dynamic + "
                  "wb_x1, fresh start and restart; acts at closed / reflecting x1 "
                  "walls)");
    } else {
      recipe::Log("<hydro>/wall_wb_onesided = false (recipe none: the old wall cell)");
    }
  }
  if (pin->DoesParameterExist("hydro","wall_wb_onesided")) {
    wall_wb_onesided = pin->GetBoolean("hydro","wall_wb_onesided");
  }
  // switch the x1 well-balanced reconstruction off above this radius (0 = never); see
  // the declaration in hydro.hpp.  Only the position-aware x1 path (spherical polar and
  // cubed sphere, where x1v is the radius) honours it -- a Cartesian mesh has no radius
  // to compare against.
  wb_rmax = pin->GetOrAddReal("hydro","wb_rmax",0.0);
  if (pin->DoesParameterExist("hydro","wb_poison_test")) wb_poison = pin->GetBoolean("hydro","wb_poison_test");
  // ...and off BELOW this radius (0 = never), same path, same caveat.
  wb_rmin = pin->GetOrAddReal("hydro","wb_rmin",0.0);
  nan_report = pin->GetOrAddBoolean("problem","nan_report",false);
  // amdef-1019: am_transport defaults on (cs, multi-D sp, fresh dynamic runs): AmModule
  {
    // (DefaultOn first: it is false at once when am_transport itself is named)
    const bool amd = hydro::AmModule::DefaultOn(ppack, pin, "hydro");
    if (amd || hydro::AmModule::Named(pin, "hydro")) {
      pam = new hydro::AmModule(ppack, pin, "hydro");
      pam->am_default = amd;
    }
  }
  if (hydro::MergeLat::Named(pin)) {
    pmlat = new hydro::MergeLat(ppack, pin);
  }
  // layoutaudit-1103 <hydro>/merge_lat_blockend = flat (the old rule: a 'narrow' merged
  // group at EVERY MeshBlock edge takes a zero slope, so the result depends on the
  // MeshBlock layout) | exchange (internal block faces exchange the one cell the group
  // PLM needs; layout-independent).  Absent: exchange on a fresh start under the default
  // recipe (recorded, so its restarts carry it), flat on a restart (bitwise the old run).
  if (pmlat != nullptr && pmlat->on) {
    std::string mb;
    if (pin->DoesParameterExist("hydro", "merge_lat_blockend")) {
      mb = pin->GetString("hydro", "merge_lat_blockend");
    } else if (recipe::NewDefaults(pin)) {
      mb = pin->GetOrAddString("hydro", "merge_lat_blockend", "exchange");
      recipe::Log("<hydro>/merge_lat_blockend = exchange (default on a fresh start)");
    } else {
      mb = "flat";
    }
    if (mb == "exchange") {
      mlat_exch = true;
    } else if (mb != "flat") {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<hydro>/merge_lat_blockend = flat | exchange" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    auto &ind = ppack->pmesh->mb_indcs;
    if (mlat_exch && (ind.nx2 < 2*pmlat->fbig || ind.nx3 < 2*pmlat->fbig)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<hydro>/merge_lat_blockend = exchange needs MeshBlock "
                << "nx2, nx3 >= 2 x the largest merge_lat_f (" << pmlat->fbig << ")"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (mlat_exch) {
      const int nfs = nhydro + nscalars + NDERIVED;
      pbval_ml = new MeshBoundaryValuesCC(ppack, pin, false);
      pbval_ml->seam_raw = true;
      pbval_ml->seam_same_panel = true;
      pbval_ml->InitializeBuffers(nfs);
      if (pseamfs != nullptr) {
        pmlfs = pseamfs;   // seam faces and internal faces never coincide: one record
      } else {
        const int n1 = ind.nx1 + 2*ind.ng, n2 = ind.nx2 + 2*ind.ng;
        const int n3 = ind.nx3 + 2*ind.ng;
        pmlfs = new DvceFaceFld5D<Real>("mlat_blockend", nmb, nfs, n3, n2, n1);
        Kokkos::realloc(pmlfs->x1f, nmb, nfs, 1, 1, 1);
        mlat_own_rec = true;
      }
      // internal-face mask from the (static) neighbour table: same panel, same level,
      // not polar; buffers 8 / 12 = x2 inner / outer face, 24 / 28 = x3
      auto &nghbr = ppack->pmb->nghbr;
      auto &mbpanel = ppack->pmb->mb_panel;
      auto &mblev = ppack->pmb->mb_lev;
      Kokkos::realloc(mlat_imask, nmb);
      auto hm = Kokkos::create_mirror_view(mlat_imask);
      const int nb_[4] = {8, 12, 24, 28};
      for (int m = 0; m < nmb; ++m) hm(m) = 0;
      for (int m = 0; m < ppack->nmb_thispack; ++m) {
        int msk = 0;
        for (int q = 0; q < 4; ++q) {
          const auto nb = nghbr.h_view(m, nb_[q]);
          if (nb.gid >= 0 && nb.panel == mbpanel.h_view(m) && nb.lev == mblev.h_view(m) &&
              nb.polar <= 0) {
            msk |= (1 << q);
          }
        }
        hm(m) = msk;
      }
      Kokkos::deep_copy(mlat_imask, hm);
    }
  }
  if (pin->GetOrAddBoolean("mesh", "polar_ring_avg", false)) {
    if (pmlat != nullptr && pmlat->on) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<mesh>/polar_ring_avg and <hydro>/merge_lat_rc cannot "
                << "be combined" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    pravg = new hydro::RingAvg(ppack, pin);
  }
  // radiation-modified CFL signal speed (hydro_newdt.cpp).  DEFAULT: on with an
  // IMPLICIT M1 transport (<rad_m1>/transport = implicit | implicit_x1, any closure),
  // whose time step has no light-speed limit; off without <rad_m1> and with the explicit
  // transport (dt already set by the light speed).  Read only when named, so the
  // parameter dump of a run without the key is unchanged.
  if (pin->DoesBlockExist("rad_m1") &&
      pin->DoesParameterExist("rad_m1","transport")) {
    std::string tr = pin->GetString("rad_m1","transport");
    rad_signal_speed = (tr.compare("implicit") == 0 || tr.compare("implicit_x1") == 0);
  }
  if (pin->DoesParameterExist("hydro","rad_signal_speed")) {
    rad_signal_speed = pin->GetBoolean("hydro","rad_signal_speed");
  }
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
  // the x1 well-balanced scheme sees the TRUE potential unless a problem generator calls
  // EnableWBEffectivePotential(): these are shallow copies, the very same allocations
  if (!use_phi_wb) {
    phicc_wb = phicc0;
    phi_wb_x1f = phi0.x1f;
  }
  if (use_wellbalance_dynamic) {
    // select well-balanced scheme assumption (no default).  Test for compatibility of options
    std::string wb_opt = pin->GetString("hydro","wb_option");
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
                << std::endl << "<hydro> wb_option = '" << wb_opt
                << "' not implemented for hydro" << std::endl;
      std::exit(EXIT_FAILURE);
    }
      // wball-1104 <hydro>/wb_dirs = x1 (default) | all (read only when named)
      if (pin->DoesParameterExist("hydro","wb_dirs")) {
        const std::string wd = pin->GetString("hydro","wb_dirs");
        if (wd.compare("all") == 0) {
          wb_all = true;
        } else if (wd.compare("x1") != 0) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<hydro>/wb_dirs must be x1 or all" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
      if (wb_all) {
        std::string bad = "";
        if (!use_wb_x1) bad = "<hydro>/wb_x1 = true";
        if (use_wb_x2) bad = "<hydro>/wb_x2 = false (the legacy x2 path)";
        if (pmy_pack->pcoord->is_general_relativistic ||
            pmy_pack->pcoord->is_special_relativistic) {
          bad = "a Newtonian mesh";
        }
        if (pmy_pack->pmesh->one_d) bad = "a 2-D or 3-D mesh";
        if (!bad.empty()) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<hydro>/wb_dirs = all needs " << bad << std::endl;
          std::exit(EXIT_FAILURE);
        }
        wb_weight = pin->GetOrAddBoolean("hydro","wb_weight",true);
        wbw_r0 = pin->GetOrAddReal("hydro","wb_w_r0",0.1);
        wbw_r1 = pin->GetOrAddReal("hydro","wb_w_r1",0.3);
        wbw_m0 = pin->GetOrAddReal("hydro","wb_w_m0",0.2);
        wbw_m1 = pin->GetOrAddReal("hydro","wb_w_m1",0.5);
        wbw_eta = pin->GetOrAddReal("hydro","wb_w_eta",1.0e-4);
        auto &ic = pmy_pack->pmesh->mb_indcs;
        const int c1 = ic.nx1 + 2*ic.ng;
        const int c2 = (ic.nx2 > 1) ? (ic.nx2 + 2*ic.ng) : 1;
        const int c3 = (ic.nx3 > 1) ? (ic.nx3 + 2*ic.ng) : 1;
        Kokkos::realloc(wbq2, nmb, 15, c3, c2, c1);
        if (c3 > 1) Kokkos::realloc(wbq3, nmb, 15, c3, c2, c1);
        Kokkos::realloc(wbw, nmb, 7, c3, c2, c1);
        Kokkos::deep_copy(wbw, 1.0);
        // wball-1104 <hydro>/wb_seam (opt-in, read only when named): the WB-aware seam
        // resample of the density and energy (bvals.hpp, MeshBoundaryValuesCC::wb_seam)
        if (pin->DoesParameterExist("hydro","wb_seam") && pin->GetBoolean("hydro","wb_seam")) {
          if (!pmy_pack->pmesh->use_cubed_sphere || pmy_pack->pmesh->multilevel) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "<hydro>/wb_seam needs a single-level cubed sphere"
                      << std::endl;
            std::exit(EXIT_FAILURE);
          }
          pbval_u->wb_seam = true;
          pbval_u->wbs_phi = phicc0;
          pbval_u->wbs_eos = peos->eos_data;
          pbval_u->wbs_opt = (wb_option == WBOption::polytropic) ? 3 : 1;
          pbval_u->wbs_etg = use_etotgrav;
          pbval_u->wbs_ien = IEN;
          if (global_variable::my_rank == 0) {
            std::cout << "hydro: wb_seam = true (seam resample of rho, E walked to the ghost "
                      << "potential, " << ((pbval_u->wbs_opt == 3) ? "local polytrope" :
                      "isothermal") << ")" << std::endl;
          }
        }
        if (global_variable::my_rank == 0) {
          std::cout << "hydro: wb_dirs = all (x2/x3 dynamic well balance), wb_weight = "
                    << (wb_weight ? "true" : "false") << std::endl;
        }
      }
      // select well-balanced scheme direction
      if (!use_wb_x1 && !use_wb_x2) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<hydro> wb_direction not set!" << std::endl;
        std::exit(EXIT_FAILURE);
      }
  }
        
  // determine if static wellbalance is enabled
  use_wellbalance_static = pin->GetOrAddBoolean("hydro","wellbalance_static",false);
  use_wellbalance_static_reconst_perturb = pin->GetOrAddBoolean("hydro","wellbalance_static_reconst",false);
  // the deviation reconstruction reads the background arrays, which only
  // wellbalance_static allocates -- the pair is one scheme, not two options
  if (use_wellbalance_static_reconst_perturb && !use_wellbalance_static) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<hydro>/wellbalance_static_reconst needs "
              << "<hydro>/wellbalance_static = true" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // allocate array of flags used with wellbalance
  if (use_wellbalance_static) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    int ncells1 = indcs.nx1 + 2*(indcs.ng);
    int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
    int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
    Kokkos::realloc(u0wb, nmb, nhydro, ncells3, ncells2, ncells1);
    Kokkos::realloc(w0wb, nmb, nhydro, ncells3, ncells2, ncells1);
    Kokkos::realloc(w0facewb.x1f, nmb, nhydro, ncells3, ncells2, ncells1+1);
    Kokkos::realloc(w0facewb.x2f, nmb, nhydro, ncells3, ncells2+1, ncells1);
    Kokkos::realloc(w0facewb.x3f, nmb, nhydro, ncells3+1, ncells2, ncells1);
    // background pressure: one scalar per cell and per face, filled once (see
    // SetWbBackgroundPressure) instead of being re-derived from the background every
    // stage
    Kokkos::realloc(pwb, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(pfacewb.x1f, nmb, ncells3, ncells2, ncells1+1);
    Kokkos::realloc(pfacewb.x2f, nmb, ncells3, ncells2+1, ncells1);
    Kokkos::realloc(pfacewb.x3f, nmb, ncells3+1, ncells2, ncells1);
  }

  // for time-evolving problems, continue to construct methods, allocate arrays
  if (evolution_t.compare("stationary") != 0) {
    // determine if FOFC is enabled
    // default false; true for the problem generators that set fofc_pgen_default
    // (hydro.hpp), except on a restart whose file lacks the key (old default kept)
    bool fofc_dflt = false;
    if (fofc_pgen_default && !pin->DoesParameterExist("hydro","fofc")) {
      if (global_variable::restart_run) {
        if (global_variable::my_rank == 0) {
          std::cout << "hydro: restart input has no hydro/fofc: keeping the old default "
                    << "(off); set it true to switch" << std::endl;
        }
      } else {
        fofc_dflt = true;
      }
    }
    use_fofc = pin->GetOrAddBoolean("hydro","fofc",fofc_dflt);
    if (use_fofc) {fofc_rep.Init(pin, "hydro");}

    // select reconstruction method (default PLM)
    std::string xorder = pin->GetOrAddString("hydro","reconstruct","plm");
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
                << std::endl << "<hydro> reconstruct = '" << xorder << "' not implemented"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // select RADIAL (x1) reconstruction on the cubed sphere / spherical polar (default
    // PLM, i.e. today's GridPiecewiseLinearX1, bitwise unchanged). `reconstruct` above
    // never reaches the x1 sweep on those two grids -- see the note on str_r1_ in
    // hydro_fluxes.cpp -- so this is a SEPARATE knob, not a fallback for it.
    {
      bool cs_or_sp = pmy_pack->pmesh->use_cubed_sphere ||
                      pmy_pack->pmesh->use_spherical_polar;
      std::string xorder1 = pin->GetOrAddString("hydro","reconstruct_x1","plm");
      auto &indcs = pmy_pack->pmesh->mb_indcs;
      if (xorder1.compare("plm") == 0) {
        recon_method_x1 = ReconstructionMethod::plm;
      } else if (xorder1.compare("ppmx") == 0) {
        // See the identical refusal + citation in mhd.cpp: Mignone (2014) [JCP 270,
        // 784, arXiv:1404.0537] gives no curvilinear Colella-Sekora extremum-preserving
        // PPM limiter, only the CW84-style one (his Eqs. 45-48); use ppm4 instead.
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "<hydro>/reconstruct_x1=ppmx is not implemented: Mignone "
          << "(2014) gives no curvilinear Colella-Sekora extremum-preserving PPM "
          << "limiter; use reconstruct_x1=ppm4 (curvilinear CW84 limiter) or plm"
          << std::endl;
        std::exit(EXIT_FAILURE);
      } else if (xorder1.compare("ppm4") == 0 ||
                 xorder1.compare("wenoz") == 0) {
        // The published, non-uniform-grid-EXACT reconstruction from A. Mignone, JCP
        // 270, 784 (2014), arXiv:1404.0537: ppm4 = his curvilinear 4th-order PPM
        // (Sec. 3.3); wenoz = his curvilinear WENO3 (Sec. 3.2, third order -- a genuine
        // order reduction from the uniform-grid 5th-order WENO-Z used elsewhere, since
        // the paper gives no higher-order curvilinear WENO). See
        // reconstruct/mignone_curvilinear.hpp for the implementation.
        if (indcs.ng < 3) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "<hydro>/reconstruct_x1=" << xorder1 << " requires at least "
            << "3 ghost zones, but <mesh>/nghost=" << indcs.ng << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (use_fofc && indcs.ng < 4) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "FOFC and <hydro>/reconstruct_x1=" << xorder1 << " require "
            << "at least 4 ghost zones, but <mesh>/nghost=" << indcs.ng << std::endl;
          std::exit(EXIT_FAILURE);
        }
        // The Kappeli & Mishra dynamic radial well-balancing is PLM-specific (see the
        // identical refusal in mhd.cpp); refuse rather than silently drop it.
        if (use_wellbalance_dynamic && use_wb_x1) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "<hydro>/reconstruct_x1=" << xorder1 << " is not implemented "
            << "with <hydro>/wellbalance_dynamic + wb_x1 (the K&M radial well-balancing "
            << "is PLM-only); use <hydro>/reconstruct_x1=plm or turn wellbalance_dynamic "
            << "off" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (xorder1.compare("ppm4") == 0) {
          recon_method_x1 = ReconstructionMethod::ppm4;
        } else {
          recon_method_x1 = ReconstructionMethod::wenoz;
        }
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<hydro>/reconstruct_x1 = '" << xorder1
                  << "' not implemented" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (cs_or_sp && global_variable::my_rank == 0) {
        std::cout << "Radial (x1) reconstruction on the "
                  << (pmy_pack->pmesh->use_cubed_sphere ? "cubed sphere" : "spherical "
                      "polar") << " grid (<hydro>/reconstruct_x1) = " << xorder1
                  << std::endl;
      }
    }

    // select Riemann solver (no default).  Test for compatibility of options
    std::string rsolver = pin->GetString("hydro","rsolver");
    // Special relativistic dynamic solvers
    if (pmy_pack->pcoord->is_special_relativistic) {
      if (evolution_t.compare("dynamic") == 0) {
        if (rsolver.compare("llf") == 0) {
          rsolver_method = Hydro_RSolver::llf_sr;
        } else if (rsolver.compare("hlle") == 0) {
          rsolver_method = Hydro_RSolver::hlle_sr;
        } else if (rsolver.compare("hllc") == 0) {
          rsolver_method = Hydro_RSolver::hllc_sr;
        // Error for anything else
        } else {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<hydro> rsolver = '" << rsolver
                    << "' not implemented for SR dynamics" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "kinematic dynamics not implemented for SR" <<std::endl;
        std::exit(EXIT_FAILURE);
      }

    // General relativistic dynamic solvers
    } else if (pmy_pack->pcoord->is_general_relativistic) {
      if (evolution_t.compare("dynamic") == 0) {
        if (rsolver.compare("llf") == 0) {
          rsolver_method = Hydro_RSolver::llf_gr;
        } else if (rsolver.compare("hlle") == 0) {
          rsolver_method = Hydro_RSolver::hlle_gr;
        // Error for anything else
        } else {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<hydro> rsolver = '" << rsolver
                    << "' not implemented for GR dynamics" << std::endl;
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
        rsolver_method = Hydro_RSolver::llf;
      // HLLE solver
      } else if (rsolver.compare("hlle") == 0) {
        rsolver_method = Hydro_RSolver::hlle;
      // HLLC solver
      } else if (rsolver.compare("hllc") == 0) {
        if (peos->eos_data.is_ideal) {
          rsolver_method = Hydro_RSolver::hllc;
        } else {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<hydro>/rsolver = hllc cannot be used with "
                    << "isothermal EOS" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      // LHLLC solver
      } else if (rsolver.compare("lhllc") == 0) {
        if (peos->eos_data.is_ideal) {
          rsolver_method = Hydro_RSolver::lhllc;
        } else {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<hydro>/rsolver = lhllc cannot be used with "
                    << "isothermal EOS" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      // HLLCLM solver
      } else if (rsolver.compare("hllclm") == 0) {
        if (peos->eos_data.is_ideal) {
          rsolver_method = Hydro_RSolver::hllclm;
        } else {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<hydro>/rsolver = hllclm cannot be used with "
                    << "isothermal EOS" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      // AUSMPUP solver
      } else if (rsolver.compare("ausmpup") == 0) {
        if (peos->eos_data.is_ideal) {
          rsolver_method = Hydro_RSolver::ausmpup;
          peos->eos_data.ausm_mcut = pin->GetOrAddReal("hydro","ausm_mcut",1.0e-13);
          peos->eos_data.ausm_mcut_p = pin->GetOrAddReal("hydro","ausm_mcut_p",1.0);
          peos->eos_data.ausm_wall_hllc = pin->GetOrAddBoolean("hydro","ausm_wall_hllc",
                                                                 false);
        } else {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<hydro>/rsolver = ausmpup cannot be used with "
                    << "isothermal EOS" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      // Roe solver
      } else if (rsolver.compare("roe") == 0) {
        // The Roe solver linearizes the flux Jacobian about a Roe-averaged state whose
        // construction assumes an ideal gas throughout (see RoeFluxAdb, which is written
        // entirely in terms of gm1). Substituting a local Gamma_1 would NOT give a valid
        // Roe average: the linearized matrix would no longer satisfy the jump condition
        // F_R - F_L = A (U_R - U_L), which is the defining property of the scheme, and
        // shock speeds would be silently wrong. A correct general-EOS Roe solver needs
        // Vinokur-Montagne averaging, which is not implemented. Refuse instead.
        if (peos->eos_data.IsGeneral()) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<hydro>/rsolver = roe is not supported with "
                    << "<hydro>/eos = general: the Roe average is only defined for an "
                    << "ideal gas. Use llf, hlle, hllc, hllclm, lhllc or ausmpup."
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
        rsolver_method = Hydro_RSolver::roe;
      // Error for anything else
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<hydro> rsolver = '" << rsolver << "' not implemented"
                  << " for dynamic problems" << std::endl;
        std::exit(EXIT_FAILURE);
      }

    // Non-relativistic kinematic solvers
    } else {
      // Advect solver
      if (rsolver.compare("advect") == 0) {
        rsolver_method = Hydro_RSolver::advect;
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<hydro> rsolver = '" << rsolver << "' not implemented"
                  << " for kinematic problems" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    // Final memory allocations
    {
      // allocate second registers, fluxes
      auto &indcs = pmy_pack->pmesh->mb_indcs;
      int ncells1 = indcs.nx1 + 2*(indcs.ng);
      int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
      int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
      Kokkos::realloc(u1,       nmb, (nhydro+nscalars), ncells3, ncells2, ncells1);
      Kokkos::realloc(uflx.x1f, nmb, (nhydro+nscalars), ncells3, ncells2, ncells1);
      Kokkos::realloc(uflx.x2f, nmb, (nhydro+nscalars), ncells3, ncells2, ncells1);
      Kokkos::realloc(uflx.x3f, nmb, (nhydro+nscalars), ncells3, ncells2, ncells1);

      // allocate array of flags used with FOFC
      if (use_fofc) {
        Kokkos::realloc(fofc,  nmb, ncells3, ncells2, ncells1);
        // duale-1112: the trial state carries the scalars under <hydro>/dual_energy
        Kokkos::realloc(utest, nmb, nhydro + (dual.on ? nscalars : 0), ncells3, ncells2,
                        ncells1);
        // per-cell FOFC flag count for the `hydro_fofc` output variable (hydro.hpp)
        Kokkos::realloc(fofc_cnt, nmb, ncells3, ncells2, ncells1);
        Kokkos::deep_copy(fofc_cnt, 0.0);
        // the general-EOS single-state LLF fallback in hydro_fofc.cpp indexes wder;
        // fail loudly here rather than read a 1x1x1x1 placeholder on the device
        if (peos->eos_data.IsGeneral() &&
            (wder.extent_int(0) != nmb || wder.extent_int(4) != ncells1)) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<hydro>/fofc=true with a general EOS, but the "
                    << "derived-variable array wder is not allocated" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
// destructor

Hydro::~Hydro() {
  if (sbf_on && global_variable::my_rank == 0) {
    std::cout << "<hydro>/split_bval_fast: energy-only refreshes " << sbf_nfast
              << ", full refreshes " << sbf_nfull << " (rank 0)" << std::endl;
  }
  if (c2pi_on && global_variable::my_rank == 0) {
    std::cout << "<hydro>/c2p_idempotent: split conversions skipped " << c2pi_nskip
              << ", run " << c2pi_nrun << " (rank 0)" << std::endl;
  }
  if (pam != nullptr) {delete pam;}
  if (pmlat != nullptr) {delete pmlat;}
  if (pravg != nullptr) {delete pravg;}
  delete peos;
  delete pbval_u;
  if (pbval_fs != nullptr) {delete pbval_fs;}
  if (pseamfs != nullptr) {delete pseamfs;}
  if (pbval_ml != nullptr) {delete pbval_ml;}
  if (mlat_own_rec && pmlfs != nullptr) {delete pmlfs;}
  if (pvisc != nullptr) {delete pvisc;}
  if (pcond != nullptr) {delete pcond;}
  if (psrc != nullptr) {delete psrc;}
}

} // namespace hydro
