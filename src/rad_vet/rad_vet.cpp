//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet.cpp
//! \brief RadiationVet constructor: reads the <rad_vet> keys, checks them and allocates
//! the module's arrays.

#include <float.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <utility>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "units/units.hpp"
#include "bvals/bvals.hpp"
#include "coordinates/coordinates.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_clean.hpp"
#include "rad_vet/rad_vet_closure.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "utils/recipe.hpp"
#include "utils/repro_sum.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {
extern bool g_vra_adv;   // vet_rotate_advance (rad_vet_vetgd.cpp)
//----------------------------------------------------------------------------------------
//! \fn RadiationVet::RadiationVet
//! \brief reads <rad_vet>, checks the configuration and allocates the arrays.
//! A key a restart file lacks keeps the default it ran with (fresh defaults: nw).

RadiationVet::RadiationVet(MeshBlockPack *ppack, ParameterInput *pin) :
    u0("vet_u0",1,1,1,1,1),
    u1("vet_u1",1,1,1,1,1),
    coarse_u0("vet_coarse_u0",1,1,1,1,1),
    uflx("vet_uflx",1,1,1,1,1),
    pmy_pack(ppack) {
  // (1) parameters
  c_light = pin->GetReal("rad_vet","c_light");
  if (!(c_light > 0.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/c_light must be positive" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  chat_over_c = pin->GetOrAddReal("rad_vet","chat_over_c",1.0);
  if (!(chat_over_c > 0.0) || chat_over_c > 1.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/chat_over_c must be in (0,1]" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  chat = chat_over_c*c_light;
  // transport: implicit only
  {std::string tr = pin->GetOrAddString("rad_vet","transport","explicit");
  if (tr.compare("implicit") != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/transport = '" << tr << "': only transport = implicit "
      << "(the clean scheme)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  }
  // ---- scheme = clean (rad_vet_clean*.cpp), the only scheme of a fresh implicit run
  const bool nw = !global_variable::restart_run;
  {const bool scd = nw;
  if (scd && !pin->DoesParameterExist("rad_vet","scheme")) {
    recipe::Log("<rad_vet>/scheme = clean (the default with transport = implicit)");
  }
  std::string sc = pin->GetOrAddString("rad_vet","scheme",scd ? "clean" : "default");
  if (sc.compare("clean") == 0) {
    CleanPreflight(pin);
    // clean_w_fse: thin-share thinness g = max(|H|/J, eps C1/1.5), eps = min(kappa_P/chi,
    // 1), w = smoothstep((g - c0)/(c1 - c0)): upwind streaming where the field is
    // anisotropic, not merely non-local
    cln_wfse = pin->GetOrAddBoolean("rad_vet","clean_w_fse",true);
    // clean_w_fse_tau: eps -> 1 - (1 - eps)(1 - exp(-tau_d)), tau_d = chi V/A_f
    cln_wfst = cln_wfse && pin->GetOrAddBoolean("rad_vet","clean_w_fse_tau",true);
    cln_c0 = pin->GetOrAddReal("rad_vet","clean_w_c0",cln_wfse ? 0.3 : 0.1);
    cln_c1 = pin->GetOrAddReal("rad_vet","clean_w_c1",cln_wfse ? 0.6 : 0.5);
    // clean_face_hctau: the bracket's ray face flux F_y interpolated in optical depth
    // between the cell centres, kept where the cells' H jump (CleanStagBracket)
    cln_hctau = pin->GetOrAddBoolean("rad_vet","clean_face_hctau",false);
    // clean_face_jse: face half ranges in the shape of the cells' half-range fluxes
    cln_jse = pin->GetOrAddBoolean("rad_vet","clean_face_jse",true);
    // clean_face_jse_cell: its factor on the donor cell's h+-, not the face value
    cln_jsec = pin->GetOrAddBoolean("rad_vet","clean_face_jse_cell",false);
    // clean_face_area: sp / cs streaming share in area form (CleanStagArea) also under
    // clean_stag_gflux
    cln_farea = pin->GetOrAddBoolean("rad_vet","clean_face_area",false);
    // clean_bc_area: an open x1 end's half-range flux c h+- E in area form, net outward
    // stream h+ + h- x A_c/A_f
    cln_bcarea = pin->GetOrAddBoolean("rad_vet","clean_bc_area",true);
    cln_tol = pin->GetOrAddReal("rad_vet","clean_tol",1.0e-10);
    cln_cnorm = pin->GetOrAddReal("rad_vet","clean_init_cnorm",1.0e-10);
    cln_lcnorm = pin->GetOrAddReal("rad_vet","clean_lin_cnorm",1.0e-10);
    // clean_maxnl, clean_maxit, clean_conv_max do not take implicit_maxit,
    // implicit_lin_maxit, implicit_tol_rel_max (a named other value prints a NOTE).
    // clean_lin_pass: a linear solve ending at clean_maxit leaves its iterate to the next
    // pass (no stall test)
    {Real ov[3] = {0.0, 0.0, 0.0};   // named implicit_maxit, _lin_maxit, _tol_rel_max
    const char *ok[3] = {"implicit_maxit", "implicit_lin_maxit", "implicit_tol_rel_max"};
    const char *ck[3] = {"clean_maxnl", "clean_maxit", "clean_conv_max"};
    for (const InputBlock &ib : pin->block) {
      if (ib.block_name != "rad_vet") continue;
      for (const InputLine &ln : ib.line) {
        if (ln.param_comment == ParameterInput::default_comment) continue;
        for (int n = 0; n < 3; ++n) {
          if (ln.param_name == ok[n]) ov[n] = pin->GetReal("rad_vet", ok[n]);
        }
      }
    }
    bool nck[3];
    for (int n = 0; n < 3; ++n) nck[n] = pin->DoesParameterExist("rad_vet", ck[n]);
    cln_maxnl = pin->GetOrAddInteger("rad_vet","clean_maxnl",200);
    cln_maxit = pin->GetOrAddInteger("rad_vet","clean_maxit",800);
    cln_bmax = pin->GetOrAddReal("rad_vet","clean_conv_max",1.0e-5);
    const Real cv[3] = {static_cast<Real>(cln_maxnl), static_cast<Real>(cln_maxit),
                        cln_bmax};
    for (int n = 0; n < 3; ++n) {
      if (global_variable::my_rank == 0 && ov[n] > 0.0 && !nck[n] && ov[n] != cv[n]) {
        std::cout << "### NOTE: <rad_vet>/" << ok[n] << " = " << ov[n] << " is not the "
                  << "clean solve's: " << ck[n] << " = " << cv[n] << " (name " << ck[n]
                  << " to change it)" << std::endl;
      }
    }
    const bool pl = (ov[1] > 0.0);
    cln_linpass = pin->GetOrAddBoolean("rad_vet","clean_lin_pass",nw || pl);}
    cln_tol_nl = pin->GetOrAddReal("rad_vet","clean_tol_nl",1.0e-11);
    cln_raymin = pin->GetOrAddReal("rad_vet","clean_ray_min",1.0e-3);
    // clean_rays = sc only
    if (pin->DoesParameterExist("rad_vet","clean_rays") &&
        pin->GetString("rad_vet","clean_rays") != "sc") {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_vet>/clean_rays = sc (the FV rays fv1 | fv | lld were "
        << "removed)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_rot_ring = off | live | jump (live where the vet_gd rays rotate with a
    // cross-fade, CleanInit): the rotating union swept as the current (+ next) set, older
    // sets as moment offsets (CleanRingMoments). clean_rot_fuse: the live ring's two sets
    // in one fused sweep (bitwise the two sweeps; a second intensity array)
    cln_rfuse = pin->GetOrAddBoolean("rad_vet","clean_rot_fuse",nw);
    if (pin->DoesParameterExist("rad_vet","clean_rot_ring")) {
      std::string rr = pin->GetString("rad_vet","clean_rot_ring");
      cln_ring = (rr.compare("live") == 0) ? 1 : ((rr.compare("jump") == 0) ? 2 : 0);
      cln_ring_named = true;
      if (cln_ring == 0 && rr.compare("off") != 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "<rad_vet>/clean_rot_ring = off | live | jump" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    cln_fillpass = pin->GetOrAddInteger("rad_vet","clean_fill_passes",8);
    // clean_fill_cs_remap: the fill's h+- of a cs seam neighbour in our face families
    cln_fillcs = pin->GetOrAddBoolean("rad_vet","clean_fill_cs_remap",false);
    cln_pred = pin->GetOrAddBoolean("rad_vet","clean_predict",true);
    // implicit_predictor (step unless named none): the stage solve starts from the
    // predicted E (CleanFirstIterate)
    cln_kpred = !(pin->DoesParameterExist("rad_vet","implicit_predictor") &&
                  !pin->IsParameterDefaulted("rad_vet","implicit_predictor") &&
                  pin->GetString("rad_vet","implicit_predictor") == "none");
    // clean_pred_stage: sdirk2's stage-2 solve starts from E_old + h K1 of this step's
    // stage 1, not of the last step's stage 2
    cln_pst = pin->GetOrAddBoolean("rad_vet","clean_pred_stage",false);
    // clean_gas_secant = false | true | stage2 | max | max_stage2: exchange slopes from
    // the chord through the last pass's gas T (max: only where steeper than the tangent)
    {const std::string gs = pin->GetOrAddString("rad_vet","clean_gas_secant","false");
    cln_gsec = (gs == "true") ? 1 : ((gs == "stage2") ? 2 : ((gs == "max") ? 3 :
               ((gs == "max_stage2") ? 4 : 0)));}
    // clean_pred_xch: stage-2 start E_old + h K, K the stage-1 slope where the exchange
    // is slow, the last stage-2 slope where stiff (weights 1 : z,
    // z = h (ch kE + 4 c kP a T^4/e))
    cln_pxc = pin->GetOrAddBoolean("rad_vet","clean_pred_xch",false);
    {std::string bf = pin->GetOrAddString("rad_vet","clean_bc","halfrange");
    cln_bdry_form = (bf.compare("marshak") == 0) ? 1 : 0;}
    // clean_layout = staggered only
    if (pin->GetOrAddString("rad_vet","clean_layout",nw ? "staggered" : "cell") !=
        "staggered") {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_vet>/clean_layout = staggered (the cell layout was removed)"
        << std::endl;
      std::exit(EXIT_FAILURE);
    }
    cln_pcr = pin->GetOrAddBoolean("rad_vet","clean_stag_pcr",true);
    // clean_stag_opac_update (default on for a fresh run, every mesh) and clean_face_opac
    // (default harmonic on a fresh run), USER 11-21, recorded; a restart file without a
    // key keeps the old default (echoed)
    auto old = [&](const char *k, const char *v) {
      if (!nw && !pin->DoesParameterExist("rad_vet", k)) {
        recipe::Log(std::string("<rad_vet>/") + k + " = " + v + " (a restart file "
                    "without the key: the old default)");
      }
    };
    old("clean_stag_opac_update", "false");
    old("clean_face_opac", "arith");
    cln_opt = pin->GetOrAddBoolean("rad_vet","clean_stag_opac_update",nw);
    // clean_stag_opac_passes: opacities of passes 0 .. N-1, then frozen
    cln_optn = pin->GetOrAddInteger("rad_vet","clean_stag_opac_passes",3);
    // clean_stag_opac_drag: a refresh rebuilds the clean_stag_vimp drag share b and Q
    // (CleanStagDrag)
    old("clean_stag_opac_drag", "false");
    cln_optdr = pin->GetOrAddBoolean("rad_vet","clean_stag_opac_drag",nw);
    {const std::string fo = pin->GetOrAddString("rad_vet","clean_face_opac",
                                                nw ? "harmonic" : "arith");
    if (fo != "arith" && fo != "harmonic") {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_face_opac = arith | harmonic"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    cln_hface = (fo == "harmonic");}
    // clean_stag_fs_line = pcr | thomas: the F solve's x1 lines on a GPU; old files'
    // clean_stag_fs_pcr (true) / clean_stag_fs_pcr_cs (cs: false) when unnamed
    {const bool cs = pmy_pack->pmesh->use_cubed_sphere;
    auto named = [&](const char *k, const bool d) {
      return pin->DoesParameterExist("rad_vet", k) ? pin->GetBoolean("rad_vet", k) : d;
    };
    if (cs) old("clean_stag_fs_line", "thomas");
    const bool pcr = named("clean_stag_fs_pcr", true) &&
                     (!cs || named("clean_stag_fs_pcr_cs", nw));
    const std::string fl = pin->GetOrAddString("rad_vet","clean_stag_fs_line",
                                               pcr ? "pcr" : "thomas");
    if (fl != "pcr" && fl != "thomas") {
      std::cout << "### FATAL ERROR: <rad_vet>/clean_stag_fs_line = pcr | thomas"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    cln_fspcr = (fl == "pcr");}
    cln_stall_n = pin->GetOrAddInteger("rad_vet","clean_bicg_stall_its",
                                      cln_linpass ? 0 : 50);
    cln_stall_f = pin->GetOrAddReal("rad_vet","clean_bicg_stall_fac",0.99);
    // face tie: a face two MeshBlocks share (or one across a periodic wrap) takes its
    // owner's F, momentum exchange, energy-flux memory and Phi (CleanStagFaceTie).
    // clean_stag_tie_gid: on every mesh the cs rule (CleanStagCsTie), owner the lower gid
    old("clean_stag_tie_gid", "false");
    cln_tgid = pin->GetOrAddBoolean("rad_vet","clean_stag_tie_gid",
                                    nw && !pmy_pack->pmesh->multilevel);
    if (cln_tgid && pmy_pack->pmesh->multilevel) {
      std::cout << "### FATAL ERROR: <rad_vet>/clean_stag_tie_gid needs a uniform mesh"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_stag_fsweeps N: line solves of the F solve up to N times
    if (pin->DoesParameterExist("rad_vet","clean_stag_fsweeps")) {
      cln_fsweeps = std::max(pin->GetInteger("rad_vet","clean_stag_fsweeps"), 1);
    }
    // clean_stag_fschur (default on where it applies): F lines exact across MeshBlocks
    // and periodic wraps (CleanStagFSchur); needs a uniform mesh and MeshBlocks of >= 3
    // cells along a line (>= 6 on longer or open lines)
    {Mesh *pm = pmy_pack->pmesh;
    // (cubed sphere: a lateral line is a ring of four panels' MeshBlocks)
    const bool cs = pm->use_cubed_sphere;
    bool ok = cln_fsweeps == 1 && !pm->multilevel, lag = false;
    for (int d = 0; d < 3; ++d) {
      // (layout_invariant_mode = tile: tiles as MeshBlocks, CleanStagFSchurTile)
      const int n = (d > 0 && repro::tile_n > 0) ? repro::tile_n :
          ((d == 0) ? pm->mb_indcs.nx1 :
                      ((d == 1) ? pm->mb_indcs.nx2 : pm->mb_indcs.nx3));
      const int nb = ((cs && d > 0) ? 4 : 1)*((d == 0) ? pm->mesh_indcs.nx1 :
                     ((d == 1) ? pm->mesh_indcs.nx2 : pm->mesh_indcs.nx3))/n;
      const bool per = (cs && d > 0) || (pm->mesh_bcs[2*d] == BoundaryFlag::periodic);
      const int nlong = (cs && per && nb%2 == 0) ? 7 : 6;   // cs: the kept segment
      const int nmin = (nb > 2 || (nb == 2 && !per)) ? nlong : ((nb == 2) ? 3 : 0);
      ok = ok && n >= nmin && (nmin == 0 || pm->mb_indcs.ng >= 3);
      lag = lag || nmin > 0;
    }
    cln_fschur = ok;
    if (pin->DoesParameterExist("rad_vet","clean_stag_fschur")) {
      cln_fschur = pin->GetBoolean("rad_vet","clean_stag_fschur");
      if (cln_fschur && !ok) {
        std::cout << "### FATAL ERROR: <rad_vet>/clean_stag_fschur needs "
                  << "clean_stag_fsweeps = 1, a uniform mesh, "
                  << "nghost >= 3 and MeshBlocks of >= 3 cells along a line of two "
                     "round a "
                  << "periodic one, >= 6 along longer or open ones" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    } else if (!ok && lag && !cs && cln_fsweeps == 1) {   // no silent lagged ends
      std::cout << "### FATAL ERROR: <rad_vet>/clean_stag_fschur cannot take this mesh "
                << "(nghost >= 3, MeshBlock size, uniform mesh): name clean_stag_fsweeps "
                << "(line ends iterated) or clean_stag_fschur = false" << std::endl;
      std::exit(EXIT_FAILURE);
    } else if (!ok && global_variable::my_rank == 0) {
      std::cout << "<rad_vet> clean_stag_fschur off: the F solve's line ends lagged "
                << "(clean_stag_fsweeps, mesh or MeshBlock size)" << std::endl;
    }
    // layout_tile: a periodic line joins whole MeshBlocks (CleanStagFSchurTile)
    for (int d = 1; d < 3 && cln_fschur && repro::tile_n > 0 && !cs; ++d) {
      const int nx = (d == 1) ? pm->mb_indcs.nx2 : pm->mb_indcs.nx3;
      const int nbm = ((d == 1) ? pm->mesh_indcs.nx2 : pm->mesh_indcs.nx3)/nx;
      if (nx > repro::tile_n && nbm%2 == 1 && nbm > 1 &&
          pm->mesh_bcs[2*d] == BoundaryFlag::periodic) {
        std::cout << "### FATAL ERROR: <mesh>/layout_tile with clean_stag_fschur: a "
                  << "periodic line of MeshBlocks of more than one tile needs an even "
                  << "number of them (or one)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    // cubed sphere without fschur iterates its seams (fsweeps 64, Thomas lines)
    if (cs && !cln_fschur && !pin->DoesParameterExist("rad_vet","clean_stag_fsweeps")) {
      cln_fsweeps = 64;
    }
    if (cs && !cln_fschur && cln_fspcr && global_variable::my_rank == 0) {
      std::cout << "<rad_vet> clean_stag_fs_line = pcr: thomas on the cubed sphere "
                << "without clean_stag_fschur" << std::endl;
    }}
    // clean_consistency: default on Cartesian, and on sp / cs with the vet_gd rays (sc)
    const bool curv = pmy_pack->pmesh->use_spherical_polar ||
                      pmy_pack->pmesh->use_cubed_sphere;
    const bool vg = pin->DoesParameterExist("rad_vet","vet_gd") &&
                    pin->GetBoolean("rad_vet","vet_gd");
    // sp / cs: clean_rays recorded; a restart file without it ran other rays: fatal
    if (curv && vg && !pin->DoesParameterExist("rad_vet","clean_rays")) {
      if (global_variable::restart_run) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "<rad_vet>/clean_rays: the restart file ran the removed FV "
                          "rays "
          << "(fv1); name clean_rays = sc" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      pin->GetOrAddString("rad_vet","clean_rays","sc");
    }
    cln_split = true;
    cln_consist = pin->GetOrAddBoolean("rad_vet","clean_consistency", !curv || vg);
    cln_maxpic = pin->GetOrAddInteger("rad_vet","clean_maxpic",12);
    // staggered: energy flux apart from the face momentum (split), face momentum of its
    // own F solve (fsolve) in the Newton loop (fwork); with clean_consistency the rays'
    // kinetic momentum flux (fkin) and face flux as flux factors (gflux); donor-cell
    // advective share (upwind). Read only when named.
    auto stag = [&](const char *key, const bool def) {
      return pin->DoesParameterExist("rad_vet", key) ? pin->GetBoolean("rad_vet", key) :
                                                      def;
    };
    cln_fsolve = cln_upw = cln_split;
    cln_fkin = cln_fwork = cln_gflux = cln_split && cln_consist;
    // clean_stag_freal: realisable sb (comoving part <= 2 |beta| hc kE E_k) and restart
    // faces <= c E of their donor
    cln_freal = stag("clean_stag_freal", cln_gflux && cln_fsolve);
    // clean_stag_upwind_central: an upwind face carries the central mean of a E within
    // the TVD bound (CleanStagAdvectCorr)
    cln_ucen = stag("clean_stag_upwind_central", false);
    // clean_stag_central: the rest as the plain central mean, the work drain and the gas
    // kick in the same Patankar weight (CleanStagMprkWeight)
    cln_ucp = stag("clean_stag_central", false);
    cln_ucen = cln_ucen || cln_ucp;
    // clean_stag_csvec_radial: cs cell vector exact for a radial field (true only)
    cln_csvr = stag("clean_stag_csvec_radial", true);
    if (!cln_csvr) {
      std::cout << "### FATAL ERROR: <rad_vet>/clean_stag_csvec_radial = false (the old "
                << "cell vector) is gone" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_stag_thinf: the bracket's ray correction in the face half ranges, realisable
    // (CleanStagBracket)
    cln_thinf = pin->GetOrAddBoolean("rad_vet","clean_stag_thinf",nw);
    // clean_stag_flim: false only (fwb + fcone bound F)
    if (pin->GetOrAddBoolean("rad_vet","clean_stag_flim",false)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_stag_flim = true was removed (only "
                                "false)"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_stag_fwb: F rows in deviation form about the rays' own state, their exact
    // steady state (CleanStagFSolve)
    old("clean_stag_fwb", "false");
    old("clean_stag_fcone", "false");
    old("clean_stag_fwb_stage", "false");
    cln_fwb = pin->GetOrAddBoolean("rad_vet","clean_stag_fwb",
                                   nw && cln_fsolve && cln_consist);
    if (cln_fwb && !(cln_fsolve && cln_consist)) {
      std::cout << "### FATAL ERROR: <rad_vet>/clean_stag_fwb needs the staggered F "
                   "solve "
                << "with clean_consistency" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_stag_fwb_stage: fwb's residual s of the stage-old state, once per stage;
    // off: of each solve's E (thin-row F = the rays' flux factor times E)
    cln_fwbst = cln_fwb && pin->GetOrAddBoolean("rad_vet","clean_stag_fwb_stage",false);
    // clean_stag_fwb_t2: fwb_stage's s of the stage-time state E_old + h K (second order
    // in time)
    old("clean_stag_fwb_t2", "false");
    cln_fwbt2 = pin->GetOrAddBoolean("rad_vet","clean_stag_fwb_t2",nw && cln_fwbst);
    if (cln_fwbt2 && cln_fwb && !cln_fwbst) {
      std::cout << "### FATAL ERROR: <rad_vet>/clean_stag_fwb_t2 needs "
                << "clean_stag_fwb_stage" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    cln_fwbt2 = cln_fwbt2 && cln_fwbst;
    // clean_stag_fwb_line / _ring: fwb's s of the x1 faces less its line mean, of the phi
    // faces (3-D sp, periodic x3) less its ring mean: no net radial force on a line, no
    // torque on a ring (CleanStagFwbProject)
    cln_fwbpj[0] = cln_fwb && stag("clean_stag_fwb_line", false);
    const bool ring_ok = pmy_pack->pmesh->use_spherical_polar &&
        pmy_pack->pmesh->three_d &&
        !pmy_pack->pmesh->multilevel &&
        pmy_pack->pmesh->mesh_bcs[BoundaryFace::inner_x3] == BoundaryFlag::periodic;
    old("clean_stag_fwb_ring", "false");
    cln_fwbpj[2] = cln_fwb && pin->GetOrAddBoolean("rad_vet","clean_stag_fwb_ring",
                                                   nw && cln_fwb && ring_ok);
    if ((cln_fwbpj[0] || cln_fwbpj[2]) && pmy_pack->pmesh->multilevel) {
      std::cout << "### FATAL ERROR: <rad_vet>/clean_stag_fwb_line, _ring need a uniform "
                << "mesh" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (cln_fwbpj[2] && !(pmy_pack->pmesh->use_spherical_polar &&
        pmy_pack->pmesh->three_d &&
        pmy_pack->pmesh->mesh_bcs[BoundaryFace::inner_x3] == BoundaryFlag::periodic)) {
      std::cout << "### FATAL ERROR: <rad_vet>/clean_stag_fwb_ring needs 3-D spherical "
                << "polar with a periodic x3" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_stag_fwb_cart (read only when named, default false): on a uniform Cartesian
    // mesh fwb's s of the faces of every direction less its line mean, as fwb_line's x1
    // faces: linear momentum exact without a ledger; no-op on curvilinear meshes
    if (cln_fwb && pin->DoesParameterExist("rad_vet", "clean_stag_fwb_cart") &&
        pin->GetBoolean("rad_vet", "clean_stag_fwb_cart") &&
        !pmy_pack->pmesh->multilevel && !pmy_pack->pmesh->use_spherical_polar &&
        !pmy_pack->pmesh->use_cubed_sphere) {
      cln_fwbpj[0] = cln_fwbpj[1] = cln_fwbpj[2] = true;
    }
    // clean_stag_fcone (default on for a fresh run with fwb, recorded; a restart file
    // without the key: off): every cell's flux vector within its realisable cone,
    // |F| <= c E, by construction: the faces scaled towards 0 (CleanStagFCone)
    cln_fcone = pin->GetOrAddBoolean("rad_vet","clean_stag_fcone",nw && cln_fwb);
    if (cln_fcone && !cln_fsolve) {
      std::cout << "### FATAL ERROR: <rad_vet>/clean_stag_fcone needs the staggered F "
                << "solve" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_stag_fcone_gas: the gas takes the momentum fcone cuts at the face (momentum
    // and L_z exact; E pays the kinetic energy)
    cln_fcgas = cln_fcone && pin->GetOrAddBoolean("rad_vet","clean_stag_fcone_gas",false);
    // clean_stag_fcone_lz (read only when named, default off): fcone takes a cell's
    // excess from its faces without L_z lever arm first (CleanLzFamily: spherical polar
    // x1, x2; cubed sphere x1; Cartesian x3), cut and counted as now; a face with an arm
    // is scaled only where those faces alone exceed c E, and its cut is the gas's at the
    // face (as fcone_gas; clean_lz_spin: a rigid rotation), so that L_z is exact (cubed
    // sphere: CleanCsTorque)
    cln_fclz = cln_fcone && pin->DoesParameterExist("rad_vet","clean_stag_fcone_lz") &&
               pin->GetBoolean("rad_vet","clean_stag_fcone_lz");
    // clean_stag_fcone_emax (read only when named; needs clean_stag_fcone): every inner
    // face's F within c max(E_L, E_R), also an opaque face's: no face F streams more than
    // its brighter cell holds, so the gas's face terms stay below the cap h chi_f E
    if (pin->DoesParameterExist("rad_vet","clean_stag_fcone_emax")) {
      cln_fcemx = pin->GetBoolean("rad_vet","clean_stag_fcone_emax");
      if (cln_fcemx && !cln_fcone) {
        std::cout << "### FATAL ERROR: <rad_vet>/clean_stag_fcone_emax needs "
                  << "clean_stag_fcone" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    // clean_stag_xdonor_freeze (0 off): from this pass on the xdonor cross flux C is held
    // (no donor/cap flip between passes)
    cln_xdfz = pin->GetOrAddInteger("rad_vet","clean_stag_xdonor_freeze",4);
    if (const int pd = pin->GetOrAddInteger("rad_vet","clean_pass_damp",0); pd != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_pass_damp = " << pd
                << " was removed (only 0)"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_stag_xdonor = lag | p1 | off: the off-diagonal stress flux as an M-matrix
    // term (CleanStagXDonor); clean_stag_xsrc: also the face sources past their explicit
    // share
    const std::string xd = pin->GetOrAddString("rad_vet","clean_stag_xdonor",
                                               nw ? "lag" : "off");
    if (xd != "off" && xd != "p1" && xd != "lag") {
      std::cout << "### FATAL ERROR <rad_vet>/clean_stag_xdonor = off | p1 | lag"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    cln_xd = (xd != "off");
    cln_xlag = (xd == "lag");
    cln_xsrc = cln_xd && pin->GetOrAddBoolean("rad_vet","clean_stag_xsrc",nw);
    // clean_stag_eicarry: a stage start's gas eint < 0 enters the gas row as is, paid by
    // its implicit exchange (CleanOldVector)
    cln_eicry = stag("clean_stag_eicarry", false);
    cln_colf = stag("clean_stag_colflux", false);
    // clean_cs_lz_abs = 2 (cubed sphere, the only form): the gas phi momentum is the
    // faces' absorption as its r and theta; the L_z it does not take is booked as the
    // radiation's (CleanLzBook)
    if (pmy_pack->pmesh->use_cubed_sphere) cln_lza = 2;
    if (pin->DoesParameterExist("rad_vet", "clean_cs_lz_abs") &&
        pin->GetInteger("rad_vet", "clean_cs_lz_abs") != 2) {
      std::cout << "### FATAL ERROR: <rad_vet>/clean_cs_lz_abs = 0 | 1 (the old cubed-"
                << "sphere gas phi momentum) is gone; 2 is the only form" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    cln_lzbook = (cln_lza == 2) || (pmy_pack->pmesh->use_spherical_polar &&
                                    pmy_pack->pmesh->three_d);
    // clean_cs_lz_tel: L_z by telescoping (gas phi momentum of CleanCsTorque, radiation
    // L_z of its face momenta), so the ledger holds the radial boundaries' torque alone
    if (cln_lza == 2 || pin->DoesParameterExist("rad_vet", "clean_cs_lz_tel")) {
      cln_lzt = pin->GetOrAddBoolean("rad_vet","clean_cs_lz_tel",nw);
    }
    // clean_cs_lz_bnd_free: no radiation stress torque through the cs radial boundaries
    if (cln_lzt) cln_lzbf = pin->GetOrAddBoolean("rad_vet","clean_cs_lz_bnd_free",nw);
    if (cln_lzt && cln_lza != 2) {
      std::cout << "### FATAL ERROR: <rad_vet>/clean_cs_lz_tel needs the cubed sphere "
                << "and clean_layout = staggered" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_lz_spin (read only when named; needs clean_stag_fcone_lz, on the cubed sphere
    // clean_cs_lz_tel): the L_z of the momentum that fcone's cut of the faces with an
    // arm and clean_wb_mcap's cap withhold is the gas's as one rigid rotation dOmega =
    // X/I, I = sum rho V R^2 (CleanWriteBack), not the face's gas: L_z exact, no kick
    if (pin->DoesParameterExist("rad_vet","clean_lz_spin")) {
      cln_lzsp = pin->GetBoolean("rad_vet","clean_lz_spin");
      if (cln_lzsp && (!cln_fclz || (pmy_pack->pmesh->use_cubed_sphere && !cln_lzt))) {
        std::cout << "### FATAL ERROR: <rad_vet>/clean_lz_spin needs clean_stag_fcone_lz "
                  << "(the cubed sphere: and clean_cs_lz_tel)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    if (cln_lzbook) {
      const bool cs = (cln_lza == 2);
      cln_lzpin = pin;
      for (int k1 = 0; k1 < 2; ++k1) {   // a restart file's non-finite value: NOTE, 0
        const std::string key = M1LzKey(pin, cs, k1 == 1);
        Real &v = (k1 == 1) ? cln_lzs1 : cln_lzl;
        v = pin->GetOrAddReal("rad_vet", key, 0.0);
        if (std::isfinite(v)) continue;
        if (global_variable::my_rank == 0) {
          std::cout << "### NOTE: <rad_vet>/" << key << " = " << v << " in the restart "
                    << "file is not finite: reset to 0 (the L_z ledger restarts from 0)"
                    << std::endl;
        }
        v = 0.0;
        pin->SetReal("rad_vet", key, 0.0);
      }
    }
    // clean_xt_cons: first iterate's gas T with G = 1 - beta^2/2, then T from the pass's
    // gas row; clean_xt_norm: the local exchange with E responding,
    // e(T) + A T^4/(1 + hc kE) = eint + hc kE E/(1 + hc kE)
    cln_xtcons = stag("clean_xt_cons", false);
    // clean_xt_exact = cons | off | x2: one A' per cell and pass (CleanRightHandSide);
    // cons: the work by which the row and the gas kinetic energy differ goes to E whole;
    // x2: also every T root from the local 2x2 of the E and gas rows
    {
      const std::string xe = pin->GetOrAddString("rad_vet","clean_xt_exact",
                                                 nw ? "cons" : "off");
      if (xe != "off" && xe != "cons" && xe != "x2") {
        std::cout << "### FATAL ERROR <rad_vet>/clean_xt_exact = off | cons | x2"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      cln_xte = (xe == "cons") ? 1 : ((xe == "x2") ? 2 : 0);
      if (cln_xte > 0) cln_xtcons = true;
    }
    if (cln_ucp && cln_xte != 1) {   // CLW_PD: the work's q
      std::cout << "### FATAL ERROR <rad_vet>/clean_stag_central needs clean_xt_exact = "
                << "cons" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    cln_xtnorm = stag("clean_xt_norm", false) || cln_xte == 2;
    // clean_stag_rayratio: the bracket reads the rays only as ratios times the state's E
    cln_rayratio = pin->GetOrAddBoolean("rad_vet","clean_stag_rayratio",false) &&
                   cln_gflux;
    // clean_ring_jmin: a ring J ratio J_s/J_b only where both J >= clean_ray_min E
    // (elsewhere the offset is invalid)
    cln_ringjmin = pin->GetOrAddBoolean("rad_vet","clean_ring_jmin",true);
    // clean_rot_ring_fill_rst: at a restart from a file without the clean state the ring
    // fill takes the older sets' offsets swept at the restart state (CleanRingMoments)
    cln_rfrst = pin->GetOrAddBoolean("rad_vet","clean_rot_ring_fill_rst",true);
    // clean_rst_mshare: at such a restart the M1 share's memory starts at the first
    // solution's steady flux (CleanStagFSolve)
    cln_rmsh = pin->GetOrAddBoolean("rad_vet","clean_rst_mshare",true);
    // clean_source_convex: the rays' vet_scatter source with kappa_P/chi capped at 1,
    // S = J + min(kappa_P/chi, 1)(B - J) between J and B
    cln_srccvx = pin->GetOrAddBoolean("rad_vet","clean_source_convex",true);
    if (cln_freal && !(cln_gflux && cln_fsolve)) {
      std::cout << "### FATAL ERROR: <rad_vet>/clean_stag_freal needs clean_consistency "
                << "(the F solve and gflux)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    cln_vimp = stag("clean_stag_vimp", cln_split);
    // cubed sphere: only with clean_stag_vimp_cs
    if (cln_vimp && pmy_pack->pmesh->use_cubed_sphere) old("clean_stag_vimp_cs", "false");
    if (cln_vimp && pmy_pack->pmesh->use_cubed_sphere &&
        !pin->GetOrAddBoolean("rad_vet","clean_stag_vimp_cs",nw)) {
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_vet> staggered on the cubed sphere: clean_stag_vimp = false"
                  << std::endl;
      }
      cln_vimp = false;
    }
    // clean_clos1_budget (read only when named; false, the default: the strict norm): a
    // pass whose iterate the next pass takes the closure of in the budget's norm, to 0.1
    // clean_conv_max, the floor of the passes after it at the budget's ceiling (a step
    // dt_retry redoes strictly, as its budget)
    cln_c1b = pin->DoesParameterExist("rad_vet","clean_clos1_budget") &&
              pin->GetBoolean("rad_vet","clean_clos1_budget");
    // clean_vel_settled (read only when named; false, the default: every pass): the
    // velocity update (its F solve) only in a pass that can end the solve, one whose E
    // change settled; the others leave v, F to the next pass (the same fixed point)
    cln_vset = pin->DoesParameterExist("rad_vet","clean_vel_settled") &&
               pin->GetBoolean("rad_vet","clean_vel_settled");
    // clean_bud_maxnorm (read only when named; false, the default: the L2 norm): the
    // budget passes' linear test and Eisenstat-Walker |r0| in the max norm over the cells,
    // the norm of the pass's own test (CleanNewtonUpdate), not an L2 sum over the cells
    cln_bmn = pin->DoesParameterExist("rad_vet","clean_bud_maxnorm") &&
              pin->GetBoolean("rad_vet","clean_bud_maxnorm");
    // clean_gas_tfloor (default true): the gas row's T bounded below by the EOS
    // tfloor, and the gas eint written back from that root (E takes the rest)
    cln_gtf = stag("clean_gas_tfloor", true);
    // clean_wb_fct: the write-back's Zalesak limit on the written faces' energy flux,
    // E > 0 by construction (CleanWbFct)
    old("clean_wb_fct", "false");
    cln_wbfct = pin->GetOrAddBoolean("rad_vet","clean_wb_fct",
                                     nw && cln_gtf && cln_xte > 0);
    if (cln_wbfct && !(cln_gtf && cln_xte > 0)) {
      std::cout << "### FATAL ERROR <rad_vet>/clean_wb_fct needs clean_gas_tfloor and "
                << "clean_xt_exact = cons | x2" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_wb_floor_fail: a write-back whose E floor takes energy its own gas does not
    // give is not admissible (stage fails; BE step to dt_retry, efloor_column); accepted
    // at the last cut with a CLEAN WARNING
    if (pin->DoesParameterExist("rad_vet","clean_wb_floor_fail")) {
      cln_wbff = pin->GetBoolean("rad_vet","clean_wb_floor_fail");
    }
    // clean_wb_exch_local (needs clean_gas_tfloor, EOS tfloor > 0): the write-back's cell
    // exchange solved exactly from the transported E (E > 0, gas >= e(rho, tfloor))
    if (pin->DoesParameterExist("rad_vet","clean_wb_exch_local")) {
      cln_wbxl = pin->GetBoolean("rad_vet","clean_wb_exch_local");
    }
    // clean_wb_mcap (read only when named; default false): the write-back gives the gas
    // the face momentum terms at the Newton passes' cap h chi_f E, the momentum whose
    // work the E row paid (CleanMomentumChange); the terms withheld are printed
    if (pin->DoesParameterExist("rad_vet","clean_wb_mcap")) {
      cln_wbmc = pin->GetBoolean("rad_vet","clean_wb_mcap");
    }
    // clean_wb_floor_diag N (diagnostic, read only when named): the terms of the first
    // N cells of each write-back whose E the floor raises (CLEAN WBFDIAG)
    if (pin->DoesParameterExist("rad_vet","clean_wb_floor_diag")) {
      cln_wbfdg = pin->GetInteger("rad_vet","clean_wb_floor_diag");
    }
    // clean_gas_newton_switch W (0 off, W in [3, 64]), _min, clean_stall_fac: a stalled
    // solve's remaining passes on the bracketed gas root without the Newton slope
    if (!nw && !pin->DoesParameterExist("rad_vet","clean_gas_newton_switch")) {
      recipe::Log("<rad_vet>/clean_gas_newton_switch = 0 (a restart file without the "
                  "key: "
                  "the old default)");
    }
    cln_gnsw = pin->GetOrAddInteger("rad_vet","clean_gas_newton_switch",nw ? 4 : 0);
    cln_gnswmin = pin->GetOrAddInteger("rad_vet","clean_gas_newton_switch_min",
                                       nw ? 6 : 20);
    cln_stfac = pin->GetOrAddReal("rad_vet","clean_stall_fac",0.5);
    if (cln_gnsw < 0 || cln_gnsw > 64 || (cln_gnsw > 0 && cln_gnsw < 3) ||
        cln_gnswmin < 1 || !(cln_stfac > 0.0 && cln_stfac <= 1.0)) {
      std::cout << "### FATAL ERROR <rad_vet>/clean_gas_newton_switch must be 0 or in "
                << "[3, 64], clean_gas_newton_switch_min >= 1, clean_stall_fac in (0, 1]"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // preconditioner: mg (CleanStagVCycle)
    {// clean_stag_mg_color = mirror | parity: the mg's line colours (CleanStagLineSolve)
    const std::string mc = pin->GetOrAddString("rad_vet","clean_stag_mg_color",
                                               nw ? "mirror" : "parity");
    if (mc != "parity" && mc != "mirror") {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_stag_mg_color = parity | mirror"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    cln_stagmg = (mc == "mirror") ? 2 : 1;
    // layout_tile N: the mg acts on fixed N x N (x2, x3) tiles as on MeshBlocks, so the
    // solve is the same for every layout of whole tiles (one MeshBlock along x1)
    if (repro::tile_n > 0) {
      const int nt = repro::tile_n;
      const auto &mi = pmy_pack->pmesh->mb_indcs;
      const bool d3 = pmy_pack->pmesh->three_d;
      if (nt < 1 || mi.nx2 % nt != 0 || (d3 && mi.nx3 % nt != 0) ||
          pmy_pack->pmesh->nmb_rootx1 != 1) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<mesh>/layout_tile N needs N dividing meshblock nx2 "
                  << "(and nx3) and one MeshBlock along x1" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      cln_mgt = {mi.nx1, nt, d3 ? nt : 1};
    }
    // clean_pc_sweep = sym | rb | rbr | rbrsym: the mg smoothing per level, a forward
    // red-black sweep, + a second red pass (rbr), + a black-red post-sweep (sym)
    const std::string sw = pin->GetOrAddString("rad_vet","clean_pc_sweep",
                                               nw ? "sym" : "rb");
    cln_pcsweep = (sw == "rbr") + 2*(sw == "sym") + 3*(sw == "rbrsym");
    if (!cln_pcsweep && sw != "rb") {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_pc_sweep = " << sw << ": rb | rbr | sym "
                << "| rbrsym (wave was removed)" << std::endl;
      std::exit(EXIT_FAILURE);
    }}
    // clean_stag_lat_sums (sp / Cartesian, with clean_stag_xdonor): the faces' lateral
    // stress from the cells' sums Y (CleanStagCsLateral), the cs form; off: 4-cell
    // central differences
    old("clean_stag_lat_sums", "false");
    cln_ysum = pin->GetOrAddBoolean("rad_vet","clean_stag_lat_sums",
                                    nw && cln_xd) &&
               pmy_pack->pmesh->multi_d && !pmy_pack->pmesh->use_cubed_sphere;
    if (cln_ysum && !cln_xd) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_stag_lat_sums: with clean_stag_xdonor = "
                << "lag | p1" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_closure: predict | iterate | extrap | extrap1 | extrap2 | eddington | sphere
    {std::string cl = pin->GetOrAddString("rad_vet","clean_closure","extrap1");
    cln_closit = (cl.compare("iterate") == 0);
    cln_extrap = (cl.compare("extrap") == 0) ? 1 : ((cl.compare("extrap2") == 0) ? 2 :
                 ((cl.compare("extrap1") == 0) ? 3 : 0));
    cln_edd = (cl.compare("eddington") == 0);
    cln_sphere = (cl.compare("sphere") == 0);
    if (!cln_closit && cln_extrap == 0 && !cln_edd && !cln_sphere &&
        cl.compare("predict") != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_vet>/clean_closure = predict | iterate | extrap | extrap1 "
                        "| "
        << "extrap2 | eddington | sphere" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (cln_sphere) cln_sph_r = pin->GetReal("rad_vet","clean_sphere_radius");}
    {const std::string nn = pin->GetOrAddString("rad_vet","clean_nl_norm","global");
    if (nn.compare("global") != 0 && nn.compare("cell") != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_vet>/clean_nl_norm = global | cell" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    cln_nl_cell = (nn.compare("cell") == 0);}
    // clean_troot_rtol: the exchange root stops at its residual noise
    cln_trtol = pin->GetOrAddReal("rad_vet","clean_troot_rtol",1.0e-13);
    if (pin->GetOrAddBoolean("rad_vet","clean_conv_est_f",false)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_conv_est_f = true was removed "
                << "(only false)"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_verbose: 0 the CLEAN cycle lines, warnings and events; 1 also the per-solve
    // diagnostics
    if (pin->DoesParameterExist("rad_vet","clean_verbose")) {
      cln_verbose = pin->GetInteger("rad_vet","clean_verbose");
    }
    // clean_conv = budget | strict. budget: every row measured in the cell's own step
    // budget S (sum of its gross terms), tolerance clean_conv_kappa x the step's
    // time-error estimate, between clean_tol_nl and clean_conv_max, BiCGStab in that norm
    // with Eisenstat-Walker forcing; strict: E, T relative to themselves, clean_tol_nl
    {const std::string cv = pin->GetOrAddString("rad_vet","clean_conv","budget");
    if (cv != "budget" && cv != "strict") {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_conv = budget | strict" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    cln_budget = (cv == "budget");}
    cln_bkap = pin->GetOrAddReal("rad_vet","clean_conv_kappa",1.0);
    // a step dt_retry redoes is solved to the strict tolerance
    if (pin->GetOrAddBoolean("rad_vet","clean_retry_budget",false)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_retry_budget = true was removed "
                << "(only false)"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (!pin->GetOrAddBoolean("rad_vet","clean_halo_direct",true)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_halo_direct = false was removed "
                << "(only true)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_init = given | steady (read only when named): steady relaxes the radiation
    // before step 1
    if (pin->DoesParameterExist("rad_vet","clean_init")) {
      const std::string ci = pin->GetString("rad_vet","clean_init");
      if (ci != "given" && ci != "steady") {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "<rad_vet>/clean_init = given | steady" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      cln_steady = (ci == "steady");
    }
    if (cln_steady) {
      // problem/he_ic_equil: its relax keys give the defaults
      const bool heq = pin->DoesParameterExist("problem","he_ic_equil") &&
                       pin->GetBoolean("problem","he_ic_equil");
      auto eqk = [&](const char *k, const Real d) {
        return (heq && pin->DoesParameterExist("problem", k)) ?
               pin->GetReal("problem", k) : d;
      };
      cln_steady_n = pin->GetOrAddInteger("rad_vet","clean_init_n",
          static_cast<int>(eqk("he_ic_equil_relax_n", heq ? 40.0 : 50.0)));
      cln_steady_dtfac = pin->GetOrAddReal("rad_vet","clean_init_dtfac",
                                           eqk("he_ic_equil_dtfac", 1.0e6));
      cln_steady_tol = pin->GetOrAddReal("rad_vet","clean_init_tol",1.0e-8);
    }
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/scheme = '" << sc << "': only scheme = clean (with "
      << "transport = implicit; a restart file without the key ran the old solver)"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }
  }
  // geometry guard: uniform Cartesian, a spherical-polar wedge or full sphere, or the
  // cubed sphere (checks in SphericalS1Check, VetLatInit); no relativistic coordinates
  sph_geom = false;
  cs_geom = false;
  {
    Mesh *pm_ = ppack->pmesh;
    std::string why;
    const bool sp_ = pm_->use_spherical_polar && !pm_->use_cubed_sphere;
    const bool cs_ = pm_->use_cubed_sphere && !pm_->use_spherical_polar;
    if (pm_->use_spherical_polar && !sp_) {why += " mesh/use_spherical_polar";}
    // full sphere (use_polar_boundary) and theta stretch allowed on sp: a polar face is a
    // physical face of area 0 (reflect is exact there)
    const bool full_ok = sp_;
    if (pm_->use_polar_boundary) {
      if (!full_ok) {
        why += " mesh/use_polar_boundary (the full sphere needs transport = implicit on "
               "spherical polar)";
      }
    } else if (sp_ && (pm_->mesh_size.x2min <= 1.0e-10 ||
                       pm_->mesh_size.x2max >= M_PI - 1.0e-10)) {
      why += " a theta range that reaches a pole without mesh/use_polar_boundary (use "
             "ix2_bc = ox2_bc = polar with use_polar_boundary, or a wedge clear of the "
             "poles)";
    }
    if (pm_->use_grid_stretch_r && !sp_ && !cs_) {why += " mesh/use_grid_stretch_r";}
    if (pm_->use_grid_stretch_r_poly && !sp_ && !cs_) {
      why += " mesh/use_grid_stretch_r_poly";
    }
    if (pm_->use_grid_stretch_theta && !full_ok) {
      why += " mesh/use_grid_stretch_theta (needs transport = implicit on spherical "
             "polar)";
    }
    if (ppack->pcoord != nullptr &&
        (ppack->pcoord->is_special_relativistic ||
         ppack->pcoord->is_general_relativistic ||
         ppack->pcoord->is_dynamical_relativistic)) {
      why += " relativistic coordinates";
    }
    if (!why.empty()) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_vet> (explicit and implicit transport) supports a "
        << "uniform Cartesian mesh and, for the implicit solve only (closure "
        << "eddington | m1 | minerbo | kershaw), a spherical-polar wedge or the full "
        << "sphere and (closure eddington only) the cubed sphere; this input sets:"
        << why << "." << std::endl
        << "The M1 kernels would run with Cartesian uniform-dx arithmetic on it."
        << std::endl;
      std::exit(EXIT_FAILURE);
    }
    sph_geom = sp_;
    if (cs_) {
      sph_geom = true;
      cs_geom = true;
    }
  }
  nstage = M1_NSTAGE;
  impl_cfl = -1.0;
  marshak_q = 0.5;
  ibc_x1min = M1_IBC_MARSHAK;
  ibc_x1max = M1_IBC_MARSHAK;
  iflux_x1min = 0.0;
  iflux_x1max = 0.0;
  iebath_x1min = 0.0;
  iebath_x1max = 0.0;
  impl_nstep = 0.0;
  impl_itsum = 0.0;
  impl_itmax = 0.0;
  impl_solver = M1_ISOLV_LINE_JACOBI;
  trans_on = false;
  trans_x3 = false;
  bicg_on = false;
  impl_prec_fb = -1;
  bcg_capped = false;
  bcg_r0rel = 0.0;
  ew_fprev = 0.0;
  ew_etaprev = 0.0;
  impl_pred = false;
  pred_ok = false;
  pred_dt = 0.0;
  impl_offdiag = M1_OD_LAGGED;
  od_now = M1_OD_LAGGED;
  impl_vimp = false;
  vimp_now = false;
  t2_enth_be = false;
  time_scheme = 0;
  t2_ok = false;
  t2_solve = 0;
  t2_fail = false;
  t2_dbg_fail = -1;
  t2_nstep = 0.0;
  t2_nbe = 0.0;
  t2_nfall = 0.0;
  t2_dtprev = 0.0;
  t2_vprev = false;
  t2_vext = true;
  t2_nclip = 0.0;
  pred2_ok = false;
  pred2_dt = 0.0;
  t2_lin_save = 0.0;
  t2_onep_s = -1.0;
  impl_tlim = M1_TLIM_NONE;
  // Anderson acceleration
  aa_nh = 0;
  aa_head = 0;
  aa_hasp = false;
  aa_fnp = -1.0;
  st_edges = true;
  impl_prec = 0;
  e_floor = pin->GetOrAddReal("rad_vet","e_floor",(FLT_MIN));

  // removed explicit face-flux options: only their old defaults are accepted
  if (pin->DoesParameterExist("rad_vet","scaled_prefactor") &&
      pin->GetReal("rad_vet","scaled_prefactor") != 20.0) {
    std::cout << "### FATAL ERROR: <rad_vet>/scaled_prefactor = "
              << pin->GetString("rad_vet","scaled_prefactor") << " was removed (only 20)"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // opacities per unit mass (code units); kappa_e and kappa_f default to kappa_p
  {std::string op = pin->GetOrAddString("rad_vet","opacity","const");
  if (op.compare("const") == 0) {
    opacity_type = M1_OPAC_CONST;
  } else if (op.compare("powerlaw") == 0) {
    opacity_type = M1_OPAC_POWERLAW;
  } else if (op.compare("user") == 0) {
    opacity_type = M1_OPAC_USER;
  } else if (op.compare("table") == 0) {
    opacity_type = M1_OPAC_TABLE;
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
  // opac_abs_kappa_s [cm^2/g]; applied at every table lookup (M1TableOpacities)
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
               opacity_type != M1_OPAC_USER && opacity_type != M1_OPAC_TABLE);
  // vet_scatter: electron scattering.
  //  (i) opacity split: table os = kappa_e, of = kappa_a = kappa_R - kappa_e (guarded,
  //      kappa_a >= vet_kappa_floor), kappa_e = vet_kappa_es (const) or
  //      sigma_T x_e/(mu m_u) (eos); other laws: kappa_f, kappa_s. rho kappa_e goes to
  //      opac(M1_OP_S); the extinction stays the total M1_OP_T; scattering exchanges no
  //      energy.
  //  (ii) closure source VetScatterMix: vet_scatter_form = ma | absorption,
  //      vet_scatter_j = relaxed | current.
  // Default on for a fresh run with ray data (vet_sc, or vet_col + vet_gd), opacity =
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
  bool vet_rays = false;      // closure = vet_sc, or vet_col + vet_gd = true
  if (pin->DoesParameterExist("rad_vet","closure")) {
    const std::string cl0 = pin->GetString("rad_vet","closure");
    vet_rays = (cl0.compare("vet_sc") == 0) ||
               (cl0.compare("vet_col") == 0 && pin->DoesParameterExist("rad_vet","vet_gd")
                && pin->GetBoolean("rad_vet","vet_gd"));
  }
  if (pin->DoesParameterExist("rad_vet","vet_scatter")) {
    vscat = pin->GetBoolean("rad_vet","vet_scatter");
  } else if (fresh && vet_rays && opacity_type == M1_OPAC_TABLE && gen_eos_tab) {
    vscat = pin->GetOrAddBoolean("rad_vet","vet_scatter",true);
  }
  if (vscat) {
    std::string sf = pin->GetOrAddString("rad_vet","vet_scatter_form","ma");
    if (sf.compare("ma") == 0) {
      vscat_form = 1;
    } else if (sf.compare("absorption") == 0) {
      vscat_form = 2;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_vet>/vet_scatter_form = '" << sf << "' not implemented "
        << "(ma | absorption)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  // vet_scatter_j: relaxed (J = E*, the end-of-step pair of B*) | current
    std::string sj = pin->GetOrAddString("rad_vet","vet_scatter_j","relaxed");
    if (sj.compare("current") == 0) {
      vscat_jrel = false;
    } else if (sj.compare("relaxed") == 0) {
      vscat_jrel = true;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_vet>/vet_scatter_j = '" << sj << "' not implemented "
        << "(current | relaxed)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  // ma + current: kappa_P/kappa_R amplifies the lag between B* and J = E^n, the Picard
  // loop does not converge
    if (vscat_form == 1 && !vscat_jrel) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_vet>/vet_scatter_form = ma with vet_scatter_j = current "
        << "is not allowed: the kappa_P/kappa_R factor of the ma form (~250-520 in an "
        << "iron bump) amplifies the lag between the relaxed end-of-step B* and the "
        << "start-of-step J = E^n, and the implicit Picard loop does not converge "
        << "(no convergence). Use vet_scatter_j = relaxed (default) or "
        << "vet_scatter_form = absorption." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::string km = "const";
    if (opacity_type == M1_OPAC_TABLE) {
      // eos where a tabulated general EOS supplies x_e
      km = pin->GetOrAddString("rad_vet","vet_scatter_kappa_e",
                               (fresh && gen_eos_tab) ? "eos" : "const");
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
      std::cout << "rad_vet: vet_scatter on, form " << sf << ", J (relaxed source) " << sj
                << (vscat_jrel ? " (default)" : "");
      if (opacity_type == M1_OPAC_TABLE) {
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
  // advect_split: advective enthalpy-flux split of the E equation
  advect_split = pin->GetOrAddBoolean("rad_vet","advect_split",false);
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
                   opacity_type == M1_OPAC_USER ||
                   opacity_type == M1_OPAC_TABLE)) {
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
    // floorcap_lz (read only when named): the velocity cap scales the momentum without
    // L_z only (phi_hat: spherical polar e3, Cartesian (-y, x)/R, the cubed sphere
    // GnomonicPhiHat; the part along it is kept, also where it alone exceeds
    // floorcap_v), so that L_z is exact
    if (pin->DoesParameterExist("rad_vet","floorcap_lz")) {
      fcap_lz = pin->GetBoolean("rad_vet","floorcap_lz");
    }
    if (pin->DoesParameterExist("rad_vet","floorcap_report")) {
      fcap_report = pin->GetInteger("rad_vet","floorcap_report");
    }
    // floorcap_mesh_vol: counters weighted by the mesh cell volume, |dp| and L_z in the
    // frame of y (CleanVol, CleanLzArm)
    {
      if (!nw && !pin->DoesParameterExist("rad_vet","floorcap_mesh_vol")) {
        recipe::Log("<rad_vet>/floorcap_mesh_vol = false (a restart file without the "
                    "key: the old default)");
      }
      fcap_mvol = pin->GetOrAddBoolean("rad_vet","floorcap_mesh_vol",nw);
    }
    if (fcap_rho > 0.0 && (!coupling || !(arad > 0.0))) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/floorcap_rho needs coupling and arad > 0"
                << std::endl;
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
    force_ref = M1_FREF_NONE;
  } else if (fr.compare("wb_arad") == 0) {
    force_ref = M1_FREF_WB_ARAD;
  } else if (fr.compare("clean") == 0) {
    force_ref = M1_FREF_CLEAN;
  } else if (fr.compare("wb_arad_md") == 0) {
    force_ref = M1_FREF_WB_MD;
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
  const bool fwb = (force_ref == M1_FREF_WB_ARAD) || (force_ref == M1_FREF_WB_MD);
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
  // time2_fref_wacc: the exact reference-work payment also for the backward-Euler steps
  // of a hesdirk2 run
  fref_wacc_t2 = pin->GetOrAddBoolean("rad_vet", "time2_fref_wacc", true);
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
  if (force_ref != M1_FREF_NONE && !(coupling && gas_feedback)) {
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

  // closure
  vet_sc = false;
  vet_full = false;
  vet_mbs = nullptr;
  vet_nguard = 0.0;
  vet_ncell = 0.0;
  vet_guard_max = 0.0;
  vet_d11_min = 1.0;
  vet_time = 0.0;
  vet_itime = 0.0;
  vet_ncall = 0.0;
  tau_closure = false;
  tau_ready = false;
  vet_col = false;
  vcol_sph = false;
  vcol_axis_flux = false;
  vcol_nc = vcol_np = vcol_nmu = vcol_every = vcol_nray = 0;
  vcol_built = false;
  vcol_time = vcol_ncall = vcol_nskip = 0.0;
  vcol_sq = false;
  vcol_qmin = vcol_qmax = 0.0;
  vcol_lc = 1;
  {std::string cl = pin->GetOrAddString("rad_vet","closure","m1");
  chi_kind = M1_CHI_LEVERMORE;
  if (cl.compare("m1") == 0) {
    eddington = false;
  } else if (cl.compare("minerbo") == 0) {
    eddington = false;
    chi_kind = M1_CHI_MINERBO;
  } else if (cl.compare("kershaw") == 0) {
    eddington = false;
    chi_kind = M1_CHI_KERSHAW;
  } else if (cl.compare("eddington") == 0) {
    eddington = true;
  } else if (cl.compare("vet_sc") == 0) {
    // the implicit multi-D solve reads (chi, n) from the SC formal solution
    // (rad_vet_vet.cpp); other uses of chi keep Levermore's
    eddington = false;
    vet_sc = true;
    VetRotInit(pin);   // the vet_gd rotation keys and defaults
    // vet_gd_rotate_every has no vet_sc clock: inert
    if (pin->DoesParameterExist("rad_vet","vet_gd_rotate_every") &&
        pin->GetInteger("rad_vet","vet_gd_rotate_every") > 0 &&
        global_variable::my_rank == 0) {
      std::cout << "### NOTE: <rad_vet>/vet_gd_rotate_every has no effect with closure = "
                << "vet_sc (the set does not rotate); name vet_gd_rotate_dt" << std::endl;
    }
    RotAvgParse(pin);  // vet_rotate_avg (shared key)
    // clean_rot_ring named: the rotating set on the clean ring (vet_gd_rotate_xfade)
    if (vsc_rot && cln_ring) {
      vgd_xf = pin->GetOrAddBoolean("rad_vet","vet_gd_rotate_xfade",false) && VgdRotOn();
      vgd_rot_slope = 0;
    }
    // vet_source_noesrc (rad_vet_vet.cpp)
    vsrc_noes = pin->GetOrAddBoolean("rad_vet","vet_source_noesrc",fresh);
  } else if (cl.compare("tau") == 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/closure = tau was removed; use closure = vet_col with "
      << "vet_gd = true (spherical polar, cubed sphere) or vet_sc (Cartesian)"
      << std::endl;
    std::exit(EXIT_FAILURE);
  } else if (cl.compare("vet_col") == 0) {
    // per-column 1-D formal solution through tau_ten (rad_vet_vetcol.cpp)
    eddington = false;
    tau_closure = true;
    vet_col = true;
    vcol_nc = pin->GetOrAddInteger("rad_vet","vet_col_ncore",8);
    vcol_np = pin->GetOrAddInteger("rad_vet","vet_col_nsub",1);
    vcol_nmu = pin->GetOrAddInteger("rad_vet","vet_col_nmu",4);
    // vet_every (alias vet_col_every)
    vcol_every = pin->GetOrAddInteger("rad_vet",
                                      M1KeyName(pin, "vet_every", {"vet_col_every"}), 1);
    // vet_col_fk_min: lower clamp of f_K (restart files without it: 1/3)
    vcol_fkmin = pin->GetOrAddReal("rad_vet","vet_col_fk_min",
                                   global_variable::restart_run ? (1.0/3.0) : 0.0);
    if (!(vcol_fkmin >= 0.0 && vcol_fkmin <= 1.0/3.0 + 1.0e-12)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/vet_col_fk_min must be in [0, 1/3]"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // vet_col_reflect_top (sp): with implicit_bc_x1max = reflect the formal solution
    // mirrors the outgoing intensity at the top face (rad_vet_vetcol.cpp)
    if (sph_geom) {
      vcol_rtop = pin->GetOrAddBoolean("rad_vet","vet_col_reflect_top",
                                       !global_variable::restart_run);
    }
    // vet_col_surface_face (sp): acts with vet_col_surface_q and implicit_marshak_face =
    // linear (rad_vet.hpp)
    if (sph_geom) {
      vcol_sqf = pin->GetOrAddBoolean("rad_vet","vet_col_surface_face",
                                      !global_variable::restart_run);
    }
    // vet_col_surface_q: outer Marshak q from the formal solution
    vcol_sq = pin->GetOrAddBoolean("rad_vet","vet_col_surface_q",
                                   !global_variable::restart_run);
    vcol_qmin = pin->GetOrAddReal("rad_vet","vet_col_surface_qmin",1.0e-3);
    vcol_qmax = pin->GetOrAddReal("rad_vet","vet_col_surface_qmax",1.0);
    {
      // vet_source (alias vet_col_source), read by vet_sc too: relaxed_esc under the new
      // default set (recipe::NewDefaults), else relaxed
      const std::string vk = M1KeyName(pin, "vet_source", {"vet_col_source"});
      if (!pin->DoesParameterExist("rad_vet", vk)) {
        const bool nd = recipe::NewDefaults(pin);
        pin->GetOrAddString("rad_vet", vk, nd ? "relaxed_esc" : "relaxed");
        recipe::Log(std::string("<rad_vet>/vet_source = ") + (nd ? "relaxed_esc (default)"
            : (recipe::GetMode(pin) == recipe::Mode::none
               ? "relaxed (recipe none: the old default)"
               : "relaxed (restart input has no key: the old default)")));
      }
      std::string vs = pin->GetString("rad_vet", vk);
    if (vs.compare("relaxed") == 0) {
      vcol_srelax = true;
    } else if (vs.compare("relaxed_esc") == 0) {
      // relaxed with the cell's photon escape
      vcol_srelax = true;
      vcol_sesc = true;
    } else if (vs.compare("gas") == 0) {
      vcol_srelax = false;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_vet>/vet_col_source = '" << vs << "' not implemented "
        << "(relaxed | relaxed_esc | gas)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    }
    if (pin->DoesParameterExist("rad_vet","vet_source_noesrc")) {
      vsrc_noes = pin->GetBoolean("rad_vet","vet_source_noesrc");
    }
    // vet_bc_diff = face | cell: inner diffusion intensity of the formal solutions at the
    // face (S + mu(-dS/dtau) with S of the face, 1.5 S_a - 0.5 S_b) or from the first
    // cell
    {
      const char *bdd = global_variable::restart_run ? "cell" : "face";
      std::string bd = pin->GetOrAddString("rad_vet","vet_bc_diff",bdd);
      if (bd.compare("face") == 0) {
        vet_bc_face = true;
      } else if (bd.compare("cell") == 0) {
        vet_bc_face = false;
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<rad_vet>/vet_bc_diff must be cell or face"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    // vet_seg = mean | split | emit | jump: mean integrates the SC segment with
    // chi = (chi_u + chi_0)/2 and S linear in tau; split integrates its two parts at the
    // face with a face S limited toward S_u at an opaque -> thin interface
    // (rad_vet_vetseg.hpp); emit limits only where S_u > S_0; jump adds the
    // opaque -> thin jump term
    {
      const char *sgd = global_variable::restart_run ? "mean" : "jump";
      std::string sg = pin->GetOrAddString("rad_vet","vet_seg",sgd);
      if (sg.compare("split") == 0) {
        vet_seg_split = true;
        vet_seg_emit = false;
      } else if (sg.compare("emit") == 0) {
        vet_seg_split = true;
        vet_seg_emit = true;
      } else if (sg.compare("jump") == 0) {
        vet_seg_split = true;
        vet_seg_emit = true;
        vet_seg_jump = true;
        vet_seg_plateau = pin->GetOrAddBoolean("rad_vet", "vet_seg_plateau", false);
      } else if (sg.compare("mean") == 0) {
        vet_seg_split = false;
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<rad_vet>/vet_seg must be mean, split, emit or jump"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    std::string ax = pin->GetOrAddString("rad_vet","vet_col_axis","radial");
    if (ax.compare("flux") == 0) {
      vcol_axis_flux = true;
    } else if (ax.compare("radial") != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_vet>/vet_col_axis = '" << ax << "' not implemented "
        << "(radial | flux)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // vet_col_lat (rad_vet_vetcol.cpp): every key read only when named
    if (pin->DoesParameterExist("rad_vet","vet_col_lat")) {
      vlat_on = pin->GetBoolean("rad_vet","vet_col_lat");
    }
    // Athena++ x2rat / x3rat (not AthenaK keys) noted for VetLatInit
    for (const char *rk : {"x2rat", "x3rat"}) {
      if (pin->DoesParameterExist("mesh", rk) && pin->GetReal("mesh", rk) != 1.0) {
        vlat_rat_keys += std::string(" mesh/") + rk + " = " + pin->GetString("mesh", rk);
      }
    }
    // vet_gd (rad_vet_vetgd.cpp): rides on the vet_col_lat interface
    if (pin->DoesParameterExist("rad_vet","vet_gd")) {
      vgd_on = pin->GetBoolean("rad_vet","vet_gd");
      if (vgd_on) {vlat_on = true;}
    }
    // default set of a fresh start that names no set key: the geodesic vet_gd_geo = 3
    // (92 directions); 2-D axisymmetric with a phi extent below 2 pi: the staggered GL
    // 9 x nphi, nphi nearest 12 with 2 pi/nphi dividing the extent (fatal if none). A
    // restart naming none keeps HEALPix nside 3.
    if (vgd_on && !pin->DoesParameterExist("rad_vet","vet_gd_gl_nmu") &&
        !pin->DoesParameterExist("rad_vet","vet_gd_ls") &&
        !pin->DoesParameterExist("rad_vet","vet_gd_geo")) {
      const auto &ms = ppack->pmesh->mesh_size;
      const double per = ms.x3max - ms.x3min;
      const bool ax2 = ppack->pmesh->use_spherical_polar && !ppack->pmesh->three_d &&
          std::abs(per - 2.0*M_PI) >= 1.0e-8;
      const bool fresh = recipe::NewDefaults(pin) && !global_variable::restart_run;
      int np = 0;
      for (int q = 4; q <= 64 && fresh && ax2; ++q) {
        const double m = per*q/(2.0*M_PI);
        if (std::round(m) >= 1.0 && std::fabs(m - std::round(m)) < 1.0e-9 &&
            (np == 0 || std::abs(q - 12) < std::abs(np - 12))) {np = q;}
      }
      if (fresh && ax2 && np == 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "<rad_vet>/vet_gd: 2-D axisymmetric, no 2 pi/nphi (nphi 4..64) "
          << "divides the phi extent " << per << ": name the set (vet_gd_gl_nmu and "
          << "vet_gd_gl_nphi)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (fresh && !ax2) {
        pin->GetOrAddInteger("rad_vet","vet_gd_geo",3);
      }
      if (fresh && ax2) {
        pin->GetOrAddInteger("rad_vet","vet_gd_gl_nmu",9);
        pin->GetOrAddInteger("rad_vet","vet_gd_gl_nphi",np);
        pin->GetOrAddBoolean("rad_vet","vet_gd_gl_stagger",true);
      }
      recipe::Log(!fresh ? std::string("<rad_vet>/vet_gd: HEALPix nside 3 (") +
          (global_variable::restart_run ? "restart input has no key)" : "recipe none)")
          : (ax2 ? "<rad_vet>/vet_gd_gl_nmu = 9, vet_gd_gl_nphi = " + std::to_string(np) +
                   ", stagger (default on 2-D axisymmetric sp: exact phi walls)"
                 : std::string("<rad_vet>/vet_gd_geo = 3 (default)")));
    }
    // a fresh vet_gd run whose set rotates and names no vet_rotate_avg averages 4 rotated
    // sets (with the time clock: the cross-fade)
    if (vgd_on && recipe::NewDefaults(pin) && !global_variable::restart_run &&
        !pin->DoesParameterExist("rad_vet","vet_rotate_avg")) {
      const bool rdt = !pin->DoesParameterExist("rad_vet","vet_gd_rotate_dt") ||
          pin->GetString("rad_vet","vet_gd_rotate_dt") == "auto" ||
          pin->GetReal("rad_vet","vet_gd_rotate_dt") > 0.0;
      const bool rev = pin->DoesParameterExist("rad_vet","vet_gd_rotate_every") &&
          pin->GetInteger("rad_vet","vet_gd_rotate_every") > 0;
      if (rdt || rev) {
        pin->GetOrAddInteger("rad_vet","vet_rotate_avg",4);
        if (rdt && !rev) {
          pin->GetOrAddBoolean("rad_vet","vet_gd_rotate_xfade",true);
        }
        recipe::Log("<rad_vet>/vet_rotate_avg = 4 (default with a rotating vet_gd set)");
      }
    }
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_gl_nmu")) {
      vgd_gl_nmu = pin->GetInteger("rad_vet","vet_gd_gl_nmu");
    }
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_gl_nphi")) {
      vgd_gl_nphi = pin->GetInteger("rad_vet","vet_gd_gl_nphi");
    }
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_gl_stagger")) {
      vgd_gl_stag = pin->GetBoolean("rad_vet","vet_gd_gl_stagger");
    }
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_ls")) {
      vgd_ls = pin->GetInteger("rad_vet","vet_gd_ls");
    }
    // vet_gd_geo = n: geodesic set (vgd_ls = -n)
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_geo")) {
      const int ng_ = pin->GetInteger("rad_vet","vet_gd_geo");
      if (ng_ > 0) {vgd_ls = -ng_;}
    }
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_rotate_every")) {
      vgd_rot = pin->GetInteger("rad_vet","vet_gd_rotate_every");
    }
    // rotation keys (VetRotInit): vet_gd_rotate_dt = <dt_rot> | auto (20 x the first
    // step's dt, recorded) on a physical-time clock; default auto on a fresh vet_gd run
    // naming neither _dt nor _every and no vet_gd_beam. Sub-keys _t0, _align, _stretch,
    // _slope.
    VetRotInit(pin);   // shared with vet_sc (rad_vet_vetgd.cpp)
    // clean sp / cs, fresh, rotating on the time clock: vet_rotate_avg = 4 with the
    // cross-fade (recorded)
    if ((sph_geom || cs_geom) && !global_variable::restart_run &&
        VgdRotOn() &&
        (vgd_rot_dt > 0.0 || vgd_rot_auto)) {
      pin->GetOrAddInteger("rad_vet","vet_rotate_avg",4);
      pin->GetOrAddBoolean("rad_vet","vet_gd_rotate_xfade",true);
    }
    // a restart with an averaged set and no cross-fade key takes the cross-fade
    if (sph_geom && global_variable::restart_run &&
        VgdRotOn() && pin->DoesParameterExist("rad_vet","vet_rotate_avg") &&
        pin->GetInteger("rad_vet","vet_rotate_avg") > 1 &&
        !pin->DoesParameterExist("rad_vet","vet_gd_rotate_xfade")) {
      pin->GetOrAddBoolean("rad_vet","vet_gd_rotate_xfade",true);
    }
    RotAvgParse(pin);   // also after vet_sc's VetRotInit
    // scheme = clean: no slope correction of a rotating set (set before the first
    // VetGdRotStep); clean_rot_ring = jump: a block's first step is backward Euler
    if (vgd_on && VgdRotOn()) vgd_rot_slope = 0;
    if (cln_ring == 2) vgd_rot_slope = 2;
    // vet_gd_rotate_xfade (read only when named)
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_rotate_xfade")) {
      // a fixed set has nothing to cross-fade: inert there
      vgd_xf = pin->GetBoolean("rad_vet","vet_gd_rotate_xfade") && VgdRotOn();
      // vet_gd_rotate_xfade_m: 0 only
      if (vgd_xf && pin->DoesParameterExist("rad_vet","vet_gd_rotate_xfade_m") &&
          pin->GetInteger("rad_vet","vet_gd_rotate_xfade_m") != 0) {
        std::cout << "### FATAL ERROR <rad_vet>/vet_gd_rotate_xfade_m = "
                  << pin->GetInteger("rad_vet","vet_gd_rotate_xfade_m")
                  << " was removed (only 0)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      // fused two-set sweep (VetGdXffOk decides where it applies)
      if (vgd_xf) {vgd_xff = true;}
      if (vgd_xf && global_variable::my_rank == 0) {
        std::cout << "<rad_vet> vet_gd_rotate_xfade = true: the closure is "
                  << "cross-faded from the current to the next ray set over each "
                  << "rotation block (two sweeps per build)" << std::endl;
      }
    }
    // the cubed sphere with clean_rays = sc: vet_gd's own sweep (intensities to cln_fvi)
    cln_cssc = pmy_pack->pmesh->use_cubed_sphere && vgd_on;
    // clean_w_face_min: false only
    if (pin->GetOrAddBoolean("rad_vet","clean_w_face_min",false)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_w_face_min = true was removed "
                << "(only false)"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // vet_rotate_avg_ema_wb keeps the ring's closure of a static state, but the cs face
    // half ranges (cln_fh) are the two swept sets'
    if (cln_cssc && vra_ema_wb && global_variable::my_rank == 0) {
      std::cout << "### NOTE <rad_vet>/vet_rotate_avg_ema_wb on the cubed sphere with "
                << "clean_rays = sc: the face half ranges follow the two swept sets, so "
                << "a static state is not kept to round-off" << std::endl;
    }
    // clean_rot_ring unnamed: live where the vet_gd rays rotate in time with a cross-fade
    // (decided before a restart file's channels are counted)
    if (!cln_ring_named && vgd_on && vra_on && g_vra_adv &&
        vgd_xf && (vgd_rot_dt > 0.0 || vgd_rot_auto) && vgd_rot_align &&
        vra_k <= M1_RA_NF) {
      cln_ring = 1;
    }
    // vet_bc_bath: the top shell's inward rays start from the incident bath
    // implicit_ebath_x1max of a Marshak outer x1
    if (vgd_on) {
      vgd_tbath = pin->GetOrAddBoolean("rad_vet", "vet_bc_bath",
                                       !global_variable::restart_run);
    }
    // vet_bc_inner / vet_bc_reflect shared with vet_sc (VetBcInit; alias vet_gd_inner).
    // vet_interp = auto | log | corner: log bilinear in ln chi, ln S; corner one segment
    // per stencil corner; auto corner where the chi contrast exceeds kVetCornerRatio
    if (vgd_on) {
      const std::string vi = pin->GetOrAddString("rad_vet","vet_interp",
                                 global_variable::restart_run ? "log" : "auto");
      if (vi.compare("log") == 0) {
        vet_interp = 1;
      } else if (vi.compare("corner") == 0) {
        vet_interp = 2;
      } else if (vi.compare("auto") == 0) {
        vet_interp = 3;
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<rad_vet>/vet_interp = " << vi << " is not available "
                  << "on vet_gd (auto | log | corner)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    if (vgd_on) {VetBcInit(pin, true);}
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_allow_clamp")) {
      vgd_allow_clamp = pin->GetBoolean("rad_vet","vet_gd_allow_clamp");
    }
    // vet_gd_band_global: default true on a fresh sp full-sphere vet_gd start
    // (recipe::NewDefaults), recorded
    {
      Mesh *pm = ppack->pmesh;
      const bool fullsph = pm->use_spherical_polar && pm->use_polar_boundary &&
      pin->GetOrAddString("mesh", "ix2_bc", "") == "polar" &&
      pin->GetOrAddString("mesh", "ox2_bc", "") == "polar" &&
      pin->GetOrAddString("mesh", "ix3_bc", "") == "periodic" &&
      std::abs(pm->mesh_size.x3max - pm->mesh_size.x3min - 2.0*M_PI) < 1.0e-8;
      if (vgd_on && fullsph && !pin->DoesParameterExist("rad_vet","vet_gd_band_global")) {
        if (recipe::NewDefaults(pin)) {
          pin->GetOrAddBoolean("rad_vet","vet_gd_band_global",true);
          recipe::Log("<rad_vet>/vet_gd_band_global = true (default with vet_gd on a "
                      "fresh sp "
                      "full-sphere start)");
        } else {
          recipe::Log("<rad_vet>/vet_gd_band_global = false (" + std::string(
              recipe::GetMode(pin) == recipe::Mode::none ? "recipe none" :
              "restart input has no key") + ": the old default)");
        }
      }
    }
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_band_global")) {
      vgd_gband = pin->GetBoolean("rad_vet","vet_gd_band_global");
    }
    // vet_gd_band_ring (read only when named, VetGdInit): band intensities lent per shell
    // from a few slabs; unnamed: on where it applies except the sp raw band
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_band_ring")) {
      vgd_ring_named = true;
      vgd_ring_par = pin->GetBoolean("rad_vet","vet_gd_band_ring");
    }
    // vet_gd_halo_compact: default 2 on a fresh vet_gd start (recipe::NewDefaults),
    // recorded; VetGdInit switches it off where it does not apply
    if (vgd_on && !pin->DoesParameterExist("rad_vet","vet_gd_halo_compact")) {
      if (recipe::NewDefaults(pin)) {
        pin->GetOrAddInteger("rad_vet","vet_gd_halo_compact",2);
        recipe::Log("<rad_vet>/vet_gd_halo_compact = 2 (default with vet_gd on a fresh "
                    "start)");
      } else {
        recipe::Log("<rad_vet>/vet_gd_halo_compact = 0 (" + std::string(
            recipe::GetMode(pin) == recipe::Mode::none ? "recipe none" :
            "restart input has no key") + ": the old default)");
      }
    }
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_halo_compact")) {
      vgd_hcomp = pin->GetInteger("rad_vet","vet_gd_halo_compact");
    }
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_twin") &&
        pin->GetBoolean("rad_vet","vet_gd_twin")) {
      std::cout << "### FATAL ERROR <rad_vet>/vet_gd_twin = true: the twin is not in "
                   "this "
                << "build" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // vet_source_noesrc: fresh default true for a plain vet_col (no vet_gd)
    if (!vgd_on && fresh && !pin->DoesParameterExist("rad_vet","vet_source_noesrc")) {
      vsrc_noes = pin->GetOrAddBoolean("rad_vet","vet_source_noesrc",true);
    }
    if (pin->DoesParameterExist("rad_vet","vet_source_noesrc_h2be")) {
      vsrc_noes_h2be = pin->GetBoolean("rad_vet","vet_source_noesrc_h2be");
    }
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_wall_interp")) {
      vgd_wint = pin->GetBoolean("rad_vet","vet_gd_wall_interp");
    }
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_wall_sym")) {
      vgd_wsym = pin->GetBoolean("rad_vet","vet_gd_wall_sym");
    }
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_theta_seam_sym")) {
      vgd_wst = pin->GetBoolean("rad_vet","vet_gd_theta_seam_sym");
    }
    // vet_gd_theta_wall: periodic | mirror (specular) | vacuum (no incoming intensity at
    // a physical theta wall) | extrap (zero gradient)
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_theta_wall")) {
      const std::string tw = pin->GetString("rad_vet","vet_gd_theta_wall");
      if (tw.compare("mirror") == 0) {
        vgd_thmir = true;
      } else if (tw.compare("vacuum") == 0) {
        vgd_thvac = 1;
      } else if (tw.compare("extrap") == 0) {
        vgd_thvac = 2;
      } else if (tw.compare("periodic") != 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<rad_vet>/vet_gd_theta_wall must be periodic, mirror,"
                  << " vacuum or extrap" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    } else if (vgd_on && pmy_pack->pmesh->use_spherical_polar &&
               !pmy_pack->pmesh->use_polar_boundary) {
      // unnamed: a reflecting theta wall mirrors, a periodic one wraps
      vgd_thmir = (pmy_pack->pmesh->mesh_bcs[BoundaryFace::inner_x2] ==
                   BoundaryFlag::reflect);
    }
    // exact halo lists where they apply (VetGdInit switches them off elsewhere)
    vgd_hx_on = true;
    // vet_gd_halo_ipc: default on CUDA without a theta wall (VetGdSetup keeps it on the
    // raw 3-D sweep only; off-node partners take MPI)
#if defined(KOKKOS_ENABLE_CUDA)
    vgd_ipc_on = vgd_on && !vgd_thmir && vgd_thvac == 0;
#endif
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_halo_ipc")) {
      vgd_ipc_on = pin->GetBoolean("rad_vet","vet_gd_halo_ipc");
    }
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_replace")) {
      vgd_replace = pin->GetBoolean("rad_vet","vet_gd_replace");
    }
    // vet_gd_tangential: full tensor by default; recipe none on cs: D_rr only
    if (vgd_on && pmy_pack->pmesh->use_cubed_sphere &&
        recipe::GetMode(pin) == recipe::Mode::none) {vgd_tan = false;}
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_tangential")) {
      vgd_tan = pin->GetBoolean("rad_vet","vet_gd_tangential");
    }
    if (vlat_on) {
      // vet_col_lat_cut_thin (VetLatCsFinish)
      vlat_cut_thin = pin->GetOrAddReal("rad_vet","vet_col_lat_cut_thin",
                                        global_variable::restart_run ? 0.0 : 1.0);
      if (pin->DoesParameterExist("rad_vet","vet_col_lat_offdiag")) {
        std::string so = pin->GetString("rad_vet","vet_col_lat_offdiag");
        if (so.compare("none") == 0 || so.compare("false") == 0) {
          vlat_odm = 0;
        } else if (so.compare("lagged") == 0 || so.compare("true") == 0) {
          vlat_odm = 1;
        } else if (so.compare("operator") == 0) {
          vlat_odm = 2;
        } else {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<rad_vet>/vet_col_lat_offdiag = '" << so
                    << "' (none | lagged | operator)" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
      // vet_col_lat_posfix = drop | lag (rad_vet.hpp vlat_posfix_lag)
      if (pin->DoesParameterExist("rad_vet","vet_col_lat_posfix")) {
        std::string sp = pin->GetString("rad_vet","vet_col_lat_posfix");
        if (sp.compare("lag") == 0) {
          vlat_posfix_lag = true;
        } else if (sp.compare("drop") != 0) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<rad_vet>/vet_col_lat_posfix = '" << sp
                    << "' (drop | lag)" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
    }
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/closure = '" << cl << "' not implemented "
      << "(m1 | minerbo | kershaw | eddington | vet_sc | tau | vet_col)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  }

  // the M1 face states read two cells beyond the face
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  if (indcs.ng < 2) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet> requires at least 2 ghost zones, but <mesh>/nghost="
      << indcs.ng << std::endl;
    std::exit(EXIT_FAILURE);
  }

  dtnew = (FLT_MAX);
  dt_sub = 0.0;
  nsub = 1;

  // allocate arrays
  int nmb = std::max((ppack->nmb_thispack), (ppack->pmesh->nmb_maxperrank));
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  Kokkos::realloc(u0, nmb, M1_NVAR, ncells3, ncells2, ncells1);
  Kokkos::realloc(u1, nmb, M1_NVAR, ncells3, ncells2, ncells1);
  // face fluxes of the explicit stage tasks (never in the implicit task list): dummies
  Kokkos::realloc(uflx.x1f, nmb, M1_NVAR, 1, 1, 1);
  Kokkos::realloc(uflx.x2f, nmb, M1_NVAR, 1, 1, 1);
  Kokkos::realloc(uflx.x3f, nmb, M1_NVAR, 1, 1, 1);
  // rho*kappa in the ghost cells too; vet_scatter: a 4th component, rho kappa_e
  Kokkos::realloc(opac, nmb, vscat ? (M1_NOPAC + 1) : M1_NOPAC, ncells3, ncells2,
                  ncells1);
  Kokkos::deep_copy(opac, 0.0);
  if (coupling) {
    Kokkos::realloc(ugas1, nmb, 4, ncells3, ncells2, ncells1);
  }
  // force_reference: a valid View to capture; the pgen replaces it (SetForceReference)
  if (force_ref != M1_FREF_NONE) {
    Kokkos::realloc(arad_ref, nmb, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(arad_ref, 0.0);
  }
  // MHD: |B|^2/2 of the current (b0) and U^n (b1) fields, every cell incl. ghosts
  if (fl_mhd) {
    Kokkos::realloc(emag0, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(emag1, nmb, ncells3, ncells2, ncells1);
  } else {
    Kokkos::realloc(emag0, 1, 1, 1, 1);
    Kokkos::realloc(emag1, 1, 1, 1, 1);
  }
  Kokkos::realloc(cnt, M1_NCNT);
  for (int n=0; n<M1_NCNT; ++n) {cnt.h_view(n) = 0.0;}
  cnt.modify_host();
  cnt.sync_device();

  if (ppack->pmesh->multilevel) {
    int nccells1 = indcs.cnx1 + 2*(indcs.ng);
    int nccells2 = (indcs.cnx2 > 1)? (indcs.cnx2 + 2*(indcs.ng)) : 1;
    int nccells3 = (indcs.cnx3 > 1)? (indcs.cnx3 + 2*(indcs.ng)) : 1;
    Kokkos::realloc(coarse_u0, nmb, M1_NVAR, nccells3, nccells2, nccells1);
  }

  // the implicit solver's parameters, checks and arrays
  ImplicitInit(pin);
  // clean sc rays with a reflecting top: the mirror reads the last build
  cln_scmem = vgd_on && sph_geom && ibc_x1max == M1_IBC_REFLECT;
  // the ray cut (CleanRayCut): sc rays on 3-D sp, one MeshBlock along x1, no reflecting
  // top
  const Mesh *pm = ppack->pmesh;
  const bool ok = vgd_on && sph_geom &&
                  pm->three_d && !pm->multilevel && pm->mesh_indcs.nx1 == indcs.nx1 &&
                  !cln_scmem;
  cln_rcut_on = ok;
  if (!ok && global_variable::my_rank == 0) {
    std::cout << "<rad_vet> clean_ray_cut: not applicable here (staggered sc rays on "
              << "3-D spherical polar or the cubed sphere, one MeshBlock along x1 only)"
              << std::endl;
  }
  if (cln_rcut_on) {
    cln_rcut_r = pin->GetOrAddReal("rad_vet","clean_ray_cut_r",0.01);
    cln_rcut_mg = pin->GetOrAddReal("rad_vet","clean_ray_cut_margin",10.0);
    cln_rcut_tau = pin->GetOrAddReal("rad_vet","clean_ray_cut_taumin",0.0);
    cln_rcut_wait = pin->GetOrAddInteger("rad_vet","clean_ray_cut_wait",10);
    // the cut and its count in a restart file (CleanRayCut)
    cln_rcut = pin->GetOrAddInteger("rad_vet","clean_ray_cut_now",-1);
    cln_rcut_n = pin->GetOrAddInteger("rad_vet","clean_ray_cut_count",0);
    pin_rcut_ = pin;
  }
  // polar_ring_avg: false only
  if (pin->GetOrAddBoolean("rad_vet", "polar_ring_avg", false)) {
    std::cout << "### FATAL ERROR <rad_vet>/polar_ring_avg = true was removed (only "
                 "false)"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if (sph_geom) {SphericalS1Check(pin);}
  if (cs_geom) {
    pin_cs_ = pin;
    CubedS1Init();
  }

  // (5) boundary buffers
  pbval_u = new MeshBoundaryValuesCC(ppack, pin, false);
  // <mesh>/cs_seam_blockend = extend also on E, F (bvals_cc.cpp's along-seam resample)
  pbval_u->cs_blockend_ext_ = CsSeamBlockendExtend(ppack, pin);
  pbval_u->InitializeBuffers(M1_NVAR);
  // (F2, F3) is the tangential pair; E and F1 (radial) are scalars
  pbval_u->SetVectorPairs(M1_NVAR, {{M1_F2, M1_F3}});

  if (global_variable::my_rank == 0) {
    std::cout << "<rad_vet>: c=" << c_light << " chat/c=" << chat_over_c
              << " advect_split=" << (advect_split ? "true" : "false")
              // the closure as named in the input
              << " closure=" << (vet_col ? "vet_col" : vet_sc ? "vet_sc" :
                                 tau_closure ? "tau" : eddington ? "eddington" :
                                 (chi_kind == M1_CHI_MINERBO) ? "minerbo" :
                                 (chi_kind == M1_CHI_KERSHAW) ? "kershaw" : "m1")
              << std::endl;
    std::cout << "         kappa_p=" << kappa_p << " kappa_e=" << kappa_e
              << " kappa_f=" << kappa_f << " kappa_s=" << kappa_s
              << " coupling=" << (coupling ? "true" : "false")
              << " gas_feedback=" << (gas_feedback ? "true" : "false")
              << " arad=" << arad << std::endl;
    std::cout << "         opacity=" << (opacity_type == M1_OPAC_TABLE ? "table" : "")
              << " force_reference="
              << (force_ref == M1_FREF_WB_ARAD ? "wb_arad" :
                  (force_ref == M1_FREF_CLEAN ? "clean" :
                  (force_ref == M1_FREF_WB_MD ? "wb_arad_md" : "none")))
              << " temp_unit_kelvin=" << otab.tunit
              << " rho_unit_cgs=" << otab.dunit << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::SetOpacityTables
//! \brief hand the module the two opacity tables (opacity = table)

void RadiationVet::SetOpacityTables(const DvceArray2D<Real> &kr,
                                    const DvceArray2D<Real> &kp,
                                    const DvceArray1D<Real> &lT,
                                    const DvceArray1D<Real> &lD,
                                    const int nT, const int nD) {
  otab.kr = kr;
  otab.kp = kp;
  otab.lT = lT;
  otab.lD = lD;
  otab.nT = nT;
  otab.nD = nD;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::SetForceReference
//! \brief the reference acceleration of force_reference = wb_arad; resolves
//! force_reference_work = auto

void RadiationVet::SetForceReference(const DvceArray4D<Real> &a) {
  arad_ref = a;
  // force_reference_work = auto: split only if the pgen gives the reference work
  if (fref_wsplit_auto) {
    fref_wsplit = fref_wsplit_ok;
    fref_wsplit_auto = false;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> force_reference_work = auto -> "
                << (fref_wsplit ? "split" : "full") << std::endl;
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::SetX1maxBathColumns
//! \brief per-column incident bath E_bath(m,k,j) of the outer x1 Marshak face; also
//! allocates ibadv_x1max (lab-frame advective energy flux through that face)

void RadiationVet::SetX1maxBathColumns(const DvceArray3D<Real> &eb) {
  iebath_hi_col = eb;
  iebath_hi_col_on = true;
  Kokkos::realloc(ibadv_x1max, eb.extent_int(0), eb.extent_int(1), eb.extent_int(2));
  Kokkos::deep_copy(ibadv_x1max, 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::SetX2FaceFlux
//! \brief the imposed x2-face fluxes (m,k,i) of implicit_bc_x2min|max = flux

void RadiationVet::SetX2FaceFlux(const DvceArray3D<Real> &flo,
                                 const DvceArray3D<Real> &fhi) {
  ifx2_lo = flo;
  ifx2_hi = fhi;
  ifx2_on = true;
}

//----------------------------------------------------------------------------------------
//! \fn RadiationVet::FrefWaccOn, FrefWaccGam0, FrefWacc
//! \brief force_reference_work = split: the accumulator of the reference work the gas
//! receives in the hydro stages (fref_wacc)

bool RadiationVet::FrefWaccOn() const {
  // time2_fref_wacc: also the backward-Euler steps of a hesdirk2 run (first step, a
  // redo; not a helper solve)
  const bool t2be = fref_wacc_t2 && (time_scheme == M1_TIME_HESDIRK2) &&
                    (!t2_ok || FrefWaccStage()) && !t2_inprobe;
  return fref_wsplit && ((time_scheme == M1_TIME_BE) || t2be);
}

bool RadiationVet::FrefWaccStage() const {
  return t2_ok;
}

Real RadiationVet::FrefWaccGam0(Real bdt, Real dt_mesh) const {
  if (FrefWaccStage()) return 0.0;   // each stage solve pays its hydro stage's work
  for (int s = 0; s < fref_hnst; ++s) {
    if (std::abs(bdt - fref_hbeta[s]*dt_mesh) <= 1.0e-12*fref_hbeta[s]*dt_mesh) {
      return fref_hgam0[s];
    }
  }
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
            << "<rad_vet>/force_reference_work = split with time_scheme = be: no hydro "
            << "stage has beta dt = " << bdt << " (dt = " << dt_mesh << "; integrators "
            << "rk1, rk2, rk3 supported)"
            << std::endl;
  std::exit(EXIT_FAILURE);
  return 0.0;
}

DvceArray4D<Real> RadiationVet::FrefWacc() {
  const int n0 = arad_ref.extent_int(0), n1 = arad_ref.extent_int(1);
  const int n2 = arad_ref.extent_int(2), n3 = arad_ref.extent_int(3);
  if (fref_wacc.extent_int(0) != n0 || fref_wacc.extent_int(1) != n1 ||
      fref_wacc.extent_int(2) != n2 || fref_wacc.extent_int(3) != n3) {
    Kokkos::realloc(fref_wacc, n0, n1, n2, n3);
    Kokkos::deep_copy(fref_wacc, 0.0);
  }
  return fref_wacc;
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationVet::RssMode
//! \brief the closure tensor <hydro>/rad_signal_speed reads (hydro_newdt.cpp)

int RadiationVet::RssRstMode() const {return RssMode() + 16;}

int RadiationVet::RssMode() const {
  // clean: its own Eddington tensor cln_s, once CleanInit ran
  if (cln_s.extent(0) > 0 && cln_init) return 5;
  if (tau_closure && tau_ten.extent(0) > 0) return 2;
  if (vet_sc && vet_cell.extent(0) > 0) return vet_full ? 4 : 3;
  if (!eddington && !tau_closure && !vet_sc) return 1;
  return 0;
}

int RadiationVet::RssRstNch(int mode) const {
  // extended restart state (mode + 16): also the other opacity components and the tau
  // tensor's lateral slots (read before a step's first rebuild)
  if (mode >= 16) {
    const int b = mode - 16;
    return RssRstNch(b) + (opac.extent_int(1) - 1)
           + ((b == 2 && tau_ten.extent_int(1) > 4) ? (tau_ten.extent_int(1) - 4) : 0);
  }
  if (mode == 2 || mode == 3) return 5;
  if (mode == 4) return 4;
  if (mode == 5) return 7;   // opac_T and the six components of cln_s's D
  return 1;
}

void RadiationVet::RssSlot(int mode, int n, DvceArray5D<Real> *&a, int &c) {
  if (mode >= 16) {   // base slots, then opac (all but T), then tau_ten 4..
    const int b = mode - 16, nb = RssRstNch(b), nop = opac.extent_int(1) - 1;
    if (n < nb) {RssSlot(b, n, a, c); return;}
    const int k = n - nb;
    if (k < nop) {a = &opac; c = (k < M1_OP_T) ? k : (k + 1); return;}
    a = &tau_ten;
    c = 4 + (k - nop);
    return;
  }
  if (n == 0) {a = &opac; c = M1_OP_T; return;}
  if (mode == 2) {a = &tau_ten; c = n - 1; return;}
  if (mode == 5) {a = &cln_s; c = CLN_D11 + n - 1; return;}
  a = &vet_cell;
  c = (mode == 3) ? (M1_VET_CHI + n - 1) : (M1_VET_D11 + n - 1);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RssRstPack
//! \brief signal-speed inputs of MeshBlocks m0 .. m0+nmb-1, channel by channel, into
//! dst (nmb, nch, k, j, i)

void RadiationVet::RssRstPack(DvceArray5D<Real> &dst, int nmb, int mode, int m0) {
  const int nch = RssRstNch(mode);
  for (int n=0; n<nch; ++n) {
    DvceArray5D<Real> *a = nullptr;
    int c = 0;
    RssSlot(mode, n, a, c);
    Kokkos::deep_copy(Kokkos::subview(dst, Kokkos::ALL, n, Kokkos::ALL, Kokkos::ALL,
                                      Kokkos::ALL),
                      Kokkos::subview(*a, std::make_pair(m0, m0 + nmb), c, Kokkos::ALL,
                                      Kokkos::ALL, Kokkos::ALL));
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RssPrime
//! \brief before the first dt with <hydro>/rad_signal_speed: put back the signal-speed
//! inputs a restart file carries, else fill opac from the current state (the closure
//! tensor starts from what the module has)

void RadiationVet::RssPrime(bool restart) {
  if (rss_stage_mode >= 16) {
    // extended record; the tensor slots only into the same closure
    const int mode = rss_stage_mode, b = mode - 16, nb = RssRstNch(b);
    const bool tens = (b == 2 && tau_closure) || (b == 5 && RssMode() == 5) ||
                      ((b == 3 || b == 4) && vet_sc && vet_cell.extent(0) > 0 &&
                       (vet_full == (b == 4)));
    if (b == 2 && tau_closure && !tau_ready) {VetColInit();}
    const int nmb = static_cast<int>(rss_stage.extent(0));
    const int nall = static_cast<int>(rss_stage.extent(1));
    const int nop = opac.extent_int(1) - 1;
    for (int n=0; n<nall; ++n) {
      const bool isop = (n == 0) || (n >= nb && n < nb + nop);
      if (!isop && !tens) continue;
      DvceArray5D<Real> *a = nullptr;
      int c = 0;
      RssSlot(mode, n, a, c);
      if (a->extent_int(1) <= c) continue;
      Kokkos::deep_copy(Kokkos::subview(*a, std::make_pair(0,nmb), c, Kokkos::ALL,
                                        Kokkos::ALL, Kokkos::ALL),
                        Kokkos::subview(rss_stage, Kokkos::ALL, n, Kokkos::ALL,
                                        Kokkos::ALL, Kokkos::ALL));
    }
    rss_stage = DvceArray5D<Real>();
    rss_stage_mode = -1;
    return;
  }
  if (rss_stage_mode >= 0) {
    const int mode = rss_stage_mode;
    // the tensor is put back only into the same kind of closure; opac always
    int want = 1;
    if (mode == 2 && tau_closure) {
      if (!tau_ready) {VetColInit();}
      want = RssRstNch(mode);
    } else if (((mode == 3 || mode == 4) && vet_sc && vet_cell.extent(0) > 0 &&
                (vet_full == (mode == 4))) || (mode == 5 && RssMode() == 5)) {
      want = RssRstNch(mode);
    } else if (mode >= 2 && global_variable::my_rank == 0) {
      std::cout << "### WARNING: <rad_vet> the restart file's signal-speed closure "
                   "tensor "
                << "does not match this run's closure; only the opacity is restored"
                << std::endl;
    }
    const int nmb = static_cast<int>(rss_stage.extent(0));
    for (int n=0; n<want; ++n) {
      DvceArray5D<Real> *a = nullptr;
      int c = 0;
      RssSlot(mode, n, a, c);
      Kokkos::deep_copy(Kokkos::subview(*a, std::make_pair(0,nmb), c, Kokkos::ALL,
                                        Kokkos::ALL, Kokkos::ALL),
                        Kokkos::subview(rss_stage, Kokkos::ALL, n, Kokkos::ALL,
                                        Kokkos::ALL, Kokkos::ALL));
    }
    rss_stage = DvceArray5D<Real>();
    rss_stage_mode = -1;
    return;
  }
  if (restart && global_variable::my_rank == 0) {
    std::cout << "### WARNING: restart file has no <hydro|mhd>/rad_signal_speed inputs; "
              << "the opacity is rebuilt from the restored state and this restart is not "
              << "bitwise." << std::endl;
  }
  (void) Opacity(nullptr, 1);
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetBcInit
//! \brief x1 boundary intensity of the formal solutions, one rule for vet_sc and vet_gd.
//!   vet_bc_inner = auto | diff | emit (alias vet_gd_inner): rays entering at a Marshak
//!     or flux inner face; diff = the bottom cell's diffusion intensity (vet_bc_diff),
//!     emit = the face's own emission (bath implicit_ebath_x1min, or black sphere 4 F/c),
//!     auto = emit at Marshak, diff at flux. Restart files without it: vet_sc auto if
//!     vet_bc_bath, else diff; vet_gd diff.
//!   vet_bc_reflect = auto | mirror | open at a reflecting x1 end: mirror = specular
//!     (H = 0 at the face); open = diffusion at x1min, vacuum at x1max; auto = diffusion
//!     at x1min, mirror at x1max on vet_sc. vet_gd cannot mirror (its global set has no
//!     specular images). Restart files without it: open.

void RadiationVet::VetBcInit(ParameterInput *pin, const bool gd) {
  const bool rst = global_variable::restart_run;
  const std::string kin = gd ? M1KeyName(pin, "vet_bc_inner", {"vet_gd_inner"})
                             : std::string("vet_bc_inner");
  const bool bnamed = pin->DoesParameterExist("rad_vet", "vet_bc_bath");
  const bool bath = bnamed ? pin->GetBoolean("rad_vet", "vet_bc_bath") : !rst;
  std::string vi;
  if (pin->DoesParameterExist("rad_vet", kin)) {
    vi = pin->GetString("rad_vet", kin);
  } else if (!gd && (rst || (bnamed && !bath))) {
    vi = bath ? "auto" : "diff";     // vet_sc: the vet_bc_bath rule
  } else if (gd && rst) {
    vi = "diff";                     // vet_gd: vet_gd_inner unset
  } else {
    vi = pin->GetOrAddString("rad_vet", "vet_bc_inner", "auto");
  }
  if (vi.compare("diff") == 0) {
    vet_bc_inner = 0;
  } else if (vi.compare("emit") == 0) {
    vet_bc_inner = 1;
  } else if (vi.compare("auto") == 0) {
    vet_bc_inner = 2;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<rad_vet>/" << kin << " = '" << vi << "' is not a choice "
              << "(auto | diff | emit)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string vr;
  if (pin->DoesParameterExist("rad_vet", "vet_bc_reflect")) {
    vr = pin->GetString("rad_vet", "vet_bc_reflect");
  } else if (rst) {
    vr = "open";
  } else {
    vr = pin->GetOrAddString("rad_vet", "vet_bc_reflect", "auto");
  }
  if (vr.compare("auto") == 0) {
    // inner end: diffusion intensity (the isotropic part of H = 0); outer end: mirror
    vet_bc_mirror = false;
    vet_bc_mirror_hi = !gd;
  } else if (vr.compare("mirror") == 0) {
    if (gd) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/vet_bc_reflect = mirror is not available on "
                << "vet_gd (its direction set is global: the specular image of a "
                << "direction is not in it); use auto or open" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    vet_bc_mirror = true;
    vet_bc_mirror_hi = true;
  } else if (vr.compare("open") == 0) {
    vet_bc_mirror = false;
    vet_bc_mirror_hi = false;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<rad_vet>/vet_bc_reflect = '" << vr << "' is not a choice "
              << "(auto | mirror | open)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::VetBcEmit
//! \brief vet_bc_inner resolves to emit at the current inner x1 face

bool RadiationVet::VetBcEmit() const {
  if (ibc_x1min == M1_IBC_OUTFLOW) {return true;}   // clean: nothing enters
  if (vet_bc_inner == 1) {
    return (ibc_x1min == M1_IBC_MARSHAK) || (ibc_x1min == M1_IBC_FLUX);
  }
  if (vet_bc_inner == 2) {return (ibc_x1min == M1_IBC_MARSHAK);}
  return false;
}

//----------------------------------------------------------------------------------------
//! \fn RadiationVet::~RadiationVet
//! \brief waits for pending table work, frees IPC buffers, prints the reports

RadiationVet::~RadiationVet() {
  if (vgd_tab_state == 1) {vgd_tabf.wait();}   // the table prefetch
  if (vgd_tab_state == 3) {vgd_tex.fence();}
#if defined(KOKKOS_ENABLE_CUDA)
  if (vgd_ipc_state == 1) {     // vet_gd_halo_ipc
    (void)cudaDeviceSynchronize();
    for (void *q : vgd_ipc_open) {if (q != nullptr) {(void)cudaIpcCloseMemHandle(q);}}
    (void)cudaFree(vgd_ipc_rbuf);
    (void)cudaFree(vgd_ipc_sig);
    (void)cudaFree(vgd_ipc_err);
  }
#endif
#if defined(KOKKOS_ENABLE_HIP)
  if (vgd_ipc_state == 1) {     // vet_gd_halo_ipc (HIP)
    (void)hipDeviceSynchronize();
    for (void *q : vgd_ipc_open) {if (q != nullptr) {(void)hipIpcCloseMemHandle(q);}}
    (void)hipFree(vgd_ipc_rbuf);
    (void)hipFree(vgd_ipc_sig);
    (void)hipFree(vgd_ipc_err);
  }
#endif
  if (vet_col) {VetColReport();}
  ReportCounters();
  ImplicitReport();
  Time2Report();      // time_scheme = hesdirk2
  if (vet_sc) {VetReport();}
  VetFree();
  delete pbval_u;
  if (pbval_d3 != nullptr) {delete pbval_d3;}
  if (pbval_sg != nullptr) {delete pbval_sg;}
  if (pbval_ea != nullptr) {delete pbval_ea;}
  if (pbval_vs != nullptr) {delete pbval_vs;}
  if (pbval_hr != nullptr) {delete pbval_hr;}
  if (pbval_gd != nullptr) {delete pbval_gd;}
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ReportCounters
//! \brief one line at the end of the run: cost and health of the implicit energy
//! exchange solve (failures: cells that hit M1_MAXIT without reaching M1_RTOL)

void RadiationVet::ReportCounters() {
  if (!coupling) return;
  cnt.modify_device();
  cnt.sync_host();
  Real nsolve = cnt.h_view(M1_CNT_NSOLVE);
#if MPI_PARALLEL_ENABLED
  {
    Real loc[M1_NCNT], glob[M1_NCNT];
    for (int n=0; n<M1_NCNT; ++n) {loc[n] = cnt.h_view(n);}
    MPI_Reduce(loc, glob, M1_NCNT, MPI_ATHENA_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
    Real lmax = cnt.h_view(M1_CNT_ITMAX), gmax;
    MPI_Reduce(&lmax, &gmax, 1, MPI_ATHENA_REAL, MPI_MAX, 0, MPI_COMM_WORLD);
    for (int n=0; n<M1_NCNT; ++n) {cnt.h_view(n) = glob[n];}
    cnt.h_view(M1_CNT_ITMAX) = gmax;
    nsolve = cnt.h_view(M1_CNT_NSOLVE);
  }
#endif
  if (global_variable::my_rank != 0) return;
  Real mean = (nsolve > 0.0) ? (cnt.h_view(M1_CNT_ITSUM)/nsolve) : 0.0;
  std::cout << "<rad_vet> implicit energy solve: cells=" << nsolve
            << " iterations mean=" << mean << " max=" << cnt.h_view(M1_CNT_ITMAX)
            << " FAILURES=" << cnt.h_view(M1_CNT_NFAIL)
            << " (no bracket: " << cnt.h_view(M1_CNT_NBRAK) << ")" << std::endl;
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationVet::SetSubsteps
//! \brief the stage dt of a mesh step: one solve per step (no sub-cycling)

int RadiationVet::SetSubsteps(Real dt_mesh) {
  nsub = 1;
  dt_sub = dt_mesh;
  return nsub;
}


//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RotAvgParse
//! \brief <rad_vet>/vet_rotate_avg (+ _weight) and vet_rotate_avg_ema*, after the
//! closure's rotation keys

void RadiationVet::RotAvgParse(ParameterInput *pin) {
  // vet_rotate_avg = K (read only when named; < 2 off): closure averaged over the last
  // K rotated sets (rad_vet_rotavg.cpp)
  if (pin->DoesParameterExist("rad_vet","vet_rotate_avg")) {
    vra_k = pin->GetInteger("rad_vet","vet_rotate_avg");
    if (vra_k > M1_RA_NF) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/vet_rotate_avg must be <= " << M1_RA_NF
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    vra_on = (vra_k >= 2) && (vgd_on || vet_sc);
    // vet_rotate_avg_weight = j (sets weighted by their J) | equal
    if (vra_on && pin->DoesParameterExist("rad_vet","vet_rotate_avg_weight")) {
      const std::string w = pin->GetString("rad_vet","vet_rotate_avg_weight");
      if (w.compare("equal") == 0) {
        vra_eq = true;
      } else if (w.compare("j") != 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<rad_vet>/vet_rotate_avg_weight must be j or equal"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    const bool rot = vgd_on ? VgdRotOn() : vsc_rot;
    if (vra_on && !rot && global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_rotate_avg = " << vra_k << ": the direction set does "
                << "not rotate (vet_gd_rotate_dt / _every 0, or pinned): inert"
                << std::endl;
    }
    if (vra_on && rot && global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_rotate_avg = " << vra_k << ": the "
                << (vgd_on ? "vet_gd" : "vet_sc") << " closure is the mean over the "
                << "last " << vra_k << " rotated sets" << std::endl;
    }
  }
  // vet_rotate_avg_ema = N (default 128 for a fresh run whose set rotates; 0 off): the
  // exponential moving average over the rotated sets
  const bool nw = !global_variable::restart_run;
  const bool erot = nw && (vgd_on || vet_sc) && (vgd_on ? VgdRotOn() : vsc_rot);
  if (erot || pin->DoesParameterExist("rad_vet","vet_rotate_avg_ema")) {
    vra_ema_n = pin->GetOrAddReal("rad_vet","vet_rotate_avg_ema",128.0);
  }
  // vet_rotate_avg_ema_wb acts on the EMA ring only
  if (pin->DoesParameterExist("rad_vet","vet_rotate_avg_ema_wb") && vra_ema_n < 1.0) {
    std::cout << "### FATAL ERROR <rad_vet>/vet_rotate_avg_ema_wb needs "
              << "vet_rotate_avg_ema >= 1" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (vra_ema_n >= 1.0 && (vgd_on || vet_sc)) {
    vra_ema = true;
    vra_k = 3;   // slot 0 the EMA, slot 1 the change-test aux, slot 2 the Delta
    if (pin->DoesParameterExist("rad_vet","vet_rotate_avg_ema_change")) {
      vra_ema_eps = pin->GetReal("rad_vet","vet_rotate_avg_ema_change");
    }
    if (pin->DoesParameterExist("rad_vet","vet_rotate_avg_ema_change_j")) {
      vra_ema_chj = pin->GetBoolean("rad_vet","vet_rotate_avg_ema_change_j");
    }
    {
      // vet_rotate_avg_ema_change_form = linear | step
      const std::string fm = pin->GetOrAddString("rad_vet",
                                                 "vet_rotate_avg_ema_change_form",
                                                 "linear");
      if (fm.compare("step") == 0) {
        vra_ema_chstep = true;
      } else if (fm.compare("linear") != 0) {
        std::cout << "### FATAL ERROR <rad_vet>/vet_rotate_avg_ema_change_form = '" << fm
                  << "' is not a choice (linear | step)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    if (pin->DoesParameterExist("rad_vet","vet_rotate_avg_ema_wb")) {
      vra_ema_wb = pin->GetBoolean("rad_vet","vet_rotate_avg_ema_wb");
    }
    if (pin->DoesParameterExist("rad_vet","vet_rotate_avg_ema_dema")) {
      vra_ema_dema = pin->GetBoolean("rad_vet","vet_rotate_avg_ema_dema");
      if (vra_ema_dema) {vra_k = 4;}   // slot 2 = e2, slot 3 the Delta
    }
    if (pin->DoesParameterExist("rad_vet","vet_rotate_avg_ema_diag")) {
      vra_ema_diag = pin->GetInteger("rad_vet","vet_rotate_avg_ema_diag");
    }
    vra_on = true;
    if (pin->DoesParameterExist("rad_vet","vet_rotate_avg_weight")) {
      vra_eq = (pin->GetString("rad_vet","vet_rotate_avg_weight").compare("equal") == 0);
    }
    if (pin->DoesParameterExist("rad_vet","vet_rotate_avg_ema_tau")) {
      vra_ema_tau = pin->GetReal("rad_vet","vet_rotate_avg_ema_tau");
    }
    if (pin->DoesParameterExist("rad_vet","vet_rotate_avg_ema_lc")) {
      vra_ema_lc = pin->GetReal("rad_vet","vet_rotate_avg_ema_lc");
    }
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_rotate_avg_ema = " << vra_ema_n << ": "
                << "the closure is an exponential moving average over the rotated sets "
                << "(a = 1/N per block" << ((vra_ema_tau > 0.0) ? ", time window tau = "
                + std::to_string(vra_ema_tau) : "") << ((vra_ema_lc > 0.0) ?
                ", window = " + std::to_string(vra_ema_lc) + " x r/c" : "")
                << "; change test eps " << vra_ema_eps
                << ")" << (vra_ema_dema ? "; lag-corrected double EMA 2 e1 - e2" : "")
                << "; vet_rotate_avg is not used"
                << std::endl;
    }
  }
}

} // namespace radvet
