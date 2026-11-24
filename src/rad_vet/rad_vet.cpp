//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet.cpp
//! \brief constructor, parameter reading and allocation for the RadiationVet class

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
// <rad_vet>/vet_gd_twin_srw (defined in rad_vet_vetgd.cpp)
extern bool g_vra_adv;   // vet_rotate_advance (rad_vet_vetgd.cpp)
//----------------------------------------------------------------------------------------
// constructor, initializes data structures and parameters

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
  // ---- MILESTONE 3a hook: which transport scheme runs.  Everything the implicit one
  // needs is read and allocated by ImplicitInit at the end of this constructor.
  // (a fresh run's default implicit is set before, recipe::CleanDefaults; a restart file
  // without the key ran explicit)
  {std::string tr = pin->GetOrAddString("rad_vet","transport","explicit");
  if (tr.compare("explicit") == 0) {
    transport = M1_TRANSPORT_EXPLICIT;
  } else if (tr.compare("implicit_x1") == 0) {
    transport = M1_TRANSPORT_IMPLICIT_X1;
  } else if (tr.compare("implicit") == 0) {
    transport = M1_TRANSPORT_IMPLICIT;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/transport = '" << tr << "' is not a valid choice "
      << "(explicit | implicit_x1 | implicit)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  }
  // ---- <rad_vet>/scheme = clean: the clean stage solve (rad_vet_clean.cpp),
  // the default of a fresh implicit run (v20.4).  The production values below are the
  // defaults of a fresh clean run; a restart file without a key keeps its old default
  const bool nw = !global_variable::restart_run;
  {const bool scd = nw && (transport == M1_TRANSPORT_IMPLICIT);
  if (scd && !pin->DoesParameterExist("rad_vet","scheme")) {
    recipe::Log("<rad_vet>/scheme = clean (the default with transport = implicit)");
  }
  std::string sc = pin->GetOrAddString("rad_vet","scheme",scd ? "clean" : "default");
  if (sc.compare("clean") == 0) {
    clean = true;
    CleanPreflight(pin);
    // clean_w_fse (default on): the thin share's thinness is the old production's
    // (implicit_thin_switch = fse): g = max(|H|/J, eps C1/1.5), eps = min(kappa_P/chi,
    // 1), w = smoothstep((g - c0)/(c1 - c0)) with c0, c1 = 0.3, 0.6 unless named.
    // Upwind streaming holds where the field is anisotropic (|H|/J), not where it is
    // merely non-local (C1 of an unresolved opaque photosphere cell, I != S across it,
    // isotropic field); eps removes the non-locality of scattering
    cln_wfse = pin->GetOrAddBoolean("rad_vet","clean_w_fse",true);
    // clean_w_fse_tau: scattering isotropises only over the cell's optical depth tau_d =
    // chi V/A_f, so eps -> 1 - (1 - eps)(1 - exp(-tau_d)) (a thin or empty cell keeps C1)
    cln_wfst = cln_wfse && pin->GetOrAddBoolean("rad_vet","clean_w_fse_tau",true);
    cln_c0 = pin->GetOrAddReal("rad_vet","clean_w_c0",cln_wfse ? 0.3 : 0.1);
    cln_c1 = pin->GetOrAddReal("rad_vet","clean_w_c1",cln_wfse ? 0.6 : 0.5);
    // clean_face_hctau (default off): the bracket's rays' face flux F_y interpolated in
    // optical depth between the cell centres (old production implicit_thin_hcface = tau,
    // CleanStagBracket), the rays' face flux kept where the cells' H jump (beams)
    cln_hctau = pin->GetOrAddBoolean("rad_vet","clean_face_hctau",false);
    // clean_face_jse (default on): the face half ranges in the shape of the cells'
    // half-range fluxes at the face (old production implicit_thin_recon = jse)
    cln_jse = pin->GetOrAddBoolean("rad_vet","clean_face_jse",true);
    // clean_face_jse_cell (default off): its factor on the donor cell's h+- (old
    // M1ThinJs), not on the reconstructed face value
    cln_jsec = pin->GetOrAddBoolean("rad_vet","clean_face_jse_cell",false);
    // clean_face_hclat (default off): the bracket's F_y on the lateral faces of sp / cs is
    // the mean of the two cells' H (old production hc with implicit_thin_hcface_gd = x1)
    cln_hclat = pin->GetOrAddBoolean("rad_vet","clean_face_hclat",false);
    // clean_face_area (default off): the sp / cs streaming share in area form
    // (CleanStagArea) also under clean_stag_gflux: its h+- are the donor cell's rays
    cln_farea = pin->GetOrAddBoolean("rad_vet","clean_face_area",false);
    // clean_bc_area (default on): an open x1 end's half-range flux c h+- E of the end
    // cell's point intensities in area form: its net outward stream h+ + h- x A_c/A_f
    // (an isotropic field none: the bath balances; old: Marshak q = H_f/J)
    cln_bcarea = pin->GetOrAddBoolean("rad_vet","clean_bc_area",true);
    // clean_ray_source_relaxed (default off): the rays' source the vet_source /
    // vet_scatter_j form (relaxed_esc: the local exchange with escape, its relaxed pair)
    cln_rsrel = pin->GetOrAddBoolean("rad_vet","clean_ray_source_relaxed",false);
    cln_tol = pin->GetOrAddReal("rad_vet","clean_tol",1.0e-10);
    cln_cnorm = pin->GetOrAddReal("rad_vet","clean_init_cnorm",1.0e-10);
    cln_lcnorm = pin->GetOrAddReal("rad_vet","clean_lin_cnorm",1.0e-10);
    // clean_maxnl (200), clean_maxit (800) and clean_conv_max (1e-5) do not take the old
    // production's implicit_maxit, implicit_lin_maxit and implicit_tol_rel_max (a named
    // other value prints a NOTE).  clean_lin_pass (default on, the old production's
    // form): a linear solve that ends at clean_maxit leaves its iterate to the next pass,
    // without the stall test (CleanSolve); a restart file without the key keeps the form
    // a named implicit_lin_maxit set
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
    // clean_rays = sc (vet_gd's own sweep; the FV rays fv1 | fv | lld are removed)
    if (pin->DoesParameterExist("rad_vet","clean_rays") &&
        pin->GetString("rad_vet","clean_rays") != "sc") {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_vet>/clean_rays = sc (the FV rays fv1 | fv | lld were "
        << "removed)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_rot_ring = off | live | jump (read only when named; live wherever the vet_gd
    // rays rotate with a cross-fade, CleanInit): a rotating union swept as the current
    // set (+ the next one), the older ones carried as moment offsets (CleanRingMoments)
    // clean_rot_fuse (default on): the live ring's two sets
    // (b, b+1) of a build in one fused sweep (vet_gd_rotate_xfade_fuse; the raw sp band
    // and the cubed sphere's compact lists: one message per partner, VetGdGbHaloX /
    // VetGdCsHaloCX), the closure of each from its own intensities: bitwise the two
    // sweeps; a second intensity array
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
    // clean_fill_cs_remap (default off): the fill's h+- of a cubed-sphere seam neighbour
    // in our face families (CleanClosureFill)
    cln_fillcs = pin->GetOrAddBoolean("rad_vet","clean_fill_cs_remap",false);
    cln_pred = pin->GetOrAddBoolean("rad_vet","clean_predict",true);
    // implicit_predictor (step unless named none): the stage solve starts from the
    // predicted E (CleanFirstIterate), as the old production's step predictor (a
    // restart file's default "none", which ImplicitInit adds, is not named)
    cln_kpred = !(pin->DoesParameterExist("rad_vet","implicit_predictor") &&
                  !pin->IsParameterDefaulted("rad_vet","implicit_predictor") &&
                  pin->GetString("rad_vet","implicit_predictor") == "none");
    // clean_pred_stage (default false): sdirk2's stage-2 solve starts from E_old + h K1 of
    // this step's stage 1 (E^n + dt K1), not of the last step's stage 2
    cln_pst = pin->GetOrAddBoolean("rad_vet","clean_pred_stage",false);
    // clean_gas_secant = false (default) | true | stage2: the exchange's slopes from the
    // chord through the last pass's gas T (CleanRightHandSide), in every solve or in the
    // stage-2 solves only; max | max_stage2: the chord only where it is steeper than the
    // tangent
    {const std::string gs = pin->GetOrAddString("rad_vet","clean_gas_secant","false");
    cln_gsec = (gs == "true") ? 1 : ((gs == "stage2") ? 2 : ((gs == "max") ? 3 :
               ((gs == "max_stage2") ? 4 : 0)));}
    // clean_pred_xch (default false): sdirk2's stage-2 solve starts from E_old + h K, K the
    // stage-1 slope of this step where the cell's exchange is slow, the last step's
    // stage-2 slope where it is stiff (weights 1 : z, z = h (ch kE + 4 c kP a T^4/e))
    cln_pxc = pin->GetOrAddBoolean("rad_vet","clean_pred_xch",false);
    {std::string bf = pin->GetOrAddString("rad_vet","clean_bc","halfrange");
    cln_bdry_form = (bf.compare("marshak") == 0) ? 1 : 0;}
    // clean_layout = staggered, the only layout (the cell layout was removed; a restart
    // file without the key ran it)
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
    // (old production implicit_opac_lag = pred, implicit_opac_lag_passes: the opacities
    // of passes 0 .. N-1, then frozen: a change of the discretisation, order 2 for N >= 2)
    cln_optn = pin->GetOrAddInteger("rad_vet","clean_stag_opac_passes",3);
    // clean_stag_opac_drag (default on for a fresh run, USER 11-21, recorded; a restart
    // file without the key keeps off, echoed): a refresh rebuilds the clean_stag_vimp
    // drag share b and Q of the faces (CleanStagDrag) with its chi
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
    // clean_stag_fs_line = pcr | thomas: the F solve's x1 lines on a GPU, on every mesh
    // (the cubed sphere: with clean_stag_fschur); default pcr on a fresh run, USER 11-21,
    // recorded; unnamed in a restart file, the old keys clean_stag_fs_pcr (true) and, on
    // the cubed sphere, clean_stag_fs_pcr_cs (false: thomas, echoed)
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
    // (the face tie: a face two MeshBlocks share, or one with itself across a periodic
    // wrap, takes its owner's F, momentum exchange, energy-flux memory and (the
    // write-back) Phi, CleanStagFaceTie)
    // clean_stag_tie_gid (default on for a fresh run on a uniform mesh, USER 11-21: one
    // rule on every mesh; recorded, a restart file without it keeps the old default):
    // spherical polar, Cartesian, the tie's T, d and Phi by the cubed sphere's rule
    // (CleanStagCsTie): the lower gid's, a block with itself across a periodic wrap its
    // low face's (F stays the low face's, its storage)
    old("clean_stag_tie_gid", "false");
    cln_tgid = pin->GetOrAddBoolean("rad_vet","clean_stag_tie_gid",
                                    nw && !pmy_pack->pmesh->multilevel);
    if (cln_tgid && pmy_pack->pmesh->multilevel) {
      std::cout << "### FATAL ERROR: <rad_vet>/clean_stag_tie_gid needs a uniform mesh"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_stag_fsweeps N (read only when named, default 1): the staggered F solve's
    // line solves up to N times, the line ends taking the neighbours' F of the last sweep
    if (pin->DoesParameterExist("rad_vet","clean_stag_fsweeps")) {
      cln_fsweeps = std::max(pin->GetInteger("rad_vet","clean_stag_fsweeps"), 1);
    }
    // clean_stag_fschur (read only when named; default on where it applies): the F
    // solve's lines exact across MeshBlocks and periodic wraps (CleanStagFSchur), in
    // place of the sweeps; on a uniform mesh, and as a segment goes out in three cells
    // at each end of a line, MeshBlocks of >= 3 cells along it (>= 6 where it goes both
    // ways: more than two round a periodic line, two or more along one with ends);
    // elsewhere the line ends lagged, as the log says
    {Mesh *pm = pmy_pack->pmesh;
    // (the cubed sphere: a lateral line is a ring of four panels' MeshBlocks)
    const bool cs = pm->use_cubed_sphere;
    bool ok = cln_fsweeps == 1 && !pm->multilevel, lag = false;
    for (int d = 0; d < 3; ++d) {
      // (layout_invariant_mode = tile: as on MeshBlocks of a tile, CleanStagFSchurTile)
      const int n = (d > 0 && repro::tile_n > 0) ? repro::tile_n :
          ((d == 0) ? pm->mb_indcs.nx1 :
                      ((d == 1) ? pm->mb_indcs.nx2 : pm->mb_indcs.nx3));
      const int nb = ((cs && d > 0) ? 4 : 1)*((d == 0) ? pm->mesh_indcs.nx1 :
                     ((d == 1) ? pm->mesh_indcs.nx2 : pm->mesh_indcs.nx3))/n;
      const bool per = (cs && d > 0) || (pm->mesh_bcs[2*d] == BoundaryFlag::periodic);
      const int nlong = (cs && per && nb%2 == 0) ? 7 : 6;   // (cs: the kept segment)
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
    } else if (!ok && lag && !cs && cln_fsweeps == 1) {   // (no silent lagged ends)
      std::cout << "### FATAL ERROR: <rad_vet>/clean_stag_fschur cannot take this mesh "
                << "(nghost >= 3, MeshBlock size, uniform mesh): name clean_stag_fsweeps "
                << "(line ends iterated) or clean_stag_fschur = false" << std::endl;
      std::exit(EXIT_FAILURE);
    } else if (!ok && global_variable::my_rank == 0) {
      std::cout << "<rad_vet> clean_stag_fschur off: the F solve's line ends lagged "
                << "(clean_stag_fsweeps, mesh or MeshBlock size)" << std::endl;
    }
    // (layout_invariant_mode = tile: the join of a periodic line in whole MeshBlocks,
    // CleanStagFSchurTile)
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
    // (the cubed sphere without it iterates its seams: clean_stag_fsweeps default 64,
    // its lines by Thomas)
    if (cs && !cln_fschur && !pin->DoesParameterExist("rad_vet","clean_stag_fsweeps")) {
      cln_fsweeps = 64;
    }
    if (cs && !cln_fschur && cln_fspcr && global_variable::my_rank == 0) {
      std::cout << "<rad_vet> clean_stag_fs_line = pcr: thomas on the cubed sphere "
                << "without clean_stag_fschur" << std::endl;
    }}
    // the consistency term is the default on Cartesian meshes and, with the vet_gd
    // rays (sc), on spherical polar and the cubed sphere (well-balanced only with them)
    const bool curv = pmy_pack->pmesh->use_spherical_polar ||
                      pmy_pack->pmesh->use_cubed_sphere;
    const bool vg = pin->DoesParameterExist("rad_vet","vet_gd") &&
                    pin->GetBoolean("rad_vet","vet_gd");
    // (sp and cs, fresh: recorded so that a restart keeps it; a restart file without
    // the key ran the removed FV rays fv1: FATAL unless the input names clean_rays)
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
    // (with the split energy flux, where the rays are not opt-in)
    cln_split = true;
    cln_consist = pin->GetOrAddBoolean("rad_vet","clean_consistency", !curv || vg);
    cln_maxpic = pin->GetOrAddInteger("rad_vet","clean_maxpic",12);
    // the staggered scheme: the energy flux apart from the face momentum (split), the
    // face momentum of its own F solve (fsolve) in the Newton loop (fwork), with
    // clean_consistency the rays' kinetic momentum flux in it (fkin) and their face flux
    // as flux factors (gflux), the donor-cell share of the advective flux (upwind).  Read
    // only when named: clean_stag_vimp, the gas velocity implicit in the face momentum;
    // clean_stag_sbimp (default off), the Ma et al. beta terms linearised in E
    auto stag = [&](const char *key, const bool def) {
      return pin->DoesParameterExist("rad_vet", key) ? pin->GetBoolean("rad_vet", key) :
                                                      def;
    };
    cln_fsolve = cln_upw = cln_split;
    cln_fkin = cln_fwork = cln_gflux = cln_split && cln_consist;
    // clean_stag_freal (default on with the F solve and gflux): realisable by
    // construction, the comoving part of sb at most 2 |beta| hc kE E_k and the faces of a
    // restart from the old scheme's file at most c E of their donor
    cln_freal = stag("clean_stag_freal", cln_gflux && cln_fsolve);
    // clean_stag_upwind_central (default off): an upwind face's fixed point carries the
    // central mean of the cells' a E within the TVD bound (CleanStagAdvectCorr), the old
    // production's cell flux a E
    cln_ucen = stag("clean_stag_upwind_central", false);
    // clean_stag_central (default off): the rest that of the plain central mean, the
    // radiation work's drain on the E rows' diagonal and the gas's kick in the same
    // Patankar weight (CleanStagAdvectCorr, CleanRightHandSide, CleanStagMprkWeight)
    cln_ucp = stag("clean_stag_central", false);
    cln_ucen = cln_ucen || cln_ucp;
    // clean_stag_csvec_radial (true, the only form): the cubed sphere's cell vector exact
    // for a radial field (CleanStagCsVector); false (the old least-squares fit) is gone
    cln_csvr = stag("clean_stag_csvec_radial", true);
    if (!cln_csvr) {
      std::cout << "### FATAL ERROR: <rad_vet>/clean_stag_csvec_radial = false (the old "
                << "cell vector) is gone" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_stag_thinf (default on): the staggered bracket's ray correction in the face
    // half ranges, realisable (CleanStagBracket), not as the gated lagged flux DG
    cln_thinf = pin->GetOrAddBoolean("rad_vet","clean_stag_thinf",nw);
    // (fwb + fcone bound F, and keep free streaming; every mesh: the face F is the
    // momentum along the face's unit normal, the cubed sphere's too, CleanCsNormA)
    if (pin->GetOrAddBoolean("rad_vet","clean_stag_flim",false)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_stag_flim = true was removed (only "
                                "false)"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_stag_fwb (default on for a fresh run with the F solve and clean_consistency,
    // recorded; a restart file without the key: off): the F rows in deviation form about
    // the rays' own state, so that it is their exact steady state (CleanStagFSolve)
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
    // clean_stag_fwb_stage (default off, recorded): fwb's residual s of the stage-old
    // state, once per stage (CleanStagFSolve); off: of each solve's E, F_ray too, so that
    // the thin rows' F is the rays' flux factor times the E they leave (old prod's
    // c (h+ E_L + h- E_R) of E^{n+1}); E_old's carried the stage slopes' O(h K) into F
    // (HG: a flux dipole at the opaque photosphere, FDIP)
    cln_fwbst = cln_fwb && pin->GetOrAddBoolean("rad_vet","clean_stag_fwb_stage",false);
    // clean_stag_fwb_t2 (default on for a fresh run with fwb_stage, USER 11-21, recorded;
    // a restart file without the key: off): fwb_stage's s of the stage-time state
    // E_old + h K, K the last slope (CleanStagFSolve): second order in time (no-op
    // without fwb; FATAL with fwb but without fwb_stage)
    old("clean_stag_fwb_t2", "false");
    cln_fwbt2 = pin->GetOrAddBoolean("rad_vet","clean_stag_fwb_t2",nw && cln_fwbst);
    if (cln_fwbt2 && cln_fwb && !cln_fwbst) {
      std::cout << "### FATAL ERROR: <rad_vet>/clean_stag_fwb_t2 needs "
                << "clean_stag_fwb_stage" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    cln_fwbt2 = cln_fwbt2 && cln_fwbst;
    // clean_stag_fwb_line (default off), clean_stag_fwb_ring (default on for a fresh run
    // with fwb on a uniform 3-D spherical polar mesh with periodic x3, USER 11-21,
    // recorded; a restart file without the key: off): fwb's s of the x1 faces
    // less its line mean, of the phi faces (3-D spherical polar, periodic x3) less its
    // ring mean, in proportion to the streaming share: no net radial force on a line, no
    // torque on a ring, momentum and L_z exact without a ledger (CleanStagFwbProject)
    cln_fwbpj[0] = cln_fwb && stag("clean_stag_fwb_line", false);
    const bool ring_ok = pmy_pack->pmesh->use_spherical_polar && pmy_pack->pmesh->three_d &&
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
        pin->GetBoolean("rad_vet", "clean_stag_fwb_cart") && !pmy_pack->pmesh->multilevel &&
        !pmy_pack->pmesh->use_spherical_polar && !pmy_pack->pmesh->use_cubed_sphere) {
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
    // clean_stag_fcone_gas (default off): the gas takes the momentum fcone cuts, at the
    // face, as its share of the face's conservative term (the cut F's absorption not
    // withheld), so gas + radiation momentum and L_z are exact; E pays the kinetic energy
    cln_fcgas = cln_fcone && pin->GetOrAddBoolean("rad_vet","clean_stag_fcone_gas",false);
    // clean_stag_thinf_held (default off, the v20.7 value, USER 11-21): the thin flux
    // factor of a donor blended by theta = min(1, max(J, P)/(clean_ray_min E_old)), P its
    // largest ray outflow: theta g(rays) + (1 - theta) h (its own half range); a
    // ray-shadowed donor is no dam
    cln_thinh = cln_thinf &&
        pin->GetOrAddBoolean("rad_vet","clean_stag_thinf_held",false);
    // clean_stag_xdonor_freeze (default 4; 0 = off): from this pass of a solve on, the
    // xdonor cross flux C is held (CleanStagXDonor): its donor side and cap no longer
    // switch from pass to pass (HG C2D: a 1-cell donor/cap flip, period-2 passes at 2e-4)
    cln_xdfz = pin->GetOrAddInteger("rad_vet","clean_stag_xdonor_freeze",4);
    if (const int pd = pin->GetOrAddInteger("rad_vet","clean_pass_damp",0); pd != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_pass_damp = " << pd
                << " was removed (only 0)"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_stag_xdonor = lag (default) | p1 | off: the off-diagonal stress flux an
    // M-matrix term (CleanStagXDonor); clean_stag_xsrc (default on): also the face
    // sources past their explicit share; clean_stag_xsweep (default off): a monotone
    // line-Jacobi sweep after each linear solve
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
    cln_xsw = cln_xd && stag("clean_stag_xsweep", false);
    // clean_stag_kpos = off (default) | faces (= true) | cell | redo: the stage increment
    // in MPRK form (CleanStagKPos): faces, every slope face flux and the cell part
    // moved; cell, the cell part and the x1 boundary only, the face fluxes explicit;
    // redo, as cell, and where the solve gives X <= 0 the stage again with those
    // donors' faces moved (CleanStagKRedo, rounds, the last all faces)
    {
      const std::string kp = !pin->DoesParameterExist("rad_vet", "clean_stag_kpos") ?
                             "off" : pin->GetString("rad_vet", "clean_stag_kpos");
      if (kp != "off" && kp != "false" && kp != "true" && kp != "faces" &&
          kp != "cell" && kp != "redo" && kp != "delta") {
        std::cout << "### FATAL ERROR <rad_vet>/clean_stag_kpos = off | faces | cell | "
                     "redo"
                  << " | delta" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      // (delta: as faces, but each face's slope flux G splits into tau (al E + ar E) of
      // the stage's own coefficients, implicit, and the rest dG = G - tau (al E_s,L +
      // ar E_s,R) moved as q X_donor: near equilibrium dG -> 0, CleanStagKDelta)
      cln_kmode = (kp == "true" || kp == "faces") ? 1 : (kp == "cell" ? 2 :
                  (kp == "redo" ? 3 : (kp == "delta" ? 4 : 0)));
      cln_kpos = (cln_kmode > 0);
    }
    cln_kreal = stag("clean_stag_kreal", false);
    // clean_stag_kpat (default off): the stage increment's gas kick in Patankar form, the
    // start's E and gas eint >= 0, the momentum rest in the faces' F (CleanStagKPat)
    cln_kpat = stag("clean_stag_kpat", false);
    // clean_stag_eicarry (default off): the stage start's gas eint < 0 enters the gas row
    // as is, paid by its implicit exchange, not taken as 0 (CleanOldVector)
    cln_eicry = stag("clean_stag_eicarry", false);
    if (cln_kpat && cln_kreal) {
      std::cout << "### FATAL ERROR <rad_vet>/clean_stag_kpat with clean_stag_kreal is not "
                << "implemented" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    cln_colf = stag("clean_stag_colflux", false);
    // clean_cs_lz_abs = 2 (the cubed sphere, the only form): the gas phi momentum is
    // the faces' absorption as its r and theta (CleanStagMomentum), and the L_z it does
    // not take is booked as the radiation's (CleanLzBook), so that the total stays
    // exact; 0 | 1 (the torque row of the radiation's cell vector, no ledger) are gone
    if (pmy_pack->pmesh->use_cubed_sphere) cln_lza = 2;
    if (pin->DoesParameterExist("rad_vet", "clean_cs_lz_abs") &&
        pin->GetInteger("rad_vet", "clean_cs_lz_abs") != 2) {
      std::cout << "### FATAL ERROR: <rad_vet>/clean_cs_lz_abs = 0 | 1 (the old cubed-"
                << "sphere gas phi momentum) is gone; 2 is the only form" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    cln_lzbook = (cln_lza == 2) || (pmy_pack->pmesh->use_spherical_polar &&
                                    pmy_pack->pmesh->three_d);
    // clean_cs_lz_tel (default on for a fresh staggered cubed-sphere run, recorded):
    // its L_z by telescoping, the gas phi momentum of CleanCsTorque and the radiation's
    // L_z that of its face momenta, so that the ledger holds the radial boundaries'
    // torque alone (force_reference = clean withholds its phi row, CleanCsTorque)
    if (cln_lza == 2 || pin->DoesParameterExist("rad_vet", "clean_cs_lz_tel")) {
      cln_lzt = pin->GetOrAddBoolean("rad_vet","clean_cs_lz_tel",nw);
    }
    // clean_cs_lz_bnd_free: no radiation stress torque through the radial boundaries on
    // the cubed sphere (CleanCsTorque; with clean_cs_lz_tel, default on for a fresh run)
    if (cln_lzt) cln_lzbf = pin->GetOrAddBoolean("rad_vet","clean_cs_lz_bnd_free",nw);
    if (cln_lzt && cln_lza != 2) {
      std::cout << "### FATAL ERROR: <rad_vet>/clean_cs_lz_tel needs the cubed sphere "
                << "and clean_layout = staggered" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_cs_lat_wb (default off): a lateral face projects its cells' lateral stress
    // sums after the rotation of the cell's radial direction into the face plane, so
    // that a spherically symmetric P gives no lateral face force (CleanStagCsLat)
    cln_latwb = stag("clean_cs_lat_wb", false);
    if (cln_latwb && cln_lza != 2) {
      std::cout << "### FATAL ERROR: <rad_vet>/clean_cs_lat_wb needs the cubed sphere "
                << "and clean_layout = staggered" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (cln_lzbook) {
      const bool cs = (cln_lza == 2);
      cln_lzpin = pin;
      for (int k1 = 0; k1 < 2; ++k1) {   // (a restart file's non-finite value: NOTE, 0)
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
    cln_xpat = stag("clean_stag_xpat", false);
    // clean_xt_cons (default off): the first iterate's gas T with A's G = 1 - beta^2/2,
    // and after the first right-hand side T from the pass's gas row (its G and beta
    // terms: one Newton update at X = E_k), so the first pass's A is the gas row's;
    // clean_xt_norm (default off, study): the local exchange with E responding,
    // e(T) + A T^4/(1 + hc kE) = eint + hc kE E/(1 + hc kE)
    cln_xtcons = stag("clean_xt_cons", false);
    // clean_xt_exact = cons (default) | off | x2: the exchange's bookkeeping with one
    // A' per cell and pass (CleanRightHandSide): cons, the first iterate's T from the
    // pass's gas row (clean_xt_cons), the work the row and the gas's kinetic energy
    // differ by to E whole (eint takes the heating only), the closure printed (CLEAN
    // XT); x2 also every T root from the local 2x2 of the E and gas rows, E + e kept
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
    if (cln_kpat && cln_xte == 0 && !cln_kpos && !cln_xpat) {   // (CLW_KF..+2: the rest)
      std::cout << "### FATAL ERROR <rad_vet>/clean_stag_kpat needs clean_xt_exact" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (cln_ucp && (cln_kpos || cln_xpat || cln_xte != 1)) {   // (CLW_PD: the work's q)
      std::cout << "### FATAL ERROR <rad_vet>/clean_stag_central needs clean_xt_exact = "
                << "cons, without clean_stag_kpos / xpat" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (cln_kpos && cln_xte == 2) {   // (x2's local 2x2 has no kpos gas term)
      std::cout << "### FATAL ERROR <rad_vet>/clean_stag_kpos with clean_xt_exact = x2 "
                << "is not implemented" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    cln_xtnorm = stag("clean_xt_norm", false) || cln_xte == 2;
    // clean_stag_rayratio (default off): the bracket reads the rays only as ratios (their
    // face flux over the donor's ray J, the model flux of E), times the state's E
    cln_rayratio = pin->GetOrAddBoolean("rad_vet","clean_stag_rayratio",false) &&
                   cln_gflux;
    // clean_ring_jmin (default on): a clean_rot_ring J ratio J_s/J_b is measured only
    // where both sets' J are at least clean_ray_min E (a ratio of two ray-shadowed J is
    // no ray pattern); elsewhere the offset is invalid (left out of the mix)
    cln_ringjmin = pin->GetOrAddBoolean("rad_vet","clean_ring_jmin",true);
    // clean_rot_ring_fill_rst (default on): the ring's fill at a restart from a file
    // without the clean state (cln_rstold: an old scheme's, whose closure averaged a
    // full window) takes the older sets' offsets swept at the restart state, as the
    // steady start (CleanRingMoments)
    cln_rfrst = pin->GetOrAddBoolean("rad_vet","clean_rot_ring_fill_rst",true);
    // clean_rst_mshare (default on): at a restart from a file without the clean state
    // the M1 share's memory starts at its steady flux of the first solution
    // (CleanStagFSolve)
    cln_rmsh = pin->GetOrAddBoolean("rad_vet","clean_rst_mshare",true);
    // clean_source_convex (default on): the rays' vet_scatter source of the clean scheme
    // (the state's own, CleanRays) with the absorption share kappa_P/chi at most 1, as
    // the thermal source without vet_scatter: S = J + min(kappa_P/chi, 1) (B - J) between
    // J and B (kappa_P > chi, a cool envelope's Planck mean above its Rosseland
    // extinction, made S = J + 750 (B - J) < 0 there: floored to 1e-300, a source hole in
    // the rays)
    cln_srccvx = pin->GetOrAddBoolean("rad_vet","clean_source_convex",true);
    if (cln_freal && !(cln_gflux && cln_fsolve)) {
      std::cout << "### FATAL ERROR: <rad_vet>/clean_stag_freal needs clean_consistency "
                << "(the F solve and gflux)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    cln_vimp = stag("clean_stag_vimp", cln_split);
    // (the cubed sphere: only with clean_stag_vimp_cs, default on for a fresh run, USER
    // 11-21, recorded; a restart file without the key keeps off, echoed)
    if (cln_vimp && pmy_pack->pmesh->use_cubed_sphere) old("clean_stag_vimp_cs", "false");
    if (cln_vimp && pmy_pack->pmesh->use_cubed_sphere &&
        !pin->GetOrAddBoolean("rad_vet","clean_stag_vimp_cs",nw)) {
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_vet> staggered on the cubed sphere: clean_stag_vimp = false"
                  << std::endl;
      }
      cln_vimp = false;
    }
    cln_sbimp = stag("clean_stag_sbimp", false) && cln_split;
    // clean_gas_tfloor (default true): the gas row's T bounded below by the EOS
    // tfloor, and the gas eint written back from that root (E takes the rest)
    cln_gtf = stag("clean_gas_tfloor", true);
    // clean_wb_fct (default on for a fresh run with clean_gas_tfloor and clean_xt_exact,
    // USER 11-21, recorded; a restart file without the key: off, echoed): the
    // write-back's Zalesak limit on the written faces' energy flux, E > 0 by construction
    // (CleanWbFct); its budget is the gas row's of clean_gas_tfloor with
    // clean_xt_exact's work
    old("clean_wb_fct", "false");
    cln_wbfct = pin->GetOrAddBoolean("rad_vet","clean_wb_fct",
                                     nw && cln_gtf && cln_xte > 0);
    if (cln_wbfct && !(cln_gtf && cln_xte > 0)) {
      std::cout << "### FATAL ERROR <rad_vet>/clean_wb_fct needs clean_gas_tfloor and "
                << "clean_xt_exact = cons | x2" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_wb_floor_fail (read only when named; default false): a write-back whose E
    // floor takes energy the cell's own gas does not give (from the radial column or
    // created, M1_POS_FLR_MV + M1_POS_FLR_UN) is not admissible: the stage fails (the
    // backward-Euler redo), a backward-Euler step fails into dt_retry (reason
    // efloor_column); at the last cut the step is accepted, counted, CLEAN WARNING
    if (pin->DoesParameterExist("rad_vet","clean_wb_floor_fail")) {
      cln_wbff = pin->GetBoolean("rad_vet","clean_wb_floor_fail");
    }
    // clean_wb_exch_local (read only when named; default false; needs clean_gas_tfloor
    // and an EOS tfloor > 0): the write-back's cell exchange solved exactly from the
    // transported E (CleanWriteBack): E > 0 and the gas >= e(rho, tfloor) by
    // construction, energy exact; the face FCT limits the transported E
    if (pin->DoesParameterExist("rad_vet","clean_wb_exch_local")) {
      cln_wbxl = pin->GetBoolean("rad_vet","clean_wb_exch_local");
    }
    // clean_wb_floor_diag N (diagnostic, read only when named): the terms of the first
    // N cells of each write-back whose E the floor raises (CLEAN WBFDIAG)
    if (pin->DoesParameterExist("rad_vet","clean_wb_floor_diag")) {
      cln_wbfdg = pin->GetInteger("rad_vet","clean_wb_floor_diag");
    }
    // clean_gas_newton_switch (default W = 4 on a fresh run, USER 11-21, recorded; 0 = off;
    // W in [3, 64]), clean_gas_newton_switch_min (default 6 on a fresh run, else 20),
    // clean_stall_fac (0.5): the old production's implicit_gas_newton_switch, a stalled
    // solve's remaining passes on the bracketed gas root without the Newton slope; a
    // restart file without the keys keeps 0 / 20
    if (!nw && !pin->DoesParameterExist("rad_vet","clean_gas_newton_switch")) {
      recipe::Log("<rad_vet>/clean_gas_newton_switch = 0 (a restart file without the "
                  "key: "
                  "the old default)");
    }
    cln_gnsw = pin->GetOrAddInteger("rad_vet","clean_gas_newton_switch",nw ? 4 : 0);
    cln_gnswmin = pin->GetOrAddInteger("rad_vet","clean_gas_newton_switch_min",nw ? 6 : 20);
    cln_stfac = pin->GetOrAddReal("rad_vet","clean_stall_fac",0.5);
    if (cln_gnsw < 0 || cln_gnsw > 64 || (cln_gnsw > 0 && cln_gnsw < 3) ||
        cln_gnswmin < 1 || !(cln_stfac > 0.0 && cln_stfac <= 1.0)) {
      std::cout << "### FATAL ERROR <rad_vet>/clean_gas_newton_switch must be 0 or in "
                << "[3, 64], clean_gas_newton_switch_min >= 1, clean_stall_fac in (0, 1]"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // the preconditioner: the old production mg (CleanStagVCycle)
    {// clean_stag_mg_color = mirror (default) | parity: the mg's line colours,
    // mirror-symmetric or (j + k) parity (CleanStagLineSolve)
    const std::string mc = pin->GetOrAddString("rad_vet","clean_stag_mg_color",
                                               nw ? "mirror" : "parity");
    if (mc != "parity" && mc != "mirror") {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_stag_mg_color = parity | mirror"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    cln_stagmg = (mc == "mirror") ? 2 : 1;
    // (layout_invariant_mode = tile: the mg acts on fixed tiles of N x N (x2, x3) cells
    // of the mesh, N the mesh's layout_tile, as on MeshBlocks of that size (line colours,
    // couplings kept, coarse aggregates inside a tile), so that its action, and the
    // solve, are the same for every MeshBlock layout whose MeshBlocks are whole tiles
    // (CleanStagLineSolve; one MeshBlock along x1))
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
    // clean_pc_sweep = sym (default of a fresh run, recorded; a restart file without the
    // key keeps rb) | rb | rbr | rbrsym | wave: the mg's smoothing on each level, a
    // forward red-black sweep, + a second red pass (rbr), + a black-red post-sweep after
    // the coarse correction (sym); wave: symmetric Gauss-Seidel over the x1 lines in
    // wavefronts from the MeshBlock edges in and back out, before and after the coarse
    // correction (a ray-dark region's donor-cell transport across x2, x3)
    // (CleanStagVCycle)
    const std::string sw = pin->GetOrAddString("rad_vet","clean_pc_sweep",
                                               nw ? "sym" : "rb");
    cln_pcsweep = (sw == "rbr") + 2*(sw == "sym") + 3*(sw == "rbrsym") + 4*(sw == "wave");
    if (!cln_pcsweep && sw != "rb") {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_pc_sweep = rb | rbr | sym | rbrsym | "
                                "wave"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }}
    // clean_stag_stencil: the operator applied from its rows (CleanStencilRow), built
    // once per pass (default, as production's implicit_op_stencil; on the cubed sphere
    // the rows without the lateral stress, which the stored coefficients of its sums add,
    // CleanStagCsLateral); not with clean_stag_sbimp
    cln_stst = stag("clean_stag_stencil", !cln_sbimp);
    if (cln_stst && cln_sbimp) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_stag_stencil: not with clean_stag_sbimp"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_stag_lat_sums: spherical polar / Cartesian, the lateral stress of the faces
    // from the cells' sums Y (CleanStagCsLateral), the cubed sphere's form (spherical
    // polar phi faces: the torque form); with clean_stag_xdonor, whose cross flux reads
    // them, and no clean_stag_sbimp (the operator keeps the central rows).  Default on
    // for a fresh run with those, USER 11-21 uniform settings, recorded; a restart file
    // without the key keeps the 4-cell central differences (echoed)
    old("clean_stag_lat_sums", "false");
    cln_ysum = pin->GetOrAddBoolean("rad_vet","clean_stag_lat_sums",
                                    nw && cln_xd && !cln_sbimp) &&
               pmy_pack->pmesh->multi_d && !pmy_pack->pmesh->use_cubed_sphere;
    if (cln_ysum && (!cln_xd || cln_sbimp)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_stag_lat_sums: with clean_stag_xdonor = "
                << "lag | p1 and without clean_stag_sbimp" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_closure (default extrap1: the rays' source from converged states; predict's
    // explicit K1 state patterns the cs photosphere, STAGFIX5)
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
    // clean_skip_time2_hist (read only when named; default true): no arrays for the
    // implicit solve's predictors (ipred, ipred2) and the vet_col extrap2 builds (x2_tt,
    // x2_hr, x2_q), which the clean scheme neither reads nor writes; their restart
    // channels are zeros and carry no extrap2 history (x2_t = -1e300)
    cln_not2h = true;
    if (pin->DoesParameterExist("rad_vet","clean_skip_time2_hist")) {
      cln_not2h = pin->GetBoolean("rad_vet","clean_skip_time2_hist");
    }
    {const std::string nn = pin->GetOrAddString("rad_vet","clean_nl_norm","global");
    if (nn.compare("global") != 0 && nn.compare("cell") != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_vet>/clean_nl_norm = global | cell" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    cln_nl_cell = (nn.compare("cell") == 0);}
    // (the exchange root stops at its residual noise, CleanExchangeTemperature)
    cln_trtol = pin->GetOrAddReal("rad_vet","clean_troot_rtol",1.0e-13);
    if (pin->GetOrAddBoolean("rad_vet","clean_conv_est_f",false)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_conv_est_f = true was removed "
                << "(only false)"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_verbose (read only when named): 0 (default) the CLEAN cycle lines, warnings
    // and events; 1 also the per-solve diagnostics (CLEAN MCAP, OLDVEC, FCONE, FWB, LZ,
    // WBFCT, KPOS, KREAL, XPAT, the write-back closure and the ray cut)
    if (pin->DoesParameterExist("rad_vet","clean_verbose")) {
      cln_verbose = pin->GetInteger("rad_vet","clean_verbose");
    }
    // clean_conv = budget (default) | strict: budget measures every row in the
    // cell's own budget for the step, S = the sum of its gross terms (old value, fluxes
    // through the faces, exchange, work, beta terms), so that a change or residual is
    // the fraction of the energy the cell exchanges in the step.  Its tolerance is
    // clean_conv_kappa (1.0) x the step's time-error estimate in the same weights,
    // between clean_tol_nl (10 clean_tol) and clean_conv_max (1e-5), and BiCGStab in
    // that norm with Eisenstat-Walker forcing; strict: E, T relative to themselves,
    // clean_tol_nl
    {const std::string cv = pin->GetOrAddString("rad_vet","clean_conv","budget");
    if (cv != "budget" && cv != "strict") {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_conv = budget | strict" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    cln_budget = (cv == "budget");}
    cln_bkap = pin->GetOrAddReal("rad_vet","clean_conv_kappa",1.0);
    // clean_conv_max (read above): the tolerance's ceiling, default 1e-5.  The MeshBlock
    // layout moves the answer to about it (a solve may stop a pass earlier or later),
    // which a step dt_retry redoes amplifies (mpirej_retry, 1 vs 2 ranks: 1e-4 9.1e-2,
    // 1e-5 9.2e-6)
    // (a step dt_retry redoes is solved to the strict tolerance)
    if (pin->GetOrAddBoolean("rad_vet","clean_retry_budget",false)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_retry_budget = true was removed "
                << "(only false)"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (!pin->GetOrAddBoolean("rad_vet","clean_halo_direct",true)) cln_hdl.state = -1;
    // env CLEAN_POISON_GHOST=1 (gate only): the ghost layers clean_halo_direct does not
    // send are NaN
    // (1: NaN; 2: +1e30; 3: -1e30, for readers that compare or take min / max)
#if VET_DEBUG
    if (std::getenv("CLEAN_POISON_GHOST") != nullptr) {
      const int pm = std::atoi(std::getenv("CLEAN_POISON_GHOST"));
      cln_hdl.poison = true;
      cln_hdl.pval = (pm == 2) ? 1.0e30 : ((pm == 3) ? -1.0e30 : cln_hdl.pval);
    }
#endif
    // clean_init = given | steady (read only when named: an input without it keeps a
    // byte-identical parameter dump).  A start out of radiative equilibrium takes steady:
    // the rays, a steady solution, run ahead of the given field's front.
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
      // (problem/he_ic_equil starts: the old module's relax, he_ic_equil_relax_n steps
      // of he_ic_equil_dtfac x dt, defaults 40 and 1e6; at 1e3-1e6 dt the gas cooled to
      // a photosphere with inward F, RSG L/L_in -0.13, and the re-march inverted rho)
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
    if (transport != M1_TRANSPORT_IMPLICIT) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_vet>/scheme = clean needs transport = implicit"
        << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else {   // (the old production solver, scheme = default, is not in this build)
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/scheme = '" << sc << "': only scheme = clean (with "
      << "transport = implicit; a restart file without the key ran the old solver)"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }
  }
  // GEOMETRY GUARD (docs/dev/rad_vet_curvilinear_design.md sect. 0).  Every M1 kernel,
  // explicit and implicit, is CARTESIAN with a UNIFORM cell: divergences divide by
  // mb_size.dx1..3 (rad_vet_update.cpp, rad_vet_fluxes.cpp, rad_vet_implicit.cpp), the
  // flux vector is Cartesian, and no metric enters.  On a spherical-polar or cubed-sphere
  // mesh the run would proceed silently with that arithmetic (on sp dx2 is Delta theta
  // in radians: the transverse transport is wrong by r and r sin(theta)), so refuse it.
  // The radial stretches are refused too although they only act inside the sp/cs
  // geometry (on a plain Cartesian mesh the hydro ignores them and so would M1), so
  // that no stretched input can reach the uniform-dx kernels; and relativistic
  // coordinates, since the M1 here is the O(v/c) flat-space system.
  //
  // STAGE S1 (tests_m1/runs_5a_sp_s1) lifts the refusal for ONE configuration: a
  // spherical-polar WEDGE that does not touch a pole, with the radial stretches allowed,
  // transport = implicit, closure = eddington and the restrictions checked by
  // SphericalS1Check after ImplicitInit; STAGE S2 (tests_m1/runs_5b_sp_s2) adds the
  // chi(f) closures m1 / minerbo / kershaw there.  Everything else on sp (the poles,
  // explicit transport, vet_sc, tau) stays refused.
  // STAGE CS0 (branch m1-cs-implicit, rt_design_1003/CS_IMPLICIT_VET.md): the cubed
  // sphere for transport = implicit, closure = eddington, time_scheme = be, bicgstab +
  // rbgs_fwd only (checked in SphericalS1Check, which runs for both meshes).  x1 = r on
  // it as on sp, so the sp radial forms (A_f/V_i, dxface, the Marshak face) carry over;
  // the transverse faces take the plain two-point form over dxface (no skew factor).
  sph_geom = false;
  cs_geom = false;
  {
    Mesh *pm_ = ppack->pmesh;
    std::string why;
    const bool sp_ = pm_->use_spherical_polar && !pm_->use_cubed_sphere;
    const bool cs_ = pm_->use_cubed_sphere && !pm_->use_spherical_polar;
    if (pm_->use_spherical_polar && !sp_) {why += " mesh/use_spherical_polar";}
    if (cs_ && transport != M1_TRANSPORT_IMPLICIT) {
      why += " mesh/use_cubed_sphere with transport != implicit";
    }
    // S1: the FULL sphere (mesh/use_polar_boundary, ix2_bc = ox2_bc = polar) and the
    // theta stretch (mesh/use_grid_stretch_theta, e.g. mesh/theta_plateau) are allowed
    // for transport = implicit on sp.  Every implicit kernel takes a polar face as a
    // PHYSICAL face (flag != block, periodic) and the pole face has area 0
    // (coordinates.cpp), so "pole = reflect at a zero-area face" is exact for the E/F
    // rows; the transverse geometry reads Coordinates (areas, volumes, dxface).  The
    // remaining restrictions (implicit_bc_x2 = auto | reflect at a polar face, no vet_sc
    // / tau / vet_col_lat / vet_gd) are checked in SphericalS1Check and VetLatInit.
    const bool full_ok = sp_ && (transport == M1_TRANSPORT_IMPLICIT);
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
  impl_halo_ipc = false;
  hi_state = 0;
  hi_seq = 0;
  hi_nchk = 0;
  hi_rbuf = nullptr;
  hi_sig = nullptr;
  hi_err = nullptr;
  hi_half = 0;
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
  // ---- milestone 3b phase B
  impl_solver = M1_ISOLV_LINE_JACOBI;
  trans_on = false;
  trans_x3 = false;
  pbval_th = nullptr;
  pbval_tq = nullptr;
  halo_shell = false;
  // ---- milestone 3b phase C
  bicg_on = false;
  pbval_kr = nullptr;
  impl_prec_fb = -1;
  bcg_capped = false;
  bcg_r0rel = 0.0;
  ew_fprev = 0.0;
  ew_etaprev = 0.0;
  impl_pred = false;
  pred_ok = false;
  pred_dt = 0.0;
  // ---- milestone 3b phase D
  impl_offdiag = M1_OD_LAGGED;
  od_now = M1_OD_LAGGED;
  impl_vimp = false;
  vimp_now = false;
  t2_enth_be = false;
  pbval_vm = nullptr;
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
  // ---- milestone 3e: Anderson acceleration of the Picard map
  aa_nh = 0;
  aa_head = 0;
  aa_hasp = false;
  aa_fnp = -1.0;
  // ---- milestone 3g: the gas-radiation energy coupling
  impl_halo_direct = false;
  halo_direct_on = false;
  impl_stencil = false;
  st_edges = true;
  impl_kfuse = 0;
  impl_prec = 0;
  // ---- end of the 3a hook
  cfl_rad = pin->GetOrAddReal("rad_vet","cfl_rad",0.4);
  e_floor = pin->GetOrAddReal("rad_vet","e_floor",(FLT_MIN));
  subcycle = pin->GetOrAddBoolean("rad_vet","subcycle",true);

  // (the explicit transport's face-flux options, removed with it: only their old
  // defaults are accepted, CleanPreflight)
  if (pin->DoesParameterExist("rad_vet","scaled_prefactor") &&
      pin->GetReal("rad_vet","scaled_prefactor") != 20.0) {
    std::cout << "### FATAL ERROR: <rad_vet>/scaled_prefactor = "
              << pin->GetString("rad_vet","scaled_prefactor") << " was removed (only 20)"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // (1b) opacities, per unit mass, in code units.  kappa_e defaults to kappa_p and
  // kappa_f to kappa_p as well, so a single kappa_p gives a grey absorbing medium.
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
  // the unit conversions of the TABLE lookup, which is in kelvin and g/cm^3 (design
  // sect. 4 / the stage-2 plan sect. 4).  T[K] = T_code*temp_unit_kelvin: the same
  // direction Conduction uses with Units::temperature_cgs(), which for the He box's
  // <units> 1/1/1 with mu = 1 is mu m_u/k_B = 1.202724e-8 K per code unit.
  otab.tunit = pin->GetOrAddReal("rad_vet","temp_unit_kelvin",1.0);
  otab.dunit = pin->GetOrAddReal("rad_vet","rho_unit_cgs",1.0);
  otab.kunit = pin->GetOrAddReal("rad_vet","kappa_unit",1.0);
  if (!(otab.tunit > 0.0) || !(otab.dunit > 0.0) || !(otab.kunit > 0.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/temp_unit_kelvin, rho_unit_cgs and kappa_unit must all "
      << "be positive" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // <rad_vet>/opac_abs_rho_max (read only when named, default 0 = off, bitwise inert): an
  // ABSORPTION MASK of the table opacity for tenuous gas.  Cells with a CODE density
  // below it get kappa_P = kappa_E = 0 (no emission, no absorption) and the transport
  // opacity <rad_vet>/opac_abs_kappa_s [cm^2/g] (electron scattering, default 0.34 = 0.2
  // (1 + X) at X = 0.7) instead of the table's Rosseland mean.  Meant for an artificial,
  // optically thin, thermally decoupled ambient (ry_per_accretor's hot 300 km/s ambient,
  // design docs/dev/accretor_rhd_design.md sect. 2.3): with an edge-held table kappa_P it
  // would cool to T_rad in ~1e-6 s and collapse.  Applied inside M1TableOpacities, i.e.
  // at EVERY table lookup of the module (opacity task, the Picard re-evaluation, the T
  // solve, time2), so the mask is the same wherever kappa is used.
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
  // vet_scatter (read only when named, default off = bitwise): ONE
  // switch for electron scattering across the module.
  //  (i) the opacity split.  opacity = table: os = kappa_e, of = kappa_a = kappa_R -
  //      kappa_e (Ma+2026 sect. 2.2.4, guarded so that kappa_a >= vet_kappa_floor), with
  //      kappa_e = vet_kappa_es (vet_scatter_kappa_e = const, default 0.34 cm^2/g) or
  //      sigma_T x_e/(mu m_u) from the general EOS (= eos).  Other opacity laws: the
  //      input's own kappa_f (absorption) and kappa_s (scattering).  rho kappa_e goes to
  //      opac(M1_OP_S); every consumer of the extinction keeps the TOTAL M1_OP_T.
  //      Scattering exchanges no energy (coherent, no Compton term): the exchange keeps
  //      kappa_P / kappa_E.
  //  (ii) the closure source of vet_sc / vet_col / vet_col_lat / vet_gd, VetScatterMix
  //      (rad_vet.hpp), form vet_scatter_form = ma (Jiang 2021 eq. 6) | absorption,
  //      J next to B* vet_scatter_j = relaxed (default) | current (absorption only;
  //      ma + current FATALs below).
  // (user 10-09): on a FRESH run with a VET closure that carries ray data
  // (closure = vet_sc, or vet_col with vet_gd = true), opacity = table and a tabulated
  // general EOS (which supplies kappa_e from its Saha x_e), vet_scatter defaults ON with
  // vet_scatter_kappa_e = eos; everything else (M1 / eddington / tau / vet_col without
  // vet_gd, other opacity laws, the ideal gas) and a restart whose file lacks the key
  // keep the old default (off).  The resolved value is recorded.
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
    // vet_scatter_j default = relaxed (J = E*, the end-of-step pair of B*).
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
    // form ma + j = current is unstable: the kappa_P/kappa_R factor of the ma form
    // (~250-520 in the BSG iron bump) amplifies the lag between the relaxed
    // (end-of-step) B* and the start-of-step J = E^n; the implicit Picard loop then
    // fails to converge (NON-CONV on the BSG, bsg_1001/VET_SCATTER.md round 2).
    if (vscat_form == 1 && !vscat_jrel) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_vet>/vet_scatter_form = ma with vet_scatter_j = current "
        << "is not allowed: the kappa_P/kappa_R factor of the ma form (~250-520 in an "
        << "iron bump) amplifies the lag between the relaxed end-of-step B* and the "
        << "start-of-step J = E^n, and the implicit Picard loop does not converge "
        << "(NON-CONV on the BSG). Use vet_scatter_j = relaxed (default) or "
        << "vet_scatter_form = absorption." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::string km = "const";
    if (opacity_type == M1_OPAC_TABLE) {
      // eos where a tabulated general EOS supplies x_e (fresh runs)
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

  // (1c) matter coupling.  It reads rho, v and the gas energy from hydro's CONSERVED
  // u0 and writes u0(IEN) and u0(IM1..3) back, so it needs a <hydro> block unless every
  // opacity vanishes (the pure-transport tests of milestone 1a).
  // m1-mhd: or an <mhd> block (the same conserved layout; IEN also carries |B|^2/2)
  const bool have_mhd = pin->DoesBlockExist("mhd");
  bool have_hydro = pin->DoesBlockExist("hydro") || have_mhd;
  if (have_mhd && pin->DoesBlockExist("hydro")) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet> with both <hydro> and <mhd> (ion-neutral) is not "
      << "supported: the module couples to ONE fluid" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // <mhd> on the cubed sphere too (EmagBuild takes the C2P's orthonormal
  // |B|^2 there; the momentum layout is the hydro one: covariant on the panel basis)
  const std::string fblk = have_mhd ? "mhd" : "hydro";
  fl_on = have_hydro;
  fl_mhd = have_mhd;
  // (1c) the advective enthalpy-flux split of the E equation (design sect. 3 "Moving
  // fluid"), default off (it was on only with the explicit thick_flux = ap_hll)
  advect_split = pin->GetOrAddBoolean("rad_vet","advect_split",false);
  coupling = pin->GetOrAddBoolean("rad_vet","coupling",!opac_zero);
  gas_feedback = pin->GetOrAddBoolean("rad_vet","gas_feedback",true);
  // debug switches (see rad_vet.hpp); every default is the production value
  dbg_gas_force = pin->GetOrAddBoolean("rad_vet","dbg_gas_force",true);
  if (!dbg_gas_force && global_variable::my_rank == 0) {
    std::cout << "### WARNING: <rad_vet> DEBUG switch dbg_gas_force = false is active"
              << " -- this is not a physical configuration" << std::endl;
  }
  // (the debug switches opac_freeze, dbg_gas_heat and dbg_gas_force_trans, removed:
  // only their production values; dbg_tensor: CleanPreflight)
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
      << "the opacity enters the thick-limit FLUX only" << std::endl;
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
  // the design's gas-only EOS requirement (sect. 1)
  if (have_hydro && pin->GetOrAddBoolean(fblk,"eos_radiation",false)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet> requires a GAS-ONLY EOS: <hydro>/eos_radiation must be "
      << "false (the radiation energy and pressure are evolved, not in the EOS)"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // radiation constant, code units: the equilibrium energy density is arad*T^4 with T
  // the EOS's own code temperature.  Required whenever there is emission/absorption.
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
  // floor-gas cap (rad_vet.hpp; read only when named)
  if (pin->DoesParameterExist("rad_vet","floorcap_rho")) {
    fcap_rho = pin->GetReal("rad_vet","floorcap_rho");
    if (pin->DoesParameterExist("rad_vet","floorcap_v")) {
      fcap_v = pin->GetReal("rad_vet","floorcap_v");
    }
    if (pin->DoesParameterExist("rad_vet","floorcap_report")) {
      fcap_report = pin->GetInteger("rad_vet","floorcap_report");
    }
    // floorcap_mesh_vol (scheme = clean; default on for a fresh run, USER 11-21,
    // recorded; a restart file without the key keeps off, echoed): the counters
    // weighted by the mesh's cell volume, |dp| and L_z of the gas momentum in the frame
    // of y on every mesh (CleanVol, CleanLzArm); the state is the same
    if (clean || pin->DoesParameterExist("rad_vet","floorcap_mesh_vol")) {
      if (clean && !nw && !pin->DoesParameterExist("rad_vet","floorcap_mesh_vol")) {
        recipe::Log("<rad_vet>/floorcap_mesh_vol = false (a restart file without the "
                    "key: the old default)");
      }
      fcap_mvol = pin->GetOrAddBoolean("rad_vet","floorcap_mesh_vol",nw);
      if (fcap_mvol && !clean) {
        std::cout << "### FATAL ERROR: <rad_vet>/floorcap_mesh_vol needs scheme = clean"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
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
  // Unit audit (units-physical-inputs): c_light, arad, temp_unit_kelvin, rho_unit_cgs and
  // kappa_unit are CODE-unit constants typed into the input, not derived from <units>.
  // Print each next to the value <units> implies (general EOS: T_unit = m_u/k_B per
  // code temperature, <units>/mu divided out; ideal: times mu) and warn on a mismatch.
  // Diagnostic only: the input values are used unchanged.
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

  // (2a) <rad_vet>/force_reference.  `wb_arad` means: an EXTERNAL well-balanced scheme
  // already delivers rho*arad_ref(z) to the gas momentum (box_convection's effective
  // potential), so the coupling must hand the gas only the RESIDUAL radiative force.
  // The reference array itself comes from the problem generator (SetForceReference);
  // without it the option is a fatal, not a silent no-op.  `clean` (scheme = clean):
  // the reference is the clean solve's own radiative force of its first coupled solve,
  // i.e. of the (steady-started) initial state; the hydro carries that state as its
  // static well-balanced background, so a discrete equilibrium of both is exact.
  {std::string fr = pin->GetOrAddString("rad_vet","force_reference","none");
  if (fr.compare("none") == 0) {
    force_ref = M1_FREF_NONE;
  } else if (fr.compare("wb_arad") == 0) {
    force_ref = M1_FREF_WB_ARAD;
  } else if (fr.compare("clean") == 0 && clean) {
    force_ref = M1_FREF_CLEAN;
  } else if (fr.compare("wb_arad_md") == 0 && clean) {
    force_ref = M1_FREF_WB_MD;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/force_reference = '" << fr << "' is not a valid choice "
      << "(none | wb_arad | clean or wb_arad_md with scheme = clean)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  }
  // force_reference_work (m1-keydefault-0927): with force_reference = wb_arad the
  // default is `auto` = split when the problem generator gives the gas the reference
  // work (it sets fref_wsplit_ok before SetForceReference: box_convection, he_star_vet),
  // else full; resolved and echoed in SetForceReference.  `full` keeps the pre-0927
  // scheme.  Without wb_arad the key is read only when named and has no effect.
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
  // time2_fref_wacc (default true, USER 10-22; false = the old payment):
  // the exact payment above also for the backward-Euler steps of a hesdirk2 run, which
  // otherwise pay v' dt rho a_ref (B deck, first step: +76 energy units created, rE jump)
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

  // (2a) the RSLA start-up check of design sect. 2
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
  tau_time = 0.0;
  tau_ncall = 0.0;
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
    // the implicit multi-D solve reads (chi, n) from the short-characteristics formal
    // solution (rad_vet_vet.cpp); every other use of chi (explicit wave speeds, the
    // 1-D branch) keeps Levermore's
    eddington = false;
    vet_sc = true;
    VetRotInit(pin);   // the vet_gd rotation keys and defaults
    RotAvgParse(pin);  // vet_rotate_avg (shared key)
    // scheme = clean with clean_rot_ring named: the rotating set on the clean ring, as
    // vet_gd's (vet_gd_rotate_xfade, default false; no slope correction)
    if (clean && vsc_rot && cln_ring) {
      vgd_xf = pin->GetOrAddBoolean("rad_vet","vet_gd_rotate_xfade",false) && VgdRotOn();
      vgd_rot_slope = 0;
    }
    // vet_source_noesrc for vet_sc too (rad_vet_vet.cpp).
    // (#11): DEFAULT true on a fresh run, as on vet_col / vet_gd (the VET source of a
    // backward-Euler solve reads the physical E^n, not E^n + dt (chat/c) esrc); a restart
    // whose file lacks the key keeps false (echoed).  Bitwise wherever esrc = 0.
    vsrc_noes = pin->GetOrAddBoolean("rad_vet","vet_source_noesrc",fresh);
  } else if (cl.compare("tau") == 0) {
    // the multi-D implicit solve reads (chi, n) from the column optical depth
    // (rad_vet_tau.cpp), in the uniaxial form; every other use of chi (explicit wave
    // speeds, the 1-D branch, the explicit coupling) keeps Levermore's, as closure = m1
    eddington = false;
    tau_closure = true;
  } else if (cl.compare("vet_col") == 0) {
    // STAGE S5 (rad_vet_vetcol.cpp): the tensor of a per-column 1-D formal solution,
    // carried by the tau closure's machinery (tau_ten, built once per step)
    eddington = false;
    tau_closure = true;
    vet_col = true;
    vcol_nc = pin->GetOrAddInteger("rad_vet","vet_col_ncore",8);
    vcol_np = pin->GetOrAddInteger("rad_vet","vet_col_nsub",1);
    vcol_nmu = pin->GetOrAddInteger("rad_vet","vet_col_nmu",4);
    // the canonical key vet_every, vet_col_every its alias
    vcol_every = pin->GetOrAddInteger("rad_vet",
                                      M1KeyName(pin, "vet_every", {"vet_col_every"}), 1);
    // vet_col_fk_min (m1-sp-order2b): the lower clamp of f_K, default 0 on every mesh (a
    // restart whose file lacks the key keeps 1/3, echoed)
    vcol_fkmin = pin->GetOrAddReal("rad_vet","vet_col_fk_min",
                                   global_variable::restart_run ? (1.0/3.0) : 0.0);
    if (!(vcol_fkmin >= 0.0 && vcol_fkmin <= 1.0/3.0 + 1.0e-12)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/vet_col_fk_min must be in [0, 1/3]"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // vet_col_reflect_top (m1-sp-order2b, sp only): with implicit_bc_x1max = reflect the
    // formal solution mirrors the outgoing intensity at the top face instead of taking a
    // vacuum top (rad_vet_vetcol.cpp).  Sp default true (a restart whose file lacks the
    // key keeps false, echoed).
    if (sph_geom) {
      vcol_rtop = pin->GetOrAddBoolean("rad_vet","vet_col_reflect_top",
                                       !global_variable::restart_run);
    }
    // vet_col_surface_face (m1-sp-order2b, sp only; rad_vet.hpp): sp default true (a
    // restart whose file lacks the key keeps false, echoed); acts only together with
    // vet_col_surface_q and implicit_marshak_face = linear
    if (sph_geom) {
      vcol_sqf = pin->GetOrAddBoolean("rad_vet","vet_col_surface_face",
                                      !global_variable::restart_run);
    }
    // vet_col_surface_q: DEFAULT true since m1-defaults2 (tests_m1/runs_5e_vetcol2:
    // T-S4 L1 2.1e-3 -> 5.2e-5 at n = 256, Milne 1.8e-3 -> 1.1e-4).  A restart whose
    // file lacks the key keeps false; the resolved value is echoed.
    vcol_sq = pin->GetOrAddBoolean("rad_vet","vet_col_surface_q",
                                   !global_variable::restart_run);
    vcol_qmin = pin->GetOrAddReal("rad_vet","vet_col_surface_qmin",1.0e-3);
    vcol_qmax = pin->GetOrAddReal("rad_vet","vet_col_surface_qmax",1.0);
    {
      // the canonical key vet_source (vet_col_source its alias), read by vet_sc too
      // (VetInit) (USER 10-20): DEFAULT relaxed_esc under the new default set
      // (utils/recipe.hpp: a fresh run, or <job>/recipe = prod1019); a restart whose
      // input names neither key, or recipe = none, keeps the old default relaxed
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
      // (opt-in): relaxed with the cell's photon escape (rad_vet_vetsrc.hpp)
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
    // vet_bc_diff (the vet_sc key, rad_vet_vet.cpp): the inner diffusion
    // intensity of the formal solutions at the FACE (face) or from the first CELL (cell):
    // vet_gd / vet_col_lat start S + mu(-dS/dtau) with S of the face, 1.5 S_a - 0.5 S_b
    // (the S their half-cell carry starts from), vet_col's core rays E + (S_f - S_0)
    // (vet_col_order2 keeps its own face start).  DEFAULT face; a restart
    // whose file lacks the key keeps cell; cell reproduces older runs bitwise.
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
    // vet_seg: mean (bitwise older runs) integrates the SC segment
    // upwind point -> cell centre with chi = (chi_u + chi_0)/2 and S linear in tau; split
    // integrates its two parts at the cell face with chi_u, chi_0 and a face S limited
    // toward S_u at an opaque -> thin interface (rad_vet_vetseg.hpp); emit applies that
    // limiter only where S_u > S_0; jump (coreL-1012) adds the opaque -> thin jump
    // term.  DEFAULT emit, jump (fresh runs); a
    // restart whose file lacks the key keeps mean (bitwise older runs); files since
    // echo their value.
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
    // vet_col_lat (m1-vetcol-lat, rad_vet_vetlat.cpp): every key read only when named, so
    // that an input without them writes the same restart and parameter dump as before
    if (pin->DoesParameterExist("rad_vet","vet_col_lat")) {
      vlat_on = pin->GetBoolean("rad_vet","vet_col_lat");
    }
    // Athena++ x2rat / x3rat (not AthenaK keys) noted for VetLatInit
    for (const char *rk : {"x2rat", "x3rat"}) {
      if (pin->DoesParameterExist("mesh", rk) && pin->GetReal("mesh", rk) != 1.0) {
        vlat_rat_keys += std::string(" mesh/") + rk + " = " + pin->GetString("mesh", rk);
      }
    }
    // vet_gd (m1-vet-gd, rad_vet_vetgd.cpp): rides on the vet_col_lat interface
    if (pin->DoesParameterExist("rad_vet","vet_gd")) {
      vgd_on = pin->GetBoolean("rad_vet","vet_gd");
      if (vgd_on) {vlat_on = true;}
    }
    // (USER 11-17): the DEFAULT set of a FRESH start that names no set key,
    // recorded for restarts; a restart whose input names none keeps HEALPix nside 3.
    // (USER: <= 120 rays, rotation as in production; RAYSET2): the geodesic
    // set vet_gd_geo = 3 (92 directions), averaged over 4 rotated sets (geodef7 below):
    // the accuracy of the fixed geo 4 at 92 rays.  A 2-D axisymmetric mesh with a phi
    // extent below 2 pi, which must rotate the set onto itself (the icosahedral set has
    // no 4-fold axis), takes the staggered product set GL 9 x nphi, nphi nearest 12 with
    // 2 pi/nphi dividing the extent (FATAL if none does); nmu odd: the staggered rings
    // pair with their equatorial mirrors (GL 10 x 12: a spurious F_theta 1.8e-8 of F_r
    // 5.5e-4 on the rhse_sp steady state)
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
    // (USER: "rotate with averaging"): a FRESH vet_gd run whose set rotates
    // (vet_gd_rotate_dt unnamed, i.e. auto, or > 0; or vet_gd_rotate_every > 0) and that
    // names no vet_rotate_avg averages 4 rotated sets (on the time clock of scheme =
    // clean with vet_gd_rotate_xfade = true: the live ring); never one rotating set alone
    if (vgd_on && recipe::NewDefaults(pin) && !global_variable::restart_run &&
        !pin->DoesParameterExist("rad_vet","vet_rotate_avg")) {
      const bool rdt = !pin->DoesParameterExist("rad_vet","vet_gd_rotate_dt") ||
          pin->GetString("rad_vet","vet_gd_rotate_dt") == "auto" ||
          pin->GetReal("rad_vet","vet_gd_rotate_dt") > 0.0;
      const bool rev = pin->DoesParameterExist("rad_vet","vet_gd_rotate_every") &&
          pin->GetInteger("rad_vet","vet_gd_rotate_every") > 0;
      if (rdt || rev) {
        pin->GetOrAddInteger("rad_vet","vet_rotate_avg",4);
        if (clean && rdt && !rev) {
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
    // <rad_vet>/vet_gd_geo = n: the geodesic set (vgd_ls = -n; needs
    // vet_gd_gl_nmu = 0 and no vet_gd_ls)
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_geo")) {
      const int ng_ = pin->GetInteger("rad_vet","vet_gd_geo");
      if (ng_ > 0) {vgd_ls = -ng_;}
    }
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_rotate_every")) {
      vgd_rot = pin->GetInteger("rad_vet","vet_gd_rotate_every");
    }
    // the rotation on a physical-time clock (absent or 0 = the cycle clock of
    // vet_gd_rotate_every, bitwise).  vet_gd_rotate_dt = <dt_rot> | auto. auto = 20 x the
    // FIRST step's dt (the dt of the first Mesh::NewTimeStep, after its global min, so
    // the same on every rank), fixed for the run and recorded in the parameters as the
    // number (VetGdAlignDt), so a restart keeps it.  DEFAULT auto on a FRESH vet_gd run
    // that names neither vet_gd_rotate_dt nor vet_gd_rotate_every (either named =
    // explicit input, kept) and has no vet_gd_beam (which excludes a
    // rotation); a restart whose file lacks both keys keeps the old default (no
    // rotation). The sub-keys are read (and echoed) with the physical clock:
    // vet_gd_rotate_t0 (0), _align (true), _stretch (0.02 on a fresh run, 0 on a restart
    // whose file lacks it), _slope (correct on a fresh run, keep on a restart whose file
    // lacks it).
    VetRotInit(pin);   // shared with vet_sc (rad_vet_vetgd.cpp)
    // clean staggered sp and cs, fresh, with a rotating set on the time clock: the set
    // time-averaged by default, vet_rotate_avg = 4 with the cross-fade (clean_rot_ring
    // live), recorded so that a restart keeps them; a restart file without the keys keeps
    // the plain rotation
    if (clean && (sph_geom || cs_geom) && !global_variable::restart_run &&
        VgdRotOn() &&
        (vgd_rot_dt > 0.0 || vgd_rot_auto)) {
      pin->GetOrAddInteger("rad_vet","vet_rotate_avg",4);
      pin->GetOrAddBoolean("rad_vet","vet_gd_rotate_xfade",true);
    }
    // (a restart of the old production's averaged set, whose file has no cross-fade:
    // the sp sc rays take it, the live ring; without it the first steps from the
    // production state take the K-set closure's jumps, which its R136 restart cannot)
    if (clean && sph_geom && global_variable::restart_run &&
        VgdRotOn() && pin->DoesParameterExist("rad_vet","vet_rotate_avg") &&
        pin->GetInteger("rad_vet","vet_rotate_avg") > 1 &&
        !pin->DoesParameterExist("rad_vet","vet_gd_rotate_xfade")) {
      pin->GetOrAddBoolean("rad_vet","vet_gd_rotate_xfade",true);
    }
    RotAvgParse(pin);   // (also after vet_sc's VetRotInit)
    // scheme = clean: no slope correction of a rotating set (the closure of CleanFvSets),
    // set before the first VetGdRotStep (a restart on a block boundary acts at once).
    // clean_rot_ring = jump: a block's first step is backward Euler from the first step
    // on (a restart on a block boundary too)
    if (clean && vgd_on && VgdRotOn()) vgd_rot_slope = 0;
    if (cln_ring == 2) vgd_rot_slope = 2;
    // vet_gd_rotate_xfade (read only when named; rad_vet.hpp)
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_rotate_xfade")) {
      // (a fixed set has nothing to cross-fade: the key is inert there)
      vgd_xf = pin->GetBoolean("rad_vet","vet_gd_rotate_xfade") && VgdRotOn();
      // (vet_gd_rotate_xfade_m > 0 and its _lean, _pen: the old scheme's sparse
      // next-set sweeps, removed)
      if (vgd_xf && pin->DoesParameterExist("rad_vet","vet_gd_rotate_xfade_m") &&
          pin->GetInteger("rad_vet","vet_gd_rotate_xfade_m") != 0) {
        std::cout << "### FATAL ERROR <rad_vet>/vet_gd_rotate_xfade_m = "
                  << pin->GetInteger("rad_vet","vet_gd_rotate_xfade_m")
                  << " was removed (only 0)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      // the fused two-set sweep (VetGdXffOk decides where it applies)
      if (vgd_xf) {vgd_xff = true;}
      if (vgd_xf && global_variable::my_rank == 0) {
        std::cout << "<rad_vet> vet_gd_rotate_xfade = true: the closure is "
                  << "cross-faded from the current to the next ray set over each "
                  << "rotation block (two sweeps per build)" << std::endl;
      }
    }
    // the cubed sphere with clean_rays = sc: vet_gd's own sweep (intensities to cln_fvi)
    cln_cssc = clean && pmy_pack->pmesh->use_cubed_sphere && vgd_on;
    // clean_thin_donor (default off): the thin share w c h E of every face from its
    // donor alone, the half range h out through the face the donor's own (its rays
    // about the face normal: the six faces partition the donor's free-streaming
    // outflow by projected area, w <= 1 keeps it below), no lagged rays' flux DG; on
    // the cubed sphere with the cells' own half ranges (cln_rfh: CLN_HO, which the
    // rotating ring's offsets carry as the closure)
    cln_tdn = pin->GetOrAddBoolean("rad_vet","clean_thin_donor",false);
    cln_rfh = cln_cssc && cln_tdn;
    // (a face's thin share w from the thicker of its two cells' C1)
    if (pin->GetOrAddBoolean("rad_vet","clean_w_face_min",false)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_w_face_min = true was removed "
                << "(only false)"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // rotwb2-1121: vet_rotate_avg_ema_wb keeps the ring's closure of a static state, but
    // the cs face half ranges (cln_fh, CleanFaceHalf) are the two swept sets' (no offset)
    if (cln_cssc && vra_ema_wb && !cln_rfh && global_variable::my_rank == 0) {
      std::cout << "### NOTE <rad_vet>/vet_rotate_avg_ema_wb on the cubed sphere with "
                << "clean_rays = sc: the face half ranges follow the two swept sets, so "
                << "a static state is not kept to round-off" << std::endl;
    }
    // clean_rot_ring not named: live wherever the vet_gd rays rotate in time
    // with a cross-fade (decided here, before a restart file's channels are counted)
    if (clean && !cln_ring_named && vgd_on && vra_on && g_vra_adv &&
        vgd_xf && (vgd_rot_dt > 0.0 || vgd_rot_auto) && vgd_rot_align &&
        vra_k <= M1_RA_NF) {
      cln_ring = 1;
    }
    // the inward rays of the top shell start from the incident bath
    // implicit_ebath_x1max of a MARSHAK outer x1 (vet_sc's vet_bc_bath, the same key and
    // default: true on a fresh run; a restart whose file lacks the key keeps the vacuum
    // top).  With implicit_ebath_x1max = 0 (production) the sweep is unchanged, bitwise.
    if (vgd_on) {
      vgd_tbath = pin->GetOrAddBoolean("rad_vet", "vet_bc_bath",
                                       !global_variable::restart_run);
    }
    // vet_gd_inner = diff | emit, vet_gd_beam* (read only when named).
    // vet_bc_inner / vet_bc_reflect, shared with vet_sc
    // (VetBcInit); vet_gd_inner is read as the alias of vet_bc_inner
    // vet_interp = auto | log | corner, shared with vet_sc (lin
    // is vet_sc's old form only).  vet_gd: log = its old bilinear in ln chi, ln S;
    // corner = one segment per stencil corner (VetGdSweep); auto = corner where the
    // stencil's chi contrast exceeds kVetCornerRatio, log elsewhere (as vet_sc).  DEFAULT
    // auto on a fresh run; a restart whose file lacks the key keeps log.
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
    // vet_path = cells | none, shared with vet_sc (VetInit): the step of a
    // ray that jumps more than one lateral cell subdivided at the theta-row / phi-layer
    // nodes where its cells have a chi (or non-thin S) contrast > kVetCornerRatio
    // (VetGdSweep).  OPT-IN: DEFAULT none (fresh and restart; bitwise 281dd686), cells
    // only when named.
    if (vgd_on) {
      const std::string vp = pin->GetOrAddString("rad_vet","vet_path", "none");
      if (vp.compare("none") == 0) {
        vet_path = 0;
      } else if (vp.compare("cells") == 0) {
        vet_path = 1;
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<rad_vet>/vet_path = " << vp << " is not a choice "
                  << "(cells | none)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    if (vgd_on) {VetBcInit(pin, true);}
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_allow_clamp")) {
      vgd_allow_clamp = pin->GetBoolean("rad_vet","vet_gd_allow_clamp");
    }
    // (read only when named, default false = bitwise): the band keyed on the global cell
    // + the cell-window clamp (layout-independent sweep; sph2.inc) prod1103d default
    // (USER 11-03): vet_gd_band_global = true on a FRESH spherical-polar FULL-SPHERE
    // start with vet_gd on (polar boundaries at both poles, x3 periodic over 2 pi;
    // recipe::NewDefaults) when the key is not named; recorded for restarts; a restart
    // whose input lacks the key keeps false (bitwise); an explicit key wins.
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
    // vet_gd_band_ring (read only when named; bitwise, VetGdInit): the band intensities
    // lent per shell from a few wide slabs; unnamed: on where it applies, except the
    // raw band of spherical polar (speed over memory); true there saves its memory
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_band_ring")) {
      vgd_ring_named = true;
      vgd_ring_par = pin->GetBoolean("rad_vet","vet_gd_band_ring");
    }
    // prod1104b default (USER 11-04): vet_gd_halo_compact = 2 on a FRESH vet_gd start
    // (recipe::NewDefaults) when not named (bitwise; 29-40 % faster on cs vet_gd runs;
    // the vet_gd setup still switches it off where it does not apply); recorded for
    // restarts; a restart whose input lacks the key keeps 0; an explicit key wins.
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
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_cs_xgroups")) {
      // speed only (bitwise any value); 1 = the plain per-shell exchange
      vgd_csxg = std::max(1, pin->GetInteger("rad_vet","vet_gd_cs_xgroups"));
    }
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_twin") &&
        pin->GetBoolean("rad_vet","vet_gd_twin")) {
      std::cout << "### FATAL ERROR <rad_vet>/vet_gd_twin = true: the twin is not in "
                   "this "
                << "build" << std::endl;
      std::exit(EXIT_FAILURE);
    }
#if VET_DEBUG
    vgd_ring_poison = (std::getenv("VGD_RING_POISON") != nullptr) &&
                      (std::atoi(std::getenv("VGD_RING_POISON")) != 0);
#endif
    // the same fresh default for a plain vet_col (no vet_gd), as
    // vet_sc.  vet_gd without the twin keeps false: there the pair is needed together
    // (NOTE-2026-10-09-viper-sp-seam-mass-leak: noesrc alone exposes the gd ray effect)
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
    // vet_gd_theta_wall (read only when named): periodic (default, the old
    // theta periodicity of the sweep) | mirror (specular theta walls) | vacuum (no
    // incoming intensity at a physical theta wall) | extrap (zero gradient); vetgd.cpp
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
    } else if (vgd_on && clean && pmy_pack->pmesh->use_spherical_polar &&
               !pmy_pack->pmesh->use_polar_boundary) {
      // (unnamed, scheme = clean): the rays' theta walls follow the mesh: a reflecting
      // wall mirrors, a periodic one wraps
      vgd_thmir = (pmy_pack->pmesh->mesh_bcs[BoundaryFace::inner_x2] ==
                   BoundaryFlag::reflect);
    }
    // the exact halo lists, made from the reads alone, where they apply (MPI halo with
    // vet_gd_halo_compact > 0, shell lists; VetGdInit switches them off elsewhere)
    vgd_hx_on = true;
    // vet_gd_halo_ipc (default below)
    // (default on for clean on CUDA without a theta wall: VetGdSetup keeps it on the raw
    // 3-D sweep only; off-node partners take the MPI path)
#if defined(KOKKOS_ENABLE_CUDA)
    vgd_ipc_on = vgd_on && clean && !vgd_thmir && vgd_thvac == 0;
#endif
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_halo_ipc")) {
      vgd_ipc_on = pin->GetBoolean("rad_vet","vet_gd_halo_ipc");
    }
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_replace")) {
      vgd_replace = pin->GetBoolean("rad_vet","vet_gd_replace");
    }
    // (USER 11-13 "make some defaults", the R136/HG/RSG key set): the full tensor is
    // the default on every mesh; recipe none keeps the old cubed sphere default false
    // (D_rr only).  Name false for D_rr only.
    if (vgd_on && pmy_pack->pmesh->use_cubed_sphere &&
        recipe::GetMode(pin) == recipe::Mode::none) {vgd_tan = false;}
    if (vgd_on && pin->DoesParameterExist("rad_vet","vet_gd_tangential")) {
      vgd_tan = pin->GetBoolean("rad_vet","vet_gd_tangential");
    }
    if (vlat_on) {
      if (pin->DoesParameterExist("rad_vet","vet_col_lat_taucut")) {
        vlat_taucut = pin->GetReal("rad_vet","vet_col_lat_taucut");
      }
      // scheme = clean: its rays sweep from the first shell (no depth cut), decided here
      // so that every build, also one before the first solve (CleanInit), uses it
      if (clean) vlat_taucut = std::numeric_limits<Real>::max();
      // vet_col_lat_cut_thin (VetLatCsFinish).  DEFAULT 1 on a fresh
      // run; a restart whose file lacks the key keeps 0 (the depth cut alone, bitwise)
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
      // vet_col_lat_posfix = drop (default, bitwise) | lag (read only when
      // named): what the positivity trigger does to the lateral D_r,lat / Q operator
      // terms (rad_vet.hpp vlat_posfix_lag)
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
      // (stage S4, USER 10-18): on the cubed sphere D_r,lat is the Cartesian
      // vector v (VetGdCsMoments) in the cs forms M1CsLat0 / M1CsLatF (VetLatOpCs /
      // VetLatStencilAddCs): the default is the sp one (operator)
    }
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/closure = '" << cl << "' not implemented "
      << "(m1 | minerbo | kershaw | eddington | vet_sc | tau | vet_col)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (vet_sc && transport != M1_TRANSPORT_IMPLICIT) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/closure = vet_sc needs <rad_vet>/transport = implicit "
      << "on a multi-D mesh" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (tau_closure && transport != M1_TRANSPORT_IMPLICIT) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/closure = tau needs <rad_vet>/transport = implicit "
      << "on a multi-D mesh" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (chi_kind != M1_CHI_LEVERMORE && transport == M1_TRANSPORT_EXPLICIT) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/closure = '" << cl << "' needs an implicit transport: "
      << "the explicit HLL wave speeds are Levermore's" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  }

  // the M1 face states read two cells beyond the face (PLM of the old <rad_vet>/
  // reconstruct, the one value left, CleanPreflight)
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  if (indcs.ng < 2) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet> requires at least 2 ghost zones, but <mesh>/nghost="
      << indcs.ng << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // (2) who owns the mesh timestep.  With sub-cycling on, the radiation CFL is NOT
  // folded into Mesh::NewTimeStep -- unless nothing else sets a timestep at all, in
  // which case the mesh step would be unbounded.
  bool other_sets_dt = (pin->DoesBlockExist("hydro") || pin->DoesBlockExist("mhd") ||
                        pin->DoesBlockExist("z4c") || pin->DoesBlockExist("particles"));
  sets_mesh_dt = (!subcycle) || (!other_sets_dt);

  // (3) PD-ARS explicit weights.  With S == 0 the tableau (design sect. 5) is
  //   U1   = U^n + dt L(U^n)
  //   U^n1 = U^n + dt/2 [L(U^n) + L(U1)]
  // i.e. Heun / SSP-RK2 written in the AthenaK gam0/gam1/beta form.
  gam0[0] = 1.0;  gam1[0] = 0.0;  beta[0] = 1.0;
  gam0[1] = 0.5;  gam1[1] = 0.5;  beta[1] = 0.5;

  dtnew = (FLT_MAX);
  dt_sub = 0.0;
  nsub = 1;

  // (4) allocate arrays
  int nmb = std::max((ppack->nmb_thispack), (ppack->pmesh->nmb_maxperrank));
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  Kokkos::realloc(u0, nmb, M1_NVAR, ncells3, ncells2, ncells1);
  Kokkos::realloc(u1, nmb, M1_NVAR, ncells3, ncells2, ncells1);
  // (memory): the face fluxes are the explicit stage's (CalculateFluxes,
  // RKUpdate, flux correction); the implicit task list never adds them
  {const bool expl = (transport < M1_TRANSPORT_IMPLICIT_X1);
  const int f3 = expl ? ncells3 : 1, f2 = expl ? ncells2 : 1, f1 = expl ? ncells1 : 1;
  Kokkos::realloc(uflx.x1f, nmb, M1_NVAR, f3, f2, f1);
  Kokkos::realloc(uflx.x2f, nmb, M1_NVAR, f3, f2, f1);
  Kokkos::realloc(uflx.x3f, nmb, M1_NVAR, f3, f2, f1);}
  // rho*kappa is needed in the GHOST cells too: the face opacity at the first and last
  // active face is the arithmetic mean over a ghost and an active cell.
  // vet_scatter: a 4th component, rho kappa_e (M1_OP_S)
  Kokkos::realloc(opac, nmb, vscat ? (M1_NOPAC + 1) : M1_NOPAC, ncells3, ncells2,
                  ncells1);
  Kokkos::deep_copy(opac, 0.0);
  if (coupling) {
    Kokkos::realloc(ugas1, nmb, 4, ncells3, ncells2, ncells1);
  }
  // force_reference: allocated (and zeroed) here so that the kernel always has a valid
  // View to capture; the pgen REPLACES it with its own through SetForceReference.
  if (force_ref != M1_FREF_NONE) {
    Kokkos::realloc(arad_ref, nmb, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(arad_ref, 0.0);
  }
  // m1-mhd: |B|^2/2 of the current (b0) and U^n (b1) fields, every cell incl. ghosts
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

  // (4b) MILESTONE 3a hook: the implicit solver's own parameters, checks and arrays.
  // Returns immediately with transport = explicit.
  ImplicitInit(pin);
  // clean sc rays with a reflecting top: its mirror reads the last build (rad_vet.hpp)
  cln_scmem = clean && vgd_on && sph_geom && ibc_x1max == M1_IBC_REFLECT;
  // the ray cut (CleanRayCut) wherever it applies: the staggered layout's sc rays on 3-D
  // spherical polar or the cubed sphere, one MeshBlock along x1, no reflecting top
  if (clean) {
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
      // (the cut and its count in a restart file, CleanRayCut)
      cln_rcut = pin->GetOrAddInteger("rad_vet","clean_ray_cut_now",-1);
      cln_rcut_n = pin->GetOrAddInteger("rad_vet","clean_ray_cut_count",0);
      pin_rcut_ = pin;
    }
  }
  // <rad_vet>/polar_ring_avg (the implicit M1 on hydro's coarse polar phi groups):
  // removed
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
  // <mesh>/cs_seam_blockend = extend also on the radiation moments (E, F; the quadratic
  // along-seam resample of bvals_cc.cpp): the shared second pass of RecvAndUnpackCC,
  // clamped to the PANEL at an interior block end (default clamp, bitwise)
  pbval_u->cs_blockend_ext_ = CsSeamBlockendExtend(ppack, pin);
  pbval_u->InitializeBuffers(M1_NVAR);
  // (F2, F3) is the tangential pair (what the default (IVY, IVZ) rule gives by luck);
  // E and F1 (radial) are scalars.  Stated explicitly, see SetVectorPairs.
  pbval_u->SetVectorPairs(M1_NVAR, {{M1_F2, M1_F3}});

  if (global_variable::my_rank == 0) {
    std::cout << "<rad_vet>: c=" << c_light << " chat/c=" << chat_over_c
              << " cfl_rad=" << cfl_rad
              << " advect_split=" << (advect_split ? "true" : "false")
              // the closure as named in the input (was eddington ? eddington : m1, which
              // printed m1 for minerbo/kershaw/vet_sc/tau/vet_col; the printout only)
              << " closure=" << (vet_col ? "vet_col" : vet_sc ? "vet_sc" :
                                 tau_closure ? "tau" : eddington ? "eddington" :
                                 (chi_kind == M1_CHI_MINERBO) ? "minerbo" :
                                 (chi_kind == M1_CHI_KERSHAW) ? "kershaw" : "m1")
              << " subcycle=" << (subcycle ? "true" : "false") << std::endl;
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
//! \brief the per-column incident bath E_bath(m,k,j) of the outer x1 Marshak face (see
//! rad_vet.hpp, iebath_hi_col); also allocates ibadv_x1max (the lab-frame advective
//! energy flux through that face, filled by the implicit face kernel)

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
//! \brief force_reference_work = split, time_scheme = be: the accumulator of the
//! reference work the gas receives in the hydro stages (rad_vet.hpp, fref_wacc)

bool RadiationVet::FrefWaccOn() const {
  // time2_fref_wacc: also the backward-Euler steps of a hesdirk2 run (no
  // valid slope: the first step, a redo; not a helper solve, not multi-rate / csplit)
  const bool t2be = fref_wacc_t2 && (time_scheme == M1_TIME_HESDIRK2) &&
                    (!t2_ok || FrefWaccStage()) && !t2_inprobe && (csplit == 0);
  return fref_wsplit && ((time_scheme == M1_TIME_BE) || t2be) &&
         (mr_every <= 1);
}

bool RadiationVet::FrefWaccStage() const {
  return clean && t2_ok;
}

Real RadiationVet::FrefWaccGam0(Real bdt, Real dt_mesh) const {
  if (FrefWaccStage()) return 0.0;   // (each stage solve pays its hydro stage's work)
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
//! \brief the closure tensor <hydro>/rad_signal_speed reads (hydro_newdt.cpp uses this)

int RadiationVet::RssRstMode() const {return RssMode() + 16;}

int RadiationVet::RssMode() const {
  // (scheme = clean: its own current Eddington tensor, cln_s, which the stage solve
  // rebuilds every step; tau_ten is not rebuilt there; not before CleanInit: a restart's
  // closure, clean_closure = iterate, precedes the cubed sphere's geometry cln_gn)
  if (clean && cln_s.extent(0) > 0 && cln_init) return 5;
  if (tau_closure && tau_ten.extent(0) > 0) return 2;
  if (vet_sc && vet_cell.extent(0) > 0) return vet_full ? 4 : 3;
  if (!eddington && !tau_closure && !vet_sc) return 1;
  return 0;
}

int RadiationVet::RssRstNch(int mode) const {
  // the extended restart state (mode + 16): also the other opacity components
  // and, with the tau closure, the tensor's lateral slots (both are read before the first
  // rebuild of a step: time2_k1_stiff reads opac P/E at the stage formation)
  if (mode >= 16) {
    const int b = mode - 16;
    return RssRstNch(b) + (opac.extent_int(1) - 1)
           + ((b == 2 && tau_ten.extent_int(1) > 4) ? (tau_ten.extent_int(1) - 4) : 0);
  }
  if (mode == 2 || mode == 3) return 5;
  if (mode == 4) return 4;
  if (mode == 5) return 7;   // (opac_T and the six components of cln_s's D)
  return 1;
}

void RadiationVet::RssSlot(int mode, int n, DvceArray5D<Real> *&a, int &c) {
  if (mode >= 16) {   // the base slots, then opac (all but T), then tau_ten 4..
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
//! \brief the signal-speed inputs of MeshBlocks m0 .. m0+nmb-1, channel by channel, into
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
//! \brief before the first dt of a run with <hydro>/rad_signal_speed: put back the
//! signal-speed inputs a restart file carries (the restarted run then computes the dt
//! of the straight run), or -- fresh start, or a file without them -- fill opac from the
//! current state so that the first steps see the radiation-modified speed, not the gas
//! speed alone (He wedge: 5.7x the safe dt).  The closure tensor then starts from what
//! the module has (E/3 or M1Chi(F/cE)).

void RadiationVet::RssPrime(bool restart) {
  if (rss_stage_mode >= 16) {
    // the extended record; the tensor slots only into the same closure
    const int mode = rss_stage_mode, b = mode - 16, nb = RssRstNch(b);
    const bool tens = (b == 2 && tau_closure) || (b == 5 && RssMode() == 5) ||
                      ((b == 3 || b == 4) && vet_sc && vet_cell.extent(0) > 0 &&
                       (vet_full == (b == 4)));
    if (b == 2 && tau_closure && !tau_ready) {TauClosureInit();}
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
      if (!tau_ready) {TauClosureInit();}
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
// destructor

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetBcInit
//! \brief the x1 boundary intensity of the formal solutions, ONE
//! rule for vet_sc and vet_gd.
//!   vet_bc_inner = auto | diff | emit (alias vet_gd_inner on vet_gd), the rays entering
//!     at a MARSHAK or FLUX inner x1 face: diff = the diffusion intensity of the bottom
//!     cell (each closure's form, vet_bc_diff); emit = the face's own emission, the
//!     isotropic bath implicit_ebath_x1min of a Marshak face or the black sphere 4 F/c of
//!     a flux face (F = implicit_flux_x1min); auto = emit at a Marshak face, diff at a
//!     flux face.  DEFAULT auto on a fresh run.  A restart whose file lacks the key keeps
//!     the closure's old rule: vet_sc auto if its vet_bc_bath is on (the bath at a
//!     Marshak inner face, m1-h2div), else diff; vet_gd diff (vet_gd_inner unset).  A
//!     fresh vet_sc input that names vet_bc_bath = false (both ends off, the old meaning)
//!     and not vet_bc_inner keeps diff too.
//!   vet_bc_reflect = auto | mirror | open, a REFLECTING x1 end: mirror = the entering
//!     rays carry the leaving intensity of the mirror ray (specular reflection, H = 0 at
//!     the face, the M1 solve's zero-flux face); open = the old vet_sc rule (diffusion
//!     intensity at x1min, VACUUM at x1max: inconsistent with the zero-flux M1 face).
//!     auto = the diffusion intensity at a reflecting x1min (consistent with H = 0; the
//!     sp/Cartesian twin pair rswt_r1000 / rswt_cart agrees better with it) and the
//!     mirror at a reflecting x1max on vet_sc; vet_gd cannot mirror (FORCED: its
//!     direction set is global, so the specular image of a direction is not in the set;
//!     the vet_col column of the hybrid mirrors at a reflecting x1min, vet_col_order2,
//!     and vet_gd refuses a reflecting x1max).  DEFAULT auto on a fresh run; a restart
//!     whose file lacks the key keeps open.  vet_sc mirror: single MeshBlock along x1
//!     (VetShortChar).

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
    vi = bath ? "auto" : "diff";     // vet_sc: the old vet_bc_bath rule
  } else if (gd && rst) {
    vi = "diff";                     // vet_gd: the old default (vet_gd_inner unset)
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
    // the inner end keeps the diffusion intensity (H = 0 of a zero-flux face is its
    // isotropic part; the rswt_cart / rswt_r1000 pair agrees 1.5x better with it than
    // with the mirror, and vet_gd cannot mirror); the outer end mirrors (its old vacuum
    // contradicts the zero-flux M1 face)
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
//! \brief vet_bc_inner resolves to emit at the current inner x1
//! face (emit: a Marshak or flux face; auto: a Marshak face)

bool RadiationVet::VetBcEmit() const {
  if (ibc_x1min == M1_IBC_OUTFLOW) {return true;}   // (scheme = clean) nothing enters
  if (vet_bc_inner == 1) {
    return (ibc_x1min == M1_IBC_MARSHAK) || (ibc_x1min == M1_IBC_FLUX);
  }
  if (vet_bc_inner == 2) {return (ibc_x1min == M1_IBC_MARSHAK);}
  return false;
}

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
  if (vgd_ipc_state == 1) {     // the HIP port
    (void)hipDeviceSynchronize();
    for (void *q : vgd_ipc_open) {if (q != nullptr) {(void)hipIpcCloseMemHandle(q);}}
    (void)hipFree(vgd_ipc_rbuf);
    (void)hipFree(vgd_ipc_sig);
    (void)hipFree(vgd_ipc_err);
  }
#endif
  if (tau_closure) {TauClosureReport();}
  ReportCounters();
  ImplicitReport();   // milestone 3a; a no-op in transport = explicit
  Time2Report();      // time_scheme = hesdirk2; a no-op otherwise
  if (vet_sc) {VetReport();}
  VetFree();
  delete pbval_u;
  if (pbval_th != nullptr) {delete pbval_th;}
  if (pbval_tq != nullptr) {delete pbval_tq;}
  if (pbval_kr != nullptr) {delete pbval_kr;}
  if (pbval_vm != nullptr) {delete pbval_vm;}
  if (pbval_d3 != nullptr) {delete pbval_d3;}
  if (pbval_sg != nullptr) {delete pbval_sg;}
  if (pbval_vl != nullptr) {delete pbval_vl;}
  if (pbval_ea != nullptr) {delete pbval_ea;}
  if (pbval_vs != nullptr) {delete pbval_vs;}
  if (pbval_hr != nullptr) {delete pbval_hr;}
  if (pbval_gd != nullptr) {delete pbval_gd;}
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ReportCounters
//! \brief one line at the end of the run with the cost and the health of the implicit
//! energy-exchange solve.  A non-zero failure count means cells left the bracketed
//! Newton at M1_MAXIT without reaching M1_RTOL and is a result to distrust.

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
//! \brief Set the substep dt and the number of substeps for a mesh step dt_mesh.
//! N_sub = ceil(dt_mesh/dt_rad) with dt_rad = cfl_rad*min(dx)/chat (= dtnew), and every
//! substep uses dt_mesh/N_sub.  When this module owns the mesh timestep there is nothing
//! to sub-cycle and N_sub = 1.

int RadiationVet::SetSubsteps(Real dt_mesh) {
  if (sets_mesh_dt || !subcycle) {
    nsub = 1;
  } else {
    nsub = 1;
    if (dtnew > 0.0 && dt_mesh > dtnew) {
      nsub = static_cast<int>(std::ceil(dt_mesh/dtnew));
      if (nsub < 1) {nsub = 1;}
    }
  }
  dt_sub = dt_mesh/static_cast<Real>(nsub);
  return nsub;
}


//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RotAvgParse
//! \brief <rad_vet>/vet_rotate_avg (+ _weight), after the rotation keys of
//! the closure (vet_gd: in the vet_col block; vet_sc: in its closure branch)

void RadiationVet::RotAvgParse(ParameterInput *pin) {
  // vet_rotate_avg = K (read only when named; absent or < 2 = off, bitwise;
  // a restart whose file lacks the key keeps the old behaviour): the closure averaged
  // over the last K rotated direction sets (rad_vet_rotavg.cpp).  The key is shared by
  // every rotating closure; vet_gd: needs vet_gd_rotate_dt or _every (inert otherwise)
  if (pin->DoesParameterExist("rad_vet","vet_rotate_avg")) {
    vra_k = pin->GetInteger("rad_vet","vet_rotate_avg");
    if (vra_k > M1_RA_NF) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/vet_rotate_avg must be <= " << M1_RA_NF
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    vra_on = (vra_k >= 2) && (vgd_on || vet_sc);
    // vet_rotate_avg_weight = j (default: the sets weighted by their formal-solution
    // J, i.e. the normalised moments of the union of the K sets) | equal
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
  // vet_rotate_avg_ema = N (default 128 + change_form step for a fresh run
  // whose direction set rotates, recorded; else read only when named; 0 = off; see
  // rad_vet.hpp): the exponential moving average over the rotated sets (one offset slot)
  const bool nw = !global_variable::restart_run;
  const bool erot = nw && (vgd_on || vet_sc) && (vgd_on ? VgdRotOn() : vsc_rot);
  if (erot || pin->DoesParameterExist("rad_vet","vet_rotate_avg_ema")) {
    vra_ema_n = pin->GetOrAddReal("rad_vet","vet_rotate_avg_ema",128.0);
  }
  // vet_rotate_avg_ema_wb acts on the EMA ring only (no silent no-op)
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
      // (linear, the fresh default again since USER 11-21: step costs the ring's time
      // order, E/F p 1.66/1.39 vs 1.95/1.85 on deck D; a restart keeps its recorded form)
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
