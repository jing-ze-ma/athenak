//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file deep_hot_jupiter.cpp
//! \brief Problem generator for the deep hot Jupiter.
//!
//! REFERENCE: Heng, Menou, Phillipps, MNRAS, 413, 2380 (2011); Deitrick, Mendonça, Schroffenegger, Grimm, Tsai, Heng, ApJS, 248, 30 (2020)
//!
//! problem/rt_use_cons (default TRUE here since 2026-09-22, as in box_convection;
//! red_giant defaults it to its mode-3 flag).  Where the two-stream sweep reads the
//! cell's thermodynamic state from.  With it FALSE the sweep reads w0, the primitives the
//! PREVIOUS ConToPrim wrote, i.e. the state at the START of the stage: the RK update,
//! the explicit source terms and the implicit radial conduction have all moved the gas
//! since, so the radiative source is computed against a state that is one whole stage
//! stale.  That lag is not harmless -- it is a delayed thermostat, and a delayed
//! thermostat is a numerical FORCING whose amplitude is linear in dt.  It was measured
//! driving the organ-pipe acoustic modes of the closed FeCZ boxes (saturated v_rms
//! strictly monotone in dt and in nothing else; see the memory note
//! rt-source-dt-forcing.md), which is why box_convection hands the sweep the post-RK
//! state and why TRUE is the recommended setting for new runs here as well.
//!
//! IT WAS FALSE UNTIL 2026-09-22, so that every run ever made from these inputs
//! reproduced bit for bit.  The default is now TRUE: the two states differ by O(dt) in
//! the radiative source and the lagged one is the wrong physics, so an input that wants
//! the old answer has to say rt_use_cons = false.  sp_mhd_prod3-type inputs, which never
//! set the line, therefore change; cs_mhd_prod3 already set it true and does not.
//! (The conserved state is not guaranteed positive after an RK stage; the sweep's
//! EintFromCons path carries its own non-positive guard -- see utils/two_stream_rt.hpp.)


#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

// tusplit-1102: the private header shared with deep_hot_jupiter_rt_{hse,src,diag}.cpp
#include "deep_hot_jupiter_rt.hpp"
#include "utils/ck_2nd.hpp"   // ck2nd-1102 (keys read only when named)

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  bool use_etotgrav = false;
  bool use_wellbalance_static = false;
  bool use_wellbalance_dynamic = false;
  const bool use_spherical_polar = pmy_mesh_->use_spherical_polar;
  // The cubed sphere needs the cell's PANEL to turn (x2,x3) into a direction; see
  // CSCellAngles.  Captured by value / as a DualArray so the device lambdas below can
  // read it without touching `this`.
  const bool use_cubed_sphere_ = pmy_mesh_->use_cubed_sphere;
  bool user_srcs = pin->GetOrAddBoolean("problem","user_srcs",false);
  if (user_srcs) user_srcs_func = SourceFunc;
  // problem/flux_hst (default true): append the nine flux columns of DhjFluxHistory to
  // every history output.  Read without recording it in the input, so that a restart
  // file written by a run that does not set the key is byte-for-byte what it was; an
  // explicit flux_hst = false leaves the history exactly as before.
  // Needs the curvilinear area/volume Views (cubed sphere or spherical polar); a
  // Cartesian box has 1x1x1x1 placeholders there, so it is skipped.
  const bool flux_hst = (pin->DoesParameterExist("problem","flux_hst")
                         ? pin->GetBoolean("problem","flux_hst") : true)
                        && (pmy_mesh_->use_cubed_sphere || use_spherical_polar);
  int fh_rst_req = -1;   // problem/flux_hst_rst: -1 unset, 0 false, 1 true
  if (flux_hst) {
    user_hist = true;
    user_hist_func = DhjFluxHistory;
    // problem/flux_hst_rst (rthst-1101): write the radiative face-flux cache these
    // columns read into the implicit-ck restart block, so that the rows a continuation
    // writes before its first full RT call equal the uninterrupted run's (they had no
    // radiative part).  Diagnostic only.  Default = ck_implicit (resolved below, once
    // ck_implicit is known; no block exists without it).  Read without recording it in
    // the input; false writes restart files byte-identical to the ones before the key.
    if (pin->DoesParameterExist("problem","flux_hst_rst")) {
      fh_rst_req = pin->GetBoolean("problem","flux_hst_rst") ? 1 : 0;
    }
    // stswall-1029: the RKG resistive energy flux through the top and the wall (see
    // Resistivity::track_x1_bflx); folded into Etot_top / Etot_bot by DhjFluxHistory
    {
      MeshBlockPack *prk = pmy_mesh_->pmb_pack;
      if (prk->pmhd != nullptr && prk->pmhd->presist != nullptr &&
          prk->pmhd->presist->use_rkg_sts) {
        prk->pmhd->presist->track_x1_bflx = true;
      }
    }
    // problem/flux_hst_rkavg (see its note); read without recording it in the input
    const bool rkavg = pin->DoesParameterExist("problem","flux_hst_rkavg")
                       ? pin->GetBoolean("problem","flux_hst_rkavg") : true;
    const std::string integ = pin->GetOrAddString("time","integrator","rk2");
    dhj_fhst_nst = 0;
    if (integ == "rk1") {
      dhj_fhst_nst = 1;  dhj_fhst_w[0] = 1.0;
    } else if (integ == "rk2") {
      dhj_fhst_nst = 2;  dhj_fhst_w[0] = 0.5;  dhj_fhst_w[1] = 0.5;
    } else if (integ == "rk3") {
      dhj_fhst_nst = 3;  dhj_fhst_w[0] = 1.0/6.0;  dhj_fhst_w[1] = 1.0/6.0;
      dhj_fhst_w[2] = 2.0/3.0;
    }
    // other integrators (imex etc.): the last stage, as before
    dhj_fhst_avg = rkavg && user_srcs && (dhj_fhst_nst > 0);
    if (dhj_fhst_avg) {
      auto &ic = pmy_mesh_->mb_indcs;
      dhj_fhst_ptr = new DvceArray4D<Real>("dhj_fhst", pmy_mesh_->pmb_pack->nmb_thispack,
                                           4, ic.nx3 + 2*ic.ng, ic.nx2 + 2*ic.ng);
    }
  }
  // problem/flux_hst_wall (default false): four more columns after the nine flux
  // columns, splitting the inner-wall fluxes by the sign of the local mass flux:
  // Mdot_bot_in / Mdot_bot_out (sum over faces with fx1(IDN,is) > 0 / < 0) and
  // Efl_bot_in / Efl_bot_out (the fluid energy flux fx1(IEN,is) A over the same faces).
  dhj_flux_hst_wall = flux_hst && pin->DoesParameterExist("problem","flux_hst_wall")
                      && pin->GetBoolean("problem","flux_hst_wall");
  // problem/flux_hst_floor (default false; read without recording it, like flux_hst):
  // three more columns after the others, Efloor [erg/s] and Mfloor [g/s], the energy and
  // mass that ConToPrim's floors added to the active cells, averaged over the interval
  // since the previous history row, and Efloor_rt, the part of Efloor added by the
  // conversion right after the operator-split radiation (the rest is what the hydro
  // stages left below a floor).  Nothing else books them: they are not fluxes, so
  // they are the part of d(tot-E)/dt that Etot_bot - Etot_top does not contain (see
  // Hydro::c2p_track).  Costs two active-cell kernels per ConToPrim; off by default.
  // relaxdrag default (USER 11-01): false on a FRESH dhj start (recipe::NewDefaults; the
  // resolved value is written into the input, so restarts continue with it); a restart
  // whose input lacks the key keeps the old default true; an explicit key wins
  if (!pin->DoesParameterExist("problem","relax_drag")) {
    if (recipe::NewDefaults(pin)) {
      pin->GetOrAddBoolean("problem","relax_drag",false);
      recipe::Log("problem/relax_drag = false (default on a fresh dhj start)");
    } else {
      recipe::Log("problem/relax_drag = true (old default: restart input or recipe "
                  "none)");
    }
  }
  if (pin->DoesParameterExist("problem","relax_drag")) {
    dhj_relax_drag = pin->GetBoolean("problem","relax_drag");
  }
  if (pin->DoesParameterExist("problem","relax_drag_rot")) {
    dhj_relax_drag_rot = pin->GetReal("problem","relax_drag_rot");
  }
  if (global_variable::my_rank == 0 && (!dhj_relax_drag || dhj_relax_drag_rot != 5.0)) {
    std::cout << "dhj: start-up Rayleigh drag: problem/relax_drag = "
              << (dhj_relax_drag ? "true" : "false") << ", relax_drag_rot = "
              << dhj_relax_drag_rot << std::endl;
  }
  dhj_flux_hst_floor = flux_hst && pin->DoesParameterExist("problem","flux_hst_floor")
                       && pin->GetBoolean("problem","flux_hst_floor");
  if (dhj_flux_hst_floor) {
    MeshBlockPack *pfl = pmy_mesh_->pmb_pack;
    if (pfl->phydro != nullptr) pfl->phydro->EnableC2PTrack();
    if (pfl->pmhd != nullptr) pfl->pmhd->EnableC2PTrack();
    dhj_floor_t0 = pmy_mesh_->time;
  }
  // problem/coriolis_am, problem/am_transport (default ON in the rotating frame, user
  // 10-01) and problem/am_hst (DIAGNOSTIC, default false): see namespace amh.  All read
  // without recording, so a run that does not rotate is untouched (also its restarts).
  {
    MeshBlockPack *pam = pmy_mesh_->pmb_pack;
    const bool curv = pmy_mesh_->use_cubed_sphere || use_spherical_polar;
    const bool rotating = user_srcs && curv && (hot_jupiter_param.omega != 0.0);
    const bool mhd = (pam->pmhd != nullptr);
    const bool wbs = mhd ? pam->pmhd->use_wellbalance_static
                         : pam->phydro->use_wellbalance_static;
    Coordinates *pco = pam->pcoord;
    // the source variants AmCell replicates: the plain cs / sp angular sources.  Any
    // other one is named here (user 10-01: explicit, never a silent default off)
    std::string bad = "";
    if (wbs) bad = mhd ? "<mhd>/wellbalance_static" : "<hydro>/wellbalance_static";
    if (pmy_mesh_->use_cubed_sphere && pco->cs_wellbalanced_src) {
      bad = "cs_wellbalanced_src";
    }
    if (use_spherical_polar && pco->sp_wellbalanced_src) bad = "sp_wellbalanced_src";
    if (use_spherical_polar && pco->sp_cart_polar_momentum) {
      bad = "sp_cart_polar_momentum";
    }
    if (use_spherical_polar && pco->sp_face_avg) bad = "sp_face_avg";
    const bool trq_ok = curv && bad.empty();
    {
      const std::string sm = pin->DoesParameterExist("problem","sponge_top_mode")
                             ? pin->GetString("problem","sponge_top_mode") : "radial";
      if (sm != "all" && sm != "radial") {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "problem/sponge_top_mode must be all or radial, not '"
                  << sm << "'" << std::endl;
        exit(EXIT_FAILURE);
      }
      dhj_sponge_top_radial = (sm == "radial");
      if (pin->DoesParameterExist("problem","sponge_top_plo_bar")) {
        dhj_sp_top_plo_bar = pin->GetReal("problem","sponge_top_plo_bar");
      }
      if (pin->DoesParameterExist("problem","sponge_top_phi_bar")) {
        dhj_sp_top_phi_bar = pin->GetReal("problem","sponge_top_phi_bar");
      }
      if (pin->DoesParameterExist("problem","sponge_top_tau")) {
        dhj_sp_top_tau = pin->GetReal("problem","sponge_top_tau");
      }
      if (global_variable::my_rank == 0) {
        std::cout << "dhj: top sponge " << dhj_sp_top_plo_bar << " .. "
                  << dhj_sp_top_phi_bar << " bar, tau " << dhj_sp_top_tau << " s"
                  << std::endl;
      }
      if (global_variable::my_rank == 0) {
        std::cout << "dhj: problem/sponge_top_mode = " << sm
                  << (pin->DoesParameterExist("problem","sponge_top_mode") ? ""
                      : " (default; 'all' = the pre-10-01 full-velocity sponge)")
                  << std::endl;
      }
    }
    amh::fix = pin->DoesParameterExist("problem","coriolis_am")
               ? pin->GetBoolean("problem","coriolis_am") : rotating;
    amh::trans = pin->DoesParameterExist("problem","am_transport")
                 ? pin->GetBoolean("problem","am_transport") : rotating;
    // amdef-1019: the generic <hydro|mhd>/am_transport (default on cs / sp) is the same
    // transport correction: a DEFAULT one yields to problem/am_transport, a named one
    // together with it is a FATAL (applied twice)
    if (amh::trans) {
      if (hydro::AmModule *pg = hydro::AmModuleOf(pam)) {
        pg->YieldToPgen(pin, "problem/am_transport");
      }
    }
    if (rotating && amh::trans && !trq_ok) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "rotating frame with " << bad << ": problem/am_transport "
                << "does not replicate that angular source, so the axial angular "
                << "momentum is NOT conserved with it.  Set problem/am_transport = false "
                << "explicitly to run anyway (L_abs then drifts)" << std::endl;
      exit(EXIT_FAILURE);
    }
    if (global_variable::my_rank == 0 && rotating && !amh::trans) {
      std::cout << "### WARNING: problem/am_transport = false in a rotating run"
                << (bad.empty() ? std::string("") : " (" + bad + ")")
                << ": the transport angular-momentum correction is OFF, L_abs drifts"
                << std::endl;
    }
    if (global_variable::my_rank == 0 && rotating && !amh::fix) {
      std::cout << "### WARNING: problem/coriolis_am = false in a rotating run: the "
                << "Coriolis torque is not continuity-consistent, L_abs drifts"
                << std::endl;
    }
    amh::on = flux_hst && pin->DoesParameterExist("problem","am_hst")
              && pin->GetBoolean("problem","am_hst");
    amh::split = amh::on && pin->DoesParameterExist("problem","am_hst_split")
                 && pin->GetBoolean("problem","am_hst_split");
    if ((amh::fix || amh::trans || amh::on) && !(curv && user_srcs)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "problem/coriolis_am, am_transport and am_hst need "
                << "user_srcs and the cubed sphere or spherical polar" << std::endl;
      exit(EXIT_FAILURE);
    }
    if ((amh::trans || amh::split) && !trq_ok) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "problem/am_transport and am_hst_split replicate the "
                << "plain angular source only: not with cs/sp_wellbalanced_src, "
                << "wellbalance_static, sp_cart_polar_momentum or sp_face_avg"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    if (amh::on) {
      if (pin->DoesParameterExist("problem","am_hst_pdeep")) {
        amh::pdeep = pin->GetReal("problem","am_hst_pdeep");
      }
      const std::string integ = pin->GetOrAddString("time","integrator","rk2");
      amh::nst = 0;
      if (integ == "rk1") {
        amh::nst = 1;  amh::w[0] = 1.0;
      } else if (integ == "rk2") {
        amh::nst = 2;  amh::w[0] = 0.5;  amh::w[1] = 0.5;
      } else if (integ == "rk3") {
        amh::nst = 3;  amh::w[0] = 1.0/6.0;  amh::w[1] = 1.0/6.0;  amh::w[2] = 2.0/3.0;
      }
      // + Ttide (tidetrq-1101) with problem/roche
      const int nfh = 9 + (dhj_flux_hst_wall ? 4 : 0) + (dhj_flux_hst_floor ? 3 : 0)
                      + (amh::split ? 12 : 8) + (hot_jupiter_param.roche ? 1 : 0)
                      + 1;  // + Tdrag (relaxdrag-1101)
      if (amh::nst == 0 || nfh > NHISTORY_VARIABLES) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "problem/am_hst needs rk1/rk2/rk3 and at most "
                  << NHISTORY_VARIABLES << " user history columns (have " << nfh << ")"
                  << " (e.g. problem/am_hst_split = false frees 4)" << std::endl;
        exit(EXIT_FAILURE);
      }
      for (int n=0; n<amh::NACC; ++n) amh::acc[n] = 0.0;
      amh::t0 = pmy_mesh_->time;
    }
    if (global_variable::my_rank == 0 && (rotating || amh::fix || amh::trans)) {
      std::cout << "dhj: angular momentum: problem/coriolis_am = "
                << (amh::fix ? "true" : "false") << ", problem/am_transport = "
                << (amh::trans ? "true" : "false") << ", problem/am_hst = "
                << (amh::on ? "true" : "false") << std::endl;
    }
  }
  // <hydro>/lhllc_x1_phi_min, <mhd>/lhlld_x1_phi_min (EOS_Data::x1_phi_min) on dhj:
  //   FRESH start (recipe::NewDefaults; prod1103, USER 11-03): 0 = the low-Mach fix on
  //     the radial faces too, hydro AND MHD (mhdab_1102: MHD state within 0.5-2 % of 1,
  //     the r1m v6 deck's value), = the global default.
  //   RESTART whose input lacks the key (or <job>/recipe = none): 1, the old dhj-only
  //     default (user 09-30: the full HLLC/HLLD velocity-jump term on x1 faces), so
  //     running chains are unchanged.
  // An input (or restart) that names the key keeps its value; GetOrAdd records the
  // chosen value so this run's restarts carry it.  To change it on a restart use
  // `-i file` with the key (the command line cannot add a key a restart file lacks).
  {
    MeshBlockPack *plm = pmy_mesh_->pmb_pack;
    const bool nd = recipe::NewDefaults(pin);
    const Real xd = nd ? 0.0 : 1.0;
    const std::string why = nd ? " (default on a fresh start)" : (recipe::GetMode(pin) ==
        recipe::Mode::none ? " (recipe none: the old dhj default)" : " (restart input "
        "has no key: the old dhj default)");
    if (plm->phydro != nullptr) {
      if (!pin->DoesParameterExist("hydro", "lhllc_x1_phi_min")) {
        recipe::Log(std::string("<hydro>/lhllc_x1_phi_min = ") + (nd ? "0" : "1") + why);
      }
      plm->phydro->peos->eos_data.x1_phi_min =
          pin->GetOrAddReal("hydro", "lhllc_x1_phi_min", xd);
    }
    if (plm->pmhd != nullptr) {
      if (!pin->DoesParameterExist("mhd", "lhlld_x1_phi_min")) {
        recipe::Log(std::string("<mhd>/lhlld_x1_phi_min = ") + (nd ? "0" : "1") + why);
      }
      plm->pmhd->peos->eos_data.x1_phi_min =
          pin->GetOrAddReal("mhd", "lhlld_x1_phi_min", xd);
    }
    if (global_variable::my_rank == 0) {
      if (plm->phydro != nullptr) {
        std::cout << "dhj: <hydro>/lhllc_x1_phi_min = "
                  << plm->phydro->peos->eos_data.x1_phi_min << std::endl;
      }
      if (plm->pmhd != nullptr) {
        std::cout << "dhj: <mhd>/lhlld_x1_phi_min = "
                  << plm->pmhd->peos->eos_data.x1_phi_min << std::endl;
      }
    }
  }
  // problem/budget_dt > 0 (DIAGNOSTIC): the per-shell ledger (EBudCk).  Read only if the
  // key exists, so a run without it is untouched (also its restart files).
  ebud::dt_out = pin->DoesParameterExist("problem","budget_dt")
                 ? pin->GetReal("problem","budget_dt") : 0.0;
  ebud::on = (ebud::dt_out > 0.0);
  if (ebud::on) {
    MeshBlockPack *peb = pmy_mesh_->pmb_pack;
    if (peb->phydro == nullptr || peb->pmhd != nullptr
        || pin->GetString("time","integrator") != "rk2") {
      std::cout << "### FATAL ERROR: problem/budget_dt needs pure hydro and rk2"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    const int n1 = pmy_mesh_->mb_indcs.nx1;
    ebud::acc = new DvceArray2D<Real>("eb_acc", ebud::NQ, n1 + 1);
    ebud::tmp = new DvceArray2D<Real>("eb_tmp", 2, n1 + 1);
    ebud::fhst = new DvceArray2D<Real>("eb_fhst", 4, n1);
    {
      auto &ic = pmy_mesh_->mb_indcs;
      ebud::cpre = new DvceArray5D<Real>("eb_cpre", peb->nmb_thispack, 2,
                                         ic.nx3 + 2*ic.ng, ic.nx2 + 2*ic.ng,
                                         ic.nx1 + 2*ic.ng);
    }
    if (!peb->phydro->c2p_track) peb->phydro->EnableC2PTrack();
  }
  // problem/wall_closed (DEFAULT TRUE, user 09-27; false = the old reservoir wall): make
  // the inner radial boundary a CLOSED wall.
  // The user boundary (HydrostaticEquilibrium) leaves the inner ghosts at the state the
  // problem generator gave them -- the initial hydrostatic column at rest -- so the is
  // face was a Riemann problem against a FIXED reservoir: a first cell that moves at
  // v_r exchanges mass with it (in at the ghost's entropy, out at its own) and the face
  // carried a net fluid energy flux of 1.5e30 erg/s into WASP-121b w1x (rot 60-63),
  // 160 x sigma T_int^4 4 pi r_in^2, at zero mean mass flux.  With wall_closed the flux
  // at the is face uses the mirror of the interior-side state (as ix1_bc = reflect
  // does): zero mass flux, zero advected energy, the wall pressure from the interior.
  // The ghosts keep their role in the reconstruction slopes and, under MHD, the field
  // boundary is unchanged.  Recorded with GetOrAdd: a restart written before the key
  // existed picks up the closed wall.
  if (pin->GetOrAddBoolean("problem","wall_closed",true)) {
    if (pmy_mesh_->pmb_pack->phydro != nullptr) {
      pmy_mesh_->pmb_pack->phydro->wall_closed_ix1 = true;
    }
    if (pmy_mesh_->pmb_pack->pmhd != nullptr) {
      auto *pmhd_w = pmy_mesh_->pmb_pack->pmhd;
      pmhd_w->wall_closed_ix1 = true;
      // problem/wall_mhd_mirror (MHD only): conducting (DEFAULT, condwall-1028, user
      // 10-29), stress_free (wallsf-1028) or legacy (bitwise 03b468843).  The ghosts of
      // this boundary hold B_t ODD and B_r even (the field is normal to the wall,
      // HydrostaticEquilibrium below), while the legacy mirror kept
      // B_t EVEN at the is face: the wall then exerted the anchored Maxwell stress
      // -B_r B_t on a free-slip flow and did the work -B_r (v_t.B_t), +1.04 % of
      // Lsw_abs into WASP-121b b100 (wallpoynt_1028/RESULTS.md), which a wall at rest
      // cannot do.  stress_free flips B_t in the mirror (zero stress, zero Poynting
      // flux; v_t stays even, free slip as in hydro), takes the wall-edge corner EMFs
      // from the interior row and zeroes the resistive Poynting flux at the wall.
      const std::string wmm =
          pin->GetOrAddString("problem", "wall_mhd_mirror", "conducting");
      // conducting (condwall-1028): a perfectly conducting wall at rest.  The edge EMFs
      // tangential to the wall (ideal and resistive) are 0, so B_r on the wall is held;
      // the Riemann mass+energy fluxes through the wall are 0 (no Poynting flux); the
      // mirror is the legacy one (B_t even, v_t even: free slip), so the anchored field
      // keeps its Maxwell stress -B_r B_t on the wall (booked as Twall in the
      // am_hst_split ledger).
      // A restart whose input records wall_mhd_mirror keeps that value.
      if (wmm == "stress_free") {
        pmhd_w->wall_stress_free = true;
        if (pmhd_w->presist != nullptr) pmhd_w->presist->wall_bt_zero = true;
      } else if (wmm == "conducting") {
        pmhd_w->wall_conducting_ix1 = true;
        if (pmhd_w->presist != nullptr) pmhd_w->presist->wall_cond_ix1 = true;
      } else if (wmm != "legacy") {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "problem/wall_mhd_mirror = '" << wmm
                  << "'; expected stress_free, conducting or legacy" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (global_variable::my_rank == 0) {
        std::cout << "  inner radial wall (MHD): wall_mhd_mirror = " << wmm << std::endl;
      }
    }
    if (global_variable::my_rank == 0) {
      std::cout << "  inner radial wall: CLOSED (problem/wall_closed = true)"
                << std::endl;
    }
  }
  // problem/photosphere_dump = <file> writes the tau = 2/3 level per band and per
  // g-point at the end of the run, using the run's own correlated-k opacity.  Needs one
  // RT evaluation first (the per-cell tables are allocated lazily), so restart and take a
  // single cycle rather than running it at nlim = 0.
  if (!pin->GetOrAddString("problem", "photosphere_dump", "").empty()) {
    pgen_final_func = DhjPhotosphereDump;
  }
  // problem/olr_dump = <file> (default off) writes the net top-of-domain longwave flux
  // per ck band and column, plus the absorbed stellar flux, at the end of the run
  // (synthetic phase curves).  Same one-cycle-restart rule as photosphere_dump.
  if (!pin->GetOrAddString("problem", "olr_dump", "").empty()) {
    pgen_final_func = DhjFinalDumps;
  }
  // read before anything restart-sensitive: the outer BC needs it on restarts too
  bc_outer_maxwell = pin->GetOrAddBoolean("problem","bc_outer_maxwell",true);
  {
    // proddef-1101 (USER): fresh-start default outflow1 (Caltech, all five planets); a
    // restart whose input lacks the key keeps the old hse
    if (!pin->DoesParameterExist("problem","top_bc")) {
      recipe::Log(recipe::NewDefaults(pin) ? "problem/top_bc = outflow1 (default)"
                  : "problem/top_bc = hse (old default: restart input or recipe none)");
    }
    const std::string tb = pin->GetOrAddString("problem","top_bc",
                                               recipe::NewDefaults(pin) ? "outflow1"
                                                                        : "hse");
    if (tb == "hse") {
      dhj_top_bc = 0;
    } else if (tb == "outflow1") {
      dhj_top_bc = 1;
    } else if (tb == "outflow0") {
      dhj_top_bc = 2;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "problem/top_bc must be hse, outflow1 or outflow0, not '"
                << tb << "'" << std::endl;
      exit(EXIT_FAILURE);
    }
    if (global_variable::my_rank == 0 && dhj_top_bc != 0) {
      std::cout << "deep_hot_jupiter_rt: problem/top_bc = " << tb
                << " (zero-gradient outer-x1 ghosts, v_r = max(v_r_ie, 0)"
                << (dhj_top_bc == 1 ? "; one-way top face" : "; NOT one-way") << ")"
                << std::endl;
    }
    if (dhj_top_bc == 1) {
      if (pmy_mesh_->pmb_pack->phydro != nullptr) {
        pmy_mesh_->pmb_pack->phydro->diode_ox1 = true;
      }
      if (pmy_mesh_->pmb_pack->pmhd != nullptr) {
        pmy_mesh_->pmb_pack->pmhd->diode_ox1 = true;
      }
    }
  }
  if (global_variable::my_rank == 0) {
    std::cout << "deep_hot_jupiter_rt: outer-x1 Maxwell-stress term in the ghost "
              << "extrapolation is " << (bc_outer_maxwell ? "ON" : "off") << std::endl;
  }
  // blocked-band RT scaling harness: sized once, before any RT call.
  // problem/rt_split: DEFAULT TRUE since 2026-09-22.  The three-kernel chain-parallel
  // path is the same arithmetic in the same block-summation order as the monolithic
  // loop, so it is bitwise identical, and rt_ck (here) / rt_grey (box_convection,
  // red_giant) force it true anyway.  The serial path is kept -- the grey picket fence
  // with rt_ck = false is the one configuration that can still reach it -- because the
  // split grey sweep sizes its column at compile time and fatals when n1 > RT_NNC (72),
  // which the monolithic path does not -- and because on THAT path the two are NOT
  // bitwise equal (7.1 % in 1-KE after 20 cycles of deep_hot_jupiter_rt_ideal_xe; see
  // two_stream_rt.hpp).  A picket-fence run that wants the old answer must say false.
  rt_split = pin->GetOrAddBoolean("problem","rt_split",true);
  rt_ck = pin->GetOrAddBoolean("problem","rt_ck",false);
  // problem/rt_layer_legacy is read in full below; peek at it here because the
  // ck_spherical / ck_beam_sph defaults must not pick a value their own guard refuses.
  const bool layer_legacy_peek =
      pin->GetOrAddBoolean("problem","rt_layer_legacy",false);
  // problem/ck_spherical and problem/ck_beam_sph: the SPHERICAL thermal two-stream and
  // the PSEUDO-SPHERICAL direct beam of the correlated-k kernel (see
  // two_stream_rt::ck_spherical / ::ck_beam_sph).  They are pure GEOMETRY, so since
  // 2026-09-22 they DEFAULT TRUE on a curvilinear mesh -- pmesh->use_spherical_polar or
  // pmesh->use_cubed_sphere -- and are REFUSED on a Cartesian one, where the
  // plane-parallel form is the only one with a meaning.  The default is additionally
  // conditioned on the two things the flags' own guards in two_stream_rt.hpp require,
  // so that turning them on by default can never turn a working run into a fatal:
  // rt_ck (they convert the correlated-k kernel alone; the grey sweep is spherical
  // already, 8bca3dfa) and !rt_layer_legacy (the staggered whole-cell ck layers were
  // never converted).  Set either line explicitly to override.
  const bool cksph_default =
      (use_spherical_polar || use_cubed_sphere_) && rt_ck && !layer_legacy_peek;
  ck_spherical = pin->GetOrAddBoolean("problem","ck_spherical",cksph_default);
  ck_beam_sph = pin->GetOrAddBoolean("problem","ck_beam_sph",cksph_default);
  // problem/ck_sweep_cache (0 off / 1 per-pass carry / 2 per-column store): keep the
  // per-(cell, chain) layer operator in the thread instead of rebuilding it in each of
  // the four column passes.  Bitwise at every setting; the switch exists so the A/B can
  // be run on one binary.  Mode 2 pays private memory for it and is worth it only when
  // there ARE four passes, i.e. under ck_spherical -- so that is the default there and
  // mode 1 is the default otherwise.  See two_stream_rt::ck_sweep_cache and
  // tests_ck_sweep_cost/README.md.
  two_stream_rt::ck_sweep_cache =
      pin->GetOrAddInteger("problem","ck_sweep_cache", ck_spherical ? 2 : 1);
  if (two_stream_rt::ck_sweep_cache < 0 || two_stream_rt::ck_sweep_cache > 2) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_sweep_cache must "
              << "be 0, 1 or 2, got " << two_stream_rt::ck_sweep_cache << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // problem/ck_sweep_form (0 four-pass / 1 tm / 2 sd): solve the spherical face
  // coupling EXACTLY in two column passes instead of estimating it with two probe
  // passes.  See two_stream_rt::ck_sweep_form for the two recurrences, and
  // tests_ck_sweep_form/README.md for what the probe's lag was worth.  The remaining
  // refusals (ck_spherical off, ck_sweep_cache != 2, ck_implicit) are stated at the
  // first RT call, where those flags are all final.
  //
  // Since 2026-09-22 this DEFAULTS TO 1 (tm) wherever it can run, i.e. under
  // ck_spherical.  The four-pass form's probe carries a mixing constant that is right
  // only to O(beta^2) -- 2.3e-2 of the face flux at the production area ratio 3.45 --
  // while tm and sd solve the coupled column exactly (they agree to 3.5e-17) and tm is
  // 1.36x faster on the GPU; the numbers are in tests_ck_sweep_form/README.md.  The
  // default is conditioned on every one of the flag's own refusals, so that flipping it
  // cannot turn a working input into a startup fatal: it falls back to 0 with
  // ck_spherical off (plane-parallel: the two-pass sweep IS the exact form there) and at
  // ck_sweep_cache != 2 (the only setting the probe-free forms are instantiated at).
  // Until defaults-0927 (2026-09-27) it also fell back to 0 under ck_implicit; tm is now
  // the default there too (the T4 linear path needs it, tm carries its own tridiagonal
  // since phase T1), so an implicit input that does not name the key moves from the
  // four-pass to the tm form.  Set the line explicitly to override.
  const bool cksweepform_set = pin->DoesParameterExist("problem","ck_sweep_form");
  const int cksweepform_default =
      (ck_spherical && two_stream_rt::ck_sweep_cache == 2) ? 1 : 0;
  two_stream_rt::ck_sweep_form =
      pin->GetOrAddInteger("problem","ck_sweep_form", cksweepform_default);
  if (two_stream_rt::ck_sweep_form < 0 || two_stream_rt::ck_sweep_form > 2) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_sweep_form must "
              << "be 0 (four-pass), 1 (tm) or 2 (sd), got "
              << two_stream_rt::ck_sweep_form << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // State the RESOLVED value at startup, because the default is now conditional: an
  // input that says nothing gets the exact form on a spherical ck run and the historical
  // four-pass one everywhere else, and which of the two it got must not have to be
  // inferred from the input file.
  if (global_variable::my_rank == 0) {
    static const char *const ckform_name[3] = {"four-pass (probe)", "tm", "sd"};
    std::cout << "### deep_hot_jupiter_rt: problem/ck_sweep_form = "
              << two_stream_rt::ck_sweep_form << " ("
              << ckform_name[two_stream_rt::ck_sweep_form] << "), "
              << (cksweepform_set ? "set explicitly" :
                  (cksweepform_default == 1 ? "the default here (ck_spherical)" :
                   "the default here (no spherical face coupling to solve, or "
                   "ck_sweep_cache != 2)")) << std::endl;
  }
  // The refusal, stated here rather than only at the first RT call, so that a Cartesian
  // input naming either flag dies at startup with the reason.
  if ((ck_spherical || ck_beam_sph) &&
      !(use_spherical_polar || use_cubed_sphere_)) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_spherical and "
              << "problem/ck_beam_sph are the SPHERICAL forms of the correlated-k "
              << "thermal two-stream and of the direct stellar beam.  They have no "
              << "meaning on a Cartesian mesh (mesh/use_spherical_polar and "
              << "mesh/use_cubed_sphere are both false), where the plane-parallel form "
              << "is the geometry.  Remove the line, or set it false." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // problem/ck_implicit: the BACKWARD-EULER correlated-k column solve (see
  // two_stream_rt::ck_implicit and utils/two_stream_column_ck.hpp).  With it on the gas
  // receives exactly the converged flux divergence -- no relaxation factor, no rt_de_max
  // clip -- at ck_impl_maxit ck sweeps per RK stage instead of one; false = the
  // semi-implicit per-cell apply.  DEFAULT TRUE for FRESH runs since defaults2-0927
  // (user decision 09-27: the production scheme), and with it the T4 set, ck_impl_arat
  // 1e30 and ck_impl_nosync, whose defaults below are conditional on it.  The default
  // is conditioned on the switch's own refusal in two_stream_rt.hpp (rt_ck,
  // rt_use_cons, !rt_layer_legacy; rt_implicit_column is refused by this pgen anyway)
  // and on ck_sweep_form != 2 (sd has no tridiagonal), so it can never turn a working
  // input into a startup fatal.  On a RESTART the default stays false: the value is
  // recorded, so a file written by an older binary carries its "ck_implicit = 0" and
  // keeps it, and a file written before the key existed ran the semi-implicit apply,
  // which it keeps too.  Set the line explicitly to override.
  {const bool ckimp_default = !restart && rt_ck && !layer_legacy_peek &&
       pin->GetOrAddBoolean("problem","rt_use_cons",true) &&
       (two_stream_rt::ck_sweep_form != 2);
  two_stream_rt::ck_implicit =
      pin->GetOrAddBoolean("problem","ck_implicit",ckimp_default);}
  // problem/flux_hst_rst: default = flux_hst && ck_implicit; without ck_implicit an
  // explicit true is a no-op (no ck restart block to carry the cache) and warns
  two_stream_rt::rt_fhst_rst = flux_hst && two_stream_rt::ck_implicit &&
                               (fh_rst_req != 0);
  if (flux_hst && fh_rst_req == 1 && !two_stream_rt::ck_implicit &&
      global_variable::my_rank == 0) {
    std::cout << "### WARNING: problem/flux_hst_rst = true needs problem/ck_implicit; "
              << "ignored (restart files carry no flux-history cache)" << std::endl;
  }
  two_stream_rt::ck_impl_tol = pin->GetOrAddReal("problem","ck_impl_tol",1.0e-8);
  two_stream_rt::ck_impl_dtol = pin->GetOrAddReal("problem","ck_impl_dtol",1.0e-8);
  two_stream_rt::ck_impl_norm_eps =
      pin->GetOrAddReal("problem","ck_impl_norm_eps",1.0e-3);
  // default changed 09-30 (user), was 8
  // prod1103c default (USER 11-03 "default become 16"): ck_impl_maxit = 16 on a FRESH
  // start (recipe::NewDefaults) with ck_implicit when the key is not named (every Caltech
  // dhj production deck: "reset 09-30 (user): original 3x key"); recorded for restarts.
  // Otherwise 24 (the 09-30 default) as before; an explicit key wins.
  if (!pin->DoesParameterExist("problem","ck_impl_maxit") && two_stream_rt::ck_implicit) {
    const bool nd = recipe::NewDefaults(pin);
    pin->GetOrAddInteger("problem","ck_impl_maxit", nd ? 16 : 24);
    recipe::Log(nd ? "problem/ck_impl_maxit = 16 (default with ck_implicit)"
                   : "problem/ck_impl_maxit = 24 (restart input or recipe none: the old "
                     "default)");
  }
  two_stream_rt::ck_impl_maxit = pin->GetOrAddInteger("problem","ck_impl_maxit",24);
  two_stream_rt::ck_impl_dtmax = pin->GetOrAddReal("problem","ck_impl_dtmax",0.5);
  two_stream_rt::ck_impl_verbose =
      pin->GetOrAddBoolean("problem","ck_impl_verbose",false);
  two_stream_rt::ck_impl_demax = pin->GetOrAddReal("problem","ck_impl_demax",0.5);
  two_stream_rt::ck_impl_refresh_kappa =
      pin->GetOrAddBoolean("problem","ck_impl_refresh_kappa",false);
  // problem/ck_impl_tau_min: the two-level split (thin cells take the bracketed per-cell
  // implicit step, thick ones the tridiagonal).  This is what makes the Newton converge
  // at 10x and 100x the production dt; see utils/two_stream_column_ck.hpp.
  two_stream_rt::ck_impl_tau_min =
      pin->GetOrAddReal("problem","ck_impl_tau_min",0.0);
  // defaults-0927: 1e30 (no two-level split) under ck_implicit -- arat 2 FAILED (9 of 10
  // ranks non-converged, README_T1_stall) and every T4/c2 run uses 1e30; 2.0 otherwise
  two_stream_rt::ck_impl_arat = pin->GetOrAddReal("problem","ck_impl_arat",
      two_stream_rt::ck_implicit ? 1.0e30 : 2.0);
  two_stream_rt::ck_impl_debug = pin->GetOrAddInteger("problem","ck_impl_debug",0);
  two_stream_rt::ck_impl_ncloc = pin->GetOrAddInteger("problem","ck_impl_ncloc",0);
  two_stream_rt::ck_impl_ncref = pin->GetOrAddInteger("problem","ck_impl_ncref",0);
  dtloc_every = pin->GetOrAddInteger("problem","dtloc_every",0);
  two_stream_rt::ck_impl_colskip =
      pin->GetOrAddBoolean("problem","ck_impl_colskip",true);
  // ck_impl_once / _jreuse / _pred / _every: read here only when named; otherwise the
  // value (prod1103b fresh default or the old one) is set and recorded below, once the
  // fused-Newton configuration is known (before ck_impl_xstep)
  const bool once_named = pin->DoesParameterExist("problem","ck_impl_once");
  two_stream_rt::ck_impl_once =
      once_named ? pin->GetBoolean("problem","ck_impl_once") : false;
  // ck-fast levers (tests_ck_implicit/README_fast.md), all default off: the stored
  // operator kept across calls (xstep = k cycles, xstep_thr = max relative T/rho change),
  // the adaptive chord Jacobian (jreuse = the contraction that triggers a rebuild), no
  // confirmation pass (pred, pred_fac), and no allocation/blocking fill in the pass path
  // (nosync).  See utils/two_stream_column_ck.hpp.  ck_impl_xstep is read further down,
  // once ck_impl_every and ck_impl_frozen_op are known (its default depends on them).
  two_stream_rt::ck_impl_xstep_thr =
      pin->GetOrAddReal("problem","ck_impl_xstep_thr",0.0);
  const bool jreuse_named = pin->DoesParameterExist("problem","ck_impl_jreuse");
  two_stream_rt::ck_impl_jreuse =
      jreuse_named ? pin->GetReal("problem","ck_impl_jreuse") : 0.0;
  two_stream_rt::ck_impl_jreuse_xc =
      pin->GetOrAddBoolean("problem","ck_impl_jreuse_xc",false);
  two_stream_rt::ck_impl_jreuse_act =
      pin->GetOrAddReal("problem","ck_impl_jreuse_act",0.25);
  const bool pred_named = pin->DoesParameterExist("problem","ck_impl_pred");
  two_stream_rt::ck_impl_pred =
      pred_named ? pin->GetBoolean("problem","ck_impl_pred") : false;
  two_stream_rt::ck_impl_pred_fac =
      pin->GetOrAddReal("problem","ck_impl_pred_fac",0.5);
  two_stream_rt::ck_impl_pred_chk =
      pin->GetOrAddBoolean("problem","ck_impl_pred_chk",false);
  // defaults-0927: on under ck_implicit (bitwise on CPU and GPU, -4.3 %/cycle)
  two_stream_rt::ck_impl_nosync =
      pin->GetOrAddBoolean("problem","ck_impl_nosync",two_stream_rt::ck_implicit);
  two_stream_rt::ck_impl_warm_step =
      pin->GetOrAddBoolean("problem","ck_impl_warm_step",false);
  two_stream_rt::ck_impl_cvkeep =
      pin->GetOrAddBoolean("problem","ck_impl_cvkeep",false);
  // ck-cadence (tests_ck_implicit/README_cadence.md), default off: the full implicit
  // call every N cycles, the per-cell linearised implicit step in between, and the
  // per-column refresh guard.  See utils/two_stream_column_ck.hpp.
  const bool every_named = pin->DoesParameterExist("problem","ck_impl_every");
  two_stream_rt::ck_impl_every =
      every_named ? pin->GetInteger("problem","ck_impl_every") : 1;
  two_stream_rt::ck_impl_every_thr =
      pin->GetOrAddReal("problem","ck_impl_every_thr",0.0);
  // ck-fast2 lever 1, default off: the deep diffusion handover (per-layer dtau of every
  // chain >= ck_dif_dtau below the chain cut).  See utils/two_stream_column_ck.hpp.
  two_stream_rt::ck_dif_dtau = pin->GetOrAddReal("problem","ck_dif_dtau",0.0);
  two_stream_rt::ck_dif_margin = pin->GetOrAddInteger("problem","ck_dif_margin",0);
  // ck-fast2 lever 3 (problem/ck_beam_par) is read below, after the T4 set it needs.
  // the storing pass without the chain kernel (see two_stream_rt::ck_store_split)
  two_stream_rt::ck_store_split = pin->GetOrAddBoolean("problem","ck_store_split",true);
  // problem/ck_impl_frozen_op: freeze the exchange operator over the Newton passes (the
  // sweep is linear in B_b at frozen opacity, so only the opacity-dependent coefficients
  // have to be rebuilt -- and they do not change); ck_impl_frozen_cof decides whether the
  // layer coefficients are stored as well as kappa rho (memory against expm1); and
  // ck_impl_warm starts the Newton from the previous call's converged increment.  See
  // utils/two_stream_column_ck.hpp.
  // defaults-0927 (2026-09-27): the T4 set -- ck_impl_frozen_op, ck_impl_lin,
  // ck_impl_fuse, ck_impl_jac_lin, ck_impl_cvsec -- DEFAULTS ON under ck_implicit
  // wherever their own refusals are met (tm sweep on a spherical ck run, sweep cache 2,
  // no ck_impl_refresh_kappa, a double-precision chain, ck_impl_debug <= 0); they are
  // reformulations of the same Newton (same fixed point; lin vs chain 4e-15; see
  // docs/dev/switch_inventory_2026-09-24.md sect. 3 and tests_ck_implicit/README_T4.md).
  // Without ck_implicit they stay off, as before.  An explicit value always wins; the
  // values are recorded, so a restart keeps what its file stored.
  const bool t4def = two_stream_rt::ck_implicit && ck_spherical &&
                     two_stream_rt::ck_sweep_form == 1 &&
                     two_stream_rt::ck_sweep_cache == 2 &&
                     !two_stream_rt::ck_impl_refresh_kappa && (RT_FP32 == 0) &&
                     two_stream_rt::ck_impl_debug <= 0;
  two_stream_rt::ck_impl_frozen_op =
      pin->GetOrAddBoolean("problem","ck_impl_frozen_op",t4def);
  two_stream_rt::ck_impl_frozen_cof =
      pin->GetOrAddBoolean("problem","ck_impl_frozen_cof",true);
  two_stream_rt::ck_impl_warm =
      pin->GetOrAddBoolean("problem","ck_impl_warm",false);
  // problem/ck_impl_lin (+ ck_impl_lin_check): the linear re-apply kernel of the tm sweep
  // for the Newton passes that do not build the Jacobian.  See
  // utils/two_stream_column_ck.hpp.  Default off; the requirements are checked at the
  // first RT call.
  two_stream_rt::ck_impl_lin = pin->GetOrAddBoolean("problem","ck_impl_lin",
      t4def && two_stream_rt::ck_impl_frozen_op && two_stream_rt::ck_impl_frozen_cof);
  two_stream_rt::ck_impl_lin_check =
      pin->GetOrAddInteger("problem","ck_impl_lin_check",0);
  two_stream_rt::ck_impl_lin_thr = pin->GetOrAddInteger("problem","ck_impl_lin_thr",1);
  if (two_stream_rt::ck_impl_lin_thr != 1 && two_stream_rt::ck_impl_lin_thr != 4) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_impl_lin_thr "
              << "must be 1 or 4, got " << two_stream_rt::ck_impl_lin_thr << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // phase T4 (tests_ck_implicit/README_T4.md): the fused residual + step kernel, the
  // Jacobian from the linear kernel's stored factorisation, the Jacobian on pass 0, and
  // light-weight launches.  All default off.  See utils/two_stream_column_ck.hpp.
  two_stream_rt::ck_impl_fuse = pin->GetOrAddBoolean("problem","ck_impl_fuse",t4def);
  two_stream_rt::ck_impl_jac_lin = pin->GetOrAddBoolean("problem","ck_impl_jac_lin",
      t4def && two_stream_rt::ck_impl_lin && two_stream_rt::ck_impl_lin_thr == 1);
  // problem/ck_sph_face: see two_stream_rt::ck_sph_face.  problem/ck_sph_dilute = true
  // is the older spelling of ck_sph_face = 1.
  // f6-default (user 10-02): on a FRESH spherical ck run (rt_ck, ck_spherical, tm sweep)
  // the default is 6, the monotone VEF ladder, with its companions ck_sph_top = 1 and
  // (via ck_vef_on below) ck_vef_every = 150, ck_impl_rowfb = true, ck_impl_rsec = 4: the
  // validated WASP-121b production combination (ckf6_mono_1001, ckf6_cs_mhd_1001).
  // Face 6 needs the production linear Newton (f6_ok below).  When a key that face 6
  // cannot run with is set (or defaults off) and ck_sph_face is not named, the default
  // falls back to 0 (the historical Eddington face) with a WARNING; the regression
  // inputs (ck_nquad = 1, the semi-implicit apply) rely on that.  On a RESTART the
  // default stays 0: the value is recorded, so a file written by a binary that knew the
  // key keeps its face, and a file written before the key existed ran face 0.  Explicit
  // keys always win.
  two_stream_rt::ck_sph_dilute = pin->GetOrAddBoolean("problem","ck_sph_dilute",false);
  {const bool f6_cand = !restart && !two_stream_rt::ck_sph_dilute && rt_ck &&
       ck_spherical && two_stream_rt::ck_sweep_form == 1 &&
       !pin->DoesParameterExist("problem","ck_sph_face");
  bool f6_def = false;
  if (f6_cand) {
    // the requirements of face 6 (as f6_ok below), ck_sph_top = 1 unless named
    const int top_named = pin->DoesParameterExist("problem","ck_sph_top") ?
                          pin->GetInteger("problem","ck_sph_top") : 1;
    const int nq_named = pin->DoesParameterExist("problem","ck_nquad") ?
                         pin->GetInteger("problem","ck_nquad") : 2;
    std::string why;
    if (!two_stream_rt::ck_implicit) why += " ck_implicit = false;";
    if (!two_stream_rt::ck_impl_frozen_op) why += " ck_impl_frozen_op = false;";
    if (!two_stream_rt::ck_impl_lin) why += " ck_impl_lin = false;";
    if (!two_stream_rt::ck_impl_jac_lin) why += " ck_impl_jac_lin = false;";
    if (two_stream_rt::ck_impl_lin_thr != 1) why += " ck_impl_lin_thr != 1;";
    if (nq_named != 2) why += " ck_nquad != 2;";
    if (two_stream_rt::ck_dif_dtau > 0.0) why += " ck_dif_dtau > 0;";
    if (top_named != 1) why += " ck_sph_top != 1;";
    f6_def = why.empty();
    if (!f6_def && global_variable::my_rank == 0) {
      std::cout << "### WARNING in deep_hot_jupiter_rt: problem/ck_sph_face is not set "
                << "and its default (6) cannot run with" << why << " falling back to "
                << "ck_sph_face = 0 (Eddington face).  Name ck_sph_face to silence this."
                << std::endl;
    }
  }
  two_stream_rt::ck_sph_face = pin->GetOrAddInteger("problem","ck_sph_face",
      f6_def ? 6 : (two_stream_rt::ck_sph_dilute ? 1 : 0));}
  // problem/ck_sph_top: see two_stream_rt::ck_sph_top; default 1 (vacuum top) with face
  // 6, which needs it, else 0 (bitwise the historical datum)
  two_stream_rt::ck_sph_top = pin->GetOrAddInteger("problem","ck_sph_top",
      two_stream_rt::ck_sph_face == 6 ? 1 : 0);
  // ck_sph_face = 6 (the M-matrix form of the VEF closure) reads the same formal
  // solution, refresh, restart state and face-5 Newton defaults
  two_stream_rt::ck_vef_on = (two_stream_rt::ck_sph_face == 5 ||
                              two_stream_rt::ck_sph_face == 6);
  // problem/ck_wall_flux_exact (wallfix-1002, see two_stream_rt::ck_wall_flux_exact):
  // the face 5 / 6 wall delivers sigma T_int^4 exactly (was 1.0718x).  Default true on a
  // fresh run; a restart whose input lacks the key keeps the old wall (one line).
  {const bool wnamed = pin->DoesParameterExist("problem","ck_wall_flux_exact");
  two_stream_rt::ck_wall_flux_exact =
      pin->GetOrAddBoolean("problem","ck_wall_flux_exact",!restart);
  if (restart && !wnamed && two_stream_rt::ck_vef_on && global_variable::my_rank == 0) {
    std::cout << "deep_hot_jupiter_rt: restart input has no problem/ck_wall_flux_exact: "
              << "keeping the old face-" << two_stream_rt::ck_sph_face << " wall "
              << "(1.07x internal flux); set it true to switch" << std::endl;
  }}
  // problem/ck_vef_every, ck_vef_ncore: the formal solution of ck_sph_face = 5.  With
  // face 5 the DEFAULTS of ck_vef_every (150), ck_impl_rowfb (true) and ck_impl_rsec (4)
  // are the measured cheapest accuracy-neutral setting (cksph_test_0930/PORT.md);
  // explicit keys override; face 0 keeps the old defaults (bitwise)
  const bool vef5 = two_stream_rt::ck_vef_on;
  two_stream_rt::ck_vef_every = pin->GetOrAddInteger("problem","ck_vef_every",
                                                     vef5 ? 150 : 1);
  two_stream_rt::ck_vef_ncore = pin->GetOrAddInteger("problem","ck_vef_ncore",8);
  two_stream_rt::ck_vef_sync = pin->GetOrAddBoolean("problem","ck_vef_sync",true);
  // problem/ck_vef_report (diagnostic, print only; two_stream_rt.hpp)
  if (pin->DoesParameterExist("problem", "ck_vef_xlim")) {
    two_stream_rt::ck_vef_xlim = pin->GetReal("problem", "ck_vef_xlim");
  }
  if (pin->DoesParameterExist("problem", "ck_f6_report")) {
    two_stream_rt::ck_f6_report = pin->GetInteger("problem", "ck_f6_report");
  }
  if (pin->DoesParameterExist("problem", "ck_vef_report")) {
    two_stream_rt::ck_vef_report = pin->GetInteger("problem", "ck_vef_report");
  }
  // face 5 runs on the production linear path (ck_impl_lin, the paired tier, its
  // Jacobian from the factorisation); faces 1-4 on the chain kernel only
  const bool vef_lin_ok = two_stream_rt::ck_impl_lin_thr == 1 &&
      pin->GetOrAddInteger("problem","ck_nquad",2) == 2 &&
      two_stream_rt::ck_dif_dtau <= 0.0;
  const bool lin_on = two_stream_rt::ck_impl_lin || two_stream_rt::ck_impl_jac_lin;
  // face 6 lives in the linear kernels only (ck_f6_build, rt_chain_ck_lin6,
  // rt_chain_ck_jlin6): the frozen-operator Newton with ck_impl_lin and
  // ck_impl_jac_lin, the paired tier, vacuum top
  const bool f6_ok = two_stream_rt::ck_implicit && two_stream_rt::ck_impl_frozen_op &&
      two_stream_rt::ck_impl_lin && two_stream_rt::ck_impl_jac_lin && vef_lin_ok &&
      two_stream_rt::ck_sph_top == 1;
  if (two_stream_rt::ck_sph_face == 6 && !f6_ok) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_sph_face = 6 needs "
              << "ck_implicit, ck_impl_frozen_op, ck_impl_lin, ck_impl_jac_lin, "
              << "ck_impl_lin_thr = 1, ck_nquad = 2, ck_dif_dtau = 0 and ck_sph_top = 1"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (two_stream_rt::ck_sph_face < 0 || two_stream_rt::ck_sph_face > 6 ||
      (two_stream_rt::ck_sph_face != 0 &&
       (!two_stream_rt::ck_spherical || two_stream_rt::ck_sweep_form != 1)) ||
      (two_stream_rt::ck_sph_face >= 1 && two_stream_rt::ck_sph_face <= 4 && lin_on) ||
      (two_stream_rt::ck_sph_face >= 5 && lin_on && !vef_lin_ok) ||
      (two_stream_rt::ck_sph_face >= 3 && two_stream_rt::ck_dif_dtau > 0.0) ||
      (two_stream_rt::ck_sph_face >= 5 && two_stream_rt::ck_sph_top == 2)) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_sph_face must be "
              << "0..6, and != 0 needs ck_spherical and ck_sweep_form = 1; 1..4 need "
              << "ck_impl_lin = false and ck_impl_jac_lin = false; 5 with ck_impl_lin "
              << "needs ck_impl_lin_thr = 1, ck_nquad = 2 and ck_dif_dtau = 0; 3 and 4 "
              << "also need ck_dif_dtau = 0, 5 ck_sph_top != 2"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // problem/ck_sph_top (prototype, default 0): see two_stream_rt::ck_sph_top
  if (two_stream_rt::ck_sph_top < 0 || two_stream_rt::ck_sph_top > 2 ||
      (two_stream_rt::ck_sph_top == 2 && lin_on)) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_sph_top must be "
              << "0..2, and 2 needs ck_impl_lin = false and ck_impl_jac_lin = false"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  two_stream_rt::ck_impl_jac0 = pin->GetOrAddBoolean("problem","ck_impl_jac0",false);
  two_stream_rt::ck_impl_jneg = pin->GetOrAddBoolean("problem","ck_impl_jneg",false);
  two_stream_rt::ck_impl_cvsec = pin->GetOrAddBoolean("problem","ck_impl_cvsec",
      t4def && two_stream_rt::ck_impl_fuse);
  // problem/ck_beam_par (ck-fast2 lever 3): the storing pass's pseudo-spherical beam in
  // its own per-target kernels (see two_stream_rt::ck_beam_par).  defaults-0927: ON
  // wherever its requirements (checked below) hold, off otherwise (ck-fast2: -25 % RT;
  // it differs from the in-kernel beam only by deposits below e^-60 of the beam).
  {const bool bpdef = two_stream_rt::ck_implicit && two_stream_rt::ck_impl_frozen_op &&
                      !two_stream_rt::ck_impl_refresh_kappa &&
                      two_stream_rt::ck_impl_lin &&
                      two_stream_rt::ck_impl_lin_thr == 1 &&
                      two_stream_rt::ck_sweep_form == 1 && ck_spherical && ck_beam_sph;
  two_stream_rt::ck_beam_par = pin->GetOrAddBoolean("problem","ck_beam_par",bpdef);}
  // problem/ck_impl_rsec: secant bound on the thick rows' diagonal (0 = off, bitwise)
  two_stream_rt::ck_impl_rsec = pin->GetOrAddReal("problem","ck_impl_rsec",
                                                  two_stream_rt::ck_vef_on ? 4.0 : 0.0);
  two_stream_rt::ck_impl_lw = pin->GetOrAddBoolean("problem","ck_impl_lw",false);
  // problem/ck_impl_glob: Newton globalisation of the fused T4 step, none | ls | ls_sub
  // (utils/two_stream_column_ck.hpp, tests_ck_implicit/README_glob.md).  none = bitwise.
  {
    const std::string gl = pin->GetOrAddString("problem","ck_impl_glob","none");
    if (gl == "none") {
      two_stream_rt::ck_impl_glob = 0;
    } else if (gl == "ls") {
      two_stream_rt::ck_impl_glob = 1;
    } else if (gl == "ls_sub") {
      two_stream_rt::ck_impl_glob = 2;
    } else {
      std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_impl_glob = " << gl
                << " (none | ls | ls_sub)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    two_stream_rt::ck_impl_ls_ntry =
        pin->GetOrAddInteger("problem","ck_impl_ls_ntry",4);
    two_stream_rt::ck_impl_ls_c = pin->GetOrAddReal("problem","ck_impl_ls_c",1.0e-4);
    two_stream_rt::ck_impl_sub_max =
        pin->GetOrAddInteger("problem","ck_impl_sub_max",32);
  }
  // problem/ck_impl_esc (cheaper ls_sub escalation) and problem/ck_impl_aa (Anderson
  // acceleration of the chord iteration, depth), both default off = bitwise
  // (utils/two_stream_column_ck.hpp, tests_ck_implicit/README_jac.md)
  two_stream_rt::ck_impl_esc = pin->GetOrAddInteger("problem","ck_impl_esc",0);
  two_stream_rt::ck_impl_esc_rho = pin->GetOrAddReal("problem","ck_impl_esc_rho",0.9);
  two_stream_rt::ck_impl_esc_extra =
      pin->GetOrAddInteger("problem","ck_impl_esc_extra",8);
  two_stream_rt::ck_impl_aa = pin->GetOrAddInteger("problem","ck_impl_aa",0);
  two_stream_rt::ck_impl_aa_rst = pin->GetOrAddBoolean("problem","ck_impl_aa_rst",true);
  // problem/ck_impl_floorbound / ck_impl_kkt_demax: the Newton as a bound-constrained
  // solve (e >= e_floor; KKT cells on the floor / on the demax bound count as converged).
  // See utils/two_stream_column_ck.hpp.  DEFAULT TRUE since defaults2-0927 on FRESH runs
  // of the fused ck_implicit Newton (ck_implicit, ck_impl_fuse, ck_impl_debug <= 0,
  // ck_impl_glob = none, ck_impl_aa = 0: the only path they exist on); false elsewhere
  // and on restarts (a file that lacks the keys keeps the unconstrained Newton).  Where
  // no cell reaches a bound both are the identity: bitwise on vs off on the CPU gates
  // (dhj_ck_spherical and the 64x8x8 gate input, 20 cycles, floor = kkt = 0).  Where
  // they bind (both WASP-121b productions, every call) the unconstrained Newton cannot
  // converge: 1x 466/500 calls NOT-CONVERGED off vs 6/500 on, 10x 500/500 vs 133/500
  // (defaults2_0927/B).  Set either line explicitly to override.
  {const bool fbdef = !restart && two_stream_rt::ck_implicit &&
       two_stream_rt::ck_impl_fuse && two_stream_rt::ck_impl_debug <= 0 &&
       two_stream_rt::ck_impl_glob == 0 && two_stream_rt::ck_impl_aa == 0;
  two_stream_rt::ck_impl_floorbound =
      pin->GetOrAddBoolean("problem","ck_impl_floorbound",fbdef);
  two_stream_rt::ck_impl_kkt_demax =
      pin->GetOrAddBoolean("problem","ck_impl_kkt_demax",fbdef);}
  two_stream_rt::ck_impl_stalldbg =
      pin->GetOrAddInteger("problem","ck_impl_stalldbg",0);
  two_stream_rt::ck_impl_negdbg = pin->GetOrAddInteger("problem","ck_impl_negdbg",0);
  two_stream_rt::ck_impl_jfd = pin->GetOrAddInteger("problem","ck_impl_jfd",0);
  two_stream_rt::ck_impl_rowfb = pin->GetOrAddBoolean("problem","ck_impl_rowfb",
                                                     two_stream_rt::ck_vef_on);
  if (two_stream_rt::ck_vef_on && global_variable::my_rank == 0) {
    std::cout << "deep_hot_jupiter_rt: ck_sph_face = " << two_stream_rt::ck_sph_face
              << " (VEF), effective: ck_sph_top = "
              << two_stream_rt::ck_sph_top << ", ck_vef_every = "
              << two_stream_rt::ck_vef_every << ", ck_vef_sync = "
              << two_stream_rt::ck_vef_sync << ", ck_impl_rowfb = "
              << two_stream_rt::ck_impl_rowfb << ", ck_impl_rsec = "
              << two_stream_rt::ck_impl_rsec << ", ck_wall_flux_exact = "
              << two_stream_rt::ck_wall_flux_exact << std::endl;
  }
  two_stream_rt::ck_impl_negpiv = pin->GetOrAddBoolean("problem","ck_impl_negpiv",false);
  two_stream_rt::ck_impl_jnet = pin->GetOrAddBoolean("problem","ck_impl_jnet",false);
  two_stream_rt::ck_impl_tol_ptop = pin->GetOrAddReal("problem","ck_impl_tol_ptop",0.0);
  two_stream_rt::ck_impl_norm_eps_top =
      pin->GetOrAddReal("problem","ck_impl_norm_eps_top",0.1);
  if ((two_stream_rt::ck_impl_rowfb || two_stream_rt::ck_impl_negpiv) &&
      !two_stream_rt::ck_impl_fuse) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_impl_rowfb and "
              << "ck_impl_negpiv live in the fused step (ck_impl_fuse)." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // problem/ck_impl_kkt_row: DEFAULT TRUE since default-flips (09-26; 1x production
  // NOT-CONVERGED 20/2113 -> 0 at the same cost).  It acts only with ck_impl_floorbound
  // or ck_impl_kkt_demax on the fused glob = 0 path; false restores the old rows.  The
  // value is recorded, so a restart keeps its stored value; a restart file written
  // before the key existed picks up true.
  two_stream_rt::ck_impl_kkt_row =
      pin->GetOrAddBoolean("problem","ck_impl_kkt_row",true);
  // problem/ck_impl_conserve: 0 = off (bitwise), 1 (default) = the reported face fluxes
  // (Lir_top, the ledger) consistent with the energy applied in the cells a ck_implicit
  // call leaves on a bound, 2 = in every cell.  The solution is unchanged either way
  // (CkConserveFlux, utils/two_stream_rt.hpp).  Read only if present, so a run without
  // it writes the same restart header.
  // default changed 09-30 (user), was 0
  // prod1103b default (USER 10-28 production key; audit 11-03): ck_impl_conserve = 2 on a
  // FRESH start (recipe::NewDefaults) with ck_implicit when the key is not named
  // (diagnostic: the solution is unchanged); recorded for the run's restarts.
  if (!pin->DoesParameterExist("problem","ck_impl_conserve") && two_stream_rt::ck_implicit &&
      recipe::NewDefaults(pin)) {
    pin->GetOrAddInteger("problem","ck_impl_conserve",2);
    recipe::Log("problem/ck_impl_conserve = 2 (default with ck_implicit)");
  }
  two_stream_rt::ck_impl_conserve =
      pin->DoesParameterExist("problem","ck_impl_conserve")
      ? pin->GetInteger("problem","ck_impl_conserve") : 1;
  if (global_variable::my_rank == 0) {
    std::cout << "deep_hot_jupiter_rt: ck_impl_maxit " << two_stream_rt::ck_impl_maxit
              << ", ck_impl_conserve " << two_stream_rt::ck_impl_conserve
              << ", sponge_top " << hot_jupiter_param.sponge_top
              << ", sponge_bottom " << hot_jupiter_param.sponge_bottom
              << " (defaults 24, 1, 1, 0 since 09-30)" << std::endl;
  }
  // problem/ck_impl_osc = N: per-cell Aitken damping of oscillating cells from pass N
  // on (utils/two_stream_column_ck.hpp; cknewton_0928).  0 = off (default, bitwise).
  two_stream_rt::ck_impl_osc = pin->GetOrAddInteger("problem","ck_impl_osc",0);
  if (two_stream_rt::ck_impl_osc < 0 || (two_stream_rt::ck_impl_osc > 0 &&
      two_stream_rt::ck_impl_aa > 0)) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_impl_osc must be "
              << ">= 0, and 0 with ck_impl_aa > 0; got " << two_stream_rt::ck_impl_osc
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if ((two_stream_rt::ck_impl_floorbound || two_stream_rt::ck_impl_kkt_demax ||
       two_stream_rt::ck_impl_rsec > 0.0 || two_stream_rt::ck_impl_osc > 0) &&
      (!two_stream_rt::ck_implicit || !two_stream_rt::ck_impl_fuse ||
       two_stream_rt::ck_impl_debug > 0 || two_stream_rt::ck_impl_glob != 0)) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_impl_floorbound, "
              << "ck_impl_kkt_demax, ck_impl_rsec and ck_impl_osc are implemented on the "
              << "fused "
              << "ck_implicit Newton only (ck_implicit, ck_impl_fuse, "
              << "ck_impl_debug <= 0, "
              << "ck_impl_glob = none)." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // problem/ck_impl_reuse_jac (0 off / 1 chord / 2 scaled chord) and problem/ck_impl_seed
  // (0 off / 1 semi-implicit / 2 exact lagged): the two phase-4 levers, both default off.
  // See tests_ck_implicit/README_phase4.md.
  two_stream_rt::ck_impl_reuse_jac =
      pin->GetOrAddInteger("problem","ck_impl_reuse_jac",0);
  two_stream_rt::ck_impl_seed = pin->GetOrAddInteger("problem","ck_impl_seed",0);
  if (two_stream_rt::ck_impl_reuse_jac < 0 || two_stream_rt::ck_impl_reuse_jac > 2) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_impl_reuse_jac "
              << "must be 0, 1 or 2, got " << two_stream_rt::ck_impl_reuse_jac
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (two_stream_rt::ck_impl_seed < 0 || two_stream_rt::ck_impl_seed > 2) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_impl_seed must be "
              << "0, 1 or 2, got " << two_stream_rt::ck_impl_seed << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (two_stream_rt::ck_impl_frozen_op && two_stream_rt::ck_impl_refresh_kappa &&
      global_variable::my_rank == 0) {
    std::cout << "### WARNING in deep_hot_jupiter_rt: problem/ck_impl_frozen_op is "
              << "IGNORED with ck_impl_refresh_kappa = true: the operator is frozen "
              << "only where the opacity is." << std::endl;
  }
  // prod1103b defaults (USER 09-24 dhj production "T4 + C2 + E4", approved for every
  // Caltech dhj production; audit 11-03): on a FRESH start (recipe::NewDefaults) with the
  // fused frozen-operator ck_implicit Newton (ck_impl_fuse, ck_impl_frozen_op, debug <= 0,
  // glob = none, aa = 0, reuse_jac = 0, colskip, no warm start, no rt_rad_force), each of
  // these keys that is not named takes the production value and is recorded (GetOrAdd):
  //   ck_impl_once = true (C2: the ck RT once per cycle, operator split after the RK
  //   stages), ck_impl_jreuse = 0.2 (C2 adaptive chord), ck_impl_pred = true (C2
  //   predictor), ck_impl_every = 4 (E4 cadence; only with ck_impl_once on).
  // Otherwise (restart input without the key, recipe none, another ck configuration) the
  // old default (false / 0 / false / 1) is recorded as before.  An explicit key wins.
  // ck_impl_xstep is NOT defaulted (audit_dhj_0928, below).
  {
    const bool nd = recipe::NewDefaults(pin);
    const bool c2ok = two_stream_rt::ck_implicit && two_stream_rt::ck_impl_fuse &&
        two_stream_rt::ck_impl_frozen_op && two_stream_rt::ck_impl_debug <= 0 &&
        two_stream_rt::ck_impl_glob == 0 && two_stream_rt::ck_impl_aa == 0 &&
        two_stream_rt::ck_impl_reuse_jac == 0 && two_stream_rt::ck_impl_colskip &&
        !two_stream_rt::ck_impl_warm && !two_stream_rt::rt_rad_force;
    const bool c2def = nd && c2ok;
    const std::string why = c2def ? " (default with the fused T4 ck_implicit)" :
        (!nd ? (recipe::GetMode(pin) == recipe::Mode::none ? " (recipe none: the old "
        "default)" : " (restart input has no key: the old default)") : "");
    const bool logit = two_stream_rt::ck_implicit && (c2def || !nd);
    if (!once_named) {
      two_stream_rt::ck_impl_once = pin->GetOrAddBoolean("problem","ck_impl_once",c2def);
      if (logit) recipe::Log(std::string("problem/ck_impl_once = ") +
                             (c2def ? "true" : "false") + why);
    }
    if (!jreuse_named) {
      two_stream_rt::ck_impl_jreuse =
          pin->GetOrAddReal("problem","ck_impl_jreuse",c2def ? 0.2 : 0.0);
      if (logit) recipe::Log(std::string("problem/ck_impl_jreuse = ") +
                             (c2def ? "0.2" : "0") + why);
    }
    if (!pred_named) {
      two_stream_rt::ck_impl_pred = pin->GetOrAddBoolean("problem","ck_impl_pred",c2def);
      if (logit) recipe::Log(std::string("problem/ck_impl_pred = ") +
                             (c2def ? "true" : "false") + why);
    }
    if (!every_named) {
      const bool e4 = c2def && two_stream_rt::ck_impl_once;
      two_stream_rt::ck_impl_every =
          pin->GetOrAddInteger("problem","ck_impl_every",e4 ? 4 : 1);
      if (logit) recipe::Log(std::string("problem/ck_impl_every = ") + (e4 ? "4" : "1") +
                             (e4 || !nd ? why : " (ck_impl_once off)"));
    }
  }
  // problem/ck_impl_xstep: keep the stored operator of ck_impl_frozen_op over k CYCLES
  // (under ck_impl_every = 4, k <= 4 is bitwise off).  DEFAULT 0 (off) since
  // audit_dhj_0928: on MI300A (WASP-121b nx1 76/74, current tip) k = 8 gains only
  // 3.2 % (1x) / 3.9 % (10x) per simulated second, and its deviation reaches 1.3-2x the
  // 1e-14-kick spread at 1e-3..1 bar and > 1 bar within the first rotation (it settles
  // into the noise by ~1 rotation), so it is a RELAXATION-PHASE option only.  Set it
  // explicitly (8) for spin-up; the value is recorded, so a restart keeps its value.
  two_stream_rt::ck_impl_xstep =
      pin->GetOrAddInteger("problem","ck_impl_xstep",0);
  if (two_stream_rt::ck_impl_xstep > 0 && !(two_stream_rt::ck_implicit &&
      two_stream_rt::ck_impl_frozen_op && !two_stream_rt::ck_impl_refresh_kappa)) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_impl_xstep keeps "
              << "the stored operator of ck_impl_frozen_op across calls; it needs "
              << "ck_implicit and ck_impl_frozen_op (and not ck_impl_refresh_kappa)."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (two_stream_rt::ck_impl_jreuse > 0.0 && (two_stream_rt::ck_impl_reuse_jac != 0 ||
      two_stream_rt::ck_impl_glob != 0)) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_impl_jreuse is the "
              << "adaptive chord; it replaces ck_impl_reuse_jac (set it 0) and is not "
              << "combined with ck_impl_glob." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (two_stream_rt::ck_impl_cvkeep && !two_stream_rt::ck_impl_cvsec) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_impl_cvkeep keeps "
              << "the secant cv of ck_impl_cvsec across calls; set ck_impl_cvsec."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (two_stream_rt::ck_impl_pred && (!two_stream_rt::ck_impl_fuse ||
      two_stream_rt::ck_impl_debug > 0 || two_stream_rt::ck_impl_glob != 0 ||
      two_stream_rt::ck_impl_aa > 0)) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_impl_pred needs "
              << "ck_impl_fuse (ck_impl_debug <= 0) and neither ck_impl_glob nor "
              << "ck_impl_aa." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (two_stream_rt::ck_impl_every > 1 && !(two_stream_rt::ck_implicit &&
      two_stream_rt::ck_impl_once && two_stream_rt::ck_impl_fuse &&
      two_stream_rt::ck_impl_colskip && two_stream_rt::ck_impl_debug <= 0 &&
      two_stream_rt::ck_impl_glob == 0 && !two_stream_rt::ck_impl_warm)) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_impl_every > 1 "
              << "needs ck_implicit, ck_impl_once, ck_impl_fuse, ck_impl_colskip, "
              << "ck_impl_debug <= 0, ck_impl_glob = 0 and ck_impl_warm off."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (two_stream_rt::ck_beam_par && !(two_stream_rt::ck_implicit &&
      two_stream_rt::ck_impl_frozen_op && !two_stream_rt::ck_impl_refresh_kappa &&
      two_stream_rt::ck_impl_lin && two_stream_rt::ck_impl_lin_thr == 1 &&
      two_stream_rt::ck_sweep_form == 1 && two_stream_rt::ck_spherical &&
      two_stream_rt::ck_beam_sph)) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_beam_par needs "
              << "ck_implicit, ck_impl_frozen_op, ck_impl_lin, ck_impl_lin_thr = 1, "
              << "ck_sweep_form = 1, ck_spherical and ck_beam_sph." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (two_stream_rt::ck_dif_dtau > 0.0 && !(two_stream_rt::ck_implicit &&
      two_stream_rt::ck_impl_frozen_op && !two_stream_rt::ck_impl_refresh_kappa &&
      two_stream_rt::ck_impl_lin && two_stream_rt::ck_impl_jac_lin &&
      two_stream_rt::ck_sweep_form == 1 && two_stream_rt::ck_spherical &&
      two_stream_rt::ck_impl_debug <= 0 && two_stream_rt::ck_dif_margin >= 0)) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_dif_dtau (the deep "
              << "diffusion handover) needs ck_implicit, ck_impl_frozen_op, ck_impl_lin, "
              << "ck_impl_jac_lin, ck_sweep_form = 1, ck_spherical, ck_impl_debug <= 0 "
              << "and ck_dif_margin >= 0." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (two_stream_rt::ck_impl_once && !two_stream_rt::ck_implicit) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_impl_once moves "
              << "the WHOLE radiation operator out of the RK stage and is only meant "
              << "for the implicit column solve; set problem/ck_implicit = true."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (two_stream_rt::ck_impl_once && two_stream_rt::rt_rad_force) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_impl_once is not "
              << "compatible with problem/rt_rad_force -- the force is applied inside "
              << "the RK stage from the sweep's arrays, which the split step no longer "
              << "refreshes there." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (two_stream_rt::ck_impl_once) {
    // once per hydro step, with the full cycle dt, after the last RK stage.  The
    // machinery is the generic operator split of hydro_tasks.cpp (user_split_func with
    // user_split_once), the same one box_convection's rt_once_per_cycle used.
    user_split_func = DhjCkRtSplit;
    user_split_once = true;
  }
  // rkgfuse-1102 problem/ck_strang (default false, read only when named): the ck RT of
  // ck_impl_once Strang-split around the RK stages, ck(dt/2) H ck(dt/2), instead of one
  // full-dt call after them (first order).  Under MHD the first half is taken by the
  // Driver after the leading RKG super-step (user_split_pre_late), so with <mhd>/
  // rkg_strang the RKG halves nest outside: R(dt/2) ck(dt/2) H ck(dt/2) R(dt/2).  The
  // ck_impl_every cadence counts cycles (two_stream_column_ck.hpp, ck_cad_half).
  if (pin->DoesParameterExist("problem", "ck_strang") &&
      pin->GetBoolean("problem", "ck_strang")) {
    if (!two_stream_rt::ck_impl_once) {
      std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_strang needs "
                << "problem/ck_impl_once = true" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (ebud::on) {
      std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_strang is not "
                << "supported with the energy-budget diagnostic (ebud, one call per "
                << "cycle)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    user_split_once = false;
    user_split_pre_late = (pmy_mesh_->pmb_pack->pmhd != nullptr);
  }
  // problem/rt_test_* (see the note at rt_test_freeze).  Default off.
  rt_test_freeze = pin->GetOrAddBoolean("problem","rt_test_freeze",false);
  if (rt_test_freeze) {
    if (pmy_mesh_->pmb_pack->pmhd != nullptr || global_variable::nranks != 1) {
      std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/rt_test_freeze is "
                << "a serial hydro-only test mode" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    rt_test_dt = pin->GetOrAddReal("problem","rt_test_dt",0.0);
    rt_test_out = pin->GetOrAddString("problem","rt_test_out","");
    rt_test_every = pin->GetOrAddInteger("problem","rt_test_every",1);
    rt_test_col[0] = pin->GetOrAddInteger("problem","rt_test_col_m",0);
    rt_test_col[1] = pin->GetOrAddInteger("problem","rt_test_col_k",-1);
    rt_test_col[2] = pin->GetOrAddInteger("problem","rt_test_col_j",-1);
    user_split_func = DhjRtFreezeSplit;
    user_split_once = true;
  }
  // ckhesdirk-1102 problem/ck_time_scheme (read only when named): split (default: the
  // operator split above, bitwise) | hesdirk2 (the ck heating as the implicit part of
  // the IMEX H-ESDIRK2 step around the Heun stages, in-stage; see DhjCkT2Stage)
  if (pin->DoesParameterExist("problem","ck_time_scheme")) {
    const std::string cts = pin->GetString("problem","ck_time_scheme");
    if (cts.compare("hesdirk2") == 0) {
      two_stream_rt::ck_t2_on = true;
    } else if (cts.compare("split") != 0) {
      std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_time_scheme = '"
                << cts << "' is not a choice (split | hesdirk2)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  if (pin->DoesParameterExist("problem","ck_count_every")) {
    dhj_ck_count_every = pin->GetInteger("problem","ck_count_every");
  }
  // ck2nd-1102: ck_impl_tol_inc, ck_cad_drift (read only when named; utils/ck_2nd.hpp)
  two_stream_rt::Ck2ndInit(pin, two_stream_rt::ck_impl_fuse &&
                           (two_stream_rt::ck_impl_debug <= 0),
                           two_stream_rt::ck_impl_glob, two_stream_rt::ck_impl_every);
  if (two_stream_rt::ck_t2_on) {
    std::string why;
    if (!two_stream_rt::ck_implicit) why += " ck_implicit = false;";
    if (pin->GetString("time","integrator").compare("rk2") != 0) why += " integrator != rk2;";
    if (rt_test_freeze) why += " rt_test_freeze;";
    if (ebud::on) why += " budget_dt;";
    if (two_stream_rt::rt_rad_force) why += " rt_rad_force;";
    if (two_stream_rt::ck_impl_warm) why += " ck_impl_warm;";
    // prod-1102: ck_strang (operator-split ck halves) and the in-stage hesdirk2 ck are
    // two different placements of the same ck heating; together the ck would be lost
    // or doubled
    if (pin->DoesParameterExist("problem", "ck_strang") &&
        pin->GetBoolean("problem", "ck_strang")) why += " ck_strang;";
    if (!why.empty()) {
      std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_time_scheme = "
                << "hesdirk2 is not wired for:" << why << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // the ck leaves the operator split: it runs inside each Heun stage (SourceFunc)
    user_split_func = nullptr;
    user_split_once = false;
    DhjCkT2Init(pin, pmy_mesh_);
  }
  two_stream_rt::rt_test_mu0 = pin->GetOrAddReal("problem","rt_test_mu0",-2.0);
  two_stream_rt::ck_test_kgrey = pin->GetOrAddReal("problem","ck_test_kgrey",0.0);
  if (two_stream_rt::ck_implicit && global_variable::my_rank == 0) {
    std::cout << "deep_hot_jupiter_rt: correlated-k source is IMPLICIT "
              << "(ck_implicit), tol " << two_stream_rt::ck_impl_tol << ", maxit "
              << two_stream_rt::ck_impl_maxit << std::endl;
  }
  // problem/rt_cell_report + problem/rt_report_every: the per-cell radiative-balance
  // report and, with it, the rt_desum line -- sum(src*dt*dx) against sum(de*dx), i.e.
  // what the sweep asked for against what the gas actually received once the
  // semi-implicit relaxation, the rt_de_max limiter and the floors have had their say.
  // red_giant.cpp and box_convection.cpp already read these; this path did not, so the
  // diagnostic was unreachable here.  Defaults are the namespace defaults, so this is
  // inert unless asked for.
  rt_cell_report = pin->GetOrAddBoolean("problem","rt_cell_report",false);
  rt_report_every = pin->GetOrAddInteger("problem","rt_report_every",100);
  // problem/rt_outer_verbose + problem/rt_apply_debug(_n): the two DIAGNOSTIC prints the
  // apply kernel already carries -- the per-200-cycle clip/rescue census and the per-cell
  // energy balance of the top rt_apply_debug_n cells of the (ck_dump_m, ck_dump_k,
  // ck_dump_j) column.  Neither was reachable from a dhj input; both default off and
  // change nothing that is applied to the gas.
  two_stream_rt::rt_outer_verbose =
      pin->GetOrAddBoolean("problem","rt_outer_verbose",false);
  two_stream_rt::rt_apply_debug = pin->GetOrAddInteger("problem","rt_apply_debug",0);
  two_stream_rt::rt_apply_debug_n =
      pin->GetOrAddInteger("problem","rt_apply_debug_n",8);
  rt_star_teff = pin->GetOrAddReal("problem","ck_star_teff",6000.0);
  // problem/albedo: observed Bond albedo in place of the Parmentier+2015 fit (see
  // two_stream_rt.hpp rt_albedo_set).  Read only when present, so an input without it
  // keeps an identical parameter dump and the fit, bit for bit.
  if (pin->DoesParameterExist("problem", "albedo")) {
    two_stream_rt::rt_albedo_set = pin->GetReal("problem", "albedo");
    if (!(two_stream_rt::rt_albedo_set >= 0.0 && two_stream_rt::rt_albedo_set < 1.0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "problem/albedo must be in [0, 1); got "
                << two_stream_rt::rt_albedo_set << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (global_variable::my_rank == 0) {
      std::cout << "deep_hot_jupiter_rt: Bond albedo fixed at problem/albedo = "
                << two_stream_rt::rt_albedo_set << " (replaces the Parmentier+2015 fit)"
                << std::endl;
    }
  }
  rt_dump_file = pin->GetOrAddString("problem","ck_dump_file","");
  two_stream_rt::rt_dump_kap = pin->GetOrAddBoolean("problem","ck_dump_kap",false);
  rt_dump_m = pin->GetOrAddInteger("problem","ck_dump_m",0);
  rt_dump_j = pin->GetOrAddInteger("problem","ck_dump_j",-1);
  rt_dump_k = pin->GetOrAddInteger("problem","ck_dump_k",-1);
  rt_ck_pcut = pin->GetOrAddReal("problem","ck_pcut_bar",10.0);
  // the grey-sweep fix, wired here so a dhj run can reproduce a pre-fix run bitwise.
  // Defaults to the FIX (false); see utils/two_stream_rt.hpp.
  // rt_layer_legacy: every layer used to span a whole cell in optical depth but carry a
  // source running between two cell CENTRES, offset from it by half a cell.
  two_stream_rt::rt_layer_legacy =
      pin->GetOrAddBoolean("problem","rt_layer_legacy",false);
  // problem/rt_ali_diag: the accelerated-Lambda diagonal in the per-cell semi-implicit
  // apply.  A PHYSICS CORRECTION, default TRUE: without it an optically thick cell is
  // handed 1/x of its own radiative source (x up to 1e5), so the thick layers of the
  // column never re-relax.  The LEGACY value is false, and reproduces a pre-fix run bit
  // for bit.  See the long note in utils/two_stream_rt.hpp.
  two_stream_rt::rt_ali_diag =
      pin->GetOrAddBoolean("problem","rt_ali_diag",true);
  rt_de_max = pin->GetOrAddReal("problem","rt_de_max",0.5);
  // problem/rt_floor_consistent: the radiative source may not cool a cell below the state
  // the EOS floors restore it to.  Default false = today's arithmetic, bit for bit; see
  // the long note in utils/two_stream_rt.hpp and tests_ck_topslab/README.md.
  two_stream_rt::rt_floor_consistent =
      pin->GetOrAddBoolean("problem","rt_floor_consistent",false);
  rt_semi_implicit = pin->GetOrAddBoolean("problem","rt_semi_implicit",true);
  // the sub-cycled local relaxation; see two_stream_rt.hpp, rt_relax_sub.  Default 1 =
  // the single-step form, bit for bit.  Read here, the one site both the from-scratch and
  // the restart path go through.
  rt_relax_sub = pin->GetOrAddInteger("problem","rt_relax_sub",1);
  rt_relax_xcrit = pin->GetOrAddReal("problem","rt_relax_xcrit",1.0);
  rt_relax_submax = pin->GetOrAddInteger("problem","rt_relax_submax",32);
  if (global_variable::my_rank == 0) {
    std::cout << "deep_hot_jupiter_rt: RT two-stream source is "
              << (rt_semi_implicit ? "SEMI-IMPLICIT" : "EXPLICIT") << std::endl;
  }
  // problem/rt_use_cons: read the RT solver's thermodynamic state from the CONSERVED
  // u0 instead of w0 (the previous stage's ConToPrim output, which the explicit source
  // terms and the implicit radial conduction have since moved).  See the long note in
  // utils/two_stream_rt.hpp.  DEFAULT TRUE since 2026-09-22 (it was false): the lagged
  // w0 state is a delayed thermostat whose forcing is linear in dt, so an input that
  // wants the old answer has to say rt_use_cons = false.  Read here, the one site both
  // the from-scratch and the restart path go through.
  rt_use_cons = pin->GetOrAddBoolean("problem","rt_use_cons",true);
  // per-cycle single-meshblock diagnostic dump. Read here, which is the one place both
  // the from-scratch and the restart path go through (CallProblemGenerator is called by
  // BOTH ProblemGenerator constructors), so the two can never drift apart -- see the
  // rot_potential restart bug for what happens when they do.
  diag_gid = pin->GetOrAddInteger("problem","diag_gid",-1);
  diag_cycle_min = pin->GetOrAddInteger("problem","diag_cycle_min",0);
  diag_cycle_max = pin->GetOrAddInteger("problem","diag_cycle_max",2147483647);
  {
    const std::string bn = pin->GetOrAddString("job","basename","athena");
    std::snprintf(diag_basename, sizeof(diag_basename), "%s", bn.c_str());
  }
  if (diag_gid >= 0 && global_variable::my_rank == 0) {
    std::cout << "deep_hot_jupiter_rt: cycle diagnostic dump for gid " << diag_gid
              << ", cycles [" << diag_cycle_min << "," << diag_cycle_max << "]"
              << std::endl;
  }
  rt_int_at_cut = pin->GetOrAddBoolean("problem","ck_int_at_cut",true);
  if (rt_ck && !rt_split) {
    // The correlated-k solver only exists inside the split path. Without this, rt_ck=true
    // with rt_split=false silently ran the GREY picket fence and looked like it worked.
    rt_split = true;
    if (global_variable::my_rank == 0) {
      std::cout << "deep_hot_jupiter_rt: problem/rt_ck implies problem/rt_split; "
                << "enabling the split RT path" << std::endl;
    }
  }
  // The EXACT implicit column (problem/rt_implicit_column = 3) is GREY-ONLY on this
  // branch: the correlated-k / picket-fence band sweep uses the centre-to-centre layer
  // branch, which the mode-3 column solve does not linearise.  This pgen also never
  // forwards the parameter to two_stream_rt, so refuse it here instead of ignoring it
  // silently.
  if (pin->GetOrAddInteger("problem","rt_implicit_column",0) != 0) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/rt_implicit_column is "
              << "not supported here.  The exact implicit column (mode 3) is GREY-ONLY "
              << "on this branch and the correlated-k / picket-fence band sweep carries "
              << "no Jacobian, so the column solve would run against a zero Jacobian."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // problem/ck_nquad: DEFAULT 2 (two-point Gauss) since defaults-0927 (2026-09-27; user
  // 09-25 "ck_nquad = 2 for dhj production": the exact diffusion limit, nquad 1 is ~10 %
  // low there; rms eint vs a tight reference 1.0e-2 vs 1.6e-2; +6 % ms/cycle at nx1 256).
  // 1 restores the diffusivity-factor form.  The value is recorded, so a restart keeps
  // what its file stored; an input that does not name it changes.
  ck_nq = pin->GetOrAddInteger("problem","ck_nquad",2);
  if (ck_nq != 1 && ck_nq != 2) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_nquad must be 1 "
              << "(diffusivity factor) or 2 (Gauss), got " << ck_nq << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // the ck table is also needed when the radiative diffusion takes its Rosseland
  // opacity from it (<hydro|mhd>/rad_kappa_src = table), even with the grey RT
  bool ck_need_table = rt_ck;
  {
    MeshBlockPack *pk = pmy_mesh_->pmb_pack;
    Conduction *pcd = (pk->pmhd != nullptr) ? pk->pmhd->pcond : pk->phydro->pcond;
    if (pcd != nullptr && pcd->rad_kappa_tab) ck_need_table = true;
    // The tau blend (<hydro|mhd>/rad_tau_hi > 0) splits each x1 face flux between the
    // conduction (weight rad_w) and the two-stream (1 - rad_w).  The conduction always
    // differences with the mesh's own face areas and volumes (spherical on a curvilinear
    // mesh); with ck_spherical = false the correlated-k deposit is plane-parallel (unit
    // areas, dx1).  The blended flux then mixes two divergences and does not telescope,
    // i.e. the blend is not conservative.  ck_spherical = true uses the same areas and
    // volumes as the conduction, and the blend is exact.
    if (rt_ck && !ck_spherical && (use_spherical_polar || use_cubed_sphere_) &&
        pcd != nullptr && pcd->rad_tau_mode) {
      std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/rt_ck = true with "
                << "problem/ck_spherical = false on a curvilinear mesh, together with "
                << "the radiative-conduction tau blend (rad_tau_hi > 0).  The ck "
                << "deposit is then plane-parallel while the conduction is spherical, "
                << "so the blended x1 face flux does not telescope (not conservative). "
                << "Set problem/ck_spherical = true (the curvilinear default), or turn "
                << "the blend off (rad_tau_hi = 0)." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  if (ck_need_table && ck_lk_ptr == nullptr) {
    read_ck_table(pin->GetOrAddString("problem","ck_table",
                                      "data/exo_fms_ck/ck/Premixed_1x_g8_11.txt"),
                  rt_ck_pcut);
    // problem/ck_interp_T: linear (default, bitwise) | pchip (cubic Hermite on exact
    // node slopes) for the band Planck fractions -- see correlated_k.hpp
    {
      const std::string it = pin->GetOrAddString("problem","ck_interp_T","linear");
      if (it == "pchip") {
        correlated_k::ck_pf_pchip = true;
      } else if (it != "linear") {
        std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_interp_T = '"
                  << it << "'; use linear or pchip." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    build_planck_fractions(rt_ck_pcut);
    read_ck_continuum(pin->GetOrAddString("problem","ck_data_dir",
                                          "data/exo_fms_ck"),
                      pin->GetOrAddString("problem","ck_swflux",
                                          "sw_band_flux_W121_11.txt"),
                      rt_star_teff,
                      // problem/ck_ce_table: composition table under ck_data_dir. Read
                      // only when set, so inputs without it keep an identical parameter
                      // dump. Pair the _hiT table with ck/Premixed_1x_g8_11_hiT.txt.
                      pin->DoesParameterExist("problem", "ck_ce_table") ?
                          pin->GetString("problem", "ck_ce_table") :
                          std::string("CE_tables/FastChem_ck_1x_int.txt"));
    ck_selftest();
    ck_rt_selftest();
    // The chain set is fixed by the table once correlated-k is on: every (band, g-point)
    // for each of the two Gauss points of the two-stream angular quadrature. Note this is
    // TWICE the 88 usually quoted for an 11 x 8 scheme -- 88 counts band x g, and this
    // kernel also carries the 2-point angular quadrature the picket fence uses.
    rt_nchain = CK_NB*CK_NG*ck_nq;
    if (global_variable::my_rank == 0) {
      std::cout << "  " << CK_NB << " bands x " << CK_NG << " g-points x " << ck_nq
                << " angular point(s) = " << rt_nchain << " column solves per cell";
      if (ck_nq == 1) {
        std::cout << " (diffusivity factor " << CK_DIFFUSIVITY << ")" << std::endl;
      } else {
        std::cout << " (2-point Gauss)" << std::endl;
      }
    }
  }
  user_bcs_func = HydrostaticEquilibrium;
  // NOTE: this function must NOT return early on a restart. Only the initial condition
  // (the "probini" kernel) and the magnetic field at the end are genuinely one-off;
  // everything in between builds BACKGROUND state that lives in memory only and is
  // therefore gone after a restart -- the gravitational potential phicc0/phi0 consumed by
  // etotgrav and by the well-balanced reconstruction, and the well-balanced reference
  // atmosphere u0wb/w0wb/w0facewb. Returning here left phi identically zero, which
  // silently changed the physics at the restart point for any run with etotgrav = true.
  // The two one-off blocks are guarded with `if (!restart)` individually instead.
  if (pmy_mesh_->one_d || pmy_mesh_->two_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "deep hot Jupiter problem generator only works in 3D" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (pin->GetInteger("mesh", "nx1") != pin->GetInteger("meshblock", "nx1")) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "deep hot Jupiter problem generator only allows one meshblock in r direction for the RT to work properly" << std::endl;
    exit(EXIT_FAILURE);
  }

  // capture variables for kernel
  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  // the cell's PANEL, which is what turns (x2,x3) into a direction on the cubed sphere
  auto &mbpanel_ = pmbp->pmb->mb_panel;
  auto &size = pmbp->pmb->mb_size;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  int &ng = indcs.ng;
  int n1m1 = indcs.nx1 + 2*ng - 1;
  int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
  int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
    
  // enrol the per-cycle diagnostic. The RT array is sized on the first RT call (which
  // happens before the first hook call); the conduction one has to be sized here,
  // because Conduction was built by AddPhysics before this pgen ran.
  if (diag_gid >= 0) {
    rt_diag = true;
    Conduction *pc = (pmbp->phydro != nullptr) ? pmbp->phydro->pcond
                                               : ((pmbp->pmhd != nullptr)
                                                  ? pmbp->pmhd->pcond : nullptr);
    if (pc != nullptr) {
      pc->EnableDiag(pmbp->nmb_thispack, n3m1+1, n2m1+1, n1m1+1);
    }
    user_cycle_func = DhjCycleDiag;
  }

  Real r0, r1;
  r0 = pmy_mesh_->mesh_size.x1min;
  r1 = pmy_mesh_->mesh_size.x1max;

  // Select either Hydro or MHD
  DvceArray5D<Real> u0_;
  DvceArray5D<Real> w0_;
    
  DvceArray4D<Real> phi0_x1f;
  DvceArray4D<Real> phi0_x2f;
  DvceArray4D<Real> phi0_x3f;
  DvceArray4D<Real> phicc0;
  DvceArray5D<Real> u0wb;
  DvceArray5D<Real> w0wb;
  DvceArray5D<Real> w0facewb_x1f;
  DvceArray5D<Real> w0facewb_x2f;
  DvceArray5D<Real> w0facewb_x3f;
    
  auto &x1v_ = pmbp->pcoord->x1v;
  auto &x1f_ = pmbp->pcoord->xx1f;
  auto &x2v_ = pmbp->pcoord->x2v;
  auto &x2f_ = pmbp->pcoord->xx2f;
  auto &x3v_ = pmbp->pcoord->x3v;
  auto &x3f_ = pmbp->pcoord->xx3f;
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &area2 = pmbp->pcoord->area.x2f;
  auto &area3 = pmbp->pcoord->area.x3f;
  auto &dxe1 = pmbp->pcoord->dxedge.x1e;
  auto &dxe2 = pmbp->pcoord->dxedge.x2e;
  auto &dxe3 = pmbp->pcoord->dxedge.x3e;

  Real gamma;
  if (pmbp->phydro != nullptr) {
    u0_ = pmbp->phydro->u0;
    w0_ = pmbp->phydro->w0;
    gamma = pmbp->phydro->peos->eos_data.gamma;
    use_etotgrav = pmbp->phydro->use_etotgrav;
    use_wellbalance_static = pmbp->phydro->use_wellbalance_static;
    use_wellbalance_dynamic = pmbp->phydro->use_wellbalance_dynamic;
    phi0_x1f = pmbp->phydro->phi0.x1f;
    phi0_x2f = pmbp->phydro->phi0.x2f;
    phi0_x3f = pmbp->phydro->phi0.x3f;
    phicc0 = pmbp->phydro->phicc0;
    u0wb = pmbp->phydro->u0wb;
    w0wb = pmbp->phydro->w0wb;
    w0facewb_x1f = pmbp->phydro->w0facewb.x1f;
    w0facewb_x2f = pmbp->phydro->w0facewb.x2f;
    w0facewb_x3f = pmbp->phydro->w0facewb.x3f;
  } else if (pmbp->pmhd != nullptr) {
    u0_ = pmbp->pmhd->u0;
    w0_ = pmbp->pmhd->w0;
    gamma = pmbp->pmhd->peos->eos_data.gamma;
    use_etotgrav = pmbp->pmhd->use_etotgrav;
    use_wellbalance_static = pmbp->pmhd->use_wellbalance_static;
    use_wellbalance_dynamic = pmbp->pmhd->use_wellbalance_dynamic;
    phi0_x1f = pmbp->pmhd->phi0.x1f;
    phi0_x2f = pmbp->pmhd->phi0.x2f;
    phi0_x3f = pmbp->pmhd->phi0.x3f;
    phicc0 = pmbp->pmhd->phicc0;
    u0wb = pmbp->pmhd->u0wb;
    w0wb = pmbp->pmhd->w0wb;
    w0facewb_x1f = pmbp->pmhd->w0facewb.x1f;
    w0facewb_x2f = pmbp->pmhd->w0facewb.x2f;
    w0facewb_x3f = pmbp->pmhd->w0facewb.x3f;
  }
  Real gm1 = gamma - 1.0;
  Real igm1 = 1.0/gm1;
  // by-value copy, capturable in the device lambdas below. Every pressure -> internal
  // energy conversion goes through it: for a general EOS e is not p/(gamma-1), and the
  // well-balanced background arrays built here feed the static scheme directly.
  auto eos = (pmbp->phydro != nullptr) ? pmbp->phydro->peos->eos_data
                                       : pmbp->pmhd->peos->eos_data;
    
//  Real Teq = 1469.0;
//  Real grav_acc = -942.0;
//  Real ap = 9.44e9;
//  Real omega = 2.06e-5;
//  Real Rgas = 4.593e7;
//  Real met = 0.0;
    
    // This problem generator (and utils/two_stream_rt.hpp) reads every dimensional
    // <problem> key as CGS without conversion -- Teq [K], grav [cm/s^2], ap [cm], omega
    // [1/s], ck_pcut_bar with 1 bar = 1e6 code pressure, bbot_gauss via sqrt(4 pi) --
    // so the code unit system must be cgs.  Refuse anything else rather than mis-scale.
    if (pmbp->punit != nullptr &&
        (pmbp->punit->length_cgs() != 1.0 || pmbp->punit->mass_cgs() != 1.0 ||
         pmbp->punit->time_cgs() != 1.0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "deep_hot_jupiter_rt reads its dimensional inputs as "
                << "cgs: <units> length_cgs, mass_cgs and time_cgs must all be 1"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    Real Teq = pin->GetReal("problem","Teq");
    Real grav_acc = -pin->GetReal("problem","grav");
    Real ap = pin->GetReal("problem","ap");
    // Gravity model: constant g (default, historical) or a point mass with the same g0 at
    // ap. Needs a RADIAL x1, which spherical polar and the CUBED SPHERE both have -- in
    // the Cartesian branch x1v is a height in a plane-parallel box and there is no r to
    // fall off with.
    const bool grav_pmass = pin->GetOrAddBoolean("problem","grav_point_mass",false);
    if (grav_pmass && !(use_spherical_polar || use_cubed_sphere_)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "problem/grav_point_mass requires a radial x1: "
                << "mesh/use_spherical_polar or mesh/use_cubed_sphere"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    // The stellar tide is NOT extended to the cubed sphere: it writes acceleration
    // components in the spherical (r, theta, phi) basis, and the gnomonic tangent pair is
    // neither that basis nor orthonormal, so the expressions would be silently wrong
    // rather than merely untested. It is deliberately off in production anyway.
    const bool tide = pin->GetOrAddBoolean("problem","stellar_tide",false);
    if (tide && !use_spherical_polar) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "problem/stellar_tide requires mesh/use_spherical_polar"
                << std::endl;
      exit(EXIT_FAILURE);
    }

    // Radial grid stretch, copied out of the Mesh so the device lambdas below capture
    // plain values. This function rebuilds radii with CellCenterX/LeftEdgeX, which are
    // UNSTRETCHED; each is passed through ApplyRStretch so the initial condition and the
    // gravitational potential land on the same grid the Coordinates arrays describe.
    const bool str_r_ = pmy_mesh_->use_grid_stretch_r;
    const bool str_rp_ = pmy_mesh_->use_grid_stretch_r_poly;
    const Real fstr_r_ = pmy_mesh_->fStretchR;
    const Real rmin_ = pmy_mesh_->mesh_size.x1min;
    const Real rmax_ = pmy_mesh_->mesh_size.x1max;
    Real cpoly_[NSTRETCH_R_POLY];
    for (int n=0; n<NSTRETCH_R_POLY; ++n) cpoly_[n] = pmy_mesh_->fStretchRPoly[n];
    Real omega = pin->GetReal("problem","omega");
    // problem/rot_potential: the centrifugal potential joins the gravitational one in
    // the column, the ghosts and the well-balanced arrays (see TotPotAt / ZEffFromPot)
    const bool rot_potential = pin->GetOrAddBoolean("problem","rot_potential",false);
    const Real om2_ = rot_potential ? SQR(omega) : 0.0;
    // roche-dhj-1031 problem/roche (default false: rgms_ = 0 and every TotPotAt below
    // is the old one, bitwise).  The host star's tide as a POTENTIAL: TidePot (companion
    // + indirect term) joins TotPotAt wherever this function builds a potential -- the
    // equipotential IC (z_eff of Phi_tot), the ghosts, phicc0 and the phi0 faces, hence
    // the well-balanced radial force (each column's own Phi) and the etotgrav energy.
    // The tangential tidal force is the explicit analytic -grad TidePot in SourceFunc
    // (or, with problem/roche_wb_tan, the hydro's own equilibrium residual: wbtan),
    // next to the centrifugal one; Coriolis unchanged.  Keys: roche_a [cm] (default
    // 0.02571 AU = WASP-121b) and G M_* from Kepler, Omega^2 a^3 = G (M_p + M_*), with
    // G M_p = grav ap^2 (point mass).  The centrifugal force stays in om2 (rot_potential,
    // required) -- not counted a second time anywhere.
    const bool roche_on = hot_jupiter_param.roche;
    Real rgms_ = 0.0, ra_ = 0.0;
    if (roche_on) {
      const bool wbx1 = (pmbp->phydro != nullptr) ? pmbp->phydro->use_wb_x1
                                                  : pmbp->pmhd->use_wb_x1;
      std::string bad = "";
      if (!rot_potential) bad = "problem/rot_potential = true";
      if (!grav_pmass) bad = "problem/grav_point_mass = true";
      if (!(use_spherical_polar || use_cubed_sphere_)) bad = "sp or the cubed sphere";
      if (!use_wellbalance_dynamic || !wbx1) bad = "wellbalance_dynamic + wb_x1";
      if (use_wellbalance_static) bad = "no wellbalance_static";
      if (!use_etotgrav) bad = "etotgrav = true";
      if (!user_srcs) bad = "problem/user_srcs = true";
      // the column-top closure tau = kappa p/|grad Phi| is wired into the band (ck)
      // solver only (rt_pre_geom fills rt_cf slot 4); the old grey paths keep g_pm
      if (!rt_ck) bad = "problem/rt_ck = true";
      if (tide) bad = "problem/stellar_tide = false (the tide is IN the Roche potential)";
      if (pin->DoesParameterExist("hydro","frame_omega") ||
          pin->DoesParameterExist("mhd","frame_omega")) {
        bad = "no <hydro|mhd>/frame_omega (the pgen applies the frame forces)";
      }
      if (!bad.empty()) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "problem/roche needs " << bad << std::endl;
        exit(EXIT_FAILURE);
      }
      ra_ = pin->GetOrAddReal("problem","roche_a",0.02571*1.495978707e13);
      rgms_ = SQR(omega)*ra_*ra_*ra_ - (-grav_acc)*ap*ap;
      // problem/roche_gms_scale (TEST ONLY, default 1): G M_* x scale, e.g. 1e-12 for
      // the no-star gate (must reproduce rot_potential to round-off); the lobe check is
      // skipped below 1e-3
      const Real gsc = pin->GetOrAddReal("problem","roche_gms_scale",1.0);
      rgms_ *= gsc;
      if (!(rgms_ > 0.0)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "problem/roche: Omega^2 a^3 <= G M_p" << std::endl;
        exit(EXIT_FAILURE);
      }
      hot_jupiter_param.roche_gms = rgms_;
      hot_jupiter_param.roche_a = ra_;
      const Real gm1r = (-grav_acc)*ap*ap;
      const Real rl1 = roche::L1Distance(gm1r, rgms_, ra_, om2_);
      // the potential at L1 (on the substellar axis: theta = pi/2, phi = 0) and the
      // highest potential on the top face (poles, theta = 0) -- the top equipotential is
      // closed (inside the planet's Roche lobe) iff the latter is below the former
      const Real phil1 = TotPotAt(grav_acc, ap, rl1, rl1 - ap, true, om2_, 0.5*M_PI, 0.0,
                                  rgms_, ra_);
      Real phtop = -1.0e300;
      for (int it=0; it<=180; ++it) {
        const Real th = M_PI*it/180.0;
        for (int ip=0; ip<=36; ++ip) {
          phtop = std::max(phtop, TotPotAt(grav_acc, ap, r1, r1 - ap, true, om2_, th,
                                           M_PI*ip/36.0, rgms_, ra_));
        }
      }
      if (global_variable::my_rank == 0) {
        std::cout << "deep_hot_jupiter_rt: problem/roche ON: a = " << ra_ << " cm, G M_* = "
                  << rgms_ << " (M_* = " << rgms_/6.6743e-8/1.98892e33 << " Msun), "
                  << "G M_p = " << gm1r << ", L1 at " << rl1 << " cm (x1max/r_L1 = "
                  << r1/rl1 << "), Phi(L1) = " << phil1 << ", max Phi on x1max = "
                  << phtop << std::endl;
      }
      // A WARNING, not a stop: the outer boundary holds the top hydrostatically, and the
      // region near L1 is outside the domain, so nothing can overflow in the model; but
      // the gas above that level is not bound to the planet (WASP-121b production grid:
      // x1max = 0.66 r_L1 on the axis, yet the lobe's POLAR radius is ~0.64 r_L1, so the
      // top cells near the poles sit just outside the Roche lobe)
      // topbc-1101: the AREA fraction of the top shell (r = x1max) on open
      // equipotentials (Phi > Phi(L1)), the guide for problem/top_bc: 0 -> hse, the
      // top is inside the lobe; > 0 -> outflow1 (no hydrostatic state above Phi(L1))
      {
        double aop = 0.0, atot = 0.0;
        const int nth = 360, nph = 180;   // phi in [0, pi]: the potential is even in phi
        for (int it=0; it<nth; ++it) {
          const Real th = M_PI*(it + 0.5)/nth;
          const double w = std::sin(th);
          for (int ip=0; ip<nph; ++ip) {
            const Real ph = M_PI*(ip + 0.5)/nph;
            atot += w;
            if (TotPotAt(grav_acc, ap, r1, r1 - ap, true, om2_, th, ph, rgms_, ra_)
                > phil1) aop += w;
          }
        }
        if (global_variable::my_rank == 0) {
          std::cout << "deep_hot_jupiter_rt: problem/roche: " << 100.0*aop/atot
                    << " % of the top shell (r = x1max) on OPEN equipotentials "
                    << "(Phi > Phi(L1)); problem/top_bc = "
                    << (dhj_top_bc == 0 ? "hse" : (dhj_top_bc == 1 ? "outflow1"
                                                                      : "outflow0"))
                    << (aop > 0.0 ? " (outflow1 recommended)" : " (hse recommended)")
                    << std::endl;
        }
      }
      if (gsc >= 1.0e-3 && !(phtop < phil1) && global_variable::my_rank == 0) {
        std::cout << "### WARNING: problem/roche: the equipotential through part of the "
                  << "domain top is OPEN (max Phi(x1max) >= Phi(L1)): those top cells lie "
                  << "outside the planet's Roche lobe" << std::endl;
      }
    }
    // roche-wbtan-1101 problem/roche_wb_tan (fresh default: proddef-1101 below; see
    // namespace wbtan): the tangential centrifugal (+ tidal, with problem/roche) force as
    // the hydro's own discrete residual of the equipotential equilibrium.  Works with
    // rot_potential alone (centrifugal only) or with roche.
    // proddef-1101 (USER): fresh-start default true with problem/roche on a single-level
    // mesh (roche itself already requires everything else roche_wb_tan needs); a restart
    // whose input lacks the key keeps the old false
    if (!pin->DoesParameterExist("problem","roche_wb_tan")) {
      const bool sl = !(pmy_mesh_->adaptive || pmy_mesh_->multilevel);
      const bool wd = roche_on && sl && recipe::NewDefaults(pin);
      if (roche_on) {
        recipe::Log(std::string("problem/roche_wb_tan = ") + (wd ? "true (default with "
                    "problem/roche)" : (!sl ? "false (needs a single-level mesh)"
                    : "false (old default: restart input or recipe none)")));
      }
      dhj_wbt.on = pin->GetOrAddBoolean("problem","roche_wb_tan",wd);
    } else {
      dhj_wbt.on = pin->GetBoolean("problem","roche_wb_tan");
    }
    if (dhj_wbt.on) {
      std::string bad = "";
      if (!rot_potential) bad = "problem/rot_potential = true";
      if (!(use_spherical_polar || use_cubed_sphere_)) bad = "sp or the cubed sphere";
      const bool wbx1t = (pmbp->phydro != nullptr) ? pmbp->phydro->use_wb_x1
                                                   : pmbp->pmhd->use_wb_x1;
      if (!use_wellbalance_dynamic || !wbx1t) {
        bad = "wellbalance_dynamic + wb_x1";
      }
      if (use_wellbalance_static) bad = "no wellbalance_static";
      if (!use_etotgrav) bad = "etotgrav = true (the work is booked by the energy flux)";
      if (!user_srcs) bad = "problem/user_srcs = true";
      if (tide) bad = "problem/stellar_tide = false";
      if (pmy_mesh_->adaptive || pmy_mesh_->multilevel) bad = "a single-level mesh";
      if (!bad.empty()) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "problem/roche_wb_tan needs " << bad << std::endl;
        exit(EXIT_FAILURE);
      }
      dhj_wbt.ready = false;
    }
    Real Rgas = pin->GetReal("problem","Rgas");
    Real met = pin->GetReal("problem","met");

    // <problem>/Rgas is the ideal-gas R/mu that fixes this atmosphere's mean molecular
    // weight. Under a general EOS composition lives in the EOS instead, so Rgas is unused
    // and the Kelvin conversion is eos.temp_cgs = (pres_cgs/dens_cgs) m_u/k_B, which
    // carries NO mu. The two therefore agree only at mu = 1, i.e. Rgas*temp_cgs = 1. That
    // is not required for a physical run -- but it is required for the general path to
    // reproduce the ideal one, which is how this pgen is verified, so say so loudly.
    if (eos.IsGeneral() && global_variable::my_rank == 0) {
      Real mu_implied = 1.0/(Rgas*eos.temp_cgs);
      if (fabs(mu_implied - 1.0) > 1.0e-4) {
        std::cout << std::endl << "### WARNING! in " << __FILE__ << " at line "
                  << __LINE__ << std::endl
                  << "<problem>/Rgas = " << Rgas << " implies mean molecular weight "
                  << mu_implied << ", but the general EOS supplies its own composition "
                  << "and is being asked for temperature directly." << std::endl
                  << "Rgas is now unused; this run will NOT reproduce the eos=ideal run. "
                  << "Set Rgas = " << 1.0/eos.temp_cgs << " for that comparison."
                  << std::endl << std::endl;
      }
    }

    // <problem>/met is [M/H] in dex and feeds the opacity fit in get_kapr(); the EOS
    // scales its metal electron donors by eos_metal_mh, which DEFAULTS to met and so
    // normally agrees automatically. It can only differ if it was set explicitly, and the
    // result would be an atmosphere opaque at one metallicity and conducting at another,
    // with the electron fraction -- hence the Ohmic resistivity -- wrong by
    // 10^(met - eos_metal_mh). Refuse rather than let that pass.
    if (eos.IsGeneral() && eos.MetalIonization()) {
      if (fabs(met - eos.MetalMetallicity()) > 1.0e-6) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "<problem>/met = " << met << " but the EOS metal donors were built "
                  << "with eos_metal_mh = " << eos.MetalMetallicity()
                  << "; both are [M/H] in dex and must agree. Remove the explicit "
                  << "eos_metal_mh and it will follow met." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    Real iap = 1.0/ap;
    
    Real grav = -grav_acc;
    Real Tirr = Teq*sqrt(2);
    Real Tint;
    get_Tint(Teq, Tint);
    Real mus = 1.0;
    
    const int N = 10000;
    // roche-dhj-1031: the column must reach the largest z_eff in the domain; the tide
    // RAISES Phi at the poles (and the terminator), past the 1.1 (r1 - r0) margin
    Real zcol_ = (r1-r0)*1.1;
    if (roche_on) {
      auto x1vh = Kokkos::create_mirror_view_and_copy(HostMemSpace(), x1v_);
      const int nimax = x1vh.extent_int(1) - 1;
      Real rt = x1vh(0, nimax);
      for (int m=0; m<pmbp->nmb_thispack; ++m) rt = std::max(rt, x1vh(m, nimax));
      // the largest z_eff of the top ghosts with and without the tide; the column grows
      // by the DIFFERENCE only, so a negligible tide (the no-star gate,
      // roche_gms_scale -> 0) leaves the rot_potential column (and IC) unchanged to
      // round-off.  (The rot_potential column itself stops at 1.1 (r1 - r0) and the top
      // ghosts read its linear extrapolation, as before.)
      Real ztop = 0.0, ztop0 = 0.0;
      for (int it=0; it<=90; ++it) {
        const Real th = 0.5*M_PI*it/90.0;
        for (int ip=0; ip<=36; ++ip) {
          ztop = std::max(ztop, ZEffFromPot(grav_acc, ap,
                              TotPotAt(grav_acc, ap, rt, rt - ap, grav_pmass, om2_, th,
                                       M_PI*ip/36.0, rgms_, ra_), grav_pmass));
          ztop0 = std::max(ztop0, ZEffFromPot(grav_acc, ap,
                               TotPotAt(grav_acc, ap, rt, rt - ap, grav_pmass, om2_, th),
                               grav_pmass));
        }
      }
      zcol_ += std::max(0.0, ztop - ztop0);
      if (global_variable::my_rank == 0) {
        std::cout << "deep_hot_jupiter_rt: problem/roche: column to z = " << zcol_
                  << " (1.1 (r1 - r0) + the tidal rise of the top ghosts' z_eff, "
                  << ztop - ztop0 << ")"
                  << std::endl;
      }
    }
    DualArray1D<Real> zarr("zarr", N);
    DualArray1D<Real> logparr("logparr", N);
    get_wb_eos_arr(eos, Rgas, grav_acc, ap, grav_pmass, N, zcol_,
                   zarr.h_view, logparr.h_view);
//    zarr.template modify<HostMemSpace>();
//    zarr.template sync<DevExeSpace>();
//    logparr.template modify<HostMemSpace>();
//    logparr.template sync<DevExeSpace>();
    zarr.modify_host();
    zarr.sync_device();
    logparr.modify_host();
    logparr.sync_device();
    
    DualArray1D<Real> Tarr_init("Tarrinit", N);
    DualArray1D<Real> lgparr_init("lgparrinit", N);
    // problem/ic_profile: a tabulated T(p) in place of the picket fence (+ 0.9 nabla_ad
    // adjustment).  Empty (default) keeps the picket-fence path bitwise.
    const std::string ic_prof = pin->GetOrAddString("problem", "ic_profile", "");
    if (ic_prof.empty()) {
      get_picket_fence_pT_arr(eos, Rgas, gamma, Tint, Tirr, met, grav, mus, N,
                              Tarr_init.h_view, lgparr_init.h_view);
    } else {
      read_ic_profile(ic_prof, N, Tarr_init.h_view, lgparr_init.h_view);
    }
    // problem/ic_adiabat_ext (default true; ic_profile only; w121hroche-1101): BELOW the
    // table's bottom row T(p) continues along the ADIABAT through that row,
    // dln T/dln p = GradAd(p, T) of the deck's own EOS (general table: H2 dissociation,
    // ionization), integrated (RK4) over ic_adiabat_ext_dex decades (default 2) on
    // Ne = 801 points uniform in log10 p.  Before, get_init_Tp(_host) continued the LAST
    // SEGMENT of the table linearly in log10 p, which the inner wall reaches whenever it
    // cuts equipotentials deeper than the table (rot_potential margins, problem/roche:
    // ~4400 bar substellar vs 2241 bar for WASP-121b).  Only lookups at log10 p > the
    // table bottom change; every lookup inside the table is bitwise (and a deck whose
    // column never leaves the table is bitwise).  false: the old linear extrapolation.
    // Beyond the extension's end: its last segment, linearly (WARNING printed).
    int ne_ext = 0;
    const int ne_alloc = 801;
    DualArray1D<Real> Tarr_ext("Tarrext", ne_alloc);
    DualArray1D<Real> lgparr_ext("lgparrext", ne_alloc);
    const bool ic_ad_ext = !ic_prof.empty() &&
                           pin->GetOrAddBoolean("problem", "ic_adiabat_ext", true);
    const Real ic_ad_dex = ic_prof.empty() ? 0.0 :
                           pin->GetOrAddReal("problem", "ic_adiabat_ext_dex", 2.0);
    if (ic_ad_ext) {
      ne_ext = ic_adiabat_extension(eos, gamma, Rgas, N, Tarr_init.h_view,
                                    lgparr_init.h_view, ic_ad_dex, ne_alloc,
                                    Tarr_ext.h_view, lgparr_ext.h_view);
      Tarr_ext.modify_host();
      Tarr_ext.sync_device();
      lgparr_ext.modify_host();
      lgparr_ext.sync_device();
    }
    
    Tarr_init.modify_host();
    Tarr_init.sync_device();
    lgparr_init.modify_host();
    lgparr_init.sync_device();
    
//    DvceArray1D<Real> zinitarr("zinitarr", N);
//    DvceArray1D<Real> logpinitarr("logpinitarr", N);
    DualArray1D<Real> zarr_init("zarrinit", N);
    DualArray1D<Real> logparr_init("logparrinit", N);
    // problem/ic_profile: extend the IC column below z = 0 to the deepest cell or inner
    // ghost (smallest radius, on the equator where rot_potential lowers z_eff most),
    // with a 10 % margin.  Picket-fence path (no ic_profile): zlow = 0, unchanged.
    Real zlow_ic = 0.0;
    if (!ic_prof.empty() && (use_spherical_polar || use_cubed_sphere_)) {
      auto x1v_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), x1v_);
      Real rg = x1v_h(0, 0);
      for (int m=0; m<pmbp->nmb_thispack; ++m) rg = std::min(rg, x1v_h(m, 0));
      Real zg = (om2_ > 0.0) ?
          ZEffFromPot(grav_acc, ap, TotPotAt(grav_acc, ap, rg, rg - ap, grav_pmass, om2_,
                                             0.5*M_PI), grav_pmass) : (rg - ap);
      // roche-dhj-1031: the lowest Phi on the deepest shell is on the substellar or the
      // antistellar point of the equator
      if (roche_on) {
        for (const Real ph : {0.0, M_PI}) {
          zg = std::min(zg, ZEffFromPot(grav_acc, ap,
                            TotPotAt(grav_acc, ap, rg, rg - ap, grav_pmass, om2_,
                                     0.5*M_PI, ph, rgms_, ra_), grav_pmass));
        }
      }
      if (zg < 0.0) zlow_ic = 1.1*zg;
      std::cout << "deep_hot_jupiter_rt: ic_profile column extended to z = " << zlow_ic
                << " cm below the 250 bar anchor (deepest ghost r = " << rg << ")"
                << std::endl;
    }
    get_init_eos_arr(eos, Rgas, grav_acc, ap, grav_pmass, Tarr_init.h_view,
                     lgparr_init.h_view, N, zcol_, zarr_init.h_view,
                     logparr_init.h_view, zlow_ic, ne_ext, Tarr_ext.h_view,
                     lgparr_ext.h_view);
    if (ne_ext > 1 && global_variable::my_rank == 0) {
      // where the column bottom sits: in the table or on the adiabat (diagnostic)
      const Real lgbot = logparr_init.h_view(0)/std::log(10.0);
      const Real lgtab = lgparr_init.h_view(N-1), lgend = lgparr_ext.h_view(ne_ext-1);
      Real tbot;
      get_init_Tp_host(N, Tarr_init.h_view, lgparr_init.h_view,
                       std::exp(logparr_init.h_view(0)), tbot, ne_ext, Tarr_ext.h_view,
                       lgparr_ext.h_view);
      std::cout << "deep_hot_jupiter_rt: IC column bottom p = " << std::pow(10.0, lgbot)/1.0e6
                << " bar, T = " << tbot << " K ("
                << ((lgbot > lgtab) ? "on the adiabat below the table bottom "
                                    : "inside the table; table bottom ")
                << std::pow(10.0, lgtab)/1.0e6 << " bar)" << std::endl;
      if (lgbot > lgend) {
        std::cout << "### WARNING: problem/ic_adiabat_ext: the IC column reaches "
                  << std::pow(10.0, lgbot)/1.0e6 << " bar, beyond the adiabat extension ("
                  << std::pow(10.0, lgend)/1.0e6 << " bar): raise problem/ic_adiabat_ext_dex"
                  << std::endl;
      }
    }
    
    zarr_init.modify_host();
    zarr_init.sync_device();
    logparr_init.modify_host();
    logparr_init.sync_device();

    // problem/remap_file (default empty = off; remap_0929): a FROM-SCRATCH start whose
    // active cells come from a file written by docs/handover/scripts/dhj_remap.py, i.e.
    // a restart radially remapped onto THIS grid (spin up on a coarse radial grid, then
    // continue on a fine one; the tool also refines the cubed sphere C(n) -> C(f n)).
    // The file is written for THIS grid and MeshBlock layout.  File: 8-byte magic,
    // int32 nmb_total, nx1, nx2, nx3; then per gid
    // int32 lx1, lx2, lx3, level, panel and double [5][nx3][nx2][nx1] = rho, v1, v2, v3
    // (= u0(IM1..3)/rho as stored, covariant on the cubed sphere), and the thermal
    // variable: p for magic "RADREMAP" (eint =
    // EintFromP(rho, p), any EOS) or the internal energy density for "RADREMPE".
    // E = eint + KE + rho*Phi_tot at the NEW cell centres (etotgrav; rot_potential
    // included), exactly as the analytic IC.  The ghosts keep the analytic column and are
    // refilled by Driver::Initialize; everything that is not u0 (ck caches, warm starts,
    // flux-history accumulators, floor bookkeeping) starts fresh as on any fresh start.
    // MHD: the file carries the hydro state only.  B is the analytic bbot field set below
    // (a hydro spin-up continued as MHD), or, with problem/remap_bfile, the remapped
    // field of an MHD restart INSTEAD of the seed.  Ignored on a restart.
    const std::string remap_file = pin->GetOrAddString("problem","remap_file","");
    const bool use_remap_ = (!restart && !remap_file.empty());
    const int rnx1 = indcs.nx1, rnx2 = indcs.nx2, rnx3 = indcs.nx3;
    DvceArray5D<Real> wrm_("remap_w", use_remap_ ? pmbp->nmb_thispack : 1, 5,
                           use_remap_ ? rnx3 : 1, use_remap_ ? rnx2 : 1,
                           use_remap_ ? rnx1 : 1);
    bool remap_eint_ = false;
    if (use_remap_) {
      std::ifstream rf(remap_file, std::ios::binary);
      char mg[8];
      std::int32_t hd[4];
      rf.read(mg, 8);
      rf.read(reinterpret_cast<char*>(hd), sizeof(hd));
      const std::string smg(mg, 8);
      remap_eint_ = (smg == "RADREMPE");
      if (!rf || (smg != "RADREMAP" && smg != "RADREMPE") ||
          hd[0] != pmy_mesh_->nmb_total || hd[1] != rnx1 || hd[2] != rnx2 ||
          hd[3] != rnx3) {
        std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/remap_file "
                  << remap_file << ": unreadable, bad magic or grid mismatch (file "
                  << hd[0] << " MeshBlocks of " << hd[1] << "x" << hd[2] << "x" << hd[3]
                  << ")" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      const std::size_t nrec = static_cast<std::size_t>(5)*rnx3*rnx2*rnx1;
      const std::streamoff recb = 5*sizeof(std::int32_t) + nrec*sizeof(double);
      auto hrm = Kokkos::create_mirror_view(wrm_);
      std::vector<double> buf(nrec);
      for (int m=0; m<pmbp->nmb_thispack; ++m) {
        const int gid = pmbp->gids + m;
        rf.seekg(8 + sizeof(hd) + gid*recb, std::ios::beg);
        std::int32_t ll[5];
        rf.read(reinterpret_cast<char*>(ll), sizeof(ll));
        rf.read(reinterpret_cast<char*>(buf.data()), nrec*sizeof(double));
        const LogicalLocation &lc = pmy_mesh_->lloc_eachmb[gid];
        if (!rf || ll[0] != lc.lx1 || ll[1] != lc.lx2 || ll[2] != lc.lx3 ||
            ll[3] != lc.level || ll[4] != lc.panel) {
          std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/remap_file "
                    << "record " << gid << " does not match this MeshBlock's "
                    << "LogicalLocation" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        std::size_t q = 0;
        for (int n=0; n<5; ++n) {
          for (int k=0; k<rnx3; ++k) {
            for (int j=0; j<rnx2; ++j) {
              for (int i=0; i<rnx1; ++i) hrm(m,n,k,j,i) = buf[q++];
            }
          }
        }
      }
      Kokkos::deep_copy(wrm_, hrm);
      if (global_variable::my_rank == 0) {
        std::cout << "deep_hot_jupiter_rt: problem/remap_file " << remap_file << " read ("
                  << pmy_mesh_->nmb_total << " MeshBlocks, nx1 = " << rnx1 << ", "
                  << (remap_eint_ ? "eint" : "p") << ")" << std::endl;
      }
    }
    const bool rm_eint_ = remap_eint_;
    auto &rm_ccell_ = pmbp->pcoord->cos_cell;   // read only on the cubed sphere

    // one-off: the initial condition. Skipped on a restart, where u0/w0 come from file.
    if (!restart) {
    par_for("probini", DevExeSpace(), 0, (pmbp->nmb_thispack-1), 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
        
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      int nx1 = indcs.nx1;
      Real &x2min = size.d_view(m).x2min;
      Real &x2max = size.d_view(m).x2max;
      int nx2 = indcs.nx2;
      Real &x3min = size.d_view(m).x3min;
      Real &x3max = size.d_view(m).x3max;
      int nx3 = indcs.nx3;
        
      Real x1v, x2v, x3v;
      if (use_spherical_polar) {
        x1v = x1v_(m,i);
        x2v = x2v_(m,j);
        x3v = x3v_(m,k);
      } else {
        x1v = CellCenterX(i-is, nx1, x1min, x1max);
        ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
        if (use_cubed_sphere_) x1v = x1v_(m,i);   // stored centroid/face, as sp
        x2v = CellCenterX(j-js, nx2, x2min, x2max);
        x3v = CellCenterX(k-ks, nx3, x3min, x3max);
      }
      Real r = x1v;
      if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
        
      Real lam, theta, phi;
      if (use_spherical_polar) {
        theta = x2v;
        lam = -x2v+M_PI/2.0;
        phi = x3v-M_PI;
      } else if (use_cubed_sphere_) {
        CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
      } else {
        theta = -(x2v*iap-M_PI/2.0);
        lam = x2v*iap;
        phi = x1v*iap;
      }
        
      Real pwb, denwb;
      // rot_potential: sample the column at the equivalent height of this cell's total
      // potential (barotropic in Phi_tot); otherwise at its geometric height
      const Real zs = (om2_ > 0.0) ?
          ZEffFromPot(grav_acc, ap,
                      TotPotAt(grav_acc, ap, r, x1v, grav_pmass,
                               om2_, theta, phi, rgms_, ra_),
                      grav_pmass) : x1v;
      get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,zs,denwb,pwb);
        
      Real p, den;
      get_init_eos(eos, Rgas, grav_acc, Tarr_init.d_view,lgparr_init.d_view,N,zarr_init.d_view,logparr_init.d_view,zs,den,p,
                   ne_ext, Tarr_ext.d_view, lgparr_ext.d_view);
//      get_init_eos_arr(lam, phi, N, zinitarr, logpinitarr);
//      get_init_eos(zinitarr,logpinitarr,x3v,lam,phi,den,p);
//      p = pwb;
//      den = denwb;

      if (!use_remap_) {
      u0_(m,IDN,k,j,i) = den;
      u0_(m,IM1,k,j,i) = 0.0;
      u0_(m,IM2,k,j,i) = 0.0;
      u0_(m,IM3,k,j,i) = 0.0;
      u0_(m,IEN,k,j,i) = EintFromP(eos, igm1, den, p);
        
      w0_(m,IDN,k,j,i) = den;
      w0_(m,IVX,k,j,i) = 0.0;
      w0_(m,IVY,k,j,i) = 0.0;   // was IVX three times; harmless, ConsToPrim rebuilt w0
      w0_(m,IVZ,k,j,i) = 0.0;
      w0_(m,IEN,k,j,i) = EintFromP(eos, igm1, den, p);
      } else {
        // problem/remap_file: active cells from the file, ghosts analytic (see above)
        Real v1 = 0.0, v2 = 0.0, v3 = 0.0, eint0;
        if (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke) {
          den = wrm_(m,0,k-ks,j-js,i-is);
          v1  = wrm_(m,1,k-ks,j-js,i-is);
          v2  = wrm_(m,2,k-ks,j-js,i-is);
          v3  = wrm_(m,3,k-ks,j-js,i-is);
          eint0 = rm_eint_ ? wrm_(m,4,k-ks,j-js,i-is)
                           : EintFromP(eos, igm1, den, wrm_(m,4,k-ks,j-js,i-is));
        } else {
          eint0 = EintFromP(eos, igm1, den, p);
        }
        // v1..3 = u0(IM1..3)/rho as stored (on the cubed sphere the COVARIANT
        // components); KE = 0.5 m_i v^i with the gnomonic metric, as ConToPrim forms it
        Real vu2 = v2, vu3 = v3;
        if (use_cubed_sphere_) {
          const Real c = rm_ccell_(m,k,j), det = 1.0 - c*c;
          vu2 = (v2 - c*v3)/det;
          vu3 = (v3 - c*v2)/det;
        }
        u0_(m,IDN,k,j,i) = den;
        u0_(m,IM1,k,j,i) = den*v1;
        u0_(m,IM2,k,j,i) = den*v2;
        u0_(m,IM3,k,j,i) = den*v3;
        u0_(m,IEN,k,j,i) = eint0 + 0.5*den*(v1*v1 + v2*vu2 + v3*vu3);
        w0_(m,IDN,k,j,i) = den;
        w0_(m,IVX,k,j,i) = v1;
        w0_(m,IVY,k,j,i) = vu2;
        w0_(m,IVZ,k,j,i) = vu3;
        w0_(m,IEN,k,j,i) = eint0;
      }
        
        Real phicc = TotPotAt(grav_acc, ap, r, x1v, grav_pmass,
                              om2_, theta, phi, rgms_, ra_);
        if (use_etotgrav) {
            u0_(m,IEN,k,j,i) += den*phicc;
        }
    });
    }  // end of !restart guard on the initial condition

    // roche-wbtan-1101: the equipotential equilibrium W_eq (rho, eint; v = 0) of EVERY
    // cell incl. the ghosts, by the IC formula above (fresh start AND restart: the probe
    // in SourceFunc runs once per run segment)
    if (dhj_wbt.on) {
      const int nmbw = pmbp->nmb_thispack;
      wbtan::Alloc(dhj_wbt, pmbp);
      auto weq_ = *dhj_wbt.weq;
      par_for("wbtan_weq", DevExeSpace(), 0, nmbw-1, 0, n3m1, 0, n2m1, 0, n1m1,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        Real x1v, x2v, x3v;
        if (use_spherical_polar) {
          x1v = x1v_(m,i);
          x2v = x2v_(m,j);
          x3v = x3v_(m,k);
        } else {
          x1v = x1v_(m,i);   // cubed sphere: stored centroid, as the IC
          x2v = CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min, size.d_view(m).x2max);
          x3v = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min, size.d_view(m).x3max);
        }
        const Real r = x1v;
        x1v -= ap;
        Real lam = 0.0, theta, phi;
        if (use_spherical_polar) {
          theta = x2v;
          phi = x3v-M_PI;
        } else {
          CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
        }
        const Real zs = ZEffFromPot(grav_acc, ap,
                                    TotPotAt(grav_acc, ap, r, x1v, grav_pmass,
                                             om2_, theta, phi, rgms_, ra_),
                                    grav_pmass);
        Real p, den;
        get_init_eos(eos, Rgas, grav_acc, Tarr_init.d_view, lgparr_init.d_view, N,
                     zarr_init.d_view, logparr_init.d_view, zs, den, p,
                     ne_ext, Tarr_ext.d_view, lgparr_ext.d_view);
        weq_(m,0,k,j,i) = den;
        weq_(m,1,k,j,i) = EintFromP(eos, igm1, den, p);
      });
      dhj_wbt.ready = false;
      // wbtan-rstg-1101: a restart takes g_t from the file when it carries it
      if (restart) wbtan::LoadFromRestart(this, dhj_wbt, pmbp, "deep_hot_jupiter_rt");
      if (global_variable::my_rank == 0) {
        std::cout << "deep_hot_jupiter_rt: problem/roche_wb_tan ON: tangential "
                  << (roche_on ? "centrifugal + tidal" : "centrifugal")
                  << " force = the hydro's own residual of the equipotential state"
                  << std::endl;
      }
    }

    // problem/seed + problem/seed_amp: a deterministic density perturbation on the
    // from-scratch initial condition.  It exists to build ENSEMBLES: the cubed-sphere
    // flow amplifies round-off to O(1) within 0.1 rotation, so distinct realizations
    // need only a 1e-10 kick, and events like the vertex dt collapse can only be studied
    // as a rate over many of them.  The random number is a hash of (seed, global cell
    // index), so it is identical for any MeshBlock decomposition or MPI rank count.
    const int ic_seed = pin->GetOrAddInteger("problem","seed",0);
    const Real seed_amp = pin->GetOrAddReal("problem","seed_amp",0.0);
    // problem/seed_restart = true: apply the same kick to the state read from a restart
    // (a chaotic-noise twin of a restarted run).  Default false: restarts untouched.
    const bool seed_rst = pin->GetOrAddBoolean("problem","seed_restart",false);
    if ((!restart || seed_rst) && seed_amp > 0.0) {
      // global index offset of each MeshBlock, from its LogicalLocation on the host
      // (panel, global i/j/k of the block's first active cell). The LEVEL is
      // deliberately NOT hashed: the root level shifts when the MeshBlock size changes,
      // which would make the perturbation decomposition-dependent.
      DualArray2D<int> gidx("ic_seed_gidx", pmbp->nmb_thispack, 4);
      for (int m=0; m<pmbp->nmb_thispack; ++m) {
        LogicalLocation &loc = pmy_mesh_->lloc_eachmb[pmbp->gids + m];
        gidx.h_view(m,0) = loc.panel;
        gidx.h_view(m,1) = loc.lx1*indcs.nx1;
        gidx.h_view(m,2) = loc.lx2*indcs.nx2;
        gidx.h_view(m,3) = loc.lx3*indcs.nx3;
      }
      gidx.modify_host();
      gidx.sync_device();

      if (global_variable::my_rank == 0) {
        std::cout << "deep_hot_jupiter_rt: IC density perturbation seed = " << ic_seed
                  << ", amplitude = " << seed_amp << std::endl;
      }
      auto u0seed = u0_;
      auto w0seed = w0_;
      auto gidx_ = gidx;
      const std::uint64_t sd = static_cast<std::uint64_t>(
                                 static_cast<std::int64_t>(ic_seed));
      const Real amp = seed_amp;
      par_for("probseed", DevExeSpace(), 0, (pmbp->nmb_thispack-1), ks, ke, js, je,
      is, ie,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        std::uint64_t h = SplitMix64Mix(sd);
        h = SplitMix64Mix(h ^ static_cast<std::uint64_t>(gidx_.d_view(m,0)));
        h = SplitMix64Mix(h ^ static_cast<std::uint64_t>(gidx_.d_view(m,1) + (i-is)));
        h = SplitMix64Mix(h ^ static_cast<std::uint64_t>(gidx_.d_view(m,2) + (j-js)));
        h = SplitMix64Mix(h ^ static_cast<std::uint64_t>(gidx_.d_view(m,3) + (k-ks)));
        // top 53 bits -> double in [0,1)
        Real rr = static_cast<Real>(h >> 11)*(1.0/9007199254740992.0);
        Real fac = 1.0 + amp*(2.0*rr - 1.0);
        u0seed(m,IDN,k,j,i) *= fac;
        u0seed(m,IM1,k,j,i) *= fac;
        u0seed(m,IM2,k,j,i) *= fac;
        u0seed(m,IM3,k,j,i) *= fac;
        w0seed(m,IDN,k,j,i) *= fac;
      });
    }

    par_for("probwb", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
        
        Real &x1min = size.d_view(m).x1min;
        Real &x1max = size.d_view(m).x1max;
        int nx1 = indcs.nx1;
        Real &x2min = size.d_view(m).x2min;
        Real &x2max = size.d_view(m).x2max;
        int nx2 = indcs.nx2;
        Real &x3min = size.d_view(m).x3min;
        Real &x3max = size.d_view(m).x3max;
        int nx3 = indcs.nx3;
        
        Real x1v, x2v, x3v;
        if (use_spherical_polar) {
          x1v = x1v_(m,i);
          x2v = x2v_(m,j);
          x3v = x3v_(m,k);
        } else {
          x1v = CellCenterX(i-is, nx1, x1min, x1max);
          ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
          if (use_cubed_sphere_) x1v = x1v_(m,i);   // stored centroid/face, as sp
          x2v = CellCenterX(j-js, nx2, x2min, x2max);
          x3v = CellCenterX(k-ks, nx3, x3min, x3max);
        }
        Real r = x1v;
        if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
          
        Real lam, theta, phi;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -x2v+M_PI/2.0;
          phi = x3v-M_PI;
        } else if (use_cubed_sphere_) {
          CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
        } else {
          theta = -(x3v*iap-M_PI/2.0);
          lam = x3v*iap;
          phi = x2v*iap;
        }
          
        Real pwb, denwb;
        const Real zs = (om2_ > 0.0) ?
            ZEffFromPot(grav_acc, ap,
                        TotPotAt(grav_acc, ap, r, x1v, grav_pmass,
                                 om2_, theta, phi, rgms_, ra_),
                        grav_pmass) : x1v;
        get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,zs,denwb,pwb);
        
        Real p, den;
        get_init_eos(eos, Rgas, grav_acc, Tarr_init.d_view,lgparr_init.d_view,N,zarr_init.d_view,logparr_init.d_view,zs,den,p,
                   ne_ext, Tarr_ext.d_view, lgparr_ext.d_view);
//        get_init_eos_arr(lam, phi, N, zinitarr, logpinitarr);
//        get_init_eos(zinitarr,logpinitarr,x1v,lam,phi,den,p);
//        p = pwb;
//        den = denwb;
        
        Real phicc = TotPotAt(grav_acc, ap, r, x1v, grav_pmass,
                              om2_, theta, phi, rgms_, ra_);
        
      if (use_etotgrav || use_wellbalance_dynamic) {
          if (use_spherical_polar) {
            x1v = x1f_(m,i);
            x2v = x2v_(m,j);
            x3v = x3v_(m,k);
          } else {
            x1v = LeftEdgeX(i-is, nx1, x1min, x1max);
            ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
            if (use_cubed_sphere_) x1v = x1f_(m,i);   // stored centroid/face, as sp
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          }
          r = x1v;
          if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
          if (use_spherical_polar) {
            theta = x2v;
            lam = -x2v+M_PI/2.0;
            phi = x3v-M_PI;
          } else if (use_cubed_sphere_) {
            CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
          } else {
            theta = -(x3v*iap-M_PI/2.0);
            lam = x3v*iap;
            phi = x2v*iap;
          }
          phicc = TotPotAt(grav_acc, ap, r, x1v, grav_pmass,
                           om2_, theta, phi, rgms_, ra_);
          phi0_x1f(m,k,j,i) = phicc;
          if (i == ie) {
              if (use_spherical_polar) {
                x1v = x1f_(m,i+1);
              } else {
                x1v = LeftEdgeX(i+1-is, nx1, x1min, x1max);
                ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
                if (use_cubed_sphere_) x1v = x1f_(m,i+1);   // stored centroid/face, as sp
              }
              r = x1v;
              if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
              if (use_spherical_polar) {
                theta = x2v;
                lam = -x2v+M_PI/2.0;
                phi = x3v-M_PI;
              } else if (use_cubed_sphere_) {
                CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
              } else {
                theta = -(x3v*iap-M_PI/2.0);
                lam = x3v*iap;
                phi = x2v*iap;
              }
              phicc = TotPotAt(grav_acc, ap, r, x1v, grav_pmass,
                               om2_, theta, phi, rgms_, ra_);
              phi0_x1f(m,k,j,i+1) = phicc;
          }
          
          if (use_spherical_polar) {
            x1v = x1v_(m,i);
            x2v = x2f_(m,j);
            x3v = x3v_(m,k);
          } else {
            x1v = CellCenterX(i-is, nx1, x1min, x1max);
            ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
            if (use_cubed_sphere_) x1v = x1v_(m,i);   // stored centroid/face, as sp
            x2v = LeftEdgeX(j-js, nx2, x2min, x2max);
            x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          }
          r = x1v;
          if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
          if (use_spherical_polar) {
            theta = x2v;
            lam = -x2v+M_PI/2.0;
            phi = x3v-M_PI;
          } else if (use_cubed_sphere_) {
            CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
          } else {
            theta = -(x3v*iap-M_PI/2.0);
            lam = x3v*iap;
            phi = x2v*iap;
          }
          phicc = TotPotAt(grav_acc, ap, r, x1v, grav_pmass,
                           om2_, theta, phi, rgms_, ra_);
          phi0_x2f(m,k,j,i) = phicc;
          if (j == je) {
              if (use_spherical_polar) {
                x2v = x2f_(m,j+1);
              } else {
                x2v = LeftEdgeX(j+1-js, nx2, x2min, x2max);
              }
              if (use_spherical_polar) {
                theta = x2v;
                lam = -x2v+M_PI/2.0;
                phi = x3v-M_PI;
              } else if (use_cubed_sphere_) {
                CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
              } else {
                theta = -(x3v*iap-M_PI/2.0);
                lam = x3v*iap;
                phi = x2v*iap;
              }
              phicc = TotPotAt(grav_acc, ap, r, x1v, grav_pmass,
                               om2_, theta, phi, rgms_, ra_);
              phi0_x2f(m,k,j+1,i) = phicc;
          }
          
          if (use_spherical_polar) {
            x1v = x1v_(m,i);
            x2v = x2v_(m,j);
            x3v = x3f_(m,k);
          } else {
            x1v = CellCenterX(i-is, nx1, x1min, x1max);
            ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
            if (use_cubed_sphere_) x1v = x1v_(m,i);   // stored centroid/face, as sp
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            x3v = LeftEdgeX(k-ks, nx3, x3min, x3max);
          }
          r = x1v;
          if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
          if (use_spherical_polar) {
            theta = x2v;
            lam = -x2v+M_PI/2.0;
            phi = x3v-M_PI;
          } else if (use_cubed_sphere_) {
            CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
          } else {
            theta = -(x3v*iap-M_PI/2.0);
            lam = x3v*iap;
            phi = x2v*iap;
          }
          phicc = TotPotAt(grav_acc, ap, r, x1v, grav_pmass,
                           om2_, theta, phi, rgms_, ra_);
          phi0_x3f(m,k,j,i) = phicc;
          if (k == ke) {
              if (use_spherical_polar) {
                x3v = x3f_(m,k+1);
              } else {
                x3v = LeftEdgeX(k+1-ks, nx3, x3min, x3max);
              }
              if (use_spherical_polar) {
                theta = x2v;
                lam = -x2v+M_PI/2.0;
                phi = x3v-M_PI;
              } else if (use_cubed_sphere_) {
                CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
              } else {
                theta = -(x3v*iap-M_PI/2.0);
                lam = x3v*iap;
                phi = x2v*iap;
              }
              phicc = TotPotAt(grav_acc, ap, r, x1v, grav_pmass,
                               om2_, theta, phi, rgms_, ra_);
              phi0_x3f(m,k+1,j,i) = phicc;
          }
      }
        if (use_wellbalance_static) {
            if (use_spherical_polar) {
              x1v = x1f_(m,i);
              x2v = x2v_(m,j);
              x3v = x3v_(m,k);
            } else {
              x1v = LeftEdgeX(i-is, nx1, x1min, x1max);
              ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
              if (use_cubed_sphere_) x1v = x1f_(m,i);   // stored centroid/face, as sp
              x2v = CellCenterX(j-js, nx2, x2min, x2max);
              x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            }
            if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
            Real denwb;
            get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
            w0facewb_x1f(m,IDN,k,j,i) = denwb;
            w0facewb_x1f(m,IM1,k,j,i) = 0.0;
            w0facewb_x1f(m,IM2,k,j,i) = 0.0;
            w0facewb_x1f(m,IM3,k,j,i) = 0.0;
            w0facewb_x1f(m,IEN,k,j,i) = EintFromP(eos, igm1, denwb, pwb);
            if (i == ie) {
                if (use_spherical_polar) {
                  x1v = x1f_(m,i+1);
                } else {
                  x1v = LeftEdgeX(i+1-is, nx1, x1min, x1max);
                  ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
                  if (use_cubed_sphere_) x1v = x1f_(m,i+1);   // as sp
                }
                if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
                get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
                w0facewb_x1f(m,IDN,k,j,i+1) = denwb;
                w0facewb_x1f(m,IM1,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IM2,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IM3,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IEN,k,j,i+1) = EintFromP(eos, igm1, denwb, pwb);
            }
            
            if (use_spherical_polar) {
              x1v = x1v_(m,i);
              x2v = x2f_(m,j);
              x3v = x3v_(m,k);
            } else {
              x1v = CellCenterX(i-is, nx1, x1min, x1max);
              ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
              if (use_cubed_sphere_) x1v = x1v_(m,i);   // stored centroid/face, as sp
              x2v = LeftEdgeX(j-js, nx2, x2min, x2max);
              x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            }
            if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
            get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
            w0facewb_x2f(m,IDN,k,j,i) = denwb;
            w0facewb_x2f(m,IM1,k,j,i) = 0.0;
            w0facewb_x2f(m,IM2,k,j,i) = 0.0;
            w0facewb_x2f(m,IM3,k,j,i) = 0.0;
            w0facewb_x2f(m,IEN,k,j,i) = EintFromP(eos, igm1, denwb, pwb);
            if (j == je) {
                w0facewb_x2f(m,IDN,k,j+1,i) = denwb;
                w0facewb_x2f(m,IM1,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IM2,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IM3,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IEN,k,j+1,i) = EintFromP(eos, igm1, denwb, pwb);
            }
            
            if (use_spherical_polar) {
              x1v = x1v_(m,i);
              x2v = x2v_(m,j);
              x3v = x3f_(m,k);
            } else {
              x1v = CellCenterX(i-is, nx1, x1min, x1max);
              ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
              if (use_cubed_sphere_) x1v = x1v_(m,i);   // stored centroid/face, as sp
              x2v = CellCenterX(j-js, nx2, x2min, x2max);
              x3v = LeftEdgeX(k-ks, nx3, x3min, x3max);
            }
            if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
            get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
            w0facewb_x3f(m,IDN,k,j,i) = denwb;
            w0facewb_x3f(m,IM1,k,j,i) = 0.0;
            w0facewb_x3f(m,IM2,k,j,i) = 0.0;
            w0facewb_x3f(m,IM3,k,j,i) = 0.0;
            w0facewb_x3f(m,IEN,k,j,i) = EintFromP(eos, igm1, denwb, pwb);
            if (k == ke) {
                w0facewb_x3f(m,IDN,k+1,j,i) = denwb;
                w0facewb_x3f(m,IM1,k+1,j,i) = 0.0;
                w0facewb_x3f(m,IM2,k+1,j,i) = 0.0;
                w0facewb_x3f(m,IM3,k+1,j,i) = 0.0;
                w0facewb_x3f(m,IEN,k+1,j,i) = EintFromP(eos, igm1, denwb, pwb);
            }
        }
    });
    if (use_etotgrav || use_wellbalance_dynamic) {
        int &ng = indcs.ng;
        int n1m1 = indcs.nx1 + 2*ng - 1;
        int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
        int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
        // roche-dhj-1031 problem/roche_rst_add (default false; read only on a restart
        // with problem/roche): the restart file was written WITHOUT the tide, so its
        // E = eint + KE + rho Phi_old; add rho TidePot so eint is unchanged.  ONE-SHOT:
        // the key is reset to false in the parameters, so restart files written by this
        // run (which carry Phi_Roche) do not add it a second time.
        bool rst_add_ = false;
        if (restart && roche_on && pin->DoesParameterExist("problem","roche_rst_add")) {
          rst_add_ = pin->GetBoolean("problem","roche_rst_add");
          pin->SetBoolean("problem","roche_rst_add",false);
          if (global_variable::my_rank == 0 && rst_add_) {
            std::cout << "deep_hot_jupiter_rt: problem/roche_rst_add: E += rho TidePot "
                      << "on the restarted state (one-shot)" << std::endl;
          }
        }
        par_for("wbgrav", DevExeSpace(), 0, (pmbp->nmb_thispack-1), 0, n3m1, 0, n2m1, 0, n1m1,
        KOKKOS_LAMBDA(int m, int k, int j, int i) {
            
            Real &x1min = size.d_view(m).x1min;
            Real &x1max = size.d_view(m).x1max;
            int nx1 = indcs.nx1;
            Real &x2min = size.d_view(m).x2min;
            Real &x2max = size.d_view(m).x2max;
            int nx2 = indcs.nx2;
            Real &x3min = size.d_view(m).x3min;
            Real &x3max = size.d_view(m).x3max;
            int nx3 = indcs.nx3;
            
            Real x1v, x2v, x3v;
            if (use_spherical_polar) {
              x1v = x1v_(m,i);
              x2v = x2v_(m,j);
              x3v = x3v_(m,k);
            } else {
              x1v = CellCenterX(i-is, nx1, x1min, x1max);
              ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
              if (use_cubed_sphere_) x1v = x1v_(m,i);   // stored centroid/face, as sp
              x2v = CellCenterX(j-js, nx2, x2min, x2max);
              x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            }
            Real r = x1v;
            if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
              
            Real lam, theta, phi;
            if (use_spherical_polar) {
              theta = x2v;
              lam = -x2v+M_PI/2.0;
              phi = x3v-M_PI;
            } else if (use_cubed_sphere_) {
              CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
            } else {
              theta = -(x3v*iap-M_PI/2.0);
              lam = x3v*iap;
              phi = x2v*iap;
            }
            
            Real phicc = TotPotAt(grav_acc, ap, r, x1v, grav_pmass,
                                  om2_, theta, phi, rgms_, ra_);
            phicc0(m,k,j,i) = phicc;
            if (rst_add_ && use_etotgrav) {
              u0_(m,IEN,k,j,i) += u0_(m,IDN,k,j,i)*
                  (phicc - TotPotAt(grav_acc, ap, r, x1v, grav_pmass, om2_, theta));
            }
            // The FACE potentials on every face, ghost faces included.  The initial
            // condition above fills only the active faces, but the well-balanced
            // background of the outermost ghost cells (which the x1 sweep reconstructs)
            // reads the ghost faces too; left at zero they made that background diverge.
            Real x1fl, x1fr;
            if (use_spherical_polar || use_cubed_sphere_) {
              x1fl = x1f_(m,i);
              x1fr = x1f_(m,i+1);
            } else {
              x1fl = LeftEdgeX(i-is, nx1, x1min, x1max);
              ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1fl);
              x1fr = LeftEdgeX(i+1-is, nx1, x1min, x1max);
              ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1fr);
            }
            Real zl = x1fl, zr = x1fr;
            if (use_spherical_polar || use_cubed_sphere_) { zl -= ap; zr -= ap; }
            phi0_x1f(m,k,j,i) = TotPotAt(grav_acc, ap, x1fl, zl, grav_pmass, om2_, theta,
                                          phi, rgms_, ra_);
            if (i == n1m1) {
              phi0_x1f(m,k,j,i+1) = TotPotAt(grav_acc, ap, x1fr, zr, grav_pmass, om2_,
                                              theta, phi, rgms_, ra_);
            }
        });
//        par_for("wbgravbc", DevExeSpace(), 0, (pmbp->nmb_thispack-1), 0, n3m1, 0, n2m1,
//        KOKKOS_LAMBDA(int m, int k, int j) {
//          if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::reflect) {
//            for (int i=0; i<ng; ++i) {
//              phicc0(m,k,j,ie+i+1) = phicc0(m,k,j,ie-i);
//              u0_(m,IDN,k,j,ie+i+1) = u0_(m,IDN,k,j,ie-i);
//              u0_(m,IEN,k,j,ie+i+1) = u0_(m,IEN,k,j,ie-i);
//              phi0_x1f(m,k,j,ie+i+1) = phicc0(m,k,j,ie);
//            }
//          }
//        });
    }
    if (use_wellbalance_static) {
        int &ng = indcs.ng;
        int n1m1 = indcs.nx1 + 2*ng - 1;
        int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
        int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
        par_for("wbcc", DevExeSpace(), 0,(pmbp->nmb_thispack-1), 0, n3m1, 0, n2m1, 0, n1m1,
        KOKKOS_LAMBDA(int m, int k, int j, int i) {
            
          Real &x1min = size.d_view(m).x1min;
          Real &x1max = size.d_view(m).x1max;
          int nx1 = indcs.nx1;
          Real &x2min = size.d_view(m).x2min;
          Real &x2max = size.d_view(m).x2max;
          int nx2 = indcs.nx2;
          Real &x3min = size.d_view(m).x3min;
          Real &x3max = size.d_view(m).x3max;
          int nx3 = indcs.nx3;
            
          Real x1v, x2v, x3v;
          if (use_spherical_polar) {
            x1v = x1v_(m,i);
            x2v = x2v_(m,j);
            x3v = x3v_(m,k);
          } else {
            x1v = CellCenterX(i-is, nx1, x1min, x1max);
            ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
            if (use_cubed_sphere_) x1v = x1v_(m,i);   // stored centroid/face, as sp
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          }
          Real r = x1v;
          if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
            
          Real pwb, denwb;
          get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
          u0wb(m,IDN,k,j,i) = denwb;
          u0wb(m,IM1,k,j,i) = 0.0;
          u0wb(m,IM2,k,j,i) = 0.0;
          u0wb(m,IM3,k,j,i) = 0.0;
          u0wb(m,IEN,k,j,i) = EintFromP(eos, igm1, denwb, pwb);
          if (use_etotgrav) {
              Real theta;
              if (use_spherical_polar) {
                theta = x2v;
              } else if (use_cubed_sphere_) {
                Real lam_, phi_;
                CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam_, phi_);
              } else {
                theta = -(x3v*iap-M_PI/2.0);
              }
              Real phicc = TotPotAt(grav_acc, ap, r, x1v, grav_pmass, om2_, theta);
              u0wb(m,IEN,k,j,i) += denwb*phicc;
          }
          w0wb(m,IDN,k,j,i) = denwb;
          w0wb(m,IM1,k,j,i) = 0.0;
          w0wb(m,IM2,k,j,i) = 0.0;
          w0wb(m,IM3,k,j,i) = 0.0;
          w0wb(m,IEN,k,j,i) = EintFromP(eos, igm1, denwb, pwb);
        });
    }

    // the cubed-sphere tangent-pair trig, for the orthonormal-frame bcc below
    auto &ccell_ic_ = pmbp->pcoord->cos_cell;
    auto &scell_ic_ = pmbp->pcoord->sin_cell;
    // initialize magnetic fields if MHD. One-off like the initial condition above: on a
    // restart b0/bcc0 are read from file.
    if (!restart && pmbp->pmhd != nullptr) {
      // Base field in CODE units (Heaviside-Lorentz, p_mag = B^2/2): problem/bbot, or
      // problem/bbot_gauss converted by ReadHotJupiterBbot (pgen.cpp), which also prints
      // it in both units.  B_r = bbot cos(theta) (r0/r)^3: bbot is the polar field at r0.
      // default changed 09-30 (user), was none (fatal if absent)
      Real bbot = hot_jupiter ? hot_jupiter_param.bbot
                              : pin->GetOrAddReal("problem","bbot",0.0);
      // problem/remap_bfile (see below) REPLACES the seed: with the key set the dipole is
      // not applied (bbot = 0 here), so the field is the remapped one only and the
      // energy of the active cells is E_file + e_mag(remap) with no add/subtract of a
      // dipole.  Read only when named, so a run without it embeds the same input.
      const std::string remap_bfile =
          pin->DoesParameterExist("problem","remap_bfile") ?
          pin->GetString("problem","remap_bfile") : std::string("");
      if (!remap_bfile.empty() && bbot != 0.0) {
        if (global_variable::my_rank == 0) {
          std::cout << "deep_hot_jupiter_rt: problem/remap_bfile is set, so the dipole "
                    << "seed (problem/bbot or bbot_gauss, bbot = " << bbot
                    << " [code]) is NOT applied; B is the remapped field only"
                    << std::endl;
        }
        bbot = 0.0;
      }
      auto &b0 = pmbp->pmhd->b0;
      auto &bcc0 = pmbp->pmhd->bcc0;
      par_for("pgen_b0", DevExeSpace(), 0,(pmbp->nmb_thispack-1),0, n3m1, 0, n2m1, 0, n1m1, //ks,ke,js,je,is,ie,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
          if (use_spherical_polar) {
            Real x1v = x1v_(m,i);
            Real x2v = x2v_(m,j);
            Real x3v = x3v_(m,k);
            Real x1fl = x1f_(m,i);
            Real x2fl = x2f_(m,j);
            Real x3fl = x3f_(m,k);
            Real x1fr = x1f_(m,i+1);
            Real x2fr = x2f_(m,j+1);
            Real x3fr = x3f_(m,k+1);
              
            Real A1 = 0.0;
            Real A2 = 0.0;
            Real A2ip = 0.0;
            Real A2kp = 0.0;
            Real A2ipkp = 0.0;
            Real A3 = 0.5*bbot*r0*sin(x2fl)/SQR(x1fl/r0);
            Real A3ip = 0.5*bbot*r0*sin(x2fl)/SQR(x1fr/r0);
            Real A3jp = 0.5*bbot*r0*sin(x2fr)/SQR(x1fl/r0);
            Real A3ipjp = 0.5*bbot*r0*sin(x2fr)/SQR(x1fr/r0);
//            Real A3 = 0.5*bbot*x1fl*sin(x2fl);
//            Real A3ip = 0.5*bbot*x1fr*sin(x2fl);
//            Real A3jp = 0.5*bbot*x1fl*sin(x2fr);
//            Real A3ipjp = 0.5*bbot*x1fr*sin(x2fr);
              
            b0.x1f(m,k,j,i) = (dxe3(m,k,j+1,i)*A3jp - dxe3(m,k,j,i)*A3)/area1(m,k,j,i) - (dxe2(m,k+1,j,i)*A2kp - dxe2(m,k,j,i)*A2)/area1(m,k,j,i);
            if ((mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar && j == js) || (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar && j == je+1)) {
              b0.x2f(m,k,j,i) = - ((x1fr*(x3fr-x3fl))*A3ip - (x1fl*(x3fr-x3fl))*A3) / (0.5*(SQR(x1fr)-SQR(x1fl))*(x3fr-x3fl));
            } else {
              b0.x2f(m,k,j,i) = - (dxe3(m,k,j,i+1)*A3ip - dxe3(m,k,j,i)*A3)/area2(m,k,j,i);
            }
            b0.x3f(m,k,j,i) = (dxe2(m,k,j,i+1)*A2ip - dxe2(m,k,j,i)*A2)/area3(m,k,j,i);
            Real b0x1fip = (dxe3(m,k,j+1,i+1)*A3ipjp - dxe3(m,k,j,i+1)*A3ip)/area1(m,k,j,i+1) - (dxe2(m,k+1,j,i+1)*A2ipkp - dxe2(m,k,j,i+1)*A2ip)/area1(m,k,j,i+1);
            Real b0x2fjp;
            if ((mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar && j == js-1) || (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar && j == je)) {
              b0x2fjp = - ((x1fr*(x3fr-x3fl))*A3ipjp - (x1fl*(x3fr-x3fl))*A3jp) / (0.5*(SQR(x1fr)-SQR(x1fl))*(x3fr-x3fl));
            } else {
              b0x2fjp = - (dxe3(m,k,j+1,i+1)*A3ipjp - dxe3(m,k,j+1,i)*A3jp)/area2(m,k,j+1,i);
            }
            Real b0x3fkp = (dxe2(m,k+1,j,i+1)*A2ipkp - dxe2(m,k+1,j,i)*A2kp)/area3(m,k+1,j,i);
            if (i==n1m1) b0.x1f(m,k,j,i+1) = b0x1fip;
            if (j==n2m1) b0.x2f(m,k,j+1,i) = b0x2fjp;
            if (k==n3m1) b0.x3f(m,k+1,j,i) = b0x3fkp;
              
            Real lw, rw;
            lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
            rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
            bcc0(m,IBX,k,j,i) = lw*b0.x1f(m,k,j,i) + rw*b0x1fip;
            lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
            rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
            bcc0(m,IBY,k,j,i) = lw*b0.x2f(m,k,j,i) + rw*b0x2fjp;
            lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
            rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
            bcc0(m,IBZ,k,j,i) = lw*b0.x3f(m,k,j,i) + rw*b0x3fkp;
              
            u0_(m,IEN,k,j,i) += 0.5*(SQR(bcc0(m,IBX,k,j,i))+SQR(bcc0(m,IBY,k,j,i))+SQR(bcc0(m,IBZ,k,j,i)));
          } else if (use_cubed_sphere_) {
            // CUBED SPHERE.  The same dipole, built the same way: each face value is the
            // circulation of A around that face's own edges divided by the face area,
            // using the very `dxedge` and `area` arrays mhd_ct.cpp consumes.  That makes
            // div B vanish to ROUND-OFF rather than to truncation -- CT preserves the
            // divergence it is handed, so a field projected face by face would keep a
            // small permanent monopole on this grid.  Follows cs_test iprob = 12.
            //
            // A = A_phi phihat with A_phi = 0.5*bbot*r0*sin(theta)/(r/r0)^2, which is the
            // spherical-polar expression above.  In Cartesian components
            // sin(theta)*phihat is (-y,x,0)/r, so no angle is ever formed and the axis
            // is not special:
            //   A = 0.5*bbot*r0/(rf/r0)^2 * (-qhat_y, qhat_x, 0).
            // A.rhat = 0 identically, so the RADIAL edge potential vanishes and only the
            // two tangential edge components are needed.  Unlike spherical polar, where
            // phihat IS the x3 direction, both panel tangents see the field here, so the
            // x2-edge component is not zero.
            const int p = mbpanel_.d_view(m);
            // The component of A on an edge's own UNIT tangent, which is what dxedge*A
            // consumes.  `along_xi` selects the x2 edge (radial face, xi CENTRE, eta
            // face) over the x3 edge (radial face, xi face, eta CENTRE).  x2/x3 are PANEL
            // coordinates on [-1,1], not angles: the angle is pi/4 times the coordinate.
            auto Aedge = [&](const int ii, const int jj, const int kk,
                             const bool along_xi) {
              const Real rf = x1f_(m,ii);
              const Real xi = 0.25*M_PI*(along_xi ? x2v_(m,jj) : x2f_(m,jj));
              const Real et = 0.25*M_PI*(along_xi ? x3f_(m,kk) : x3v_(m,kk));
              Real qh[3], e1[3], e2[3];
              cubed_sphere::PanelToCart(p, xi, et, qh);
              cubed_sphere::PanelTangents(p, xi, et, e1, e2);
              const Real aphi = 0.5*bbot*r0/SQR(rf/r0);
              const Real ax = -aphi*qh[1];
              const Real ay =  aphi*qh[0];
              const Real *t = along_xi ? e1 : e2;
              const Real tn = sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);
              return (ax*t[0] + ay*t[1])/tn;
            };
            b0.x1f(m,k,j,i) =
                (dxe3(m,k,j+1,i)*Aedge(i,j+1,k,false) - dxe3(m,k,j,i)*Aedge(i,j,k,false)
               - dxe2(m,k+1,j,i)*Aedge(i,j,k+1,true) + dxe2(m,k,j,i)*Aedge(i,j,k,true))
                /area1(m,k,j,i);
            b0.x2f(m,k,j,i) =
                -(dxe3(m,k,j,i+1)*Aedge(i+1,j,k,false)
                - dxe3(m,k,j,i)*Aedge(i,j,k,false))/area2(m,k,j,i);
            b0.x3f(m,k,j,i) =
                (dxe2(m,k,j,i+1)*Aedge(i+1,j,k,true)
               - dxe2(m,k,j,i)*Aedge(i,j,k,true))/area3(m,k,j,i);
            // the far faces of the last cell in each direction, which no cell owns
            Real b0x1fip = (dxe3(m,k,j+1,i+1)*Aedge(i+1,j+1,k,false)
                          - dxe3(m,k,j,i+1)*Aedge(i+1,j,k,false)
                          - dxe2(m,k+1,j,i+1)*Aedge(i+1,j,k+1,true)
                          + dxe2(m,k,j,i+1)*Aedge(i+1,j,k,true))/area1(m,k,j,i+1);
            Real b0x2fjp = -(dxe3(m,k,j+1,i+1)*Aedge(i+1,j+1,k,false)
                           - dxe3(m,k,j+1,i)*Aedge(i,j+1,k,false))/area2(m,k,j+1,i);
            Real b0x3fkp = (dxe2(m,k+1,j,i+1)*Aedge(i+1,j,k+1,true)
                          - dxe2(m,k+1,j,i)*Aedge(i,j,k+1,true))/area3(m,k+1,j,i);
            if (i==n1m1) b0.x1f(m,k,j,i+1) = b0x1fip;
            if (j==n2m1) b0.x2f(m,k,j+1,i) = b0x2fjp;
            if (k==n3m1) b0.x3f(m,k+1,j,i) = b0x3fkp;

            // bcc must be formed the way ConsToPrim forms it on THIS grid, or the
            // magnetic energy added here and the one subtracted on the first C2P
            // disagree.  On the cubed sphere that is GnomonicEquiangleRaiseVelMHD: the
            // radial component interpolated to the centroid x1v, and the tangential pair
            // rotated into the ORTHONORMAL frame, B.e_xi = (b_xi + c b_eta)/s with the
            // eta slot unchanged.  The plain averages used before are components on the
            // non-orthogonal tangent pair, whose sum of squares differs from |B|^2 by
            // (c^2 (b2^2 + b3^2) + 2 c b2 b3)/s^2 -- a few percent of the internal
            // energy in the top cells, where the gas energy is ~1 erg/cm^3.
            {
              const Real lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
              const Real rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
              bcc0(m,IBX,k,j,i) = lw*b0.x1f(m,k,j,i) + rw*b0x1fip;
            }
            {
              const Real cg = ccell_ic_(m,k,j), sg = scell_ic_(m,k,j);
              const Real by_n = 0.5*(b0.x2f(m,k,j,i) + b0x2fjp);
              const Real bz_n = 0.5*(b0.x3f(m,k,j,i) + b0x3fkp);
              bcc0(m,IBY,k,j,i) = (by_n + cg*bz_n)/sg;
              bcc0(m,IBZ,k,j,i) = bz_n;
            }

            u0_(m,IEN,k,j,i) += 0.5*(SQR(bcc0(m,IBX,k,j,i))
                                   + SQR(bcc0(m,IBY,k,j,i))
                                   + SQR(bcc0(m,IBZ,k,j,i)));
          }
      });

      // problem/remap_bfile (dhj-remap-mhd, default empty = off): the ACTIVE face fields
      // of a remapped MHD restart, given as edge line integrals of a vector potential
      // in the radial gauge (written by docs/handover/scripts/dhj_remap.py next to
      // remap_file).  b = curl/area with the code's own areas, exactly as the dipole
      // above, so div B vanishes to round-off on this grid.  File: 8-byte magic
      // "RADREMB1", int32 nmb_total, nx1, nx2, nx3; then per gid int32 lx1, lx2, lx3,
      // level, panel and double L2[nx3+1][nx2][nx1+1] (the x2 edges, along xi at the
      // eta faces) and L3[nx3][nx2+1][nx1+1] (the x3 edges, along eta at the xi faces).
      // The dipole seed is skipped when this key is set (bbot = 0 above), so the swap
      // of magnetic energy below removes zero; ghosts hold B = 0 until the boundary
      // exchange and the physical boundaries fill them.  Cubed sphere only.
      if (!remap_bfile.empty()) {
        if (!use_cubed_sphere_ || !use_remap_) {
          std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/remap_bfile "
                    << "needs the cubed sphere and problem/remap_file" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        const int bn1 = indcs.nx1, bn2 = indcs.nx2, bn3 = indcs.nx3;
        DvceArray4D<Real> l2_("remap_l2", pmbp->nmb_thispack, bn3+1, bn2, bn1+1);
        DvceArray4D<Real> l3_("remap_l3", pmbp->nmb_thispack, bn3, bn2+1, bn1+1);
        std::ifstream bf(remap_bfile, std::ios::binary);
        char mg[8];
        std::int32_t hd[4];
        bf.read(mg, 8);
        bf.read(reinterpret_cast<char*>(hd), sizeof(hd));
        if (!bf || std::string(mg, 8) != "RADREMB1" || hd[0] != pmy_mesh_->nmb_total ||
            hd[1] != bn1 || hd[2] != bn2 || hd[3] != bn3) {
          std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/remap_bfile "
                    << remap_bfile << ": unreadable, bad magic or grid mismatch"
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
        const std::size_t n2r = static_cast<std::size_t>(bn3+1)*bn2*(bn1+1);
        const std::size_t n3r = static_cast<std::size_t>(bn3)*(bn2+1)*(bn1+1);
        const std::streamoff recb = 5*sizeof(std::int32_t) + (n2r + n3r)*sizeof(double);
        auto hl2 = Kokkos::create_mirror_view(l2_);
        auto hl3 = Kokkos::create_mirror_view(l3_);
        std::vector<double> buf(n2r + n3r);
        for (int m=0; m<pmbp->nmb_thispack; ++m) {
          const int gid = pmbp->gids + m;
          bf.seekg(8 + sizeof(hd) + gid*recb, std::ios::beg);
          std::int32_t ll[5];
          bf.read(reinterpret_cast<char*>(ll), sizeof(ll));
          bf.read(reinterpret_cast<char*>(buf.data()), (n2r + n3r)*sizeof(double));
          const LogicalLocation &lc = pmy_mesh_->lloc_eachmb[gid];
          if (!bf || ll[0] != lc.lx1 || ll[1] != lc.lx2 || ll[2] != lc.lx3 ||
              ll[3] != lc.level || ll[4] != lc.panel) {
            std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/remap_bfile "
                      << "record " << gid << " does not match this MeshBlock's "
                      << "LogicalLocation" << std::endl;
            std::exit(EXIT_FAILURE);
          }
          std::size_t q = 0;
          for (int k=0; k<=bn3; ++k) {
            for (int j=0; j<bn2; ++j) {
              for (int i=0; i<=bn1; ++i) hl2(m,k,j,i) = buf[q++];
            }
          }
          for (int k=0; k<bn3; ++k) {
            for (int j=0; j<=bn2; ++j) {
              for (int i=0; i<=bn1; ++i) hl3(m,k,j,i) = buf[q++];
            }
          }
        }
        Kokkos::deep_copy(l2_, hl2);
        Kokkos::deep_copy(l3_, hl3);
        // 1. the faces of every active cell (the far faces by the last cell)
        par_for("pgen_b0_remap", DevExeSpace(), 0, (pmbp->nmb_thispack-1), ks, ke, js, je,
                is, ie, KOKKOS_LAMBDA(int m, int k, int j, int i) {
          const int a = i - is, b = j - js, c = k - ks;
          u0_(m,IEN,k,j,i) -= 0.5*(SQR(bcc0(m,IBX,k,j,i)) + SQR(bcc0(m,IBY,k,j,i))
                                 + SQR(bcc0(m,IBZ,k,j,i)));
          b0.x1f(m,k,j,i) = (l3_(m,c,b+1,a) - l3_(m,c,b,a) - l2_(m,c+1,b,a)
                           + l2_(m,c,b,a))/area1(m,k,j,i);
          b0.x2f(m,k,j,i) = -(l3_(m,c,b,a+1) - l3_(m,c,b,a))/area2(m,k,j,i);
          b0.x3f(m,k,j,i) = (l2_(m,c,b,a+1) - l2_(m,c,b,a))/area3(m,k,j,i);
          if (i == ie) {
            b0.x1f(m,k,j,i+1) = (l3_(m,c,b+1,a+1) - l3_(m,c,b,a+1) - l2_(m,c+1,b,a+1)
                               + l2_(m,c,b,a+1))/area1(m,k,j,i+1);
          }
          if (j == je) {
            b0.x2f(m,k,j+1,i) = -(l3_(m,c,b+1,a+1) - l3_(m,c,b+1,a))/area2(m,k,j+1,i);
          }
          if (k == ke) {
            b0.x3f(m,k+1,j,i) = (l2_(m,c+1,b,a+1) - l2_(m,c+1,b,a))/area3(m,k+1,j,i);
          }
        });
        // 2. bcc and the magnetic energy of the active cells, as the dipole branch
        par_for("pgen_bcc_remap", DevExeSpace(), 0, (pmbp->nmb_thispack-1), ks, ke, js,
                je, is, ie, KOKKOS_LAMBDA(int m, int k, int j, int i) {
          const Real lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
          const Real rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
          bcc0(m,IBX,k,j,i) = lw*b0.x1f(m,k,j,i) + rw*b0.x1f(m,k,j,i+1);
          const Real cg = ccell_ic_(m,k,j), sg = scell_ic_(m,k,j);
          const Real by_n = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
          const Real bz_n = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
          bcc0(m,IBY,k,j,i) = (by_n + cg*bz_n)/sg;
          bcc0(m,IBZ,k,j,i) = bz_n;
          u0_(m,IEN,k,j,i) += 0.5*(SQR(bcc0(m,IBX,k,j,i)) + SQR(bcc0(m,IBY,k,j,i))
                                 + SQR(bcc0(m,IBZ,k,j,i)));
        });
        if (global_variable::my_rank == 0) {
          std::cout << "deep_hot_jupiter_rt: problem/remap_bfile " << remap_bfile
                    << " read (active face fields replaced)" << std::endl;
        }
      }

      // GATE ON THE CONSTRUCTION, printed once.  Two numbers, because div B alone is
      // satisfied by a ZERO field -- which is exactly what the cubed sphere silently had
      // before this branch existed.  B is the discrete curl of A on the very edges CT
      // uses, so sum(area*B) over a cell must vanish to ROUND-OFF at t=0, and |B|max must
      // be of order bbot.
      {
        auto b1 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), b0.x1f);
        auto b2 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), b0.x2f);
        auto b3 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), b0.x3f);
        auto a1 = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                      pmbp->pcoord->area.x1f);
        auto a2 = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                      pmbp->pcoord->area.x2f);
        auto a3 = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                      pmbp->pcoord->area.x3f);
        Real dvb = 0.0, bscale = 0.0;
        for (int m=0; m<pmbp->nmb_thispack; ++m) {
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je; ++j) {
              for (int i=is; i<=ie; ++i) {
                Real fl = b1(m,k,j,i+1)*a1(m,k,j,i+1) - b1(m,k,j,i)*a1(m,k,j,i)
                        + b2(m,k,j+1,i)*a2(m,k,j+1,i) - b2(m,k,j,i)*a2(m,k,j,i)
                        + b3(m,k+1,j,i)*a3(m,k+1,j,i) - b3(m,k,j,i)*a3(m,k,j,i);
                Real amx = fmax(a1(m,k,j,i+1), fmax(a2(m,k,j+1,i), a3(m,k+1,j,i)));
                Real bmg = fmax(fabs(b1(m,k,j,i)),
                                fmax(fabs(b2(m,k,j,i)), fabs(b3(m,k,j,i))));
                bscale = fmax(bscale, bmg);
                dvb = fmax(dvb, fabs(fl)/fmax(bmg*amx, 1.0e-300));
              }
            }
          }
        }
#if MPI_PARALLEL_ENABLED
        Real dd_[2] = {dvb, bscale};
        MPI_Allreduce(MPI_IN_PLACE, dd_, 2, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
        dvb = dd_[0];  bscale = dd_[1];
#endif
        if (global_variable::my_rank == 0) {
          std::cout << "deep_hot_jupiter_rt: initial B -- max |sum area*B|"
                    << " / (|B| max area) = " << dvb << ", max |B| = " << bscale
                    << " [code] = " << bscale*std::sqrt(4.0*M_PI) << " G (bbot = "
                    << bbot << " [code] = " << bbot*std::sqrt(4.0*M_PI) << " G)"
                    << std::endl;
        }
      }
    }

    // On a restart u0 and b0 are restored exactly, ghost zones included, but bcc0 is a
    // DERIVED array that the restart file does not carry, so it starts at zero. The
    // inner-x1 user boundary swaps the magnetic energy of its ghost cells by subtracting
    // 0.5*bcc0^2 and adding it back from the current face fields; with bcc0 still zero it
    // subtracts nothing and so double counts the magnetic energy on the very first step.
    // Fill bcc0 from the restored face fields here, which is exactly the state the
    // initialisation above leaves on a fresh start. The interpolation must match the one
    // ConsToPrim uses, hence the split on use_spherical_polar.
    if (restart && pmbp->pmhd != nullptr) {
      auto &b0 = pmbp->pmhd->b0;
      auto &bcc0 = pmbp->pmhd->bcc0;
      par_for("pgen_bcc_restart", DevExeSpace(), 0, (pmbp->nmb_thispack-1),
              0, n3m1, 0, n2m1, 0, n1m1,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        if (use_spherical_polar) {
          Real lw, rw;
          lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
          rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
          bcc0(m,IBX,k,j,i) = lw*b0.x1f(m,k,j,i) + rw*b0.x1f(m,k,j,i+1);
          lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
          rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
          bcc0(m,IBY,k,j,i) = lw*b0.x2f(m,k,j,i) + rw*b0.x2f(m,k,j+1,i);
          lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
          rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
          bcc0(m,IBZ,k,j,i) = lw*b0.x3f(m,k,j,i) + rw*b0.x3f(m,k+1,j,i);
        } else if (use_cubed_sphere_) {
          // as GnomonicEquiangleRaiseVelMHD: centroid-weighted radial, orthonormal frame
          Real lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
          Real rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
          bcc0(m,IBX,k,j,i) = lw*b0.x1f(m,k,j,i) + rw*b0.x1f(m,k,j,i+1);
          const Real cg = ccell_ic_(m,k,j), sg = scell_ic_(m,k,j);
          const Real by_n = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
          const Real bz_n = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
          bcc0(m,IBY,k,j,i) = (by_n + cg*bz_n)/sg;
          bcc0(m,IBZ,k,j,i) = bz_n;
        } else {
          bcc0(m,IBX,k,j,i) = 0.5*(b0.x1f(m,k,j,i) + b0.x1f(m,k,j,i+1));
          bcc0(m,IBY,k,j,i) = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
          bcc0(m,IBZ,k,j,i) = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
        }
      });
    }

  // problem/rt_test_rho / rt_test_T / rt_test_dT / rt_test_nwave: the test initial
  // state (see the note at rt_test_freeze).  Default off; from scratch only.
  {
    const Real trho = pin->GetOrAddReal("problem","rt_test_rho",0.0);
    const Real tTk = pin->GetOrAddReal("problem","rt_test_T",0.0);
    const Real tdT = pin->GetOrAddReal("problem","rt_test_dT",0.0);
    const Real tnw = pin->GetOrAddReal("problem","rt_test_nwave",0.0);
    // rt_test_rstep > 0 (test only): rho = rt_test_rho_top above r = rt_test_rstep (the
    // sharp-photosphere check of tests_ck_sph/closure_0930); needs rt_test_rho > 0.
    const Real trst = pin->GetOrAddReal("problem","rt_test_rstep",0.0);
    const Real trtop = pin->GetOrAddReal("problem","rt_test_rho_top",0.0);
    if (!restart && (trho > 0.0 || tTk > 0.0 || tdT != 0.0)) {
      if (pmbp->phydro == nullptr || use_cubed_sphere_) {
        std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/rt_test_rho, _T "
                  << "and _dT need a hydro run on a spherical-polar mesh" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      auto u0t = pmbp->phydro->u0;
      auto eost = eos;
      const bool etgt = use_etotgrav;
      auto phit = phicc0;
      auto &x1vt = pmbp->pcoord->x1v;
      const Real r0t = pmy_mesh_->mesh_size.x1min;
      const Real r1t = pmy_mesh_->mesh_size.x1max;
      const Real tcode = (tTk > 0.0) ? tTk/eost.temp_cgs : 0.0;
      par_for("pgen_rt_test_init", DevExeSpace(), 0, (pmbp->nmb_thispack-1),
              ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        const Real rho0 = u0t(m,IDN,k,j,i);
        const Real ph = etgt ? phit(m,k,j,i) : 0.0;
        const Real eold = EintFromCons(u0t, m, k, j, i, 0.0, false, etgt, ph);
        Real ke0 = u0t(m,IEN,k,j,i) - eold - (etgt ? rho0*ph : 0.0);
        Real rho = rho0;
        if (trho > 0.0) {
          rho = (trst > 0.0 && x1vt(m,i) > trst) ? trtop : trho;
          u0t(m,IDN,k,j,i) = rho;
          u0t(m,IM1,k,j,i) = 0.0;
          u0t(m,IM2,k,j,i) = 0.0;
          u0t(m,IM3,k,j,i) = 0.0;
          ke0 = 0.0;
        }
        Real t = (tcode > 0.0) ? tcode : eost.Temperature(rho0, eold);
        const Real xi = (x1vt(m,i) - r0t)/(r1t - r0t);
        t *= 1.0 + tdT*((tnw > 0.0) ? sin(tnw*M_PI*xi) : 1.0);
        u0t(m,IEN,k,j,i) = eost.EnergyFromTemperature(rho, t) + ke0
                         + (etgt ? rho*ph : 0.0);
      });
      if (global_variable::my_rank == 0) {
        std::cout << "deep_hot_jupiter_rt: rt_test initial state: rho " << trho
                  << " T " << tTk << " dT " << tdT << " nwave " << tnw << std::endl;
      }
    }
  }

  // the correlated-k Rosseland table for the radiative diffusion (no-op unless
  // <mhd>/rad_kappa_src = table)
  ck_build_rosseland_table((pmbp->pmhd != nullptr) ? pmbp->pmhd->pcond
                                                  : pmbp->phydro->pcond);
  {
    MeshBlockPack *pmbp_ = pmy_mesh_->pmb_pack;
    Conduction *pc = (pmbp_->pmhd != nullptr) ? pmbp_->pmhd->pcond : pmbp_->phydro->pcond;
    if (pc != nullptr && pc->iso_cond_type.compare("radiative") == 0) {
      // the tau blend hands the deep column over between the diffusion and the
      // correlated-k two-stream; the grey scheme has no cut and would double count
      if (pc->rad_tau_mode && !rt_ck) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "rad_tau_hi > 0 (the optical-depth blend) needs "
                  << "problem/rt_ck = true" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      // the internal flux enters on exactly ONE path: through the conduction wall
      // (rad_flux_inner, blend on) or as the two-stream's intensity at the cut
      // (ck_int_at_cut, blend off).  Both at once would inject it twice.
      if (!pc->rad_tau_mode && rt_int_at_cut && pc->rad_flux_inner != 0.0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "rad_flux_inner != 0 without the optical-depth blend "
                  << "injects the internal flux twice: set problem/ck_int_at_cut = false "
                  << "or rad_tau_hi > 0" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      // rad_flux_inner < 0: the planet's internal flux sigma T_int^4 through the wall
      if (pc->rad_flux_inner < 0.0) {
        Real Tint;
        get_Tint(hot_jupiter_param.Teq, Tint);
        const Real fint = 5.670374419e-5*SQR(SQR(Tint));
        pc->rad_flux_inner = fint/(pmbp_->punit->pressure_cgs()
                                   *pmbp_->punit->velocity_cgs());
        if (global_variable::my_rank == 0) {
          std::cout << "  radiative diffusion: wall flux sigma T_int^4 = " << fint
                    << " erg/s/cm^2 (T_int = " << Tint << " K)" << std::endl;
        }
      }
    }
  }

  // problem/rt_test_freeze: the frozen state is the state the problem generator (or the
  // restart) leaves behind
  if (rt_test_freeze) {
    auto &u0f = pmy_mesh_->pmb_pack->phydro->u0;
    rt_frz_ptr = new DvceArray5D<Real>("rt_frz", u0f.extent(0), u0f.extent(1),
                                       u0f.extent(2), u0f.extent(3), u0f.extent(4));
    Kokkos::deep_copy(*rt_frz_ptr, u0f);
    // and the temperature cache the RT takes its EOS guess from (TGuess(wtemp)): it is
    // restored with the state, so that a horizontally uniform state stays bitwise
    // uniform instead of inheriting the per-column guesses of the discarded RK stages
    auto &wt = pmy_mesh_->pmb_pack->phydro->wtemp;
    rt_frzT_ptr = new DvceArray4D<Real>("rt_frzT", wt.extent(0), wt.extent(1),
                                        wt.extent(2), wt.extent(3));
    Kokkos::deep_copy(*rt_frzT_ptr, wt);
  }

  return;
}
