//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_params.cpp
//! \brief keys of the scheme (called by the constructor), and the checks of the keys an
//! input names.
//!
//!  called from: rad_vet.cpp

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_args.hpp"
#include "utils/recipe.hpp"
#include "utils/repro_sum.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CheckKeys
//! \brief checks the <rad_vet> / <time> keys an input names under scheme = clean: a
//! replaced key names its clean counterpart, unread keys are listed once, and a key that
//! would change the physics stops the run unless off.

void RadiationVet::CheckKeys(ParameterInput *pin) {
  // block/key (* ends a prefix): nullptr read, "" no effect, "!" fatal unless off, "x"
  // removed (fatal when named), "=v" fatal unless v (old restart value), "+v" too, with a
  // NOTE where a restart file of an older code has no such key; else clean key
  static const char *tab[][2] = {
    {"rad_vet/implicit_bc_advect", ""}, {"rad_vet/implicit_bc_*", nullptr},
    // inner-flux time hooks: removed
    {"rad_vet/implicit_flux_x1min_mod_*", "x"},
    {"rad_vet/implicit_flux_x1min_ramp_*", "x"},
    {"rad_vet/implicit_flux_x1min_table*", "x"},
    {"rad_vet/implicit_flux_x1m*", nullptr},
    {"rad_vet/implicit_ebath_*", nullptr}, {"rad_vet/implicit_beam*", nullptr},
    {"rad_vet/implicit_cfl", nullptr},
    {"rad_vet/implicit_mr_*", nullptr}, {"rad_vet/implicit_det_trace", nullptr},
    {"rad_vet/implicit_tol", "clean_tol_nl"}, {"rad_vet/implicit_lin_tol", "clean_tol"},
    // read by ImplicitInit; a value other than the clean key's prints a NOTE
    {"rad_vet/implicit_maxit", nullptr},
    {"rad_vet/implicit_tol_rel_max", nullptr},
    {"rad_vet/implicit_lin_maxit", nullptr},
    {"rad_vet/implicit_closure_lag", "clean_closure"},
    {"rad_vet/time2_vet", "clean_closure"}, {"rad_vet/time2_vet_col", "clean_closure"},
    {"rad_vet/time2_vet_sc", "clean_closure"},
    {"rad_vet/time2_vet_extrap", "clean_closure"},
    {"rad_vet/advect_split", ""}, {"rad_vet/time2_one_pass_safety", ""},
    {"rad_vet/time2_lin_tol*", "clean_tol"},
    {"rad_vet/implicit_ghost_x1min", "!"}, {"rad_vet/ic_rad*", "clean_init"},
    {"rad_vet/implicit_predictor", "=step"}, {"rad_vet/implicit_predictor_tcap", "=1.1"},
    {"rad_vet/implicit_resid_fatal", nullptr},
    // the old scheme's predictor order: the layout of its restart files (Time2RstNchWant)
    {"rad_vet/implicit_predictor_order", nullptr},
    {"rad_vet/implicit_hold", nullptr}, {"rad_vet/implicit_pos_floor", nullptr},
    {"rad_vet/implicit_enthalpy", "=plm"}, {"rad_vet/implicit_pos_gas", "+true"},
    {"rad_vet/implicit_pos_gas_frac", "=1e-3"},
    {"rad_vet/implicit_pos_gas_cons", "+true"},
    {"rad_vet/implicit_pos_floor_cons", "=true"},
    {"rad_vet/implicit_pos_floor_gas_keep", "+1"},
    {"rad_vet/implicit_am_cs", "=true"}, {"rad_vet/be_legacy", "=false"},
    {"rad_vet/implicit_*", ""},
    // hesdirk2 stage formation (Time2FormStage / Time2Stash) and dt_retry: read;
    // ImplicitSolve's acceptance, posfix and KE cap: none
    {"rad_vet/time2_k1_reset", "!"},   // K1 reset of held cells: removed
    {"rad_vet/time2_k1_redist", "=0"}, {"rad_vet/time2_k1_stiff", "=0"},
    {"rad_vet/time2_face_fct", "+true"},
    {"time/dt_retry_max", ""},
    // a rejected stage: always the backward-Euler redo of the cycle at its dt
    {"rad_vet/time2_reject", ""}, {"time/dt_retry_t2_be", "=0"},
    {"time/dt_retry_fac_t2", ""}, {"time/dt_retry_t2be_now", ""},
    {"rad_vet/time2_accept_fallback", ""}, {"rad_vet/time2_posfix*", ""},
    {"rad_vet/time2_inc_kecap", ""}, {"rad_vet/time2_pred_cap", ""},
    // implicit solve's position limiters, Zalesak round cap and x2 vet fields: none
    {"rad_vet/time2_adm_poslim", ""}, {"rad_vet/time2_dtr_poslim", ""},
    {"rad_vet/zalesak_cap_rounds", ""},
    {"rad_vet/time2_vet_x2_fields", ""},
    // cell-layout keys: removed
    {"rad_vet/clean_alpha", "x"}, {"rad_vet/clean_recon", "x"},
    {"rad_vet/clean_face_model", "x"}, {"rad_vet/clean_consistency_weight", "x"},
    {"rad_vet/clean_psi_eps", "x"}, {"rad_vet/clean_pos_*", "x"},
    {"rad_vet/clean_diag_cells", "x"}, {"rad_vet/clean_maxpos", "x"},
    {"rad_vet/clean_vfreeze", "x"}, {"rad_vet/clean_gas_exact", "x"},
    {"rad_vet/clean_yres", "x"}, {"rad_vet/clean_nl_efloor", "x"},
    {"rad_vet/clean_pc", "x"}, {"rad_vet/clean_pc_sweeps", "x"},
    {"rad_vet/clean_pc_factor*", "x"}, {"rad_vet/clean_pc_reuse", "x"},
    {"rad_vet/clean_dissip_power", "x"},
    {"rad_vet/clean_line_teams", "x"}, {"rad_vet/clean_line_solver", "x"},
    {"rad_vet/clean_line_chunk", "x"},
    {"rad_vet/clean_fv_*", "x"},   // FV rays: removed
    // explicit transport's face-flux and source options: only their old defaults
    {"rad_vet/thick_flux", "=none"}, {"rad_vet/ap_form", "=alpha2"},
    {"rad_vet/reconstruct", "=plm"}, {"rad_vet/split_vel", "=recon"},
    {"rad_vet/source_form", "=full"}, {"rad_vet/f_source", "=cell"},
    {"rad_vet/dbg_tensor", "=none"}, {"rad_vet/dbg_tensor_tilt", ""},
    // the one value left
    {"rad_vet/scheme", "=clean"}, {"rad_vet/clean_layout", "=staggered"},
    {"rad_vet/clean_rays", "=sc"}, {"rad_vet/clean_halo_direct", "=true"},
    {"rad_vet/clean_cs_lz_abs", "=2"}, {"rad_vet/clean_stag_csvec_radial", "=true"},
    {"rad_vet/clean_stag_flim", "=false"}, {"rad_vet/clean_pass_damp", "=0"},
    {"rad_vet/clean_conv_est_f", "=false"}, {"rad_vet/clean_retry_budget", "=true"},
    {"rad_vet/vet_gd_twin", "=false"},
    // study keys: removed
    {"rad_vet/clean_stag_colflux", "=false"}, {"rad_vet/clean_stag_fwb_cart", "=false"},
    {"rad_vet/clean_stag_fwb_line", "=false"},
    {"rad_vet/clean_stag_fwb_stage", "=false"}, {"rad_vet/clean_stag_fwb_t2", "=false"},
    {"rad_vet/clean_stag_rayratio", "=false"}, {"rad_vet/clean_vel_settled", "=false"},
    {"rad_vet/clean_xt_cons", "=false"}, {"rad_vet/clean_xt_norm", "=false"},
    {"rad_vet/dbg_energy_tally", "=false"},
    {"rad_vet/doppler_face", "=donor"}, {"rad_vet/vet_gd_rotate_stretch", "=0.02"},
    {"rad_vet/vet_rotate_avg_ema_dema", "=false"},
    {"rad_vet/vet_rotate_avg_ema_diag", "=0"}, {"rad_vet/vet_rotate_avg_ema_lc", "=0"},
    {"rad_vet/vet_rotate_avg_ema_tau", "=0"},
    // linear solve and passes: constants
    {"rad_vet/clean_bicg_stall_its", "=0"}, {"rad_vet/clean_bicg_stall_fac", "=0.99"},
    {"rad_vet/clean_init_cnorm", "=1e-10"}, {"rad_vet/clean_lin_cnorm", "=1e-10"},
    {"rad_vet/clean_lin_pass", "+true"}, {"rad_vet/clean_maxpic", "=12"},
    {"rad_vet/clean_nl_norm", "=global"}, {"rad_vet/clean_troot_rtol", "=1e-13"},
    {"rad_vet/clean_gas_newton_switch", "+4"},
    {"rad_vet/clean_gas_newton_switch_min", "+6"}, {"rad_vet/clean_stall_fac", "=0.5"},
    {"rad_vet/clean_gas_secant", "+max_stage2"}, {"rad_vet/clean_pred_stage", "=false"},
    {"rad_vet/clean_pred_xch", "+true"}, {"rad_vet/clean_predict", "=true"},
    // closure: constants
    {"rad_vet/clean_fill_passes", "=8"}, {"rad_vet/clean_fill_cs_remap", "=false"},
    {"rad_vet/clean_ring_jmin", "=true"},
    {"rad_vet/clean_ray_cut_margin", "=10"}, {"rad_vet/clean_ray_cut_taumin", "=0"},
    {"rad_vet/polar_ring_avg", "=false"}, {"rad_vet/clean_w_face_min", "=false"},
    // face coefficients and bracket: constants
    {"rad_vet/clean_bc", "=halfrange"}, {"rad_vet/clean_bc_area", "=true"},
    {"rad_vet/clean_face_area", "=false"}, {"rad_vet/clean_face_hctau", "+true"},
    {"rad_vet/clean_face_jse", "=true"}, {"rad_vet/clean_face_jse_cell", "=false"},
    {"rad_vet/clean_face_opac", "+harmonic"}, {"rad_vet/clean_w_fse", "=true"},
    {"rad_vet/clean_w_fse_tau", "=true"}, {"rad_vet/clean_source_convex", "=true"},
    {"rad_vet/clean_cs_lz_bnd_free", "+true"}, {"rad_vet/clean_stag_eicarry", "+true"},
    {"rad_vet/clean_wb_floor_fail", "=false"}, {"rad_vet/clean_rst_mshare", "=true"},
    {"rad_vet/clean_rot_ring_fill_rst", "=true"},
    // staggered operator: constants; where a rule decides (fwb, freal with
    // clean_consistency; fschur, fsweeps by the mesh; fwb_ring by the x3 wrap) the rule
    {"rad_vet/clean_stag_fwb", nullptr}, {"rad_vet/clean_stag_fwb_ring", nullptr},
    {"rad_vet/clean_stag_freal", nullptr},
    {"rad_vet/clean_stag_fschur", nullptr}, {"rad_vet/clean_stag_fsweeps", nullptr},
    {"rad_vet/clean_gas_tfloor", "=true"}, {"rad_vet/clean_wb_fct", "+true"},
    {"rad_vet/clean_xt_exact", "+cons"}, {"rad_vet/clean_stag_central", "+true"},
    {"rad_vet/clean_stag_upwind_central", ""},
    {"rad_vet/clean_stag_fs_line", "+pcr"},
    {"rad_vet/clean_stag_pcr", "=true"}, {"rad_vet/clean_stag_opac_update", "+true"},
    {"rad_vet/clean_stag_opac_drag", "+true"}, {"rad_vet/clean_stag_lat_sums", "+true"},
    {"rad_vet/clean_stag_mg_color", "+mirror"}, {"rad_vet/clean_stag_tie_gid", "+true"},
    {"rad_vet/clean_stag_xdonor", "+lag"}, {"rad_vet/clean_stag_xsrc", "+true"},
    {"rad_vet/clean_stag_vimp", "=true"}, {"rad_vet/clean_stag_vimp_cs", "+true"},
    // write-back, time scheme and restart writer: constants
    {"rad_vet/time2_cfl_post", "=true"}, {"rad_vet/time2_enth_be", ""},
    {"rad_vet/time2_fref_wacc", "=true"}, {"rad_vet/time2_gtime_mom", "=true"},
    {"rad_vet/zalesak_active", "+true"}, {"rad_vet/zalesak_halo_skip", "=true"},
    {"rad_vet/restart_pack_per_mb", "=true"}, {"rad_vet/cs_seam_avg", ""},
    // rays and rotation: constants; rules (recipe::Rule): the rotation axis (3d on the
    // cubed sphere), the cross-fade and the advancing window (the ring), ema_wb (curved
    // meshes), the global band (full sphere), the slope at a set change (a Cartesian set)
    {"rad_vet/vet_gd_rotate_axis", nullptr}, {"rad_vet/vet_gd_rotate_xfade", nullptr},
    {"rad_vet/vet_rotate_advance", nullptr}, {"rad_vet/vet_rotate_avg_ema_wb", nullptr},
    {"rad_vet/vet_gd_band_global", nullptr}, {"rad_vet/vet_gd_rotate_slope", nullptr},
    {"rad_vet/vet_gd_rotate_slope_frac", nullptr},
    {"rad_vet/vet_gd_rotate_align", "=true"}, {"rad_vet/vet_gd_rotate_t0", "=0"},
    {"rad_vet/vet_gd_rotate_xfade_m", "=0"}, {"rad_vet/vet_rotate_avg_weight", "=j"},
    {"rad_vet/vet_rotate_avg_ema_change", "=0.05"},
    {"rad_vet/vet_rotate_avg_ema_change_j", "=true"},
    {"rad_vet/vet_gd_halo_compact", "+2"}, {"rad_vet/vet_interp", "+auto"},
    {"rad_vet/vet_seg", "+jump"}, {"rad_vet/vet_seg_plateau", "=false"},
    {"rad_vet/vet_scatter_form", "=ma"}, {"rad_vet/vet_scatter_j", "=relaxed"},
    {"rad_vet/vet_source_noesrc_h2be", "=true"}, {"rad_vet/vet_tensor", "+full"},
    {"rad_vet/vet_sc_kernel", "!"},
    {"rad_vet/vet_col_lat_cut_thin", "+1"},
    // early files' names of clean_stag_fs_line
    {"rad_vet/clean_stag_fs_pcr", "=true"}, {"rad_vet/clean_stag_fs_pcr_cs", "=true"},
    // old solver's diagnostics, multi-rate, coupling split, K1 FCT: removed
    {"rad_vet/dbg_energy_audit", "!"}, {"rad_vet/dbg_cs_metric_selftest", "!"},
    {"rad_vet/implicit_det_trace", "!"}, {"rad_vet/dbg_hydro_off", "!"},
    {"rad_vet/time2_k1_fct", "!"}, {"rad_vet/implicit_mr_every", "!"},
    {"rad_vet/coupling_split", "!"},
    // Zalesak fixed-point and local-round variants: only off
    {"rad_vet/zalesak_converge", "!"}, {"rad_vet/zalesak_local_rounds", "!"},
    {"rad_vet/zalesak_colsweep", "!"},
    // old step history: never kept, its restart channels are skipped
    {"rad_vet/clean_skip_time2_hist", "=true"},
    // hesdirk2 start probe: only be
    {"rad_vet/time2_start", "=be"},
    // extrapolated vet_sc builds: only 1
    {"rad_vet/vet_sc_every", "=1"},
    // restart fingerprints (debug hook): removed
    {"rad_vet/dbg_rst_fp*", "x"},
    // explicit transport's sub-cycling: removed
    {"rad_vet/cfl_rad", ""}, {"rad_vet/subcycle", ""},
    // vet_gd twin and symmetric reference: inert here (with vet_gd, a twin named true
    // stops in rad_vet.cpp)
    {"rad_vet/vet_gd_twin*", ""},
    {"rad_vet/vet_gd_symref*", ""},
    {"rad_vet/vet_gd_async", "!"}, {"rad_vet/time2_stage_theta", "!"},
    {"rad_vet/time2_stage_local_be", "!"},
    // study keys: removed
    {"rad_vet/vet_gd_smooth", "!"}, {"rad_vet/vet_dump", "=none"},
    {"rad_vet/vet_dump_every", ""}, {"rad_vet/vet_gd_halo_cache_mb", "!"},
    {"rad_vet/vet_gd_beam", "!"}, {"rad_vet/vet_gd_beam_*", "x"},
    {"rad_vet/vet_gd_beam2_*", "x"}, {"rad_vet/vet_gd_mms", "!"},
    {"rad_vet/vet_gd_hr_split", "=none"}, {"rad_vet/vet_gd_hr_split_*", ""},
    {"rad_vet/vet_gd_hr_latfit", "=-1"}, {"rad_vet/vet_gd_hr_latfit_*", ""},
    {"rad_vet/sp_face_sync", "!"}, {"rad_vet/sp_face_sync_relax", ""},
    // the column formal solution, its lateral correction and
    // the relaxed source of lagged builds: no effect; only the values the rays had
    {"rad_vet/vet_col_ncore", ""}, {"rad_vet/vet_col_nsub", ""},
    {"rad_vet/vet_col_nmu", ""}, {"rad_vet/vet_col_fk_min", ""},
    {"rad_vet/vet_col_reflect_top", ""}, {"rad_vet/vet_col_axis", "=radial"},
    {"rad_vet/vet_col_surface_*", ""},
    {"rad_vet/vet_col_lat", ""}, {"rad_vet/vet_col_lat_offdiag", ""},
    {"rad_vet/vet_col_lat_posfix", ""}, {"rad_vet/vet_every", "=1"},
    {"rad_vet/vet_col_every", "=1"}, {"rad_vet/vet_source", ""},
    {"rad_vet/vet_col_source", ""}, {"rad_vet/vet_gd_replace", ""},
    {"rad_vet/vet_gd_tangential", ""}, {"rad_vet/vet_gd_rotate_slope_reuse", ""},
    {"rad_vet/time2_vet_source*", ""}, {"rad_vet/time2_vet_col_source*", ""},
    {"rad_vet/cs_vertex_mono*", ""},
    // old solver's diagnostics and vet_col build keys: no effect
    {"rad_vet/report_newton_fb", ""}, {"rad_vet/dbg_t2_admiss", ""},
    {"rad_vet/dbg_opac_part", ""}, {"rad_vet/doppler_vel_dbg", ""},
    {"rad_vet/dbg_trans_memory", ""}, {"rad_vet/dbg_cell_lo", ""},
    {"rad_vet/dbg_cell_hi", ""}, {"rad_vet/vet_col_lat_posfix_force", ""},
    {"rad_vet/vet_col_lat_dump", ""}, {"rad_vet/vet_col_lat_dump_every", ""},
    {"rad_vet/polar_ring_avg_lump", "!"}, {"rad_vet/polar_ring_avg_gz", "!"},
    {"rad_vet/time2_vstage", ""}, {"rad_vet/vet_col_order2", ""},
    {"rad_vet/vet_col_team", ""}, {"rad_vet/vet_col_team_size", ""},
    {"rad_vet/vet_col_chunk", ""}, {"rad_vet/vet_col_spread", ""},
    // vet_sc on several MeshBlocks: ray kernel, angle sweep and lag removed
    {"rad_vet/vet_mb_kernel", "x"}, {"rad_vet/vet_mb_tblock*", "x"},
    {"rad_vet/vet_mb_mom_fuse", "x"}, {"rad_vet/vet_mb_mom_batch", "x"},
    {"rad_vet/vet_mb_halo", "x"}, {"rad_vet/vet_mb_overlap", "x"},
    {"rad_vet/vet_mb_angles", "x"}, {"rad_vet/vet_mb_lag", "x"},
    {"rad_vet/vet_mb_agroup", "=1"}, {"rad_vet/vet_mb_agg", "=true"},
    // vet_sc Milne diagnostic: removed
    {"rad_vet/vet_milne", "!"},
    // vet_path = cells: removed (decks and restart files store none)
    {"rad_vet/vet_path", "=none"},
    // vet_seg_cross: removed
    {"rad_vet/vet_seg_cross", "!"},
    // vet_col_lat depth cut: the clean rays always start at the first shell
    {"rad_vet/vet_col_lat_taucut", ""},
    // vet_gd upwind profile shift and dense shell launch: removed
    {"rad_vet/vet_gd_upwind", "=bilinear"}, {"rad_vet/vet_gd_upwind_*", ""},
    {"rad_vet/vet_gd_shell_list", "=true"},
    // cs sweep's pipelined direction groups: removed
    {"rad_vet/vet_gd_cs_xgroups", "=1"},
    // inward-only band ring slabs: removed
    {"rad_vet/vet_gd_band_ring_inward", "!"},
    // staggered increment and cross-flux study forms: removed
    {"rad_vet/clean_stag_kpos", "!"}, {"rad_vet/clean_stag_kpat", "=false"},
    {"rad_vet/clean_stag_kreal", "!"}, {"rad_vet/clean_stag_xsweep", "!"},
    {"rad_vet/clean_stag_xpat", "!"}, {"rad_vet/clean_thin_donor", "!"},
    {"rad_vet/clean_stag_stencil", "=true"}, {"rad_vet/clean_stag_sbimp", "!"},
    {"rad_vet/clean_lz_bnd_free", "!"}, {"rad_vet/clean_stag_thinf_held", "!"},
    {"rad_vet/clean_face_hclat", "!"},
    {"rad_vet/clean_ray_source_relaxed", "!"},
    // removed, named by decks or restart files: seam clamp, pass print, fixed set's band
    // reach (covered by the bound's band); halo depth cap and ray cut only at kept values
    {"rad_vet/cs_cross_mono", ""}, {"rad_vet/clean_nl_print", ""},
    {"rad_vet/vet_gd_reach", ""}, {"rad_vet/clean_halo_depth", "=0"},
    {"rad_vet/clean_ray_cut", "=true"}};
  const bool say = (global_variable::my_rank == 0);
  std::string none;
  std::vector<std::string> drop;
  for (InputBlock &ib : pin->block) {
    if (ib.block_name != "rad_vet" && ib.block_name != "time") continue;
    for (InputLine &ln : ib.line) {
      const std::string k = ib.block_name + "/" + ln.param_name, &v = ln.param_value;
      for (const auto &t : tab) {
        const std::string p(t[0]);
        const size_t n = p.size() - (p.back() == '*');
        if (k.compare(0, n, p, 0, n) != 0 || (n == p.size() && k.size() != n)) continue;
        if (t[1] == nullptr) break;
        // handled here: never reported as unread (a default-echo line: no message)
        ln.accessed = true;
        std::string c(t[1]);
        if (c[0] == '+') c[0] = '=';
        // (a removed key's value counts also where an older code recorded its default)
        if (ln.param_comment == ParameterInput::default_comment &&
            !(c[0] == '=' && global_variable::restart_run)) break;
        if (c == "x") {
          RadVetFatal("<" + k + "> belonged to a removed path (the cell layout, the FV "
                      "rays, a debug hook removed or a retired study key); remove it "
                      "from the input");
        }
        if (c[0] == '=' && !ParameterInput::SameValue(v, c.substr(1))) {
          RadVetFatal("<" + k + "> = " + v + ": the key was removed, the code runs " +
                      c.substr(1) + " (name " + ln.param_name + "=" + c.substr(1) +
                      " to continue so)");
        }
        if (c == "!" && v != "false" && v != "0" && v != "off") {
          RadVetFatal("<" + k + "> = " + v + " has no counterpart in scheme = clean");
        }
        if (c.empty() || c == "!") none += " " + k;   // (accepted value)
        if (c != "!" && c[0] != '=' && k.compare(0, 17, "rad_vet/implicit_") == 0) {
          drop.push_back(ln.param_name);
        }
        if (say && !c.empty() && c != "!" && c[0] != '=') {
          std::cout << "### WARNING (scheme = clean): <" << k << "> has no effect, the "
                    << "clean solve reads <rad_vet>/" << c << std::endl;
        }
        break;
      }
    }
  }
  if (say && !none.empty()) {
    std::cout << "### WARNING (scheme = clean): no effect:" << none << std::endl;
  }
  if (say && global_variable::restart_run &&
      !pin->DoesParameterExist("rad_vet", "rst_version")) {
    std::string mv;
    for (const auto &t : tab) {
      const std::string k = t[0], b = k.substr(0, k.find('/'));
      if (t[1] == nullptr || t[1][0] != '+' ||
          pin->DoesParameterExist(b, k.substr(b.size() + 1))) continue;
      mv += " " + k.substr(b.size() + 1) + "=" + (t[1] + 1);
    }
    if (!mv.empty()) {
      std::cout << "### NOTE: the restart file is of an older code and holds no value of "
                << "these removed keys; where they apply the code runs" << mv
                << std::endl;
    }
  }
  // implicit_* keys the clean solve does not read leave a fresh run's input, so
  // ImplicitInit parses its defaults; a restart keeps its file's keys (its layout)
  for (InputBlock &ib : pin->block) {
    if (ib.block_name != "rad_vet" || global_variable::restart_run) continue;
    ib.line.remove_if([&](const InputLine &ln) {
      return std::find(drop.begin(), drop.end(), ln.param_name) != drop.end();
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ReadSchemeKeys
//! \brief keys of the implicit scheme

void RadiationVet::ReadSchemeKeys(ParameterInput *pin) {
  // ---- scheme = clean (rad_vet_tasks.cpp), the only scheme of a fresh implicit run
  const bool nw = !global_variable::restart_run;
  // the fresh default of a scheme key, on a restart without it too (fix)
  auto fix = [&](const char *k, const char *v, const char *o) {
    return nw || recipe::Fix(pin, "rad_vet", k, v, o);
  };
  {
    CheckKeys(pin);
    // rst_version, recorded: a restart file without it is of an older code, whose value of
    // a key it does not hold is not the only one left; name it to continue with that
    if (!nw && !pin->DoesParameterExist("rad_vet","rst_version")) {
      auto has = [&](const char *k) {return pin->DoesParameterExist("rad_vet", k);};
      const bool curved = pmy_pack->pmesh->use_spherical_polar ||
                          pmy_pack->pmesh->use_cubed_sphere;
      const bool crv = curved && VetCurvedRays(pin);
      const bool rot = has("vet_gd_rotate_dt") &&
                       (pin->GetString("rad_vet","vet_gd_rotate_dt") == "auto" ||
                        pin->GetReal("rad_vet","vet_gd_rotate_dt") > 0.0);
      const bool be = !has("time_scheme") ||
                      pin->GetString("rad_vet","time_scheme") == "be";
      const std::vector<std::tuple<const char*, const char*, bool>> old = {
          {"scheme", "clean", true}, {"clean_layout", "staggered", true},
          {"clean_rays", "sc", crv}, {"be_legacy", "false", be},
          {"implicit_enthalpy", "plm", true},
          {"vet_gd_rotate_stretch", "0.02", rot},
          {"vet_gd_rotate_slope", "correct", !curved && rot},
          {"vet_rotate_advance", "true", crv && has("vet_rotate_avg") &&
           pin->GetInteger("rad_vet","vet_rotate_avg") > 1},
          {"vet_rotate_avg_ema_wb", "true", crv && has("vet_rotate_avg_ema") &&
           pin->GetReal("rad_vet","vet_rotate_avg_ema") >= 1.0}};
      for (const auto &[k, v, need] : old) {
        if (!need || has(k)) continue;
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "<rad_vet>: the restart file does not hold " << k << " and is of "
          << "an older code, which ran another value; the code runs " << v << " (name "
          << k << "=" << v << " to continue so; inputs/stars/restart_old_files.txt)"
          << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    pin->GetOrAddInteger("rad_vet","rst_version",1);
    // thin-share thinness g = max(|H|/J, eps C1/1.5), eps = min(kappa_P/chi,
    // 1), w = smoothstep((g - c0)/(c1 - c0)): upwind streaming where the field is
    // anisotropic, not merely non-local
    // eps -> 1 - (1 - eps)(1 - exp(-tau_d)), tau_d = chi V/A_f
    cln_c0 = pin->GetOrAddReal("rad_vet","clean_w_c0",0.3);
    cln_c1 = pin->GetOrAddReal("rad_vet","clean_w_c1",0.6);
    cln_tol = pin->GetOrAddReal("rad_vet","clean_tol",1.0e-10);
    // clean_maxnl, clean_maxit, clean_conv_max do not take implicit_maxit,
    // implicit_lin_maxit, implicit_tol_rel_max (a named other value prints a NOTE).
    // A linear solve ending at clean_maxit leaves its iterate to the next pass.
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
    recipe::Old(pin, "rad_vet", "clean_maxit", "800");
    cln_maxit = pin->GetOrAddInteger("rad_vet","clean_maxit",nw ? 200 : 800);
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
    }
    cln_tol_nl = pin->GetOrAddReal("rad_vet","clean_tol_nl",1.0e-11);
    cln_raymin = pin->GetOrAddReal("rad_vet","clean_ray_min",1.0e-3);
    // clean_rot_ring = off | live | jump (live where the vet_gd rays rotate with a
    // cross-fade, InitArrays): the rotating union swept as the current (+ next) set,
    // older sets as moment offsets (RingMoments). clean_rot_fuse: the live ring's two
    // sets in one fused sweep (bitwise the two sweeps; a second intensity array)
    cln_rfuse = pin->GetOrAddBoolean("rad_vet","clean_rot_fuse",
                                     fix("clean_rot_fuse", "true", "false"));
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
    // clean_stag_opac_passes: opacities of passes 0 .. N-1, then frozen
    cln_optn = pin->GetOrAddInteger("rad_vet","clean_stag_opac_passes",3);
    // fschur where it applies: F lines exact across MeshBlocks
    // and periodic wraps (FluxSchur); needs a uniform mesh and MeshBlocks of >= 3
    // cells along a line (>= 6 on longer or open lines)
    {Mesh *pm = pmy_pack->pmesh;
    // (cubed sphere: a lateral line is a ring of four panels' MeshBlocks)
    const bool cs = pm->use_cubed_sphere;
    bool ok = !pm->multilevel, lag = false;
    for (int d = 0; d < 3; ++d) {
      // (layout_invariant_mode = tile: tiles as MeshBlocks, FluxSchurTile)
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
    if (!ok && lag && !cs) {   // no silent lagged ends
      std::cout << "### FATAL ERROR: <rad_vet>: the F solve's lines need nghost >= 3 and "
                << "MeshBlocks of >= 3 cells along a line of two round a periodic one, "
                << ">= 6 along longer or open ones" << std::endl;
      std::exit(EXIT_FAILURE);
    } else if (!ok && global_variable::my_rank == 0) {
      std::cout << "<rad_vet> the F solve's line ends are lagged (mesh or MeshBlock size)"
                << std::endl;
    }
    // layout_tile: a periodic line joins whole MeshBlocks (FluxSchurTile)
    for (int d = 1; d < 3 && cln_fschur && repro::tile_n > 0 && !cs; ++d) {
      const int nx = (d == 1) ? pm->mb_indcs.nx2 : pm->mb_indcs.nx3;
      const int nbm = ((d == 1) ? pm->mesh_indcs.nx2 : pm->mesh_indcs.nx3)/nx;
      if (nx > repro::tile_n && nbm%2 == 1 && nbm > 1 &&
          pm->mesh_bcs[2*d] == BoundaryFlag::periodic) {
        std::cout << "### FATAL ERROR: <mesh>/layout_tile: a "
                  << "periodic line of MeshBlocks of more than one tile needs an even "
                  << "number of them (or one)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    // cubed sphere without fschur iterates its seams (Thomas lines)
    if (cs && !cln_fschur) cln_fsweeps = 64;
    recipe::Rule(pin, "rad_vet", "clean_stag_fschur", cln_fschur ? "true" : "false");
    recipe::Rule(pin, "rad_vet", "clean_stag_fsweeps", std::to_string(cln_fsweeps));}
    // clean_consistency: default on Cartesian, and on sp / cs with the vet_gd rays (sc)
    const bool curv = pmy_pack->pmesh->use_spherical_polar ||
                      pmy_pack->pmesh->use_cubed_sphere;
    const bool vg = VetCurvedRays(pin);
    cln_split = true;
    cln_consist = pin->GetOrAddBoolean("rad_vet","clean_consistency", !curv || vg);
    // staggered: energy flux apart from the face momentum (split), face momentum of its
    // own F solve (fsolve) in the Newton loop (fwork); with consistency the rays'
    // kinetic momentum flux (fkin) and face flux as flux factors (gflux); donor-cell
    // advective share (upwind). Read only when named.
    auto stag = [&](const char *key, const bool def) {
      return pin->DoesParameterExist("rad_vet", key) ? pin->GetBoolean("rad_vet", key) :
                                                      def;
    };
    cln_fsolve = cln_upw = cln_split;
    cln_fkin = cln_fwork = cln_gflux = cln_split && cln_consist;
    for (const char *k : {"clean_stag_fwb", "clean_stag_freal"}) {
      recipe::Rule(pin, "rad_vet", k, cln_consist ? "true" : "false");
    }
    // freal: realisable sb (comoving part <= 2 |beta| hc kE E_k) and restart
    // faces <= c E of their donor
    cln_freal = cln_gflux;
    // clean_stag_thinf: the bracket's ray correction in the face half ranges, realisable
    // (FaceBracket)
    cln_thinf = pin->GetOrAddBoolean("rad_vet","clean_stag_thinf",
                                     fix("clean_stag_thinf", "true", "false"));
    // fwb: F rows in deviation form about the rays' own state, their exact steady state
    // (FluxSolve); fwb_ring: its s of the phi faces (3-D sp, periodic x3) less its
    // ring mean: no torque on a ring (FluxBalanceProject)
    cln_fwb = cln_consist;
    cln_fwbring = cln_fwb && pmy_pack->pmesh->use_spherical_polar &&
        pmy_pack->pmesh->three_d &&
        pmy_pack->pmesh->mesh_bcs[BoundaryFace::inner_x3] == BoundaryFlag::periodic;
    recipe::Rule(pin, "rad_vet", "clean_stag_fwb_ring", cln_fwbring ? "true" : "false");
    // clean_stag_fcone (default on with fwb, recorded): every cell's flux vector within
    // its realisable cone, |F| <= c E, by construction: the faces scaled towards 0
    // (FluxCone)
    cln_fcone = pin->GetOrAddBoolean("rad_vet","clean_stag_fcone",cln_fwb &&
                                     fix("clean_stag_fcone", "true", "false"));
    if (cln_fcone && !cln_fsolve) {
      std::cout << "### FATAL ERROR: <rad_vet>/clean_stag_fcone needs the staggered F "
                << "solve" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_stag_fcone_lz (default on, not on Cartesian meshes): fcone takes a cell's
    // excess from its faces without L_z lever arm first (LzFamily: spherical polar
    // x1, x2; cubed sphere x1; Cartesian x3), cut and counted as now; a face with an arm
    // is scaled only where those faces alone exceed c E, and its cut is the gas's at the
    // face (as fcone_gas; clean_lz_spin: a rigid rotation), so that L_z is exact (cubed
    // sphere: CsTorque)
    // clean_stag_fcone_gas (default false): the gas takes the momentum fcone cuts at each
    // inner face (the cut F's absorption is not withheld), so gas + radiation momentum
    // and L_z are exact; E pays the kinetic energy (the write-back's residual), which
    // a radiation-poor cell cannot (test floor_momentum's exch0 box blows up)
    cln_fcgas = cln_fcone && pin->GetOrAddBoolean("rad_vet","clean_stag_fcone_gas",false);
    cln_fclz = cln_fcone && pin->GetOrAddBoolean("rad_vet","clean_stag_fcone_lz", curv &&
                                fix("clean_stag_fcone_lz", "true", "false"));
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
    // clean_cs_lz_abs = 2 (cubed sphere, the only form): the gas phi momentum is the
    // faces' absorption as its r and theta; the L_z it does not take is booked as the
    // radiation's (LzBook)
    if (pmy_pack->pmesh->use_cubed_sphere) cln_lza = 2;
    cln_lzbook = (cln_lza == 2) || (pmy_pack->pmesh->use_spherical_polar &&
                                    pmy_pack->pmesh->three_d);
    // clean_cs_lz_tel: L_z by telescoping (gas phi momentum of CsTorque, radiation
    // L_z of its face momenta), so the ledger holds the radial boundaries' torque alone
    if (cln_lza == 2 || pin->DoesParameterExist("rad_vet", "clean_cs_lz_tel")) {
      cln_lzt = pin->GetOrAddBoolean("rad_vet","clean_cs_lz_tel",
                                     fix("clean_cs_lz_tel", "true", "false"));
    }
    // clean_cs_lat_wb: a spherically symmetric P gives the cs lateral faces no force
    // (CsLat).  Default on for the cubed sphere, also on a restart without it
    if (cln_lza == 2) {
      cln_latwb = pin->GetOrAddBoolean("rad_vet","clean_cs_lat_wb",
                                       fix("clean_cs_lat_wb", "true", "false"));
      recipe::OldStored(pin, "rad_vet", "clean_cs_lat_wb", "false", "true");
    } else if (stag("clean_cs_lat_wb", false)) {
      std::cout << "### FATAL ERROR: <rad_vet>/clean_cs_lat_wb needs the cubed sphere "
                << "and clean_layout = staggered" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (cln_lzt && cln_lza != 2) {
      std::cout << "### FATAL ERROR: <rad_vet>/clean_cs_lz_tel needs the cubed sphere "
                << "and clean_layout = staggered" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_lz_spin (default on with clean_stag_fcone_lz, not on Cartesian meshes; needs
    // clean_stag_fcone_lz, on the cubed sphere clean_cs_lz_tel): the L_z of the momentum
    // that fcone's cut of the faces with an arm and clean_wb_mcap's cap withhold is the
    // gas's as one rigid rotation dOmega = X/I, I = sum rho V R^2 (WriteBack), not the
    // face's gas: L_z exact, no kick
    const bool lzs_ok = cln_fclz && (!pmy_pack->pmesh->use_cubed_sphere || cln_lzt);
    cln_lzsp = pin->GetOrAddBoolean("rad_vet","clean_lz_spin", curv && lzs_ok &&
                                    fix("clean_lz_spin", "true", "false"));
    if (cln_lzsp && !lzs_ok) {
      std::cout << "### FATAL ERROR: <rad_vet>/clean_lz_spin needs clean_stag_fcone_lz "
                << "(the cubed sphere: and clean_cs_lz_tel)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (cln_lzbook) {
      const bool cs = (cln_lza == 2);
      cln_lzpin = pin;
      for (int k1 = 0; k1 < 2; ++k1) {   // a restart file's non-finite value: NOTE, 0
        const std::string key = VetLzKey(pin, cs, k1 == 1);
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
    // clean_clos1_budget (default; false: the strict norm): a pass whose iterate the next
    // pass takes the closure of in the budget's norm, to 0.1 clean_conv_max, the floor of
    // the passes after it at the budget's ceiling (a step dt_retry redoes strictly, as
    // its budget)
    cln_c1b = pin->GetOrAddBoolean("rad_vet","clean_clos1_budget",
                                   fix("clean_clos1_budget", "true", "false"));
    // clean_bud_maxnorm (default; false: the L2 norm): the budget passes' linear test and
    // Eisenstat-Walker |r0| in the max norm over the cells, the norm of the pass's own
    // test (NewtonUpdate)
    cln_bmn = pin->GetOrAddBoolean("rad_vet","clean_bud_maxnorm",
                                   fix("clean_bud_maxnorm", "true", "false"));
    // clean_wb_exch_local (default): the write-back's cell exchange
    // solved exactly from the transported E (E > 0, gas >= e(rho, tfloor))
    cln_wbxl = pin->GetOrAddBoolean("rad_vet","clean_wb_exch_local",
                                   fix("clean_wb_exch_local", "true", "false"));
    // clean_wb_mcap (read only when named; default false): the write-back gives the gas
    // the face momentum terms at the Newton passes' cap h chi_f E, the momentum whose
    // work the E row paid (FaceMomentumChange); the terms withheld are printed
    if (pin->DoesParameterExist("rad_vet","clean_wb_mcap")) {
      cln_wbmc = pin->GetBoolean("rad_vet","clean_wb_mcap");
    }
    // clean_wb_floor_diag N (diagnostic, read only when named): the terms of the first
    // N cells of each write-back whose E the floor raises (CLEAN WBFDIAG)
    if (pin->DoesParameterExist("rad_vet","clean_wb_floor_diag")) {
      cln_wbfdg = pin->GetInteger("rad_vet","clean_wb_floor_diag");
    }
    // clean_dbg_inf C (test, read only when named): the first stage-1 write-back of cycle
    // C gives rank 0's first active cell a non-finite E iterate (the test of the
    // non-finite rejection, tests/vet_clean inf_reject)
    if (pin->DoesParameterExist("rad_vet","clean_dbg_inf")) {
      cln_dbginf = pin->GetInteger("rad_vet","clean_dbg_inf");
    }
    // preconditioner: mg (VCycle)
    {
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
    // clean_stag_mg_levels N (default 2): the mg's coarse levels at
    // most, each halving the tile's x2 and x3 where even (>= 4); 0: the sweep alone
    cln_mgnl = pin->GetOrAddInteger("rad_vet","clean_stag_mg_levels",cln_mgnl);
    if (cln_mgnl < 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_stag_mg_levels N >= 0" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // clean_pc_sweep = sym (the default, recorded) | rb | rbr | rbrsym: the mg's
    // smoothing on each level, a forward red-black sweep, + a second red pass (rbr), + a
    // black-red post-sweep after the coarse correction (sym) (VCycle)
    const std::string sw = pin->GetOrAddString("rad_vet","clean_pc_sweep",
        fix("clean_pc_sweep", "sym", "rb") ? "sym" : "rb");
    cln_pcsweep = (sw == "rbr") + 2*(sw == "sym") + 3*(sw == "rbrsym");
    if (!cln_pcsweep && sw != "rb") {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/clean_pc_sweep = " << sw << ": rb | rbr | sym "
                << "| rbrsym (wave was removed)" << std::endl;
      std::exit(EXIT_FAILURE);
    }}
    // sp / Cartesian: the faces' lateral stress from the cells' sums Y
    // (CsLateralStress), the cs form
    cln_ysum = pmy_pack->pmesh->multi_d && !pmy_pack->pmesh->use_cubed_sphere;
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
    // clean_verbose (read only when named): 0 (default) the CLEAN cycle lines, warnings
    // and events; 1 also the per-solve diagnostics (CLEAN MCAP, OLDVEC, FCONE, FWB, LZ,
    // WBFCT, KPOS, KREAL, XPAT, the write-back closure and the ray cut); 2 also the
    // CLEAN PASS line of every pass (its changes de, dv against tol_nl)
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
    // clean_init = given | steady: steady relaxes the radiation
    // before step 1
    {
      const std::string ci = pin->GetOrAddString("rad_vet","clean_init","given");
      if (ci != "given" && ci != "steady") {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "<rad_vet>/clean_init = given | steady" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      cln_steady = (ci == "steady");
    }
    if (cln_steady) {
      // problem/star_ic_equil: its relax keys give the defaults
      const bool heq = pin->DoesParameterExist("problem","star_ic_equil") &&
                       pin->GetBoolean("problem","star_ic_equil");
      auto eqk = [&](const char *k, const Real d) {
        return (heq && pin->DoesParameterExist("problem", k)) ?
               pin->GetReal("problem", k) : d;
      };
      cln_steady_n = pin->GetOrAddInteger("rad_vet","clean_init_n",
          static_cast<int>(eqk("star_ic_equil_relax_n", heq ? 40.0 : 50.0)));
      cln_steady_dtfac = pin->GetOrAddReal("rad_vet","clean_init_dtfac",
                                           eqk("star_ic_equil_dtfac", 1.0e6));
      cln_steady_tol = pin->GetOrAddReal("rad_vet","clean_init_tol",1.0e-8);
    }
  }
}

}  // namespace radvet
