# hegiant-1013: every key changed vs the 10-09 He giant deck (hegiant_fresh897_raven_eft.athinput 9c89872a) -> FIN3 deck
# (Viper deck; Raven adds vet_gd_twin_lowmem = true; Caltech = Viper with @FILES@ paths). Paths excluded.

## CHANGED (7)
| block/key | old | new | note |
|---|---|---|---|
| hydro/eos_logt_min | 3.5 | 3.4 | hegiant-1013: covers tfloor 3e3 K (as E1); was 3.5 |
| hydro/tfloor_kelvin | 4.0e3 | 3.0e3 | hegiant-1013: E1 rule (user 10-08): below T_rad at x1max (IC T_RE(200 Rsun) ~3.5e3 K); was 4.0e3 |
| problem/he_sponge_r0 | 1.384443e13 | 1.25226e13 | hegiant-1013: E1 rule r0 = 0.9 x1max (180 Rsun); was 199 Rsun (0.995 x1max) |
| rad_m1/implicit_blend_xthin | 30 | 0 | final_1012 |
| rad_m1/implicit_opac_newton | false | true |  |
| rad_m1/implicit_vimp | false | true | final_1012 |
| rad_m1/time_scheme | be | hesdirk2 | final_1012 |

## ADDED (named; 181) -- FIN3 <rad_m1> keys from AG Car E1 final + E1 IC/diagnostic keys + binary defaults named
| block/key | value | note |
|---|---|---|
| hydro/fofc_report | 512 | hegiant-1013: diagnostic only (as E1; bg check counts) |
| hydro/fofc_report_dt | 3.6e3 | hegiant-1013: diagnostic only (as E1) |
| problem/he_bc_inner | wall | hegiant-1013: named (default, as E1) |
| problem/he_esrc_const | 0 | hegiant-1013: named at the binary default (E1 PAR_DUMP list) |
| problem/he_gravity | wb | hegiant-1013: named (default; he_ic_equil needs it) |
| problem/he_ic_equil | false | hegiant-1013: DI IC OFF: bg check F (equil on) re-marched atmosphere falls back, 1.03-1.10 R_ph drho 0.41 vs 0 |
| problem/he_ic_equil_dtfac | 1.0 | hegiant-1013: as E1 |
| problem/he_ic_equil_last_march | true | hegiant-1013: as E1 (v6 f1La) |
| problem/he_ic_equil_maxit | 1 | hegiant-1013: as E1 |
| problem/he_ic_equil_relax_n | 40 | hegiant-1013: as E1 |
| problem/he_ic_equil_relax_tol | 1.0e-4 | hegiant-1013: as E1 |
| problem/he_ic_equil_tol | 1.0e-9 | hegiant-1013: as E1 |
| problem/he_ic_rad | file | hegiant-1013: NOT E1's sph1d: smoke 12179235 FATAL (sph1d HSE loop: Gamma>0.9 clamped, scaffold carries 78-99  |
| problem/he_ic_rad_atm_post | false | hegiant-1013: superseded by he_ic_equil (exclusive in code) |
| problem/he_seed_k | 4 | hegiant-1013: named at the binary default (E1 PAR_DUMP list) |
| problem/he_seed_kdist | box | hegiant-1013: named at the binary default (E1 PAR_DUMP list) |
| problem/he_seed_rng | 1234 | hegiant-1013: named at the binary default (E1 PAR_DUMP list) |
| problem/he_wind_ic | 0 | hegiant-1013: named at the binary default (E1 PAR_DUMP list) |
| problem/hot_jupiter | 0 | hegiant-1013: named at the binary default (E1 PAR_DUMP list) |
| problem/mlt_closure | frozen | hegiant-1013: named as E1 (= mlt_flux_frozen) |
| problem/mlt_closure_taper | 0.01 | hegiant-1013: named at the binary default (E1 PAR_DUMP list) |
| problem/nan_report | 0 | hegiant-1013: named at the binary default (E1 PAR_DUMP list) |
| rad_m1/advect_split | 0 |  |
| rad_m1/ap_form | alpha2 |  |
| rad_m1/cfl_rad | 0.4 |  |
| rad_m1/dbg_energy_tally | 0 |  |
| rad_m1/dbg_gas_force | 1 |  |
| rad_m1/dbg_gas_force_trans | 1 |  |
| rad_m1/dbg_gas_heat | 1 |  |
| rad_m1/dbg_opac_part | 0 |  |
| rad_m1/dbg_tensor | none |  |
| rad_m1/dbg_tensor_tilt | 0.3 |  |
| rad_m1/dbg_trans_memory | 1 |  |
| rad_m1/f_source | cell |  |
| rad_m1/ic_rad | fs_relax | prodfinal-1011 |
| rad_m1/ic_rad_fs_iter | 3 |  |
| rad_m1/ic_rad_relax_dtfac | 1 |  |
| rad_m1/ic_rad_relax_n | 40 |  |
| rad_m1/ic_rad_relax_tol | 0.0001 |  |
| rad_m1/implicit_accel | none |  |
| rad_m1/implicit_allow_multid | 0 |  |
| rad_m1/implicit_anderson_beta | 1 |  |
| rad_m1/implicit_anderson_m | 5 |  |
| rad_m1/implicit_anderson_start | 1 |  |
| rad_m1/implicit_bc_advect | 1 |  |
| rad_m1/implicit_bcg_fallback | line |  |
| rad_m1/implicit_bcg_max_restarts | 2 |  |
| rad_m1/implicit_bcg_sync | 1 |  |
| rad_m1/implicit_blend_fmode | max |  |
| rad_m1/implicit_blend_mode | flux |  |
| rad_m1/implicit_blend_r0 | 1.5 |  |
| rad_m1/implicit_blend_xthin_mode | all |  |
| rad_m1/implicit_blend_xthin_wmin | 0 |  |
| rad_m1/implicit_bmom_half | 1 |  |
| rad_m1/implicit_cfl | -1 |  |
| rad_m1/implicit_closure_relax | 1 |  |
| rad_m1/implicit_closure_relax_thin | 0 |  |
| rad_m1/implicit_conv_est | 1 |  |
| rad_m1/implicit_cs2_hll | drag | final_1012 |
| rad_m1/implicit_cs2_shift | false | final_1012 |
| rad_m1/implicit_ebath_x1max | 0 |  |
| rad_m1/implicit_ebath_x1min | 0 |  |
| rad_m1/implicit_eos_cache | 0 |  |
| rad_m1/implicit_eos_cache_check | 1 |  |
| rad_m1/implicit_eos_cache_nt | 2 |  |
| rad_m1/implicit_f0_init | auto |  |
| rad_m1/implicit_face_weight | distance |  |
| rad_m1/implicit_fast_kernels | 1 |  |
| rad_m1/implicit_gas_newton_switch_min | 20 |  |
| rad_m1/implicit_gas_src | conservative2 | final_1012 |
| rad_m1/implicit_halo_ipc | 0 |  |
| rad_m1/implicit_halo_ovl_faces | 0 |  |
| rad_m1/implicit_hr_excess | none | prodfinal-1011 |
| rad_m1/implicit_hr_face | sf |  |
| rad_m1/implicit_hr_face_gd | none | FIN3 (prodswitch-1013, fin3-1012) |
| rad_m1/implicit_hr_model | 0 |  |
| rad_m1/implicit_krylov_pipe | 0 |  |
| rad_m1/implicit_lin_ew_gamma | 0.9 |  |
| rad_m1/implicit_line_solver | pcr |  |
| rad_m1/implicit_lres_test | 0 |  |
| rad_m1/implicit_marshak_face | linear |  |
| rad_m1/implicit_mg_halo | 1 |  |
| rad_m1/implicit_mg_levels | 3 |  |
| rad_m1/implicit_nofb_work | 1 |  |
| rad_m1/implicit_op_team_red | 0 |  |
| rad_m1/implicit_opac_newton_guard | 0.5 |  |
| rad_m1/implicit_opac_newton_slope_max | 3 | opacity-cliff damping (H/He recombination; He giant fix 90a1d462), user 10-08 |
| rad_m1/implicit_partition | none |  |
| rad_m1/implicit_pcr_check | 0 |  |
| rad_m1/implicit_pcr_team | 0 |  |
| rad_m1/implicit_picard_log | 0 |  |
| rad_m1/implicit_precond_float | 0 |  |
| rad_m1/implicit_predictor_frozen | false | FIN3: UN #10 dropped (frozen predictor costs Cartesian T-ladder order; inert on sp). coordinator 2026-10-13, r |
| rad_m1/implicit_predictor_order | 2 |  |
| rad_m1/implicit_rad_inertia | 1 |  |
| rad_m1/implicit_recon | dc |  |
| rad_m1/implicit_recon_lag | step |  |
| rad_m1/implicit_recon_npass | -1 |  |
| rad_m1/implicit_recon_w | -1 |  |
| rad_m1/implicit_res_floor | 0 |  |
| rad_m1/implicit_src_stable | 1 |  |
| rad_m1/implicit_stall_fac | 0.5 |  |
| rad_m1/implicit_thin_bq | lin | final_1012 |
| rad_m1/implicit_thin_bq_gd | cell | final_1012 |
| rad_m1/implicit_thin_corr | hc | final_1012 |
| rad_m1/implicit_thin_fsa | 0 |  |
| rad_m1/implicit_thin_fsj | 0 |  |
| rad_m1/implicit_thin_h0 | 0.3 |  |
| rad_m1/implicit_thin_h1 | 0.6 |  |
| rad_m1/implicit_thin_hcface | tau | final_1012 |
| rad_m1/implicit_thin_hcface_gd | x1 | final_1012 |
| rad_m1/implicit_thin_hcsf | j | final_1012 |
| rad_m1/implicit_thin_js_face | hcface | jsfix 8bfeae8c, d3bf42d7 final method (user 2026-10-12) |
| rad_m1/implicit_thin_js_gd | x1 | final_1012 |
| rad_m1/implicit_thin_kn0 | 1 |  |
| rad_m1/implicit_thin_knp | 4 |  |
| rad_m1/implicit_thin_recon | jse | final_1012 |
| rad_m1/implicit_thin_switch | fse | final_1012 |
| rad_m1/implicit_trans_fmax | 1 |  |
| rad_m1/implicit_trans_limit | none |  |
| rad_m1/implicit_tsolve_opac_start | 20 |  |
| rad_m1/implicit_vimp_bc | true | final_1012 |
| rad_m1/implicit_vimp_fold | 0 |  |
| rad_m1/kappa_e | 0 |  |
| rad_m1/kappa_f | 0 |  |
| rad_m1/kappa_p | 0 |  |
| rad_m1/kappa_s | 0 |  |
| rad_m1/opac_a | 0 |  |
| rad_m1/opac_b | 0 |  |
| rad_m1/opac_freeze | 0 |  |
| rad_m1/reconstruct | plm |  |
| rad_m1/rho_ref | 1 |  |
| rad_m1/rsla_force | 0 |  |
| rad_m1/rsla_vmax | -1 |  |
| rad_m1/rsla_warn | 0.1 |  |
| rad_m1/scaled_prefactor | 20 |  |
| rad_m1/source_form | full |  |
| rad_m1/split_vel | recon |  |
| rad_m1/subcycle | 1 |  |
| rad_m1/t_ref | 1 |  |
| rad_m1/thick_flux | none |  |
| rad_m1/time2_lin_tol_fac | 1 |  |
| rad_m1/time2_one_pass_safety | 30 |  |
| rad_m1/time2_pred_cap | 0.5 | final_1012 |
| rad_m1/time2_start | be | final_1012 |
| rad_m1/time2_vet_col_source | blend | final_1012 |
| rad_m1/time2_vet_sc | predict | final_1012 |
| rad_m1/time2_vstage | true | final_1012 |
| rad_m1/vet_bc_diff | face | final_1012 |
| rad_m1/vet_bc_inner | auto | FIN3 (prodswitch-1013, fin3-1012) |
| rad_m1/vet_bc_reflect | auto | FIN3 (prodswitch-1013, fin3-1012) |
| rad_m1/vet_col_axis | radial |  |
| rad_m1/vet_col_chunk | 0 |  |
| rad_m1/vet_col_fk_min | 0 |  |
| rad_m1/vet_col_lat_cut_thin | 1.0 | FIN3 (prodswitch-1013, fin3-1012) |
| rad_m1/vet_col_ncore | 8 |  |
| rad_m1/vet_col_nmu | 4 |  |
| rad_m1/vet_col_nsub | 1 |  |
| rad_m1/vet_col_order2 | 1 |  |
| rad_m1/vet_col_reflect_top | 1 |  |
| rad_m1/vet_col_source | relaxed |  |
| rad_m1/vet_col_surface_face | 1 |  |
| rad_m1/vet_col_surface_q | 1 |  |
| rad_m1/vet_col_surface_qmax | 1 |  |
| rad_m1/vet_col_surface_qmin | 0.001 |  |
| rad_m1/vet_col_team | 1 |  |
| rad_m1/vet_col_team_size | 0 |  |
| rad_m1/vet_gd_halo_ipc | false | hegiant-1013: layout key: CUDA IPC halo is single-node CUDA only (Raven runs 2 nodes; HIP no) |
| rad_m1/vet_gd_rotate_align | true | final_1012 |
| rad_m1/vet_gd_rotate_dt | 4.0e2 | hegiant-1013: ~20 dt at dt ~20 s (BSG 1.75e3 = 20 x 87.5 s, E1 1.6e4 = 22 x 727 s) |
| rad_m1/vet_gd_rotate_slope | correct | d3bf42d7 final method (user 2026-10-12) |
| rad_m1/vet_gd_rotate_stretch | 0.02 | one method: = BSG and fresh default (user 2026-10-12) |
| rad_m1/vet_gd_rotate_t0 | 0.0 | final_1012 |
| rad_m1/vet_gd_symref | col | final_1012 |
| rad_m1/vet_gd_twin_state | trans | FIN3 (prodswitch-1013, fin3-1012) |
| rad_m1/vet_gd_twin_weight | none | FIN3 (prodswitch-1013, fin3-1012) |
| rad_m1/vet_interp | auto | FIN3 (prodswitch-1013, fin3-1012) |
| rad_m1/vet_kappa_es | 0.34 |  |
| rad_m1/vet_kappa_floor | 1e-05 |  |
| rad_m1/vet_rotate_avg | 4 | FIN3 (prodswitch-1013, fin3-1012) |
| rad_m1/vet_seg | jump | final_1012 |

## REMOVED (0)
