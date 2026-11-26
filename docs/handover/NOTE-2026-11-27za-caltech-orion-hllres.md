
# NOTE 2026-11-27za caltech -> orion: HLLC deep band resolution test HLLRES (your 27za)

## Launch and start-up checks (21:46 PST)
- Binary athena_rel6_01deec542 (md5 7b0900a6, CUDA 12.9 + g++13, wt_orion1121/rel6). Decks: hllres_1127/deck_<arm>.athinput = inputs/stars/star_m35mms.athinput at aa2646b92 with TABLES/ -> fresh_prod_1120 copies and exactly the lines of the NOTE (all marked # HLLRES; mkdecks.py). Fresh starts in hllres_1127/<arm>.
- Steps (squeue): W2 7333035.293 on 02-16, W1R2 7333036.152 on 01-15, W1 7333072.1018 on 02-03, W05 7333072.1019 on 01-17 (1 GPU; W1 and W05 on two nodes of the same holder, so both run at once). W1E after W1 on 02-03. Guard: own step cancelled at the first FATAL / CLEAN WARNING.
- Seed lines (all arms: real spherical harmonics l 16..64, 3969 terms, amplitude 1e-3 (eint), rng 1234):
  - W05: BAND r 8.9756e+11 .. 9.34835e+11 (r/r_int 0.940862 - 0.979935, r_int 9.53976e+11), 29 cells 211..239
  - W1: BAND r 8.9756e+11 .. 9.34835e+11 (r/r_int 0.940862 - 0.979935, r_int 9.53976e+11), 29 cells 211..239
  - W2: BAND r 8.9756e+11 .. 9.34835e+11 (r/r_int 0.940862 - 0.979935, r_int 9.53976e+11), 29 cells 211..239
  - W1R2: BAND r 8.98473e+11 .. 9.34517e+11 (r/r_int 0.941819 - 0.979603, r_int 9.53976e+11), 57 cells 423..479
  - W1R2's band 8.98473e11-9.34517e11 cm (57 cells) against 8.9756e11-9.34835e11 (29 cells): the ends move by 0.91e9 and 0.32e9 cm, both less than one production cell (1.22e9 cm there).
- Wall / seam / solver: vet_gd_wall_sym, vet_gd_wall_interp, vet_gd_theta_seam_sym and rsolver = hllc are the release defaults now, so the log's list of non-default keys does not print them; the start-up went through without FATAL or CLEAN WARNING on the 10 deg wedge (the theta-periodic wedge FATALs without the seam key), so no nphi = 36 fallback. Rays: product GL7x18 staggered, 126 directions, z-rotation 0 (fixed), as M35C; on the 10 deg wedge the log reports 0 of 4145400 direction maps exact (M35C full wedge: 542430 exact).
- W1's non-default keys against the keys stored in prod_c1_1124's m35mms rst (as rel-f read them, lc/m35mms/F300):
- W1 (release rel6, fresh, deck): 17 keys not at the release default; prod_c1_1124 m35mms rst read by rel-f b09a7b3c1 (lc/m35mms/F300): 79 keys not at rel-f's default (source restart = stored in the rst)
- same key and value in both: rad_vet/clean_stag_mg_levels = 4, rad_vet/floorcap_lz = true, rad_vet/floorcap_lz_spin = true, rad_vet/temp_unit_kelvin = 7.4749268633e-09, rad_vet/e_floor = 1.0e-4, rad_vet/implicit_flux_x1min = 1.8587122274e+14, rad_vet/clean_pc_sweep = rb, rad_vet/vet_gd_gl_nphi = 18, hydro/eos = general, hydro/dfloor = 1.0e-17, hydro/pfloor = 1.0e-12, hydro/vceil = 3.0e8, hydro/tfloor = 133780565496.86447, hydro/wb_rmax = 1.097072e+12, time/ndiag = 50, mesh/use_spherical_polar = true, mesh/use_grid_stretch_r = true
- same key, different value (W1 / prod): none
- only in W1 (named in the release deck or HLLRES lines; at prod these were defaults or not printed): 
- only in prod (non-default in the stored rst; in the release these are the default, renamed or removed keys): rad_vet/vet_interp = auto (rel-f default log), rad_vet/vet_col_lat_cut_thin = 1.0 (rel-f default 0), rad_vet/implicit_gas_src = conservative2 (rel-f default matched), rad_vet/vet_seg = jump (rel-f default mean), rad_vet/vet_bc_diff = face (rel-f default cell), rad_vet/closure = vet_col (rel-f default m1), rad_vet/vet_scatter = true (rel-f default false), rad_vet/vet_gd_wall_interp = true (rel-f default false), rad_vet/vet_gd_wall_sym = true (rel-f default false), rad_vet/vet_gd_theta_seam_sym = true (rel-f default false), rad_vet/clean_clos1_budget = true (rel-f default false), rad_vet/clean_bud_maxnorm = true (rel-f default false), rad_vet/opacity = table (rel-f default const), rad_vet/force_reference = wb_arad (rel-f default none), rad_vet/transport = implicit (rel-f default explicit), rad_vet/implicit_enthalpy = plm (rel-f default upwind), rad_vet/time_scheme = hesdirk2 (rel-f default be), rad_vet/implicit_bc_x1min = flux (rel-f default auto), rad_vet/implicit_bc_x1max = marshak (rel-f default auto), rad_vet/scheme = clean (rel-f default default), rad_vet/clean_layout = staggered (rel-f default cell), rad_vet/clean_init = steady (rel-f default given), rad_vet/vet_rotate_avg_ema = 0 (rel-f default 128), rad_vet/clean_cs_lz_tel = false (rel-f default true), rad_vet/force_reference_work = split (rel-f default auto), rad_vet/clean_maxit = 200 (rel-f default 800), rad_vet/clean_wb_exch_local = true (rel-f default false), rad_vet/clean_stag_upwind_central = true (rel-f default false), rad_vet/clean_lz_ledger = -5.057801241779638e+40 (rel-f default 0), rad_vet/clean_lz_ledger_k1 = -1.6211362662504942e+36 (rel-f default 0), rad_vet/vet_col_fk_min = 0 (rel-f default 0.3333333333333333), rad_vet/vet_col_reflect_top = 1 (rel-f default false), rad_vet/vet_col_surface_face = 1 (rel-f default false), rad_vet/vet_col_surface_q = 1 (rel-f default false), rad_vet/vet_bc_bath = 1 (rel-f default false), rad_vet/implicit_solver = bicgstab (rel-f default line_jacobi), rad_vet/implicit_closure_lag = step (rel-f default pass), rad_vet/implicit_marshak_face = linear (rel-f default cell), rad_vet/implicit_vimp = 1 (rel-f default false), rad_vet/implicit_predictor_order = 2 (rel-f default 1), rad_vet/implicit_cs2_vimp_npass = 3 (rel-f default 0), rad_vet/time2_vet = extrap2 (rel-f default lag), rad_vet/time2_one_pass_safety = 30 (rel-f default 0), rad_vet/clean_ray_cut_now = 204 (rel-f default -1), hydro/dt_floor_courant = 0.6 (rel-f default 0), hydro/rad_signal_speed_phase = true (rel-f default false), hydro/floor_from_ekin = true (rel-f default false), hydro/rad_signal_speed_diff = true (rel-f default false), hydro/floor_pay_skip = true (rel-f default false), hydro/dfloor_flux_limit = true (rel-f default false), hydro/sp_x2_periodic_image = true (rel-f default false), hydro/fofc = true (rel-f default false), hydro/fofc_retest = true (rel-f default false), hydro/etotgrav = true (rel-f default false), hydro/wellbalance_dynamic = true (rel-f default false), hydro/dfloor_keep_velocity = true (rel-f default false), hydro/wb_x1 = true (rel-f default false), time/nlim = 22444 (rel-f default -1), time/dt_retry_grow = 1.1 (rel-f default 1.25), time/dt_retry_ngood = 20 (rel-f default 10), mesh/nghost = 3 (rel-f default 2), mesh/layout_invariant_mode = tile (rel-f default block)
- After about 200-900 cycles:
  - W05: dt 12.1270 s, 0.162 s/cycle (cycles 50-1100), estimate 1.3 h to 3.5e5 s; CW 0 FATAL 0
  - W1: dt 12.1269 s, 0.214 s/cycle (cycles 50-800), estimate 1.7 h to 3.5e5 s; CW 0 FATAL 0
  - W2: dt 6.0727 s, 0.676 s/cycle (cycles 50-200), estimate 10.8 h to 3.5e5 s; CW 0 FATAL 0
  - W1R2: dt 6.0556 s, 0.330 s/cycle (cycles 50-450), estimate 5.3 h to 3.5e5 s; CW 0 FATAL 0
  - W2's dt halves (6.07 s): the lateral cell sets it; W1R2's dt 6.06 s (radial). W1R2 is below 12 h, so its tlim stays 3.5e5 s. W2 at about 11 h ends inside its holder (7333035, to 11-26 23:16).

## Readings at the W1 / W05 end (W2 and W1R2 still running)
- Run health:
  - W05: rc 0, wall 5750 s, CW 0, FATAL 0, DT_RETRY 0; dt 12.13 s.
  - W1: rc 0, wall 8662 s, CW 0, FATAL 0, DT_RETRY 0; dt 12.13 s; 0.30 s/cycle late in the run.
  - W1E (fresh HLLE counterpart of W1) started at 00:10 on 02-03.
- Full tables (items 1-5, every 1e4 s, A/B windows) and item 6 are in docs/handover/caltech_results/HLLRES/w1w05/. Figures F1-F4 (idort_style) are there too: W05 and W1 against prod_c1 and M35C.

| at 0.74 R_ph,IC | W05 (64, 10 deg) | M35C (128, 20 deg) | W1 (128, 10 deg) | prod_c1 (256, 20 deg) |
|---|---|---|---|---|
| lateral cell vs production | 2x | 2x | 1x | 1x |
| peak k20 at 1.8e5 / 2.6e5 / 3.46e5 s | 14 / 14 / 14 | 15 / 14 / 14 | 26 / 26 / 26 | 26 / 25 / 23 |
| cells per wavelength at 3.46e5 | 9.1 | 9.1 | 9.8 | 11.1 |
| band rms at 3.46e5 [km/s] (peak-centred band) | 8.4e-2 | 8.4e-2 | 1.6e-2 | 2.0e-2 |
| e-fold 1-2e5 / 2-2.6e5 / 2.6-3.5e5 s [1e4 s] | 2.11 / 3.80 / 6.43 | 2.22 / 3.49 / 5.31 | 2.54 / 2.53 / 5.60 | 2.44 / 2.55 / 5.40 |
| ln-residual rms of the fits | 0.26-0.38 | 0.07-0.09 | 0.07-0.13 | 0.05-0.08 |
| band axis share (counting share) | 0.98-0.99 (0.28) | 0.64-0.87 (0.13) | 0.90-0.97 (0.14) | 0.55-0.65 (0.08) |
| radial: peak r, FWHM at 3.46e5 | 0.754, 0.023 R (23 cells) | 0.747, 0.054 R (54 cells) | 0.735, 0.047 R (47 cells) | 0.747, 0.047 R (48 cells) |
| node spacing of the strongest coefficient at 3.46e5 [R_ph] | (1 node) | 0.031 | 0.015 / 0.019 / 0.021 | 0.015 / 0.018 / 0.021 |
| hst 2-3.5e5 s: KE_int / KEr_int [1e39 erg/sr], L_top/L_in | 11.7 / 5.19, 0.9996 | - | 12.6 / 5.57, 1.0015 | 11.8 / 5.16, 1.0003 |

- **W1 vs prod_c1 (the wedge-change control)** passes both of your bounds:
  - growth e-fold within 4 % in every window (bound 25 %);
  - peak k20 26 against 23-26 (bound 20 %: within 13 % at 3.46e5 s, equal at 1.8e5 s);
  - FWHM 47 vs 48 cells; node spacing equal to 0.001 R;
  - amplitude 0.8x at 3.46e5 s;
  - the axis share is higher (0.97 vs 0.65), from the 10 deg periodic box (its counting share is 0.14 vs 0.08).
  - So the narrow-wedge ladder is usable.
- **Lateral (W05 vs W1, so far):** halving the lateral resolution halves the peak k20 (26 -> 14), the same 9-10 cells per wavelength. That is the lateral grid-artefact signature of your fixed reading. W05 reproduces M35C (14, 9.1 cpw, same amplitude) on a different wedge.
  - The growth rate barely moves (e-fold 2.1 vs 2.5e4 s early). The amplitude at equal t is 5x larger on the coarse grid.
  - The verdict waits for W2 (does the peak move to k20 ~50?).
- **Radial:** W05's band is half as wide in R as W1's (23 vs 47 cells at the same radial grid) and has one node. So the radial structure also follows the lateral scale. W1R2 (radial halving) will say whether FWHM is fixed in cells or in R.
- **Item 5 (frequency):** unreliable. The strongest coefficient changes mode between dumps, so only 9 (W1) and 7 (prod) same-mode dumps after 2e5 s remain, with gaps. Phase residuals are 1.4-1.6 rad; |omega|/N ~0.5-0.7 with large error. W05 gives 1.81e-4 rad/s = 1.13 N (68 dumps, residual 0.7 rad). The 1800 s cadence aliases above 1.75e-3 rad/s.
- F4: the v_r maps at 0.74 show the checkerboard of axis-aligned modes, one with the W1 cell and one with twice the cell for W05. Colour ranges +-0.072 (W05) and +-0.042 km/s (W1).

## Readings at the W1R2 / W1E end (W2 still running, ends about 09:00 PST)
- Run health:
  - W1R2: rc 0, wall 24804 s, CW 0, FATAL 0, DT_RETRY 0; dt 6.06 s.
  - W1E (fresh HLLE, W1's grid): rc 0, wall 8647 s, CW 0, FATAL 0, DT_RETRY 0.
- Full tables and figures F1-F4: docs/handover/caltech_results/HLLRES/w1r2_w1e/.

| at 0.74 R_ph,IC | W1 (1x radial) | W1R2 (0.5x radial cell) | W1E (HLLE) | prod_c1 |
|---|---|---|---|---|
| peak k20 at 1.8e5 / 3.46e5 s | 26 / 26 | 26 / 22 | 8 / 8 (no band) | 26 / 23 |
| cells per wavelength at 3.46e5 | 9.8 | 11.6 | - | 11.1 |
| band rms at 3.46e5 [km/s] | 1.6e-2 | 6.9e-2 (peak-band) | band absent: B window 4e-6 | 2.0e-2 |
| e-fold 1-2e5 / 2-2.6e5 / 2.6-3.5e5 [1e4 s] | 2.54 / 2.53 / 5.60 | 1.94 / 2.73 / 4.41 | no growth | 2.44 / 2.55 / 5.40 |
| radial FWHM at 2.6e5 s | 0.043 R = 44 cells | 0.047 R = 95 cells | - | 0.043 R = 44 cells |
| radial FWHM at 3.46e5 s | 0.047 R = 47 cells | 0.028 R = 56 cells | - | 0.047 R = 48 cells |
| node spacing [R_ph] | 0.015-0.021 | 0.013-0.017 | - | 0.015-0.021 |
| hst KE_int / KEr_int [1e39 erg/sr], L_top/L_in | 12.6 / 5.57, 1.0015 | 12.4 / 5.41, 1.0001 | 11.7 / 5.05, 1.0012 | 11.8 / 5.16, 1.0003 |

- **Radial (W1 vs W1R2):** at 2.6e5 s the band has the same width in R (0.043 vs 0.047 R_ph) with twice the cells in W1R2, and the node spacing in R is the same (0.013-0.021). By your reading that is **physical (fixed in R), not set by the radial grid**. At 3.46e5 W1R2's band is narrower (0.028 R) and 2-4x stronger: halving the radial cell makes the band grow faster (e-fold 1.94 vs 2.54e4 s early) rather than vanish.
- **Lateral cells per wavelength at radial refinement:** k20 stays 26 to 2.6e5 s and drifts to 22 at 3.46e5 (11.6 cpw). The lateral scale still follows the lateral grid, as W05 vs W1 showed.
- **W1E (fresh HLLE on W1's grid):** no band at 0.74 at all. The peak sits at k20 8 with rms about 1e-5 km/s and no growth (e-folds negative or 1.7e5 s); the convective/shell power is at r 0.85 R_ph. This confirms on a fresh start what R6 (HLLE) showed on the restart: the band needs HLLC.
- **Frequency (item 5):** W1R2 keeps one strongest mode (11, 0) through 70 dumps after 2e5 s: |omega| = 1.86e-4 rad/s = 1.16 N, phase residual 0.25 rad. W05 gives 1.13 N likewise. So the band oscillates at about N at 0.74 (alias caveat: 1800 s cadence, Nyquist 1.75e-3 rad/s).
- W2 (lateral halving) will decide the lateral question; final readings when it ends.
