# NOTE 2026-11-24zzzzj4 caltech -> orion: production-restart test of rel-f b09a7b3c1 (your 24zzzzj4)

## m35zams (11-24 22:24): rst m35zams_m35zams.00003.rst (hard link, c50816), rel-f 80d16d39bd1f47bd3f22b5711c660519 against athena_c4speed_d979d5268; own keys: rad_vet/time2_gtime_slope=true rad_vet/clean_stag_thinf=true rad_vet/clean_stag_lat_sums=true rad_vet/clean_stag_tie_gid=true rad_vet/clean_stag_opac_drag=true rad_vet/vet_gd_allow_clamp=false hydro/vceil=3.0e8 rad_vet/clean_pc_sweep=rb rad_vet/clean_stag_mg_levels=2 time/cfl_number=0.3 rad_vet/clean_conv_max=1e-5
- R1: 2026-11-24 22:21:32 rc=0 (srun) ended: Terminating on cycle limit wall=25 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- R50: 2026-11-24 22:22:07 rc=0 (srun) ended: Terminating on cycle limit wall=34 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- R25a: 2026-11-24 22:22:46 rc=0 (srun) ended: Terminating on cycle limit wall=38 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- R25b: 2026-11-24 22:23:28 rc=0 (srun) ended: Terminating on cycle limit wall=31 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- C50: 2026-11-24 22:24:08 rc=0 (srun) ended: Terminating on cycle limit wall=40 s CW 0 FATAL 0 inf 0 DT_RETRY 0
#### NOTE lines of the restart (rel-f, 1 step)
- ### NOTE: <problem>/he_* (hs_*) keys are now named star_*; read so: he_bc_inner he_bc_linhao he_bc_outer he_bc_outer_r2 he_esrc_const he_gm he_gravity he_grid_dump he_ic_balance he_ic_cols he_ic_eint_from_t he_ic_file he_ic_profile he_ic_rad he_ic_rad_atm_post
- ### NOTE: <hydro>/fofc_flag_exchange = true (restart file without it: the fix applies; name fofc_flag_exchange=false to keep the old behaviour)
- ### NOTE: <rad_vet>/time2_gtime_pos = true (restart file without it: the fix applies; name time2_gtime_pos=false to keep the old behaviour)
- NOTE lines: 3; "the fix applies": 2; "added from the command line": 0; "old behaviour": 2; FATAL 0, CLEAN WARNING 0; clean_init lines: []
#### 50 straight against 25 + rst + 25
- user.hst: R50 50 rows, R25a 25, R25b 25; rows not printed-identical to R50 at the same t: 0 (first t -) -> BITWISE (printed)
- hydro.hst: R50 50 rows, R25a 25, R25b 25; rows not printed-identical to R50 at the same t: 0 (first t -) -> BITWISE (printed)
- final hydro_w dump, payload md5: R50 ced72921829d, R25b ced72921829d -> IDENTICAL
- final vet dump, payload md5: R50 afc7bad49c89, R25b afc7bad49c89 -> IDENTICAL
#### against the running binary (C50, same rst, same keys, 50 cycles)
- row 0: t 1.76869032e+05 / 1.76869032e+05, dt 3.480511 / 3.480511, L_top rel 3.21e-12, KE_int rel 1.11e-15, M_tot rel 0.00e+00, Etot rel 8.88e-16
- row 1: t 1.76872513e+05 / 1.76872513e+05, dt 3.480511 / 3.480511, L_top rel 1.82e-10, KE_int rel -1.16e-13, M_tot rel -6.66e-16, Etot rel 1.78e-15
- row 2: t 1.76875993e+05 / 1.76875993e+05, dt 3.480511 / 3.480511, L_top rel -1.94e-10, KE_int rel -1.14e-13, M_tot rel -1.78e-15, Etot rel -1.11e-16
- row 3: t 1.76879474e+05 / 1.76879474e+05, dt 3.480511 / 3.480511, L_top rel 1.09e-10, KE_int rel -6.92e-14, M_tot rel 2.22e-15, Etot rel -4.44e-16
- row 4: t 1.76882954e+05 / 1.76882954e+05, dt 3.480511 / 3.480511, L_top rel -1.57e-10, KE_int rel -4.35e-14, M_tot rel 1.78e-15, Etot rel 1.33e-15
- row 9: t 1.76900357e+05 / 1.76900357e+05, dt 3.480511 / 3.480511, L_top rel 3.36e-09, KE_int rel -2.17e-12, M_tot rel -9.99e-16, Etot rel -2.89e-15
- row 24: t 1.76952565e+05 / 1.76952565e+05, dt 3.480511 / 3.480511, L_top rel -6.55e-07, KE_int rel 6.13e-09, M_tot rel 4.44e-16, Etot rel -7.57e-14
- row 49: t 1.77039577e+05 / 1.77039577e+05, dt 3.480511 / 3.480511, L_top rel 9.07e-07, KE_int rel 1.45e-08, M_tot rel 2.22e-16, Etot rel -2.47e-13
- R50: s/cyc (elapsed prints) 0.2937 (c50820-50860); its sum per cycle 35.8, solves 4.48; CW 0 FATAL 0 inf 0 DT_RETRY 0
- C50: s/cyc (elapsed prints) 0.2925 (c50820-50860); its sum per cycle 35.7, solves 4.44; CW 0 FATAL 0 inf 0 DT_RETRY 0
#### END-row residuals over the 50 cycles (R_E = dEtot - dt (L_in - L_top + P_esrc - P_sponge - E_hflx); R_M = dM_tot + dt M_hflx; bsg: R_L = dLz_tot + dt Lz_hflx)
- R50: R_E / Etot RMS 9.96e-16 max 2.55e-15 sum -2.12e-16; R_M / M_tot RMS 1.09e-15 max 2.23e-15 sum -6.20e-16
- C50: R_E / Etot RMS 1.22e-15 max 2.98e-15 sum 6.74e-16; R_M / M_tot RMS 1.29e-15 max 2.98e-15 sum -8.68e-16
