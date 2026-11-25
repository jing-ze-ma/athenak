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

## m35mms (11-24 22:30): rst m35mms_m35mms.00004.rst (hard link, c22144), rel-f 80d16d39bd1f47bd3f22b5711c660519 against athena_c4speed_d979d5268; own keys: rad_vet/time2_gtime_slope=true rad_vet/clean_stag_thinf=true rad_vet/clean_stag_lat_sums=true rad_vet/clean_stag_tie_gid=true rad_vet/clean_stag_opac_drag=true rad_vet/vet_gd_allow_clamp=false hydro/vceil=3.0e8 rad_vet/clean_pc_sweep=rb rad_vet/clean_stag_mg_levels=4 time/cfl_number=0.3 rad_vet/clean_conv_max=1e-5
- R1: 2026-11-24 22:25:12 rc=0 (srun) ended: Terminating on cycle limit wall=32 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- R50: 2026-11-24 22:26:41 rc=0 (srun) ended: Terminating on cycle limit wall=89 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- R25a: 2026-11-24 22:28:14 rc=0 (srun) ended: Terminating on cycle limit wall=93 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- R25b: 2026-11-24 22:29:17 rc=0 (srun) ended: Terminating on cycle limit wall=63 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- C50: 2026-11-24 22:30:45 rc=0 (srun) ended: Terminating on cycle limit wall=88 s CW 0 FATAL 0 inf 0 DT_RETRY 0
#### NOTE lines of the restart (rel-f, 1 step)
- ### NOTE: <problem>/he_* (hs_*) keys are now named star_*; read so: he_bc_inner he_bc_linhao he_bc_outer he_bc_outer_r2 he_esrc_const he_gm he_gravity he_grid_dump he_ic_balance he_ic_cols he_ic_eint_from_t he_ic_file he_ic_profile he_ic_rad he_ic_rad_atm_post
- ### NOTE: <hydro>/fofc_flag_exchange = true (restart file without it: the fix applies; name fofc_flag_exchange=false to keep the old behaviour)
- ### NOTE: <rad_vet>/time2_gtime_pos = true (restart file without it: the fix applies; name time2_gtime_pos=false to keep the old behaviour)
- NOTE lines: 3; "the fix applies": 2; "added from the command line": 0; "old behaviour": 2; FATAL 0, CLEAN WARNING 0; clean_init lines: []
#### 50 straight against 25 + rst + 25
- user.hst: R50 50 rows, R25a 25, R25b 25; rows not printed-identical to R50 at the same t: 0 (first t -) -> BITWISE (printed)
- hydro.hst: R50 50 rows, R25a 25, R25b 25; rows not printed-identical to R50 at the same t: 0 (first t -) -> BITWISE (printed)
- final hydro_w dump, payload md5: R50 0c015e3cac61, R25b 0c015e3cac61 -> IDENTICAL
- final vet dump, payload md5: R50 3f010f7b9af1, R25b 3f010f7b9af1 -> IDENTICAL
#### against the running binary (C50, same rst, same keys, 50 cycles)
- row 0: t 2.59212479e+05 / 2.59212479e+05, dt 12.006718 / 12.006718, L_top rel -1.68e-14, KE_int rel 3.11e-15, M_tot rel -1.55e-15, Etot rel -1.44e-15
- row 1: t 2.59224485e+05 / 2.59224485e+05, dt 12.006725 / 12.006725, L_top rel -1.37e-13, KE_int rel -3.00e-15, M_tot rel 2.00e-15, Etot rel 4.44e-16
- row 2: t 2.59236492e+05 / 2.59236492e+05, dt 12.006732 / 12.006732, L_top rel 2.39e-13, KE_int rel 2.22e-16, M_tot rel 3.11e-15, Etot rel -1.22e-15
- row 3: t 2.59248499e+05 / 2.59248499e+05, dt 12.006738 / 12.006738, L_top rel -1.94e-13, KE_int rel 1.49e-14, M_tot rel 1.11e-15, Etot rel 4.44e-16
- row 4: t 2.59260506e+05 / 2.59260506e+05, dt 12.006745 / 12.006745, L_top rel -3.53e-13, KE_int rel 2.51e-14, M_tot rel -2.22e-15, Etot rel 2.66e-15
- row 9: t 2.59320540e+05 / 2.59320540e+05, dt 12.006780 / 12.006780, L_top rel 7.97e-13, KE_int rel 1.22e-14, M_tot rel 3.77e-15, Etot rel 3.11e-15
- row 24: t 2.59500642e+05 / 2.59500642e+05, dt 12.006886 / 12.006886, L_top rel 7.58e-07, KE_int rel -7.18e-11, M_tot rel 0.00e+00, Etot rel -2.39e-13
- row 49: t 2.59800816e+05 / 2.59800816e+05, dt 12.007021 / 12.007021, L_top rel -9.81e-07, KE_int rel -1.41e-09, M_tot rel 2.22e-16, Etot rel 4.20e-12
- R50: s/cyc (elapsed prints) -; its sum per cycle 61.3, solves 9.26; CW 0 FATAL 0 inf 0 DT_RETRY 0
- C50: s/cyc (elapsed prints) -; its sum per cycle 60.8, solves 9.16; CW 0 FATAL 0 inf 0 DT_RETRY 0
#### END-row residuals over the 50 cycles (R_E = dEtot - dt (L_in - L_top + P_esrc - P_sponge - E_hflx); R_M = dM_tot + dt M_hflx; bsg: R_L = dLz_tot + dt Lz_hflx)
- R50: R_E / Etot RMS 2.29e-15 max 6.75e-15 sum -6.07e-15; R_M / M_tot RMS 2.53e-15 max 4.72e-15 sum -3.78e-16
- C50: R_E / Etot RMS 2.00e-15 max 4.42e-15 sum -1.20e-15; R_M / M_tot RMS 2.56e-15 max 6.61e-15 sum -2.08e-15

## bsg (11-24 22:35): rst bsg_bsg3d.00008.rst (hard link, c29127), rel-f 80d16d39bd1f47bd3f22b5711c660519 against athena_sdorder3_95813cf85; own keys: none
- R1: 2026-11-24 22:30:53 rc=0 (srun) ended: Terminating on cycle limit wall=31 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- R50: 2026-11-24 22:32:00 rc=0 (srun) ended: Terminating on cycle limit wall=67 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- R25a: 2026-11-24 22:33:09 rc=0 (srun) ended: Terminating on cycle limit wall=69 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- R25b: 2026-11-24 22:33:57 rc=0 (srun) ended: Terminating on cycle limit wall=48 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- C50: 2026-11-24 22:35:30 rc=0 (srun) ended: Terminating on cycle limit wall=93 s CW 0 FATAL 0 inf 0 DT_RETRY 0
#### NOTE lines of the restart (rel-f, 1 step)
- ### NOTE: <problem>/he_* (hs_*) keys are now named star_*; read so: he_bc_inner he_bc_linhao he_bc_outer he_bc_outer_r2 he_esrc_const he_gm he_gravity he_grid_dump he_ic_balance he_ic_cols he_ic_eint_from_t he_ic_file he_ic_profile he_ic_profile_amp he_ic_prof
- ### NOTE: <hydro>/merge_lat_mom = radial (restart file without it: the fix applies; name merge_lat_mom=rigid to keep the old behaviour)
- ### NOTE: <hydro>/fofc_flag_exchange = true (restart file without it: the fix applies; name fofc_flag_exchange=false to keep the old behaviour)
- ### NOTE: <rad_vet>/time2_gtime_pos = true (restart file without it: the fix applies; name time2_gtime_pos=false to keep the old behaviour)
- NOTE lines: 4; "the fix applies": 3; "added from the command line": 0; "old behaviour": 3; FATAL 0, CLEAN WARNING 0; clean_init lines: []
#### 50 straight against 25 + rst + 25
- user.hst: R50 50 rows, R25a 25, R25b 25; rows not printed-identical to R50 at the same t: 0 (first t -) -> BITWISE (printed)
- hydro.hst: R50 50 rows, R25a 25, R25b 25; rows not printed-identical to R50 at the same t: 0 (first t -) -> BITWISE (printed)
- final hydro_w dump, payload md5: R50 497fb089d38f, R25b 497fb089d38f -> IDENTICAL
- final vet dump, payload md5: R50 5632e73547e4, R25b 5632e73547e4 -> IDENTICAL
#### against the running binary (C50, same rst, same keys, 50 cycles)
- row 0: t 4.83873577e+06 / 4.83873577e+06, dt 169.427068 / 169.427068, L_top rel 9.40e-12, KE_int rel -3.04e-04, M_tot rel -4.44e-16, Etot rel 8.88e-16
- row 1: t 4.83890520e+06 / 4.83890520e+06, dt 169.427508 / 169.427523, L_top rel 1.38e-11, KE_int rel -3.30e-04, M_tot rel -4.00e-15, Etot rel 6.66e-16
- row 2: t 4.83907463e+06 / 4.83907463e+06, dt 169.427970 / 169.427985, L_top rel 2.23e-11, KE_int rel -3.55e-04, M_tot rel -3.22e-15, Etot rel 1.55e-15
- row 3: t 4.83924406e+06 / 4.83924406e+06, dt 169.428441 / 169.428456, L_top rel 1.61e-11, KE_int rel -3.80e-04, M_tot rel 3.77e-15, Etot rel 8.88e-16
- row 4: t 4.83941349e+06 / 4.83941349e+06, dt 169.428920 / 169.428935, L_top rel -1.38e-11, KE_int rel -4.05e-04, M_tot rel 1.33e-15, Etot rel 6.44e-15
- row 9: t 4.84026010e+06 / 4.84026010e+06, dt 169.079814 / 169.079814, L_top rel 2.43e-11, KE_int rel -5.20e-04, M_tot rel -2.44e-15, Etot rel 3.33e-15
- row 24: t 4.84278201e+06 / 4.84278201e+06, dt 167.299694 / 167.299693, L_top rel 5.31e-07, KE_int rel -7.97e-04, M_tot rel -2.00e-15, Etot rel -1.59e-13
- row 49: t 4.84691279e+06 / 4.84691292e+06, dt 163.504413 / 163.517619, L_top rel -1.81e-05, KE_int rel -1.12e-03, M_tot rel 1.55e-15, Etot rel 1.05e-11
- R50: s/cyc (elapsed prints) -; its sum per cycle 32.8, solves 6.84; CW 0 FATAL 0 inf 0 DT_RETRY 0
- C50: s/cyc (elapsed prints) -; its sum per cycle 33.0, solves 6.88; CW 0 FATAL 0 inf 0 DT_RETRY 0
#### END-row residuals over the 50 cycles (R_E = dEtot - dt (L_in - L_top + P_esrc - P_sponge - E_hflx); R_M = dM_tot + dt M_hflx; bsg: R_L = dLz_tot + dt Lz_hflx)
- R50: R_E / Etot RMS 2.86e-15 max 8.59e-15 sum -8.58e-17; R_M / M_tot RMS 2.59e-15 max 7.14e-15 sum 7.76e-16; R_L RMS 2.438e+31 max 5.471e+31 sum -8.090e+32
- C50: R_E / Etot RMS 3.00e-15 max 8.51e-15 sum 1.35e-15; R_M / M_tot RMS 2.53e-15 max 5.28e-15 sum -1.24e-15; R_L RMS 3.549e+33 max 1.111e+34 sum 8.744e+34
