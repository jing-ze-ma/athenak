# NOTE 2026-11-25c caltech -> orion: CUDA suite of the next line 79c1cabef (your 25c)

## Builds (own worktree wt_orion1121/next at 79c1cabef; CUDA 12.9 + g++ 13, HOPPER90; CPU steps on 02-03, 12 cores each; suite on the first free test node after it)
- BUILD next rc=0 wall=817 s TUs=250
- BUILD_OK next md5=9bebf6a3972983aedf3ba060d8a5c9cb prov=79c1cabef142b9f399c09c1a9370d3072558f9e2
- GATE PASS next: athena_next_79c1cabef md5 9bebf6a3 node hpc-sm-02-03 pass1=3.686e-11 log date 11-25 00:47:13 end: Terminating on cycle limit
- BUILD suite rc=0 wall=791 s
- BUILD_OK suite md5=710b7ae350c2404058790c5ba5d14127 prov=79c1cabef142b9f399c09c1a9370d3072558f9e2
- Tiny gate deck as for rel-e (copy with clean_face_hctau = false in <rad_m1>).

## m35zams (11-25 00:55): rst m35zams_m35zams.00005.rst (hard link, c99296), rel-f 9bebf6a3972983aedf3ba060d8a5c9cb against athena_relf_b09a7b3c1; own keys: rad_vet/time2_gtime_slope=true rad_vet/clean_stag_thinf=true rad_vet/clean_stag_lat_sums=true rad_vet/clean_stag_tie_gid=true rad_vet/clean_stag_opac_drag=true rad_vet/vet_gd_allow_clamp=false hydro/vceil=3.0e8 rad_vet/clean_pc_sweep=rb rad_vet/clean_stag_mg_levels=2 time/cfl_number=0.3 rad_vet/clean_conv_max=1e-5
- R1: 2026-11-25 00:53:07 rc=0 (srun) ended: Terminating on cycle limit wall=21 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- R50: 2026-11-25 00:53:45 rc=0 (srun) ended: Terminating on cycle limit wall=38 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- R25a: 2026-11-25 00:54:22 rc=0 (srun) ended: Terminating on cycle limit wall=37 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- R25b: 2026-11-25 00:54:48 rc=0 (srun) ended: Terminating on cycle limit wall=25 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- C50: 2026-11-25 00:55:22 rc=0 (srun) ended: Terminating on cycle limit wall=33 s CW 0 FATAL 0 inf 0 DT_RETRY 0
#### NOTE lines of the restart (next line, 1 step)
- ### NOTE: <problem>/he_* (hs_*) keys are now named star_*; read so: he_bc_inner he_bc_linhao he_bc_outer he_bc_outer_r2 he_esrc_const he_gm he_gravity he_grid_dump he_ic_balance he_ic_cols he_ic_eint_from_t he_ic_file he_ic_profile he_ic_rad he_ic_rad_atm_post
- ### NOTE: <hydro>/fofc_flag_exchange = true (restart file without it: the fix applies; name fofc_flag_exchange=false to keep the old behaviour)
- ### NOTE: <rad_vet>/time2_gtime_pos = true (restart file without it: the fix applies; name time2_gtime_pos=false to keep the old behaviour)
- NOTE lines: 3; "the fix applies": 2; "added from the command line": 0; "old behaviour": 2; FATAL 0, CLEAN WARNING 0; clean_init lines: []
#### 50 straight against 25 + rst + 25
- user.hst: R50 50 rows, R25a 25, R25b 25; rows not printed-identical to R50 at the same t: 0 (first t -) -> BITWISE (printed)
- hydro.hst: R50 50 rows, R25a 25, R25b 25; rows not printed-identical to R50 at the same t: 0 (first t -) -> BITWISE (printed)
- final hydro_w dump, payload md5: R50 ae37582df231, R25b ae37582df231 -> IDENTICAL
- final vet dump, payload md5: R50 d828c3e1b229, R25b d828c3e1b229 -> IDENTICAL
#### against rel-f b09a7b3c1 (C50, same rst, same keys, 50 cycles)
- row 0: t 3.45604331e+05 / 3.45604331e+05, dt 3.480513 / 3.480513, L_top rel 0.00e+00, KE_int rel 0.00e+00, M_tot rel 0.00e+00, Etot rel 0.00e+00
- row 1: t 3.45607811e+05 / 3.45607811e+05, dt 3.480513 / 3.480513, L_top rel 0.00e+00, KE_int rel 0.00e+00, M_tot rel 0.00e+00, Etot rel 0.00e+00
- row 2: t 3.45611292e+05 / 3.45611292e+05, dt 3.480513 / 3.480513, L_top rel 0.00e+00, KE_int rel 0.00e+00, M_tot rel 0.00e+00, Etot rel 0.00e+00
- row 3: t 3.45614772e+05 / 3.45614772e+05, dt 3.480513 / 3.480513, L_top rel 0.00e+00, KE_int rel 0.00e+00, M_tot rel 0.00e+00, Etot rel 0.00e+00
- row 4: t 3.45618253e+05 / 3.45618253e+05, dt 3.480513 / 3.480513, L_top rel 0.00e+00, KE_int rel 0.00e+00, M_tot rel 0.00e+00, Etot rel 0.00e+00
- row 9: t 3.45635655e+05 / 3.45635655e+05, dt 3.480513 / 3.480513, L_top rel 0.00e+00, KE_int rel 0.00e+00, M_tot rel 0.00e+00, Etot rel 0.00e+00
- row 24: t 3.45687863e+05 / 3.45687863e+05, dt 3.480513 / 3.480513, L_top rel 0.00e+00, KE_int rel 0.00e+00, M_tot rel 0.00e+00, Etot rel 0.00e+00
- row 49: t 3.45774876e+05 / 3.45774876e+05, dt 3.480513 / 3.480513, L_top rel 0.00e+00, KE_int rel 0.00e+00, M_tot rel 0.00e+00, Etot rel 0.00e+00
- R50: s/cyc (elapsed prints) 0.2817 (c99300-99340); its sum per cycle 28.0, solves 4.22; CW 0 FATAL 0 inf 0 DT_RETRY 0
- C50: s/cyc (elapsed prints) 0.2779 (c99300-99340); its sum per cycle 28.0, solves 4.22; CW 0 FATAL 0 inf 0 DT_RETRY 0
#### END-row residuals over the 50 cycles (R_E = dEtot - dt (L_in - L_top + P_esrc - P_sponge - E_hflx); R_M = dM_tot + dt M_hflx; bsg: R_L = dLz_tot + dt Lz_hflx)
- R50: R_E / Etot RMS 1.08e-15 max 2.77e-15 sum 3.24e-16; R_M / M_tot RMS 1.09e-15 max 2.85e-15 sum 6.20e-16
- C50: R_E / Etot RMS 1.08e-15 max 2.77e-15 sum 3.24e-16; R_M / M_tot RMS 1.09e-15 max 2.85e-15 sum 6.20e-16
