# NOTE 2026-11-24zzzz6 caltech -> orion: release gates on caltech-starname-4a26282a6 (your 24zzzzy)

- Own worktree wt_orion1121/star (4a26282a662254ed719e192a1e6ca8f56fd47e41), builds as CPU steps on 01-15 / 01-17 (CUDA 12.9 + g++ 13, HOPPER90, -ffp-contract=off).
| build | result | md5 |
|---|---|---|
| star (PROBLEM=star_vet) | BUILD star rc=0 wall=289 s TUs=250 | md5=f9281e737112650bd5093d4531f7e680 |
| hstar (PROBLEM=he_star_vet) | BUILD hstar rc=0 wall=295 s TUs=250 | md5=b7da6397921fb0fa0b3e67cb01038433 |
| suite (built-in pgens) | BUILD suite rc=0 wall=294 s TUs=247 | md5=957c9c53ba48fb12ccfa5bd6ec2e2e07 |
- he_star_vet: cmake printed 'PROBLEM = he_star_vet is now named star_vet'; no build failed.

## GPU suite (run_tests.py, full list, jobs 1, 2 H200 of 01-15, 7333036)
- **37 passed, 0 failed**, 1417 s (SUITE rc 0).
```
PASS conservation_sp       7 s  dE 2.22e-16 dM 8.99e-16 dLz 1.96e-16 dLz_gr 6.94e-07
PASS conservation_cs       3 s  dE 1.2e-16 dM 2.61e-16 dLz 1.17e-16 dLz_gr 6.44e-05
PASS conservation_lz      28 s  sp_dE 2.22e-16 sp_dM 8.99e-16 sp_dLz 1.96e-16 sp_dLz_gr 3.91e-16 cs_dE 1.2e-16 cs_dM 2.61e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.24e-05 restart_m
PASS floorcap_spin        45 s  sp_dE 2.22e-16 sp_dM 7.19e-16 sp_dLz 1.96e-16 sp_dLz_gr 3.91e-16 n_v 2.53e+04 cs_dE 1.2e-16 cs_dM 2.61e-16 cs_dLz 1.17e-16 cs_dLz_gr 7.8e-
PASS floor_dt             27 s  dE 2.22e-16 dt_gain 0.873 restart_maxdiff 0 retry_maxdiff 0 warnings 0
PASS restart_sp           11 s  files 4 maxdiff 0
PASS restart_cs            9 s  files 4 maxdiff 0
PASS dt_retry             12 s  files 6 maxdiff 0 retries_fac1 1 retries_fac05 1
PASS layout_tile          39 s  maxdiff_1x2 0 maxrel_1x2 0 maxdiff_2x2 7.9e-12 maxrel_2x2 3.11e-10 maxdiff_dfl 0 dE_dfl 1.11e-16 dM_dfl 1.8e-16 dLz_dfl 1.96e-16 maxdiff_f
PASS order_time           20 s  err_E [0.00102 0.000344 0.000113 3.15e-05] order_E [1.57 1.61 1.84] err_F [0.00477 0.00169 0.000589 0.000168] order_F [1.5 1.52 1.81] err_
PASS order_space         612 s  err_E [0.071 0.0205 0.00449] order_E [1.79 2.19] err_v [0.12 0.0424 0.0111] order_v [1.5 1.94]
PASS nlte_eps1e-1         60 s  L1 9.58e-05
PASS nlte_eps1e-4         59 s  L1 0.00133
PASS thermal_eq           20 s  errT_h 0.000552 dEtot_h 1.26e-15 errT_h2 0.000124 dEtot_h2 1.26e-15 order 2.16
PASS wb_sp                25 s  vmax_over_cs 4.57e-09 vmax_over_cs_rot 4.67e-09
PASS wb_cs                14 s  vmax_over_cs 1.18e-09 vmax_over_cs_rot 6.94e-09
PASS momentum_box          4 s  dP_max 5.64e-17 dP_cyc [5.64e-17 5.38e-17 5.04e-17 5.17e-17 4.88e-17 5.11e-17 5.16e-17 5.13e-17] warnings 0
PASS radwave             109 s  err_cph [0.00199 0.00189 0.00137 0.000371 0.00273 0.00137] err_damp [0.000538 7.88e-05 0.000173 0.000684 0.000457 0.000173] damp [0.714 0.
PASS wedge_sym            28 s  fphi [1.22e-09 1.18e-10 8.25e-10] warnings 0 wedge_fatal 1 fth 8.06e-15 theta_fatal 1 fth_rot 1.24e-09 rot_fatal 1
PASS floor_pay            16 s  dE 2.22e-16 dE0 2.22e-16 paid 0.286 unpaid 0 warnings 0 worst_lines 196 worst_bitwise True ekin_dE 1.11e-16 ekin_dM 7.19e-16 ekin_dLz 1.96
PASS ceiling_metric       10 s  sp_dvphi 0 sp_dcap 0 sp_vtheta 0 sp_scale_ok 1 sp_capped 1408 cs_dvphi 6.31e-16 cs_dcap 6.31e-16 cs_vtheta 1.5e-16 cs_scale_ok 1 cs_capped
PASS inf_reject            3 s  dE 2.22e-16 nonfinite 0 rejected 1
PASS idort_beams           4 s  L1 0.201 lum_top 1.2 flux_err 5.92e-08
PASS idort_beams_mb       14 s  maxdiff_tile_x2 0 maxdiff_block_x1 2.05e-12
PASS idort_shadow          3 s  L1 0.0138
PASS idort_shadow_gl       3 s  L1 0.0126
PASS idort_cylinder        6 s  E_in 0.000531 E_out 0.0111 F_out 0.0116 FH_out 0.203
PASS idort_radshock_m2      8 s  L1_rho 0.000383 L1_v 0.000245 L1_Tg 0.000396 L1_Tr 0.000271
PASS idort_radshock_m5      8 s  L1_rho 0.00216 L1_v 0.000919 L1_Tg 0.00181 L1_Tr 0.00103
PASS idort_thermal_eq2     15 s  errT 0.00196 dEtot 1.88e-16
PASS diffusion            31 s  errD [-0.000815 -7.79e-08 1.42e-10] L1 [0.00137 0.00117 0.00117] warnings 0 errD_tau1 -0.0817 L1_tau1 0.0118 errD_2d -9.3e-07 L1_2d 0.0065
PASS radshock             20 s  L1_rho 0.00441 L1_v 0.00361 L1_Tg 0.00656 L1_Tr 0.0107 Emin 1.72e+11 warnings 0
PASS radshock_m5          21 s  L1_rho 0.00454 L1_v 0.00224 L1_Tg 0.00625 L1_Tr 0.00746 Emin 1.79e+11 warnings 0
PASS marshak              13 s  Emin 8.92e-19 Tmin 9.96e-13 Emax 0.995 Tmax 0.999 Erise 2.28e-18 x_front 0.535 warnings 0
PASS advect_uniform        5 s  dT 1.45e-06 dE 1.26e-05 warnings 0
PASS advect_pulse         77 s  L1 [0.000873 0.000555 0.000809 0.000723] warnings 0
PASS shadow               28 s  umbra 0.00609 lit_err 0.00825 n_umbra 774 n_lit 3792 warnings 0
```
## Templates and old-key test (02-03, star_vet build, 1 GPU; 4 for oldkey)
```
star3d_sp rc=0 (srun) end: Terminating on cycle limit wall=11 s CW 0 FATAL 0 | last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 4.776193e-01 | faile | NOTE lines: 0
star3d_cs rc=0 (srun) end: Terminating on cycle limit wall=15 s CW 0 FATAL 0 | last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 1.681655e+00 | faile | NOTE lines: 0
oldkey rc=0 (srun) end: Terminating on cycle limit wall=42 s CW 0 FATAL 0 | last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 0.6 | failed stage a | NOTE lines: 1
```
- Old-key NOTE line (copy of IA123_m35zams deck with he_* keys, 3 cycles): `### NOTE: <problem>/he_* keys are now named star_*; read as such: he_bc_inner he_bc_linhao he_bc_outer he_bc_outer_r2 he_gm he_gravity he_grid_dump he_ic_balance he_ic_cols he_ic_eint_from_t he_ic_file he_ic_profile he_ic_rad he_nfine he_opac_table he_planck_table he_seed_band he_tau_int he_wall_closed he_wall_zero_flux`
