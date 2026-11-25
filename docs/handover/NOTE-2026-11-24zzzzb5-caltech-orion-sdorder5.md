# NOTE 2026-11-24zzzzb5 caltech -> orion: sdorder5 (10bf5513c): CUDA builds, GPU suite, old-restart check (your 24zzzzb5)

## 1. Builds (own worktree wt_orion1121/sdorder5, CUDA 12.9 + g++ 13, HOPPER90, 02-09)
- BUILD sdorder5 rc=0 wall=318 s TUs=250
- BUILD_OK sdorder5 md5=9bd078772d1dae6c179a6bf6107f688d prov=10bf5513c8b40fca1f0d9b17fc06ece51a53fdc5
- GATE PASS sdorder5: athena_sdorder5_10bf5513c md5 9bd07877 node hpc-sm-02-09 pass1=3.686e-11 log date 11-24 18:03:12 end: Terminating on cycle limit
- BUILD suite rc=0 wall=315 s
- BUILD_OK suite md5=4a68156564e775e69dbcfd7440f06f97 prov=10bf5513c8b40fca1f0d9b17fc06ece51a53fdc5

## 2. GPU suite (full list, jobs 1, 2 H200 of 02-09): 37 passed, 1 failed, 1392 s (38 expected)
```
PASS conservation_sp       6 s  dE 1.11e-16 dM 8.99e-16 dLz 3.91e-16 dLz_gr 6.94e-07
PASS conservation_cs       3 s  dE 1.2e-16 dM 2.61e-16 dLz 1.17e-16 dLz_gr 6.44e-05
PASS conservation_lz      28 s  sp_dE 2.22e-16 sp_dM 7.19e-16 sp_dLz 1.96e-16 sp_dLz_gr 3.91e-16 cs_dE 1.2e-16 cs_dM 2.61e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.24e-05 restart_maxdiff 0 retry_maxdi
PASS floorcap_spin        45 s  sp_dE 2.22e-16 sp_dM 5.4e-16 sp_dLz 1.96e-16 sp_dLz_gr 3.91e-16 n_v 2.59e+04 cs_dE 1.2e-16 cs_dM 1.3e-16 cs_dLz 1.17e-16 cs_dLz_gr 7.8e-05 nospin_sp_dLz_gr 0.
PASS floor_dt             27 s  dE 1.11e-16 dt_gain 0.873 restart_maxdiff 0 retry_maxdiff 0 warnings 0
PASS restart_sp           11 s  files 4 maxdiff 0
PASS restart_cs            8 s  files 4 maxdiff 0
PASS dt_retry             12 s  files 6 maxdiff 0 retries_fac1 1 retries_fac05 1
PASS layout_tile          39 s  maxdiff_1x2 0 maxrel_1x2 0 maxdiff_2x2 9.09e-13 maxrel_2x2 7.48e-11 maxdiff_dfl 0 dE_dfl 2.22e-16 dM_dfl 1.8e-16 dLz_dfl 1.96e-16 maxdiff_fofc 0 dE_fofc 2.22e
PASS order_time           20 s  err_E [0.00102 0.000344 0.000113 3.15e-05] order_E [1.57 1.61 1.84] err_F [0.00475 0.00168 0.000588 0.000168] order_F [1.5 1.52 1.81] err_T [0.000167 6.61e-05
PASS order_stiff          24 s  err_E [4.31e-06 1.13e-06 2.9e-07] order_E [1.94 1.96] err_v [0.00449 0.00118 0.000305] order_v [1.92 1.96]
PASS order_space         573 s  err_E [0.071 0.0203 0.00444] order_E [1.81 2.19] err_v [0.12 0.0422 0.011] order_v [1.5 1.94]
PASS nlte_eps1e-1         59 s  L1 9.58e-05
PASS nlte_eps1e-4         60 s  L1 0.00133
PASS thermal_eq           20 s  errT_h 0.000552 dEtot_h 1.26e-15 errT_h2 0.000124 dEtot_h2 1.26e-15 order 2.16
PASS wb_sp                25 s  vmax_over_cs 4.57e-09 vmax_over_cs_rot 4.67e-09
PASS wb_cs                15 s  vmax_over_cs 1.18e-09 vmax_over_cs_rot 6.94e-09
PASS momentum_box          4 s  dP_max 5.64e-17 dP_cyc [5.64e-17 5.27e-17 5.23e-17 5.42e-17 5.45e-17 5.21e-17 5.35e-17 5.25e-17] warnings 0
FAIL radwave             102 s  err_cph [0.00303 0.00189 0.00137 0.000364 0.00277 0.00137] err_damp [0.00235 0.000105 4.01e-05 0.000208 0.000545 4.02e-05] damp [0.711 0.164 0.0151 0.000302 0
PASS wedge_sym            29 s  fphi [3.11e-09 9.79e-12 3.11e-10] warnings 0 wedge_fatal 1 fth 6.91e-15 theta_fatal 1 fth_rot 1.56e-09 rot_fatal 1
PASS floor_pay            15 s  dE 2.22e-16 dE0 2.22e-16 paid 0.285 unpaid 0 warnings 0 worst_lines 195 worst_bitwise True ekin_dE 1.11e-16 ekin_dM 5.4e-16 ekin_dLz 3.91e-16
PASS ceiling_metric       10 s  sp_dvphi 0 sp_dcap 0 sp_vtheta 0 sp_scale_ok 1 sp_capped 1408 cs_dvphi 6.31e-16 cs_dcap 6.31e-16 cs_vtheta 1.5e-16 cs_scale_ok 1 cs_capped 2928
PASS inf_reject            4 s  dE 2.22e-16 nonfinite 0 rejected 1
PASS idort_beams           5 s  L1 0.777 lum_top 1.79 flux_err 2.66e-07
PASS idort_beams_mb       16 s  maxdiff_tile_x2 0 maxdiff_block_x1 4.32e-12
PASS idort_shadow          3 s  L1 0.0138
PASS idort_shadow_gl       4 s  L1 0.0126
PASS idort_cylinder        7 s  E_in 0.000531 E_out 0.0111 F_out 0.0116 FH_out 0.203
PASS idort_radshock_m2      9 s  L1_rho 0.000356 L1_v 0.000213 L1_Tg 0.000385 L1_Tr 0.00026
PASS idort_radshock_m5      8 s  L1_rho 0.0017 L1_v 0.000818 L1_Tg 0.00162 L1_Tr 0.00094
PASS idort_thermal_eq2     15 s  errT 0.00196 dEtot 1.88e-16
PASS diffusion            30 s  errD [-0.000815 -7.79e-08 1.42e-10] L1 [0.00137 0.00117 0.00117] warnings 0 errD_tau1 -0.0817 L1_tau1 0.0118 errD_2d -9.3e-07 L1_2d 0.00659 th_min 1.15e-17 th
PASS radshock             19 s  L1_rho 0.00438 L1_v 0.00358 L1_Tg 0.00659 L1_Tr 0.0107 Emin 1.72e+11 warnings 0
PASS radshock_m5          20 s  L1_rho 0.00391 L1_v 0.00212 L1_Tg 0.0062 L1_Tr 0.00744 Emin 1.79e+11 warnings 0
PASS marshak              12 s  Emin 8.92e-19 Tmin 9.96e-13 Emax 0.995 Tmax 0.999 Erise 8.89e-19 x_front 0.535 warnings 0
PASS advect_uniform        4 s  dT 2.73e-06 dE 2.35e-05 warnings 0
PASS advect_pulse         76 s  L1 [0.000873 0.000575 0.00082 0.000775] warnings 0
PASS shadow               27 s  umbra 0.00609 lit_err 0.00825 n_umbra 774 n_lit 3792 warnings 0
```

## 3. Old-restart check: t9l62.00015, 5 cycles (nlim 5284), T2 args (those of Compton C0: IA123_t9l62_restart.args + dt_report=10, outputs off), time2_gtime_slope not named, 01-17, 4 ranks, T2 7d7c4a71 then this binary
- T2: 2026-11-24 18:36:18 rc=0 (srun) ended: Terminating on cycle limit wall=82 s CW 0 FATAL 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- SD5: 2026-11-24 18:37:39 rc=0 (srun) ended: Terminating on cycle limit wall=81 s CW 0 FATAL 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- time2_gtime_slope lines in the new binary's log: `<rad_vet> time2_gtime_slope = false (not in the restart file)|`

**t9l62.hydro.hst**: rows 5 / 5, columns 10; file byte-identical: True
- columns whose printed digits are identical in every row: 10 of 10

**t9l62.user.hst**: rows 5 / 5, columns 32; file byte-identical: True
- columns whose printed digits are identical in every row: 32 of 32
