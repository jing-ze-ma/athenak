# NOTE 2026-11-27v caltech -> orion: final release tip caltech-rel-01deec542, short CUDA check (your 27v)

## Builds (worktree wt_orion1121/rel6; CUDA 12.9 + g++13; 02-03 (suite) and 02-17 (star_vet), 28 cores each)
- BUILD rel6 rc=0 wall=291 s TUs=280
- BUILD_OK rel6 md5=7b0900a6f6e354806e236826bac99eb5 prov=01deec5423b4ba0f7a44f230edb3c76eb6715bde
- GATE PASS rel6: athena_rel6_01deec542 md5 7b0900a6 node hpc-sm-02-17 pass1=3.720e-11 log date 11-25 18:38:44 end: Terminating on cycle limit
- BUILD suite rc=0 wall=286 s
- BUILD_OK suite md5=8cc98fb37b17c4ed497cc84692f8327e prov=01deec5423b4ba0f7a44f230edb3c76eb6715bde
- Tiny gate deck as before (clean_face_hctau = false in <rad_m1>).

## Builds and gate
- BUILD rel6 rc=0 wall=291 s TUs=280
- BUILD_OK rel6 md5=7b0900a6f6e354806e236826bac99eb5 prov=01deec5423b4ba0f7a44f230edb3c76eb6715bde
- GATE PASS rel6: athena_rel6_01deec542 md5 7b0900a6 node hpc-sm-02-17 pass1=3.720e-11 log date 11-25 18:38:44 end: Terminating on cycle limit
- BUILD suite rc=0 wall=286 s
- BUILD_OK suite md5=8cc98fb37b17c4ed497cc84692f8327e prov=01deec5423b4ba0f7a44f230edb3c76eb6715bde
- gate deck: b1_1124/gate_relt.athinput (the old one without clean_face_hctau); the gate script now expects 3.720e-11

## Quick list, rel6 01deec542 against rel82 82231a337 (same node 02-03, one H200 each, side by side), and idort_shadow_gd
- 25 passed, 0 failed, 524 s
- 25 passed, 0 failed, 521 s
- 1 passed, 0 failed, 6 s
- rel6 quick: 25 lines, PASS 25, FAIL 0 [], SKIP 0; rel82 quick: 25 lines, PASS 25
- common 25, equal 25, different 0 []; only rel6 [], only rel82 []
- idort_shadow_gd (rel6): `PASS idort_shadow_gd       6 s  L1 0.0123 nray 162`
```
PASS conservation_sp       6 s  dE 3.33e-16 dM 8.09e-16 dLz 0 dLz_gr 5.92e-07
PASS conservation_cs       3 s  dE 1.2e-16 dM 1.3e-16 dLz 0 dLz_gr 6.48e-05
PASS restart_sp           10 s  files 4 maxdiff 0
PASS dt_retry             11 s  files 6 maxdiff 0 retries_fac1 1 retries_fac05 1
PASS layout_tile          45 s  maxdiff_1x2 0 maxrel_1x2 0 maxdiff_2x2 9.09e-13 maxrel_2x2 1.31e-10 maxdiff_dfl 0 dE_dfl 2.22e-16 dM_dfl 1.8e-16 dLz_dfl 1.96e-1
PASS layout_cs            13 s  warnings 0 maxdiff 0 dE 1.2e-16 dM 1.3e-16 dLz 0 maxdiff_dark 0 dE_dark 1.2e-16 dM_dark 1.3e-16 dLz_dark 0
PASS fofc_layout_cs       30 s  maxdiff_mb22 0 maxdiff_r2 0 maxdiff_r2mb22 0 fo_fofc 9.44e+04 maxdiff_restart 0 maxdiff_retry 0 dE 2.16e-16 dLz 0 dM 4.31e-16
PASS order_time           16 s  err_E [0.000113 3.14e-05] order_E [1.84] err_F [0.000586 0.000167] order_F [1.81] err_T [2.44e-05 6.89e-06] order_T [1.82]
PASS nlte_eps1e-1         57 s  L1 9.58e-05
PASS mhd_equil            11 s  dEtot 3.33e-16 dB 0 fatal 3 moved 0.0214 warnings 0
PASS mhd_layout           75 s  cs_r2 1.78e-15 cs_mb22 2e-15 cs_ranks 0 cs_restart 0 cs_dE 2.13e-16 cs_dLz 0 cs_divb 3.13e-15 cs_fofc 88947 sp_r2 2.78e-17 sp_mb
PASS wb_sp                29 s  vmax_over_cs 6.08e-09 vmax_over_cs_rot 6.1e-09
PASS eddington_cs         19 s  dE 0 dM 1.3e-16 dLz 0 dLz_gr 6.48e-05 warnings 0 restart_maxdiff 0 nogd 2.84e-14
PASS ray_clamp             5 s  fatal_fixed 1 fatal_ring 1
PASS momentum_box          5 s  dP_max 8.1e-17 dP_cyc [7.3e-17 8.04e-17 7.79e-17 8.02e-17 7.77e-17 7.92e-17 8.1e-17 7.96e-17] warnings 0
PASS wedge_sym            32 s  fphi [1.23e-09 2.48e-10 1.31e-09] warnings 0 wedge_fatal 1 fth 1.24e-16 theta_fatal 1 fth_rot 2.67e-10 rot_fatal 1
PASS large_step            4 s  KE_r 0.000809 rho_min 87.4 dEgas 0.000318
PASS floor_pay            28 s  dE 3.33e-16 dE0 1.67e-16 paid 0.153 unpaid 0 warnings 0 worst_lines 185 worst_bitwise True ekin_floor_dE 3.33e-16 ekin_floor_dM 
PASS floor_momentum       59 s  box_ekin [3.27e-16 2.33e-16 8.95e-17 2.06e-16 3.11e-17 9.74e-30 1.87e-16 1.9e-30] box_both [2.83e-16 2.73e-16 3.98e-16 2.94e-16 
PASS ceiling_metric       11 s  sp_dvphi 0 sp_dcap 0 sp_vtheta 0 sp_scale_ok 1 sp_capped 1408 cs_dvphi 6.31e-16 cs_dcap 6.31e-16 cs_vtheta 1.5e-16 cs_scale_ok 1
PASS cmdline_keys         11 s  same 0 note 1 typo 1 block 1
PASS idort_shadow          4 s  L1 0.0138
PASS idort_shadow_gl       4 s  L1 0.0126
PASS diffusion            30 s  errD [-0.000815 -7.79e-08 1.42e-10] L1 [0.00137 0.00117 0.00117] warnings 0 errD_tau1 -0.0817 L1_tau1 0.0118 errD_2d -9.29e-07 L
PASS advect_uniform        6 s  dT 5.38e-10 dE 2.68e-08 warnings 0
```
