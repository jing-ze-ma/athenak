# NOTE 2026-11-26t caltech -> orion: CUDA build and suite of caltech-next5-b5ade8753 (your 26t)

## Builds (worktree wt_orion1121/next5; CUDA 12.9 + g++13; 02-03, 14 cores each)
- BUILD next5 rc=0 wall=728 s TUs=250
- BUILD_OK next5 md5=706572b0a64184be8d0cd61ce95132dd prov=b5ade8753a53223605b2d52c4a47d62479066c1a
- GATE PASS next5: athena_next5_b5ade8753 md5 706572b0 node hpc-sm-02-03 pass1=3.686e-11 log date 11-25 13:06:13 end: Terminating on cycle limit
- BUILD suite rc=0 wall=737 s
- BUILD_OK suite md5=a185d1011a283438546e0b8a3aa0b319 prov=b5ade8753a53223605b2d52c4a47d62479066c1a
- Tiny gate deck as before (clean_face_hctau = false in <rad_m1>).

## 3. GPU suite, full list in 4 disjoint --only lists, one H200 of the node named in the log each (a 2-rank test shares its one GPU), jobs 1, --star-exe; srun rc of the 4: 0 0 0 0): 1 passed, 0 failed, 986 s;1 passed, 0 failed, 120 s;28 passed, 0 failed, 1020 s;28 passed, 0 failed, 1011 s;1 passed, 0 failed, 19 s;; test lines 59, SKIP lines 0 (59 tests of this branch)
- CLEAN WARNING: 1 lines in the kept run directories, by test:       1 inf_reject; (inf_reject injects one on purpose).
- star_smoke (dev against the CPU reference, bound 2e-4, nothing changed): `PASS star_smoke            3 s  dev 0.00011 warnings 0 worst L_intE`
```
PASS order_space         986 s  err_E [0.0907 0.0208 0.00434] order_E [2.12 2.26] err_v [0.106 0.0379 0.00979] order_v [1.48 1.95]
PASS star_cs_ring        120 s  maxdiff_r3 0 warnings 0 arena 1 fused 1
PASS conservation_sp      11 s  dE 3.33e-16 dM 8.09e-16 dLz 0 dLz_gr 5.92e-07
PASS conservation_cs       5 s  dE 1.2e-16 dM 1.3e-16 dLz 0 dLz_gr 6.48e-05
PASS floor_dt             32 s  dE 3.33e-16 dt_gain 0.791 restart_maxdiff 0 retry_maxdiff 0 warnings 0
PASS restart_fixes        88 s  warnings 0 notes 6 stored 6 maxdiff_lat_wb 0 moved_lat_wb 0.00177 maxdiff_fofc 0 moved_fofc 0.0296 maxdiff_merge 0 moved_merge 0.00787 maxdiff_blockend 0 moved_blockend
PASS restart_sp           12 s  files 4 maxdiff 0
PASS dt_retry             11 s  files 6 maxdiff 0 retries_fac1 1 retries_fac05 1
PASS layout_tile          45 s  maxdiff_1x2 0 maxrel_1x2 0 maxdiff_2x2 9.09e-13 maxrel_2x2 1.31e-10 maxdiff_dfl 0 dE_dfl 2.22e-16 dM_dfl 1.8e-16 dLz_dfl 1.96e-16 maxdiff_fofc 0 dE_fofc 1.11e-16 dM_fof
PASS order_time           24 s  err_E [0.00102 0.000342 0.000113 3.14e-05] order_E [1.57 1.6 1.84] err_F [0.00475 0.00168 0.000586 0.000167] order_F [1.5 1.52 1.81] err_T [0.000166 6.56e-05 2.44e-05 6
PASS nlte_eps1e-4         58 s  L1 0.00133
PASS mhd_layout           69 s  cs_r2 1.78e-15 cs_mb22 2e-15 cs_ranks 0 cs_restart 0 cs_dE 2.13e-16 cs_dLz 0 cs_divb 3.13e-15 cs_fofc 88947 sp_r2 2.78e-17 sp_mb22 2.78e-17 sp_ranks 0 sp_restart 0 sp_d
PASS mhd_radwave         204 s  err_cph [0.000754 0.000658 0.000224 0.000141 0.000293 0.000474 0.000905 0.000836 0.000684 0.000536 0.001 0.00145 0.00126 0.000342] err_damp [0.000288 0.000294 0.000219 
PASS mhd_wb               90 s  static [6.68e-09 5.52e-10 1.68e-08] field [6.77e-09 5.87e-10 2.62e-08] geneos [6.75e-09 1.68e-08] dyn [2.6e-08 2.97e-08 2.6e-08] warnings 0 monopole 4.37e-06 restart_ma
PASS mhd_vet              93 s  warnings 0 cs_layout 4.52e-15 cs_ranks 0 cs_retry 0 cs_divb 1.02e-14 cs_restart 0 cs_dM 2.61e-16 cs_dLz 0 cs_dE 1.11e-16 cs_moved 0.443 sp_layout 0 sp_ranks 0 sp_retry 
PASS wb_cs                15 s  vmax_over_cs 5.11e-10 vmax_over_cs_rot 7.4e-09
PASS eddington_cs         17 s  dE 0 dM 1.3e-16 dLz 0 dLz_gr 6.48e-05 warnings 0 restart_maxdiff 0 nogd 2.84e-14
PASS ray_clamp             5 s  fatal_fixed 1 fatal_ring 1
PASS floor_pay            28 s  dE 3.33e-16 dE0 1.67e-16 paid 0.153 unpaid 0 warnings 0 worst_lines 185 worst_bitwise True ekin_floor_dE 3.33e-16 ekin_floor_dM 1.08e-15 ekin_cart_dM 2.22e-16 ekin_dE 3
PASS floor_momentum       61 s  box_ekin [3.27e-16 2.33e-16 8.95e-17 2.06e-16 3.11e-17 9.74e-30 1.87e-16 1.9e-30] box_both [2.83e-16 2.73e-16 3.98e-16 2.94e-16 2.57e-16 2.15e-16 4.28e-17 2.14e-17] box
PASS ceiling_metric       10 s  sp_dvphi 0 sp_dcap 0 sp_vtheta 0 sp_scale_ok 1 sp_capped 1408 cs_dvphi 6.31e-16 cs_dcap 6.31e-16 cs_vtheta 1.5e-16 cs_scale_ok 1 cs_capped 2928
PASS inf_reject            3 s  dE 3.33e-16 nonfinite 0 rejected 1
PASS wb_star              11 s  warnings 0 prof_in_vr 0.00226 prof_in_lat 4.57e-09 rigid_in_vr 0.00225 rigid_in_lat 0.000213 eq_in_vr 9.41e-05 eq_in_lat 2.13e-09
PASS idort_beams_mb       23 s  maxdiff_tile_x2 0 maxdiff_block_x1 1.13e-10
PASS idort_shadow          3 s  L1 0.0138
PASS idort_cylinder        8 s  E_in 0.000533 E_out 0.0115 F_out 0.0114 FH_out 0.204
PASS idort_radshock_m2     11 s  L1_rho 0.000356 L1_v 0.000213 L1_Tg 0.000385 L1_Tr 0.00026
PASS diffusion            30 s  errD [-0.000815 -7.79e-08 1.42e-10] L1 [0.00137 0.00117 0.00117] warnings 0 errD_tau1 -0.0817 L1_tau1 0.0118 errD_2d -9.29e-07 L1_2d 0.00659 th_min 1.15e-17 th_max -0.0
PASS radshock_m5          25 s  L1_rho 0.00386 L1_v 0.00212 L1_Tg 0.0062 L1_Tr 0.00744 Emin 1.79e+11 warnings 0
PASS shadow               27 s  umbra 0.00609 lit_err 0.00825 n_umbra 774 n_lit 3792 warnings 0
PASS conservation_lz      37 s  sp_dE 3.33e-16 sp_dM 9.44e-16 sp_dLz 0 sp_dLz_gr 1.47e-16 cs_dE 1.2e-16 cs_dM 2.61e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.8e-05 restart_maxdiff 0 retry_maxdiff 0 warnings 0
PASS floorcap_spin        52 s  sp_dE 3.33e-16 sp_dM 8.09e-16 sp_dLz 0 sp_dLz_gr 2.94e-16 n_v 6.03e+04 cs_dE 1.2e-16 cs_dM 2.61e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.3e-05 nospin_sp_dLz_gr 0.366 restart_max
PASS merge_lat            63 s  sp_dM 4.88e-16 sp_dLz 9.88e-16 sp_chk 1.09e-15 sp_amin 0 sp_dE 2.25e-16 sp_vlat 8.31e-16 sp_rigid_dLz 5.01e-05 cs_dM 9.41e-16 cs_dLz 6.49e-16 cs_chk 1.39e-15 cs_amin 0 
PASS restart_cs           11 s  files 4 maxdiff 0
PASS layout_cs            14 s  warnings 0 maxdiff 0 dE 1.2e-16 dM 1.3e-16 dLz 0 maxdiff_dark 0 dE_dark 1.2e-16 dM_dark 1.3e-16 dLz_dark 0
PASS fofc_layout_cs       29 s  maxdiff_mb22 0 maxdiff_r2 0 maxdiff_r2mb22 0 fo_fofc 9.44e+04 maxdiff_restart 0 maxdiff_retry 0 dE 2.16e-16 dLz 0 dM 4.31e-16
PASS order_stiff          26 s  err_E [4.31e-06 1.13e-06 2.9e-07] order_E [1.94 1.96] err_v [0.00449 0.00118 0.000305] order_v [1.92 1.96]
PASS nlte_eps1e-1         57 s  L1 9.58e-05
PASS thermal_eq           22 s  errT_h 0.000552 dEtot_h 1.12e-15 errT_h2 0.000124 dEtot_h2 7e-16 order 2.16
PASS mhd_equil            10 s  dEtot 3.33e-16 dB 0 fatal 3 moved 0.0214 warnings 0
PASS mhd_geneos          193 s  warnings 0 cs_layout 7.51e-14 cs_ranks 0 cs_retry 0 cs_divb 1.22e-14 cs_restart 0 cs_dM 2.61e-16 cs_dLz 0 cs_dE 1.11e-16 cs_moved 0.00409 sp_layout 0 sp_ranks 0 sp_retr
PASS wb_sp                29 s  vmax_over_cs 6.08e-09 vmax_over_cs_rot 6.1e-09
PASS shell_force_cs        3 s  lat 5.2e-05 rad 5.56e-05
PASS ray_cut_slab         78 s  warnings 0 maxdiff_nofuse 0 maxdiff_low 0 maxdiff_retry 0 shared [low retry shared] said_low 1 said_shared 0 fallbacks [0 1] moves 1
PASS momentum_box          5 s  dP_max 8.1e-17 dP_cyc [7.3e-17 8.04e-17 7.79e-17 8.02e-17 7.77e-17 7.92e-17 8.1e-17 7.96e-17] warnings 0
PASS radwave             143 s  err_cph [0.00149 0.00147 0.000928 0.000239 0.00226 0.000931] err_damp [0.000865 0.000225 1.63e-05 7.97e-05 0.000655 1.63e-05] damp [0.712 0.164 0.0151 0.000174 0.164 0.
PASS wedge_sym            32 s  fphi [1.23e-09 2.48e-10 1.31e-09] warnings 0 wedge_fatal 1 fth 1.24e-16 theta_fatal 1 fth_rot 2.67e-10 rot_fatal 1
PASS cmdline_keys         12 s  same 0 note 1 typo 1 block 1
PASS star_smoke            3 s  dev 0.00011 warnings 0 worst L_intE
PASS mlt_time             21 s  warnings 0 err [0.008 0.00109 0.000155] order 2.81 err_start [0.0714 0.0369 0.0188] order_start 0.976
PASS idort_beams          12 s  L1 0.177 lum_top 1.19 flux_err 5.97e-08 F_nocone 0.683
PASS idort_shadow_gl       4 s  L1 0.0126
PASS idort_radshock_m5     11 s  L1_rho 0.00171 L1_v 0.000822 L1_Tg 0.00161 L1_Tr 0.000938
PASS idort_thermal_eq2     15 s  errT 0.00196 dEtot 1.88e-16
PASS radshock             24 s  L1_rho 0.00438 L1_v 0.00358 L1_Tg 0.00659 L1_Tr 0.0107 Emin 1.72e+11 warnings 0
PASS marshak              13 s  Emin 8.92e-19 Tmin 9.96e-13 Emax 0.995 Tmax 0.999 Erise 3.07e-18 x_front 0.535 warnings 0
PASS advect_uniform        6 s  dT 5.38e-10 dE 2.68e-08 warnings 0
PASS advect_pulse         86 s  L1 [0.000872 0.000574 0.000821 0.000776] warnings 0
PASS mpi_ranks            19 s  warnings 0 fcone_lz 7 files 15 maxdiff_sp 0 maxdiff_cs 0 maxdiff_pay 0
```
### Printed values against next3 c3f133b44 (26f suite, common tests)
- next5: 59 lines, PASS 59, FAIL 0: []

- against next3 c3f133b44: 58 common tests, printed values equal in 51, different in 7, only in rel-e ['cmdline_keys']
  - order_time
    - next: PASS err_E [0.00102 0.000342 0.000113 3.14e-05] order_E [1.57 1.6 1.84] err_F [0.00475 0.00168 0.000586 0.000167] order_F [1.5 1.52 1.81] err_T [0.000166 6.56e-05 2.44e-05 6.89e-06] order_T [1.34 1.43 1.82]
    - other: PASS err_E [0.00102 0.000344 0.000113 3.15e-05] order_E [1.57 1.61 1.84] err_F [0.00475 0.00168 0.000588 0.000168] order_F [1.5 1.52 1.81] err_T [0.000167 6.61e-05 2.46e-05 6.91e-06] order_T [1.34 1.43 1.83]
  - conservation_lz
    - next: PASS sp_dE 3.33e-16 sp_dM 9.44e-16 sp_dLz 0 sp_dLz_gr 1.47e-16 cs_dE 1.2e-16 cs_dM 2.61e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.8e-05 restart_maxdiff 0 retry_maxdiff 0 warnings 0
    - other: PASS sp_dE 3.33e-16 sp_dM 9.44e-16 sp_dLz 0 sp_dLz_gr 1.47e-16 cs_dE 1.2e-16 cs_dM 2.61e-16 cs_dLz 0 cs_dLz_gr 6.54e-05 restart_maxdiff 0 retry_maxdiff 0 warnings 0
  - floorcap_spin
    - next: PASS sp_dE 3.33e-16 sp_dM 8.09e-16 sp_dLz 0 sp_dLz_gr 2.94e-16 n_v 6.03e+04 cs_dE 1.2e-16 cs_dM 2.61e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.3e-05 nospin_sp_dLz_gr 0.366 restart_maxdiff 0 retry_maxdiff 0 layout_maxdiff 0 warnings 0
    - other: PASS sp_dE 3.33e-16 sp_dM 8.09e-16 sp_dLz 1.47e-16 sp_dLz_gr 1.47e-16 n_v 5.85e+04 cs_dE 1.2e-16 cs_dM 1.3e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.3e-05 nospin_sp_dLz_gr 0.366 restart_maxdiff 0 retry_maxdiff 0 layout_maxdiff 0 warnings 0
  - mhd_equil
    - next: PASS dEtot 3.33e-16 dB 0 fatal 3 moved 0.0214 warnings 0
    - other: PASS dEtot 2.22e-16 dB 0 fatal 3 moved 0.0214 warnings 0
  - wedge_sym
    - next: PASS fphi [1.23e-09 2.48e-10 1.31e-09] warnings 0 wedge_fatal 1 fth 1.24e-16 theta_fatal 1 fth_rot 2.67e-10 rot_fatal 1
    - other: PASS fphi [1.23e-09 2.48e-10 1.31e-09] warnings 0 wedge_fatal 1 fth 1.24e-16 theta_fatal 1 fth_rot 1.56e-09 rot_fatal 1
  - star_smoke
    - next: PASS dev 0.00011 warnings 0 worst L_intE
    - other: PASS dev 4.31e-06 warnings 0 worst L_bot
  - mpi_ranks
    - next: PASS warnings 0 fcone_lz 7 files 15 maxdiff_sp 0 maxdiff_cs 0 maxdiff_pay 0
    - other: PASS warnings 0 fcone_lz 8 files 15 maxdiff_sp 0 maxdiff_cs 0 maxdiff_pay 0
