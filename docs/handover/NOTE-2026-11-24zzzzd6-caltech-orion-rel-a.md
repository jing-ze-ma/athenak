# NOTE 2026-11-24zzzzd6 caltech -> orion: GPU gate of rel-a c880e2703 (your 24zzzzd6)

## 1. Builds (own worktree wt_orion1121/rela at c880e2703, CUDA 12.9 + g++ 13, HOPPER90, CPU steps on 02-03)
- BUILD rela rc=0 wall=718 s TUs=250
- BUILD_OK rela md5=5cfc3be547ad477553ec5db0e9c07fca prov=c880e2703605969324c5e1695fae47fbfb1509fc
- GATE PASS rela: athena_rela_c880e2703 md5 5cfc3be5 node hpc-sm-02-03 pass1=3.686e-11 log date 11-24 18:42:42 end: Terminating on cycle limit
- BUILD suite rc=0 wall=713 s
- BUILD_OK suite md5=51bcde3e842c725d5e0397982afda953 prov=c880e2703605969324c5e1695fae47fbfb1509fc
- star_vet build follows; the suite starts on the built_in binary.

## 1. Builds so far (11-24 18:50, own worktree wt_orion1121/rela at c880e2703, CUDA 12.9 + g++ 13, HOPPER90, CPU steps on 02-03)
- BUILD rela rc=0 wall=718 s TUs=250
- BUILD_OK rela md5=5cfc3be547ad477553ec5db0e9c07fca prov=c880e2703605969324c5e1695fae47fbfb1509fc
- GATE PASS rela: athena_rela_c880e2703 md5 5cfc3be5 node hpc-sm-02-03 pass1=3.686e-11 log date 11-24 18:42:42 end: Terminating on cycle limit
- BUILD suite rc=0 wall=713 s
- BUILD_OK suite md5=51bcde3e842c725d5e0397982afda953 prov=c880e2703605969324c5e1695fae47fbfb1509fc
- **star_vet: not built.** -D PROBLEM=star_vet fails at CMake; the branch has src/pgen/he_star_vet.cpp and no star_vet source or alias. Nothing edited. The error:
```
45:CMake Error at src/CMakeLists.txt:263 (target_sources):
46-  Cannot find source file:
47-    pgen/star_vet.cpp
48-  Tried extensions .c .C .c++ .cc .cpp .cxx .cu .mpp .m .M .mm .ixx .cppm
49-  .ccm .cxxm .c++m .h .hh .h++ .hm .hpp .hxx .in .txx .f .F .for .f77 .f90
50-  .f95 .f03 .hip .ispc
51--- Generating done (2.6s)
52-CMake Generate step failed.  Build files cannot be regenerated correctly.
-- Generating done (2.6s)
CMake Generate step failed.  Build files cannot be regenerated correctly.
```

## 3. Old-rst check, M35ZAMS (4 ranks, 01-17)
- wr8: 95813cf85 restarts the T2 stop rst (m35zams.stop1124b, c22155) with time2_gtime_slope = false and clean_stag_thinf = false (my reading of "keys false"; the other candidate keys as in production) for 8 cycles and writes its rst (c22163). base16: 95813cf85 from that rst for 16 cycles. rel8: rel-a from the same rst for 8 cycles, writing its own rst; rel8b: rel-a from that for 8 more. Compared by printed t.
- wr8: 2026-11-24 18:50:22 rc=0 (srun) ended: Terminating on cycle limit wall=35 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- base16: 2026-11-24 18:50:51 rc=0 (srun) ended: Terminating on cycle limit wall=29 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- rel8: 2026-11-24 18:51:30 rc=0 (srun) ended: Terminating on cycle limit wall=39 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- rel8b: 2026-11-24 18:51:57 rc=0 (srun) ended: Terminating on cycle limit wall=27 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- user.hst rel-a first 8: 8 rows, printed-identical to the base 16-cycle run at the same t: 8; first different or missing t: none
- user.hst rel-a second 8: 8 rows, printed-identical to the base 16-cycle run at the same t: 8; first different or missing t: none
- user.hst: base 16 rows; BITWISE (printed)
- hydro.hst rel-a first 8: 8 rows, printed-identical to the base 16-cycle run at the same t: 8; first different or missing t: none
- hydro.hst rel-a second 8: 8 rows, printed-identical to the base 16-cycle run at the same t: 8; first different or missing t: none
- hydro.hst: base 16 rows; BITWISE (printed)

## star_vet build
- BUILD star rc=3 wall=20 s
- 1:srun: error: hpc-sm-02-03: task 0: Exited with exit code 3
- BUILD_FAIL
- cmake: CMake Error at src/CMakeLists.txt:263 (target_sources):
- cmake:   .ccm .cxxm .c++m .h .hh .h++ .hm .hpp .hxx .in .txx .f .F .for .f77 .f90
- cmake:   .f95 .f03 .hip .ispc
- cmake: 
- cmake: 
- cmake: -- Generating done (2.6s)
- cmake: CMake Generate step failed.  Build files cannot be regenerated correctly.

## 2. GPU suite, full list (2 H200 of 02-03, jobs 1; the srun exit code was lost when the script was split, the last line of the suite is the verdict): 37 passed, 1 failed, 1401 s; test lines 38 (38 expected)
- CLEAN WARNING: 1 lines in 1 files of the kept run directories (38 directories).
```
PASS conservation_sp       6 s  dE 1.11e-16 dM 8.99e-16 dLz 3.91e-16 dLz_gr 6.94e-07
PASS conservation_cs       3 s  dE 1.2e-16 dM 2.61e-16 dLz 1.17e-16 dLz_gr 6.44e-05
PASS conservation_lz      28 s  sp_dE 2.22e-16 sp_dM 7.19e-16 sp_dLz 1.96e-16 sp_dLz_gr 3.91e-16 cs_dE 1.2e-16 cs_dM 2.61e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.24e-05 restart_maxdiff 0 retry_maxdiff 0 warni
PASS floorcap_spin        45 s  sp_dE 2.22e-16 sp_dM 5.4e-16 sp_dLz 1.96e-16 sp_dLz_gr 3.91e-16 n_v 2.59e+04 cs_dE 1.2e-16 cs_dM 1.3e-16 cs_dLz 1.17e-16 cs_dLz_gr 7.8e-05 nospin_sp_dLz_gr 0.288 restar
PASS floor_dt             27 s  dE 1.11e-16 dt_gain 0.873 restart_maxdiff 0 retry_maxdiff 0 warnings 0
PASS restart_sp           12 s  files 4 maxdiff 0
PASS restart_cs           10 s  files 4 maxdiff 0
PASS dt_retry             12 s  files 6 maxdiff 0 retries_fac1 1 retries_fac05 1
PASS layout_tile          39 s  maxdiff_1x2 0 maxrel_1x2 0 maxdiff_2x2 9.09e-13 maxrel_2x2 7.48e-11 maxdiff_dfl 0 dE_dfl 2.22e-16 dM_dfl 1.8e-16 dLz_dfl 1.96e-16 maxdiff_fofc 0 dE_fofc 2.22e-16 dM_fof
PASS order_time           20 s  err_E [0.00102 0.000344 0.000113 3.15e-05] order_E [1.57 1.61 1.84] err_F [0.00475 0.00168 0.000588 0.000168] order_F [1.5 1.52 1.81] err_T [0.000167 6.61e-05 2.46e-05 
PASS order_stiff          24 s  err_E [4.31e-06 1.13e-06 2.9e-07] order_E [1.94 1.96] err_v [0.00449 0.00118 0.000305] order_v [1.92 1.96]
PASS order_space         574 s  err_E [0.071 0.0203 0.00444] order_E [1.81 2.19] err_v [0.12 0.0422 0.011] order_v [1.5 1.94]
PASS nlte_eps1e-1         59 s  L1 9.58e-05
PASS nlte_eps1e-4         58 s  L1 0.00133
PASS thermal_eq           20 s  errT_h 0.000552 dEtot_h 1.26e-15 errT_h2 0.000124 dEtot_h2 1.26e-15 order 2.16
PASS wb_sp                25 s  vmax_over_cs 4.57e-09 vmax_over_cs_rot 4.67e-09
PASS wb_cs                14 s  vmax_over_cs 1.18e-09 vmax_over_cs_rot 6.93e-09
PASS momentum_box          4 s  dP_max 5.64e-17 dP_cyc [5.64e-17 5.27e-17 5.23e-17 5.42e-17 5.45e-17 5.21e-17 5.35e-17 5.25e-17] warnings 0
FAIL radwave             103 s  err_cph [0.00303 0.00189 0.00137 0.000364 0.00277 0.00137] err_damp [0.00235 0.000105 4.01e-05 0.000208 0.000545 4.02e-05] damp [0.711 0.164 0.0151 0.000302 0.164 0.015
PASS wedge_sym            29 s  fphi [3.11e-09 9.79e-12 3.11e-10] warnings 0 wedge_fatal 1 fth 6.91e-15 theta_fatal 1 fth_rot 1.56e-09 rot_fatal 1
PASS floor_pay            16 s  dE 2.22e-16 dE0 2.22e-16 paid 0.285 unpaid 0 warnings 0 worst_lines 195 worst_bitwise True ekin_dE 1.11e-16 ekin_dM 5.4e-16 ekin_dLz 3.91e-16
PASS ceiling_metric       10 s  sp_dvphi 0 sp_dcap 0 sp_vtheta 0 sp_scale_ok 1 sp_capped 1408 cs_dvphi 6.31e-16 cs_dcap 6.31e-16 cs_vtheta 1.5e-16 cs_scale_ok 1 cs_capped 2928
PASS inf_reject            4 s  dE 2.22e-16 nonfinite 0 rejected 1
PASS idort_beams           5 s  L1 0.177 lum_top 1.19 flux_err 5.97e-08
PASS idort_beams_mb       21 s  maxdiff_tile_x2 0 maxdiff_block_x1 7.61e-12
PASS idort_shadow          3 s  L1 0.0138
PASS idort_shadow_gl       3 s  L1 0.0126
PASS idort_cylinder        6 s  E_in 0.000531 E_out 0.0111 F_out 0.0116 FH_out 0.203
PASS idort_radshock_m2      9 s  L1_rho 0.000356 L1_v 0.000213 L1_Tg 0.000385 L1_Tr 0.00026
PASS idort_radshock_m5      8 s  L1_rho 0.0017 L1_v 0.000818 L1_Tg 0.00162 L1_Tr 0.00094
PASS idort_thermal_eq2     16 s  errT 0.00196 dEtot 1.88e-16
PASS diffusion            29 s  errD [-0.000815 -7.79e-08 1.42e-10] L1 [0.00137 0.00117 0.00117] warnings 0 errD_tau1 -0.0817 L1_tau1 0.0118 errD_2d -9.3e-07 L1_2d 0.00659 th_min 1.15e-17 th_max -0.03
PASS radshock             19 s  L1_rho 0.00438 L1_v 0.00358 L1_Tg 0.00659 L1_Tr 0.0107 Emin 1.72e+11 warnings 0
PASS radshock_m5          20 s  L1_rho 0.00391 L1_v 0.00212 L1_Tg 0.0062 L1_Tr 0.00744 Emin 1.79e+11 warnings 0
PASS marshak              12 s  Emin 8.92e-19 Tmin 9.96e-13 Emax 0.995 Tmax 0.999 Erise 8.89e-19 x_front 0.535 warnings 0
PASS advect_uniform        5 s  dT 2.73e-06 dE 2.35e-05 warnings 0
PASS advect_pulse         75 s  L1 [0.000873 0.000575 0.00082 0.000775] warnings 0
PASS shadow               28 s  umbra 0.00609 lit_err 0.00825 n_umbra 774 n_lit 3792 warnings 0
--
37 passed, 1 failed, 1401 s
srun: error: hpc-sm-02-03: task 0: Exited with exit code 1
```

### Values against the earlier GPU suites (printed digits). The pins are the bounds in run_tests.py: a test outside its pin is a FAIL line above; the lists below show every printed value that moved, inside the pin or not.
- rel-a: 38 lines, PASS 37, FAIL 1: ['radwave']

- against starname 4a26282a6 (01-15, 37/37): 37 common tests, printed values equal in 16, different in 21, only in rel-a ['order_stiff']
  - conservation_sp
    - rel-a: PASS dE 1.11e-16 dM 8.99e-16 dLz 3.91e-16 dLz_gr 6.94e-07
    - other: PASS dE 2.22e-16 dM 8.99e-16 dLz 1.96e-16 dLz_gr 6.94e-07
  - conservation_lz
    - rel-a: PASS sp_dE 2.22e-16 sp_dM 7.19e-16 sp_dLz 1.96e-16 sp_dLz_gr 3.91e-16 cs_dE 1.2e-16 cs_dM 2.61e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.24e-05 restart_maxdiff 0 retry_maxdiff 0 warnings 0
    - other: PASS sp_dE 2.22e-16 sp_dM 8.99e-16 sp_dLz 1.96e-16 sp_dLz_gr 3.91e-16 cs_dE 1.2e-16 cs_dM 2.61e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.24e-05 restart_maxdiff 0 retry_maxdiff 0 warnings 0
  - floorcap_spin
    - rel-a: PASS sp_dE 2.22e-16 sp_dM 5.4e-16 sp_dLz 1.96e-16 sp_dLz_gr 3.91e-16 n_v 2.59e+04 cs_dE 1.2e-16 cs_dM 1.3e-16 cs_dLz 1.17e-16 cs_dLz_gr 7.8e-05 nospin_sp_dLz_gr 0.288 restart_maxdiff 0 retry_maxdiff 0 layout_maxdiff 0 warnings 0
    - other: PASS sp_dE 2.22e-16 sp_dM 7.19e-16 sp_dLz 1.96e-16 sp_dLz_gr 3.91e-16 n_v 2.53e+04 cs_dE 1.2e-16 cs_dM 2.61e-16 cs_dLz 1.17e-16 cs_dLz_gr 7.8e-05 nospin_sp_dLz_gr 0.287 restart_maxdiff 0 retry_maxdiff 0 layout_maxdiff 0 warnings 0
  - floor_dt
    - rel-a: PASS dE 1.11e-16 dt_gain 0.873 restart_maxdiff 0 retry_maxdiff 0 warnings 0
    - other: PASS dE 2.22e-16 dt_gain 0.873 restart_maxdiff 0 retry_maxdiff 0 warnings 0
  - layout_tile
    - rel-a: PASS maxdiff_1x2 0 maxrel_1x2 0 maxdiff_2x2 9.09e-13 maxrel_2x2 7.48e-11 maxdiff_dfl 0 dE_dfl 2.22e-16 dM_dfl 1.8e-16 dLz_dfl 1.96e-16 maxdiff_fofc 0 dE_fofc 2.22e-16 dM_fofc 0.000635 dLz_fofc 1.96e-16 fo_fofc 4.3e+04
    - other: PASS maxdiff_1x2 0 maxrel_1x2 0 maxdiff_2x2 7.9e-12 maxrel_2x2 3.11e-10 maxdiff_dfl 0 dE_dfl 1.11e-16 dM_dfl 1.8e-16 dLz_dfl 1.96e-16 maxdiff_fofc 0 dE_fofc 1.11e-16 dM_fofc 0.000636 dLz_fofc 1.96e-16 fo_fofc 4.3e+04
  - order_time
    - rel-a: PASS err_E [0.00102 0.000344 0.000113 3.15e-05] order_E [1.57 1.61 1.84] err_F [0.00475 0.00168 0.000588 0.000168] order_F [1.5 1.52 1.81] err_T [0.000167 6.61e-05 2.46e-05 6.91e-06] order_T [1.34 1.43 1.83]
    - other: PASS err_E [0.00102 0.000344 0.000113 3.15e-05] order_E [1.57 1.61 1.84] err_F [0.00477 0.00169 0.000589 0.000168] order_F [1.5 1.52 1.81] err_T [0.000167 6.61e-05 2.45e-05 6.9e-06] order_T [1.34 1.43 1.83]
  - order_space
    - rel-a: PASS err_E [0.071 0.0203 0.00444] order_E [1.81 2.19] err_v [0.12 0.0422 0.011] order_v [1.5 1.94]
    - other: PASS err_E [0.071 0.0205 0.00449] order_E [1.79 2.19] err_v [0.12 0.0424 0.0111] order_v [1.5 1.94]
  - wb_cs
    - rel-a: PASS vmax_over_cs 1.18e-09 vmax_over_cs_rot 6.93e-09
    - other: PASS vmax_over_cs 1.18e-09 vmax_over_cs_rot 6.94e-09
  - momentum_box
    - rel-a: PASS dP_max 5.64e-17 dP_cyc [5.64e-17 5.27e-17 5.23e-17 5.42e-17 5.45e-17 5.21e-17 5.35e-17 5.25e-17] warnings 0
    - other: PASS dP_max 5.64e-17 dP_cyc [5.64e-17 5.38e-17 5.04e-17 5.17e-17 4.88e-17 5.11e-17 5.16e-17 5.13e-17] warnings 0
  - radwave
    - rel-a: FAIL err_cph [0.00303 0.00189 0.00137 0.000364 0.00277 0.00137] err_damp [0.00235 0.000105 4.01e-05 0.000208 0.000545 4.02e-05] damp [0.711 0.164 0.0151 0.000302 0.164 0.0151] warnings 0 dt_gain 2.77 dt_gain_phase [2.15 1]  <- err_cph [0.00303 0.00189 0.00137 0.000364 0.00277 0.00137] > 0.003
    - other: PASS err_cph [0.00199 0.00189 0.00137 0.000371 0.00273 0.00137] err_damp [0.000538 7.88e-05 0.000173 0.000684 0.000457 0.000173] damp [0.714 0.164 0.0152 0.000778 0.164 0.0152] warnings 0 dt_gain 2.77 dt_gain_phase [2.15 1]
  - wedge_sym
    - rel-a: PASS fphi [3.11e-09 9.79e-12 3.11e-10] warnings 0 wedge_fatal 1 fth 6.91e-15 theta_fatal 1 fth_rot 1.56e-09 rot_fatal 1
    - other: PASS fphi [1.22e-09 1.18e-10 8.25e-10] warnings 0 wedge_fatal 1 fth 8.06e-15 theta_fatal 1 fth_rot 1.24e-09 rot_fatal 1
  - floor_pay
    - rel-a: PASS dE 2.22e-16 dE0 2.22e-16 paid 0.285 unpaid 0 warnings 0 worst_lines 195 worst_bitwise True ekin_dE 1.11e-16 ekin_dM 5.4e-16 ekin_dLz 3.91e-16
    - other: PASS dE 2.22e-16 dE0 2.22e-16 paid 0.286 unpaid 0 warnings 0 worst_lines 196 worst_bitwise True ekin_dE 1.11e-16 ekin_dM 7.19e-16 ekin_dLz 1.96e-16
  - idort_beams
    - rel-a: PASS L1 0.177 lum_top 1.19 flux_err 5.97e-08
    - other: PASS L1 0.201 lum_top 1.2 flux_err 5.92e-08
  - idort_beams_mb
    - rel-a: PASS maxdiff_tile_x2 0 maxdiff_block_x1 7.61e-12
    - other: PASS maxdiff_tile_x2 0 maxdiff_block_x1 2.05e-12
  - idort_radshock_m2
    - rel-a: PASS L1_rho 0.000356 L1_v 0.000213 L1_Tg 0.000385 L1_Tr 0.00026
    - other: PASS L1_rho 0.000383 L1_v 0.000245 L1_Tg 0.000396 L1_Tr 0.000271
  - idort_radshock_m5
    - rel-a: PASS L1_rho 0.0017 L1_v 0.000818 L1_Tg 0.00162 L1_Tr 0.00094
    - other: PASS L1_rho 0.00216 L1_v 0.000919 L1_Tg 0.00181 L1_Tr 0.00103
  - radshock
    - rel-a: PASS L1_rho 0.00438 L1_v 0.00358 L1_Tg 0.00659 L1_Tr 0.0107 Emin 1.72e+11 warnings 0
    - other: PASS L1_rho 0.00441 L1_v 0.00361 L1_Tg 0.00656 L1_Tr 0.0107 Emin 1.72e+11 warnings 0
  - radshock_m5
    - rel-a: PASS L1_rho 0.00391 L1_v 0.00212 L1_Tg 0.0062 L1_Tr 0.00744 Emin 1.79e+11 warnings 0
    - other: PASS L1_rho 0.00454 L1_v 0.00224 L1_Tg 0.00625 L1_Tr 0.00746 Emin 1.79e+11 warnings 0
  - marshak
    - rel-a: PASS Emin 8.92e-19 Tmin 9.96e-13 Emax 0.995 Tmax 0.999 Erise 8.89e-19 x_front 0.535 warnings 0
    - other: PASS Emin 8.92e-19 Tmin 9.96e-13 Emax 0.995 Tmax 0.999 Erise 2.28e-18 x_front 0.535 warnings 0
  - advect_uniform
    - rel-a: PASS dT 2.73e-06 dE 2.35e-05 warnings 0
    - other: PASS dT 1.45e-06 dE 1.26e-05 warnings 0
  - advect_pulse
    - rel-a: PASS L1 [0.000873 0.000575 0.00082 0.000775] warnings 0
    - other: PASS L1 [0.000873 0.000555 0.000809 0.000723] warnings 0

- against sdorder5 10bf5513c (02-09): 38 common tests, printed values equal in 35, different in 3, only in rel-a []
  - wb_cs
    - rel-a: PASS vmax_over_cs 1.18e-09 vmax_over_cs_rot 6.93e-09
    - other: PASS vmax_over_cs 1.18e-09 vmax_over_cs_rot 6.94e-09
  - idort_beams
    - rel-a: PASS L1 0.177 lum_top 1.19 flux_err 5.97e-08
    - other: PASS L1 0.777 lum_top 1.79 flux_err 2.66e-07
  - idort_beams_mb
    - rel-a: PASS maxdiff_tile_x2 0 maxdiff_block_x1 7.61e-12
    - other: PASS maxdiff_tile_x2 0 maxdiff_block_x1 4.32e-12
