# NOTE 2026-11-24zzzzg8 caltech -> orion: GPU gate of rel-e 08c48ac10 (your 24zzzzg7 / g8)

## 1. Builds (own worktree wt_orion1121/rele at 08c48ac10, CUDA 12.9 + g++ 13, HOPPER90, CPU steps on 02-03 and 02-09 next to the untimed bsg arms): star_vet (PROBLEM = star_vet) and built_in
- BUILD rele rc=0 wall=405 s TUs=250
- BUILD_OK rele md5=8737e5b0b52061a57b49680e3f873d70 prov=08c48ac101b98ee287577ad7a6a6f5959aa63537
- GATE PASS rele: athena_rele_08c48ac10 md5 8737e5b0 node hpc-sm-02-09 pass1=3.686e-11 log date 11-24 20:49:53 end: Terminating on cycle limit
- BUILD suite rc=0 wall=406 s
- BUILD_OK suite md5=7732727852a6d57fabb0dfb9fff4830c prov=08c48ac101b98ee287577ad7a6a6f5959aa63537
- The tiny gate deck copy has one added line, clean_face_hctau = false (as for rel-b).

## 3. GPU suite, full list in 4 disjoint --only lists, one H200 of 02-03 each (a 2-rank test shares its one GPU), jobs 1, --star-exe; srun rc of the 4: 0 0 1 0): 1 passed, 0 failed, 976 s;14 passed, 0 failed, 387 s;15 passed, 1 failed, 418 s;15 passed, 0 failed, 407 s;; test lines 46, SKIP lines 0 (46 tests of this branch; the count expected is in its docs)
- CLEAN WARNING: 1 lines in the kept run directories, by test:       1 inf_reject; (inf_reject injects one on purpose).
- star_smoke (dev against the CPU reference, bound 2e-4, nothing changed): `PASS star_smoke            3 s  dev 3.69e-06 warnings 0 worst L_bot`
```
PASS order_space         976 s  err_E [0.0907 0.0208 0.00434] order_E [2.12 2.26] err_v [0.106 0.0379 0.00979] order_v [1.48 1.95]
PASS conservation_cs       6 s  dE 1.2e-16 dM 1.3e-16 dLz 1.17e-16 dLz_gr 6.48e-05
PASS conservation_lz      35 s  sp_dE 3.33e-16 sp_dM 9.44e-16 sp_dLz 2.94e-16 sp_dLz_gr 1.47e-16 cs_dE 1.2e-16 cs_dM 2.61e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.54e-05 restart_maxdiff 0 retry_maxdiff 0 warni
PASS floor_dt             34 s  dE 3.33e-16 dt_gain 0.791 restart_maxdiff 0 retry_maxdiff 0 warnings 0
PASS restart_sp           13 s  files 4 maxdiff 0
PASS fofc_layout_cs       35 s  maxdiff_mb22 0 maxdiff_r2 0 maxdiff_r2mb22 0 fo_fofc 9.44e+04 maxdiff_restart 0 maxdiff_retry 0 dE 2.16e-16 dLz 0 dM 4.31e-16
PASS order_stiff          30 s  err_E [4.31e-06 1.13e-06 2.9e-07] order_E [1.94 1.96] err_v [0.00449 0.00118 0.000305] order_v [1.92 1.96]
PASS thermal_eq           21 s  errT_h 0.000552 dEtot_h 1.12e-15 errT_h2 0.000124 dEtot_h2 7e-16 order 2.16
PASS radwave             143 s  err_cph [0.00149 0.00147 0.000928 0.000239 0.00226 0.000931] err_damp [0.000865 0.000225 1.63e-05 7.97e-05 0.000655 1.63e-05] damp [0.712 0.164 0.0151 0.000174 0.164 0.
PASS floor_pay            13 s  dE 3.33e-16 dE0 1.67e-16 paid 0.153 unpaid 0 warnings 0 worst_lines 185 worst_bitwise True ekin_dE 3.33e-16 ekin_dM 5.4e-16 ekin_dLz 1.47e-16
PASS star_smoke            3 s  dev 3.69e-06 warnings 0 worst L_bot
PASS idort_radshock_m2     10 s  L1_rho 0.000356 L1_v 0.000213 L1_Tg 0.000385 L1_Tr 0.00026
PASS idort_radshock_m5     11 s  L1_rho 0.00171 L1_v 0.000822 L1_Tg 0.00161 L1_Tr 0.000938
PASS radshock_m5          26 s  L1_rho 0.00386 L1_v 0.00212 L1_Tg 0.0062 L1_Tr 0.00744 Emin 1.79e+11 warnings 0
PASS advect_uniform        6 s  dT 5.38e-10 dE 2.68e-08 warnings 0
PASS floorcap_spin        50 s  sp_dE 3.33e-16 sp_dM 8.09e-16 sp_dLz 1.47e-16 sp_dLz_gr 1.47e-16 n_v 5.85e+04 cs_dE 1.2e-16 cs_dM 1.3e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.3e-05 nospin_sp_dLz_gr 0.366 resta
FAIL mpi_ranks             5 s    <- ERROR sp2: return code 1 (see /resnick/groups/carnegie_poc/jingze/b1_1124/runs_suitere2/mpi_ranks/sp2/run.log)
PASS restart_fixes        69 s  warnings 0 notes 5 stored 5 maxdiff_lat_wb 0 moved_lat_wb 0.00092 maxdiff_fofc 0 moved_fofc 0.00448 maxdiff_merge 0 moved_merge 0.00787 maxdiff_gtime 0 moved_gtime 6.85
PASS restart_cs           12 s  files 4 maxdiff 0
PASS dt_retry             13 s  files 6 maxdiff 0 retries_fac1 1 retries_fac05 1
PASS momentum_box          4 s  dP_max 8.1e-17 dP_cyc [7.3e-17 8.04e-17 7.79e-17 8.02e-17 7.77e-17 7.92e-17 8.1e-17 7.96e-17] warnings 0
PASS floor_momentum       55 s  box_ekin [3.27e-16 2.33e-16 8.95e-17 2.06e-16 3.11e-17 9.74e-30 1.87e-16 1.9e-30] box_both [2.83e-16 2.73e-16 3.98e-16 2.94e-16 2.57e-16 2.15e-16 4.28e-17 2.14e-17] box
PASS inf_reject            3 s  dE 3.33e-16 nonfinite 0 rejected 1
PASS wb_star               5 s  v_in 0.00309 v_all 0.672 lat_mid 1.31e-11 warnings 0
PASS idort_shadow          3 s  L1 0.0138
PASS idort_thermal_eq2     15 s  errT 0.00196 dEtot 1.88e-16
PASS diffusion            31 s  errD [-0.000815 -7.79e-08 1.42e-10] L1 [0.00137 0.00117 0.00117] warnings 0 errD_tau1 -0.0817 L1_tau1 0.0118 errD_2d -9.29e-07 L1_2d 0.00659 th_min 1.15e-17 th_max -0.0
PASS radshock             25 s  L1_rho 0.00438 L1_v 0.00358 L1_Tg 0.00659 L1_Tr 0.0107 Emin 1.72e+11 warnings 0
PASS marshak              12 s  Emin 8.92e-19 Tmin 9.96e-13 Emax 0.995 Tmax 0.999 Erise 3.07e-18 x_front 0.535 warnings 0
PASS advect_pulse         87 s  L1 [0.000872 0.000574 0.000821 0.000776] warnings 0
PASS shadow               28 s  umbra 0.00609 lit_err 0.00825 n_umbra 774 n_lit 3792 warnings 0
PASS conservation_sp       8 s  dE 3.33e-16 dM 8.09e-16 dLz 1.47e-16 dLz_gr 5.92e-07
PASS merge_lat            74 s  sp_dM 4.88e-16 sp_dLz 9.88e-16 sp_chk 1.09e-15 sp_amin 0 sp_dE 2.25e-16 sp_vlat 8.31e-16 sp_rigid_dLz 5.01e-05 cs_dM 9.41e-16 cs_dLz 6.49e-16 cs_chk 1.39e-15 cs_amin 0 
PASS layout_tile          53 s  maxdiff_1x2 0 maxrel_1x2 0 maxdiff_2x2 9.09e-13 maxrel_2x2 1.31e-10 maxdiff_dfl 0 dE_dfl 2.22e-16 dM_dfl 1.8e-16 dLz_dfl 3.92e-16 maxdiff_fofc 0 dE_fofc 1.11e-16 dM_fof
PASS order_time           24 s  err_E [0.00102 0.000344 0.000113 3.15e-05] order_E [1.57 1.61 1.84] err_F [0.00475 0.00168 0.000588 0.000168] order_F [1.5 1.52 1.81] err_T [0.000167 6.61e-05 2.46e-05 
PASS nlte_eps1e-1         59 s  L1 9.58e-05
PASS nlte_eps1e-4         57 s  L1 0.00133
PASS wb_sp                29 s  vmax_over_cs 6.08e-09 vmax_over_cs_rot 6.1e-09
PASS wb_cs                17 s  vmax_over_cs 5.11e-10 vmax_over_cs_rot 7.4e-09
PASS shell_force_cs        4 s  lat 5.2e-05 rad 5.56e-05
PASS wedge_sym            32 s  fphi [1.23e-09 2.48e-10 1.31e-09] warnings 0 wedge_fatal 1 fth 1.24e-16 theta_fatal 1 fth_rot 1.56e-09 rot_fatal 1
PASS ceiling_metric       10 s  sp_dvphi 0 sp_dcap 0 sp_vtheta 0 sp_scale_ok 1 sp_capped 1408 cs_dvphi 6.31e-16 cs_dcap 6.31e-16 cs_vtheta 1.5e-16 cs_scale_ok 1 cs_capped 2928
PASS idort_beams           6 s  L1 0.177 lum_top 1.19 flux_err 5.97e-08
PASS idort_beams_mb       23 s  maxdiff_tile_x2 0 maxdiff_block_x1 7.61e-12
PASS idort_shadow_gl       4 s  L1 0.0126
PASS idort_cylinder        8 s  E_in 0.000533 E_out 0.0115 F_out 0.0114 FH_out 0.204
```
### Printed values against rel-b (this suite's common tests)
- rel-e: 46 lines, PASS 45, FAIL 1: ['mpi_ranks']

- against rel-b c8a0fc099 (02-03): 40 common tests, printed values equal in 35, different in 5, only in rel-e ['fofc_layout_cs', 'mpi_ranks', 'restart_fixes', 'floor_momentum', 'merge_lat', 'shell_force_cs']
  - conservation_cs
    - rel-e: PASS dE 1.2e-16 dM 1.3e-16 dLz 1.17e-16 dLz_gr 6.48e-05
    - other: PASS dE 1.2e-16 dM 2.61e-16 dLz 1.17e-16 dLz_gr 6.48e-05
  - conservation_lz
    - rel-e: PASS sp_dE 3.33e-16 sp_dM 9.44e-16 sp_dLz 2.94e-16 sp_dLz_gr 1.47e-16 cs_dE 1.2e-16 cs_dM 2.61e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.54e-05 restart_maxdiff 0 retry_maxdiff 0 warnings 0
    - other: PASS sp_dE 3.33e-16 sp_dM 9.44e-16 sp_dLz 2.94e-16 sp_dLz_gr 1.47e-16 cs_dE 2.4e-16 cs_dM 2.61e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.54e-05 restart_maxdiff 0 retry_maxdiff 0 warnings 0
  - star_smoke
    - rel-e: PASS dev 3.69e-06 warnings 0 worst L_bot
    - other: PASS dev 4.1e-06 warnings 0 worst L_bot
  - wb_star
    - rel-e: PASS v_in 0.00309 v_all 0.672 lat_mid 1.31e-11 warnings 0
    - other: PASS v_in 0.00309 v_all 0.672 lat_mid 1.33e-11 warnings 0
  - wb_cs
    - rel-e: PASS vmax_over_cs 5.11e-10 vmax_over_cs_rot 7.4e-09
    - other: PASS vmax_over_cs 1.36e-09 vmax_over_cs_rot 6.93e-09

## mpi_ranks FAIL, first reading (11-24 21:11)
- 45 of 46 pass. mpi_ranks fails in 5 s with return code 1 at its 2-rank run (sp2), before athena starts:
```
--------------------------------------------------------------------------
An ORTE daemon has unexpectedly failed after launch and before
communicating back to mpirun. This could be caused by a number
of factors, including an inability to create a connection back
to mpirun due to a lack of common network interfaces and/or no
route found between them. Please check network connectivity
(including firewalls and network routing requirements).
--------------------------------------------------------------------------
```
- That is mpirun trying to start a daemon from inside my srun step (Open MPI's launcher under Slurm, 3-node holder); athena never ran, so the test has **no result yet** (it is not a code verdict either way). Nothing edited. I rerun this one test now on 2 H200 of 02-03 with the Slurm variables unset in the step's shell (a local mpirun).

## mpi_ranks rerun (11-24 21:13): PASS mpi_ranks            19 s  warnings 0 fcone_lz 8 files 15 maxdiff_sp 0 maxdiff_cs 0 maxdiff_pay 0
- 02-03, 2 H200 (one per rank), local mpirun (Slurm variables unset, --oversubscribe, plm isolated), same binary 77327278 built_in / 8737e5b0 star_vet, nothing edited; 1 passed, 0 failed, 19 s. **With it the rel-e 08c48ac10 GPU suite is 46 of 46** (45 from the 4-list run, mpi_ranks from this rerun).
