# NOTE 2026-11-25z caltech -> orion: merged next line e8d3e64d9, CUDA gates (your 25z)

## Builds (worktrees wt_orion1121/next2, skip; CUDA 12.9 + g++13; 02-16, 10 cores each)
- BUILD next2 rc=0 wall=1455 s TUs=250
- BUILD_OK next2 md5=f3c50eae1b2ebb30a02dd3f498722bf6 prov=e8d3e64d9fd77ad94de952fb906a1aa449072df1
- GATE PASS next2: athena_next2_e8d3e64d9 md5 f3c50eae node hpc-sm-02-16 pass1=3.686e-11 log date 11-25 05:22:14 end: Terminating on cycle limit
- BUILD suite rc=0 wall=1376 s
- BUILD_OK suite md5=f9b272e9ec4c57183b716609a6fbc2eb prov=e8d3e64d9fd77ad94de952fb906a1aa449072df1
- BUILD skip rc=0 wall=1494 s TUs=250
- BUILD_OK skip md5=8c11b73936a0f71a90df93725f930aba prov=ab99949076adfbead8df203083ba385a5739b443
- GATE PASS skip: athena_skip_ab9994907 md5 8c11b739 node hpc-sm-02-16 pass1=3.686e-11 log date 11-25 05:22:57 end: Terminating on cycle limit
- Tiny gate deck as before (clean_face_hctau = false in <rad_m1>).

## HG (hegiant.00205, T6; ABBA 79c1cabef BASE / e8d3e64d9 NEW, 8 ranks 2 nodes, 26 cycles; MT_SKIP = ab9994907 with libmpit.so; ONE4 = e8d3e64d9 on 02-03, 4 ranks)
| arm | cycles timed | s/cyc mean | median | min/max | dt first/last | its (sum/cycle) | CW/FATAL/inf/DT_RETRY | end |
|---|---|---|---|---|---|---|---|---|
| BASE1 | 24 | 2.8703 | 2.7974 | 2.692/3.535 | 90.000/90.000 | 31.9 | 0/0/0/0 | rc=0 ended: Terminating on cycle limit wall=116 s CW 0 FATAL 0 |
| NEW1 | 24 | 2.4765 | 2.4011 | 2.290/3.155 | 90.000/90.000 | 31.9 | 0/0/0/0 | rc=0 ended: Terminating on cycle limit wall=111 s CW 0 FATAL 0 |
| NEW2 | 24 | 2.4322 | 2.3598 | 2.252/3.116 | 90.000/90.000 | 31.9 | 0/0/0/0 | rc=0 ended: Terminating on cycle limit wall=106 s CW 0 FATAL 0 |
| BASE2 | 24 | 2.8495 | 2.7739 | 2.667/3.568 | 90.000/90.000 | 31.9 | 0/0/0/0 | rc=0 ended: Terminating on cycle limit wall=116 s CW 0 FATAL 0 |
- base (both arms) mean 2.8599 median 2.7861; new mean 2.4544 median 2.3923; new / base mean 0.8582, median 0.8587; base arms 2.8703 vs 2.8495, new arms 2.4765 vs 2.4322 (mean)
- NEW1 user.hst against BASE1: 26 / 26 rows, printed-bitwise
- NEW1 hydro.hst against BASE1: 26 / 26 rows, printed-bitwise
- NEW2 user.hst against BASE1: 26 / 26 rows, printed-bitwise
- NEW2 hydro.hst against BASE1: 26 / 26 rows, printed-bitwise
- BASE2 user.hst against BASE1: 26 / 26 rows, printed-bitwise
- BASE2 hydro.hst against BASE1: 26 / 26 rows, printed-bitwise

| arm | s/cyc mean / median | peak MiB per GPU (min / max) | user.hst / hydro.hst vs BASE1 | CW/FATAL | end |
|---|---|---|---|---|---|
| ABBA/BASE1 | 2.8703 / 2.7974 | 5 / 76139 (10 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=116 s CW 0  |
| ABBA/NEW1 | 2.4765 / 2.4011 | 53409 / 53839 (4 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=111 s CW 0  |
| ABBA/NEW2 | 2.4322 / 2.3598 | 1465 / 53633 (9 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=106 s CW 0  |
| ABBA/BASE2 | 2.8495 / 2.7739 | 60295 / 76139 (8 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=116 s CW 0  |
| MT_SKIP | 2.3880 / 2.3137 | 14591 / 53549 (9 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=106 s CW 0  |
| ONE4 | 4.0572 / 3.9414 | 104835 / 105343 (4 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=156 s CW 0  |
- MT_SKIP: Waitall cln_hd_copy 234.4 ms per cycle (215.5/249.8 min/max rank, 535 calls), all Waitall 357.4 ms per cycle; the ab9994907 and the base timer: compare with T6r (Waitall cln_hd_copy 584 ms, 644 calls per cycle, 79c1cabef)
- MT_SKIP: 2026-11-25 05:32:41 rc=0 ended: Terminating on cycle limit wall=106 s CW 0 FATAL 0
- ONE4: 2026-11-25 05:35:17 rc=0 ended: Terminating on cycle limit wall=156 s CW 0 FATAL 0

## HG: memory table corrected, and the suite
- The memory minima printed in the HG table above are wrong (two nodes append to one mem.log on NFS, samples get lost or garbled). The peak (max) is right. Parsed again with a strict host / index pattern; minima below are still lower bounds where samples were lost:
| arm | s/cyc mean / median | peak MiB per GPU (min / max) | user.hst / hydro.hst vs BASE1 | CW/FATAL | end |
|---|---|---|---|---|---|
| ABBA/BASE1 | 2.8703 / 2.7974 | 639 / 76139 (8 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=116 s CW 0  |
| ABBA/NEW1 | 2.4765 / 2.4011 | 53409 / 53839 (4 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=111 s CW 0  |
| ABBA/NEW2 | 2.4322 / 2.3598 | 1465 / 53633 (8 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=106 s CW 0  |
| ABBA/BASE2 | 2.8495 / 2.7739 | 60295 / 76139 (8 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=116 s CW 0  |
| MT_SKIP | 2.3880 / 2.3137 | 52633 / 53549 (8 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=106 s CW 0  |
| ONE4 | 4.0572 / 3.9414 | 104835 / 105343 (4 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=156 s CW 0  |
- Read: peak per GPU 79c1cabef 76.1 GB (1000 MiB), e8d3e64d9 and ab9994907 53.8 GB, one node 105.3 GB. Speed (T6, ABBA): e8d3e64d9 / 79c1cabef 0.858 (2.454 vs 2.860 s/cycle mean; base arms 2.870, 2.850, new 2.477, 2.432). The ab9994907 arm with the MPI timer: 2.388 s/cycle, Waitall cln_hd_copy 234 ms per cycle in 535 calls (79c1cabef T6r: 584 ms, 644 calls), so the fill's exchange counts fell by 17 %. ONE4: 4.057 s/cycle (median 3.94), hst bitwise with the 8-GPU base.
- The e8d3e64d9 suite: my first run used the 79c1cabef tests (wrong worktree variable), so it had no star_cs_ring / ray_cut_slab; stopped after 25 of the old tests had passed, restarted 05:36 with the e8d3e64d9 tests. The logs of the wrong run are kept under b1_1124/suiten2_oldtests and runs_suiten2*_oldtests.

## 3. GPU suite, full list in 4 disjoint --only lists, one H200 of the node named in the log each (a 2-rank test shares its one GPU), jobs 1, --star-exe; srun rc of the 4: 0 0 0 0): 1 passed, 0 failed, 987 s;1 passed, 0 failed, 122 s;26 passed, 0 failed, 784 s;26 passed, 0 failed, 770 s;1 passed, 0 failed, 19 s;; test lines 55, SKIP lines 0 (55 tests of this branch)
- CLEAN WARNING: 1 lines in the kept run directories, by test:       1 inf_reject; (inf_reject injects one on purpose).
- star_smoke (dev against the CPU reference, bound 2e-4, nothing changed): `PASS star_smoke            5 s  dev 4.31e-06 warnings 0 worst L_bot`
```
PASS order_space         987 s  err_E [0.0907 0.0208 0.00434] order_E [2.12 2.26] err_v [0.106 0.0379 0.00979] order_v [1.48 1.95]
PASS star_cs_ring        122 s  maxdiff_r3 0 warnings 0 arena 1 fused 1
PASS conservation_lz      36 s  sp_dE 3.33e-16 sp_dM 9.44e-16 sp_dLz 0 sp_dLz_gr 1.47e-16 cs_dE 1.2e-16 cs_dM 2.61e-16 cs_dLz 0 cs_dLz_gr 6.54e-05 restart_maxdiff 0 retry_maxdiff 0 warnings 0
PASS floor_dt             39 s  dE 3.33e-16 dt_gain 0.791 restart_maxdiff 0 retry_maxdiff 0 warnings 0
PASS merge_lat            66 s  sp_dM 4.88e-16 sp_dLz 9.88e-16 sp_chk 1.09e-15 sp_amin 0 sp_dE 2.25e-16 sp_vlat 8.31e-16 sp_rigid_dLz 5.01e-05 cs_dM 9.41e-16 cs_dLz 6.49e-16 cs_chk 1.39e-15 cs_amin 0 
PASS restart_cs           12 s  files 4 maxdiff 0
PASS layout_tile          47 s  maxdiff_1x2 0 maxrel_1x2 0 maxdiff_2x2 9.09e-13 maxrel_2x2 1.31e-10 maxdiff_dfl 0 dE_dfl 2.22e-16 dM_dfl 1.8e-16 dLz_dfl 1.96e-16 maxdiff_fofc 0 dE_fofc 1.11e-16 dM_fof
PASS layout_cs            12 s  warnings 0 maxdiff 0 dE 1.2e-16 dM 1.3e-16 dLz 0 maxdiff_dark 0 dE_dark 1.2e-16 dM_dark 1.3e-16 dLz_dark 0
PASS order_time           21 s  err_E [0.00102 0.000344 0.000113 3.15e-05] order_E [1.57 1.61 1.84] err_F [0.00475 0.00168 0.000588 0.000168] order_F [1.5 1.52 1.81] err_T [0.000167 6.61e-05 2.46e-05 
PASS order_stiff          25 s  err_E [4.31e-06 1.13e-06 2.9e-07] order_E [1.94 1.96] err_v [0.00449 0.00118 0.000305] order_v [1.92 1.96]
PASS nlte_eps1e-4         57 s  L1 0.00133
PASS thermal_eq           21 s  errT_h 0.000552 dEtot_h 1.12e-15 errT_h2 0.000124 dEtot_h2 7e-16 order 2.16
PASS mhd_equil            11 s  dEtot 2.22e-16 dB 0 fatal 3 moved 0.0214 warnings 0
PASS mhd_layout           70 s  cs_r2 1.78e-15 cs_mb22 2e-15 cs_ranks 0 cs_restart 0 cs_dE 2.13e-16 cs_dLz 0 cs_divb 3.13e-15 cs_fofc 88947 sp_r2 2.78e-17 sp_mb22 2.78e-17 sp_ranks 0 sp_restart 0 sp_d
PASS eddington_cs         21 s  dE 0 dM 1.3e-16 dLz 0 dLz_gr 6.48e-05 warnings 0 restart_maxdiff 0 nogd 2.84e-14
PASS ray_clamp             6 s  fatal_fixed 1 fatal_ring 1
PASS ray_cut_slab         86 s  warnings 0 maxdiff_nofuse 0 maxdiff_low 0 maxdiff_retry 0 shared [low retry shared] said_low 1 said_shared 0 fallbacks [0 1] moves 1
PASS radwave             144 s  err_cph [0.00149 0.00147 0.000928 0.000239 0.00226 0.000931] err_damp [0.000865 0.000225 1.63e-05 7.97e-05 0.000655 1.63e-05] damp [0.712 0.164 0.0151 0.000174 0.164 0.
PASS floor_pay            27 s  dE 3.33e-16 dE0 1.67e-16 paid 0.153 unpaid 0 warnings 0 worst_lines 185 worst_bitwise True ekin_floor_dE 3.33e-16 ekin_floor_dM 1.08e-15 ekin_cart_dM 2.22e-16 ekin_dE 3
PASS ceiling_metric       10 s  sp_dvphi 0 sp_dcap 0 sp_vtheta 0 sp_scale_ok 1 sp_capped 1408 cs_dvphi 6.31e-16 cs_dcap 6.31e-16 cs_vtheta 1.5e-16 cs_scale_ok 1 cs_capped 2928
PASS inf_reject            3 s  dE 3.33e-16 nonfinite 0 rejected 1
PASS idort_beams_mb       22 s  maxdiff_tile_x2 0 maxdiff_block_x1 1.13e-10
PASS idort_shadow          3 s  L1 0.0138
PASS idort_shadow_gl       4 s  L1 0.0126
PASS idort_cylinder        9 s  E_in 0.000533 E_out 0.0115 F_out 0.0114 FH_out 0.204
PASS idort_radshock_m5     10 s  L1_rho 0.00171 L1_v 0.000822 L1_Tg 0.00161 L1_Tr 0.000938
PASS idort_thermal_eq2     15 s  errT 0.00196 dEtot 1.88e-16
PASS advect_uniform        5 s  dT 5.38e-10 dE 2.68e-08 warnings 0
PASS conservation_sp      11 s  dE 3.33e-16 dM 8.09e-16 dLz 0 dLz_gr 5.92e-07
PASS conservation_cs       4 s  dE 1.2e-16 dM 1.3e-16 dLz 0 dLz_gr 6.48e-05
PASS floorcap_spin        48 s  sp_dE 3.33e-16 sp_dM 8.09e-16 sp_dLz 1.47e-16 sp_dLz_gr 1.47e-16 n_v 5.85e+04 cs_dE 1.2e-16 cs_dM 1.3e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.3e-05 nospin_sp_dLz_gr 0.366 resta
PASS restart_fixes        86 s  warnings 0 notes 6 stored 6 maxdiff_lat_wb 0 moved_lat_wb 0.00177 maxdiff_fofc 0 moved_fofc 0.0296 maxdiff_merge 0 moved_merge 0.00787 maxdiff_blockend 0 moved_blockend
PASS restart_sp           11 s  files 4 maxdiff 0
PASS dt_retry             11 s  files 6 maxdiff 0 retries_fac1 1 retries_fac05 1
PASS fofc_layout_cs       29 s  maxdiff_mb22 0 maxdiff_r2 0 maxdiff_r2mb22 0 fo_fofc 9.44e+04 maxdiff_restart 0 maxdiff_retry 0 dE 2.16e-16 dLz 0 dM 4.31e-16
PASS nlte_eps1e-1         58 s  L1 9.58e-05
PASS mhd_vet              99 s  warnings 0 cs_layout 4.52e-15 cs_ranks 0 cs_retry 0 cs_divb 1.02e-14 cs_restart 0 cs_dM 2.61e-16 cs_dLz 0 cs_dE 1.11e-16 cs_moved 0.443 sp_layout 0 sp_ranks 0 sp_retry 
PASS wb_sp                31 s  vmax_over_cs 6.08e-09 vmax_over_cs_rot 6.1e-09
PASS wb_cs                16 s  vmax_over_cs 5.11e-10 vmax_over_cs_rot 7.4e-09
PASS shell_force_cs        4 s  lat 5.2e-05 rad 5.56e-05
PASS momentum_box          5 s  dP_max 8.1e-17 dP_cyc [7.3e-17 8.04e-17 7.79e-17 8.02e-17 7.77e-17 7.92e-17 8.1e-17 7.96e-17] warnings 0
PASS wedge_sym            34 s  fphi [1.23e-09 2.48e-10 1.31e-09] warnings 0 wedge_fatal 1 fth 1.24e-16 theta_fatal 1 fth_rot 1.56e-09 rot_fatal 1
PASS floor_momentum       58 s  box_ekin [3.27e-16 2.33e-16 8.95e-17 2.06e-16 3.11e-17 9.74e-30 1.87e-16 1.9e-30] box_both [2.83e-16 2.73e-16 3.98e-16 2.94e-16 2.57e-16 2.15e-16 4.28e-17 2.14e-17] box
PASS star_smoke            5 s  dev 4.31e-06 warnings 0 worst L_bot
PASS wb_star              11 s  warnings 0 prof_in_vr 0.00226 prof_in_lat 4.57e-09 rigid_in_vr 0.00225 rigid_in_lat 0.000213 eq_in_vr 9.41e-05 eq_in_lat 2.13e-09
PASS mlt_time             20 s  warnings 0 err [0.008 0.00109 0.000155] order 2.81 err_start [0.0714 0.0369 0.0188] order_start 0.976
PASS idort_beams          11 s  L1 0.177 lum_top 1.19 flux_err 5.97e-08 F_nocone 0.683
PASS idort_radshock_m2     11 s  L1_rho 0.000356 L1_v 0.000213 L1_Tg 0.000385 L1_Tr 0.00026
PASS diffusion            31 s  errD [-0.000815 -7.79e-08 1.42e-10] L1 [0.00137 0.00117 0.00117] warnings 0 errD_tau1 -0.0817 L1_tau1 0.0118 errD_2d -9.29e-07 L1_2d 0.00659 th_min 1.15e-17 th_max -0.0
PASS radshock             25 s  L1_rho 0.00438 L1_v 0.00358 L1_Tg 0.00659 L1_Tr 0.0107 Emin 1.72e+11 warnings 0
PASS radshock_m5          25 s  L1_rho 0.00386 L1_v 0.00212 L1_Tg 0.0062 L1_Tr 0.00744 Emin 1.79e+11 warnings 0
PASS marshak              12 s  Emin 8.92e-19 Tmin 9.96e-13 Emax 0.995 Tmax 0.999 Erise 3.07e-18 x_front 0.535 warnings 0
PASS advect_pulse         86 s  L1 [0.000872 0.000574 0.000821 0.000776] warnings 0
PASS shadow               29 s  umbra 0.00609 lit_err 0.00825 n_umbra 774 n_lit 3792 warnings 0
PASS mpi_ranks            19 s  warnings 0 fcone_lz 8 files 15 maxdiff_sp 0 maxdiff_cs 0 maxdiff_pay 0
```
### Printed values against next 79c1cabef (25c suite, common tests)
- next2: 55 lines, PASS 55, FAIL 0: []

- against next 79c1cabef: 53 common tests, printed values equal in 53, different in 0, only in rel-e ['star_cs_ring', 'ray_cut_slab']

## m35zams (11-25 05:57): rst m35zams_m35zams.00007.rst (hard link, c148944), next2 f3c50eae1b2ebb30a02dd3f498722bf6 against athena_next_79c1cabef; own keys: rad_vet/time2_gtime_slope=true rad_vet/clean_stag_thinf=true rad_vet/clean_stag_lat_sums=true rad_vet/clean_stag_tie_gid=true rad_vet/clean_stag_opac_drag=true rad_vet/vet_gd_allow_clamp=false hydro/vceil=3.0e8 rad_vet/clean_pc_sweep=rb rad_vet/clean_stag_mg_levels=2 time/cfl_number=0.3 rad_vet/clean_conv_max=1e-5
- R1: 2026-11-25 05:54:41 rc=0 (srun) ended: Terminating on cycle limit wall=24 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- R50: 2026-11-25 05:55:13 rc=0 (srun) ended: Terminating on cycle limit wall=32 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- R25a: 2026-11-25 05:55:51 rc=0 (srun) ended: Terminating on cycle limit wall=38 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- R25b: 2026-11-25 05:56:21 rc=0 (srun) ended: Terminating on cycle limit wall=30 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- C50: 2026-11-25 05:56:56 rc=0 (srun) ended: Terminating on cycle limit wall=34 s CW 0 FATAL 0 inf 0 DT_RETRY 0
#### NOTE lines of the restart (next line, 1 step)
- ### NOTE: <problem>/he_* (hs_*) keys are now named star_*; read so: he_bc_inner he_bc_linhao he_bc_outer he_bc_outer_r2 he_esrc_const he_gm he_gravity he_grid_dump he_ic_balance he_ic_cols he_ic_eint_from_t he_ic_file he_ic_profile he_ic_rad he_ic_rad_atm_post
- ### NOTE: <hydro>/fofc_flag_exchange = true (restart file without it: the fix applies; name fofc_flag_exchange=false to keep the old behaviour)
- ### NOTE: <rad_vet>/time2_gtime_pos = true (restart file without it: the fix applies; name time2_gtime_pos=false to keep the old behaviour)
- NOTE lines: 3; "the fix applies": 2; "added from the command line": 0; "old behaviour": 2; FATAL 0, CLEAN WARNING 0; clean_init lines: []
#### 50 straight against 25 + rst + 25
- user.hst: R50 50 rows, R25a 25, R25b 25; rows not printed-identical to R50 at the same t: 0 (first t -) -> BITWISE (printed)
- hydro.hst: R50 50 rows, R25a 25, R25b 25; rows not printed-identical to R50 at the same t: 0 (first t -) -> BITWISE (printed)
- final hydro_w dump, payload md5: R50 9dcb7af7497b, R25b 9dcb7af7497b -> IDENTICAL
- final vet dump, payload md5: R50 1a2812d0dcdd, R25b 1a2812d0dcdd -> IDENTICAL
#### against next 79c1cabef (C50, same rst, same keys, 50 cycles)
- final hydro_w dump, payload md5: next2 R50 9dcb7af7497b, 79c1cabef C50 9dcb7af7497b -> IDENTICAL
- final vet dump, payload md5: next2 R50 1a2812d0dcdd, 79c1cabef C50 1a2812d0dcdd -> IDENTICAL
- user.hst next2 R50 against 79c1cabef C50: 50 rows, printed-bitwise
- hydro.hst next2 R50 against 79c1cabef C50: 50 rows, printed-bitwise
- row 0: t 5.18405012e+05 / 5.18405012e+05, dt 3.480520 / 3.480520, L_top rel 0.00e+00, KE_int rel 0.00e+00, M_tot rel 0.00e+00, Etot rel 0.00e+00
- row 1: t 5.18408492e+05 / 5.18408492e+05, dt 3.480520 / 3.480520, L_top rel 0.00e+00, KE_int rel 0.00e+00, M_tot rel 0.00e+00, Etot rel 0.00e+00
- row 2: t 5.18411973e+05 / 5.18411973e+05, dt 3.480520 / 3.480520, L_top rel 0.00e+00, KE_int rel 0.00e+00, M_tot rel 0.00e+00, Etot rel 0.00e+00
- row 3: t 5.18415453e+05 / 5.18415453e+05, dt 3.480520 / 3.480520, L_top rel 0.00e+00, KE_int rel 0.00e+00, M_tot rel 0.00e+00, Etot rel 0.00e+00
- row 4: t 5.18418934e+05 / 5.18418934e+05, dt 3.480520 / 3.480520, L_top rel 0.00e+00, KE_int rel 0.00e+00, M_tot rel 0.00e+00, Etot rel 0.00e+00
- row 9: t 5.18436336e+05 / 5.18436336e+05, dt 3.480520 / 3.480520, L_top rel 0.00e+00, KE_int rel 0.00e+00, M_tot rel 0.00e+00, Etot rel 0.00e+00
- row 24: t 5.18488544e+05 / 5.18488544e+05, dt 3.480520 / 3.480520, L_top rel 0.00e+00, KE_int rel 0.00e+00, M_tot rel 0.00e+00, Etot rel 0.00e+00
- row 49: t 5.18575557e+05 / 5.18575557e+05, dt 3.480520 / 3.480520, L_top rel 0.00e+00, KE_int rel 0.00e+00, M_tot rel 0.00e+00, Etot rel 0.00e+00
- R50: s/cyc (elapsed prints) 0.2716 (c148950-148990); its sum per cycle 28.6, solves 4.30; CW 0 FATAL 0 inf 0 DT_RETRY 0
- C50: s/cyc (elapsed prints) 0.2827 (c148950-148990); its sum per cycle 28.6, solves 4.30; CW 0 FATAL 0 inf 0 DT_RETRY 0
#### END-row residuals over the 50 cycles (R_E = dEtot - dt (L_in - L_top + P_esrc - P_sponge - E_hflx); R_M = dM_tot + dt M_hflx; bsg: R_L = dLz_tot + dt Lz_hflx)
- R50: R_E / Etot RMS 1.28e-15 max 3.13e-15 sum 1.29e-15; R_M / M_tot RMS 1.14e-15 max 2.48e-15 sum 2.48e-16
- C50: R_E / Etot RMS 1.28e-15 max 3.13e-15 sum 1.29e-15; R_M / M_tot RMS 1.14e-15 max 2.48e-15 sum 2.48e-16
