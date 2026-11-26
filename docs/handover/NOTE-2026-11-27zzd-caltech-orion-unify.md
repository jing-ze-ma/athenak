# NOTE 2026-11-27zzd caltech -> orion: sdirk2 with <mhd> on the GPU (your 27zzd)

## Builds (worktree wt_orion1121/un at 0b480264c; CUDA 12.9 + g++13, HOPPER90, -ffp-contract=off; 02-03)
- BUILD un rc=0 wall=614 s TUs=280
- BUILD_OK un md5=f39250ae3f7e2f956bf05b4b39295443 prov=0b480264c45736b23f1b3617cf27c9e652065ae2
- GATE PASS un: athena_un_0b480264c md5 f39250ae node hpc-sm-02-03 pass1=3.720e-11 log date 11-26 01:44:31 end: Terminating on cycle limit
- BUILD suite rc=0 wall=608 s
- BUILD_OK suite md5=62764ea88484162f37fa92b8826a0fd9 prov=0b480264c45736b23f1b3617cf27c9e652065ae2

## The six MHD tests on 1 H200 (02-17, packed next to non-timing arms; run_tests.py --only mhd_equil,mhd_vet,mhd_geneos,mhd_radwave,mhd_wb,order_mhd --jobs 4 --threads 1, mpirun --oversubscribe)
### sdirk2 (default): 6 passed, 0 failed, 763 s
```
PASS mhd_equil            17 s  dEtot 2.22e-16 dB 0 fatal 3 moved 0.0215 warnings 0
PASS mhd_geneos          328 s  warnings 0 cs_layout 6.78e-14 cs_ranks 0 cs_retry 0 cs_divb 9.56e-15 cs_restart 0 cs_dM 2.61e-16 cs_dLz 0 cs_dE 1.11e-16 cs_moved 0.00406 sp_layout 0 sp_ranks 0 sp_retry 0 sp_divb 1.14e-15 sp_restart 0 sp_dM 2.7e-16 sp_dLz 0 sp_dE 1.36e-16 sp_moved 0.107 ct_layout 0 ct_ranks 0 ct_retry 0 ct_divb 2.04e-14 ct_restart 0 ct_dM 0 ct_dE 0 ct_moved 0.0361
PASS mhd_vet             209 s  warnings 0 cs_layout 3.8e-15 cs_ranks 0 cs_retry 0 cs_divb 1.36e-14 cs_restart 0 cs_dM 2.61e-16 cs_dLz 0 cs_dE 2.23e-16 cs_moved 0.443 sp_layout 0 sp_ranks 0 sp_retry 0 sp_divb 1.29e-15 sp_restart 0 sp_dM 2.7e-16 sp_dLz 0 sp_dE 1.36e-16 sp_moved 3.5e-05 ct_layout 0 ct_ranks 0 ct_retry 0 ct_divb 1.73e-14 ct_restart 0 ct_dM 2.22e-16 ct_dE 0 ct_moved 5.41e-05
PASS order_mhd           643 s  err_E [4.23e-09 9.6e-10 3.37e-10] order_E [2.14 1.51] err_v [0.00823 0.00207 0.000509] order_v [1.99 2.03] err_B [8.79e-06 2.21e-06 5.43e-07] order_B [1.99 2.02] err_E_stiff [6.01e-06 1.53e-06 3.73e-07] order_E_stiff [1.97 2.03] err_v_stiff [0.00675 0.00172 0.00042] order_v_stiff [1.97 2.04] err_B_stiff [6.46e-06 1.63e-06 3.95e-07] order_B_stiff [1.99 2.05] err_E_ni [4.08e-06 1.03e-06 2.59e-07] order_E_ni [1.99 1.99] err_v_ni [0.00509 0.00129 0.000328] order_v_ni [1.98 1.98] err_B_ni [3.72e-06 9.41e-07 2.38e-07] order_B_ni [1.98 1.98]
PASS mhd_wb              642 s  static [6.66e-09 5.55e-10 1.68e-08] field [6.75e-09 5.87e-10 2.6e-08] geneos [6.74e-09 1.68e-08] dyn [2.59e-08 2.96e-08 2.59e-08] warnings 0 monopole 4.37e-06 restart_maxdiff 0
PASS mhd_radwave         763 s  err_cph [0.000755 0.000658 0.000215 0.00014 0.000292 0.000474 0.000905 0.000837 0.000684 0.000496 0.000995 0.00145 0.00118 0.000325] err_damp [0.000289 0.000294 0.000223 0.000263 0.000118 0.000642 5.62e-05 0.000505 0.000289 0.00175 0.00118 0.00108 0.00237 0.00325] ord_damp [3.06 3.06 3.03 3.04 0.717 2.42 4.61 0.483 3.06 3.02 2.95 3 3.02 2.95] warnings 0 nonideal_cph [0.000172 0.000794] nonideal_damp [0.000197 0.000179] geneos [0.000212 0.000223]
```
### hesdirk2 (--keys rad_vet/time2_tableau=hesdirk2): 5 passed, 1 failed, 735 s
```
PASS mhd_equil             9 s  dEtot 3.33e-16 dB 0 fatal 3 moved 0.0214 warnings 0
PASS mhd_geneos          299 s  warnings 0 cs_layout 7.51e-14 cs_ranks 0 cs_retry 0 cs_divb 1.22e-14 cs_restart 0 cs_dM 2.61e-16 cs_dLz 0 cs_dE 1.11e-16 cs_moved 0.00409 sp_layout 0 sp_ranks 0 sp_retry 0 sp_divb 1.35e-15 sp_restart 0 sp_dM 2.7e-16 sp_dLz 0 sp_dE 1.36e-16 sp_moved 0.107 ct_layout 0 ct_ranks 0 ct_retry 0 ct_divb 1.69e-14 ct_restart 0 ct_dM 2.22e-16 ct_dE 0 ct_moved 0.0361
PASS mhd_vet             211 s  warnings 0 cs_layout 4.52e-15 cs_ranks 0 cs_retry 0 cs_divb 1.02e-14 cs_restart 0 cs_dM 2.61e-16 cs_dLz 0 cs_dE 1.11e-16 cs_moved 0.443 sp_layout 0 sp_ranks 0 sp_retry 0 sp_divb 1.31e-15 sp_restart 0 sp_dM 4.05e-16 sp_dLz 0 sp_dE 1.36e-16 sp_moved 3.41e-05 ct_layout 0 ct_ranks 0 ct_retry 0 ct_divb 1.42e-14 ct_restart 0 ct_dM 1.11e-16 ct_dE 0 ct_moved 5.68e-05
FAIL order_mhd           604 s  err_E [1.31e-09 1.82e-09 3.45e-10] order_E [-0.474 2.4] err_v [0.00826 0.00209 0.000515] order_v [1.98 2.02] err_B [8.82e-06 2.23e-06 5.53e-07] order_B [1.99 2.01] err_E_stiff [5.5e-06 1.39e-06 3.43e-07] order_E_stiff [1.98 2.02] err_v_stiff [0.00617 0.00156 0.000382] order_v_stiff [1.98 2.03] err_B_stiff [5.85e-06 1.46e-06 3.55e-07] order_B_stiff [2 2.04] err_E_ni [3.71e-06 9.23e-07 2.33e-07] order_E_ni [2.01 1.99] err_v_ni [0.00462 0.00116 0.000293] order_v_ni [2 1.98] err_B_ni [3.36e-06 8.33e-07 2.1e-07] order_B_ni [2.01 1.99]  <- err_E_stiff [5.5e-06 1.39e-06 3.43e-07] not within 0.05 of 3.73e-07; err_v_stiff [0.00617 0.00156 0.000382] not within 0.05 of 0.00042; err_B_stiff [5.85e-06 1.46e-06 3.55e-07] not within 0.05 of 3.95e-07; err_E_ni [3.71e-06 9.23e-07 2.33e-07] not within 0.05 of 2.59e-07; err_v_ni [0.00462 0.00116 0.000293] not within 0.05 of 
PASS mhd_wb              615 s  static [6.68e-09 5.52e-10 1.68e-08] field [6.77e-09 5.87e-10 2.6e-08] geneos [6.75e-09 1.68e-08] dyn [2.59e-08 2.96e-08 2.59e-08] warnings 0 monopole 4.37e-06 restart_maxdiff 0
PASS mhd_radwave         735 s  err_cph [0.000754 0.000658 0.000224 0.000141 0.000293 0.000474 0.000905 0.000836 0.000684 0.000536 0.001 0.00145 0.00126 0.000342] err_damp [0.000288 0.000294 0.000219 0.000262 0.000118 0.000642 5.52e-05 0.000506 0.000289 0.00173 0.00118 0.00108 0.00232 0.00324] ord_damp [3.06 3.06 3.04 3.04 0.713 2.42 4.63 0.468 3.06 3.02 2.96 3 3.03 2.96] warnings 0 nonideal_cph [0.000173 0.000794] nonideal_damp [0.000197 0.000179] geneos [0.000221 0.000219]
```
- CLEAN WARNING: 0 in all run directories (both tableaus).
- **order_mhd, sdirk2, GPU against your CPU references:**
  - finest err_v 5.09e-4, err_E_stiff 3.73e-7, err_E_ni 2.59e-7: GPU minus CPU = 0 at the printed 3 digits;
  - finest orders v 2.03 / B 2.02 (thin); E 2.03 / v 2.04 / B 2.05 (stiff); 1.99 / 1.98 / 1.98 (Ohm + ambipolar): all equal to yours.
  - mhd_vet / mhd_geneos divb <= 2.0e-14, dE <= 2.2e-16, restart 0.
  - mhd_radwave err_cph max 1.45e-3, err_damp max 3.25e-3: equal to yours.
- **hesdirk2:** the only FAIL is order_mhd.
  - The stiff errors differ from the sdirk2 reference values the test compares against (err_E_stiff 3.43e-7 vs 3.73e-7, err_v_stiff 3.82e-4 vs 4.20e-4).
  - The thin E error is non-monotone at the coarse level (err_E 1.31e-9, 1.82e-9, 3.45e-10; order -0.47 / 2.40).
  - All other orders are 1.98-2.04.
  - That reads as order_mhd's references being sdirk2 values, not a hesdirk2 regression. Your call.
- sd.json and he.json: docs/handover/caltech_results/UNIFY/.
