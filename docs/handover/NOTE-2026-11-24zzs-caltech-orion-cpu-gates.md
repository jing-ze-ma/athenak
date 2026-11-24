# NOTE 2026-11-24zzs caltech -> orion: CPU gates A (idort, 29046c3f2) and B (compton f8e55fee7) on 28 idle cores of hpc-sm-01-15 (your 24zzq)

- CPU step: srun --overlap --gres=none -c 28 on 7333036 / hpc-sm-01-15 beside the M35ZAMS production; g++ 13.2 + hpcx mpicxx, OpenMP + MPI, -ffp-contract=off, no PROBLEM (built-in pgens). Build dirs KEPT: wt_orion1121/wt_cpu/build_cpu (29046c3f2) and wt_orion1121/wt_cpu_c/build_cpu (f8e55fee7).
- M35ZAMS s/cyc guard (stop if > 2 % slower): baseline 0.29827; monitor: 
```
/resnick/software9/external/hpcx/2.17.1/ompi/bin/mpicxx
g++ (GCC) 11.5.0 20240719 (Red Hat 11.5.0-2)
BUILD wt_cpu rc=0 65 s 29046c3f2 md5 961ed1195355cc9dadcf19e762a9c943
BUILD wt_cpu_c rc=0 65 s f8e55fee7 md5 db6ffddb793a25a76b361b836eb4f9ee
== Gate A: idort (29046c3f2)
clean VET test suite: 0 tests, full, 6 at a time x 4 threads, binary /resnick/groups/carnegie_poc/jingze/wt_orion1121/wt_cpu/build_cpu/src/athena
0 passed, 0 failed, 0 s
== Gate A: --quick
clean VET test suite: 16 tests, quick, 6 at a time x 4 threads, binary /resnick/groups/carnegie_poc/jingze/wt_orion1121/wt_cpu/build_cpu/src/athena
PASS conservation_cs       5 s  dE 8.4e-16 dM 3.91e-16 dLz 1.17e-16 dLz_gr 6.44e-05
PASS conservation_sp       6 s  dE 5.56e-16 dM 5.4e-16 dLz 0 dLz_gr 6.94e-07
PASS restart_sp           13 s  files 4 maxdiff 0
PASS dt_retry             17 s  files 6 maxdiff 0 retries_fac1 1 retries_fac05 1
PASS momentum_box          6 s  dP_max 5.87e-17 dP_cyc [5.73e-17 5.87e-17 5.75e-17 5.87e-17 5.43e-17 5.47e-17 5.28e-17 5.49e-17] warnings 0
PASS wb_sp                20 s  vmax_over_cs 4.57e-09 vmax_over_cs_rot 4.67e-09
PASS order_time           27 s  err_E [0.000113 3.15e-05] order_E [1.84] err_F [0.000589 0.000168] order_F [1.81] err_T [2.45e-05 6.9e-06] order_T [1.83]
PASS nlte_eps1e-1         24 s  L1 9.58e-05
PASS idort_shadow          2 s  L1 0.0138
PASS ceiling_metric        7 s  sp_dvphi 0 sp_dcap 0 sp_vtheta 0 sp_scale_ok 1 sp_capped 1408 cs_dvphi 6.31e-16 cs_dcap 6.31e-16 cs_vtheta 1.5e-16 cs_scale_ok 1 cs_capped 2928
PASS idort_shadow_gl       3 s  L1 0.0126
PASS advect_uniform        3 s  dT 1.45e-06 dE 1.26e-05 warnings 0
PASS floor_pay            25 s  dE 7.78e-16 dE0 4.45e-16 paid 0.286 unpaid 0 warnings 0 worst_lines 197 worst_bitwise True ekin_dE 6.67e-16 ekin_dM 7.19e-16 ekin_dLz 1.96e-16
PASS wedge_sym            28 s  fphi [1.58e-09 1.95e-11 1.87e-12] warnings 0 wedge_fatal 1 fth 1.6e-14 theta_fatal 1 fth_rot 1.65e-09 rot_fatal 1
PASS diffusion            22 s  errD [-0.000815 -7.79e-08 1.42e-10] L1 [0.00137 0.00117 0.00117] warnings 0 errD_tau1 -0.0817 L1_tau1 0.0118 errD_2d -9.3e-07 L1_2d 0.00659 th_min 1.15e-17 th_max -0.0392
PASS layout_tile          77 s  maxdiff_1x2 0 maxrel_1x2 0 maxdiff_2x2 2.27e-13 maxrel_2x2 7.81e-14 maxdiff_dfl 0 dE_dfl 1e-15 dM_dfl 3.6e-16 dLz_dfl 3.92e-16 maxdiff_fofc 0 dE_fofc 1e-15 dM_fofc 0.000636 dLz_fofc 0 fo_f
16 passed, 0 failed, 77 s
== Gate B: compton branch f8e55fee7, key off
clean VET test suite: 4 tests, full, 6 at a time x 4 threads, binary /resnick/groups/carnegie_poc/jingze/wt_orion1121/wt_cpu_c/build_cpu/src/athena
PASS compton              10 s  err [1.67e-15 4.44e-16] dEtot [2.8e-16 4.55e-16] overshoot [0 0]
PASS wb_sp                23 s  vmax_over_cs 4.57e-09 vmax_over_cs_rot 4.67e-09
PASS thermal_eq           26 s  errT_h 0.000552 dEtot_h 1.26e-15 errT_h2 0.000124 dEtot_h2 1.26e-15 order 2.16
PASS wb_cs               240 s  vmax_over_cs 1.18e-09 vmax_over_cs_rot 6.94e-09
4 passed, 0 failed, 240 s
== Gate B: key on (--keys rad_vet/compton=true)
clean VET test suite: 4 tests, full, 6 at a time x 4 threads, binary /resnick/groups/carnegie_poc/jingze/wt_orion1121/wt_cpu_c/build_cpu/src/athena
FAIL wb_sp                 1 s    <- ERROR run: return code 1 (see /resnick/groups/carnegie_poc/jingze/cpugates_1124/runs_c_on/wb_sp/run/run.log)
FAIL wb_cs                 1 s    <- ERROR run: return code 1 (see /resnick/groups/carnegie_poc/jingze/cpugates_1124/runs_c_on/wb_cs/run/run.log)
PASS compton              10 s  err [1.67e-15 4.44e-16] dEtot [2.8e-16 4.55e-16] overshoot [0 0]
PASS thermal_eq           26 s  errT_h 0.000552 dEtot_h 1.26e-15 errT_h2 0.000124 dEtot_h2 1.26e-15 order 2.16
2 passed, 2 failed, 26 s
INNER DONE
```

## Gate A idort rerun with test names (your 24zzw; same build_cpu 29046c3f2, no rebuild, 28 cores of 01-15, 11-24 12:42)
- M35ZAMS s/cyc guard: baseline 0.32305; after 0.32708; monitor: 
```
29046c3f2
clean VET test suite: 8 tests, full, 6 at a time x 4 threads, binary /resnick/groups/carnegie_poc/jingze/wt_orion1121/wt_cpu/build_cpu/src/athena
PASS idort_shadow          3 s  L1 0.0138
PASS idort_shadow_gl       4 s  L1 0.0126
PASS idort_beams          13 s  L1 0.777 lum_top 1.79 flux_err 2.66e-07
PASS idort_radshock_m2     15 s  L1_rho 0.000383 L1_v 0.000245 L1_Tg 0.000396 L1_Tr 0.000271
PASS idort_radshock_m5     14 s  L1_rho 0.00216 L1_v 0.000919 L1_Tg 0.00181 L1_Tr 0.00103
PASS idort_cylinder       23 s  E_in 0.000531 E_out 0.0111 F_out 0.0116 FH_out 0.203
PASS idort_thermal_eq2     23 s  errT 0.00196 dEtot 1.88e-16
PASS idort_beams_mb       40 s  maxdiff_tile_x2 0 maxdiff_block_x1 2.05e-11
8 passed, 0 failed, 40 s
INNER DONE
```
