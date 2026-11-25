# NOTE 2026-11-24zzzze1 caltech -> orion: GPU gate of rel-b c8a0fc099 (your 24zzzze1)

## 1. Builds (own worktree wt_orion1121/relb at c8a0fc099, CUDA 12.9 + g++ 13, HOPPER90, CPU steps on 02-09): star_vet (PROBLEM = star_vet) and built_in
- BUILD relb rc=0 wall=696 s TUs=250
- BUILD_OK relb md5=a52cb2f731ec218b4819b4883c8e28b8 prov=c8a0fc09973ccb2127e715736ed678d6873d21b0
- GATE PASS relb: athena_relb_c8a0fc099 md5 a52cb2f7 node hpc-sm-02-09 pass1=3.686e-11 log date 11-24 19:08:28 end: Terminating on cycle limit
- BUILD suite rc=0 wall=689 s
- BUILD_OK suite md5=ed110fa0d014eaa50052bd3a0cc24251 prov=c8a0fc09973ccb2127e715736ed678d6873d21b0
- **The tiny gate deck is a copy with one line added, clean_face_hctau = false in its <rad_m1> block** (gate_nols_relb.athinput; pass 1 = 3.686e-11 wanted).

## 1. Builds (own worktree wt_orion1121/relb at c8a0fc099, CUDA 12.9 + g++ 13, HOPPER90, CPU steps on 02-09): star_vet (PROBLEM = star_vet) and built_in
- BUILD relb rc=0 wall=696 s TUs=250
- BUILD_OK relb md5=a52cb2f731ec218b4819b4883c8e28b8 prov=c8a0fc09973ccb2127e715736ed678d6873d21b0
- GATE PASS relb: athena_relb_c8a0fc099 md5 a52cb2f7 node hpc-sm-02-09 pass1=3.686e-11 log date 11-24 19:08:28 end: Terminating on cycle limit
- BUILD suite rc=0 wall=689 s
- BUILD_OK suite md5=ed110fa0d014eaa50052bd3a0cc24251 prov=c8a0fc09973ccb2127e715736ed678d6873d21b0
- **The tiny gate deck is a copy with one line added, clean_face_hctau = false in its <rad_m1> block** (gate_nols_relb.athinput; pass 1 = 3.686e-11 wanted).

## 4. Old-rst check, M35ZAMS (4 ranks, 02-09): rst written by 95813cf85 (my zzzzd6 wr8, c22163, slope and thinf false stored), restarted with NO physics key named (keys absent; outputs off and nlim only). base16: 95813cf85, 16 cycles; rel8 + rel8b: rel-b star_vet, 8 + 8
- base16: 2026-11-24 19:11:33 rc=0 (srun) ended: Terminating on cycle limit wall=29 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- rel8: 2026-11-24 19:12:49 rc=0 (srun) ended: Terminating on cycle limit wall=43 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- rel8b: 2026-11-24 19:15:04 rc=0 (srun) ended: Terminating on cycle limit wall=135 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- user.hst rel-b first 8: 8 rows, printed-identical to the base 16-cycle run at the same t: 8; first different or missing t: none
- user.hst rel-b second 8: 8 rows, printed-identical to the base 16-cycle run at the same t: 8; first different or missing t: none
- user.hst: base 16 rows; BITWISE (printed)
- hydro.hst rel-b first 8: 8 rows, printed-identical to the base 16-cycle run at the same t: 8; first different or missing t: none
- hydro.hst rel-b second 8: 8 rows, printed-identical to the base 16-cycle run at the same t: 8; first different or missing t: none
- hydro.hst: base 16 rows; BITWISE (printed)

## 3. GPU suite, full list in 4 disjoint --only lists, one H200 of 02-03 each (a 2-rank test shares its one GPU), jobs 1, --star-exe; srun rc of the 4: 0 0 0 0): 1 passed, 0 failed, 979 s;13 passed, 0 failed, 352 s;13 passed, 0 failed, 353 s;13 passed, 0 failed, 313 s;; test lines 40, SKIP lines 0 (40 expected)
- CLEAN WARNING: 1 lines in the kept run directories, by test:       1 inf_reject; (inf_reject injects one on purpose).
- star_smoke (dev against the CPU reference, bound 2e-4, nothing changed): `PASS star_smoke            4 s  dev 4.1e-06 warnings 0 worst L_bot`
```
PASS order_space         979 s  err_E [0.0907 0.0208 0.00434] order_E [2.12 2.26] err_v [0.106 0.0379 0.00979] order_v [1.48 1.95]
PASS restart_sp           18 s  files 4 maxdiff 0
PASS restart_cs           16 s  files 4 maxdiff 0
PASS order_stiff          30 s  err_E [4.31e-06 1.13e-06 2.9e-07] order_E [1.94 1.96] err_v [0.00449 0.00118 0.000305] order_v [1.92 1.96]
PASS thermal_eq           21 s  errT_h 0.000552 dEtot_h 1.12e-15 errT_h2 0.000124 dEtot_h2 7e-16 order 2.16
PASS radwave             157 s  err_cph [0.00149 0.00147 0.000928 0.000239 0.00226 0.000931] err_damp [0.000865 0.000225 1.63e-05 7.97e-05 0.000655 1.63e-05] damp [0.712 0.164 0.0151 0.000174 0.164 0.
PASS wedge_sym            36 s  fphi [1.23e-09 2.48e-10 1.31e-09] warnings 0 wedge_fatal 1 fth 1.24e-16 theta_fatal 1 fth_rot 1.56e-09 rot_fatal 1
PASS floor_pay            14 s  dE 3.33e-16 dE0 1.67e-16 paid 0.153 unpaid 0 warnings 0 worst_lines 185 worst_bitwise True ekin_dE 3.33e-16 ekin_dM 5.4e-16 ekin_dLz 1.47e-16
PASS star_smoke            4 s  dev 4.1e-06 warnings 0 worst L_bot
PASS idort_beams           6 s  L1 0.177 lum_top 1.19 flux_err 5.97e-08
PASS idort_shadow_gl       5 s  L1 0.0126
PASS idort_radshock_m5     10 s  L1_rho 0.00171 L1_v 0.000822 L1_Tg 0.00161 L1_Tr 0.000938
PASS advect_uniform        7 s  dT 5.38e-10 dE 2.68e-08 warnings 0
PASS shadow               28 s  umbra 0.00609 lit_err 0.00825 n_umbra 774 n_lit 3792 warnings 0
PASS conservation_sp       9 s  dE 3.33e-16 dM 8.09e-16 dLz 1.47e-16 dLz_gr 5.92e-07
PASS conservation_cs       5 s  dE 1.2e-16 dM 2.61e-16 dLz 1.17e-16 dLz_gr 6.48e-05
PASS conservation_lz      36 s  sp_dE 3.33e-16 sp_dM 9.44e-16 sp_dLz 2.94e-16 sp_dLz_gr 1.47e-16 cs_dE 2.4e-16 cs_dM 2.61e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.54e-05 restart_maxdiff 0 retry_maxdiff 0 warni
PASS floorcap_spin        42 s  sp_dE 3.33e-16 sp_dM 8.09e-16 sp_dLz 1.47e-16 sp_dLz_gr 1.47e-16 n_v 5.85e+04 cs_dE 1.2e-16 cs_dM 1.3e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.3e-05 nospin_sp_dLz_gr 0.366 resta
PASS dt_retry             10 s  files 6 maxdiff 0 retries_fac1 1 retries_fac05 1
PASS layout_tile          41 s  maxdiff_1x2 0 maxrel_1x2 0 maxdiff_2x2 9.09e-13 maxrel_2x2 1.31e-10 maxdiff_dfl 0 dE_dfl 2.22e-16 dM_dfl 1.8e-16 dLz_dfl 3.92e-16 maxdiff_fofc 0 dE_fofc 1.11e-16 dM_fof
PASS wb_sp                29 s  vmax_over_cs 6.08e-09 vmax_over_cs_rot 6.1e-09
PASS wb_cs                18 s  vmax_over_cs 1.36e-09 vmax_over_cs_rot 6.93e-09
PASS momentum_box          6 s  dP_max 8.1e-17 dP_cyc [7.3e-17 8.04e-17 7.79e-17 8.02e-17 7.77e-17 7.92e-17 8.1e-17 7.96e-17] warnings 0
PASS ceiling_metric       14 s  sp_dvphi 0 sp_dcap 0 sp_vtheta 0 sp_scale_ok 1 sp_capped 1408 cs_dvphi 6.31e-16 cs_dcap 6.31e-16 cs_vtheta 1.5e-16 cs_scale_ok 1 cs_capped 2928
PASS idort_beams_mb       27 s  maxdiff_tile_x2 0 maxdiff_block_x1 7.61e-12
PASS radshock_m5          27 s  L1_rho 0.00386 L1_v 0.00212 L1_Tg 0.0062 L1_Tr 0.00744 Emin 1.79e+11 warnings 0
PASS advect_pulse         89 s  L1 [0.000872 0.000574 0.000821 0.000776] warnings 0
PASS floor_dt             40 s  dE 3.33e-16 dt_gain 0.791 restart_maxdiff 0 retry_maxdiff 0 warnings 0
PASS order_time           25 s  err_E [0.00102 0.000344 0.000113 3.15e-05] order_E [1.57 1.61 1.84] err_F [0.00475 0.00168 0.000588 0.000168] order_F [1.5 1.52 1.81] err_T [0.000167 6.61e-05 2.46e-05 
PASS nlte_eps1e-1         58 s  L1 9.58e-05
PASS nlte_eps1e-4         58 s  L1 0.00133
PASS inf_reject            3 s  dE 3.33e-16 nonfinite 0 rejected 1
PASS wb_star               7 s  v_in 0.00309 v_all 0.672 lat_mid 1.33e-11 warnings 0
PASS idort_shadow          5 s  L1 0.0138
PASS idort_cylinder       11 s  E_in 0.000533 E_out 0.0115 F_out 0.0114 FH_out 0.204
PASS idort_radshock_m2     12 s  L1_rho 0.000356 L1_v 0.000213 L1_Tg 0.000385 L1_Tr 0.00026
PASS idort_thermal_eq2     20 s  errT 0.00196 dEtot 1.88e-16
PASS diffusion            35 s  errD [-0.000815 -7.79e-08 1.42e-10] L1 [0.00137 0.00117 0.00117] warnings 0 errD_tau1 -0.0817 L1_tau1 0.0118 errD_2d -9.29e-07 L1_2d 0.00659 th_min 1.15e-17 th_max -0.0
PASS radshock             26 s  L1_rho 0.00438 L1_v 0.00358 L1_Tg 0.00659 L1_Tr 0.0107 Emin 1.72e+11 warnings 0
PASS marshak              13 s  Emin 8.92e-19 Tmin 9.96e-13 Emax 0.995 Tmax 0.999 Erise 3.07e-18 x_front 0.535 warnings 0
```
### Printed values against rel-a (this suite's common tests)
- rel-b: 40 lines, PASS 40, FAIL 0: []

- against rel-a c880e2703 (02-03): 38 common tests, printed values equal in 16, different in 22, only in rel-b ['star_smoke', 'wb_star']
  - order_space
    - rel-b: PASS err_E [0.0907 0.0208 0.00434] order_E [2.12 2.26] err_v [0.106 0.0379 0.00979] order_v [1.48 1.95]
    - other: PASS err_E [0.071 0.0203 0.00444] order_E [1.81 2.19] err_v [0.12 0.0422 0.011] order_v [1.5 1.94]
  - thermal_eq
    - rel-b: PASS errT_h 0.000552 dEtot_h 1.12e-15 errT_h2 0.000124 dEtot_h2 7e-16 order 2.16
    - other: PASS errT_h 0.000552 dEtot_h 1.26e-15 errT_h2 0.000124 dEtot_h2 1.26e-15 order 2.16
  - radwave
    - rel-b: PASS err_cph [0.00149 0.00147 0.000928 0.000239 0.00226 0.000931] err_damp [0.000865 0.000225 1.63e-05 7.97e-05 0.000655 1.63e-05] damp [0.712 0.164 0.0151 0.000174 0.164 0.0151] warnings 0 dt_gain 2.73 dt_gain_phase [2.16 1.06]
    - other: FAIL err_cph [0.00303 0.00189 0.00137 0.000364 0.00277 0.00137] err_damp [0.00235 0.000105 4.01e-05 0.000208 0.000545 4.02e-05] damp [0.711 0.164 0.0151 0.000302 0.164 0.0151] warnings 0 dt_gain 2.77 dt_gain_phase [2.15 1]  <- err_cph [0.00303 0.00189 0.00137 0.000364 0.00277 0.00137] > 0.003
  - wedge_sym
    - rel-b: PASS fphi [1.23e-09 2.48e-10 1.31e-09] warnings 0 wedge_fatal 1 fth 1.24e-16 theta_fatal 1 fth_rot 1.56e-09 rot_fatal 1
    - other: PASS fphi [3.11e-09 9.79e-12 3.11e-10] warnings 0 wedge_fatal 1 fth 6.91e-15 theta_fatal 1 fth_rot 1.56e-09 rot_fatal 1
  - floor_pay
    - rel-b: PASS dE 3.33e-16 dE0 1.67e-16 paid 0.153 unpaid 0 warnings 0 worst_lines 185 worst_bitwise True ekin_dE 3.33e-16 ekin_dM 5.4e-16 ekin_dLz 1.47e-16
    - other: PASS dE 2.22e-16 dE0 2.22e-16 paid 0.285 unpaid 0 warnings 0 worst_lines 195 worst_bitwise True ekin_dE 1.11e-16 ekin_dM 5.4e-16 ekin_dLz 3.91e-16
  - idort_radshock_m5
    - rel-b: PASS L1_rho 0.00171 L1_v 0.000822 L1_Tg 0.00161 L1_Tr 0.000938
    - other: PASS L1_rho 0.0017 L1_v 0.000818 L1_Tg 0.00162 L1_Tr 0.00094
  - advect_uniform
    - rel-b: PASS dT 5.38e-10 dE 2.68e-08 warnings 0
    - other: PASS dT 2.73e-06 dE 2.35e-05 warnings 0
  - conservation_sp
    - rel-b: PASS dE 3.33e-16 dM 8.09e-16 dLz 1.47e-16 dLz_gr 5.92e-07
    - other: PASS dE 1.11e-16 dM 8.99e-16 dLz 3.91e-16 dLz_gr 6.94e-07
  - conservation_cs
    - rel-b: PASS dE 1.2e-16 dM 2.61e-16 dLz 1.17e-16 dLz_gr 6.48e-05
    - other: PASS dE 1.2e-16 dM 2.61e-16 dLz 1.17e-16 dLz_gr 6.44e-05
  - conservation_lz
    - rel-b: PASS sp_dE 3.33e-16 sp_dM 9.44e-16 sp_dLz 2.94e-16 sp_dLz_gr 1.47e-16 cs_dE 2.4e-16 cs_dM 2.61e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.54e-05 restart_maxdiff 0 retry_maxdiff 0 warnings 0
    - other: PASS sp_dE 2.22e-16 sp_dM 7.19e-16 sp_dLz 1.96e-16 sp_dLz_gr 3.91e-16 cs_dE 1.2e-16 cs_dM 2.61e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.24e-05 restart_maxdiff 0 retry_maxdiff 0 warnings 0
  - floorcap_spin
    - rel-b: PASS sp_dE 3.33e-16 sp_dM 8.09e-16 sp_dLz 1.47e-16 sp_dLz_gr 1.47e-16 n_v 5.85e+04 cs_dE 1.2e-16 cs_dM 1.3e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.3e-05 nospin_sp_dLz_gr 0.366 restart_maxdiff 0 retry_maxdiff 0 layout_maxdiff 0 warnings 0
    - other: PASS sp_dE 2.22e-16 sp_dM 5.4e-16 sp_dLz 1.96e-16 sp_dLz_gr 3.91e-16 n_v 2.59e+04 cs_dE 1.2e-16 cs_dM 1.3e-16 cs_dLz 1.17e-16 cs_dLz_gr 7.8e-05 nospin_sp_dLz_gr 0.288 restart_maxdiff 0 retry_maxdiff 0 layout_maxdiff 0 warnings 0
  - layout_tile
    - rel-b: PASS maxdiff_1x2 0 maxrel_1x2 0 maxdiff_2x2 9.09e-13 maxrel_2x2 1.31e-10 maxdiff_dfl 0 dE_dfl 2.22e-16 dM_dfl 1.8e-16 dLz_dfl 3.92e-16 maxdiff_fofc 0 dE_fofc 1.11e-16 dM_fofc 0.000635 dLz_fofc 1.96e-16 fo_fofc 4.3e+04
    - other: PASS maxdiff_1x2 0 maxrel_1x2 0 maxdiff_2x2 9.09e-13 maxrel_2x2 7.48e-11 maxdiff_dfl 0 dE_dfl 2.22e-16 dM_dfl 1.8e-16 dLz_dfl 1.96e-16 maxdiff_fofc 0 dE_fofc 2.22e-16 dM_fofc 0.000635 dLz_fofc 1.96e-16 fo_fofc 4.3e+04
  - wb_sp
    - rel-b: PASS vmax_over_cs 6.08e-09 vmax_over_cs_rot 6.1e-09
    - other: PASS vmax_over_cs 4.57e-09 vmax_over_cs_rot 4.67e-09
  - wb_cs
    - rel-b: PASS vmax_over_cs 1.36e-09 vmax_over_cs_rot 6.93e-09
    - other: PASS vmax_over_cs 1.18e-09 vmax_over_cs_rot 6.93e-09
  - momentum_box
    - rel-b: PASS dP_max 8.1e-17 dP_cyc [7.3e-17 8.04e-17 7.79e-17 8.02e-17 7.77e-17 7.92e-17 8.1e-17 7.96e-17] warnings 0
    - other: PASS dP_max 5.64e-17 dP_cyc [5.64e-17 5.27e-17 5.23e-17 5.42e-17 5.45e-17 5.21e-17 5.35e-17 5.25e-17] warnings 0
  - radshock_m5
    - rel-b: PASS L1_rho 0.00386 L1_v 0.00212 L1_Tg 0.0062 L1_Tr 0.00744 Emin 1.79e+11 warnings 0
    - other: PASS L1_rho 0.00391 L1_v 0.00212 L1_Tg 0.0062 L1_Tr 0.00744 Emin 1.79e+11 warnings 0
  - advect_pulse
    - rel-b: PASS L1 [0.000872 0.000574 0.000821 0.000776] warnings 0
    - other: PASS L1 [0.000873 0.000575 0.00082 0.000775] warnings 0
  - floor_dt
    - rel-b: PASS dE 3.33e-16 dt_gain 0.791 restart_maxdiff 0 retry_maxdiff 0 warnings 0
    - other: PASS dE 1.11e-16 dt_gain 0.873 restart_maxdiff 0 retry_maxdiff 0 warnings 0
  - inf_reject
    - rel-b: PASS dE 3.33e-16 nonfinite 0 rejected 1
    - other: PASS dE 2.22e-16 nonfinite 0 rejected 1
  - idort_cylinder
    - rel-b: PASS E_in 0.000533 E_out 0.0115 F_out 0.0114 FH_out 0.204
    - other: PASS E_in 0.000531 E_out 0.0111 F_out 0.0116 FH_out 0.203
  - diffusion
    - rel-b: PASS errD [-0.000815 -7.79e-08 1.42e-10] L1 [0.00137 0.00117 0.00117] warnings 0 errD_tau1 -0.0817 L1_tau1 0.0118 errD_2d -9.29e-07 L1_2d 0.00659 th_min 1.15e-17 th_max -0.0392
    - other: PASS errD [-0.000815 -7.79e-08 1.42e-10] L1 [0.00137 0.00117 0.00117] warnings 0 errD_tau1 -0.0817 L1_tau1 0.0118 errD_2d -9.3e-07 L1_2d 0.00659 th_min 1.15e-17 th_max -0.0392
  - marshak
    - rel-b: PASS Emin 8.92e-19 Tmin 9.96e-13 Emax 0.995 Tmax 0.999 Erise 3.07e-18 x_front 0.535 warnings 0
    - other: PASS Emin 8.92e-19 Tmin 9.96e-13 Emax 0.995 Tmax 0.999 Erise 8.89e-19 x_front 0.535 warnings 0
