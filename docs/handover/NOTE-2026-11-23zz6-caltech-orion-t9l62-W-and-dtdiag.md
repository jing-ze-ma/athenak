# NOTE 2026-11-23zz6 caltech -> orion: t9l62 arm W (23zz4) + dt-limiter arms T0/T0d/T1/T1d/T2d/LAY/LAY0 (23zz5)

- W: binary caltech-t9wk-floort-2b1d35094 -> wt_orion1121/bin/athena_t9wkfl_2b1d35094 (md5 in the gate file caltech_results/H200CHK/gate_t9wkfl_2b1d35094.txt).
  Keys: IA121_t9l62_restart.args + clean_wb_floor_diag=8 + clean_wb_mcap + clean_stag_fcone_emax + hydro/floor_from_ekin=true, nlim 8050, from rst_keep/t9l62.00012.
- T/LAY arms: athena_m35fast_13679216a, env ATHENA_DT_REPORT=5, IA121_t9l62_restart.args + the arm key, rst/bin dt 1e30, job/final_restart=false.
  rst and nlim (absolute): T0/T0d prod_fresh_1121b/t9l62/rst_p1_seg/t9l62.00037.rst (00031 no longer exists; 00037 = chunk-5 final rst, cycle 30085, t ~1.574e6), nlim 30115;
  T1/T1d rst_keep/00010 (cycle 4514, t 4.31e5), nlim 4544; T2d rst_keep/00012 (cycle 6457), nlim 6487; LAY/LAY0 rst_keep/00010, nlim 4614.
- One arm per node (t9l62 needs ~104 GB/GPU on 4 H200, so two per node do not fit). Test stop rule 23zk inside each run script (t9q_1123/run_T.sh).
- W/SB continuation with ATHENA_DT_REPORT=10: skipped. No rst is kept, so there is no state to continue from.
- Files: caltech_results/IA121_T9L62/q_<arm>/.

## IA11 production t9l62 output write times (23zz5 "report back")
- rst (88.1 GB): `<restart> rst/t9l62.00010.rst: all (file closed) written at 92.3 s`; 00011 at 85.9 s.
  Per-part times for 00011: cell/face vars 58.4 s, hesdirk2 slope 83.3 s, w0 84.1 s, rad_signal_speed 85.5 s, face flux 85.9 s.
- The bins print no timing line. Estimated from the file mtimes of one dump (11:17, output index 12/3): m1_face (1.07 GB) closed at 26.81 s, m1_vet (1.61 GB) +1.77 s,
  m1 (1.88 GB) +2.08 s, hydro_w (1.34 GB) +1.54 s. So each bin takes ~1.5-2 s; m1_face's own start is unknown.
- (11-23 21:03) status: WB [cycle=7550, t 4.85395e+05 dt 12.2, CW 28, fb 16, DT_RETRY 12, pos_floor 1]; SB [cycle=7000, t 4.79148e+05 dt 9.31, CW 3, fb 2, DT_RETRY 2, pos_floor 0];
- (11-23 21:33) status: WB [cycle=8000, t 4.89076e+05 dt 3.16, CW 48, fb 28, DT_RETRY 21, pos_floor 1]; SB [cycle=7650, t 4.86040e+05 dt 11.9, CW 3, fb 2, DT_RETRY 2, pos_floor 0]; W [cycle=6750, t 4.75700e+05 dt 4.76, CW 38, fb 20, DT_RETRY 21, pos_floor 0];
- (11-23 21:38) status:

## T0 (11-23 21:38)
- 2026-11-23 21:37:34 rc=0 last: elapsed=7.001133e+01 cycle=30100 time=1.579165e+06 dt=1.173417e+02 
- FATAL 0 
- CW 0, stage fallbacks 0, DT_RETRY 0, pos_floor created 0; s/cycle n/a; dt min 95.29 at t 1.580408e+06; skip_de per cycle: max abs 1.238e+17 at c30105; 0 cycles > 1e41
- first DTREP: `DTREP cycle=30085 time=1.57756e+06 cfl*dt_hyd=107.368 cfl*dt_gasonly=187.657 rank=2 gid=5 kji=34,3,146 dir=1 r=4.5059e+12 rho=2.89929e-11 T=1.67374e+13 pg=485.265 G1=1.66667 cgas=5.28164e+06 v=-3.96447e+07 ceff=1.32938e+08 h=3.08832e+10 tauT_h=0.595814 Erad=2.56419e+06 Ddd=0.333333 cdt_kP=5083.75 frz=-1 sc=1 cell_dt=178.947 cell_dt_gas=687.417 theta_deg=90.1758 phi=1.76408`

## LAY (11-23 21:43)
- 2026-11-23 21:41:55 rc=0 last:  
- FATAL 1 ### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/v203/src/main.cpp: a key that places the cells differs from the restart file,     meshblock/nx2: restart file 128 -> now 256 
- CW 0, stage fallbacks 0, DT_RETRY 0, pos_floor created 0; s/cycle n/a; dt min ; skip_de per cycle: max abs 0.000e+00 at c; 0 cycles > 1e41

## T0d (11-23 21:43)
- 2026-11-23 21:40:47 rc=0 last: elapsed=7.561569e+01 cycle=30100 time=1.580358e+06 dt=1.927452e+02 
- FATAL 0 
- CW 0, stage fallbacks 0, DT_RETRY 0, pos_floor created 0; s/cycle n/a; dt min 192.7 at t 1.580551e+06; skip_de per cycle: max abs 1.983e+17 at c30114; 0 cycles > 1e41
- first DTREP: `DTREP cycle=30085 time=1.57756e+06 cfl*dt_hyd=187.496 cfl*dt_gasonly=187.657 rank=3 gid=6 kji=86,77,144 dir=1 r=4.44456e+12 rho=6.81208e-13 T=5.16062e+12 pg=3.51546 G1=1.66667 cgas=2.93275e+06 v=-9.44663e+07 ceff=3.01646e+06 h=3.04627e+10 tauT_h=0.00803132 Erad=23168.4 Ddd=0.333333 cdt_kP=595.969 frz=-1 sc=1 cell_dt=312.493 cell_dt_gas=312.762 theta_deg=71.1915 phi=2.86854`
- (11-23 21:49) status:

## LAY0 (11-23 21:53)
- 2026-11-23 21:49:08 rc=0 last: elapsed=3.117645e+02 cycle=4600 time=4.379579e+05 dt=6.829116e+01 
- FATAL 0 
- CW 0, stage fallbacks 0, DT_RETRY 0, pos_floor created 0; s/cycle 3.072 (c4550-4600); dt min 68.07 at t 4.389123e+05; skip_de per cycle: max abs 1.792e+24 at c4604; 0 cycles > 1e41
- first DTREP: `DTREP cycle=4515 time=432085 cfl*dt_hyd=70.0239 cfl*dt_gasonly=87.2318 rank=3 gid=7 kji=22,16,49 dir=1 r=2.31764e+12 rho=1.69288e-07 T=4.74913e+13 pg=8.03971e+06 G1=1.66667 cgas=8.89675e+06 v=-4.54801e+07 ceff=2.25751e+07 h=1.58849e+10 tauT_h=1291.34 Erad=1.6397e+08 Ddd=0.333333 cdt_kP=1.6296e+06 frz=-1 sc=1 cell_dt=233.413 cell_dt_gas=292.127 theta_deg=94.7461 phi=2.47584`
- (11-23 21:54) status: T0a [cycle=30150, t 1.58487e+06 dt 104, CW 0, fb 0, DT_RETRY 0, pos_floor 0];
- (11-23 21:54) status: T0a [cycle=30150, t 1.58508e+06 dt 98.3, CW 0, fb 0, DT_RETRY 0, pos_floor 0];
- (11-23 21:54) status: SB [cycle=8000, t 4.89336e+05 dt 3.25, CW 12, fb 5, DT_RETRY 2, pos_floor 2]; W [cycle=7150, t 4.79349e+05 dt 9.96, CW 38, fb 20, DT_RETRY 21, pos_floor 0];
- (11-23 21:55) **reply to 23zz12:** the empty 'status:' lines (21:38, 21:49) came from the LAY and 23zz11 posters, which had no running arm yet; they are fixed. W was in the main status line (21:33).
  **W now:** cycle 7150, t 4.792994e5, dt 9.79, Ltop/Lin 2.510, Lbot/Lin 0.287, CW 38, stage fallbacks 20, DT_RETRY 21, pos_floor created 0.
  floor_audit [update] skip_de per cycle: max abs 1.964e43 at c6626; **705 of 717 cycles > 1e41** (the fix does not bring skip_de down to ~1e38).
  It has NOT yet passed t 4.84e5 or 4.89e5 (at 3.5-4 s/cycle and dt ~10, ~4.84e5 comes in ~25 min).

## (11-23 22:05) reply to 23zz14: arm W evidence (W still running to 8050; now cycle 7350, t 4.8109e5, dt 11.7)
- run.log start line: `2026-11-23 21:03:51 start W BIN=athena_t9wkfl_f18bae9fa md5=30affa062badfd6025da515d33b7fc11 NODE=hpc-sm-02-17 HOLDER=7333073 ATHENA_DT_REPORT=unset RST=/resnick/groups/carnegie_poc/jingze/prod_ia11_1 ...`. The full key line is in caltech_results/IA121_T9L62/arm_W/run.log. It contains **hydro/floor_from_ekin=true**, rad_m1/clean_wb_floor_diag=8, clean_wb_mcap=true and clean_stag_fcone_emax=true. ATHENA_DT_REPORT is unset.
- Binary: **athena_t9wkfl_f18bae9fa md5 30affa062badfd6025da515d33b7fc11** (caltech-t9wk-floort-f18bae9fa, not 2b1d35094). Gate: `GATE PASS: athena_t9wkfl_f18bae9fa md5 30affa06 node hpc-sm-02-17 pass1=3.686e-11 (ref 3.686e-11, 3 passes, F/c 2.9006e+42)`.
- Echo in the log: line 24 `### NOTE: parameter 'floor_from_ekin' in block 'hydro' is absent from the restart's input: added from the command line (= true)`.
  The <hydro> startup line (`duale-1112: dual_energy 0 ... floor_pay rad_column fofc_retest 0 floor_audit 1`) does not print floor_from_ekin.
  The key is NOT in the 'never read' list, which has only job/final_restart, so the code read it.
- floor_audit lines (all three kinds: [update], [refresh], pay) for 3 cycles with [update] skip_de > 1e42 (6475, 6476, plus 6626, the max): caltech_results/IA121_T9L62/arm_W/floor_audit_c6475_6476_6626.txt. Inline:
```
<hydro> floor_audit cycle 6475 [update] cell-c2p 2483027968 e_c<=0 83 kedom 97 sub_tfloor 83 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 83 (kedom 83) floor_de 9.953766291e+40 <r> 2.443026818e+12 <r>_de 2.31967354e+12 skip(dfloor/vceil) 631575 skip_de 1.093827798e+42 skip_dm 0.8681114721 above_table_top 0
<hydro> floor_audit cycle 6475 [refresh] cell-c2p 3690987520 e_c<=0 7 kedom 132 sub_tfloor 12185 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 7 (kedom 7) floor_de 1.262945755e+38 <r> 2.419203524e+12 <r>_de 2.304716509e+12 skip(dfloor/vceil) 1071941 skip_de 1.183732641e+24 skip_dm 0 above_table_top 0
<hydro> floor_audit cycle 6475 pay: paid 1.192340414e+42 unpaid 0 given_to_rad 4.772493289e+11 cells 41508 | fofc_retest: retested 0 unflagged 0
<hydro> floor_audit cycle 6476 [update] cell-c2p 2684354560 e_c<=0 88 kedom 102 sub_tfloor 88 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 88 (kedom 88) floor_de 9.957881639e+40 <r> 2.443933997e+12 <r>_de 2.319704448e+12 skip(dfloor/vceil) 682767 skip_de 1.094978764e+42 skip_dm 0.8681114721 above_table_top 0
<hydro> floor_audit cycle 6476 [refresh] cell-c2p 3892314112 e_c<=0 8 kedom 136 sub_tfloor 12186 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 8 (kedom 8) floor_de 1.272969016e+38 <r> 2.431343035e+12 <r>_de 2.306382656e+12 skip(dfloor/vceil) 1110391 skip_de 1.183732641e+24 skip_dm 0 above_table_top 0
<hydro> floor_audit cycle 6476 pay: paid 1.193491755e+42 unpaid 0 given_to_rad 4.772493289e+11 cells 43525 | fofc_retest: retested 0 unflagged 0
<hydro> floor_audit cycle 6626 [update] cell-c2p 2.281701376e+10 e_c<=0 361 kedom 464 sub_tfloor 361 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 361 (kedom 361) floor_de 3.335832758e+41 <r> 2.443100868e+12 <r>_de 2.322210763e+12 skip(dfloor/vceil) 5781299 skip_de 2.806536063e+42 skip_dm 3.034017315 above_table_top 0
<hydro> floor_audit cycle 6626 [refresh] cell-c2p 3.395708518e+10 e_c<=0 13 kedom 534 sub_tfloor 12191 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 13 (kedom 13) floor_de 1.443068133e+38 <r> 2.414168565e+12 <r>_de 2.309179213e+12 skip(dfloor/vceil) 9213837 skip_de 2.955580135e+24 skip_dm 0 above_table_top 0
<hydro> floor_audit cycle 6626 pay: paid 3.138571693e+42 unpaid 0 given_to_rad 1.667968658e+12 cells 335770 | fofc_retest: retested 0 unflagged 0
<hydro> floor_audit cycle 6626 [update] cell-c2p 2.281701376e+10 e_c<=0 361 kedom 464 sub_tfloor 361 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 361 (kedom 361) floor_de 3.334763685e+41 <r> 2.44324273e+12 <r>_de 2.322204376e+12 skip(dfloor/vceil) 5781236 skip_de 2.806368674e+42 skip_dm 3.034017315 above_table_top 0
<hydro> floor_audit cycle 6626 [refresh] cell-c2p 3.395708518e+10 e_c<=0 13 kedom 534 sub_tfloor 12191 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 13 (kedom 13) floor_de 1.443068133e+38 <r> 2.414168565e+12 <r>_de 2.309179213e+12 skip(dfloor/vceil) 9213837 skip_de 2.955580135e+24 skip_dm 0 above_table_top 0
<hydro> floor_audit cycle 6626 pay: paid 3.138571693e+42 unpaid 0 given_to_rad 1.667968658e+12 cells 335770 | fofc_retest: retested 0 unflagged 0
<hydro> floor_audit cycle 6626 [update] cell-c2p 2.281701376e+10 e_c<=0 361 kedom 463 sub_tfloor 361 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 361 (kedom 361) floor_de 3.333424319e+41 <r> 2.44324273e+12 <r>_de 2.322153564e+12 skip(dfloor/vceil) 5781195 skip_de 2.806368674e+42 skip_dm 3.034017315 above_table_top 0
<hydro> floor_audit cycle 6626 [refresh] cell-c2p 3.395708518e+10 e_c<=0 13 kedom 534 sub_tfloor 12191 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 13 (kedom 13) floor_de 1.443068133e+38 <r> 2.414168565e+12 <r>_de 2.309179213e+12 skip(dfloor/vceil) 9213837 skip_de 2.955580135e+24 skip_dm 0 above_table_top 0
<hydro> floor_audit cycle 6626 pay: paid 3.138571693e+42 unpaid 0 given_to_rad 1.667968658e+12 cells 335770 | fofc_retest: retested 0 unflagged 0
<hydro> floor_audit cycle 6626 [update] cell-c2p 2.281701376e+10 e_c<=0 360 kedom 462 sub_tfloor 360 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 360 (kedom 360) floor_de 3.332342914e+41 <r> 2.443321357e+12 <r>_de 2.322116715e+12 skip(dfloor/vceil) 5781184 skip_de 2.806368674e+42 skip_dm 3.034017315 above_table_top 0
<hydro> floor_audit cycle 6626 [refresh] cell-c2p 3.395708518e+10 e_c<=0 13 kedom 534 sub_tfloor 12191 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 13 (kedom 13) floor_de 1.443068133e+38 <r> 2.414168565e+12 <r>_de 2.309179213e+12 skip(dfloor/vceil) 9213837 skip_de 2.955580135e+24 skip_dm 0 above_table_top 0
<hydro> floor_audit cycle 6626 pay: paid 3.138571693e+42 unpaid 0 given_to_rad 1.667968658e+12 cells 335770 | fofc_retest: retested 0 unflagged 0
<hydro> floor_audit cycle 6626 [update] cell-c2p 2.281701376e+10 e_c<=0 359 kedom 461 sub_tfloor 359 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 359 (kedom 359) floor_de 3.332214486e+41 <r> 2.443260672e+12 <r>_de 2.3221096e+12 skip(dfloor/vceil) 5781169 skip_de 2.806368674e+42 skip_dm 3.034017315 above_table_top 0
<hydro> floor_audit cycle 6626 [refresh] cell-c2p 3.395708518e+10 e_c<=0 13 kedom 534 sub_tfloor 12191 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 13 (kedom 13) floor_de 1.443068133e+38 <r> 2.414168565e+12 <r>_de 2.309179213e+12 skip(dfloor/vceil) 9213837 skip_de 2.955580135e+24 skip_dm 0 above_table_top 0
<hydro> floor_audit cycle 6626 pay: paid 3.138571693e+42 unpaid 0 given_to_rad 1.667968658e+12 cells 335770 | fofc_retest: retested 0 unflagged 0
<hydro> floor_audit cycle 6626 [update] cell-c2p 2.281701376e+10 e_c<=0 358 kedom 461 sub_tfloor 358 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 358 (kedom 358) floor_de 3.332110559e+41 <r> 2.443056596e+12 <r>_de 2.322103543e+12 skip(dfloor/vceil) 5781161 skip_de 2.806368674e+42 skip_dm 3.034017315 above_table_top 0
<hydro> floor_audit cycle 6626 [refresh] cell-c2p 3.395708518e+10 e_c<=0 13 kedom 534 sub_tfloor 12191 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 13 (kedom 13) floor_de 1.443068133e+38 <r> 2.414168565e+12 <r>_de 2.309179213e+12 skip(dfloor/vceil) 9213837 skip_de 2.955580135e+24 skip_dm 0 above_table_top 0
<hydro> floor_audit cycle 6626 pay: paid 3.138571693e+42 unpaid 0 given_to_rad 1.667968658e+12 cells 335770 | fofc_retest: retested 0 unflagged 0
<hydro> floor_audit cycle 6626 [update] cell-c2p 2.281701376e+10 e_c<=0 358 kedom 460 sub_tfloor 358 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 358 (kedom 358) floor_de 3.332110559e+41 <r> 2.443056596e+12 <r>_de 2.322103543e+12 skip(dfloor/vceil) 5781158 skip_de 2.806368674e+42 skip_dm 3.034017315 above_table_top 0
<hydro> floor_audit cycle 6626 [refresh] cell-c2p 3.395708518e+10 e_c<=0 13 kedom 534 sub_tfloor 12191 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 13 (kedom 13) floor_de 1.443068133e+38 <r> 2.414168565e+12 <r>_de 2.309179213e+12 skip(dfloor/vceil) 9213837 skip_de 2.955580135e+24 skip_dm 0 above_table_top 0
<hydro> floor_audit cycle 6626 pay: paid 3.138571693e+42 unpaid 0 given_to_rad 1.667968658e+12 cells 335770 | fofc_retest: retested 0 unflagged 0
```
- (11-23 22:24) status: T0a [cycle=30650, t 1.63707e+06 dt 82.3, CW 0, fb 0, DT_RETRY 0, pos_floor 0]; T0da [cycle=30500, t 1.65337e+06 dt 193, CW 2, fb 1, DT_RETRY 0, pos_floor 1]; T0h [cycle=30550, t 1.59182e+06 dt 29.8, CW 0, fb 0, DT_RETRY 0, pos_floor 0];
- (11-23 22:24) status: W [cycle=7650, t 4.82933e+05 dt 5.78, CW 59, fb 29, DT_RETRY 27, pos_floor 2];

## W (11-23 22:35)
- 2026-11-23 22:31:05 rc=0 last: elapsed=5.033018e+03 cycle=7700 time=4.831070e+05 dt=5.585144e+00 
- FATAL 0 
- CW 76, stage fallbacks 33, DT_RETRY 29, pos_floor created 5; s/cycle 3.774 (c6500-7700); dt min 0.02278 at t 4.752742e+05; skip_de per cycle: max abs 2.534e+43 at c7711; 1267 cycles > 1e41

## (11-23 22:35) W stopped (23zz18) + skip_de correction (23zz19)
- **W stopped** at 22:31 (scancel of the test step, per 23zz18). Last: cycle 7700, t 4.831070e5, dt 5.59. Its final section is above; its log + hst are in q_W/. The skip_de line in W's final section above uses the OLD cumulative method: ignore it, see below.
- **skip_de_per_cycle.txt fixed** (skipde.sh in caltech_results/IA121_T9L62/): now the increment of the cumulative [update] skip_de between the LAST [update] line of consecutive cycles (the accepted pass; redos are not summed). Regenerated for arm_WA, arm_WB, arm_SB, q_W, arm_S, arm_L:
  - WA: 1053 cycles, 6 > 1e41, max 1.756e+42 at c6470; 47 negative increments (DT_RETRY restores / sign of skip_de)
  - WB: 1592 cycles, 8 > 1e41, max 1.216e+42 at c6467; 113 negative increments (DT_RETRY restores / sign of skip_de)
  - SB: 1592 cycles, 16 > 1e41, max 8.816e+41 at c6463; 111 negative increments (DT_RETRY restores / sign of skip_de)
  - W: 1278 cycles, 18 > 1e41, max 1.206e+42 at c7711; 0 negative increments (DT_RETRY restores / sign of skip_de)
  - S: 392 cycles, 6 > 1e41, max 8.816e+41 at c6463; 154 negative increments (DT_RETRY restores / sign of skip_de)
  - L: 392 cycles, 5 > 1e41, max 1.749e+42 at c6470; 138 negative increments (DT_RETRY restores / sign of skip_de)
  These agree with T9FLOOR's table to within 1-2 cycles (WB 8 vs 9, WA 6 vs 8, SB 16 vs 16). The difference is probably whether a redone cycle's last line or its first accepted line is taken.
