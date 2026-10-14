# NOTE Caltech (DR, e1div-1014 part 3): fin4-dtretry build + full-size H200 gate PASS; D1 running

## Build
- origin/fin4-dtretry = 5b7e1bc5 confirmed. Incremental nofma build (build_inc_nofma.sh, from the d505d21f tree, 197 objects):
  `builds/bin/athena_gpu_he_star_m1_5b7e1bc5_nofma` md5 **6252e58e5192ce65ff1f1f627e1ce90c**.
- Reference R binary d505d21f md5 1440c05320a7bef88c4ad960019b9840.

## Gate (job 4491313, hpc-sm-01-05, 1 node x 2 H200, fresh, nlim 100, all rc=0)
Deck e1_1013/agcA_E1_final_wbrmax.athinput (md5 d29e4615); F22/F88 copies made with the TASK awk (md5 9dd1428f / 359dea09).
Keys: `time/nlim=100 time/ndiag=1 output2/dt=1.5e4 output3/dt=1.5e4 output4/dt=1e30 output5/dt=1e30`. Dir e1div_dtr/gate/{R,N,F22,F88}.

| compare | hst rows | bins (hydro_w + m1) | verdict |
|---|---|---|---|
| R vs N   | 8+8 equal | 12/12 equal | CMPRUNS BITWISE |
| R vs F22 | 8+8 equal | 12/12 equal | CMPRUNS BITWISE |
| R vs F88 | 8+8 equal | 12/12 equal | CMPRUNS BITWISE |

- F22: one `DT_RETRY cycle=22 t=1.600000e+04 dt 7.272727e+02 -> 7.272727e+02 reason=forced resid=0 cut 1`; F88 the same at cycle 88 (t=6.4e4).
  No `recovered` line printed for the forced fac-1.0 retry (only the event line); end lines `100 checkpoints (~0.59 s, 4544 MB on rank 0), 1 retries, cap -1`.
- N: `dt_retry on: fac 0.5 max 8 ngood 10 grow 1.25 views minimal cap -1`; 0 retries; checkpoint 1237 members / 333 Views (15 copied) 4544 MB on rank 0
  (largest m1.vra_d 2248 MB, m1.tau_ten 425 MB, m1.t2k1 340 MB).
- All arms: NON-CONVERGED 0, Picard mean 3.07 max 6, no FATAL/DIVERGED, floor clips 0.

## Memory (nvidia-smi peak, 500 ms sampling, per GPU = per rank)
R 55165 MiB; N 59903 MiB; F22 59885 MiB; F88 60413 MiB. **N - R = +4738 MiB/rank (~4.6 GB)**, i.e. about 2x the expected ~2.3 GB/rank,
consistent with the logged 4544 MB checkpoint on rank 0. Fits H200 141 GB with wide margin.

## s/cycle (elapsed, cycles 10-100)
R 0.5686; N 0.5704 (+0.3 %); F22 0.5780; F88 0.5830 (one extra cycle each, as expected). Init ~40 s/arm; wall 104-105 s/arm.

## D1 (target: E1 infall from rst 00022 with dt_retry)
- Job **4492078** (hpc-sm-02-04, started 2026-10-13 22:29:41 PDT, -t 06:00:00, athena -t 05:40:00, --exclude=hpc-sm-01-09).
- Dir e1div_rst/D1, rst/agcar3d.00022.rst hard-linked from R0/rst. Script e1div_rst/D1.sh = R0's exact keys (K + time/tlim=1.26e7 +
  output2/dt=2.5e4 output3/dt=2.5e4) with the 5b7e1bc5 binary, then `time/dt_retry=true time/tlim=1.30e7` last (overrides R0's 1.26e7).
- D1 report (DT_RETRY lines, Picard, dt trace, prof_*) follows in a later NOTE.

## Part 2
- R0/R1 (part 2 restart arms) were cancelled before starting (user, to free GPUs for the Linhao BSG). Their dirs e1div_rst/R0, R1 are kept
  (R0/rst/agcar3d.00022.rst also feeds D1). No part-2 reports.

## D1 result (job 4492078, hpc-sm-02-04, 2 H200; 2026-10-13 22:29:41 -> 2026-10-14 00:59:46, elapsed 02:30:05)
- **rc 0, COMPLETED, reached tlim**: `Terminating on time limit`, end `time=1.300000e+07 cycle=18676` (start rst 00022: t 1.10e7, cycle 15125;
  3551 cycles). Dumps bin/ hydro_w + m1 00111-00190 (80 pairs), rst 00023. No FATAL anywhere in run.log.
- **dt_retry end line**: `<time> dt_retry: 3551 checkpoints (2.104529e+01 s, 4.926486e+03 MB on rank 0), 116 retries, cap -1.000000e+00`
  (checkpoint overhead 21 s of ~8960 s = 0.23 %).
- **DT_RETRY lines: 199** = 116 events + 83 `recovered`. Reasons: 106 NON-CONVERGED, 10 DIVERGED. Cuts: 113 `cut 1`, 3 `cut 2`; never more
  than 2 cuts in a cycle (max 8 allowed). Several events fall inside a still-active cap (hence 83 recovered < 116). First event cycle 15198
  (t 1.1053e7); events cluster in cycles ~15900-18450 (t 1.155e7-1.284e7), last at 18661.
- **FATAL/DIVERGED (10 lines, all DT_RETRY events, all `cut 1`, all recovered by the retry; no FATAL)**:
```
DT_RETRY cycle=17521 t=1.241818e+07 dt 7.272727e+02 -> 3.636364e+02 reason=DIVERGED resid=2.945800e-02 cut 1
DT_RETRY cycle=17684 t=1.248928e+07 dt 6.717236e+02 -> 3.358618e+02 reason=DIVERGED resid=2.671058e-01 cut 1
DT_RETRY cycle=17876 t=1.256812e+07 dt 6.568749e+02 -> 3.284374e+02 reason=DIVERGED resid=1.824288e-01 cut 1
DT_RETRY cycle=17935 t=1.259345e+07 dt 7.272727e+02 -> 3.636364e+02 reason=DIVERGED resid=4.308236e-02 cut 1
DT_RETRY cycle=17953 t=1.259788e+07 dt 5.621945e+02 -> 2.810972e+02 reason=DIVERGED resid=1.433166e-02 cut 1
DT_RETRY cycle=17998 t=1.261164e+07 dt 7.272727e+02 -> 3.636364e+02 reason=DIVERGED resid=4.687301e-02 cut 1
DT_RETRY cycle=18023 t=1.261937e+07 dt 6.613023e+02 -> 3.306511e+02 reason=DIVERGED resid=1.130167e+00 cut 1
DT_RETRY cycle=18056 t=1.263334e+07 dt 6.660633e+02 -> 3.330317e+02 reason=DIVERGED resid=1.151805e-02 cut 1
DT_RETRY cycle=18161 t=1.267782e+07 dt 7.272727e+02 -> 3.636364e+02 reason=DIVERGED resid=8.422926e-02 cut 1
DT_RETRY cycle=18197 t=1.269891e+07 dt 7.272727e+02 -> 3.636364e+02 reason=DIVERGED resid=5.878629e-01 cut 1
```
  "DIVERGED" here = Picard still at 200 passes with resid 1e-2..1.1 (e.g. cycle 17521: `Picard NON-CONVERGED after 200 passes: resid=2.945800e-02`,
  NC-local rank 1 worst cell (0,17,46,348), rank 0 (0,66,47,506)); it is the dt_retry reason class, not a crash. All 10 are in t 1.242e7-1.270e7.
- **Non-converged / Picard**: 116 `Picard NON-CONVERGED after 200 passes` lines (= the 116 retry events; each one retried). End summary
  `implicit transport: solves=7.233000e+03 Picard iterations mean=3.570054e+01 max=2.000000e+02 NON-CONVERGED=0.000000e+00` (no
  unrecovered NC). hesdirk2: stage steps 2321, BE steps 1230 (stage fallbacks 1230); `NOT ADMISSIBLE` lines 1349. Worst-cell x1 index i (rank 0, 116 NC lines), per 10: 320s 3, 330s 7, 340s 19, 350s 13, 360s 7,
  370s 6, 380s 4, 390s 9, 400s 10, 410s 7, 420s 13, 430-530s 18 (mostly i 330-430; active radial index 0..535). Floor clips: E<=e_floor 76, written-back eint<=0 1947 cell-solves.
- **s/cycle** 2.531 (elapsed 70.96 s @ cycle 15150 -> 8929.95 s @ 18650; includes retries). Mean dt 2.0e6/3551 = 563 s (nominal CFL 727).
- **First 20 DT_RETRY lines**:
```
DT_RETRY cycle=15198 t=1.105309e+07 dt 7.272727e+02 -> 3.636364e+02 reason=NON-CONVERGED resid=4.071500e-07 cut 1
DT_RETRY recovered cycle=15212 t=1.105881e+07 cap 8.877841e+02 >= CFL dt 7.220607e+02: cap off
DT_RETRY cycle=15237 t=1.107709e+07 dt 7.272727e+02 -> 3.636364e+02 reason=NON-CONVERGED resid=6.366008e-06 cut 1
DT_RETRY recovered cycle=15251 t=1.108275e+07 cap 8.877841e+02 >= CFL dt 7.219867e+02: cap off
DT_RETRY cycle=15507 t=1.126836e+07 dt 7.272727e+02 -> 3.636364e+02 reason=NON-CONVERGED resid=1.635584e-07 cut 1
DT_RETRY recovered cycle=15521 t=1.127403e+07 cap 8.877841e+02 >= CFL dt 7.215909e+02: cap off
DT_RETRY cycle=15660 t=1.137455e+07 dt 7.272727e+02 -> 3.636364e+02 reason=NON-CONVERGED resid=2.484161e-08 cut 1
DT_RETRY recovered cycle=15674 t=1.138025e+07 cap 8.877841e+02 >= CFL dt 7.216125e+02: cap off
DT_RETRY cycle=15739 t=1.142764e+07 dt 7.272727e+02 -> 3.636364e+02 reason=NON-CONVERGED resid=5.517157e-07 cut 1
DT_RETRY recovered cycle=15753 t=1.143331e+07 cap 8.877841e+02 >= CFL dt 7.216959e+02: cap off
DT_RETRY cycle=15757 t=1.143599e+07 dt 6.688622e+02 -> 3.344311e+02 reason=NON-CONVERGED resid=1.117437e-08 cut 1
DT_RETRY recovered cycle=15771 t=1.144116e+07 cap 8.164822e+02 >= CFL dt 7.217147e+02: cap off
DT_RETRY cycle=15810 t=1.146909e+07 dt 7.272727e+02 -> 3.636364e+02 reason=NON-CONVERGED resid=1.728598e-07 cut 1
DT_RETRY cycle=15810 t=1.146909e+07 dt 3.636364e+02 -> 1.818182e+02 reason=NON-CONVERGED resid=3.090606e-08 cut 2
DT_RETRY recovered cycle=15825 t=1.147200e+07 cap 5.548651e+02 >= CFL dt 4.545455e+02: cap off
DT_RETRY cycle=15941 t=1.155564e+07 dt 7.272727e+02 -> 3.636364e+02 reason=NON-CONVERGED resid=1.610490e-04 cut 1
DT_RETRY recovered cycle=15955 t=1.156131e+07 cap 8.877841e+02 >= CFL dt 7.220844e+02: cap off
DT_RETRY cycle=15973 t=1.157382e+07 dt 7.272727e+02 -> 3.636364e+02 reason=NON-CONVERGED resid=6.532432e-08 cut 1
DT_RETRY recovered cycle=15987 t=1.157948e+07 cap 8.877841e+02 >= CFL dt 7.220963e+02: cap off
DT_RETRY cycle=15994 t=1.158400e+07 dt 7.272727e+02 -> 3.636364e+02 reason=NON-CONVERGED resid=9.411054e-04 cut 1
```
- **Last 20 DT_RETRY lines**:
```
DT_RETRY cycle=18132 t=1.266168e+07 dt 7.368421e+02 -> 3.684211e+02 reason=NON-CONVERGED resid=3.172795e-04 cut 1
DT_RETRY recovered cycle=18146 t=1.266742e+07 cap 8.994655e+02 >= CFL dt 7.231789e+02: cap off
DT_RETRY cycle=18161 t=1.267782e+07 dt 7.272727e+02 -> 3.636364e+02 reason=DIVERGED resid=8.422926e-02 cut 1
DT_RETRY recovered cycle=18175 t=1.268348e+07 cap 8.877841e+02 >= CFL dt 7.232292e+02: cap off
DT_RETRY cycle=18197 t=1.269891e+07 dt 7.272727e+02 -> 3.636364e+02 reason=DIVERGED resid=5.878629e-01 cut 1
DT_RETRY cycle=18210 t=1.270364e+07 dt 3.636364e+02 -> 1.818182e+02 reason=NON-CONVERGED resid=1.055507e-07 cut 1
DT_RETRY recovered cycle=18227 t=1.270821e+07 cap 8.669767e+02 >= CFL dt 7.232867e+02: cap off
DT_RETRY cycle=18249 t=1.272436e+07 dt 7.272727e+02 -> 3.636364e+02 reason=NON-CONVERGED resid=2.370015e-05 cut 1
DT_RETRY recovered cycle=18263 t=1.273003e+07 cap 8.877841e+02 >= CFL dt 7.233346e+02: cap off
DT_RETRY cycle=18270 t=1.273467e+07 dt 6.634307e+02 -> 3.317153e+02 reason=NON-CONVERGED resid=1.371836e-07 cut 1
DT_RETRY recovered cycle=18284 t=1.273989e+07 cap 8.098519e+02 >= CFL dt 7.233580e+02: cap off
DT_RETRY cycle=18355 t=1.279127e+07 dt 7.272727e+02 -> 3.636364e+02 reason=NON-CONVERGED resid=2.193632e-05 cut 1
DT_RETRY cycle=18367 t=1.279570e+07 dt 5.371901e+02 -> 2.685950e+02 reason=NON-CONVERGED resid=3.773606e-08 cut 1
DT_RETRY recovered cycle=18382 t=1.280000e+07 cap 8.196870e+02 >= CFL dt 6.714876e+02: cap off
DT_RETRY cycle=18426 t=1.283200e+07 dt 7.272727e+02 -> 3.636364e+02 reason=NON-CONVERGED resid=3.000592e-06 cut 1
DT_RETRY recovered cycle=18440 t=1.283772e+07 cap 8.877841e+02 >= CFL dt 7.233142e+02: cap off
DT_RETRY cycle=18610 t=1.296145e+07 dt 7.272727e+02 -> 3.636364e+02 reason=NON-CONVERGED resid=1.572155e-08 cut 1
DT_RETRY recovered cycle=18624 t=1.296715e+07 cap 8.877841e+02 >= CFL dt 7.231540e+02: cap off
DT_RETRY cycle=18661 t=1.299418e+07 dt 7.272727e+02 -> 3.636364e+02 reason=NON-CONVERGED resid=2.297981e-07 cut 1
DT_RETRY recovered cycle=18675 t=1.299986e+07 cap 8.877841e+02 >= CFL dt 7.231313e+02: cap off
```
- **dt every 100 cycles** (`cycle t dt elapsed_s`):
```
15200 1.105382e+07 3.636364e+02 163.1
15300 1.111782e+07 7.272727e+02 357.9
15400 1.119055e+07 7.272727e+02 542.8
15500 1.126327e+07 7.272727e+02 754.2
15600 1.133091e+07 7.272727e+02 963.3
15700 1.139927e+07 7.272727e+02 1174.1
15800 1.146182e+07 7.272727e+02 1409.9
15900 1.152582e+07 7.272727e+02 1659.6
16000 1.158618e+07 3.636364e+02 1909.6
16100 1.164545e+07 3.636364e+02 2155.5
16200 1.170541e+07 7.325032e+02 2398.8
16300 1.175933e+07 3.358586e+02 2669.3
16400 1.182182e+07 7.272727e+02 2938.3
16500 1.187131e+07 1.715572e+02 3185.7
16600 1.192764e+07 3.636364e+02 3460.6
16700 1.198618e+07 3.636364e+02 3745.1
16800 1.203745e+07 3.636364e+02 4004.8
16900 1.208800e+07 7.272727e+02 4261.0
17000 1.214545e+07 7.272727e+02 4571.2
17100 1.220393e+07 6.785087e+02 4868.8
17200 1.226343e+07 7.145322e+02 5179.2
17300 1.231055e+07 7.272727e+02 5425.1
17400 1.235642e+07 7.237495e+02 5671.1
17500 1.240291e+07 7.272727e+02 5932.0
17600 1.245052e+07 3.643892e+02 6206.6
17700 1.249547e+07 5.289823e+02 6481.8
17800 1.252851e+07 3.521127e+02 6744.5
17900 1.257572e+07 2.830285e+02 7033.0
18000 1.261236e+07 3.636364e+02 7290.6
18100 1.264811e+07 2.190469e+02 7561.1
18200 1.270000e+07 3.636364e+02 7862.6
18300 1.275129e+07 7.120984e+02 8099.6
18400 1.281308e+07 7.301587e+02 8344.2
18500 1.288145e+07 7.272727e+02 8601.1
18600 1.295418e+07 7.272727e+02 8825.2
```
- **e1div_prof.py** (CPU job 4500422, expansion, 8 procs, EXIT 0): 80 hydro_w/m1 pairs (cycles 015160-018676) ->
  docs/handover/e1div-1014/caltech_out_dtretry/prof_<cycle>.npz (80 files). The 25 slice_*.npz (147 MB total with profs) are kept on Caltech
  in e1div_rst/D1_prof/out, not pushed.
