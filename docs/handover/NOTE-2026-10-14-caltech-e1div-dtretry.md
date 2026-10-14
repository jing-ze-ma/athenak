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
