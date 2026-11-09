# TASK 2026-11-09 viper -> ORION: raycheap2 -- restart fix gates (rstabe-1109) + convergence reconciliation (CPU only)

Answer with `docs/handover/NOTE-2026-11-xx-orion-raycheap2.md` on this branch. Follow-up to your
NOTE-2026-11-09-orion-raycheap (thank you: G3 found a real restart bug, and the p = 0.84 needs reconciling).

## 1. Code and builds (CPU MPI, PROBLEM none; your recipe from the last NOTE is fine)
- FIX  = `fork/rstabe-1109` head `4e9f602c8` = `fork/prod-1112-cand` `f34ada236` + one commit: the hesdirk2
  `time2_auto_be` back-off state (skip-run length, steps left) is written to the restart file (record
  `kM1AutoBeRstMagic`, only while a skip run is active) and restored on restart (log line
  `<rad_m1> rstabe-1109: time2_auto_be skip run restored (len L, left N)`). No key; a file without the record
  restarts in stage mode exactly as before.
- BASE = `fork/prod-1112-cand` `f34ada236`.
- MOD  = `fork/raycheap2-1109` head `9fbd2d3c6` (your NEW of the last TASK + inputs): needed for section 3 only
  (the test-only key rad_m1/implicit_flux_x1min_mod_amp|period). Inputs for everything: `tests_raycheap/` of MOD.
- The predictor: when `pred_ok` is false at the dump (it is after a fallback), the original run also cold-starts
  the next step, so no predictor record is needed; the gate below checks this.

## 2. Restart-fix gates (1 rank unless noted; 40 cycles; `tests_raycheap/cmpbin.py A B`)
Deck H = `tests_raycheap/G_cs_lo_shd_ema16_h2.athinput` (hesdirk2, no mod keys; your runs showed stage fallbacks at
cycles 1, 3, 6, 11, 20 and a 16-step back-off from cycle 21). Deck P = `G_cs_lo_shd_geo2_k8e2.athinput`.

| gate | runs | PASS |
|---|---|---|
| R1 restart (FIX) | H: 40 cycles straight vs `time/nlim=20` + `-r <its rst> time/nlim=40` | m1 4/4, m1_lat 12/12 bitwise; the restarted log has the "rstabe-1109 ... restored" line |
| R2 restart (BASE) | the same pair with BASE | expected to DIFFER (documents the bug; report the first fallback cycle after 20) |
| R3 keys absent | FIX vs BASE, 40 cycles, on H and on P | bitwise 4/4 + 12/12 on both |
| R4 ranks (FIX) | H: 4 ranks vs 1 rank, 40 cycles | bitwise 4/4 + 12/12 |

All runs rc=0, 0 FATAL. Report the fallback cycle lists (NOT ADMISSIBLE / BE counts) of R1's three runs.

## 3. Convergence reconciliation (MOD binary, 24 ranks, the 4 levels c = 0.15, 0.075, 0.0375, 0.01875)
`tests_raycheap/conv.py <c0.15> <c0.075> <c0.0375> <c0.01875>` for each set:
- C1 `Z_Q64_c*.athinput` as they are (the arm that gave p 1.5 / 2.0 on Raven CPU: d 8.7e-5, 3.0e-5, 7.5e-6).
- C2 `Z_Q128T8_c*.athinput` with `rad_m1/implicit_tol=1.0e-13 rad_m1/implicit_lin_tol=1.0e-14` on the command line
  (your p = 0.84 set, tolerance-floor hypothesis).
- C3 `Z_Q64_c*.athinput` with the same tighter tolerances.
Report d, p, t_end, wall, and the stage/BE/fallback counts per run. No pass criterion: we need the numbers to
decide whether the d ~ 3e-6 floor is the solver tolerance and whether Q64 reproduces across machines.
(viper is running, on Raven CPU, a stronger-signal variant: amplitude 0.5 with the change test off, implicit_tol 1e-13.)

## 4. NOTE back
Build shas/compiler, the R1-R4 summary lines, the fallback lists, the C1-C3 conv.py outputs with the per-run
counts, wall times. Do not change code; if R1 fails, report the first differing cycle and the abe state printed.
