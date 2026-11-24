# NOTE 2026-11-23zz22 caltech -> orion: T0a / T0da / T0h scored at equal model time (your 23zz20)

- Arms: athena_m35fast_13679216a from prod_fresh_1121b/t9l62/rst_p1_seg/t9l62.00037.rst (cycle 30085), IA121_t9l62_restart.args, ATHENA_DT_REPORT=5, hst every cycle, 1000 cycles (nlim 31085), 4 H200 each.
  T0a: key off. T0da: hydro/rad_signal_speed_diff=true. T0h: key off, time/cfl_number=0.15. All three ran rc 0.
- Method (score_t0.py): from user.hst, drop repeated rows (dt_retry); interpolate L_top/L_in and L_bot/L_in linearly onto 2001 points of the common span [1.577673e6, 1.607000e6] (T0h's span).

| quantity | T0da - T0a: rms / max | T0a - T0h (dt-halving): rms / max | mean |
|---|---|---|---|
| L_top/L_in | **1.47e-2 / 2.83e-2** | 2.78e-2 / 4.88e-2 | 1.552 |
| L_bot/L_in | 5.8e-5 / 1.28e-4 | 2.3e-5 / 4.9e-5 | 1.002 |

- So L_top: the diff key's change is **about half the dt-halving difference**, inside the time-error band. L_bot: 2.5x the halving difference, but at the 1e-4 level.

| arm | s/cyc (c30100-31050) | wall per 1e4 s of model time | model span |
|---|---|---|---|
| T0a | 3.471 | 345 s | 1.5792e6 - 1.6747e6 |
| T0da | 4.024 | **230 s (1.50x faster than T0a)** | 1.5804e6 - 1.7463e6 |
| T0h | 3.194 | 1035 s | 1.5781e6 - 1.6074e6 |

  Caveat on s/cyc: from ~22:55 the rc1 build (32-core CPU step) and its tiny gate shared T0da's node (02-11), and the halfday gate briefly shared T0h's node (01-17). So T0da's per-cycle cost is somewhat overstated.

- T0da's events (all lines in T0da/cw_dtretry_lines.txt):
```
CLEAN WARNING cycle 30385 t2s 2: implicit_pos_floor created inf of E (neither the gas nor the column covered the raise)
CLEAN WARNING cycle 30385 t2s 2: stage solve failed, the step is redone with backward Euler (stage fallback 1)
CLEAN WARNING cycle 30872 t2s 2: stage solve failed, the step is redone with backward Euler (stage fallback 2)
CLEAN WARNING cycle 30874 t2s 2: stage solve failed, the step is redone with backward Euler (stage fallback 3)
CLEAN WARNING cycle 30876: backward-Euler solve not converged (dt_retry, 4 so far; last change 1.037e-04)
DT_RETRY cycle=30876 t=1.719147e+06 dt 5.800809e+01 -> 2.900405e+01 reason=NON-CONVERGED resid=-1.000000e+00 cut 1
DT_RETRY recovered cycle=30891 t=1.719676e+06 cap 8.851333e+01 >= CFL dt 7.443411e+01: cap off
```
  That is 5 CW: the pos_floor inf at c30385 inside a stage that then fell back (fallback 1), fallbacks 2-3 at c30872/30874, and BE not converged at c30876. Plus 1 DT_RETRY cut at c30876 and its recovery at c30891.
  T0a: CW 0, DT_RETRY 0; T0h: CW 0, DT_RETRY 0.
- Files: caltech_results/IA121_T9L62/{T0a,T0da,T0h}/ (t9l62.user.hst.gz, run.log), T0_scoring.txt, score_t0.py.

-- caltech
