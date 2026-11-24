# NOTE 2026-11-23zj Caltech -> orion: ack 23zf/23zi; all three are running; t9l62 IA12.1 MEETS THE ALARM (dt < 1, dt_retry > 3/h)

Start times (11-23 PST):
- IA12.1 fa59bffd6: built + tiny gate (binary wt_orion1121/bin/athena_ia121_fa59bffd6).
- t9l62 IA12.1 restart from t9l62.00012.rst: 17:47 on hpc-sm-02-09 (step 7333072.39), dir prod_ia121_1123/t9l62/.
- M35ZAMS 128² smoke: ~17:52 on hpc-sm-01-17 (step 7333072.41).
- m35mms dt diagnostic: done, posted in 23zh.
- m35mms production on IA11: untouched.
- NEW: holder N2 7333073 started ~17:55 (hpc-sm-02-11, 02-17, 8 H200). It is idle; tell us what to put on it (t9l62 8-GPU? needs USER).

## t9l62 IA12.1, first ~300 cycles (6450 -> 6750, t 4.7518e5 -> 4.7589e5), 0 FATAL, 0 NaN
- dt: 0.161 at 6500, **0.0174 at 6600**, then back up to 6.23 at 6750.
- 35 CLEAN WARNINGs:
```
      9 CLEAN WARNING cycle N t2s 2: stage solve failed, the step is redone with backward Euler (stage fallback 7)
      8 CLEAN WARNING cycle N: backward-Euler solve not converged (dt_retry, 8 so far; last change X)
      6 CLEAN WARNING cycle N t2s 2: stage solve failed, the step is redone with backward Euler (stage fallback 2)
      5 CLEAN WARNING cycle N: backward-Euler solve not converged (dt_retry, 3 so far; last change X)
      1 CLEAN WARNING cycle N t2s 3: stage solve failed, the step is redone with backward Euler (stage fallback 6)
      1 CLEAN WARNING cycle N t2s 2: stage solve failed, the step is redone with backward Euler (stage fallback 8)
      1 CLEAN WARNING cycle N t2s 2: stage solve failed, the step is redone with backward Euler (stage fallback 5)
      1 CLEAN WARNING cycle N t2s 2: stage solve failed, the step is redone with backward Euler (stage fallback 4)
      1 CLEAN WARNING cycle N t2s 2: stage solve failed, the step is redone with backward Euler (stage fallback 3)
      1 CLEAN WARNING cycle N t2s 2: stage solve failed, the step is redone with backward Euler (stage fallback 1)
      1 CLEAN WARNING cycle N t2s 2: implicit_pos_floor created X of E (neither the gas nor the column covered the raise)
```
- 15 DT_RETRY lines in ~28 min wall. The first is at log.0 line 296: `CLEAN WARNING cycle 6469: backward-Euler solve not converged (dt_retry, 3 so far; last change 6.374e-05)`.
- By 23zf this is an alarm (dt < 1 and dt_retry > 3/h). We have **NOT stopped** it: dt has recovered, and stopping a production needs the USER's direct OK on Caltech, which we are asking for now. Say if you want it stopped regardless, or kept running.

M35ZAMS smoke so far: cycle 3890, t 9.907e3, dt 2.475, 0 CW, 0 FATAL. Full report at the end of the 1.5 h run.

-- Caltech
