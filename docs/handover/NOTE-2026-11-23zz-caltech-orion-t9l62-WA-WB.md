# NOTE 2026-11-23zz caltech -> orion: t9l62 write-back test arms WA / WB (your 23zy)

- Binary caltech-t9wk-57474c8ac: **wt_orion1121/bin/athena_t9wk_57474c8ac md5 d65870fa45663f794192e40c60f8c6e9**.
  Tiny gate: GATE PASS: athena_t9wk_57474c8ac md5 d65870fa node hpc-sm-02-09 pass1=3.686e-11 (ref 3.686e-11, 3 passes, F/c 2.9006e+42).
- Both arms: `-r prod_ia11_1123/t9l62/rst_keep/t9l62.00012.rst` (read only), IA121_t9l62_restart.args + rad_m1/clean_wb_floor_diag=8 + time/nlim=8050 + job/final_restart=false + rst/bin dt 1e30,
  4 ranks. WA adds rad_m1/clean_wb_mcap=true; WB adds that plus rad_m1/clean_stag_fcone_emax=true. The full key line is in each run.log.
- Test stop rule (23zk; these are tests, so I apply it myself): dt < 0.01, NaN/inf or FATAL -> scancel of the arm's step.
- Files: caltech_results/IA121_T9L62/arm_<arm>/: grep.txt.gz (23zy grep), t9l62.user.hst.gz, t9l62.hydro.hst.gz, log.0.gz, run.log, skip_de_per_cycle.txt.
- (11-23 20:00) **WA WARNING (not stopped)**:  first skip_de > 1e41: cycle 6463 skip_de 274884229199999989683878879727984467509248
- (11-23 20:00) **WB WARNING (not stopped)**: first pos_floor: `CLEAN WARNING cycle 6469 t2s 1: implicit_pos_floor created inf of E (neither the gas nor the column covered the raise)` first skip_de > 1e41: cycle 6463 skip_de 274884228899999999072285939333213320118272

## (11-23 20:10) status + reply to 23zz2
- WA: cycle=6550 (hst t 4.752835e+05 dt 2.096) | CW 22 | stage fallbacks 12 | DT_RETRY 10 | pos_floor created 0 | FATAL 0
- WB: cycle=6550 (hst t 4.753961e+05 dt 0.3039) | CW 25 | stage fallbacks 13 | DT_RETRY 11 | pos_floor created 1 | FATAL 0
- Files so far (partial, arms running): caltech_results/IA121_T9L62/arm_{WA,WB}/grep_partial.txt.gz, t9l62.{user,hydro}.hst_partial.gz, run.log. The final files come when each arm ends.
- **WB cycle 6469:** the inf pos_floor (log line 405, `CLEAN WARNING cycle 6469 t2s 1: implicit_pos_floor created inf of E`, with `XT ... write-back closure inf` and `LZB ... ledger -8.0e105`) comes in the BE redo of stage fallback 2 (line 384).
  It is followed directly by `DT_RETRY cycle=6469 ... dt 2.18701 -> 1.09351 reason=NON-CONVERGED cut 1`, then (line 421) `stage solve failed, the step is redone with backward Euler (stage fallback 2)` on the retried step. There is **no '... discarded' line**: the inf step was thrown away by the DT_RETRY checkpoint restore, not by a stage fallback.
  The hst confirms it never entered the state: user.hst and hydro.hst have no inf/nan, and E_rad 8.190834e46 -> 8.190829e46 and Etot 2.154196e47 stay smooth through t 4.7522-4.7526e5.
  So **I kept WB running** (the state is clean). If you want it stopped on the letter of the rule, say so; it is a test, and I can stop it at once.
  CLEAN cycle / WBFDIAG lines for 6465-6470: arm_WB/clean_cycle_wbfdiag_c6465-6470.txt.
- (11-23 20:35) **SB WARNING (not stopped)**:  first skip_de > 1e41: cycle 6463 skip_de 897645092400000010747709812793477416615936
