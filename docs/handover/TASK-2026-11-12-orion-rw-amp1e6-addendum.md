# TASK 2026-11-12 addendum: 6 more grid cases move to Orion (viper standby taken by R136 production)

Same TASK, builds, runner and keys as TASK-2026-11-12-orion-rw-amp1e6.md. At 00:30 CET 11-12 the viper standby became the
R136 production holder; viper's CPU steps there were stopped. Everything of the amp-1e-6 grid is done on viper EXCEPT
the nine thin tau 0.01 cases below (six N256, three N128 killed at ~80 %); they now run ONLY on Orion (viper does not run them; no duplicates):

| file (docs/handover/scripts/rw_amp1e6/) | cases | binary |
|---|---|---|
| orion_grid256_PROD.txt | P 0.0333, 0.333, 33.3 at tau 0.01, N 256, nt 512, tol 1e-12 / lin 1e-13 / cnorm 1e-13 | 7d233ac30 |
| orion_grid256_FINAL.txt | same three, FINAL keys | 84cd20888 |
| orion_grid128_FINAL.txt | P 0.0333, 0.333 at tau 0.01, N 128, nt 256, grid tol, FINAL keys | 84cd20888 |
| orion_tt128_FINAL.txt | P 0.0333 tau 0.01 N 128, tol 1e-13 / lin 1e-14 / cnorm 1e-14, FINAL keys (RW_RUNS=runs_ttFINAL) | 84cd20888 |

Run them in the SAME run dirs as the grid arm (RW_RUNS=$PWD/runs_PROD / runs_FINAL), 6 cores in parallel (~1.5 h per case
on one viper core); P 33.3 is info only and FATALs at amp 1e-5 and 1e-6 on viper at lower N (a FATAL is a result, keep
dt_retry off). Together with the 8 N256 tolerance-check cases of the TASK: 17 cases, 17 cores (N128 ~25 min each).

```
RW_EXE=<wt_rwa6_prod>/build_none/src/athena RW_RUNS=$PWD/runs_PROD  python3 drive.py run orion_grid256_PROD.txt  --np 3 > runs_PROD.out 2>&1 &
RW_EXE=<wt_rwa6_fin>/build_none/src/athena  RW_RUNS=$PWD/runs_FINAL python3 drive.py run orion_grid256_FINAL.txt --np 3 > runs_FINAL.out 2>&1 &
```

```
RW_EXE=<wt_rwa6_fin>/build_none/src/athena RW_RUNS=$PWD/runs_FINAL   python3 drive.py run orion_grid128_FINAL.txt --np 2 > runs_FINAL128.out 2>&1 &
RW_EXE=<wt_rwa6_fin>/build_none/src/athena RW_RUNS=$PWD/runs_ttFINAL python3 drive.py run orion_tt128_FINAL.txt  --np 1 > runs_ttFINAL128.out 2>&1 &
```
Before running, `cp -r n128_ref/runs_* .` (viper's finished results.json of the tau 0.01 cells, N 16..128, and the
viper tolerance-check cases); drive.py skips cases with a results.json, and analyze.py then prints the full 128 -> 256
orders for these cells.

Pass: P 0.0333 / 0.333 cases end "ok". Return in the NOTE, per case: L1 rho, L1 E raw, L1 E mean-removed, mean E offset
(analyze.py rows PROD / FINAL, N 256 column), and attach each case's results.json (small); viper forms the 128 -> 256 orders
from the analyze.py table.
