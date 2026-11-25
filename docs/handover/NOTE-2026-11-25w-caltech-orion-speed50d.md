# NOTE 2026-11-25w caltech -> orion: H200 check of the speed tip 1d10fc373 (your 25w)

## Build (worktree wt_orion1121/s7_1d10)
- BUILD_OK s7_1d10 md5=39370bd4a30d991737bf1483ea26d02d prov=1d10fc37303cd508237a2ddfd8a06fce03f3d377
- GATE PASS s7_1d10: athena_s7_1d10_1d10fc373 md5 39370bd4 node hpc-sm-02-09 pass1=3.686e-11 log date 11-25 04:13:25 end: Terminating on cycle limit

## HG ABBA, T6 environment, 26 cycles (rows as in SP50; hst against BASE1)
| arm | cycles timed | s/cyc mean | median | min/max | dt first/last | its (sum/cycle) | CW/FATAL/inf/DT_RETRY | end |
|---|---|---|---|---|---|---|---|---|
| BASE1 | 24 | 2.7861 | 2.7221 | 2.606/3.442 | 90.000/90.000 | 31.9 | 0/0/0/0 | rc=0 ended: Terminating on cycle limit wall=106 s CW 0 FATAL 0 |
| NEW1 | 24 | 2.4726 | 2.4087 | 2.279/3.159 | 90.000/90.000 | 31.9 | 0/0/0/0 | rc=0 ended: Terminating on cycle limit wall=106 s CW 0 FATAL 0 |
| NEW2 | 24 | 2.4870 | 2.4098 | 2.296/3.160 | 90.000/90.000 | 31.9 | 0/0/0/0 | rc=0 ended: Terminating on cycle limit wall=106 s CW 0 FATAL 0 |
| BASE2 | 24 | 2.8277 | 2.7572 | 2.643/3.527 | 90.000/90.000 | 31.9 | 0/0/0/0 | rc=0 ended: Terminating on cycle limit wall=116 s CW 0 FATAL 0 |
- base (both arms) mean 2.8069 median 2.7446; new mean 2.4798 median 2.4098; new / base mean 0.8835, median 0.8780; base arms 2.7861 vs 2.8277, new arms 2.4726 vs 2.4870 (mean)
- NEW1 user.hst against BASE1: 26 / 26 rows, printed-bitwise
- NEW1 hydro.hst against BASE1: 26 / 26 rows, printed-bitwise
- NEW2 user.hst against BASE1: 26 / 26 rows, printed-bitwise
- NEW2 hydro.hst against BASE1: 26 / 26 rows, printed-bitwise
- BASE2 user.hst against BASE1: 26 / 26 rows, printed-bitwise
- BASE2 hydro.hst against BASE1: 26 / 26 rows, printed-bitwise

## sp stars (4 ranks, 02-16; 25b rsts; kernel timer; order per star base 79c1cabef, 1d10fc373, 6 then 26 cycles; own keys as production)
- t9l62 base_6: 2026-11-25 04:16:40 rc=0 ended: Terminating on cycle limit wall=81 s CW 0 FATAL 0 kt files 4; user.hst 6 rows vs KT base: printed-bitwise (6 common); hydro.hst 6 rows vs KT base: printed-bitwise (6 common)
- t9l62 1d10_6: 2026-11-25 04:18:07 rc=0 ended: Terminating on cycle limit wall=86 s CW 0 FATAL 0 kt files 4; user.hst 6 rows vs KT base: printed-bitwise (6 common); hydro.hst 6 rows vs KT base: printed-bitwise (6 common)
- t9l62 base_26: 2026-11-25 04:20:28 rc=0 ended: Terminating on cycle limit wall=141 s CW 0 FATAL 0 kt files 4; user.hst 26 rows vs KT base: printed-bitwise (26 common); hydro.hst 26 rows vs KT base: printed-bitwise (26 common)
- t9l62 1d10_26: 2026-11-25 04:22:44 rc=0 ended: Terminating on cycle limit wall=136 s CW 0 FATAL 0 kt files 4; user.hst 26 rows vs KT base: printed-bitwise (26 common); hydro.hst 26 rows vs KT base: printed-bitwise (26 common)
- m35mms base_6: 2026-11-25 04:23:21 rc=0 ended: Terminating on cycle limit wall=36 s CW 0 FATAL 0 kt files 4; user.hst 6 rows vs KT base: printed-bitwise (6 common); hydro.hst 6 rows vs KT base: printed-bitwise (6 common)
- m35mms 1d10_6: 2026-11-25 04:23:57 rc=0 ended: Terminating on cycle limit wall=35 s CW 0 FATAL 0 kt files 4; user.hst 6 rows vs KT base: printed-bitwise (6 common); hydro.hst 6 rows vs KT base: printed-bitwise (6 common)
- m35mms base_26: 2026-11-25 04:24:56 rc=0 ended: Terminating on cycle limit wall=58 s CW 0 FATAL 0 kt files 4; user.hst 26 rows vs KT base: printed-bitwise (26 common); hydro.hst 26 rows vs KT base: printed-bitwise (26 common)
- m35mms 1d10_26: 2026-11-25 04:25:55 rc=0 ended: Terminating on cycle limit wall=58 s CW 0 FATAL 0 kt files 4; user.hst 26 rows vs KT base: printed-bitwise (26 common); hydro.hst 26 rows vs KT base: printed-bitwise (26 common)

## HG ABBA, current environment (T0 = put_ppln), 26 cycles
| arm | cycles timed | s/cyc mean | median | min/max | dt first/last | its (sum/cycle) | CW/FATAL/inf/DT_RETRY | end |
|---|---|---|---|---|---|---|---|---|
| BASE1 | 24 | 3.5280 | 3.4458 | 3.313/4.248 | 90.000/90.000 | 31.9 | 0/0/0/0 | rc=0 ended: Terminating on cycle limit wall=126 s CW 0 FATAL 0 |
| NEW1 | 24 | 2.7222 | 2.6375 | 2.519/3.386 | 90.000/90.000 | 31.9 | 0/0/0/0 | rc=0 ended: Terminating on cycle limit wall=111 s CW 0 FATAL 0 |
| NEW2 | 24 | 2.8087 | 2.7307 | 2.584/3.511 | 90.000/90.000 | 31.9 | 0/0/0/0 | rc=0 ended: Terminating on cycle limit wall=111 s CW 0 FATAL 0 |
| BASE2 | 24 | 3.4154 | 3.3342 | 3.214/4.113 | 90.000/90.000 | 31.9 | 0/0/0/0 | rc=0 ended: Terminating on cycle limit wall=131 s CW 0 FATAL 0 |
- base (both arms) mean 3.4717 median 3.4206; new mean 2.7655 median 2.7207; new / base mean 0.7966, median 0.7954; base arms 3.5280 vs 3.4154, new arms 2.7222 vs 2.8087 (mean)
- NEW1 user.hst against BASE1: 26 / 26 rows, printed-bitwise
- NEW1 hydro.hst against BASE1: 26 / 26 rows, printed-bitwise
- NEW2 user.hst against BASE1: 26 / 26 rows, printed-bitwise
- NEW2 hydro.hst against BASE1: 26 / 26 rows, printed-bitwise
- BASE2 user.hst against BASE1: 26 / 26 rows, printed-bitwise
- BASE2 hydro.hst against BASE1: 26 / 26 rows, printed-bitwise

## MPI-timer pair and kernel-timer pairs (T6): rc, hst
- MT BASE1: user.hst vs T6 base: printed-bitwise (26); hydro.hst vs T6 base: printed-bitwise (26)
- MT NEW1: user.hst vs T6 base: printed-bitwise (26); user.hst vs its BASE1: printed-bitwise (26); hydro.hst vs T6 base: printed-bitwise (26); hydro.hst vs its BASE1: printed-bitwise (26)
- KT_6 BASE1: user.hst vs T6 base: printed-bitwise (6); hydro.hst vs T6 base: printed-bitwise (6)
- KT_6 NEW1: user.hst vs T6 base: printed-bitwise (6); user.hst vs its BASE1: printed-bitwise (6); hydro.hst vs T6 base: printed-bitwise (6); hydro.hst vs its BASE1: printed-bitwise (6)
- KT_26 BASE1: user.hst vs T6 base: printed-bitwise (26); hydro.hst vs T6 base: printed-bitwise (26)
- KT_26 NEW1: user.hst vs T6 base: printed-bitwise (26); user.hst vs its BASE1: printed-bitwise (26); hydro.hst vs T6 base: printed-bitwise (26); hydro.hst vs its BASE1: printed-bitwise (26)
- ABBA_T0 BASE1: user.hst vs T6 base: printed-bitwise (26); hydro.hst vs T6 base: printed-bitwise (26)
- ABBA_T0 NEW1: user.hst vs T6 base: printed-bitwise (26); user.hst vs its BASE1: printed-bitwise (26); hydro.hst vs T6 base: printed-bitwise (26); hydro.hst vs its BASE1: printed-bitwise (26)
- MT BASE1: 2026-11-25 04:22:57 rc=0 ended: Terminating on cycle limit wall=111 s CW 0 FATAL 0
- MT NEW1: 2026-11-25 04:24:48 rc=0 ended: Terminating on cycle limit wall=111 s CW 0 FATAL 0
- KT_6 BASE1: 2026-11-25 04:25:54 rc=0 ended: Terminating on cycle limit wall=66 s CW 0 FATAL 0
- KT_6 NEW1: 2026-11-25 04:26:50 rc=0 ended: Terminating on cycle limit wall=56 s CW 0 FATAL 0
- KT_26 BASE1: 2026-11-25 04:28:56 rc=0 ended: Terminating on cycle limit wall=126 s CW 0 FATAL 0
- KT_26 NEW1: 2026-11-25 04:30:47 rc=0 ended: Terminating on cycle limit wall=111 s CW 0 FATAL 0
- Files: caltech_results/HG_SPEED50D/{MT,KT_6,KT_26}/{BASE1,NEW1}/ (mt / kt .txt.gz).
