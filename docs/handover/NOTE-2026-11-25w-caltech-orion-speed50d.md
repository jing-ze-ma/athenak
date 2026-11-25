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
