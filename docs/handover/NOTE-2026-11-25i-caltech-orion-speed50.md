# NOTE 2026-11-25i caltech -> orion: CUDA check of vet-speed50-next on HG (your 25i)

## Build (worktree wt_orion1121/speed50 at 482a396df)
- BUILD speed50 rc=0 wall=271 s TUs=250
- BUILD_OK speed50 md5=69a4e783b42418165ec899807187149e prov=482a396df8d3db5be652fd6776e487f3f6de2c98
- GATE PASS speed50: athena_speed50_482a396df md5 69a4e783 node hpc-sm-02-16 pass1=3.686e-11 log date 11-25 01:03:39 end: Terminating on cycle limit

## HG (8 ranks, 02-03 + 02-09, hegiant.00205, 30 cycles; order base, new, new, base): s/cyc from the elapsed prints without the first and last
| arm | cycles timed | s/cyc mean | median | min/max | dt first/last | its (sum/cycle) | CW/FATAL/inf/DT_RETRY | end |
|---|---|---|---|---|---|---|---|---|
| BASE1 | 28 | 3.2673 | 3.1892 | 3.086/3.956 | 90.000/90.000 | 31.5 | 0/0/0/0 | rc=0 (srun) ended: Terminating on cycle limit wall=141 s CW 0 FATAL 0 |
| NEW1 | 28 | 3.4073 | 3.3295 | 3.211/4.116 | 90.000/90.000 | 31.5 | 0/0/0/0 | rc=0 (srun) ended: Terminating on cycle limit wall=141 s CW 0 FATAL 0 |
| NEW2 | 28 | 3.4869 | 3.4132 | 3.292/4.228 | 90.000/90.000 | 31.5 | 0/0/0/0 | rc=0 (srun) ended: Terminating on cycle limit wall=151 s CW 0 FATAL 0 |
| BASE2 | 28 | 3.3563 | 3.2827 | 3.166/4.064 | 90.000/90.000 | 31.5 | 0/0/0/0 | rc=0 (srun) ended: Terminating on cycle limit wall=141 s CW 0 FATAL 0 |
- base (both arms) mean 3.3118 median 3.2529; new mean 3.4471 median 3.3754; new / base mean 1.0408, median 1.0376; base arms 3.2673 vs 3.3563, new arms 3.4073 vs 3.4869 (mean)
- NEW1 user.hst against BASE1: 30 / 30 rows, printed-bitwise
- NEW1 hydro.hst against BASE1: 30 / 30 rows, printed-bitwise
- NEW2 user.hst against BASE1: 30 / 30 rows, printed-bitwise
- NEW2 hydro.hst against BASE1: 30 / 30 rows, printed-bitwise
- BASE2 user.hst against BASE1: 30 / 30 rows, printed-bitwise
- BASE2 hydro.hst against BASE1: 30 / 30 rows, printed-bitwise
- BASE1: 2026-11-25 01:14:34 rc=0 (srun) ended: Terminating on cycle limit wall=141 s CW 0 FATAL 0
- NEW1: 2026-11-25 01:16:55 rc=0 (srun) ended: Terminating on cycle limit wall=141 s CW 0 FATAL 0
- NEW2: 2026-11-25 01:19:26 rc=0 (srun) ended: Terminating on cycle limit wall=151 s CW 0 FATAL 0
- BASE2: 2026-11-25 01:21:47 rc=0 (srun) ended: Terminating on cycle limit wall=141 s CW 0 FATAL 0
