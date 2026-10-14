# NOTE 2026-10-14 Caltech E1S (TASK-2026-10-14-caltech-e1s + addendum b7e68b9f)

Caltech test-runner. Dir `/resnick/groups/carnegie_poc/jingze/e1s_1014/{scripts,probe,gate,prod}`.

## 1. Build / deck md5 (all match the TASK)
- BIN `/resnick/home/jingze/ATHENAK/builds/bin/athena_gpu_he_star_m1_5b7e1bc5_nofma` md5 6252e58e5192ce65ff1f1f627e1ce90c (5b7e1bc5, fin4-dtretry)
- deck template `agcA_E1_fin4dtr_s3.athinput.in` 9b1d3d2d; `e1s_1014/agcA_E1_fin4dtr_s3.athinput` (@FILES@ = e1_1013/files, a-w) **d32cbb98581a19f6dd1b8760bd608300**
- files: IC e81deb51, rosseland b28e97c5, planck 493f9d01
- scripts = e1caltech-1013 set (md5 as MD5SUMS) + hook e1s-1014/e1chain_hook.sh 8114774b + pulse_ana.py 34d74804.
  One local change to `link_e1.sh` (97a2239a -> 039337aa) for the addendum: after the MARGIN line,
  `MARGIN=${E1C_MARGIN:-$MARGIN}` (default unchanged; production run.cfg sets E1C_MARGIN=1200 -> athena -t = walltime - 20 min).
- LAUNCH one-word scripts: `e1s_1014/launch2.sh`, `launch4.sh` = `exec srun --mpi=pmix -n N -c 8 --cpu-bind=cores e1_1013/e1wrap.sh "$@"`.

## 2. Scaling probe (fresh, nlim 100, ndiag 1, bins/log/rst off; job 4513176, debug QoS, hpc-sm-02-10)
| arm | layout | rc | FATAL/NaN/DIVERGED/NC/NOT ADM/DT_RETRY | dt (s) | Picard mean/max | init (s) | s/cycle 10-100 | peak GB/GPU | ipc | dt_retry |
|---|---|---|---|---|---|---|---|---|---|---|
| G2 | 1x2 H200, 536x64x64, 4 MB (2/GPU) | 0 | 0/0/0/0/0/0 | 727.27 | 3.056 / 6 | 39.1 | **0.538** | 58.5 | ON | on (fac 0.5 max 8 ngood 10 grow 1.25) |
| G4 | 1x4 H200, 536x64x64, 4 MB (1/GPU) | 0 | 0/0/0/0/0/0 | 727.27 | 3.056 / 6 | 35.4 | **0.313** | 30.4 | ON | on (same) |
| G8 | 2x4, 536x64x32 | skipped: queue estimate > 6 h (normal 20 h, debug 7 h) | | | | | | | | |
| G16 | 4x4, 536x32x32 | skipped: queue estimate > 6 h (normal 28 h, debug 20 h) | | | | | | | | |
- G2 vs ref 0.570 s/cycle, 59.9 GB: consistent. G4 speed-up 1.72x on 2x the GPUs.
- hst last row (t 7.05e4 s, 8 rows) G4 vs G2 max rel diff: hydro.hst 4.2e-7 (col 4); user.hst 2.8e-2 in col 25 P_esrc
  (-5.5e21 vs L 2.8e38: a small difference of large terms, 2e-17 of L), next col 13 4.3e-7, col 21 v1sq_wall 2.2e-7, L_top (col 3) 3.9e-12.
  Not bitwise (rank count differs; same MeshBlocks), round-off grown over 100 cycles.

## 3. Cluster limits (raw in e1s_1014/probe/limits.txt) and queue estimates (sbatch --test-only, 05:47-05:53 PDT)
- partition gpu MaxTime 14-00:00:00; 18 H200 nodes x 4 (01-03 drained, 01-04/02-05 reserved, 02-10 planned; exclude 01-09).
- QoS: normal MaxWall 7 d; long MaxWall 14 d, MaxJobsPU 1; debug 30 min, MaxTRES node=10, MaxJobsPU 1, Priority 10000.
  Assoc carnegie_poc: QOS debug,long,normal; no per-assoc limits. FairShare 0.309.
- Start estimates (identical for -t 24 h / 48 h / 7 d normal / 14 d long: no backfill penalty seen for long -t):
  1x2 H200 now; 1x4 H200 13:03 (+7.2 h); 2x4 2026-10-15 01:54 (+20 h); 4x4 10-15 10:10 (+28 h); debug 2x4 12:52, 4x4 10-15 01:43.

## 4. Chosen layout
Run time = 33,300 cycles x s/cycle x 1.3:
- G2: 33,300 x 0.538 x 1.3 = 23.3 ks = 6.5 h; queue 0 -> **~6.5 h**
- G4: 33,300 x 0.313 x 1.3 = 13.5 ks = 3.8 h; queue 7.2 h -> ~11 h
- G8/G16: queue 20-28 h alone exceeds G2's whole run.
**Choice: G2 = 1 node x 2 H200, deck meshblock 536x64x64 (no second deck), LAUNCH launch2.sh, QoS long, -t 14-00:00:00**
(longest allowed; start estimate identical to 7 d normal; one link covers the run, 2 afterany backups at the same walltime).
