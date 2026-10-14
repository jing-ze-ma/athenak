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

### Production deck (coordinator relay of USER 10-14: raise tlim)
**tlim raised to 1e9 per user 10-14; the run is limited by walltime links.** Production deck
`e1s_1014/agcA_E1_fin4dtr_s3_tlim1e9.athinput` (a-w) md5 **680dadb0a706715b5c9187211e0b431c**; diff vs d32cbb98 (only line):
```
50c50
< tlim                     = 2.4192e+07   # agcmodel-1012
---
> tlim                     = 1.0e+09      # e1s-1014: USER 10-14 tlim raised to 1e9; the run is limited by walltime links
```
Meshblock unchanged (G2 = deck 536x64x64). run.cfg TLIM=1e9. The gate below ran on d32cbb98 (submitted before the tlim relay);
tlim only sets the stop time, so the restart/IPC bitwise result carries over.

## 5. Gate (job 4513204, 1x2 H200 hpc-sm-01-05, launch2.sh, deck d32cbb98)
Smoke (`e1s_1014/gate/smoke`): C fresh 60 / L 30+30 hook restart / Coff ipc=false; rc 0, FATAL 0 in all 4 runs; Picard C 3.09/6.
```
hst agcar3d.hydro.hst: C rows 11, chain rows 11, missing 0 [], extra 0 [], differing 0
hst agcar3d.user.hst: C rows 11, chain rows 11, missing 0 [], extra 0 [], differing 0
bin: C dumps 8, chain dumps 8, extra in chain 0 []
CHAINCMP BITWISE
CMPLAST BITWISE
HST agcar3d.hydro.hst BITWISE
HST agcar3d.user.hst BITWISE
IPC-OFF (Coff vs C): CMPLAST BITWISE, HST hydro BITWISE, HST user BITWISE
GATE PASS 2026-10-14_05:59:27   (-> gate/smoke/VERIFIED)
```

### Background check BG (same job; fresh via link_e1.sh, E1C_NCYC=1400, production cadence; t = 1.018e6 s)
rc 0; FATAL 0, NaN 0, DIVERGED 0, DT_RETRY 0, NON-CONVERGED 0, NOT ADMISSIBLE 0; fofc_report at-floor=0 other=0 in all 282 reports;
Picard mean 3.00 max 6; dt 727.27; s/cycle (10-1400) 0.504. `gate/BG_pulse.txt` (pulse_ana.py run):
| t (s) | vrV | vrM | vr105 | vr110 | d105 | d110 | KEe | rsurf |
|---|---|---|---|---|---|---|---|---|
| 1.0e5 | -0.042 | 0.015 | -0.060 | -0.055 | 0.9992 | 0.9982 | 1.1e36 | 1.009 |
| 5.0e5 | -0.026 | 0.109 | -0.060 | -0.055 | 1.009 | 0.990 | 9.1e36 | 1.009 |
| 9.0e5 | -0.066 | -0.045 | -0.060 | -0.055 | 1.010 | 0.988 | 1.4e37 | 1.009 |
| 1.0e6 | -0.103 | -0.043 | -0.060 | -0.055 | 1.010 | 0.988 | 2.4e37 | 1.009 |
PU S3_268 ref at 5e5 / 9e5: vrV 0.076 / -0.039, vrM 0.58 / -0.19, rsurf 1.009 / 1.011 -> same magnitudes (|v| <~ 0.1 km/s), no infall.
- shell-mean <v_r> above 1.1 R_ph: max |v_r| 0.093 km/s (same value from the first dump on, not growing); sponge (r > 7.29e13, 6 cells)
  max |v_r| 0.069 km/s; <rho>/<rho>(0) above 1.1 R_ph within 0.9949-1.003.
- user.hst (102 rows): L_top/L0 0.99987-1.00007 (last 1.00007); Mdot_top 0 -> 9.6e14 (>= 0, outward); Min_top 0 throughout;
  Picard 3.0-3.6; w_mlt 1.0; P_sponge 0 -> 9.9e29.
**BG PASS** (quiet, like the S3 reference). `prod/GO_e1` written.

## 6. STARTED (production, `e1s_1014/prod`)
- **tlim raised to 1e9 per user 10-14; the run is limited by walltime links.** Deck 680dadb0 (tlim 1e9), run.cfg TLIM=1e9, KEYS empty,
  STOP_NC=1, E1C_GATE=1, E1C_VERIFIED=gate/smoke/VERIFIED, LAUNCH=e1s_1014/launch2.sh, E1C_DIR=e1s_1014/scripts, E1C_MARGIN=1200.
- Links `agcE1s` (link_e1.sh + e1s-1014 hook), 1 node x 2 H200, `-A carnegie_poc -p gpu --qos=long --exclude=hpc-sm-01-09
  --ntasks=2 -c 8 --gres=gpu:nvidia_h200:2 --mem=240G -t 14-00:00:00` (longest allowed: partition MaxTime 14 d = long QoS MaxWall;
  normal is 7 d): **4513330** (link 1) -> 4513331 -> 4513332 (afterany, same walltime).
- Link 1 RUNNING since 2026-10-14 06:12:20 PDT on hpc-sm-01-05; fresh start (prod/rst was empty), athena `-t 335:39:57`
  (= 14 d - 20 min). md5 checked in job.4513330.out (bin 6252e58e, deck 680dadb0, IC/tables as section 1).
- First s/cycle: **0.464** (cycles 50-1200, t 8.7e5 s at 06:22); dt 727.27; FATAL 0, DT_RETRY 0, fofc at-floor 0 so far.
  At this rate 2.4192e7 s (33,300 cycles) is reached ~4.3-5.6 h after start; the run then continues toward tlim 1e9 until the walltime.
- Evolution reports (TASK section 6) are not part of this runner's job.
