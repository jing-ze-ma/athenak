# TASK raven -> Caltech: E1 FRESH with seed 1e-3 on FIN4 + dt_retry (5b7e1bc5); GPU scaling probe; long links

From the Raven coordinator (worker E1S, e1s-1014). Files: `docs/handover/e1s-1014/`. Report NOTE: `docs/handover/NOTE-2026-10-14-caltech-e1s.md`.

USER (10-14, verbatim): "can we start e1 fresh with a smaller perturbation with the optimal scheme on caltech? tell it to probe how many
nodes can it use and potentially run long jobs"

Why: PU (pulse-1014, Raven reduced 3-D E1) found that the big E1 excursion/infall (your E1 diverged at 1.24e7 s; D1 needed 116 dt
retries to 1.30e7 s) is driven by the 1e-2 seed on the full-strength MLT scaffold. With `problem/he_seed = 1.0e-3` there is no infall
(v_r > 0 in 1.0-1.2 R_ph throughout, edge 1.077 R_ph, 0 dt retries), and convection still develops (onset ~7e6 s, lateral KE 1e37 ->
3e41 by 9e6 s). This TASK is that run, full size, from t = 0.

## Rules
- Do NOT touch running/queued production: BSG (bsgLIw4 4490079 and any BSG link) and He giant (hgF4w 4489939 and links). If H200s are
  short, queue behind them. Do not cancel anything. Do not reuse or modify `e1_1013/` or `e1div_*` dirs (read-only refs).
- New dir: `/resnick/groups/carnegie_poc/jingze/e1s_1014/{scripts,probe,gate,prod}`.
- fin4-mem is NOT adopted: do not use any fin4-mem binary.
- Never put an ABSENT key on the command line (AthenaK FATALs: PU's first submit died with `time/dt_retry=true` on a fresh command line).
  Every key you override below already exists in the deck.

## 1. Binary + deck (record both md5 in the NOTE)
- Binary: your `builds/bin/athena_gpu_he_star_m1_5b7e1bc5_nofma`, md5 must be **6252e58e5192ce65ff1f1f627e1ce90c** (fin4-dtretry 5b7e1bc5).
- Files dir `<F>` = your `e1_1013/files` (IC e81deb51, rosseland b28e97c5, planck 493f9d01; check md5).
- Deck: `e1s-1014/agcA_E1_fin4dtr_s3.athinput.in` (md5 9b1d3d2d):
  ```
  sed "s#@FILES@#<F>#g" agcA_E1_fin4dtr_s3.athinput.in > e1s_1014/agcA_E1_fin4dtr_s3.athinput; chmod a-w ...
  ```
  With `<F>` = `/resnick/groups/carnegie_poc/jingze/e1_1013/files` the md5 is **d32cbb98581a19f6dd1b8760bd608300** (that path
  reproduces your old deck d29e4615 from its template, so it should match).
- What it is: the Viper fresh FIN4 + dt_retry E1 deck `agcA_E1_fin4dtr.athinput` (md5 21c55a7d, copy in `e1s-1014/` for reference;
  DR drift-1013 "VIPER E1 BACKUP with dt_retry") with exactly these differences (`diff` of the two files):
  | line | Viper 21c55a7d | this deck |
  |---|---|---|
  | 109 (`<rad_m1>`) | (absent; default false) | `vet_gd_halo_ipc = true` (as your E1/BSG; node-local CUDA IPC, bitwise; MPI fallback off-node) |
  | 351-353 (`<problem>`) | `/viper/ptmp2/.../agc/files/...` | `@FILES@/...` (IC, rosseland, planck: same files) |
  | 358 (`<problem>`) | `he_seed = 1.0e-2` | **`he_seed = 1.0e-3`** |
  Everything else is byte-identical, including the meshblock line (536 x 64 x 64; layout section below).
- vs your old E1 deck d29e4615 (FIN3, 62115d17 template), the only other differences are the FIN4 + dt_retry keys, all NAMED in the
  file: `<time> dt_retry = true`, `<rad_m1> time2_enth_be = true`, `implicit_cs2_vimp = true`, `implicit_cs2_vimp_npass = 3`
  (plus one comment). On a fresh start these are also the 5b7e1bc5 defaults; naming them keeps them on after every restart.
- Unchanged and checked: wb_rmax 2.9715e13; DI IC (he_ic_equil = true etc.); MLT scaffold ramp start 4.032e6 s, time 2.016e7 s;
  hst 1e4 s, hydro_w + m1 bins 1e5 s, log 1e5 s, rst 5e5 s; tlim 2.4192e7; mesh 536 x 128 x 128.

## 2. Scripts (one dir `e1s_1014/scripts`)
- Copy all of `e1caltech-1013/` (link_e1.sh, rst_pick.py, prune_rst.sh, smoke_e1.sh, smoke_gate.sh, chaincmp.py, cmplast.py,
  bin_convert.py) and REPLACE the hook with `e1s-1014/e1chain_hook.sh` (MD5SUMS in e1s-1014/).
- Hook change (needed with 5b7e1bc5): a RECOVERED retry prints `DT_RETRY ... reason=DIVERGED ...` (D1 had 10 such lines). The old hook
  STOPs on any `DIVERGED` in the log, i.e. it would halt the chain after a link that was fine. The new hook ignores lines containing
  `DT_RETRY` for the DIVERGED rule only. Unchanged: rc != 0, any `FATAL` (an exhausted retry is `### FATAL ERROR: <time> dt_retry:
  cycle ... failed`), `nan` in the last 60 lines, end-summary `NON-CONVERGED=` > 0 (accepted steps only) -> STOP.
  (Viper avoided this with 053bd996, which only renames the reason words; numerics = 5b7e1bc5. You stay on 5b7e1bc5.)
- LAUNCH: use a ONE-WORD launch script (the e1_1013/e1launch.sh lesson: smoke_e1.sh loses quotes of a multi-word LAUNCH), e.g.
  `e1s_1014/launchN.sh` = `exec srun --mpi=pmix -n $N -c 8 --cpu-bind=cores e1wrap.sh "$@"` (e1wrap.sh = module env +
  `CUDA_VISIBLE_DEVICES=$SLURM_LOCALID`), one rank per GPU, OMP_NUM_THREADS=1. Multi-node: same, with `-N <nodes> --ntasks-per-node=4`.

## 3. Scaling probe (fresh, short; one job per node count; nothing in prod/)
Mesh 536 x 128 x 128. Split laterally only (full radial columns, as every E1 run so far); one rank per GPU:
| arm | nodes x H200 | meshblock (nx1 x nx2 x nx3) | MeshBlocks | per GPU |
|---|---|---|---|---|
| G2  | 1 x 2 | 536 x 64 x 64 (deck) | 4  | 2 |
| G4  | 1 x 4 | 536 x 64 x 64 (deck) | 4  | 1 |
| G8  | 2 x 4 | 536 x 64 x 32 | 8  | 1 |
| G16 | 4 x 4 | 536 x 32 x 32 | 16 | 1 |
- Override with `meshblock/nx2=.. meshblock/nx3=..` (keys exist). If a lateral split FATALs at init, report the message and try the radial
  split instead (G8r 268 x 64 x 64, G16r 268 x 64 x 32); say which one you used.
- Each arm: fresh, `time/nlim=100 time/ndiag=1 output2/dt=1e30 output3/dt=1e30 output4/dt=1e30 output5/dt=1e30`, -t 0:30:00.
  Run G2 and G4 in the same 1-node job if convenient. Do the G8/G16 jobs only if the queue estimate (below) is under ~6 h; otherwise
  report the estimate and skip.
- Report per arm: rc; FATAL / NaN / DIVERGED / NON-CONVERGED / NOT ADMISSIBLE / DT_RETRY counts; dt (expect 727.27 s); Picard mean/max
  (fresh G2 ref: mean 3.07 max 6); init time; **s/cycle over cycles 10-100** (G2 ref: 0.570 s with dt_retry on 2 H200); peak memory
  per GPU (nvidia-smi sampled during the run, e.g. `nvidia-smi --query-gpu=memory.used --format=csv -lms 500` in the background on
  every node; G2 ref 59.9 GB/GPU incl. the 4.5 GB/rank dt_retry checkpoint); the `vet_gd_halo_ipc:` line (ON / MPI fallback) and the
  `dt_retry on:` line; hst last-row max rel diff vs G2 (not bitwise across layouts; expect round-off ~1e-9).
- Cluster limits (paste the raw output, briefly summarised):
  `scontrol show partition gpu`; `sacctmgr -P show qos format=Name,MaxWall,MaxTRESPU,MaxTRESPerJob,MaxJobsPU,MaxSubmitPU,GrpTRES,Priority`;
  `sacctmgr -P show assoc user=$USER format=Account,Partition,QOS,MaxJobs,MaxTRES,MaxWall,GrpTRES`; `sshare -U`; H200 nodes
  idle/mixed/allocated (`sinfo -p gpu -N -o '%N %G %t'`, h200 only). Our earlier snapshot: partition gpu MaxTime 14 d, QOS
  debug/long/normal, 18 H200 nodes (exclude hpc-sm-01-09).
- Queue: `sbatch --test-only` start estimates for 1 / 2 / 4 nodes x 4 H200 (and 1 x 2 H200) at -t 24 h, 48 h and the longest wall the
  QoS allows (normal and long QoS). Note that long -t can delay backfill: report both.
- Choose the layout with the best time-to-solution = queue wait + run time, run time = 33,300 cycles (2.4192e7 / 727 s) x s/cycle x
  ~1.3 (later Picard growth, outputs). Ref: 2 H200 at 0.57 s/cycle = ~5.3 h + margin, so ONE long link may cover the whole run on any
  layout. Prefer the fewest links with the longest walltime allowed (user: "potentially run long jobs"), but not a layout whose queue
  wait exceeds its speed gain. State the choice and the arithmetic in the NOTE. The layout is then FIXED for the whole chain
  (restarts are bitwise only for the same ranks/GPU type/binary).

## 4. Gate on the chosen layout (one job, -t ~1:30)
1. Production deck = the section-1 deck, and if the chosen meshblock is not 64 x 64, a second read-only copy with ONLY the `<meshblock>`
   nx2/nx3 lines edited (record its md5 and the diff). Production uses that file; no meshblock override on the command line.
2. Smoke + restart-bitwise: `smoke_e1.sh` from `e1s_1014/scripts` with the chosen launch and deck: C (fresh 60) vs L (30 + 30,
   hook restart) and Coff (ipc false); then `smoke_gate.sh`: CHAINCMP / CMPLAST / every hst row BITWISE; Coff BITWISE. GATE PASS writes
   VERIFIED. Any non-bitwise result: STOP and report (do not start production).
3. Background check (same job or right after): arm BG, fresh through link_e1.sh with `E1C_NCYC=1400` (t ~1.0e6 s), production
   cadence (bins 1e5). Report and check:
   - `python3 pulse_ana.py run <BG dir> BG_pulse.npz > BG_pulse.txt` (`docs/handover/pulse-1014/pulse_ana.py` + bin_convert.py beside it);
     PU S3_268 reference at 5e5 / 9e5 s: vrV 0.076 / -0.039 km/s, vrM 0.58 / -0.19, rsurf 1.009 / 1.011 R_ph.
   - user.hst (cols: 6 L_top, 17 Mdot_top, 20 Min_top, 22 Picard, 23 w_mlt, 26 P_sponge): L_top/L0 ~1.0; no systematic inflow at the top.
   - shell-mean <v_r> above 1.1 R_ph and in the sponge (r > 7.29e13 cm): |v_r| < ~1 km/s, no inward-moving floor/atmosphere layer;
     <rho>(r)/<rho>(t=0) above 1.1 R_ph quiet; fofc_report at-floor 0; NaN/FATAL 0.
   PASS = all quiet like the S3 reference. Anything else (infall, growing outer velocities, floor activity): STOP and report.
4. Write `prod/GO_e1` yourself (text: "USER 10-14 e1s request, TASK-2026-10-14-caltech-e1s, gate PASS <date>") only after 2 and 3 pass.

## 5. Production
- `prod/run.cfg`: BIN 5b7e1bc5 (md5 6252e58e), IN = production deck, KEYS empty (no test keys), TLIM=2.4192e7, STOP_NC=1, E1C_GATE=1,
  E1C_VERIFIED=<gate smoke>/VERIFIED, LAUNCH = the chosen one-word launch, E1C_DIR = e1s_1014/scripts. Fresh t = 0 (prod/rst empty).
- Links: `link_e1.sh` via the patched hook, the longest walltime allowed for the chosen layout/QoS (state it), first link + 2 backup
  links `afterany`, name `agcE1s`. Each link restarts from the newest VALID rst (rst_pick); it never starts fresh once rst exist; on
  rc != 0 / FATAL / NaN / (non-retry) DIVERGED / NON-CONVERGED the hook writes STOP and the chain halts. Do not restart after a STOP:
  report and wait.
- Push a STARTED note (jobs, layout, walltime, start estimate, first s/cycle) as soon as link 1 runs.

## 6. Early evolution report (push to the NOTE as each time is passed; table per time)
At t ~5e6, 7e6, 9e6, 1.25e7 s (and the end): from `pulse_ana.py run prod` (bins) + user.hst: vrV, vrM, vr105, KEe, rsurf, L_top/L0
(col 6 / row 0), KE_int (col 11), w_mlt (col 23), Picard mean, DT_RETRY count so far, s/cycle, dt. Push `prod_pulse.txt` and
`agcar3d.user.hst` to `docs/handover/e1s-1014/caltech_out/`. Compare against:
| t (s) | run | vrV | vrM | vr105 | KEe | rsurf | L_top/L0 | KE_int |
|---|---|---|---|---|---|---|---|---|
| 5e6    | PU S3_268 (seed 1e-3, reduced) | 0.42 | 1.83 | -0.14 | 1.0e40 | 1.011 | 0.994 | 5.9e42 |
| 7e6    | PU S3_268 | 0.29 | 3.91 | -0.89 | 9.5e40 | 1.025 | 0.960 | 2.4e43 |
| 9e6    | PU S3_268 | 1.13 | 3.77 | 3.39  | 2.0e41 | 1.056 | 0.909 | 6.8e43 |
| 1.25e7 | PU S3_268 | 0.03 | 0.93 | 0.39  | 1.2e41 | 1.072 | 0.794 | 1.1e44 |
| 5e6    | your old E1 (seed 1e-2) | 5.76 | 11.75 | 11.95 | 3.7e41 | 1.129 | | |
| 7e6    | your old E1 | 1.27 | -0.18 | -1.05 | 3.4e40 | 1.227 | | |
| 9e6    | your old E1 | -20.9 | -20.7 | -25.5 | 7.1e41 | 1.173 | | |
| 1.24e7 | your old E1 | -45.9 | -26.4 | -50.6 | 3.5e41 | 1.074 | | |
(velocities km/s, 1.0-1.2 R_ph; rsurf = outermost r/R_ph with <rho> > 1e-12; S3_268 is 268 x 32 x 32, so expect the trend, not the numbers.)
Flag at once: any infall (vrV, vrM < -5 km/s), rsurf > 1.15, DT_RETRY bursts, or Picard mean > 20.

## NOTE skeleton (`NOTE-2026-10-14-caltech-e1s.md`)
Build/deck md5 | scaling table (arm, layout, s/cycle, GB/GPU, Picard, counts) | limits + queue estimates | chosen layout + arithmetic |
gate (BITWISE lines, BG check) | STARTED (jobs, walltime) | evolution tables.
