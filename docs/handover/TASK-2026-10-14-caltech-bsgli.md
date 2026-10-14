# TASK raven -> Caltech: BSG PRODUCTION with Linhao's IC + inner BC on 1 node x 2 H200 (fin4-linhao-bc @ 8e018272), smoke now, chain AUTO after the gate

From the Raven coordinator (worker LC, bsglic-1014). USER decision: queue the new Linhao-IC BSG production on Caltech as a THIRD copy,
next to Raven bsgLI (31143829-31) and Viper bsgLIv (12190480-82). Both of those are already released.
This is a PRODUCTION run that the user asked for, like your BSG fin3 run. It is not a test task, so the "skip GPU tasks more than ~1 d out"
rule does not apply. Report the queue estimate at once (section 6, QUEUED).

## !!! NO GO FILE: the user has APPROVED this production !!!
**Caltech may release the production chain AUTOMATICALLY once its own smoke gate passes. No GO file is needed, and none will be pushed.**
- Link 1 is submitted with `--dependency=afterok:<smoke job>` (section 4).
- The smoke job runs `smoke_gate.sh` itself at its end:
  - GATE PASS -> it writes `<SMOKE_DIR>/VERIFIED`, exits 0, and link 1 starts when the scheduler allows;
  - GATE FAIL -> it exits 1, so link 1 never starts.
- The hook's safety net: a link refuses to run without `<SMOKE_DIR>/VERIFIED` (it writes `NOGO` and exits 1). Never rely on it.
  Never hand-write VERIFIED.
- If the gate FAILs: cancel the queued links (`scancel`), STOP and report (section 6). Do not rerun with changed keys, except for
  the halo_ipc case in section 2.

Background: NOTE `/raven/ptmp/jinma/bsglinhao_1014/NOTE.md` on Raven (LI). Summary: prod P FIN3 BSG (your bsgP_seed_fin3 run) with
- the IC replaced by L. Ma's bsg_star.cpp initial state:
  - `problem/he_ic_linhao`, his MESA profile as supplied (rmax 52.49 Rsun);
  - the 5 % rho perturbation instead of the 1e-6 T seed;
- his exact inner BC, `problem/he_bc_linhao = true` (XB): bd_data hydro ghosts + ghost-bath radiation.

FIN4 is the method base (fin4-1014, your FIN4 suite NOTE 1217fa43).

## 1. Build (report md5)
- Fork branch **`fin4-linhao-bc-8e018272` @ `8e018272`** (8e0182723da2...). It is a frozen snapshot.
  - Do NOT build the moving branch `fin4-linhao-bc`: its head is now bb735538, a later diagnostic commit that Raven/Viper production does not run.
  - 8e018272 = FIN4 + Linhao IC + 37798491 (the exit-abort fix) + the L_in history column (diagnostic only).
  - Raven A100 build md5 0f431a61, Viper MI300A bacea0a8.
- Same toolchain and flags as your 9972a56d / dad00933 / FIN4 builds: `-D PROBLEM=he_star_m1 -D Athena_ENABLE_MPI=ON
  -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_HOPPER90=On -D CMAKE_CXX_FLAGS=-ffp-contract=off`, default `Athena_M1_VGD_NMAX`.
- Build in a compute job (or incrementally from your FIN4 build dir), as you usually do. Name it `athena_gpu_he_star_m1_8e018272_nofma`.
- Report its md5 and `git rev-parse HEAD`.

## 2. Files and deck (`docs/handover/bsgli-caltech-1014/`, check MD5SUMS)
- One directory `<F>` with:
  - your 10-09 BSG files, already on Caltech (the deck still reads them for the ghost column, gravity and opacities):
    `ic_ma2026_f0.txt` (246a666e), `rosseland_ma2026_x0.7_z0.008_padded.txt` (450bc0c1), `planck_ma2026_x0.7_z0.008_padded.txt` (e18ca805);
  - `linhao_profile_full.txt` = gunzip of `docs/handover/bsgli-val-1014/linhao_profile_full.txt.gz`. That file is already on this branch
    (17 MB, shipped for DeltaAI; gz md5 661f76d7). The unpacked md5 MUST be **41edc225**.
  - Nothing else is shipped or read: the EOS is the ideal gas, and he_grid_dump writes into the run dir.
  - Linhao's raw `MESA_profile_combined_cgs.txt` is NOT needed: the deck names `linhao_profile_full.txt`.
- Deck: `bsgli-caltech-1014/bsgLI_bc_caltech.athinput.in` (md5 532fc7af):
  ```
  sed "s#@FILES@#<F>#g" bsgLI_bc_caltech.athinput.in > bsgLI_bc_caltech.athinput
  ```
  Change no other key. Make the deck read-only and record its md5 in the NOTE. It is Raven's production deck
  `prod/LI/bsgLI_full_fin3_bc.athinput` (md5 d37e861a; the one Raven bsgLI and Viper bsgLIv run) with only these changes:
  - the 4 file paths are templated (`he_ic_linhao_file`, `he_ic_file`, `he_opac_table`, `he_planck_table`);
  - your single-node BSG keys, as in your bsgP_seed_fin3 deck:
    - `vet_gd_halo_ipc = true`;
    - `vet_gd_twin_lowmem = true`.
  It is the DeltaAI validation deck `bsgli-val-1014/bsgLI_bc_deltaai.athinput.in` with only comment lines changed.
  - It keeps everything else unchanged: the FIN3 keys, he_ic_linhao (amp 0.05, rho_thr 1e-11, rng 0), `he_bc_linhao = true`,
    tlim 4.96e6, cfl 0.3, outputs (hst 1e3 s, bins 4.32e4 s, rst 8.64e4 s) and meshblock 256 x 64 x 64.
- **halo_ipc and lowmem: kept only if bitwise.**
  - `vet_gd_halo_ipc`: the smoke checks it on H200 (arm Coff, part of the gate).
    - If Coff is not BITWISE, or you get FATAL "halo spin timed out", the gate FAILs.
    - Then set `vet_gd_halo_ipc = false` in YOUR deck and redo the whole smoke with that deck in a NEW smoke dir.
    - Say so in the NOTE, and keep it false for the whole chain. This is the only allowed deck edit.
  - `vet_gd_twin_lowmem`:
    - DeltaAI mem2 showed it BITWISE vs false, and your 9972a56d BSG ran it.
    - false does not fit at 8 blocks/GPU (~140 GB/GPU), so no smoke arm tests it.
- Scripts, kept together in one dir (python3 + numpy):
  - `link_bsgli.sh`: the production link. It is your E1 link `link_e1.sh` with this hook, and it reads `RUN_DIR/run.cfg` at job start.
  - `bsgli_hook.sh`: `e1chain_hook.sh` with E1C_ -> BLC_. Its only change: the gate needs ONLY the smoke VERIFIED, with no GO file.
  - `rst_pick.py`, `prune_rst.sh` (9509a75c), `chaincmp.py`, `cmplast.py`, `bin_convert.py`: unchanged copies of e1caltech-1013.
  - `smoke_bsgli.sh` + `smoke_gate.sh`: the smoke and the DATA gate (section 3).
  - The smoke writes LAUNCH into run.cfg with `printf %q`, so a multi-word LAUNCH works. This fixes your E1 smoke 4452278 issue.

## Layout (fixed for the WHOLE chain)
- **1 node x 2 H200, 2 MPI ranks x 8 MeshBlocks (256 x 64 x 64) each, 1 rank per GPU, OMP_NUM_THREADS=1**. This is your BSG fin3 layout.
  Expect ~86-89 GB per H200 (your 9972a56d BSG: 85.6-88.7 GB).
- Restarts are bitwise only for the same ranks x threads x GPU type x binary. Never change the rank count, the GPU type, LAUNCH or
  the binary inside the chain. The smoke must use exactly the production launch.
- Launch: the one your BSG and E1 production used, i.e. `e1_1013/e1launch.sh`
  (= `exec srun --mpi=pmix -n 2 -c 8 --cpu-bind=cores e1wrap.sh "$@"`, with `CUDA_VISIBLE_DEVICES=$SLURM_LOCALID`).
  - Set it as `LAUNCH` for the smoke (5th argument) AND in the production run.cfg.
  - Check the run.log line `vgdrot-1011 vet_gd_halo_ipc: ON (node-local CUDA IPC ...)`.

## 3. Smoke (now): one job, 1 node x 2 H200, 1 h 15 min wall
```
S=$(sbatch --parsable --export=ALL,BLC_DIR=<abs path of bsgli-caltech-1014> <your BSG flags for 1 node x 2 H200> --ntasks=2 \
      -t 01:15:00 -J bsgLIsmk bsgli-caltech-1014/smoke_bsgli.sh <SMOKE_DIR> <BIN> <BIN_MD5> <F-deck path> <LAUNCH>)
```
It is the recipe of your E1 smoke, and of Raven 31143531 and Viper 12190373/74 on the LI deck (both BITWISE). Everything runs
through `link_bsgli.sh`. The test-only keys all exist in the deck:
- ndiag 1, hst 200 s;
- hydro_w/m1 bins 1000 s, m1_vet/m1_face bins 2500 s;
- log off; rst as the deck (t = 0 + link-end).

Arms:
- **C (fresh-start smoke)**: fresh, 60 cycles (dt 87.5 s -> t ~5250 s). Check and report:
  - rc 0, 0 FATAL / NaN / DIVERGED, NON-CONVERGED 0;
  - the `he_star_m1: he_ic_linhao: ...` line and the `he_bc_linhao ghost ...` lines;
  - dt (expected 87.5 s after the first cycle);
  - Picard mean / max;
  - the halo_ipc line, the init time;
  - the GPU peak memory (the smoke samples nvidia-smi into `<SMOKE_DIR>/gpumem.csv`);
  - s/cycle over cycles 10-60;
  - the last user.hst row.
- **L (restart BITWISE: C 60 vs L 30+30)**:
  - link 1 is fresh to cycle 30 and ends with the link-end rst;
  - link 2 runs the SAME link script; the hook picks the newest valid rst and continues to 60.
  - L must equal C BITWISE.
- **Cadence check** (chaincmp.py): every hst DATA row of C and every scheduled dump of C must be in the chain once, with none missing,
  none extra, and every one bitwise.
- **Coff**: C with `rad_m1/vet_gd_halo_ipc=false`. Its last dumps and hst must be BITWISE vs C.
- **Smoke GATE = smoke_gate.sh** (run by the smoke job at its end on the compute node; the login node runs nothing heavy).
  - It compares DATA only, never whole bin files:
    - hst rows without the '#' header lines;
    - each bin variable's data array, read with bin_convert (not the header);
    - the logs (FATAL / NaN / DIVERGED / NON-CONVERGED in all arms; the he_ic_linhao and he_bc_linhao lines in C).
  - It requires >= 4 compared dumps.
  - Your 10-13 whole-bin false positive cannot happen.
  - It writes `<SMOKE_DIR>/VERIFIED` only on GATE PASS. On FAIL, report the output of `gate_chaincmp.txt` / `gate_cmplast.txt` /
    `gate_coff.txt` plus the last 40 lines of each run.log.
- After the smoke, run the prune test on the smoke's L/rst, as for your BSG and E1: hard links in a dummy dir, keep 3. It runs on the
  login node and is light.
  - success deletes only the oldest;
  - rc != 0, STOP, FATAL, a truncated newest file and an unreadable newest file delete nothing.

## 4. Production chain (submit NOW, right after the smoke; it releases itself on the smoke's PASS)
- Run dir `<RUN>` (fresh, empty). `<RUN>/run.cfg`:
  ```
  BLC_DIR=<abs path of bsgli-caltech-1014>     # under sbatch $0 is a spooled copy, so the scripts' dir must be named
  BIN=<.../athena_gpu_he_star_m1_8e018272_nofma>
  BIN_MD5=<md5>
  IN=<F-deck path>
  KEYS=""
  TLIM=4.96e6
  STOP_NC=1
  BLC_GATE=1
  BLC_VERIFIED=<SMOKE_DIR>/VERIFIED
  LAUNCH=<the same launch as the smoke, e.g. /abs/path/e1_1013/e1launch.sh>
  ```
- Links: 12 h each, 1 node x 2 H200, `--ntasks=2`, ONE job name for all links (the hook's GUARD uses it). Submit link 1 with afterok on
  the smoke, and links 2..5 with afterany on the previous link:
  ```
  j=$(sbatch --parsable <flags> --ntasks=2 -t 12:00:00 -J bsgLIc --dependency=afterok:$S bsgli-caltech-1014/link_bsgli.sh <RUN>)
  for k in 2 3 4 5; do j=$(sbatch --parsable <flags> --ntasks=2 -t 12:00:00 -J bsgLIc --dependency=afterany:$j bsgli-caltech-1014/link_bsgli.sh <RUN>); done
  ```
- What each link does (the hook):
  - **DONE** file -> exit 0. **STOP** file -> exit 1. **NOGO** (no VERIFIED) -> exit 1.
  - **GUARD**: if an older link of the same name is RUNNING -> exit 0.
  - **Restart source**: the newest VALID rst in `<RUN>/rst` (rst_pick.py: header, the M1ROTAV1 history block, the leading size).
    - Fresh `-i` only if `<RUN>/rst` has no rst at all.
    - If rst files exist but none is valid -> STOP (never a fresh start over a run).
  - **Restart keys**: every `outputN/last_time` from the rst time (the output cadence is kept), plus `time/tlim=4.96e6`.
  - **STOP after the link**:
    - srun rc != 0;
    - FATAL;
    - NaN in the last 60 lines;
    - `DIVERGED`;
    - **NON-CONVERGED > 0** (STOP_NC=1).

    A STOP writes `<RUN>/STOP` with the reason; nothing is pruned or chained after it (later links exit at once).
  - **Prune** (keep the newest 3 rst) after a successful link only, with prune_rst.sh; log in `<RUN>/rst_prune.log`.
    A BSG rst is ~9.2 GB (your link 1), so about 28 GB are kept.
  - **DONE** when the newest valid rst has t >= 4.96e6.
- **How many links (estimate):**
  - dt = 87.5 s (your link 1: 538 of 539 printed cycles), so 4.96e6 / 87.5 = **~56,700 cycles**.
  - Speed: your 9972a56d BSG on the same layout ran **1.57 s/cycle** on average (1.2 early -> 2.0-2.3 late). FIN4 costs **+14 %**,
    so expect **~1.79 s/cycle** on average (1.37 early -> ~2.6 late).
  - So ~28 h of stepping on average. At ~41,000 s of stepping per 12 h link, that is **3 links**. If the late ~2.6 s/cycle holds
    from early on, it is **4 links**. Linhao's IC starts convection at ~4 d (KE ~40x the old IC), so Picard may rise sooner than in
    your old-IC run.
  - Queue **5 links** for margin. The rest exit DONE at once and cost nothing.
  - After each link, from your measured s/cycle and dt:
    `links_left = ceil( (4.96e6 - t) / dt_mean * s_per_cycle / 41000 )`.
    If that exceeds the links still queued, add links with afterany on the last one.
- Never change BIN, the deck, the rank count or LAUNCH inside the chain. If something must change, STOP and report first.

## 5. Priority vs your other jobs, and Do not
- **Your other jobs continue unchanged and keep their places:**
  - the He giant FIN3 chain (4454171-76);
  - BSG fin3 link 2 (4434962; let it finish, as the 10-14 Raven NOTE says);
  - the E1 restart arms R0 4483862 / R1 4483863 (DR's diagnosis).
  - Do not cancel, hold, re-nice or reorder any of them for bsgLIc.
  - The smoke and the bsgLIc links queue normally beside them. Fairshare will delay bsgLIc; that is expected.
  - If you must choose what to submit first, the order is:
    1. the running jobs;
    2. the E1 R0/R1 arms;
    3. the He giant chain;
    4. the bsgLIc smoke + chain.
- Do not touch any other run's dirs or rst (keep all old BSG and E1 data).
- Do not edit the deck beyond the @FILES@ sed. The halo_ipc fallback in section 2 is the only exception, and you must report it.
- Do not start production without VERIFIED. Do not hand-write VERIFIED.

## 6. Report: docs/handover/NOTE-2026-10-14-caltech-bsgli.md (push to bsg-files-1009)
- **QUEUED** (as soon as you read this):
  - the binary md5 + `git rev-parse HEAD`;
  - the deck md5 after the sed, and the profile md5 (41edc225);
  - the smoke job id and the link job ids;
  - the `sbatch --test-only` start estimate for 1 node x 2 H200 12 h.
- **SMOKE**:
  - the gate verdict (CHAINCMP / CMPLAST / Coff / HST lines, GATE PASS or FAIL);
  - the section 3 numbers: s/cycle over cycles 10-60, GPU peak memory per H200, init time, Picard, dt;
  - the halo_ipc line, the he_ic_linhao / he_bc_linhao lines, the prune test.
- **STARTED**: when link 1 is RUNNING (job id, node, fresh t = 0).
- **Per link**:
  - t, cycle, dt, s/cycle;
  - the `implicit transport: ... Picard ... NON-CONVERGED=` line;
  - the FATAL / NaN / DIVERGED counts, `NOT ADMISSIBLE`, and `fofc_report` at-floor > 0;
  - the newest rst, the rst_prune.log lines, and links_left.
- **Interim at t = 1e6, 2e6, 3e6, 4e6 s and at the end**: push `bsg3d.user.hst` + `bsg3d.hydro.hst` (they are small).
  Keep every bin and the last rst on disk for Raven.
