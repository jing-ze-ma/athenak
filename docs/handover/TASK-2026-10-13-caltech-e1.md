# TASK raven -> Caltech: AG Car E1 PRODUCTION on 1 node x 2 H200 (dic-prod-1013 @ dad00933), smoke now, production on GO

From the Raven coordinator (worker EC, e1caltech-1013). USER decision: the Raven copy of AG Car E1 (agcmE1 31106649) moves to
Caltech, because the Raven gpu queue now estimates about 20 h for it. Viper keeps its own E1 copy (agcvE1 12164005).
This is a PRODUCTION copy that the user asked for, like your BSG fin3 run, not a test task: the "skip GPU tasks more than ~1 d out" rule does
not apply. Do report the queue estimate at once (section 6, QUEUED), because DeltaAI is the fallback if Caltech cannot start within a few hours.
Your running BSG chain (4434961 -> 4434962 -> 4434963) has priority: do not touch it. E1 needs its own 2 H200.

Background: TASK-2026-10-12-caltech-agcmodel-E1E0.md (the E1 physics: deep inner wall 83 Rsun, 536 x 128 x 128, long scaffold ramp,
tlim 2.4192e7 s). Since then E1 got the FIN3 method (as your BSG), DI's discrete-equilibrium IC (he_ic_equil, v6 "f1La"), a denser
output cadence, and DR's wb_rmax line (below).

## !!! WAIT FOR GO !!!
**Do NOT start (or release) the production chain until the file `docs/handover/GO-e1-caltech.txt` exists on bsg-files-1009.**
The coordinator pushes it only after DR's restart-bitwise chain tests of the wb_rmax deck pass on Raven (31133340) and Viper (12181730).
- Sections 1-3 (build, files, smoke) do NOT wait: do them now.
- Submit the production links HELD (`sbatch --hold`, section 4), so they are queued but cannot start.
- Poll the branch (your 15-min poller). When GO-e1-caltech.txt appears AND your smoke gate wrote VERIFIED: copy the GO file into the
  run dir as `GO_e1`, then `scontrol release` link 1 (and the later links).
- If the GO file instead says NO-GO or names a different deck md5: do not release; follow what it says.
- The link script refuses to run without GO_e1 + VERIFIED (it writes NOGO and exits 1). That is a safety net only; never rely on it.

## 1. Build (report md5)
- Fork branch **`dic-prod-1013` @ `dad00933`** = your production 9972a56d + 17a22820 + dad00933. The two new commits are DI's opt-in
  problem/he_ic_equil (default off is bitwise to 9972a56d on A100 and MI300A). wb_rmax is already in 9972a56d (hydro/wb_rmax).
- Same toolchain and flags as your 9972a56d build (incremental from that build dir is fine): `-D PROBLEM=he_star_m1
  -D Athena_ENABLE_MPI=ON -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_HOPPER90=On -D CMAKE_CXX_FLAGS=-ffp-contract=off`,
  default `Athena_M1_VGD_NMAX`.
- Name `athena_gpu_he_star_m1_dad00933_nofma`. Report its md5 (and `git rev-parse HEAD` of the source = dad00933acc5...).

## 2. Files and deck
All in this branch; check `docs/handover/e1caltech-1013/MD5SUMS`.
- Put these three files in ONE directory `<F>` (you may already have them from the 10-12 E1/E0 TASK; check the md5):
  - `docs/handover/agcmodel-1012/caltech/ic_agcar_A_ge_deep70.txt` (md5 e81deb51, 1 MB);
  - `docs/handover/agcarA_files/rosseland_ext2_gs98_x0.36_z0.02.txt` (b28e97c5);
  - `docs/handover/agcarA_files/planck_ext2orion_gs98_x0.36_z0.02.txt` (493f9d01).
  No other file is read (the EOS is the built-in table; he_grid_dump writes agcar3d.x1grid.txt into the run dir).
- Deck: `docs/handover/e1caltech-1013/agcA_E1_final_wbrmax.athinput.in` (md5 62115d17):
  ```
  sed "s#@FILES@#<F>#g" agcA_E1_final_wbrmax.athinput.in > agcA_E1_final_wbrmax.athinput
  ```
  Change no other key. Make the deck read-only and record its md5 in the NOTE.
  - It is the Raven FINAL E1 deck (`agcmodel_1012/inp/agcA_E1_final.athinput`, md5 2d407b51, what Raven/Viper production runs)
    PLUS one `<hydro>` line, the user-approved DR fix (drift-1013), after `wb_cache_every`:
    `wb_rmax = 2.9715e13` (x1 well-balancing only up to 1.1 R_ph; above it the unsupported floor layer gets plain reconstruction).
    With the Raven paths put back, the file is byte-identical to DR's candidate deck (md5 16b0bba7).
  - Only the 3 file paths are templated. Everything else is unchanged: FIN3 keys (vet_rotate_avg 4 etc.), DI IC keys (he_ic_equil,
    v6 f1La), hst 1e4 s, hydro_w/m1 bins 1e5 s, rst 5e5 s, tlim 2.4192e7, meshblock 536 x 64 x 64.
  - `vet_gd_halo_ipc = true` is KEPT. It is CUDA-generic (node-local CUDA IPC; it falls back to MPI by itself if a partner is off-node)
    and your BSG production already runs it on 2 H200. It is speed-only and bitwise (A100: ipc on vs off BITWISE, Raven 31120050).
    The smoke checks it again on H200 (arm Coff). If it fails on H200 (FATAL "halo spin timed out", or Coff not BITWISE), set it to
    false in your deck, say so in the NOTE, and keep it false for the whole chain.
- Scripts in `docs/handover/e1caltech-1013/` (keep them together in one dir; python3 + numpy):
  - `link_e1.sh`: the production link. It is Raven's agcmE1 link script (rc_raven.sh) plus the chain hook, ported. It reads `RUN_DIR/run.cfg`
    at job start.
  - `e1chain_hook.sh`: Raven's e1chain_raven.sh (validated BITWISE on A100 + MI300A). It chooses the restart source, adds the STOP rules,
    prunes (keep 3), writes DONE at tlim, and holds the GO gate.
  - `rst_pick.py`: the newest VALID rst, plus the outputN/last_time keys that keep the output cadence (as your rst_info.py).
  - `prune_rst.sh`: the same file (9509a75c) as `prodswitch-1013/prune_rst.sh`, which your BSG already uses.
  - `smoke_e1.sh`, `smoke_gate.sh`, `chaincmp.py`, `cmplast.py`, `bin_convert.py`: smoke and DATA gate (section 3).
  - If you prefer your own BSG link script, port the hook's rules into it exactly (sections 4-5), and run the smoke through THAT script.

## Layout (fixed for the WHOLE chain)
- **1 node x 2 H200, 2 MPI ranks x 2 MeshBlocks (536 x 64 x 64) each, 1 rank per GPU, OMP_NUM_THREADS=1.** It is the same split as Viper
  (2 MI300A x 2 blocks, 55.7 GB per GPU), so about 56 GB of 141 GB per H200 is expected.
- Restarts are bitwise only for the same ranks x threads x GPU type. Never change the rank count, the GPU type or the binary inside the chain.
  The smoke must use exactly the production launch.
- The launch: `link_e1.sh` defaults to `srun -n 2 --cpu-bind=cores` (Raven: KOKKOS_MAP_DEVICE_ID_BY=mpi_rank with all GPUs visible).
  If your BSG production uses a different launch (e.g. the `CUDA_VISIBLE_DEVICES=$SLURM_LOCALID` wrapper), set the same one in run.cfg
  as `LAUNCH="..."` (everything before the binary), for the smoke AND production. Check the run.log line
  `vgdrot-1011 vet_gd_halo_ipc: ON (node-local CUDA IPC ...)`. "not possible ... MPI path" is also correct (bitwise), only slower: report it.

## 3. Smoke (now, no GO needed): one job, 1 node x 2 H200, ~1 h wall
```
sbatch --export=ALL,E1C_DIR=<abs path of e1caltech-1013> <your BSG flags for 1 node x 2 H200> --ntasks=2 -t 01:00:00 -J agcE1smk \
  e1caltech-1013/smoke_e1.sh <SMOKE_DIR> <BIN> <BIN_MD5> <F-deck path> ["<LAUNCH>"]
bash e1caltech-1013/smoke_gate.sh <SMOKE_DIR>        # after the job; login node is fine
```
This is PS's recipe (Raven 31120638 and Viper 12172704-06, both BITWISE on dad00933 + this deck without wb_rmax). Everything runs through
`link_e1.sh`. Test-only keys, which all exist in the deck: ndiag 1, hst 4e3 s, bins 1.5e4 s, log output off.
- **C (fresh-start smoke)**: fresh, 60 cycles. Check and report:
  - rc 0, 0 FATAL / NaN / DIVERGED, NON-CONVERGED 0;
  - the he_ic_equil lines (1 round; done);
  - dt (Raven: 727.27 s on every cycle);
  - Picard mean / max (Raven A100: 2.42 / 6);
  - the halo_ipc line, the init time, the GPU memory;
  - s/cycle over cycles 10-60 (Raven 4 A100: ~0.36-0.37);
  - the last user.hst row (Raven compares it with its A100 C).
- **L (restart BITWISE: C 60 vs L 30+30)**:
  - link 1 is fresh to cycle 30 and ends with the link-end rst;
  - link 2 runs the SAME link script, the hook picks the newest valid rst (cycle 30), and it continues to 60.
  - L must equal C BITWISE.
- **Cadence check**: in the same comparison (chaincmp.py), every hst DATA row of C (11 per file) and every scheduled dump (8) must be
  in the chain once: none missing, none extra, every one bitwise. The restart is from the link-end rst, which is the production case
  (the old Finalize schedule bug lost one row and one dump exactly there).
- **Coff**: C with `rad_m1/vet_gd_halo_ipc=false`. Its last dumps and hst must be BITWISE vs C.
- **Smoke GATE = smoke_gate.sh**. It compares DATA, never whole bin files:
  - hst rows without the '#' header lines;
  - each bin variable's data array, read with bin_convert (not the header);
  - the logs: FATAL / NaN / DIVERGED / NON-CONVERGED, and the he_ic_equil line;
  - it requires >= 4 compared dumps, so an empty comparison cannot pass.

  Your BSG false positive (whole-bin cmp flagging the nlim / file-counter / offset header metadata) cannot happen with it. It writes
  `<SMOKE_DIR>/VERIFIED` only on GATE PASS.
  - The production run.cfg points `E1C_VERIFIED` at that file.
  - Do not hand-write VERIFIED. If the gate FAILs, STOP and report: the output of gate_chaincmp.txt / gate_cmplast.txt, plus the
    last 40 lines of each run.log.
- Also run the prune test on the smoke's L/rst (hard links in a dummy dir, as you did for BSG; keep 3):
  - success: deletes only the oldest;
  - rc != 0, STOP, FATAL, truncated newest, unreadable newest: deletes nothing.

## 4. Production chain (queue HELD now; release only on GO, see the top)
- Run dir `<RUN>` (fresh, empty). `<RUN>/run.cfg`:
  ```
  E1C_DIR=<abs path of e1caltech-1013>     # under sbatch $0 is a spooled copy, so the scripts' dir must be named
  BIN=<.../athena_gpu_he_star_m1_dad00933_nofma>
  BIN_MD5=<md5>
  IN=<F-deck path>
  KEYS=""
  TLIM=2.4192e7
  STOP_NC=1
  E1C_GATE=1
  E1C_VERIFIED=<SMOKE_DIR>/VERIFIED
  # LAUNCH="<same launch as the smoke>"   (only if you set one there)
  ```
- Links: 12 h each, 1 node x 2 H200, `--ntasks=2`, ONE job name for all links (the hook's GUARD uses it). Use link 1 `--hold`,
  then links 2..N with `--dependency=afterany:<previous>`:
  ```
  j=$(sbatch --parsable --hold <flags> --ntasks=2 -t 12:00:00 -J agcE1c e1caltech-1013/link_e1.sh <RUN>)
  for k in 2 3; do j=$(sbatch --parsable <flags> --ntasks=2 -t 12:00:00 -J agcE1c --dependency=afterany:$j e1caltech-1013/link_e1.sh <RUN>); done
  ```
- What each link does (the hook):
  - **DONE** file -> exit 0. **STOP** file -> exit 1. **NOGO** (no GO_e1 or no VERIFIED) -> exit 1.
  - **GUARD**: if an older link of the same name is RUNNING -> exit 0.
  - **Restart source**: the newest VALID rst in `<RUN>/rst` (rst_pick.py: header, history block, leading size). Fresh `-i` only if
    `<RUN>/rst` has no rst at all. If rst files exist but none is valid -> STOP (never a fresh start over a run).
  - **Restart keys**: every `outputN/last_time` from the rst time (the output cadence is kept), plus `time/tlim=2.4192e7`.
  - **STOP after the link**:
    - srun rc != 0;
    - FATAL;
    - NaN in the last 60 lines;
    - `DIVERGED`;
    - **NON-CONVERGED > 0** (STOP_NC=1, as your BSG fin3 links).

    A STOP writes `<RUN>/STOP` with the reason; nothing is pruned or chained after it.
  - **Prune** (keep the newest 3 rst) after a successful link only: prune_rst.sh, as your BSG section 4b; log in `<RUN>/rst_prune.log`.
    An E1 rst is ~5.1 GB, so about 15 GB are kept on disk. Output bins: about 242 dumps per stem x ~170 MB = ~80 GB for the run.
  - **DONE** when the newest valid rst has t >= tlim.
- **How many links (estimate):**
  - dt = 727.27 s (Raven; it may shrink late in the ramp). So 2.4192e7 / 727.27 = **~33,300 cycles**.
  - Expected 2 H200 speed: ~0.4-0.6 s/cycle (4 A100 0.36, 2 MI300A 0.63).
  - So ~3.7-5.6 h of stepping, plus output and IC: **1 link if the speed and dt hold**.
  - Queue **3 links** for margin. The rest exit DONE at once and cost nothing.
  - After link 1, from your measured s/cycle and dt:
    `links_left = ceil( (2.4192e7 - t) / dt_mean * s_per_cycle / 41000 )`.
    If that exceeds the links still queued, add links with afterany on the last one.
- Never change BIN, the deck, the rank count or LAUNCH inside the chain. If something must change, STOP and report first.

## 5. Do not
- Do not touch your BSG chain or any other run.
- Do not edit the deck beyond the @FILES@ sed (the halo_ipc fallback in section 2 is the only exception, and it must be reported).
- Do not start production without GO_e1 + VERIFIED.

## 6. Report: docs/handover/NOTE-2026-10-13-caltech-e1.md (push to bsg-files-1009)
- **QUEUED** (as soon as you read this): binary md5, deck md5, smoke job id, held link ids, and the `sbatch --test-only` start estimate
  for 1 node x 2 H200 12 h. The coordinator decides on DeltaAI from this line.
- **SMOKE**: the section 3 numbers, the gate output (CHAINCMP / CMPLAST / HST lines, GATE PASS/FAIL), Coff verdict, halo_ipc line,
  the prune test.
- **STARTED**: when link 1 is RUNNING (job id, node, fresh t = 0). Raven then cancels its own copy.
- **Per link**: t, cycle, dt, s/cycle, the `implicit transport: ... Picard ... NON-CONVERGED=` line, FATAL / NaN counts, the counts of
  `fofc_report` at-floor > 0 and `NOT ADMISSIBLE`, the newest rst, the rst_prune.log lines, and links_left.
- **Interim at t = 4e6, 8e6, 1.2e7, 1.6e7 s and at the end**: push agcar3d.user.hst + agcar3d.hydro.hst (they are small). Keep every bin
  and the last rst on disk for Raven.
