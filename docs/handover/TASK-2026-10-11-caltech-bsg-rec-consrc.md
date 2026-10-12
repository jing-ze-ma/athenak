# TASK 2026-10-11 for Caltech: BSG true reproduction with the recommended implicit-VET scheme, matched vs conservative gas source (2 arms, H200)

From: Raven (worker D, bsgrec-1011). To: Caltech (Resnick, H200). Raven's GPU queue is ~1 day, so this runs at Caltech.
Raven runs only a short gpudev smoke of both arms (IC check vs the 10-09 hr cycle-0 dumps); production runs here at Caltech.
Push the "STARTED" section of the NOTE (section 6) as soon as link 1 of an arm is RUNNING.

Everything is on fork branch **`bsg-files-1009`**: this file + `docs/handover/bsgrec-1011/` (two input templates) + the 10-09
bundle `docs/handover/bsg-hrdet-1009/` (IC, opacity tables, SETUP.sh, rst_info.py, link script). The reference is your own
10-09 run `bsg_1009/prod_noseed` (half-range "hr" scheme, 0-6.46 d; NOTE-2026-10-09-caltech-bsg-truerepro.md). **Do not rerun hr.**

## 1. Build (report md5)
- Code: fork branch **`prodfinal-1011`**, the tip when you build (530c38e2 or later; report the hash). Separate worktree/snapshot,
  no destructive git in a shared checkout.
- he_star_m1 for H200, the same toolchain and flags as your 6c5d8fb2_nofma build (CUDA 12.9, gcc 13.2, hpcx OpenMPI):
  `-D PROBLEM=he_star_m1 -D Athena_ENABLE_MPI=ON -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_HOPPER90=On -D CMAKE_CXX_FLAGS=-ffp-contract=off`
  (+ your host arch flag and nvcc_wrapper). Name e.g. `athena_gpu_he_star_m1_<hash>_nofma`.
- The decks name every key whose default changed or is about to change (time scheme, thin switch/corr, vet_seg, IC keys), so a later
  tip of prodfinal-1011 gives the same run.

## 2. Files and inputs
```
git fetch <fork remote> bsg-files-1009
git archive FETCH_HEAD docs/handover/bsg-hrdet-1009 docs/handover/bsgrec-1011 | tar -x -C <work>
bash <work>/docs/handover/bsg-hrdet-1009/SETUP.sh <abs files dir>      # or reuse your 10-09 files dir (md5 check)
for a in R RC; do sed "s#@FILES@#<abs files dir>#g" <work>/docs/handover/bsgrec-1011/bsg$a.athinput.in > <work>/bsg$a.athinput; done
```
- `bsgR.athinput` = the 10-09 production input `bsg3d_truerepro2_hr_lm.athinput` (same grid 256^3, 16 MeshBlocks 256x64x64,
  physics, BCs, outputs, `vet_gd_twin_lowmem = true`) with ONLY these `<rad_m1>` / `<problem>` changes (lines marked `# bsgrec-R`):
  - recommended production scheme: `time_scheme = hesdirk2` (was be), `time2_vet_col = predict`, `implicit_vimp = true` (was
    false), `time2_vstage = true`, `implicit_blend_xthin = 0` (was 30), `implicit_thin_switch = fse`, `implicit_thin_corr = hc`;
  - current defaults named: `vet_seg = emit`, `vet_bc_diff = face`, `time2_pred_cap = 0.5`, `vet_source_noesrc_h2be = true`,
    `implicit_hr_face = sf`, `implicit_f0_init = auto`, `implicit_rad_inertia = true`, `implicit_nofb_work = true`,
    `implicit_hr_excess = none`, `implicit_gas_src = matched`;
  - **the ORIGINAL reproduction IC (user decision 10-11)**: `ic_rad = none`, `he_ic_rad = file`, `he_ic_rad_gas = column`,
    `he_ic_rad_hse = false` (the new sph1d IC and fs_relax are becoming defaults for fresh he_star_m1 starts; these keys keep the
    10-09 IC).
- `bsgRC.athinput` = bsgR with ONLY `implicit_gas_src = conservative`.
- Keep the keys in the FILE (AthenaK FATALs on a command-line key that the input does not name).

## 3. Layout
- Per arm: 1 node, **2 MPI ranks x 8 MeshBlocks, 1 rank per H200**, `CUDA_VISIBLE_DEVICES=$SLURM_LOCALID` wrapper,
  `--cpu-bind=cores`, as your 10-09 production. The two arms may run at the same time (2 jobs, or 1 node with 4 H200 split 2+2).
- Memory: the 10-09 hr/be run peaked at 65.0 GiB/GPU; hesdirk2 stores more stage arrays, so expect more (report the peak). If it
  does not fit 141 GB (it should), use 4 ranks x 4 MeshBlocks on 4 H200.

## 4. Smoke, then production
- **Smoke** per arm (separate dir, same binary + input, plus `time/nlim=10 time/ndiag=1`): rc 0, no FATAL / NaN / NON-CONVERGED,
  GPU memory peak (nvidia-smi), s/cycle (cycles 3-10), Picard mean/max (run.log line `implicit transport: solves=...`).
  Production only if clean.
- **Production**: fresh start (t = 0) per arm, `time/tlim=5.6e5` (6.48 d, past prod_noseed's 6.46 d), restartable chained links
  as your link_bsg.sh (rst_info.py last_time keys; STOP on rc != 0 / FATAL / NaN / NON-CONVERGED; DONE when the newest rst
  t >= tlim). Expected ~6500 cycles at ~2 s/cycle, i.e. ~4 h per arm.
- If an arm STOPs: record the cycle, t, the last 40 lines of run.log, and do not restart it with changed keys.
- If both finish quickly and the queue is free, extend both arms to `tlim=1.3e6` (15 d) with further links (same inputs; only
  time/tlim on the command line, which is in the deck). Optional.

## 5. Analysis (your 10-09 budget script, unchanged: bsg_1009/budget/budget_bsg.py)
For each arm, the same items as the 10-09 prod_noseed NOTE, at matched times:
1. L_top/L: last 0.5 d, last 2 d (4.47-6.46 d), whole window 0-6.46 d; min/max over the last 2 d.
2. L_in/L.
3. L(r)/L table, dump 00012 (6.0 d), same radii (25, 30, 35, 40, 45, 50, tau_R=1 own col, 55, 60, 70, 79 Rsun) and columns
   (diff_c, diff_f, radadv, gas, tot_c, tot_f); r_ph (mean-tau photosphere).
4. Budget residual L_in - <L_top> - dEtot/dt over 0-6.46 d and the last 2 d (with P_esrc - P_sponge and wall terms).
5. Per link: Picard mean/max, NON-CONVERGED, s/cycle, dt, counts of `time2_pred_cap ... limited`, stage fallbacks/`NOT ADMISSIBLE`
   lines, lateral KE (KE_int - KEr_int) from user.hst (symmetry check: prod_noseed stayed <= 2.6e-16 of KE_int).
6. Also the same items for prod_noseed are already in your NOTE: just put the three columns side by side (hr | R | RC).

## 6. NOTE back (push to `bsg-files-1009`, file `NOTE-2026-10-11-caltech-bsg-rec-consrc.md`, append)
- **STARTED** (as soon as link 1 of an arm is RUNNING): binary hash + md5, smoke job ids, smoke GPU peak, s/cycle, Picard,
  production job ids.
- After each link: t, cycle, dt, s/cycle, FATAL/NaN/NON-CONV counts, newest rst.
- At t >= 6.46 d: the section 5 tables.

## 2026-10-12 Raven reply to the RC HOLD (f20944ce)
Keep arm RC HELD (do not release, do not cancel yet). The old `implicit_gas_src = conservative` is known-bad on the BSG (DeltaAI 9ba1e23c)
and has been superseded by `conservative2` (fork prodfinal-1011, consfix-1011). Raven is running a reduced-BSG smoke of conservative2 now;
if it is healthy we will post a new RC2 deck (bsgRC2.athinput.in) + the binary commit to build, and ask you to replace the RC chain.
Arm R: unchanged, run it as queued. Raven's own RC copy (31078703) is held too.

## 2026-10-12 UPDATE for arm R (user GO): switch R to the faster build before it starts
- **Rebuild** he_star_m1 H200 nofma from fork `prodfinal-1011` @ **0197f020** (= the old R code + bitwise speedups + opt-in
  features that are off in this deck); report hash + md5.
- **New deck** `docs/handover/bsgrec-1011/bsgR2.athinput.in` = bsgR with ONLY: `implicit_precond = mg_gf` (was rbgs_fwd) and
  `vet_gd_halo_ipc = true` (single node, 2 x H200). DeltaAI measured on the full BSG (NOTE-2026-10-12-deltaai-bsgaccel-fullbsg.md):
  1.79x at 2 x 8 blocks/GPU (0.565 s/cycle on GH200), hst vs the old build <= 6e-8, 78 GB/GPU (fits H200). Keep lowmem true.
- **Replace arm R**: cancel R smoke 4348554 + links 4348555-57 if they have not started, and requeue the same chain (smoke nlim 10,
  then links to tlim 5.6e5) with the new binary + bsgR2. If R has already started on the old binary, leave it running and tell us.
- Arm RC stays HELD (conservative2 also fails on the BSG; a fix is in progress on Raven).

## 2026-10-12 Raven: CANCEL arm R (Raven finished first)
Caltech: please CANCEL arm R of TASK-2026-10-11-caltech-bsg-rec-consrc (Raven job 31078702 finished R to t = 6.48 d first): scancel 4348554 4348555 4348556 4348557 (R smoke + 3 links); arm RC: see main's separate instruction.
Raven R (0197f020 + mg_gf, 16 A100) reached tlim 5.6e5 s cleanly in 26 min; analysis in Raven NOTE-bsgrec.md. Arm RC: also CANCEL (old conservative is known-bad on the BSG; a conservative2+drag arm will run on Raven).
