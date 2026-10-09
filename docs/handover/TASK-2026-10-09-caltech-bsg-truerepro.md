# URGENT TASK for Caltech: stop the He giant N897, run the BSG true reproduction (2 x H200)

**User decision 10-09: "queue BSG on Caltech; we can cancel the He giant work".** From: viper. To: Caltech
(Resnick, H200 141 GB, CUDA sm_90). Everything this task needs is on the fork branch **`bsg-files-1009`** (this file +
`docs/handover/bsg-hrdet-1009/`); the code is fork branch **`mem-1009`**. Expected wall ~22 h (57k cycles x ~1.4 s).

## 1. Stop the He giant N897 cleanly (do not delete anything)
- Run dir `/resnick/groups/carnegie_poc/jingze/hegiant_1007/n897` (chain 4242926-29, binary 7adb12d3 since link 2).
- If a link is RUNNING: wait for its next rst if it is near, otherwise use the newest rst already written. Then
  `scancel` the running link and **all PENDING links** of the chain (`squeue -u $USER` fresh in the same command;
  `scancel --state=PENDING` for the queued ones).
- Record in the NOTE (section 6): last t, cycle, rst path, job ids cancelled. Keep the run dir, rsts, outputs.
  The user may resume it later from that rst.

## 2. Build (report md5)
- Code: `git fetch <fork remote> mem-1009` (github jing-ze-ma/athenak) -> commit **6c5d8fb2c** (= half-range VET 98835d99 + `vet_gd_twin_lowmem`).
  Use a separate worktree/snapshot of that commit; no destructive git in a shared checkout.
- Target he_star_m1 for H200, same toolchain as your N897 build (CUDA 12.9, gcc 13.2, hpcx OpenMPI):
  `-D PROBLEM=he_star_m1 -D Athena_ENABLE_MPI=ON -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_HOPPER90=On
  -D CMAKE_CXX_FLAGS=-ffp-contract=off` (+ your host arch flag and nvcc_wrapper as for f3a66907/7adb12d3). The
  `-ffp-contract=off` is host code only, as in Delta's `build_delta.sh he_gpu_nofma`; device code keeps nvcc defaults.
  Incremental in your existing he_star_m1 build dir is fine. Name e.g. `athena_gpu_he_star_m1_6c5d8fb2`; report md5.

## 3. Files and input
```
git fetch <fork remote> bsg-files-1009
git archive FETCH_HEAD docs/handover/bsg-hrdet-1009 | tar -x -C <work>
bash <work>/docs/handover/bsg-hrdet-1009/SETUP.sh <abs files dir>   # gunzip IC + tables, md5 check, writes inputs
```
- **Production input = `<files dir>/bsg3d_truerepro2_hr_lm.athinput`** = viper's production input
  `bsg3d_truerepro2_hr.athinput` (md5 2d81553e034d3286af8de8065dd55a6c) with only (a) the 3 file paths -> @FILES@ and
  (b) one added line `vet_gd_twin_lowmem = true` in `<rad_m1>` after `vet_gd_twin_fuse`. Keep it in the FILE:
  AthenaK FATALs on a command-line key that is not in the input (DeltaAI job 3348542).
- Grid 256^3 (meshblock 256x64x64 -> 16 MeshBlocks), tlim **4.96e6 s** (57.4 d), outputs as on viper (hst 1e3 s,
  4 bin 0.5 d, log 1e4 s, rst 1 d). Do not change physics keys; cfl 0.3.
- Helpers in the bundle: `rst_info.py RST` prints "t ncycle tlim" and the `outputN/last_time=...` keys to pass on a
  restart (so outputs do not repeat); `truerepro2_viper_link.sh.txt` = viper's link script (FRESH if no rst, else
  newest rst + last_time keys; DONE at tlim; STOP on rc != 0 / FATAL / NaN; md5 check of BIN) to adapt.

## 4. Layout
- 1 node, **2 MPI ranks x 8 MeshBlocks, 1 rank per H200**, `CUDA_VISIBLE_DEVICES=$SLURM_LOCALID` per local rank
  (wrapper), `--cpu-bind=cores`. DeltaAI measured exactly this layout with this binary commit and the lowmem key
  (NOTE-2026-10-09-deltaai-bsg-mem2.md on this branch, job 3348593): **65.0 GiB/GPU peak, 1.37 s/cycle on 2 GH200**,
  bitwise identical to the non-lowmem run. Fits 141 GB with margin.

## 5. Smoke, then production
- **Smoke** (separate dir, same binary + input + command-line keys as production, plus `time/nlim=10`): rc 0, no
  FATAL / NaN / NON-CONVERGED in the log, GPU memory peak (nvidia-smi), s/cycle (cycles 3-10). Production only if clean.
- **Production**: fresh start (t = 0) from the input, restartable chained links at your wall limit (dependency
  afterany, each link restarts from the newest rst with `rst_info.py` last_time keys), stop rule: rc != 0, FATAL,
  NaN or NON-CONVERGED -> STOP file, no further links. DONE when the newest rst t >= 4.96e6.

## 6. Notes back (push to `bsg-files-1009`, file `NOTE-2026-10-09-caltech-bsg-truerepro.md`, append)
- **Immediately when production STARTS**: start time, job ids of the chain, smoke job id + its s/cycle and GPU peak,
  binary md5, He giant stop record (section 1). **viper then cancels its pending viper/Raven BSG chains: first start
  wins**, so push this as soon as link 1 is RUNNING.
- After each link: t, cycle, dt, s/cycle, FATAL/NaN counts, newest rst.
