# TASK for FREYA (MPA cluster, MPCDF): BSG arm 2 production (same run as viper, Caltech, DeltaAI, Raven)

**CANCELLED 10-02 06:30 (user): BSG runs on viper; see NOTE-2026-10-02-bsg-drop-remote.md. Do not start this task.**

Read first: NOTE-2026-10-02-sync-bsg.md (common setup and reporting), NOTE-2026-10-02-bsg-code-ready.md (commit and
gate numbers), TASK-2026-10-01-caltech-bsg.md (full background, sections 1-5). This TASK only adds the FREYA specifics.
FREYA cannot be reached from viper (ssh blocked), so this runs from a session started on FREYA itself.

## 1. Machine (MPCDF docs, MPA-FREYA)
- Login freya01-04.bc.mpcdf.mpg.de. GPU partition **p.gpu.ampere**: 11 nodes x 4 NVIDIA A100-PCIE-40GB (sm_80);
  also p.gpu (P100/V100: do NOT use). Test partition p.test (30 min). Max 24 h. Scratch /freya/ptmp.
- Check `sinfo -p p.gpu.ampere`, `module avail gcc cuda openmpi`, and the GPU gres name (`scontrol show node <ampere node>`
  -> Gres=gpu:a100:4 or similar) before writing job scripts.

## 2. Code and input (identical on all sites)
- Build fork rt-integration at **30bf6c03** exactly (git checkout 30bf6c03; kokkos submodule). Problem he_star_m1:
  `cmake -D PROBLEM=he_star_m1 -D Athena_ENABLE_MPI=ON -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_AMPERE80=On
   -D CMAKE_CXX_COMPILER=$PWD/kokkos/bin/nvcc_wrapper -D CMAKE_BUILD_TYPE=Release` (+ a CPU build for the column gate).
  Record both md5s. (Raven, the other A100 site, uses the same flags; see its NOTE when it appears.)
- Input: docs/handover/bsg_1001_bundle (SETUP.sh <run dir>) -> bsg3d_arm2.athinput with exactly two changes:
  `<output1> dt = 1000.0` and `implicit_eos_cache` removed (comment it out).

## 3. Gates (p.test or a short p.gpu.ampere job)
- 4a column bsg_col_arm2 CPU vs 1 GPU: dt_end 1.243014e+02, L_top/L_in end 0.999466 range [0.999374, 0.999467],
  0 FATAL/nan, CPU vs GPU <= 1e-6 in hst; Picard mean report-only (2.994 or 3.000).
- 4b 3-D smoke bsg3d_arm2 on **1 node x 4 A100** (4 MeshBlocks, 1 per GPU = production layout), 20 cycles:
  dt ~88.18 s, Picard ~3, 0 FATAL/nan; log shows implicit_src_stable 1, implicit_face_weight distance,
  he_bc_hse_flux face, implicit_det_reduce 1. **Check GPU memory fits 40 GB** (nvidia-smi during the run or Kokkos
  memory report); if it does not fit, STOP and report (do not change the mesh layout without the user).

## 4. Production
- Guarded jobs as on viper (bsg_1001/DUAL_1002.md on viper, summarised): one run dir on /freya/ptmp; every job uses the
  same script: (1) if another RUNNING job of the same name exists (squeue), exit 0 (older job id wins a tie);
  (2) if DONE exists or newest rst t >= tlim, exit 0; (3) if STOP exists, exit 1; (4) fresh start if no rst, else
  restart from the newest rst with every <outputN>/last_time = floor(t/dt)*dt (t is 232 bytes after <par_end> in the
  rst header); code wall limit = job limit - 20 min; write STOP on FATAL/NaN. No scancel inside jobs.
- Submit 2 x 24 h + 2 x 6 h jobs of that script on p.gpu.ampere, 1 node, 4 GPUs, 4 MPI ranks (1 GPU each),
  tlim 4.96e6 s (57.4 d). One node at ~1.3 s/cycle (A100 guess) needs ~20-30 h for ~56,000 cycles: 2 links likely.

## 5. Report (push to the fork, never force; fetch first)
NOTE-2026-10-0x-freya-bsg.md: build md5s, gate numbers, GPU memory, s/cycle, job ids, queue wait (squeue --start),
then the progress NOTEs listed in NOTE-2026-10-02-sync-bsg.md (5/10/20/30 d). Nobody cancels another site's run.
