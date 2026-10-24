# TASK viper -> caltech (USER): run the Arepo2-IDORT GPU tests on the Caltech H200 node

**Why.** The Raven GPU queues are slow, and the user wants the single-node Arepo2-IDORT GPU tests to run on Caltech.
This continues NOTE-2026-10-24-raven-arepo2-idort-gpu.md (8c0086de6).

**No Arepo source goes through git.** Arepo2-IDORT is a private Bitbucket repo. Do not commit or push Arepo source,
diffs or patches to this GitHub fork, and do not push to Bitbucket. All Arepo branches are local only.

## What the user copies to Caltech (bundle, 192 MB)
Raven: `/raven/ptmp/jinma/arepo2_gpu/caltech_bundle_1024/`

| file | md5 | content |
|---|---|---|
| arepo_idort_gpu_1024.bundle | 905ec3cb6ad069aa1e491658956b3da8 | git bundle, complete history, 5 branches (below) |
| data/g25_inputs.tgz | 2348f0db65cc644446240951bd518235 | g25 ICs (g25_red, g25_full), tables, param_long_bd_raven.txt, run_g25.sh, tools/g25cmp.py, README_C.md |
| data/recipes.tgz | 817c23901e06d776b408b359311455ab | env.sh, build1.sh, HYPRE 2.32.0 tarball, design/, RESULTS.md, port/STAGE0.md, BENCH1.md, gpu1/ (PROGRESS.md, build/gate/time scripts, cfg, tcases), nltree/RESULTS_NLTREE.md + scripts |
| data/submodules.tgz | b4dc14b5539273f9b25622b8310ca851 | lib/safeintegral, lib/git-format-staged (public GitHub) |

Branches:
- `Arepo2-IDORT-dev` = `raven-gpu-1024` = bd36de23f3: the base.
  - Steps A and B (the GPU recipe works) need no source change.
- `idort-gpu-s0` = 0929dfe290: stage 0 (level diagnostics, FP pin, dump mode).
- `idort-gpu-s1` = 6cef4d1d69: **the island-resident GPU IDORT solve** (`IDORT_GPU`, `IdortGpu 1`).
  - It is bitwise equal to the CPU path on all gates.
- `nltree-g25` = 1108335a71: `FORCETREE_NODELOCAL(_GPU)` with `MULTIPLE_NODE_SOFTENING`, so g25 gravity can run on the GPU tree.
  - The force gate is still pending on Raven.

Unpack with `git clone -b idort-gpu-s1 arepo_idort_gpu_1024.bundle arepo2`, then `tar xzf submodules.tgz -C arepo2`.

## Build
Follow `env.sh`, `build1.sh` and `gpu1/build_g1.sh`. They are Raven scripts with absolute `/raven` paths; adapt the paths.
- **HYPRE:** build 2.32.0 CPU-only from the tarball.
- **Config:** use `gpu1/cfg/*` (CPU arm and GPU arm). Change the CUDA arch from sm_80 to **sm_90**.
- **Raven build flags:** `-fmad=false` for nvcc and `-ffp-contract=off` on the host. Keep both: bitwise CPU = GPU depends on them.
- **Layout:**
  - `NumIslandsPerNode` = number of GPUs (4);
  - ranks per node divisible by 4 (uneven islands hang);
  - `CPUsUsedPerGPU` = ranks / 4.

## Tests, in order (1 node x 4 H200)
1. **Gates**, CPU vs GPU, from `gpu1/gate_g1.sh` with the cases in `gpu1/tcases`.
   - Every gate must end `ALL SAME`.
   - Report any difference with the file and column.
2. **g25_red timing:** `gpu1/time_g1.sh`, CPU arm and GPU arm, same node, interleaved.
   - Report per full sync and per top step (the phase table format in `gpu1/PROGRESS.md`).
   - Raven reference (A100): CPU 7.19 s, GPU 3.90 s per full sync.
3. **nltree-g25:**
   - build the `examples/nltree_mns_3d` and `nltree_mns_3d_cuda` regression tests and run them (CHECKS / VERIFY);
   - then g25_red with the full GPU recipe (`VORONOI_CELL_GPU` + `FORCETREE_NODELOCAL_GPU`) against the CPU build;
   - see `nltree/RESULTS_NLTREE.md`.

## Answer first (NOTE caltech->viper)
- Node: CPU cores, GPUs, MPI and CUDA versions, and whether CUDA IPC works between the GPUs (the island halos use it).
- Whether IDORT tests can share the node with run12 or must wait for it. run12 has priority; do not disturb it.
- Then the results of 1-3. Do not change any defaults or preconditioner settings without the user.
