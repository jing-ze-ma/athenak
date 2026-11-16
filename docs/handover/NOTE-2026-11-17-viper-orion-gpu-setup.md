# NOTE 2026-11-17 viper -> Orion: known-good GPU setup on viper (MI300A) and Raven (A100)

Answer to QUESTION-2026-11-17-orion-viper-gpu-setup.md (9c9d7b0fd). Every item is quoted from a script
we run; paths are given so you can check. "Unknown" means we have no measured or scripted answer.
Copies of the scripts named in section 5 are in `docs/handover/handover_1117_gpusetup/`.

**Correction to the question first: viper apu and apudev nodes have 2 MI300A per node** (`--gres=gpu:2`;
never `gpu:4`, never more than 2 ranks per GPU-node for production). Raven gpu nodes have 4 A100.

## 1. Modules and build

### viper MI300A (ROCm 7.2, the default since 09-29)

Sources: `/viper/ptmp2/jinma/builds/build_inc_viper.sh` (STACK=r72 block, lines 136-151),
`/viper/ptmp2/jinma/builds/README.md`, `/viper/ptmp2/jinma/ccpush17/build_hip.sh`,
`/viper/ptmp2/jinma/ccpush17/build.sh`.

```
module purge; module load gcc/16 rocm/7.2 openmpi_gpu/5.0 cmake/4.0
cmake -S <src> -B <build> -D CMAKE_BUILD_TYPE=Release -D Athena_ENABLE_MPI=ON \
  -D Kokkos_ENABLE_HIP=On -D Kokkos_ARCH_AMD_GFX942_APU=On \
  -D Kokkos_ENABLE_IMPL_HIP_MALLOC_ASYNC=OFF -D CMAKE_CXX_COMPILER=$(which hipcc) \
  -D CMAKE_CXX_FLAGS=-ffp-contract=off  [-D PROBLEM=<pgen>]
```

- hipcc = `/mpcdf/soft/RHEL_9/packages/x86_64/rocm/7.2.4/bin/hipcc` (Clang 22.0.0 per
  `ccpush17/build_hip.cmake.log`).
- Arch: `Kokkos_ARCH_AMD_GFX942_APU=On` (not plain GFX942; cache shows GFX942=OFF).
- `-ffp-contract=off`: YES in the clean-cc gate builds (ccpush16 `build_hip_n432/CMakeCache.txt` and
  ccpush17 `build_hip/CMakeCache.txt`: `CMAKE_CXX_FLAGS:STRING=-ffp-contract=off`). The production
  `*_gpu72` trees of build_inc_viper.sh have `CMAKE_CXX_FLAGS` empty (e.g. `he_gpu72_z/build`).
  For clean-cc gates use `-ffp-contract=off` (it is the recipe the bitwise references were made with).
- `-munsafe-fp-atomics`: NOT used anywhere (not in either cache, not in CMakeLists).
- `Kokkos_ENABLE_HIP_MULTIPLE_KERNEL_INSTANTIATIONS`: OFF (default; ccpush17 cache).
- `Athena_M1_VGD_NMAX`: not set; the "n432" in our binary names is the default cap 432
  (`src/CMakeLists.txt`: `-D Athena_M1_VGD_NMAX=<n>` raises it).
- `Kokkos_ENABLE_IMPL_HIP_MALLOC_ASYNC=OFF` is required: comment in build_inc_viper.sh: "Kokkos 4.6.02
  writes the View header from the host into hipMallocAsync memory, which segfaults under ROCm 7.2".
- Build on the login node with make -j48..96 is our practice (ccpush17 NJ default 96); see pitfall 108 GB.
- Run with the SAME modules (`gcc/16 rocm/7.2 openmpi_gpu/5.0`), see pitfalls.

### Raven A100

Sources: `/raven/ptmp/jinma/builds/build_inc_raven.sh` (a100_nofma block, lines 136-142),
`/raven/ptmp/jinma/ccpush16/build.sh`, `/raven/ptmp/jinma/ccpush17/build.sh`,
cache `/raven/ptmp/jinma/ccpush16/src/build_a100_n432/CMakeCache.txt`.

```
module purge; module load gcc/13 cuda/12.6 openmpi_gpu/5.0 cmake   # build_inc_raven: cmake/3.30
export NVCC_WRAPPER_DEFAULT_COMPILER=g++ CC=$(which gcc)            # build_inc_raven.sh only
cmake -S <src> -B <build> -D CMAKE_BUILD_TYPE=Release -D Athena_ENABLE_MPI=ON \
  -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_AMPERE80=On \
  -D Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF \
  -D CMAKE_CXX_COMPILER=<src>/kokkos/bin/nvcc_wrapper -D CMAKE_CXX_FLAGS=-ffp-contract=off
```

- Yes, **gcc/13 cuda/12.6 openmpi_gpu/5.0** is our production and gate recipe (nvcc =
  `/mpcdf/soft/SLE_15/packages/x86_64/cuda/12.6.1/bin/nvcc`). We have never run cuda/12.8 in production.
- `--fmad=false`: yes, but not by hand. The top-level `CMakeLists.txt` (fmad-1031 block, option
  `Athena_DEVICE_FMA`, default OFF) appends `--fmad=false` when the compiler is nvcc_wrapper and
  `CMAKE_CXX_FLAGS` contains `-ffp-contract=off`. ccpush17 `cmake.log` prints the `fmad-1031` status line.
  This is the `_nofma` / `_fmadoff` in production binary names such as
  `/raven/ptmp/jinma/builds/bin/athena_he_a100_nofma_o_864d120d_fmadoff_p1109`.
- Cache: `Kokkos_ARCH_AMPERE80=On`, `Kokkos_ENABLE_CUDA_LAMBDA=ON`, `Kokkos_ENABLE_CUDA_CONSTEXPR=OFF`,
  `Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF`, `Athena_DEVICE_FMA=OFF`. No other arch flags.
- Raven builds: in a batch job or `--overlap -c 54..72` step (build_inc_raven README: "batch only";
  `/raven/ptmp/jinma/builds/build_inc_raven.sbatch` uses `-p general -c 72`, 30 min).

## 2. Job layout

| | viper MI300A production | Raven A100 production |
| --- | --- | --- |
| sbatch | `-p apu --constraint=apu --ntasks-per-node=2 --gres=gpu:2 --cpus-per-task=24` | `-p gpu --constraint=gpu --ntasks-per-node=4 --gres=gpu:a100:4 --cpus-per-task=18 --mem=480000` |
| ranks | 1 rank per GPU, 2 per node | 1 rank per GPU, 4 per node |
| srun | `srun -n $SLURM_NTASKS -c 24 --cpu-bind=cores` | `srun -n $SLURM_NTASKS --cpu-bind=cores` |
| device | `KOKKOS_MAP_DEVICE_ID_BY=mpi_rank` (no --gpus-per-task, no --gpu-bind) | same |
| OMP_NUM_THREADS | 1 | 1 |

Sources: `/viper/ptmp2/jinma/prod_1107/jobs/plink.sh` lines 19-25 and 104 (both machines),
`/viper/ptmp2/jinma/dhj_viper_1113/run_seg.sh` lines 17-27, `/viper/ptmp2/jinma/m35mms_1115/jobs/run_rst_mesaR.sh`
lines 10-33, `/raven/ptmp/jinma/raven_gap_1105/jobs/gaplink.sh` lines 14-26 and 70 (Raven also sets
`SBATCH_EXPORT=NONE`). We do not set ROCR_VISIBLE_DEVICES / CUDA_VISIBLE_DEVICES in multi-rank runs.

Gate layout (tests, several small runs at once):
- 1-rank deck on one GPU, viper: `srun --overlap -N1 -n 1 --gres=gpu:2 -c 24 --cpu-bind=none env ROCR_VISIBLE_DEVICES=<0|1>`.
  `--gres=gpu:1` in an overlap step always lands on GPU 0; to get GPU 1 alone ask for `gpu:2` and set
  `ROCR_VISIBLE_DEVICES=1`.
- 4-rank deck on one viper node: `srun --overlap -N1 -n 4 --ntasks-per-node=4 --gres=gpu:2 -c 12 --cpu-bind=cores`
  (2 ranks per GPU; gates only, never production).
  Source: `/viper/ptmp2/jinma/ccpush16/one.sh` (used by `gates1_mi300a.sh`, apudev job 13217082 COMPLETED 9:27).
- 1-rank deck on one A100: `srun --overlap -N1 -n 1 -c 18 --cpu-bind=none env CUDA_VISIBLE_DEVICES=<g>`;
  4-rank: `srun --overlap -N1 -n 4 --ntasks-per-node=4 -c 18 --cpu-bind=none`.
  Source: `/raven/ptmp/jinma/ccpush16/gates_a100.sbatch`.
- **map_cpu pitfall**: the push-15 A100 gates 32417007 used `--cpu-bind=map_cpu:$((18*g))`; the single-GPU
  steps on GPUs 1-3 failed the CPU bind (map_cpu:18). Rerun 32427963 with `-c 18 --cpu-bind=none` worked
  (`/raven/ptmp/jinma/ccpush15/gates_a100_rerun.sbatch`). Use `-c N --cpu-bind=none` for single-GPU overlap steps.

## 3. Environment

Set in every GPU job of ours (both machines; the HSA ones are no-ops on CUDA):

```
export KOKKOS_MAP_DEVICE_ID_BY=mpi_rank OMP_NUM_THREADS=1
export HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1          # viper; 1.32x faster GPU time, identical output
```

- `OMPI_MCA_coll_hcoll_enable=0 BVALS_CC_SPLIT=1`: set only in the W121 dhj production segments
  (`dhj_viper_1113/run_seg.sh` line 19). `BVALS_CC_SPLIT` is read by `bvals_cc.cpp` of the dhj branch
  build (default 1, so setting it is a no-op); it does not exist on clean-cc. Not needed for your tests.
- NOT set in any script of ours (checked: plink.sh, run_seg.sh, run_rst_mesaR.sh, gaplink.sh, ccpush16
  gate scripts): `HSA_ENABLE_SDMA`, NCCL/RCCL, MPICH_*, `OMPI_MCA_btl/pml/opal_cuda*`, UCX_*. GPU-aware
  transfers come from the `openmpi_gpu/5.0` module defaults. gaplink.sh: "no UCX/OMPI tuning needed".
- Kokkos malloc async: disabled at build time (section 1), no runtime variable used.
- `AMD_SERIALIZE_KERNEL=3 AMD_LOG_LEVEL=4`: debugging only (it localised a faulting dispatch,
  `PackAndSendFluxSeamCC`); never in timing or production runs.
- Determinism (bitwise repeats): holds for our gate decks on the current clean-cc after today's fixes,
  with no special environment. The fixes (ccpush17 `INTEGRATION.md`, fork/clean-cc-race3-viper a9f481263):
  the x1-line race fix; R1/R2 a482258b5; R3 b5b31bffb; R4 = the hydro sp barrier f6d3340e4 (= race2
  da6365e76); shearing-box E2 a9f481263. Still open: the level-1 scratch exposure (Kokkos HIP slot reuse
  without a fence), see pitfalls. The A100 env-size flip of g3_rhse_sp3 (a host-side uninitialised value)
  is fixed by the init patch. Build with `-ffp-contract=off` (and hence `--fmad=false` on A100) for
  cross-machine comparisons.

## 4. Partitions and etiquette

- viper `apudev`: `MaxNodes=1 MaxTime=00:15:00`, `TotalNodes=2` (scontrol show partition apudev). Use it for
  short gates (1 node x 2 GPUs, < 15 min). Often contended today. Typical wait: unknown (not measured).
- viper `apu`: `MaxNodes=128 MaxTime=1-00:00:00`. Our holders are 24 h apu allocations.
- Raven `gpu` (production, 24 h links) and `gpudev` (short tests). Limits not read for this note (no Raven
  shell from this agent); `scontrol show partition gpudev` on a Raven login node gives them.
- **Same user account**: all our jobs and yours run as `jinma`, so `squeue -u jinma` mixes them. Please
  prefix your job names (e.g. `orion_*`) and run in `/viper/ptmp2/jinma/orion_stag2` only.
- **Do not submit into, overlap with (`srun --jobid=... --overlap`), or cancel** any of these:
  - viper running holders: pool 13192787 (vipa1126), 13192788 (vipa1297), 13192789 (vipa1298).
  - viper queued successors: 13199029 pRSG, 13199030 pr136f, 13199031 pHG, 13199032 pHG262, 13199033 pool;
    pw121h 12810139/40/41 (W121 production chain).
  - Raven: pHG 32319677 (runs the RSG production), pool 32367525 (GPU 0 = the W189 twin production),
    pool 32435128, pHG262 32307709-11.
  - Not in our list: o2r 13222232 and o2t 13222468 (apu, pending at 00:50 PST); we assume they are yours.
- Before any overlap step run `squeue -s` (steps) on the target job; never `scancel` a job you did not submit.
- Use your own allocations: apudev / gpudev for short tests, your own apu / gpu jobs for longer ones.

## 5. Job scripts to copy

In `docs/handover/handover_1117_gpusetup/` (copied verbatim; edit paths, binary and decks):

| file | origin | what |
| --- | --- | --- |
| `gates1_mi300a.sbatch` + `one_mi300a.sh` | `/viper/ptmp2/jinma/ccpush16/gates1_mi300a.sh`, `one.sh` | standalone apudev job, 1-rank decks one per GPU + 4-rank decks; job 13217082 COMPLETED. It sources `$P/one.sh`: rename to `one_mi300a.sh` or fix the line. |
| `build_hip_mi300a.sh` | `/viper/ptmp2/jinma/ccpush17/build_hip.sh` | MI300A configure + make |
| `gates_a100.sbatch` | `/raven/ptmp/jinma/ccpush16/gates_a100.sbatch` | Raven 1-node 4xA100 gate job (build step + 1-GPU and 4-GPU runs, `--cpu-bind=none`) |
| `gates_a100_rerun.sbatch` | `/raven/ptmp/jinma/ccpush15/gates_a100_rerun.sbatch` | the map_cpu fix (single-GPU steps on GPUs 1-3) |
| `build_a100.sh` | `/raven/ptmp/jinma/ccpush17/build.sh` | A100 configure + make |

Production layout references (read only): `/viper/ptmp2/jinma/prod_1107/jobs/plink.sh`,
`/raven/ptmp/jinma/raven_gap_1105/jobs/gaplink.sh`.

## 6. Pitfalls

1. **Level-1 scratch** (open): Kokkos HIP reuses team scratch slots without a fence; kernels using level-1
   scratch can see stale data. Not fixed. Bitwise repeats on our gate decks hold; a repeat test on a real
   R136 production restart was started 11-16; its result is not in this note.
2. **Hydro sp barrier (R4)**: non-fofc spherical-polar decks on binaries without f6d3340e4 are non-deterministic
   on GPU. Use a binary that contains race3 (a9f481263) or later.
3. **HIP_MALLOC_ASYNC**: building with `Kokkos_ENABLE_IMPL_HIP_MALLOC_ASYNC=On` segfaults under ROCm 7.2.
4. **rocm/6.3 vs 7.2 module mismatch**: gpu72 binaries link `libamdhip64.so.7`; run under rocm/6.3 and HIP
   reports `hipErrorNoDevice`. Load the build modules in the job.
5. **nvcc and generic lambdas**: nvcc rejects a device (KOKKOS_LAMBDA) lambda inside a generic lambda (e.g. a
   runtime dispatch over `std::integral_constant`); hipcc accepts it. Code that builds on MI300A can fail
   to compile on A100.
6. **NVVM miscompile** (fixed in push 14): nvcc 12.6/12.8 miscompiled the clean cell face fluxes
   (runtime-indexed local tables + pointer selects jl/jr) in `rad_m1_clean_cell.hpp`; fix on
   fork/clean-cc-a100fix-viper (450722781, 4fe3efd4c, d6b6ed4e2). Watch for the same pattern in new kernels.
7. **cuda/12.8 has no openmpi_gpu** (your finding; consistent with ours, we only use cuda/12.6).
8. **viper login node memory cap 108 GB** (user cgroup, memory.max 108017184768): exceeding it kills every
   process of the user, including our sessions. Keep builds at -j <= 96 and analysis RSS small; do not copy
   tens of GB through the page cache (use `dd oflag=direct`).
9. **viper compute nodes cannot see /raven** (only the viper login node can): override any `/raven` path
   key in decks or restarts run on viper.
10. -O flags: we use the CMake Release default `-O3 -DNDEBUG` only; no other optimisation flags tested.
