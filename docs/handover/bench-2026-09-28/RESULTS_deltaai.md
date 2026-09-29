# bench-2026-09-28: DeltaAI results (NCSA DeltaAI, NVIDIA GH200)

Commit 11c9a5be (`11c9a5be7d5cb8fc00a144f48f20897282df92d9`), Kokkos d8e9af03 (4.6.02). CUDA 12.9 (cudatoolkit/25.5_12.9),
Kokkos HOPPER90 + ARMV9_GRACE, nvcc_wrapper with host compiler Cray `CC` (PrgEnv-gnu, gcc-native/14),
cray-mpich/9.0.1 with GPU-aware MPI (libmpi_gtl_cuda), Release, MPI, `-D PROBLEM=deep_hot_jupiter_rt`.
Binary md5 **8e578218635a834a0c5449fb23ac6852** (`builds/bin/athena_dhj_gpu_11c9a5be7d5c`).
Built from a clean `git checkout --detach 11c9a5be` in a persistent worktree (`scripts/build_inc_deltaai.sh`,
status checked clean before make), not from `git archive`: same sources.

Fresh start, 2000 cycles, timing window cycles 1000-2000, median of 8-cycle windows (ana_bench.py). 1 MPI rank per
GPU, 1 node. Launch (the AthenaK docs' DeltaAI recipe): `MPICH_GPU_SUPPORT_ENABLED=1 SLURM_CPU_BIND=cores
KOKKOS_MAP_DEVICE_ID_BY=mpi_rank srun -n <n> -c 16 --cpu-bind=cores athena ...` with all GPUs of the job visible to
every rank. Rank r runs on GPU r and cores 72r..72r+15 (GPU r's own Grace, per `nvidia-smi topo -m`); checked
by sampling `nvidia-smi --query-compute-apps` during every job. Jobs 3256140 (1 GPU, gh046), 3256141 (2 GPUs,
gh010), 3256179 (4 GPUs, gh137), partition ghx4-interactive; 50-cycle smoke (ck_impl_verbose) plus 2 interleaved
repeats each; FATAL 0 and NOT-CONVERGED 0 in every run.

| machine | GPU | GPUs (ranks) | 1x ms/cycle | 1x wall s / sim s | 10x ms/cycle | 10x wall s / sim s |
|---|---|---|---|---|---|---|
| DeltaAI | GH200 | 1 (1) | 18.20 / 18.19 | 1.30e-3 / 1.30e-3 | 22.68 / 22.69 | 1.88e-3 / 1.88e-3 |
| DeltaAI | GH200 | 2 (2) | 11.76 / 11.66 | 0.84e-3 / 0.84e-3 | 14.85 / 14.88 | 1.26e-3 / 1.26e-3 |
| DeltaAI | GH200 | 4 (4, 1 node) | 8.10 / 8.13 | 0.58e-3 / 0.58e-3 | 10.53 / 10.53 | 0.88e-3 / 0.88e-3 |

Values are repeat 1 / repeat 2 (as in RESULTS_viper.md). Median (min-max) of the repeats, in the template's form:

| machine | GPU | ranks | 1x ms/cycle | 1x wall/sim-s | 10x ms/cycle | 10x wall/sim-s | commit | compiler / modules |
|---|---|---|---|---|---|---|---|---|
| DeltaAI | GH200 | 1 | 18.19 (18.19-18.20) | 1.30e-3 (1.30-1.30) | 22.68 (22.68-22.69) | 1.88e-3 (1.88-1.88) | 11c9a5be | PrgEnv-gnu gcc-native/14 cudatoolkit/25.5_12.9 cray-mpich/9.0.1 |
| DeltaAI | GH200 | 2 | 11.71 (11.66-11.76) | 0.84e-3 (0.84-0.84) | 14.87 (14.85-14.88) | 1.26e-3 (1.26-1.26) | 11c9a5be | same |
| DeltaAI | GH200 | 4 | 8.12 (8.10-8.13) | 0.58e-3 (0.58-0.58) | 10.53 (10.53-10.53) | 0.88e-3 (0.88-0.88) | 11c9a5be | same |

Scaling 1 -> 2 -> 4 GPUs (per cycle): 1x 1.55x / 2.24x, 10x 1.53x / 2.15x.
Window dt: 1x 13.969 s at every rank count (viper: 13.892 s); 10x 12.045 / 11.860 / 11.913 s at 1 / 2 / 4 ranks
(viper: 12.076 / 11.991 / 11.974 s). dt differs slightly from viper, so compare wall per simulated second.
One rotation (1.101535e5 s) at 4 GPUs: 1x ~64 s wall, 10x ~97 s wall (fresh-start dt).

## Against viper (MI300A, RESULTS_viper.md), ratio viper / DeltaAI (>1 = DeltaAI faster)

| GPUs | 1x ms/cycle | 1x wall/sim-s | 10x ms/cycle | 10x wall/sim-s |
|---|---|---|---|---|
| 1 | 1.43 | 1.42 | 1.50 | 1.50 |
| 2 | 1.64 | 1.64 | 1.69 | 1.66 |
| 4 (viper: 2 nodes x 2) | 1.57 | 1.58 | 1.58 | 1.57 |

## Logs

- Job outputs: `/u/jma20/ATHENAK/runs/bench_w121bench.{3256140,3256141,3256179}.out`
- Run dirs: `/work/nvme/bivj/jma20/bench_0928/{smoke_,}n<N>_<arm>[_r<k>]_<job>/bench.log` (+ run.athinput,
  dhj.hydro.hst); GPU samples `gpus_n<N>_<job>.csv`.
- Job script: `bench_deltaai.sub` (this directory). `run_bench.sh` here has the `deltaai)` header filled
  (modules, env, launcher, data dirs); the copy on rt-integration is unchanged.

## Notes

- Do not launch with `--gpus-per-task=1 --gpu-bind=closest` on DeltaAI: Slurm's GPU cgroups block CUDA IPC, and
  GPU-aware Cray MPICH then hangs at the first halo exchange (or, with MPICH_GPU_IPC_ENABLED=0, runs 2.8x slower
  at 2 GPUs and 3.3x slower at 4). Tested but no gain at 4 GPUs: `--constraint=xpmem`,
  MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=0, MPICH_GPU_IPC_THRESHOLD=1 (all within +-0.3 %);
  MPICH_GPU_MAX_NUM_STREAMS=4 segfaults. Details: DELTAAI_FACTS.md on DeltaAI.
- A first attempt (jobs 3256102/03/07) stopped at the smoke: `ckdata10/{cia,ray,sw_flux}` were symlinks to
  viper paths. Repointed to the local copies of the same repo tables (`exo_fms_ck/{cia,ray,sw_flux}`); no timed
  run used the broken data.
