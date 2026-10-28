# NOTE viper (Raven) -> Delta + deltaai (USER): answer to 868d9520: Raven 31511100 IDORT timing and HEALTHTEST lines, and the launch/pinning

New files under `docs/handover/raven-g25full-31511100/`:
- `{g2m,c1m}/run_idort_timing_healthtest.log`: the whole run, `grep -E 'IDORT: timing|HEALTHTEST'`; 136 and 134 lines.
- `run_full.sh.txt`: the exact job script for both arms.

## 1. The 5 full IDORT solves (all-active steps 0, 32, 64, 96, 128), from those lines

| arm | step t (max), s | reductions, s | operator work max/mean |
|---|---|---|---|
| c1m | 14.43 / 16.83 / 17.25 / 14.39 / 14.05 | 1.29 / 0.77 / 0.73 / 0.57 / 0.55 | **1.25 / 1.20 / 1.26 / 1.20 / 1.19** |
| g2m (GPU IDORT) | 4.79 / 6.96 / 6.87 / 6.10 / 5.48 | 0.04 / 0.05 / 0.04 / 0.03 / 0.03 | 1.15 / 1.20 / 1.11 / 1.14 / 1.19 |

- Raven's full-solve operator imbalance is 1.19-1.26, against Delta's 1.56-1.79.
- The partial (few-cell) solves show max/mean around 40-50 on Raven too. That is expected, since few ranks hold active cells.

## 2. HEALTHTEST (Raven, first lines of each arm)
- CPU performance variation: g2m 20.8 %, worst task 62/node 0. c1m 120 %, worst task 201/node 2; the second test is 22 %, worst task 161.
  So Raven also has noisy single tasks, varying from run to run.
- Full hypercube: about 860 MB/s per pair, 94 % variation. Internode cube: about 40 GB/s, 2-3 % variation.
- Intranode cube (c1m, node 0): 2.0 GB/s per pair, 109 % variation.

## 3. Launch and pinning (both arms, the same allocation)
- `#SBATCH --ntasks-per-node=72 --cpus-per-task=1 --gres=gpu:a100:4 --mem=480000`, 4 nodes = 288 ranks.
- `srun --cpu-bind=cores bash ./mon.sh`. `mon.sh` execs `./Arepo param.txt`; local rank 0 also forks a `free`/`nvidia-smi` sampler every 5 s, which shares rank 0's core.
- **No core is left free**: 72 ranks on 72 physical cores per node (2 x 36, Xeon Platinum 8360Y). SMT is present (2 threads/core, 144 logical), but with cpus-per-task 1 + cpu-bind=cores, each rank gets one core.
- `OMP_NUM_THREADS=1`. Modules: gcc/13, cuda/12.6, **openmpi_gpu/5.0**, hdf5-serial/1.14.1, hwloc/2.12. No MPI env overrides.
- Env: `AREPO_NGBNL_GRAPH_POOL=160`, `AREPO_IDORT_PHASES=1`.
- GPU arms get `NumIslandsPerNode 4, CPUsUsedPerGPU 1, MeshCandidateSource 1, MeshCellRetirement 1, GpuMaxMemSize 16000, MaxMemSizeNode 30000`, plus `IdortGpu 1` when the build has IDORT_GPU. CPU arms get `MaxMemSizeNode 16000`.

## 4. Remarks (facts only; for Delta to weigh)
- Delta runs 64 ranks per node (16 per island) and Raven 72 per node (18 per island). Raven uses OpenMPI 5, Delta Cray MPICH.
- Raven's reductions per full solve (0.55-1.29 s) are 8-18x smaller than Delta's 10.2 s, at similar operator imbalance. So imbalance alone may not explain Delta's t_reduce. The MPI library or interconnect path is the other candidate.

No Arepo source or diffs in this push.
