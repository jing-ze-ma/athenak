# bench-2026-09-29-hebox: DeltaAI results (NCSA DeltaAI, NVIDIA GH200 120GB)

Commit 401875f0 (fork/rt-integration), `scripts/build_inc_deltaai.sh box_gpu 401875f0` (PROBLEM=box_convection,
Release, MPI, CUDA, `Kokkos_ARCH_HOPPER90`, `Kokkos_ARCH_ARMV9_GRACE`, nvcc_wrapper with CC), binary
`/u/jma20/ATHENAK/builds/bin/athena_box_gpu_401875f01beb`, md5 c646cb8f433543133ea885c3d9a71249.
Stack: **default modules: PrgEnv-gnu, gcc-native/14, cudatoolkit/25.5_12.9 (nvcc 12.9.41), cray-mpich/9.0.1,
craype-accel-nvidia90**, `MPICH_GPU_SUPPORT_ENABLED=1`. Launch = the AthenaK docs' DeltaAI recipe (every rank sees
all GPUs, `KOKKOS_MAP_DEVICE_ID_BY=mpi_rank`, `srun -n N -c 16 --cpu-bind=cores`); the `deltaai)` header of
`run_bench_hebox.sh` on this branch. Package (input, ana script) unchanged from rt-integration b3a3f53f.

**One job, 3258635** (ghx4-interactive, 1 node gh053, 4 GPUs held for the whole job; 13 min 43 s, ~1.8 SU):
smoke first, then 2 repeats of 1, 2, 4 GPUs interleaved (n1 n2 n4 n1 n2 n4), 1 MPI rank per GPU; a 1- or 2-rank run
uses GPUs 0..n-1 of the same node. Script `bench_deltaai_hebox.sub` (this directory).
Logs: `/work/nvme/bivj/jma20/bench_hebox_0929/n<N>_r<r>_3258635/bench.log`, job output
`/u/jma20/ATHENAK/runs/hebox_heboxbench.3258635.out`. Every run rc 0, FATAL 0, nan 0, NON-CONVERGED 0, fallbacks 0.
Fresh start, hesdirk2, cfl 0.3, 1200 cycles, window cycles 400-1200 (`ana_bench_hebox.py` default).

| machine | GPU | GPUs (ranks) | ms/cycle r1 / r2 | wall s / sim s r1 / r2 | dt | Picard/solve | Krylov it/linear solve | NC / fb | stack |
|---|---|---|---|---|---|---|---|---|---|
| DeltaAI | GH200 | 1 (1) | 114.77 / 114.84 | 0.712 / 0.713 | 0.16111 | 6.476 | 6.55 | 0 / 0 | PrgEnv-gnu gcc-native/14 cudatoolkit/25.5_12.9 cray-mpich/9.0.1 |
| DeltaAI | GH200 | 2 (2) | 96.96 / 96.91 | 0.602 / 0.602 | 0.16111 | 6.476 | 6.55 | 0 / 0 | same |
| DeltaAI | GH200 | 4 (4, 1 node) | 81.37 / 81.68 | 0.505 / 0.507 | 0.16111 | 6.476 | 6.55 | 0 / 0 | same |

- ms/cycle = median of the 10-cycle windows; the whole-window means are within 0.3 % of it.
- dt, Picard passes per solve (6.476), BiCGStab iterations (42.44 per transport solve, 6.55 per linear solve) are
  **identical to viper** (RESULTS_viper.md) at all three counts.
- Window check: cycles 100-1200 give 115.97 / 97.96 / 82.24 ms (r1), within 1.1 % of 400-1200 (viper: also 1.1 %).
- Repeats agree to 0.4 % or better.
- Scaling per cycle 1 -> 2 -> 4 GPUs: 1.18x / 1.41x (viper 1.20x / 1.47x).
- **Versus viper MI300A (ROCm 7.2), r1:** GH200 is 2.7 % slower at 1 GPU (114.77 vs 111.78), 4.1 % at 2 (96.96 vs
  93.17), 8.1 % at 4 (81.37 vs 75.26; viper's 4 GPUs are 2 nodes x 2). Unlike WASP-121b (bench-2026-09-28, GH200
  1.4-1.7x faster than MI300A), the GH200 has no edge on this implicit-M1 box. Not investigated (timing only).
- Wall per run incl. start-up (~6 s to cycle 0): 2.5 min at 1 GPU, 2.1 min at 2, 1.8-1.9 min at 4.
- Smoke (same job, 2 GPUs, 200 cycles, ndiag 1, Picard log): clean; 104.02 ms/cycle (cycles 50-200; viper 102.1),
  Picard 6.759 per solve (max 10; viper 6.76), 6.75 BiCGStab iterations per linear solve, dt 0.16125.
