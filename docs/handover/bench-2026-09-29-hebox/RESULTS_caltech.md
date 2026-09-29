# bench-2026-09-29-hebox: Caltech results (Resnick HPC, NVIDIA H200)

Commit 401875f0 (fork/rt-integration), `build_inc.sh hebox gpu 401875f0 box_convection` (CUDA HOPPER90,
Release, MPI), binary `builds/athena_hebox_gpu`, md5 572caa21782d215bbdf1f514a71ddcbe. Stack:
gcc/13.2 cuda/12.9 hpcx/2.17.1 ompi. Fresh start, hesdirk2, cfl 0.3, 1200 cycles, window cycles
400-1200 (ana_bench_hebox.py), 1 MPI rank per GPU (`srun --mpi=pmix`, one GPU visible per rank,
8 cores per rank), 2 repeats per job (r1 / r2). Jobs: 3615097 (1 GPU), 3615098 (2 GPUs), both on
hpc-sm-02-04. Logs: `/resnick/scratch/jingze/bench_hebox_0929/timing/n<N>_r<r>_<job>/bench.log`
(Caltech scratch, purged after 14 days). Every run rc 0, FATAL 0, nan 0, NON-CONVERGED 0, fallbacks 0.
Smoke: job 3615096 (2 GPUs), SMOKE_OK.

| machine | GPU | GPUs (ranks) | ms/cycle r1 / r2 | wall s / sim s r1 / r2 | dt | Picard/solve | Krylov it/linear solve | NC / fb | stack |
|---|---|---|---|---|---|---|---|---|---|
| Caltech | H200 | 1 (1) | 111.32 / 111.19 | 0.691 / 0.690 | 0.16111 | 6.476 | 6.55 | 0 / 0 | gcc/13.2 cuda/12.9 hpcx/2.17.1 |
| Caltech | H200 | 2 (2) | 105.77 / 105.75 | 0.657 / 0.656 | 0.16111 | 6.476 | 6.55 | 0 / 0 | gcc/13.2 cuda/12.9 hpcx/2.17.1 |
| Caltech | H200 | 4 (4) | not measured (job 3615099 cancelled by the user 09-29) | | | | | | |

- Physics identical to viper: same dt, Picard passes per solve (6.476), 42.44 BiCGStab iterations per
  transport solve and 6.55 Krylov iterations per linear solve at 1 and 2 ranks.
- **1 GPU: H200 = MI300A** (111.3 vs 111.8 ms/cycle). This box does not show the ~1.4-1.7x H200 lead
  of the WASP-121b benchmark (bench-2026-09-28).
- **2 GPUs scale poorly on H200: 1.05x** (viper 1.20x); 105.8 vs viper 93.2 ms/cycle (H200 13 % slower).
  The box is small (4 blocks of 84 x 52 x 52) and the implicit M1 solve is latency / global-reduction
  bound, so the extra rank mostly adds communication.
- The 4-GPU row was dropped: it would only extend a strong-scaling curve that is already flat, and
  commit 401875f0 predates the m1-perf merge (rt-integration 81c676c1, Newton opacity default, He box
  ~-60 %), so the absolute numbers are no longer current.
