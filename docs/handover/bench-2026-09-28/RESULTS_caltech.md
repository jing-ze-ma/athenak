# bench-2026-09-28: Caltech results (Resnick HPC, NVIDIA H200)

Commit 11c9a5be (git archive build via build_inc.sh, CUDA HOPPER90, gcc/13.2 cuda/12.9 hpcx/2.17.1 ompi,
Release, MPI), binary `builds/athena_bench_gpu` md5 d52495867f7a2eac5932676e5b99f37b. Fresh start,
2000 cycles, timing window cycles 1000-2000, median of 8-cycle windows (ana_bench.py). 1 MPI rank per GPU,
`srun --mpi=pmix`, one GPU visible per rank (`--gpus-per-task=1`), 8 cores per rank. Jobs 3603571 (1 GPU),
3603572 (2 GPUs), both on hpc-sm-01-11; 2 interleaved repeats each; FATAL 0 and NOT-CONVERGED 0 in every run.
Logs: `/resnick/scratch/jingze/bench_0928/runs` (Caltech scratch, purged after 14 days).

| machine | GPU | GPUs (ranks) | 1x ms/cycle | 1x wall s / sim s | 10x ms/cycle | 10x wall s / sim s |
|---|---|---|---|---|---|---|
| Caltech | H200 | 1 (1) | 17.99 / 17.98 | 1.29e-3 / 1.29e-3 | 22.51 / 22.51 | 1.88e-3 / 1.88e-3 |
| Caltech | H200 | 2 (2) | 11.33 / 11.33 | 0.81e-3 / 0.81e-3 | 14.38 / 14.38 | 1.19e-3 / 1.19e-3 |
| Caltech | H200 | 4 (4, 1 node) | ~8.45 (est., see below) | ~0.61e-3 (est.) | not measured | not measured |

Scaling 1 -> 2 GPUs: 1x 1.59x, 10x 1.57x (per cycle). Window dt: 1x 13.905 s, 10x 12.010 / 11.992 s at
1 / 2 ranks. H200 vs viper MI300A at equal GPU count: 1x 1.44x / 1.70x, 10x 1.52x / 1.74x faster (1 / 2 GPUs).

## 4 GPUs: launch probe only (full 4-GPU job 3603573 cancelled by the user 09-29)

Job 3603625 (hpc-sm-02-04, 1x arm only, 600 cycles, 2 interleaved repeats, window cycles 300-600):

| variant | launch | 1 GPU ms/cycle | 4 GPUs ms/cycle |
|---|---|---|---|
| A | as the table above (one GPU visible per rank) | 17.58 / 17.59 | 8.25 / 8.27 |
| B | all 4 GPUs visible, Kokkos maps device by local rank | - | 7.71 / 7.72 |
| C | B + NUMA pinning of each rank to its GPU's socket | 17.61 / 17.61 | 7.72 / 7.71 |

The 300-600 window is still in the start-up transient (dt 15.33 s, and the same 1-GPU run gives 17.60
there vs 17.99 in 1000-2000), so the 4-GPU row above is scaled by 17.99 / 17.58: **~8.45 ms/cycle (A)**,
~7.9 (B), i.e. 1 -> 4 GPUs 2.13x (A) / 2.28x (B). All GPUs visible per rank (B) is -6.7 %; NUMA pinning adds
nothing. The 10x arm was not run on 4 GPUs.
