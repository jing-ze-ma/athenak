# NOTE 2026-10-28 Delta -> viper (Raven) + deltaai (USER): like-for-like g25_full rerun 23269096 done. g2m 39.9 / c1m 73.4 s per top step (Raven 20.9 / 28.5). Imbalance explains only ~1/5 of Delta's extra MPI time: correction to a7c528ba

Follows a7c528ba and Raven's 7b1496db / f6500fc2. No Arepo source or diffs in git.

## Run
- **23269096** (interactive; the regular copy 23269089 was cancelled while PENDING). 4 x gpuA100x4, 256 ranks (64/node, `srun -c 1 --cpu-bind=cores`), Cray MPICH.
- Code: hsmlfix 2d3967b136 (builds_hf). Raven's TimeMax 29102.3 (4 top steps). `AREPO_IDORT_PHASES=1`.
- g2m: IdortGpu 1, GpuMaxMemSize 16000; then c1m. Both rc 0. Elapsed 12:47, 6.82 GPU-h.
- **g2m ran IDORT on the GPU this time:** 5 full solves (1 init + 4 top steps). The GPU memory high-water mark was 1.5 GB per task.

## 1. s / top step (Raven definition: cpu.txt total over 32 steps, mean of top steps 2-4)

| arm | Delta 23269096 | top steps 1-4 | Raven 31511100 | Delta / Raven |
|---|---|---|---|---|
| g2m (GPU IDORT) | **39.86** | 49.6 / 45.8 / 35.2 / 38.5 | 20.94 | 1.9x |
| c1m (CPU node-local) | **73.39** | 99.9 / 72.3 / 73.4 / 74.5 | 28.45 | 2.6x |

g2m cpu.txt sections, s / top step (top steps 2-4), Delta vs Raven:

| section | Delta | Raven |
|---|---|---|
| misc (IDORT inside) | 15.4 | 10.3 |
| voronoi | 12.0 | 5.1 |
| cellrounds | 6.2 | 2.5 |
| treegrav | 5.3 | 2.9 |
| crhalo | 5.1 | 1.8 |
| i/o | 3.7 | < 0.8 |

So the gap is **not only IDORT**: mesh, gravity and halo work are each about 2x Raven too.

## 2. Per full IDORT solve (mean over the full solves, step > 5 s; from idort_solver.txt)

| arm | t_step | matvecs | comm | orth | reduce | sweep | transport |
|---|---|---|---|---|---|---|---|
| Delta g2m | 10.12 | 95 | 3.26 | 0.63 | 0.47 | 2.10 | - |
| Raven g2m | 6.35 | 99 | 1.75 | 0.20 | 0.04 | 2.26 | - |
| Delta c1m | **45.71** | 92 | 13.53 | 12.83 | 11.63 | 7.73 | 3.39 |
| Raven c1m | 15.39 | 93 | 2.59 | 2.12 | 0.78 | 5.70 | 1.51 |

- The iteration counts are equal, so the difference is cost per iteration, not convergence.
- On the GPU, the sweep costs the same on both machines. The g2m gap is comm (+1.5 s) and reduce/orth (+0.4 s), and it is mostly outside IDORT (§1).
- c1m: the CPU compute (sweep + transport + lin + block) is 12.3 s on Delta vs 8.2 s on Raven (1.5x). halo + orth is 26.4 s vs 4.7 s (5.6x).

## 3. Imbalance and AREPO_IDORT_PHASES: correction to a7c528ba
- Full-solve operator work max/mean:
  - **g2m**: Delta 1.11-1.26, Raven 1.11-1.20. The same.
  - **c1m**: Delta 1.41-1.56 (this run), Raven 1.19-1.26.
- **a7c528ba overestimated the waiting.** It applied 1 - 1/imbalance to the *whole* solve time, which already includes the MPI time.
  - The waiting comes from the spread of the compute, so the right estimate is (max/mean - 1) x the mean compute time.
  - c1m: Delta ≈ **6.3 s**, Raven ≈ 1.8 s, per full solve.
  - Delta's extra halo + orth time is ~21.7 s, of which imbalance explains **~4.5 s (~1/5)**, not 2/3.
- g2m shows the same thing more clearly: equal imbalance, but comm is 1.9x and reductions are ~10x Raven.
- **The PHASES lines cannot decide it.** Every full-solve line has `solve` max/mean = 1.00 (for example g2m 12819/12783 ms, c1m 51131/51084 ms).
  - The solve ends in collectives, so every rank leaves at the same time, and the wait is hidden inside the phase.
  - The phases with a spread (fill, xinit, store; max/mean 1.3-2) are all < 0.3 s.
- **What remains is a per-call MPI cost on Delta, and it grows when the CPU is busy:**
  - Reductions cost 11.6 s per solve in c1m but 0.47 s in g2m, with the same 256 ranks, nodes, network and solver iterations.
  - The difference is that in c1m all 64 cores per node compute between collectives; in g2m the cores mostly wait on the GPU.
  - Point-to-point bandwidth is not the problem: HEALTHTEST hypercube is 908-934 MB/s per pair on Delta vs 863 MB/s on Raven.
  - The likely causes are per-iteration jitter on fully subscribed cores (HEALTHTEST CPU variation 19 % in g2m and 47 % in c1m in this run; the slowest task differs between the two arms) and/or Cray MPICH collective or progress cost. One run cannot separate them.

## 4. Possible next tests on Delta (not submitted; the USER decides)
1. c1m with fewer ranks per node (e.g. 56 or 48, leaving cores free for the OS and MPI progress), same 4 nodes. If reduce/orth drops by much more than the compute rises, it is jitter/progress, not the code.
2. The same c1m with `MPICH_ASYNC_PROGRESS=1`, or with a Delta OpenMPI build if one exists, to separate the MPI library from the hardware.
3. For g2m, the 2x gap outside IDORT (voronoi, crhalo, cellrounds) matters more than the IDORT gap. Test 1 also tests it.

Files: Delta `/work/nvme/bivj/jma20/delta_1008/arepo2/runs/23269096/{g2m,c1m}` (run.log, output/cpu.txt, output/idort_solver.txt), log `logs/a2g25rv.23269096.out`.
