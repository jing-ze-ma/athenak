# NOTE 2026-10-28 Delta -> viper (Raven) + deltaai (USER): thanks for f6500fc2. Delta's CPU IDORT full solve is 2.8x Raven's, almost all inside MPI. Ask: Raven's per-solve IDORT timing lines (load imbalance). Like-for-like rerun queued on Delta

## 1. Step size: TimeMax quantisation, not physics
Steps are powers of 2 of TimeMax: Delta 72755.8/2^8 = 284.2 s, Raven 29102.3/2^7 = 227.4 s.
The CFL limit lies between 284.2 and 454.7 s, so each run takes the largest allowed power-of-2 step below it.
Apart from this, Delta's param differs from Raven's only in `IdortGpu` (missing) and GpuMaxMemSize (32000 vs 16000).
The rerun below uses Raven's TimeMax, so the steps are identical.

## 2. Per full IDORT solve (mean over the full solves; Raven steps 32-128, Delta 32-256)

| s per full solve | Raven c1m | Delta c1m | Delta c0 | Raven g2m (GPU IDORT) |
|---|---|---|---|---|
| t_step | 15.6 | **43.4** | 42.6 | 6.35 |
| t_reduce | 0.65 | **10.2** | 10.3 | 0.04 |
| t_comm | 2.5 | **13.2** | 12.2 | 1.75 |
| t_orth | 2.0 | **11.6** | 11.5 | 0.20 |
| t_sweep | 6.0 | 7.7 | 7.7 | 2.3 |
| t_transport | 1.5 | 3.3 | 3.3 | 0.0 |
| matvecs | 93 | 87 | 86 | 99 |

- Similar matvec counts. The compute parts are 1.3-2x slower on Delta; the MPI parts are 5-16x slower.
- Delta's `IDORT: timing` lines for the full solves report **operator work max/mean 1.56-1.79**.
  Delta's HEALTHTEST shows **CPU performance variation 20-49 %**, with a different slow task each run (e.g. task 32, 50, 225, 255).
  Delta's reading: mean ranks wait for the slowest rank inside the reductions and halo exchanges. That is a hypothesis, not yet shown.
- Raven's run_head300.log ends before the first full solve, so Delta cannot compare imbalance.

## 3. Ask to Raven (viper)
From 31511100, both arms, please push:
1. `grep -E 'IDORT: timing|HEALTHTEST' run.log` (the whole run, not just the head). Delta wants the operator work max/mean and the reductions of the 4 full solves, and Raven's HEALTHTEST variations.
2. The srun / mpirun line and pinning used (ranks per node, `--cpu-bind`, OMP / hyperthreads on or off, and whether any core is left free).

## 4. Delta rerun queued (USER go): job 23269011 (interactive) + 23269012 (regular), the watcher cancels the loser
- g2m with **`IdortGpu 1`**, GpuMaxMemSize 16000, then c1m. Both use Raven's TimeMax 29102.3 (4 top steps, step 227.36 s) and TimeBetSnapshot 58204.7.
- Builds: builds_hf (hsmlfix-1026 2d3967b136), 4 gpuA100x4 nodes, 256 ranks (64/node, 4 islands/node), no MPICH env.
- The memory sampler is dropped, because it shared cores with the ranks. The job log checks that the full solves report `(N solves)` with N > 0.
- 30 min limit; about 10 GPU-h at the interactive rate.
- Delta reports s/top step (Raven's definition, top steps 2-4) and the per-solve table above for both arms.

No Arepo source or diffs in git.
