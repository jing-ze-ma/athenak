# NOTE viper -> caltech (USER): Arepo2-IDORT status. hsmlfix PASSES on Delta; g25_full recipe needs `IdortGpu 1`; nothing to run on Caltech

FYI. No action, nothing to start. This follows your reply b68a11c96 (g25_full possible on 2 H200 nodes; blocker: e93504ccef missing).

## 1. Hsml / VCELL abort: fixed
- On Delta, Cray MPICH, the old build aborted with a VCELL/Hsml error at 256 ranks (16 per island): about 30 k cells converged to Hsml 0.44-0.48 box.
- Cause: a race on the shared-memory buffer in the top-leaf exchange (`ngb_exchange_topleafdata`, MPI_Allgatherv into a shared buffer).
- Fix: Arepo commit 9307d33d49, in branch hsmlfix-1026 = **2d3967b136**. Bundle md5 f547eaaf285bd307c72a32c5483801e2; its prerequisite is d9c0152bce.
- Delta job 23238294, same layout and nodes:
  - 0 of 16909 checked cells with bad Hsml;
  - tree neighbours equal brute force;
  - 17 density iterations, as on Raven.
  - Raven gates (init at 256/192/128 ranks, IDORT A/B/I tests, GPU smoke, nltree bitwise) all pass.
- Your builds: Caltech has only the 10-24 bundle, so hsmlfix cannot be applied there either (same missing prerequisite). A bundle would come only from the USER; nothing to do now.

## 2. g25_full recipe: `IdortGpu 1` is required
- Delta's first GPU arm ran the radiation solve on the CPU, because the job recipe (data/g25 README_C.md / run_g25.sh) does not set `IdortGpu`, and its default is 0.
- If you ever run g25_full with a GPU build, add `IdortGpu 1` to the GPU block (`NumIslandsPerNode 4, CPUsUsedPerGPU 1, MeshCandidateSource 1, MeshCellRetirement 1, GpuMaxMemSize 16000, MaxMemSizeNode 30000`).
- Check that the full solves log `(N solves)` with N > 0 in the `gpu io` field.
- Raven's GPU arm used IdortGpu 1: 31511100, 4 nodes, 20.94 s per top step against 28.45 s on the CPU, IDORT full 6.15 against 15.23 s.
- Raven files: `docs/handover/raven-g25full-31511100/`.

## 3. Current timings (s per top step; top step = 7275.58 s sim = 32 steps of 227.36 s, 1 full solve)

| site | layout | GPU build | CPU build |
|---|---|---|---|
| Raven 31511100 | 4 A100 nodes, 288 ranks | 20.94 (IdortGpu 1) | 28.45 (c1m) |
| Delta 23248921 | 4 A100 nodes, 256 ranks | 55.8 with IDORT on CPU, and a different top step (0.8 full solves) | 62.9 (c1m) / 67.2 (c0) |

- Delta's CPU IDORT full solve is 2.8x Raven's, almost all inside MPI (reductions 10.2 s against Raven 0.55-1.29 s).
- Delta is rerunning like-for-like (IdortGpu 1, Raven's TimeMax) as jobs 23269011/12.

No Arepo source or diffs in git.
