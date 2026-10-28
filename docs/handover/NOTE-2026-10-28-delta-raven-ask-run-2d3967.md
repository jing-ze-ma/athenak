# NOTE 2026-10-28 Delta -> viper (Raven) (USER): please run g25_full g2m + c1m with hsmlfix 2d3967b136 on Raven. Delta and Raven ran different code, so the 2x gap in 69b9b4f4 may be the code, not Delta

Follows 69b9b4f4. No Arepo source or diffs in git.

## Why
- **The code differs.** Delta's like-for-like run 23269096 used hsmlfix **2d3967b136**. Raven's 31511100 used **e93504ccef**. 2d3967b136 is 21 commits ahead of e935 (e935 is an ancestor).
- **Two of those commits change synchronisation that runs on every step:**
  - **c59fd6214**: every barrier on a shared-memory window becomes fence + MPI_Barrier + fence. This covers NodeBarrier, ngbnl_barrier, the VCELL wave island barriers, the IDORT staging-window barriers and parallel_sort.
  - **9307d33d4**: the top-leaf exchange (ngb tree, force tree, SIDM) now gathers privately, and the island master publishes the result.
  - The other 19 are IDORT speedups (marked "exact") and diagnostics.
- **The slow sections match what those commits touch.** In g2m on Delta vs Raven (s / top step): voronoi 12.0 vs 5.1, crhalo 5.1 vs 1.8, cellrounds 6.2 vs 2.5, treegrav 5.3 vs 2.9. In c1m, reductions take 11.6 vs 0.8 s per full solve.
- **Everything else is the same** or explains only ~10 %:
  - compiler and flags: gcc 13, `-march=native -O3 -flto`;
  - parameters, step size and Config: only Raven's OUTPUT_ACCELERATION differs;
  - per-core speed: HEALTHTEST 0.74 s on Delta vs 0.80 s on Raven;
  - pair bandwidth: 920 vs 863 MB/s;
  - GPU sweep: 2.1 vs 2.3 s;
  - rank count: 256 vs 288.
- **Delta can't run e935 to check from its side.** e935 aborts at init on Delta at 256 ranks (the wall-cell Hsml bug; V2 in cdebc9cf), and Delta's 64-core nodes can't take Raven's 72 per node.

## Ask (one job, same nodes, so node and day variation cancel)
1. Build **g2m and c1m from 2d3967b136** (hsmlfix-1026; the bundle you sent, arepo_hsmlfix_1026.bundle), with the same Config.sh and Makefile settings as 31511100.
2. Run them with the 31511100 recipe: run_full.sh, 4 nodes, 288 ranks (72/node, `--cpu-bind=cores`), OpenMPI 5. Use g2m IdortGpu 1 + GpuMaxMemSize 16000, TimeMax 29102.3 (4 top steps), and the same env (AREPO_IDORT_PHASES=1 is fine).
3. **If the budget allows, in the same job:** **e93504ccef g2m + c1m** as a control on the same nodes (the 31511100 builds).
4. Report, per arm:
   - s / top step (your definition: mean of top steps 2-4);
   - the section breakdown (IDORT full, mesh, gravity, IDORT partial, rest);
   - and put `cpu.txt`, `idort_solver.txt` and the `IDORT: timing|HEALTHTEST` lines under docs/handover/raven-g25full-<job>/.

## How we read it
- **2d3967 on Raven ≈ e935 on Raven (~21 / ~28 s):** the code is not the cause. The gap is Delta's platform (Cray MPICH / Slingshot, or the 64-core single-socket node with c1m's single 64-rank island). Delta's test 1 (below) and an allreduce/barrier latency microbenchmark then split it.
- **2d3967 on Raven clearly slower than e935:** the regression is in the code, most likely the fenced barriers (c59fd6214) or the top-leaf exchange (9307d33d4). A Raven run with c59fd6214 reverted would confirm it.

## Delta meanwhile
Test 1 is queued: **23293073** (interactive) + **23293074** (regular), queue-both.
- It runs c1m with 2d3967 at 48 / 64 / 56 ranks per node on the same 4 nodes.
- The idle cores are spread evenly over the four 16-core groups.
- It shows whether leaving cores idle shrinks Delta's reduction, halo and orthogonalisation time.
