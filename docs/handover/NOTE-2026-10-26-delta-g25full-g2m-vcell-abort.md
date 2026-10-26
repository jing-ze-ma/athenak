# NOTE 2026-10-26 Delta -> viper + deltaai (USER, URGENT): g25_full g2m on d9c0152bce ABORTS in VCELL at the first mesh build; DeltaAI: hold job3 until resolved

**What happened:** Delta job 23219344, the USER-approved run of the 43f11c08 ask.
- **Setup:**
  - 4 gpuA100x4 nodes, 256 ranks, 64 per node, 4 islands per node (16 ranks per island);
  - idort-gpu-nltree-2 d9c0152bce, DeltaAI's cfg2/g2m Config.sh;
  - built with gcc 13 + CUDA 12.9 at sm_80. CUDA 12.9 cannot parse gcc 14.2.1's headers.
- **Outcome:** g2m read the IC (6,092,886 gas + 1) and finished the density loop (26 iterations). Then it hit TERMINATE at
  `VCELL: create mesh (direct cell-by-cell construction)`, 250 s in.

**The aborts:** 8 ranks, all in `vcell_check_rcut_req()` (voronoi_cell.cc:136), all at **attempt 0**. Each cell's search radius was 0.50-0.52 x BoxSize.

| task | ID | pos (cm) | rcut / box |
|---|---|---|---|
| 82 | 6064612 | 1.92e15, 1.67e16, 1.07e16 | 0.512 |
| 93 | 6089337 | 1.56e16, 7.39e15, 1.72e16 | 0.517 |
| 90 | 6073648 | 6.84e15, 1.34e16, 1.72e16 | 0.517 |
| 94 | 6090331 | 1.61e16, 7.39e15, 8.21e14 | 0.501 |
| 83 | 6085818 | 1.40e16, 2.74e14, 2.74e14 | 0.517 |
| 88 | 6063825 | 1.92e15, 3.56e15, 2.74e14 | 0.517 |
| 86 | 6071218 | 5.75e15, 7.39e15, 8.21e14 | 0.512 |
| 92 | 6064177 | 1.92e15, 9.58e15, 2.74e14 | 0.517 |

- **Where the cells are:** all on the outer faces, edges and corners of the box (BoxSize 1.75155e16), and all among the highest IDs in the IC.
- **Periodicity:** `PeriodicBoundariesOn 0`. The half-box bound guards against second periodic images, so this box may not need it at all.
- **Against Raven:** the same IC ran on Raven at 288 ranks with e93504ccef. The check is older than that (since VORONOI_CELL, 08-15 to 09-26).
  Differences here:
  - the d9c0152bce IDORT commits;
  - 256 ranks instead of 288;
  - gcc 13.

**The CPU arms are running** in the same job. c1m and c0 results follow in the job NOTE.

**DeltaAI:** your job3 runs the same g2m on the same IC and will likely hit this. Hold it until viper says go.

**Next on Delta (USER asked to try a fix):** find why attempt 0's rcut exceeds half the box for these cells. Any change stays
local to Delta. No Arepo source or diffs in git; it will be described in words in the follow-up NOTE.
