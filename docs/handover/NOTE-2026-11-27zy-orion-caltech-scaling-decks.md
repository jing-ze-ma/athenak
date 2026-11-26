# NOTE 2026-11-27zy orion -> Caltech: addendum to 27zm, scaling decks that run at 1/4/8/16 GPUs (VETPAR)

From agent VETPAR (orion), for 27zm (VETPAR + RAYUNIFY baseline).

Why: the bsg production rst holds 4 MeshBlocks (256x128x128). A restart reads its own MeshBlock layout and the code has no re-partition (`BuildTreeFromRestart` stops when <mesh>/<meshblock> differ, `src/mesh/build_tree.cpp:484`), so 27zm could not run 8/16 GPUs. These arms are FRESH starts of the bsg template instead. Every arm has the same mesh and the same `layout_tile`. The tile-local preconditioner, the exact sums and the F-solve folds are therefore the same at every GPU count, so the Krylov counts and the CLEAN lines should be identical across the arms. Please check that: it is part of the result.

## Binary
Your REL `athena_rel6_01deec542` (= aa2646b92), as 27zm.

## Deck
The release template `inputs/stars/star_bsg.athinput` (aa2646b92). TABLES/ point to your copies. Outputs off, as in 27zm: set output2/3/5/6/7 `dt = 1e30`. Keep the log and hst; `time/ncycle_out = 1`.

Command-line keys for every arm:
`mesh/layout_tile=16 time/nlim=120 rad_vet/clean_verbose=1 time/ncycle_out=1`
plus the MeshBlock of the arm:

| arm | GPUs (nodes x GPUs) | ranks | `meshblock/nx1 nx2 nx3` | MeshBlocks | per GPU |
|---|---|---|---|---|---|
| S1 | 1 (1x1) | 1 | 256 256 256 | 1 | 1 |
| S4 | 4 (1x4) | 4 | 256 128 128 | 4 | 1 |
| S8 | 8 (2x4) | 8 | 256 128 64 | 8 | 1 |
| S16 | 16 (4x4) | 16 | 256 64 64 | 16 | 1 |

- One MeshBlock along x1 in all arms: tile mode requires it. Tile 16 divides every lateral MeshBlock and is the trunk's fresh default (e2a0f36c5).
- S1 may not fit in one H200's memory (16.8 M cells, all rays on one GPU). If it fails with out-of-memory, say so and skip it; S4 is then the base.
- Use the same node type for all arms and give each arm its nodes to itself (timing, 27zo).
- 120 cycles get past the start transient (the clean_init steady solve and the first ray builds). Measure over **cycles 61-120**.

## Report per arm (cycles 61-120, as 27zm asks)
1. s/cycle (elapsed lines) and zone-cycles/s per GPU.
2. Rays:
   - the `vet_gd` fenced seconds split: source, sweeps, of which per-shell halo+walls, of which MPI wait;
   - builds per cycle.
3. VET solve:
   - Krylov iterations per cycle and solves per stage (CLEAN `it K (sum T) ... solves S`);
   - the wall of the solve (clean stage steps).
4. Hydro + rest = total - rays - solve.
5. **Invariance check:** are the CLEAN lines (K, T, S, r1, r2, guard) identical across S1/S4/S8/S16, cycle by cycle? Yes, or the first cycle and field that differ.
6. One table and one stacked bar per GPU count. Use mesa2.mplstyle + idort_style.py.
7. If time allows: nsys rank 0 of S4 and S16, cycles 61-80. Report:
   - GPU idle and the number of gaps;
   - `cudaStreamSynchronize`, `cudaDeviceSynchronize`, H2D copies;
   - MPI_Allreduce / Allgather / Waitall counts;
   - Isend bytes per cycle.

Budget: about 1.5 node-hours. Delete dumps; no rst kept.
