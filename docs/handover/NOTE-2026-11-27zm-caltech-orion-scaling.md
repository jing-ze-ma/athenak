
# NOTE 2026-11-27zm caltech -> orion: strong-scaling breakdown of a bsg cycle on H200 (your 27zm + 27zy)

- Binary: rel6 01deec542. Every arm had its nodes to itself. Cadence key time/ndiag=1 (time/ncycle_out is not a key in this code: on a restart it is ignored, on a fresh start it FATALs as an unread command-line key). Outputs off; dumps deleted.

## Series 1: fresh ladder per 27zy
- Setup: bsg release template, layout_tile 16, nlim 120, cycles 61-120 measured.
- Layouts: S1 = 1 MeshBlock 256^3; S4 = 4 x 256x128x128 (02-11); S8 = 8 x 256x128x64 (02-11 + 02-17, 7333073).
- S16 is impossible: it needs 4 nodes in one allocation, and our largest holder has 3. S12 is not a legal split (256 x 256 does not cut into 12 equal tile-16 blocks). So the ladder stops at S8. S1 fits on one H200 (no OOM).
| arm | GPUs | s/cycle (61-120) | zone-cycles/s | rays: source / sweeps (halo+walls, MPI wait) s/cycle | builds/cycle | cut shell | VET solve: clean stage steps s/cycle | Krylov it/cycle | stage - rays | total - stage (hydro + rest) |
|---|---|---|---|---|---|---|---|---|---|---|
| S1 | 1 | 2.8006 | 5.991e+06 | 0.0031 / 0.1867 (0.0000, 0.0000) | 1.1 | 148 | 2.1367 | 41.6 | 1.9469 | 0.6638 |
| S4 | 4 | 1.2637 | 1.328e+07 | 0.0204 / 0.0622 (0.0025, 0.0169) | 1.1 | 148 | 0.7308 | 41.6 | 0.6482 | 0.5329 |
| S8 | 8 | 0.8808 | 1.905e+07 | 0.0257 / 0.1316 (0.1306, 0.0988) | 1.1 | 148 | 0.5925 | 41.6 | 0.4352 | 0.2883 |
- Column notes:
  - s/cycle is over cycles 61-120 (elapsed lines).
  - The fenced ray seconds, the clean stage steps and the Krylov iterations are run totals / 120 cycles, start-up included; the code prints them once at the end. So "total - stage" mixes the two windows.
  - zone-cycles/s is for the whole mesh, not per GPU. Per GPU: 6.0e6 / 3.3e6 / 2.4e6.
- **Invariance (item 5): yes.** The CLEAN lines (it K, sum T, r1, r2, guard, solves S) of S1, S4 and S8 are identical cycle by cycle: 242/242 lines in every pair. The Krylov iterations are 41.6 per cycle in all three.
- Reading:
  - S1 -> S4 -> S8: 2.80 -> 1.26 -> 0.88 s/cycle (speed-up 2.2x at 4 GPUs, 3.2x at 8).
  - The VET solve (stage minus rays) scales 1.95 -> 0.65 -> 0.44 s.
  - The ray sweep grows from 4 to 8 GPUs: 0.062 -> 0.132 s per cycle. That is almost all per-shell halo+walls (0.003 -> 0.131 s; MPI wait 0.017 -> 0.099 s) once the band crosses the node boundary.
  - At 8 GPUs the rays are 18 % of the cycle against 7 % at 4. They are not yet dominant, but they grow while everything else shrinks.

## Series 2: production rst (bsg 00009, 4 MeshBlocks 256x128x128), cycles 5-20 of 20, ndiag 1
- A separate series from series 1. The rst cannot restart with a different MeshBlock size (FATAL `meshblock/nx2: restart file 128 -> now 64`), so only 1 and 4 GPUs.
| arm | GPUs | s/cycle (5-20) | zone-cycles/s | rays: source / sweeps (halo+walls, MPI wait) s/cycle | builds/cycle | cut shell | VET solve: clean stage steps s/cycle | Krylov it/cycle | stage - rays | total - stage (hydro + rest) |
|---|---|---|---|---|---|---|---|---|---|---|
| G1 | 1 | 2.3753 | 7.063e+06 | 0.0045 / 0.2453 (0.0000, 0.0000) | 1.0 | n/a | 2.2307 | 33.9 | 1.9809 | 0.1446 |
| G4 | 4 | 0.9228 | 1.818e+07 | 0.0270 / 0.0817 (0.0028, 0.0244) | 1.0 | n/a | 0.7903 | 33.9 | 0.6816 | 0.1325 |
- Developed convection: 33.9 Krylov iterations per cycle. The ray share at 4 GPUs is 0.11 of 0.92 s/cycle.

- Plots (mesa2 + idort_style): docs/handover/caltech_results/SCALING/scaling_S.png (series 1), docs/handover/caltech_results/SCALING/scaling_G.png (series 2). The relevant log lines of each arm are in docs/handover/caltech_results/SCALING/<arm>_lines.txt.
- Not used: my earlier fresh ladder (tile 64, cycles 5-20, F1/F4/F8) is superseded by series 1. Its F8 also overlapped another step of mine once.
- nsys of S4/S16 (item 7): not run. S16 is impossible here, and the remaining time went to the queue.

## Item 7: nsys of S4, rank 0
- Capture: delay 79 s, duration 15 s, set from S4's elapsed lines.
- By the run's own elapsed lines that window holds about 17 cycles (68-84). The per-cycle values below divide by 17; approximate.
- nsys ended the run at the end of the capture (so no Terminating line).
cycles in the capture window (elapsed 79-94 s): 17 (cycles 68-84)
| per cycle (rank 0, nsys capture window) | S4 |
|---|---|
| GPU timeline span [s] | 0.8805 |
| GPU idle [s] | 0.3685 |
| idle gaps > 10 us | 1491.2 |
| cudaStreamSynchronize | 776.0 |
| cudaDeviceSynchronize | 247.2 |
| cudaEventSynchronize | 985.4 |
| cudaMemcpyAsync calls | 316.3 |
| H2D copies | 1159.6 |
| MPI_Allreduce | 125.0 |
| MPI_Allgather | 144.8 |
| MPI_Waitall | 1593.2 |
| MPI_Isend bytes | 2864892265.4 |
- S16 is impossible here (largest holder 3 nodes). The VETPAR 27zt nsys pair on the rst (R4n) gives the same picture: about 0.4-0.5 s of GPU idle per cycle on rank 0 at 4 GPUs.
