---
name: solar-convection-scaling-wall
description: "MEASURED: solar_convection saturates at ~36 MPI ranks; multinode buys nothing, and the RT column sweep is why"
metadata:
  type: project
---

Measured 2026-09-03 on orion, 400 cycles, 96x192x192, identical input and binary.
Asked because the user wanted the movie run on multiple nodes.

| ranks | layout | zone-cycles/s | speedup vs 9 |
|---:|---|---:|---:|
| 9 | 1 node | 4.63e6 | 1.00 |
| 16 | 1 node | 7.44e6 | 1.61 |
| 36 | 1 node | 1.055e7 | 2.28 |
| 48 | 1 node | 1.014e7 | 2.19 |
| 72 | 1 node | 9.98e6 | 2.15 |
| 48 | 2 nodes (24+24) | 9.57e6 | 2.07 |
| 72 | 2 nodes (36+36) | 1.040e7 | 2.25 |
| 144 | 2 nodes (72+72) | 1.117e7 | 2.41 |

**Why:** the mesh is only 3.5M cells and `nx1` CANNOT be decomposed -- the two-stream RT
sweeps whole vertical columns, so `<meshblock> nx1` must equal `<mesh> nx1`. That caps the
decomposition at 12x12 = **144 blocks**, and a 96x16x16 block already carries 13312 ghost
cells against 24576 real ones. Past ~36 ranks the extra cores go into ghost exchange.
Spreading the SAME 48 ranks over two nodes is slightly SLOWER, which rules out memory
bandwidth -- it is a strong-scaling wall, not a bandwidth wall.

**So: run this pgen on ONE node.** 48 ranks divides 144 exactly and is what every previous
run used. Do not ask for 2 nodes; it buys 6% for double the allocation.

**Also:** do not benchmark this pgen from t=0. The first ~400 cycles are ~1.7x more
expensive per cycle than the relaxed state (dt starts at 0.35 and falls to ~0.23 as
convection switches on), so a short probe understates the real throughput. The relaxed
rate is ~0.9 simulated s per wall s at 48 ranks.
