# NOTE viper -> caltech (interim reply to TASK 8513c5b56): Arepo2-IDORT g25_full on 4 Raven GPU nodes

**Run:** Raven job 31511100, 4 nodes (288 ranks, 16 A100), g25_full IC (6.09 M cells), refinement block removed (as run_c.sh),
4 top steps, top steps 2-4 reported. Branch idort-gpu-nltree e93504ccef (the gated combined branch, incl. the g25
direct-summation fix). Both arms in the same allocation, same build process. The 8-node job (31511099, 576 ranks, your requested
size) is still queued, est. start in ~17 h.

| per top step (s) | GPU recipe (g2m) | CPU, same branch (c1m) |
|---|---|---|
| **whole top step** | **20.94** (19.31-22.45) | 28.45 (26.59-29.77) |
| IDORT full sync | 6.15 | 15.23 |
| IDORT partial steps (CPU) | 2.01 | 1.98 |
| mesh | 5.12 | 5.08 |
| gravity | 2.91 | 1.47 |
| rest | 4.74 | 4.69 |
| imbalance max/mean | 9.3 | 6.1 |
| peak node mem / GPU mem per node | 135 GB / 117 GB | 117 GB / - |

**Reading:**
- GPU is 1.36x the CPU build on the same nodes; the IDORT full sync is 2.5x faster on the GPU.
- Against your Orion figure of 35.6 s per top step on 560 cores: 20.9 s on 288 cores + 16 A100, i.e. ~1.7x faster with about
  half the cores. Caveat: no refinement in our run, and Orion's cell count over its run may differ.
- Gravity is 2x slower in the GPU arm (2.91 vs 1.47 s) with a higher imbalance (9.3 vs 6.1). Not diagnosed yet.
- Not in this run: the O(active) partial-step work (local branch idort-oactive; on g25_red 1 node it cuts the partial steps
  by 19-22 % and the whole top step by 7-10 %, bitwise) and the NodeBarrier fence fix (no cost on Raven). The merge of both is
  being gated on Raven now.
- The classic-CPU arm (c0) was not in this 4-node job; the 8-node job has it.
