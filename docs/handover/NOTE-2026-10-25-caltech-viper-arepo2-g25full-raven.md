# TASK caltech -> viper (USER 10-25): g25_full scaling test on Raven, 8 GPU nodes

No Arepo source in git (as before). Follows fa979b38.

## Why Raven
Caltech g25_full (1 H200 node, 32 ranks): GPU recipe 83.4 s per top step vs pure CPU 234.6 s; ~55 s of the 83 s is
CPU-bound partial-step work, O(NumGas) per rank. Orion reference: 35.6 s per top step on 560 cores. A fair comparison needs
~560 CPU ranks plus the GPUs, i.e. several full nodes. Caltech cannot provide that now: our held allocations are 32-core GPU
slices, and the 1-node full job could not start before 10-27, so the USER cancelled it (5436311). Raven fits: 8 GPU nodes =
576 ranks + 32 A100.

## Request
- Code: the combined branch `idort-gpu-nltree` (idort-gpu-s1 + nltree-g25 + the g25 direct-summation fix if real), merged
  with our returned branches (bundle `caltech_return_1025`, see fa979b38) once gated.
- Problem: g25_full (not g25_red), same snapshot/parameters as the Caltech 32-rank numbers.
- Layout: 8 Raven GPU nodes, 72 ranks per node = 576 ranks, 4 A100 per node shared by the node's ranks.
- Arms (same build options otherwise):
  1. GPU recipe (GPU IDORT + GPU mesh + GPU node-local tree);
  2. CPU classic tree;
  3. CPU node-local tree.
- Length: enough top steps for a stable mean (at least 3 top steps after the first).
- Report per arm: seconds per top step (mean, min, max), split into IDORT, mesh, gravity tree, CPU partial steps (the rest),
  plus peak memory per node. Also the per-rank load imbalance of the partial steps if available.
- If memory allows, one extra GPU-recipe arm on 4 nodes (288 ranks), for scaling.

## Reply
A NOTE on bsg-files-1009 with the numbers, the commit of the code used, and the job scripts' paths on Raven.
