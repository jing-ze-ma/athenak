# NOTE 2026-11-25y orion -> Caltech: memory of the new tip on HG, and HG on ONE node

caltech-mem50-9b9d05f90 = 79c1cabef + the memory commits (bitwise on CPU and MI300A: R136 33.4 -> 19.5 GB per GPU). Please build star_vet (tiny gate 3.686e-11) and, after 25w and around 25x on 02-03 / 02-09 (or any free one-node slot for part 2):
1. HG as before (8 ranks, 2 nodes, T6 environment), 26 cycles, 79c1cabef against 9b9d05f90: hst printed-bitwise, nvidia-smi peak memory per GPU (poll every 2 s) and the number at the end, s/cyc.
2. **HG on ONE node, 4 ranks (6 MeshBlocks per GPU)** with 9b9d05f90: 26 cycles of HGN's command (it ran out of memory with the old arrays: 141 GB). If it fits: peak memory per GPU, s/cyc (median), MPI wait; hst bitwise with part 1's. If it does not: the OOM line and the last memory sample. On one node the inter-node transport disappears, so the s/cyc is the interesting number. The old 79c1cabef on 4 GPUs need not be tried (it ran out of memory).
3. hst / rst: no production, no rst kept; stop by 08:00 PST.

-- orion (coordinator)
