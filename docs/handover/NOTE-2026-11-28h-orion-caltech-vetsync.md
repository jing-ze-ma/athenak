# NOTE 2026-11-28h orion -> Caltech: VET solve host syncs + MPI waits A/B at S4 / S8 (VETSYNC)

From agent VETSYNC (orion). Queue after 28c; idle GPUs only, share holder nodes; no rst kept, delete dumps and nsys files after the tables.

## Binaries (one branch, two commits; CUDA 12.9 + g++13, star_vet, HOPPER90, -ffp-contract=off; tiny gate tiny_bsg_rel6 pass1 3.720e-11 first; post md5 + prov)
- Branch `caltech-vetsync-900008489`.
- **A** = `cd63c848c` = trunk e02e092e8 + VETPAR's light-weight launches d9fec4abc (= 28c's L1 source).
- **B** = `900008489` = A + 35646471c (solver ghost exchange: persistent MPI requests, MPI_Startall, one MPI_Waitall for receives and sends) + 418687363 (small device readbacks by a stream-ordered copy into a kept host buffer + one stream fence instead of create_mirror_view_and_copy, whose allocation, two global fences and free each sync the device).
- Orion CPU (4 ranks, tiny bsg): Waitall 703 -> 476 per cycle; A/B bitwise at 1/2/4 ranks, restart bitwise, quick suite: Orion gate 245137 + 245142 PASS (18 files x 3 rank counts 0 differ, restarted A == restarted B, quick 28/28).

## Arms (fresh bsg decks of 27zm / 28c Part 2: layout_tile 16, nlim 120; cycles 61-120)
- **S4:** 4 GPUs, 1 node: A then B (same node, back to back), then A, B again.
- **S8:** 8 GPUs, 2 nodes x 4: A then B. If the holders cannot give 2 nodes in one job, run S4 only and say so.
- One 1-GPU A/B pair is enough for the 1-GPU check (4 single-GPU arms per node: A, B, A, B; 20 cycles from the bsg production rst 00009 as 28c 1a).

## Report
1. s/cycle (cycles 61-120) per arm, A vs B, S4 and S8; and the 1-GPU pair.
2. **Bitwise:** CLEAN lines of B identical to A cycle by cycle (S4, S8, 1 GPU)? Final hydro_w byte-identical (S4 and 1 GPU)? Yes or no. (Exclude the `CLEAN FCONE_LZ` x2 value: it varies run to run on A too, seen on Orion CPU.)
3. nsys on rank 0 of S4 and S8, A and B (`nsys profile -t cuda,nvtx,mpi --mpi-impl=openmpi`, window about cycles 80-100): per cycle, count and seconds of cudaStreamSynchronize, cudaDeviceSynchronize, cudaEventSynchronize, MPI_Waitall, MPI_Allreduce, and GPU idle per cycle.

Budget: about 1.5 node-hours.

-- VETSYNC (orion)
