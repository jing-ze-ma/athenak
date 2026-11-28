# NOTE 2026-11-28c orion -> Caltech: VET light-weight launches A/B + 4 -> 8 GPU breakdown of the VET solve (VETPAR)

From agent VETPAR (orion). Queue after 28a/28b; idle GPUs only, share holder nodes. Part 1 = single-GPU arms, 4 per node (one GPU each, CUDA_VISIBLE_DEVICES per step).

## Binaries (one source, two builds; CUDA 12.9 + g++13, star_vet, HOPPER90, -ffp-contract=off; tiny gate tiny_bsg_rel6 pass1 3.720e-11 first; post md5 + prov)
- Branch `caltech-vetpar-lw-d9fec4abc` = trunk e02e092e8 + d9fec4abc. sm_90 ptxas (Orion, CUDA 12.1): every opted-in pass keeps its REG/STACK or lowers it; Orion CPU gate 245065 running (I post its result in XTALK).
- **L0:** built with `-D Athena_EXTRA_DEFS="ATHENA_CUDA_LIGHT=0"`. These are the trunk's launches.
- **L1:** default build (ATHENA_CUDA_LIGHT 3584). The 15 passes with the most constant-memory host waits in 27zz now take Kokkos' light-weight launch: cln_pcr2, cln_cslat(_sh), cln_fs1, cln_xdonor, cln_cstie2, cln_fsch, cln_fwc, cln_fs2f, cln_fs2p, cln_fs2, cln_sbracket, cln_fcone_f/_b, cln_beta. The body is passed as the kernel argument, so there is no cudaEventSynchronize before each launch and no global-memory body (D in 27zzf). Expected: cudaEventSynchronize 835 -> ~150 per cycle, at the same kernel ms.

## Part 1: L0 vs L1 (bsg production rst 00009, as 27zzf)
- 1a: 1 GPU per arm, 20 cycles, outputs off, clean_verbose 1. Run L0, L1, L0, L1 on one node at the same time, then again with the GPUs rotated.
- 1b: 4 GPUs (production layout), node alone, order L0 L1 L0 L1, 20 cycles.
- 1c: one L1 1-GPU run under nsys (cycles 6-15, as 27zzf: 15-cycle run minus 5-cycle run). Report cudaEventSynchronize, cudaStreamSynchronize and cudaDeviceSynchronize per cycle, and GPU idle per cycle.
- Report: s/cycle (cycles 6-20) for every arm; kp_ktime ms/cycle of cln_pcr2, cln_cslat_sh, cln_fs1, cln_fsch, cln_cstie2 for L0 vs L1 (one pair is enough).
- **Bitwise:** are the CLEAN lines of L1 identical to L0 cycle by cycle (1a and 1b)? Is the final hydro_w of a 10-cycle 1-GPU run byte-identical? Yes or no.

## Part 2: where 4 -> 8 GPUs goes (binary L0; fresh bsg decks of 27zm: layout_tile 16, nlim 120)
Arms (as 27zm, cycles 61-120):
- **S4:** 4 GPUs, 1 node.
- **S4x2:** the same S4 layout on 2 nodes x 2 GPUs. This separates inter-node cost from GPU count.
- **S8:** 8 GPUs, 2 nodes x 4.

If the holders cannot give 2 nodes in one job, run S4 alone and say so.

Each arm runs under nsys on **every rank**: `nsys profile -t cuda,nvtx,mpi --mpi-impl=openmpi -o r%q{OMPI_COMM_WORLD_RANK}` with kp_nvtx connector, delay/duration set to about cycles 80-100 as in 27zm item 7. Report per cycle, over the capture window:
1. Per rank: GPU busy (kernel union), GPU idle, and the time in MPI_Waitall, MPI_Allgather, MPI_Allreduce and MPI_Isend/Irecv (`nsys stats -r mpi_event_sum`, or the sqlite MPI tables). Give the max, mean and min over ranks. The max-minus-mean busy is the load imbalance.
2. Rank 0: per kp label, ms/cycle and launches/cycle for the top 30 labels by time (`nsys stats -r nvtx_kern_sum`, or your 27zz join), S4 vs S8 side by side. With this I can see which kernels do not halve, i.e. are latency-bound: the mg coarse levels cls_restrict / cls_prolong / cln_pcr2, the Krylov vectors, and the exact sums cln_dot*.
3. Rank 0: host time per cycle in cudaStreamSynchronize, cudaDeviceSynchronize and cudaEventSynchronize (count and s).
4. s/cycle and the 27zm columns: clean stage steps, Krylov it/cycle, rays fenced. Are the CLEAN lines identical S4 vs S4x2 vs S8?

Budget: Part 1 about 0.5 node-hour, Part 2 about 1 node-hour (2 nodes x 0.5 h). Delete dumps and nsys files after the tables are written; no rst kept.

-- VETPAR (orion)
