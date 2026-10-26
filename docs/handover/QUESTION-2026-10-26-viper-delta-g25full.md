# QUESTION viper -> Delta (NCSA Delta, not DeltaAI) (USER): can you run Arepo2-IDORT g25_full? Reply first, do not start

The user asks whether **Delta** can run the Arepo2-IDORT g25_full timing (6.09 M cells), the measurement of TASK 8513c5b56.
Raven's 4-node backup is queued (job 31546612); DeltaAI has a GO for a 2-node GH200 job (TASK of 10-26); Caltech is full.

**Why Delta may fit:** gpuA100x4 nodes are the same layout as Raven's GPU nodes (4 x A100-40GB per node), and your 10-08
bring-up saw short jobs backfill in seconds. Raven's 4-node result (NOTE 01508663c): GPU 20.94 vs CPU 28.45 s per top step.

**What is already on /work (shared with DeltaAI):**
- `/work/nvme/bivj/jma20/arepo2_gpu_1025/`: DeltaAI's tree with the 10-24 Arepo bundle, the g25 inputs (g25_inputs.tgz md5
  2348f0db65cc644446240951bd518235, unpacked) and in `bundle/` the new incremental bundle
  `arepo_idort_gpu_nltree_2_1026_from1024.bundle` (md5 ee4103cfe2926091b085b0ed62165ac3; head idort-gpu-nltree-2 d9c0152bce,
  prerequisites 1108335a7 + 6cef4d1d6 from the 10-24 bundle). Read only: clone into your own dir (e.g. delta_1008/arepo2/).
  Do NOT build in or write to DeltaAI's tree (aarch64 builds there).
- No Arepo source or diffs in git, ever.

**The run (if yes):** 4 nodes gpuA100x4 (256 ranks, 64 per node, 4 islands per node), 3 arms one after another in one
allocation: GPU recipe g2m, CPU same branch c1m, CPU classic c0m; g25_full IC, refinement block removed, TimeMax = 4 top steps,
report steps 2-4 (param as Raven comb/run_full.sh; DeltaAI can point you to its copy). Memory: Raven peaked at 135 GB host and
117 GB GPU per node at 4 nodes (72 ranks per node), which fits 250 GB / 4 x 40 GB.

**Please reply with a NOTE:**
1. Can you run it (Arepo2 build on Delta: x86 + CUDA sm_80, cray-mpich, HYPRE/GMP/GSL/HDF5 available?), on 4 nodes?
2. Queue estimate now (sbatch --test-only) for 4 nodes x 1 h on gpuA100x4, and the cost against the 99 GPU-h balance
   (~16 GPU-h for 4 nodes x 1 h).
3. Anything blocking. Do not start before the user confirms.
