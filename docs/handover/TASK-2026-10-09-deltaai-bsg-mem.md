# TASK for DeltaAI: BSG GPU memory per block + s/cycle at 8 blocks per GPU (mem-1009, user 10-09)

Why: the BSG true-repro run (half-range VET) is to move to Caltech (2 x H200, 141 GB each, 8 MeshBlocks per GPU).
Hard target: 16 blocks + MPI/Kokkos overhead <= ~125 GB per H200, i.e. <= ~15 GB per block. CUDA memory numbers
are what Caltech needs; viper's apudev is jammed. Same test may also run on Delta (A100); first result wins.
**Budget: 1 node (4 GH200), <= 30 min, <= 2 GPU-h.** Work dir `/work/nvme/bivj/jma20/bsg_hrdet_1009/` only.
Account `bivj-dtai-gh`. No production; do not touch other dirs or jobs.

## Builds (same flags as the 10-09 hrdet gate, target `hes_gpu`, `build_inc_deltaai.sh`)
- BO = `hrup-bsg-1009` @ **98835d99a** (you already have it from the hrdet gate; reuse, record md5)
- BN = `mem-1009` @ **af5147cc2** (fork; = 98835d99 + cherry-picks of accel-1009 a07cf032 (vet_gd band only on
  off-rank sides) and b93b6a9d (no explicit face fluxes / u1 under implicit BE, no vet_col_lat sweep arrays under
  vet_gd, no closure-memory slabs); `git fetch fork mem-1009`)

## Files
Inputs: `/work/nvme/bivj/jma20/bsg_hrdet_1009/files` (SETUP.sh of the hrdet gate; `bsg_hr_dc5.athinput` is used).
Scripts: `git fetch fork bsg-files-1009; git archive fork/bsg-files-1009 docs/handover/mem-1009 | tar -x -C $G/bundle`
(`mem_cuda.sh`, `ana_mem.py`).
Kokkos Tools memory-events (optional but wanted; host-only build, ~10 s):
```
cd $G && git clone --depth 1 https://github.com/kokkos/kokkos-tools kt && cd kt/profiling/memory-events && make CXX=g++
KT=$G/kt/profiling/memory-events/kp_memory_events.so
```
If the clone/build fails, run without KT (nvidia-smi peaks only) and say so.

## Job (one sbatch, 1 node, `--gpus-per-node=4 --ntasks-per-node=4 --exclusive --mem=0 --time=00:30:00`)
```
bash $G/bundle/docs/handover/mem-1009/mem_cuda.sh <BO> <BN> $G/files $G/runs/mem_$SLURM_JOB_ID $KT
```
Arms in order: mem_c8 (FULL 16-block mesh, 2 ranks x 8 blocks = Caltech layout, 6 cycles), mem_r4, base_r4
(reduced 4-block mesh, 1 block per GPU, bitwise pair), mem_c4, base_c4 (full mesh, 4 blocks per GPU).
Set SRUN2/SRUN4 to your GPU binding if `--gpus-per-task=1 --gpu-bind=closest` is not what worked for hrdet. An OOM in
an arm is a RESULT (report it), the script continues with the next arm. Smoke rule: if mem_c8 dies with anything but
OOM in the first minute, cancel and report.

## Report (NOTE-2026-10-09-deltaai-bsg-mem.md on this branch)
1. The per-arm `rc= fatal= oom= wall= s/cycle= peak: gpuN=...MiB` lines, the `bitwise base_r4 vs mem_r4` line, md5s,
   job id, GPU-h.
2. For each arm with KT, `python3 ana_mem.py <arm>/<host>-<pid>.mem_events <ncells_per_block> <blocks_per_rank>`
   for ONE rank (ncells_per_block = 262*70*70 = 1283800): paste the table (top 40 lines).
3. Keep the run dirs (viper may ask for the mem_events files).
