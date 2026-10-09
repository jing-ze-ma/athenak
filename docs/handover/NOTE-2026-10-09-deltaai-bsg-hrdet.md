# NOTE DeltaAI: BSG half-range determinism gate -- JOB STARTED (results to follow in this file)

Answer to TASK-2026-10-09-deltaai-bsg-hrdet (bsg-files-1009 ba268b1d).

- **Status: job 3347471 (nodes gh[049,082]) RUNNING since 2026-10-09 09:18 CDT (16:18 CEST)** (DeltaAI ghx4-interactive), ~25-50 min wall expected.
- **Layout (B): 2 nodes x 8 ranks = 16 ranks, 2 MeshBlocks per GH200.** (A) not possible soon: ghx4-interactive
  MaxNodes=2, and a 4-node ghx4 job tests to start on 10-14. Binding: `srun -n 16 --ntasks-per-node=8 -c 36
  --cpu-bind=cores` + wrapper `CUDA_VISIBLE_DEVICES=SLURM_LOCALID/2` (local ranks 2g,2g+1 on Grace socket g and
  GH200 g; Kokkos then sees 1 device per rank). `--exclusive --mem=0 --time=01:00:00` (1 h instead of 30 min,
  because 2 blocks/GPU may run more slowly; only the time used is charged).
- Builds (build_inc_deltaai.sh hes_gpu, PROBLEM=he_star_m1, HOPPER90+ARMV9_GRACE, MPI on, default FMA, Cray CC host):
  - BO `athena_hes_gpu_d385de405fa8` (truerepro2-1009 d385de40) md5 `a66383005fed09b3b4f9409f74f13ab1`
  - BN `athena_hes_gpu_98835d99a16f` (hrup-bsg-1009 98835d99a) md5 `c9c6d16705d07dc6e97de6c808d02182`
- Files: SETUP.sh 3 x OK, SETUP_OK (`/work/nvme/bivj/jma20/bsg_hrdet_1009/files`).
- Gate script = bundle gate_hrdet.sh unchanged except one line after arm 1: abort if bit_old rc != 0 or FATAL (smoke
  rule). GPU memory is sampled every 20 s on both nodes during the whole job (gpumem.log).
