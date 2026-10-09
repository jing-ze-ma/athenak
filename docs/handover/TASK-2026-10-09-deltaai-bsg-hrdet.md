# TASK for DeltaAI: BSG half-range determinism gate -- URGENT (user 10-09)

The user wants this gate wherever it runs first: the same gate is queued on viper (gate8n.sh), Raven (gate16.sh)
and Delta (rad-beam-1008 `TASK-2026-10-09-delta-beam-batch3.md`). **Push a NOTE the moment your job starts**
(`docs/handover/NOTE-2026-10-09-deltaai-bsg-hrdet.md` on this branch, `bsg-files-1009`), and the results when done.
Budget ~8 GPU-h. Work dir **`/work/nvme/bivj/jma20/bsg_hrdet_1009/`** only (do not touch `delta_1008/`, which is the
Delta runner's, nor the older bsg_1001/bsg_1002 dirs). Account `bivj-dtai-gh`.

## What is tested
Branch `hrup-bsg-1009` @ **98835d99a** (half-range fix) must be (b) bitwise identical to the reference
`truerepro2-1009` @ **d385de40** with all new keys off (3 cycles), and (a)+(d) deterministic: two identical
30-cycle runs of the half-range input give bitwise identical dumps (bins every 5 cycles). BSG 3-D (Ma+2026 true
repro), 256^3, 16 MeshBlocks of 256x64x64, 1 block per rank = **16 ranks**.

## Steps
1. **Builds** (two binaries, same flags): `git fetch fork hrup-bsg-1009 truerepro2-1009 bsg-files-1009`; build
   98835d99a and d385de40 with `build_inc_deltaai.sh` target `hes_gpu` (PROBLEM=he_star_m1, HOPPER90 + ARMV9_GRACE,
   MPI on, default FMA; stack as TASK-2026-10-07-deltaai-hegiant sect. 2). Record both md5s.
2. **Files** (the 10-01 bundle used TOPS tables; this needs the ma2026 IC/tables, 12 MB gz here):
   ```
   G=/work/nvme/bivj/jma20/bsg_hrdet_1009; mkdir -p $G/bundle
   git archive fork/bsg-files-1009 docs/handover/bsg-hrdet-1009 | tar -x -C $G/bundle
   bash $G/bundle/docs/handover/bsg-hrdet-1009/SETUP.sh $G/files      # 3 x OK, SETUP_OK
   ```
3. **Job** (one sbatch; `gate_hrdet.sh` runs the 4 arms in sequence, ~25 min wall on 16 GPUs):
   `bash $G/bundle/docs/handover/bsg-hrdet-1009/gate_hrdet.sh <d385de40 binary> <98835d99 binary> $G/files $G/runs/$SLURM_JOB_ID`
   with `module reset` stack, `-A bivj-dtai-gh`, `--time=00:30:00`, `--mem=0`, `--exclusive`. Layout choice:
   - **(A) 4 nodes x 4 GH200, 4 ranks per node**, 1 rank per GPU, using the GPU binding of your 10-01 BSG gate /
     10-08 hegiant smoke (`--gpus-per-task=1 --gpu-bind=closest`, SRUN env of gate_hrdet.sh accordingly).
     Preferred if ghx4-interactive (or ghx4) allows 4 nodes and starts soon (check `scontrol show partition`
     MaxNodes and `sbatch --test-only`).
   - **(B) 2 nodes x 8 ranks**, 2 ranks (2 blocks) per GH200 (`--ntasks-per-node=8 --gpus-per-node=4`, no
     per-task binding; KOKKOS_MAP_DEVICE_ID_BY=mpi_rank maps local rank mod 4). Memory estimate: viper MI300A
     (unified) used ~66 GB per 2 blocks, i.e. ~33 GB per block -> 2 blocks ~66 GB of the 96 GB HBM, host share
     well under the node's LPDDR. Use (B) if 4 nodes cannot start within ~1 h. Same 16-rank decomposition, so the
     gate is equivalent. Check `nvidia-smi --query-gpu=memory.used --format=csv` once during `bit_old`.
   Smoke rule: arm 1 `bit_old` (3 cycles) is the smoke; if rc != 0 or FATAL (or CUDA OOM), cancel and report.
4. **Metrics** (all in the NOTE): per arm the `rc= fatal= diverged= nonconv= wall=` line; NaN check of the hst;
   s/cycle (wall between cycle 1 and the last cycle / cycles); bytewise `cmp` of every `bin/*.bin` and `*.hst`
   for bit_old vs bit_new and detH_a vs detH_b, plus
   `python3 $G/bundle/docs/handover/bsg-hrdet-1009/bitcmp.py A B <athenak>/vis/python` (prints BITWISE_PASS or the
   differing variables with max |diff|); number of files compared per pair; binary md5s; job id, layout (A/B),
   start/elapsed (relative), GPU-h. Verdict line first: `bit: PASS|FAIL  det: PASS|FAIL`.
5. Keep the run dirs (viper may ask for a bin). No production, nothing else queued.
