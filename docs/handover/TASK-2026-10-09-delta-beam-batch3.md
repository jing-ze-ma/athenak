# TASK for the Delta runner: beam batch3 = BSG half-range determinism gate -- URGENT, PROCESS FIRST

**URGENT (user 10-09): run this batch before any other waiting spb2 / beam / csfloor batch.** It is a 4-node job
(not the usual 1-node arm), so it is an exception to the runner's "1 node x 4 A100" rule; the user approved it.
Budget for this batch: ~8 GPU-h (charge it to the beam line of the ledger). Same gate is queued on viper
(`/viper/ptmp2/jinma/hrdet_1009/gate8n.sh`) and Raven (`gate16.sh`); the user wants whichever finishes first.

## What is tested
Branch `hrup-bsg-1009` @ **98835d99a** (the half-range fix) must be (b) bitwise identical to the reference
`truerepro2-1009` @ **d385de40** with all new keys off, and (a)+(d) deterministic: two identical runs of the
half-range input give bitwise identical dumps. BSG 3-D (Ma+2026 true repro), 256^3, 16 MeshBlocks of 256x64x64,
1 block per GPU = 16 ranks.

## Steps
1. **Builds (both shas, two binaries).** `git fetch origin hrup-bsg-1009 truerepro2-1009` (both on the fork), then
   `bash build_delta.sh he_gpu 98835d99a` and
   `bash build_delta.sh he_gpu d385de40` (your runner copy of docs/handover/delta-2026-10-08/build_delta.sh; default FMA target `he_gpu`; same flags
   for both, which is what the bitwise arm needs). The second build is incremental in the same build dir; the
   binaries land in `bin/athena_he_gpu_98835d99` and `bin/athena_he_gpu_d385de40`. Record both md5s.
2. **Files** (not on Delta before: the 10-01 DeltaAI BSG bundle used other TOPS tables). Branch
   **`bsg-files-1009`** (orphan, files only, 12 MB gz):
   ```
   git fetch origin bsg-files-1009
   mkdir -p /work/nvme/bivj/jma20/delta_1008/bsg_hrdet/bundle
   git archive FETCH_HEAD docs/handover/bsg-hrdet-1009 | tar -x -C /work/nvme/bivj/jma20/delta_1008/bsg_hrdet/bundle
   bash /work/nvme/bivj/jma20/delta_1008/bsg_hrdet/bundle/docs/handover/bsg-hrdet-1009/SETUP.sh \
        /work/nvme/bivj/jma20/delta_1008/bsg_hrdet/files        # must print 3 x OK and SETUP_OK
   ```
   This writes `bsg3d_truerepro2.athinput` (keys-off reference input) and `bsg_hr_dc5.athinput` (half-range input)
   with absolute paths to the gunzipped IC (`ic_ma2026_f0.txt`, md5 246a666e) and tables (Rosseland 450bc0c1,
   Planck e18ca805).
3. **Job** (one sbatch; 4 arms run in sequence inside it by `gate_hrdet.sh`):
   ```
   #!/bin/bash
   #SBATCH -J bsghrdet
   #SBATCH -A bivj-delta-gpu
   #SBATCH -p gpuA100x4
   #SBATCH --nodes=4
   #SBATCH --ntasks-per-node=4
   #SBATCH --cpus-per-task=16
   #SBATCH --gpus-per-node=4
   #SBATCH --exclusive
   #SBATCH --mem=0
   #SBATCH --time=00:25:00
   #SBATCH -o /work/nvme/bivj/jma20/delta_1008/bsg_hrdet/j.%j.out
   module unload cudatoolkit; module load cuda/12.9
   G=/work/nvme/bivj/jma20/delta_1008/bsg_hrdet
   bash $G/bundle/docs/handover/bsg-hrdet-1009/gate_hrdet.sh \
     /work/nvme/bivj/jma20/delta_1008/bin/athena_he_gpu_d385de40 \
     /work/nvme/bivj/jma20/delta_1008/bin/athena_he_gpu_98835d99 \
     $G/files $G/runs/$SLURM_JOB_ID
   ```
   Memory: 4 blocks per node. On viper (MI300A, unified) 2 blocks used ~66 GB, i.e. ~33 GB per block, which fits
   one A100-40GB per block and 4 x 33 = 132 GB of the 256 GB node RAM. If an arm dies of HOST-RAM OOM, resubmit
   on 8 nodes (`--nodes=8 --ntasks-per-node=2 --gpus-per-node=2`, still 16 ranks, 1 block per GPU, ~13 GPU-h:
   over budget, so ask first). If it dies of CUDA out-of-memory (device), do not retry: report it.
   Smoke rule: the first arm `bit_old` (3 cycles) is the smoke; if it shows rc != 0 or FATAL, cancel the job and
   report (the arms run sequentially, so nothing else is wasted).
4. **Metrics** (put all in the NOTE):
   - per arm the `rc= fatal= diverged= nonconv= wall=` line of the job output, plus
     `grep -ci NONCONV`, `grep -ci DIVERGED`, NaN check (`grep -ci nan` in the hst columns) for each `out.log`;
   - s/cycle per arm (from the `cycle=` diag lines: wall between cycle 1 and the last cycle / cycles);
   - bytewise: `cmp` every file under `bin/` and every `*.hst` between bit_old and bit_new, and between detH_a and
     detH_b (`for f in $(cd A && find . -name '*.bin' -o -name '*.hst'); do cmp -s A/$f B/$f || echo DIFF $f; done`),
     and also `python3 bundle/docs/handover/bsg-hrdet-1009/bitcmp.py A B <athenak>/vis/python` (prints
     `BITWISE_PASS` or the differing variables with max |diff|). Note: bin headers carry no wall-clock, so `cmp` is
     expected to agree with bitcmp.py; if only `cmp` differs, say which bytes.
   - number of bins compared per pair (bit: dumps at t >= 262 s or the final; detH: cycles 0,5,...,30).
   - binaries' md5, job id, nodes, start/elapsed (relative), GPU-h.
5. **NOTE back** on this branch (`rad-beam-1008`) as `docs/handover/NOTE-2026-10-09-delta-beam-batch3.md`, with the
   LEDGER line, pushed to origin. Verdict line first: `bit: PASS|FAIL  det: PASS|FAIL`.
   Do not delete the run dirs (viper may ask for one bin).
