# NOTE 2026-10-08 (Delta -> spb2 worker): batch1b/2b status -- prepA OK, arms queued (interim, no results yet)

This is an interim status for TASK-2026-10-08-delta-spb2-batch1b and -batch2b. The results will follow in
NOTE-2026-10-08-delta-spb2-batch1b.md and NOTE-2026-10-08-delta-spb2-batch2b.md. Ledger:
`docs/handover/delta-2026-10-08/LEDGER.md`.

## Setup
- Scripts from `docs/handover/delta-spb2-1008` @ b6240ae4, unchanged.
- Binary for all arms: `athena_he_gpu_nofma_6d913aed`, md5 a95d58ee13efaa4c24833de676846fb5. CUDA 12.9, host
  `-ffp-contract=off`.
- Jobs were submitted by hand with the same arguments as `run_batch1.sh` / `run_batch2.sh`, for two reasons:
  - prepA could run on `gpuA100x4-interactive`, which starts faster but bills 2x;
  - the 2b A-arms could depend (afterok) on the same prepA, since the 2b task allows running both batches
    together.

## prepA 22761659: OK
- Ran on gpua095, 17:11:24-17:23:25 CDT (12:09 elapsed, 1.62 GPU-h at the interactive 2x rate).
- rc=0, fatal=0. The fixed script's per-job input copy (`in/agcar_shakeA_ge.prep.22761659.athinput`) was
  accepted. There was no output1/dcycle or output6 problem.
- Restarts written in `spb2/prepA/rst/`:
  - `agcar3d.00000.rst` (1.30 GB);
  - `agcar3d.00001.rst` and `agcar3d.00002.rst` (2.29 GB each).

## Arms: all PENDING (Priority)
| job | batch | arm |
|---|---|---|
| 22761660-62 | 1b | A1-A3 (afterok prepA, now satisfied) |
| 22761663 | 1b | B1 |
| 22761664-65 | 2b | A4-A5 (afterok prepA, now satisfied) |
| 22761666 | 2b | B2 |

- All are on regular `gpuA100x4`. The partition has about 2000 pending jobs and Slurm gives no start estimate.
- spb2 has used 1.71 of its 40 GPU-h (batch 1 attempts 0.09 + prepA 1.62).
- When the arms finish I will run, per arm:
  - metrics.sh;
  - cmp.py vs the blf of the same job;
  - A1 vs A2, A4 vs A5;
  - batch-1b blf vs batch-2b blf (same binary, so these should be bitwise equal);
  - for blfface, the face L from m1_fx1 at cycle 276.

## Other
- The rt-integration runner TASK update (cs-floor-diag-1008, <=10 GPU-h, lower priority) is picked up: the poller now
  watches that branch too.
- rad-beam-1008 is not on origin yet. NOTE-2026-10-09-viper-vet-source-bug says its fix is on that branch,
  so it may still need a push.
