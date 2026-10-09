# NOTE 2026-10-09 Caltech: GPU capacity (answers TASK-2026-10-09-caltech-capacity; nothing submitted, BSG untouched)

Snapshot 10-09 ~10:10 PDT (Resnick HPC, account carnegie_poc).
## 1. What the account can use
- One shared GPU partition `gpu` (MaxTime 14 d; QOS debug/long/normal). Nodes, 4 GPUs each:
  - **18 x H200 141 GB** (hpc-sm-01-xx, hpc-sm-02-xx; one drained; hpc-sm-01-09 excluded by us: bad node),
  - **4 x H100 80 GB** (hpc-33-13/16/19/22),
  - 4 x L40S 48 GB (poor FP64: not useful for our double-precision runs), ~45 x P100 16 GB, 2 x V100 (old/small).
- Fair share: carnegie_poc NormShares 0.074, effective usage 0.102; our user FairShare 0.25 (we are above share -> lower priority).
- Queue: 144 pending / 169 running jobs in `gpu` at the snapshot.
## 2. Free now: none of H200/H100
- Every H200 and H100 GPU is allocated (BSG holds 2 of hpc-sm-02-11; the other 2 there are another user's). Free GPUs
  exist only on P100 (16 GB: an AG Car B block needs ~20 GB, so no) and L40S (FP64-poor).
- `sbatch --test-only` start estimates (1 node, -t 12 h and 24 h give the same):
  | request | estimated start |
  |---|---|
  | 2 x H200 | 2026-10-10 14:23 PDT (~28 h) |
  | 4 x H200 | 2026-10-10 15:18 PDT |
  | 2 x H100 | 2026-10-11 00:41 PDT |
  | 4 x H100 | 2026-10-10 22:28 PDT |
- Caveat from experience: these estimates are pessimistic. Backfill started our recent jobs much sooner (BSG 12 h link
  1 started 20 s after submission at 09:57; 15-30 min smokes typically within minutes). A 12 h 2-H200 job submitted now
  might start within hours or not until tomorrow; shorter links (4-6 h) backfill more easily.
## 3. Suggestion
- AG Car B (4 blocks x ~20 GB): fits 1 x 2 H200 (2 blocks per GPU, ~40 GB/GPU) or 4 H100; submit with short links.
- He giant fresh start: N897 needs >= 2 H200 (1 H200 OOMs at a 119 GiB allocation); the N445 fresh start ran on 1 H200.
- Running both next to BSG means 4-6 more H200 GPUs; expect queue waits as above. Waiting for viper's task.
