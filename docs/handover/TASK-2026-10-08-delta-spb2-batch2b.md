# TASK 2026-10-08 (Delta): sp-blend2 cost arms, batch 2b (implicit_lin_scaled; re-issue of batch 2)

Supersedes TASK-2026-10-08-delta-spb2-batch2. Same arms, but run with the fixed scripts of batch 1b (this
commit). Prerequisite: batch 1b prepA has written `$B/spb2/prepA/rst`.

```bash
cd <worktree of sp-blend2-1008 at this commit>/docs/handover/delta-spb2-1008
bash run_batch2.sh /work/nvme/bivj/jma20/delta_1008/bin/athena_he_gpu_nofma_6d913aed
```
- 3 jobs, about 2.7 GPU-h:
  - A4 and A5: arms blf, sblf, mg, smg, smg5 in two orders, restarting from prepA, cycles 236 -> 276;
  - B2: arms sblf, smg, blf, mg, from t = 0 to cycle 40.
- Keys, metrics and report as in TASK-2026-10-08-delta-spb2-batch2.md.
- Report in NOTE-2026-10-08-delta-spb2-batch2b.md.
- If batch 1b and 2b fit, they may run concurrently: the A4/A5 arms only need the prepA restart.
