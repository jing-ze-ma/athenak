# TASK 2026-10-08 (Delta): sp-blend2 cost arms, batch 1b (re-issue of batch 1 with fixed scripts)

This supersedes TASK-2026-10-08-delta-spb2-batch1 (stopped: NOTE-2026-10-08-delta-spb2-batch1.md). Only the
scripts in `docs/handover/delta-spb2-1008/` changed (this commit).

## What was fixed
- spb2_delta.sh writes a per-job COPY of `$B/agcar_files/agcar_shake{A,B}_ge.athinput` into `$B/spb2/in/`. The
  copy contains every key that the command line overrides:
  - `dcycle` in output1-3;
  - `dcycle` in output5 (the input's rst block; only one rst block is allowed);
  - an output6 = m1_face bin block;
  - `implicit_blend`, `implicit_lin_scaled`, `implicit_mg_levels`, `implicit_precond_float` in `<rad_m1>`.
- The script exits with a non-zero code when a run fails (rc or FATAL), so `afterok` holds the dependent arms.
  prep also fails when no restart was written.
- `srun ... < /dev/null`, so srun cannot eat the arms file.
- vis/python/bin_convert.py is taken from 03590601 (it ignores '#' comments in the header), which cmp.py needs.
- **Tested** on the viper login node with the exact script. Only W/IN0 were re-pointed and srun was shimmed
  to one CPU process, with a 480x8x8 mesh via arm keys (CPU binary he_cpu 6d913aed):
  - prepA (nlim 2) wrote rst;
  - A arms blf / sblf / blfface (restart, nlim 4) ran, with the m1_face bin written;
  - B arms blf / smg (fresh, nlim 2) ran;
  - all rc 0, no FATAL; metrics.sh and cmp.py (A and CASE=B) ran on the outputs.
  - Test dir: viper /viper/ptmp2/jinma/spblend2_1008/dtest.

## Build
Use the binary already built for batch 2: **`athena_he_gpu_nofma_6d913aed`** (md5 a95d58ee...).
With `implicit_lin_scaled` off it is the same arithmetic as f5b26846, and using it for every arm keeps one
binary. Code is unchanged since 6d913aed; this commit changes docs and vis/python only.

## Run
```bash
cd <worktree of sp-blend2-1008 at this commit>/docs/handover/delta-spb2-1008
bash run_batch1.sh /work/nvme/bivj/jma20/delta_1008/bin/athena_he_gpu_nofma_6d913aed
```
- Arm table, metrics and report: as in TASK-2026-10-08-delta-spb2-batch1.md.
  - The blfface arm now sets `output6/dcycle=276`; output6 is the m1_face block.
  - For blfface, report the face L from `m1_fx1` in `*.m1_face.*.bin` at cycle 276.
- About 4.5 GPU-h.
- Report in NOTE-2026-10-08-delta-spb2-batch1b.md.
- Then run batch 2b (TASK-2026-10-08-delta-spb2-batch2b.md), which needs prepA.
