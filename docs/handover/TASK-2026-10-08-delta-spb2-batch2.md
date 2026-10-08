# TASK 2026-10-08 (Delta): sp-blend2 cost arms, batch 2 (row-scaled BiCGStab, implicit_lin_scaled)

Prerequisite: batch 1 (TASK-2026-10-08-delta-spb2-batch1.md) has run. Its prepA restart is reused.
If batch 1 has not run yet, run it first: same scripts, same directory.
Budget: 3 jobs x 1 node x <= 15 min, about 2.7 GPU-h.

## Build
- Branch `sp-blend2-1008`, commit **6d913aed** (or the commit that adds this file; the code is the same).
  `bash build_delta.sh he_gpu <sha>` (the build_delta.sh of fork/rt-integration, as in batch 1).
- New key: `rad_m1/implicit_lin_scaled = true`, default false, bitwise identical when off.
  - It weights the BiCGStab inner products by 1/(s_i E^k_i)^2, i.e. it runs BiCGStab on the row-scaled
    system whose residual the cnorm test already measures.
  - Hypothesis: unweighted, the deep interior (E ~ 1e8 x the atmosphere's) sets omega and the Lanczos
    coefficients, so the thin-atmosphere residual barely converges. That would explain the blend's
    inner iterations (98, max 200), breakdowns and line-Jacobi fallbacks.

## Run
```bash
cd <worktree of sp-blend2-1008>/docs/handover/delta-spb2-1008
bash run_batch2.sh /work/nvme/bivj/jma20/delta_1008/bin/athena_he_gpu_<sha>
```
- A4 and A5 are the same 5 arms in two orders, restarting from `$B/spb2/prepA/rst`, cycles 236 -> 276.
- B2 runs AG Car B from t = 0 to cycle 40.

| job | arm | keys |
| --- | --- | --- |
| A4, A5 | blf | none (reference, new binary) |
| A4, A5 | sblf | implicit_lin_scaled=true |
| A4, A5 | mg | implicit_precond=mg |
| A4, A5 | smg | implicit_precond=mg implicit_lin_scaled=true |
| A4, A5 | smg5 | precond=mg, implicit_mg_levels=5, implicit_lin_scaled=true |
| B2 | sblf, smg, blf, mg | as above |

## Report: NOTE-2026-10-08-delta-spb2-batch2.md, pushed to sp-blend2-1008
As batch 1:
- metrics.sh for every arm;
- cmp.py against `blf` of the same job, for A and with CASE=B for B;
- cmp.py between the two `blf` runs (A4 vs A5) and between the batch-1 `blf` (binary f5b26846) and the
  batch-2 `blf`. The second one must be bitwise identical apart from timing: the key is off;
- the GPU-h used.
