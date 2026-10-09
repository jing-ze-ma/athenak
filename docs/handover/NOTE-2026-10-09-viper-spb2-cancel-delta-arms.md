# NOTE 2026-10-09 (viper spb2 worker -> Delta runner): CANCEL the batch 1b/2b arms

**Action on Delta:**
```bash
squeue -u $USER -h -o "%i %j %T"          # fresh check first
scancel --state=PENDING 22761660 22761661 22761662 22761663 22761664 22761665 22761666
```
- Cancel only the PENDING ones. If an arm is already RUNNING, let it finish and report it in a short NOTE, with no
  further analysis needed.
- Do not run TASK-2026-10-08-delta-spb2-batch1b / -batch2b any further, and do not resubmit them. Both are closed.
- Log the cancellation in `docs/handover/delta-2026-10-08/LEDGER.md`: spb2 used 1.71 GPU-h, the prepA restart
  included.
- Keep the remaining GPU-hours for new batches. The spb2 worker will queue a new Delta batch only if it tests
  something Raven cannot answer sooner.
- `$B/spb2/prepA/rst` can stay; a later batch may reuse it.

**Why:** Raven packs 31012413 and 31012857 (4 x A100, the same AG Car A/B states) answered every arm of 1b/2b:
- `implicit_lin_scaled` is rejected;
- `implicit_precond = mg` is not robust on B (Picard NON-CONVERGED, interior error);
- `implicit_bcg_fallback = best` (7c49ee4b) is adopted: B is 9.5 % cheaper and identical to 1e-7.

Details are in /viper/ptmp2/jinma/spblend2_1008/SPBLEND2.md.
