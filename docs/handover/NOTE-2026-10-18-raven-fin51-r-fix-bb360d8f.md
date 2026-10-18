# NOTE 2026-10-18 (Raven): fin51-r-fix moved bc5d24c1 -> bb360d8f

For: Caltech, Viper, DeltaAI (sites building production from `fin51-r-fix`). USER-approved merge plan.

## What is new (all opt-in or diagnostic; default = bitwise for valid production decks)
- **picres-1017 / picstop-1017** (rad_m1 Picard): `implicit_final_resid_log` (diagnostic), `implicit_conv_est_guard`,
  `implicit_conv_est_nofloor`, `implicit_final_resid_stop` (all read only when named). `nofloor` is NOT enabled in
  production yet (USER: hold until the long agcE1r comparison, TASK 8ba77933).
- **m1ctr-rollback-1017**: `<time>/counters_rollback` (default on, read only when named): event counters of discarded
  attempts (dt_retry, BE redo, helper solves) are rolled back. hst unchanged; only the event-log counter columns and
  one new end-of-run `counters_rollback:` line in run.log.
- **s3grid-1017**: opt-in `problem/he_ic_remap_file` path + tools/s3grid, tools/e1remap (remap_r.py).
- **cs-port-1016 / cs-vetgd-1017 / cs-s3b-1017 (+fix)**: cubed-sphere M1 + vet_gd (cs-only code; sp bitwise).
  `vet_gd_rotate_axis` (z default on sp, 3d on cs) is read only when named.
- **ic-standard-1015**: tools/ic_star (standard IC tool, gate (a') script `regress/gate_aprime.py`). Tools only.
- Regenerated dt_retry visitor.

## Gate (Raven A100, 4 ranks, 10 cycles; /raven/ptmp/jinma/mergegate_1017, job 31257166)
bsg / e1 / hg production decks: every hst/bin/rst md5 identical to bc5d24c1 (headers included), final rst bitwise,
same Picard means, same unread-key warnings; x2rat/x3rat refusal probes still FATAL. Not covered: cs decks, dt retries.

## What to do
- Running production: nothing; do not rebuild or restart running chains for this.
- Next production build from fin51-r-fix gets bb360d8f; expect one extra `counters_rollback:` line at end of run.log.
