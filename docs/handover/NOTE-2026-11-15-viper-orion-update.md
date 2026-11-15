# NOTE 2026-11-15 viper -> Orion: update after the clean-scheme handover

- Restyle checks: the user moved the restyle CPU bitwise comparisons (branch restyle-1114, merged per module into
  clean-1114) to Raven CPU nodes. The earlier plan note said they would come to you as A/B TASKs; they will NOT,
  unless a TASK-*-orion-* file for them appears on this branch. Your focus is the clean scheme only.
- Clean scheme: unchanged — you own fork/clean-cc (head at handover dc67dba7b); see
  HANDOVER-2026-11-15-clean-scheme-viper-to-orion.md. Restyle modules 1-2 are in clean-1114 (c9d1a898e); module 3
  (well-balanced: wb_x2, wb_seam, wb_rmin, wb_poison_test removed; cs format cleanup) follows; merge clean-1114 into
  clean-cc when convenient.
- GPU checks: viper's apudev partition (MI300A, 15-min jobs) is available now; ask in a NOTE when you want a GPU
  build/validation of clean-cc (e.g. a commit + a run list), and viper runs it there.
- Please acknowledge with a short NOTE when you have read the handover, and send results as
  NOTE-<date>-orion-viper-<topic>.md.
