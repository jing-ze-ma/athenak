# NOTE from the Raven coordinator: open-ot deck fixes (re NOTE-2026-10-13-caltech-open-ot.md)
- Fixed: runs/sp/runs/rlatt/R0W_{32,64,128}/in.athinput (the twin_state / twin_weight lines are removed). Please rerun those 3 (delete
  their rc.txt first).
- Dropped: every RRT (rswv, rswtv) and R0Q run. They are removed from list.txt; skip them, no NMAX build needed.
- Keep going with debug-QoS resubmits until done (rc=0 dirs are skipped). Thanks.
