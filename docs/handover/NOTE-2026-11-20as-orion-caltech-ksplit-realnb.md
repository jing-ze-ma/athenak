# NOTE 2026-11-20as orion -> caltech (USER): CONV3, KSPLIT diagnostic (new build ksplit) + one more LW bin script

1. New build **B_ksplit** = caltech-runs/ksplit, commit  (v20l plus a print only). Arms: `arms_ksplit.txt`: KS1_hg (6 cycles), then
   KS001_hg (5 cycles). Medium priority. Put the counts of 'CLEAN KSPLIT' / 'CLEAN OLDVEC' lines in the result line, and send run.log.gz.
2. LW bins: please run `python3 docs/handover/caltech_runs_1121/real_nb.py <run dir>` on LWr1170_hg and LWr260_hg (same usage as
   real_diag.py), and put the output in `caltech_results/LW_real_nb.txt`. If their bins are already deleted, run the two reruns LWs1170_hg /
   LWs260_hg from `arms_ksplit.txt` (B_v20f, KEEPBIN=1, ~10 min each), use them instead, then delete the bins. Low priority.

-- orion (for the USER)
