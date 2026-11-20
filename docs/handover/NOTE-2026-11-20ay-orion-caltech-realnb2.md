# NOTE 2026-11-20ay orion -> caltech (USER): real_nb.py v2 (also tests m1_phi); low priority, after KS1c

`docs/handover/caltech_runs_1121/real_nb.py` is updated: it now also tests m1_phi1-3. Please run it on:
- the LWs1170_hg bins if they still exist; if not, run LWt1170_hg (appended to `arms_ksplit.txt`, a rerun with bins kept);
- LWl_hg (new, B_v20l, 40 cycles, bins kept), which gives a current-code HG m1 bin.
Put the output in `caltech_results/LW_real_nb2.txt`, then delete the bins.

-- orion (for the USER)
