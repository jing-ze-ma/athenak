# NOTE 2026-11-20cb orion -> caltech (USER): BEAMFIX arms, HEAD OF THE QUEUE

New build **B_bfx** = caltech-runs/beamfix, commit 3a67711b9ec01fcd8afeb7e2fbd29765f890ce6a (he_star_m1, H200). Its src is v201e plus one key, rad_m1/clean_stag_dg_form; the default 0 is v20.1 bitwise.
Arms are in `docs/handover/caltech_runs_1121/arms_beamfix.txt`: BFX4/BFX0 for RSG (NR=4, 1 node) and HG (NR=8, 2 nodes, the L1_hg layout), 300 cycles each. The key is given on the command line (no deck change needed).
Run them **first** when slots free up. Do not preempt anything that is running. If 2-node slots are slow to free, start the RSG pair first.
Provenance: BFX0_* hst should equal L1_*_off rows 0-299 bitwise. Please say whether it does.
Results go to caltech_results/<TAG>/ as usual (info, res, hst, run.log.gz). We need the CW kinds, the dt_retry lines and the "energy created" E-floor lines from the log.

-- orion (for the USER)
