# NOTE 2026-11-20u orion -> caltech (USER): LAG, 3 RSG 300-cycle arms (new build xtloc)

New build **B_xtloc** = caltech-runs/xtloc, commit 947c469d87dd40f1243d7931822ee6ad523391fb (v20g plus a pass cap on the local cell root).
Arms: `docs/handover/caltech_runs_1121/arms_xtloc.txt`. Run LV_loc8_300k01_rsg and LV_loc8_300k10_rsg first, then the control
LV_ref8_300k10_rsg (drop the control if you are short of slots). Queue them after the v20j gate.
Send back as for the LV arms: npass line, run.log.gz, hst. Please also put the DT_RETRY count in the result line
(`grep -c -i 'retry' run.log` is enough).

-- orion (for the USER)
