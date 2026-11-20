# NOTE 2026-11-20ct orion -> caltech (USER): corrected BFXC arms (xcorner=2), HEAD of both queues

The arms that used xcorner=3 are commented out. The replacements all run on **B_bfxc20** (it holds forms 19 and 20) with **xcorner=2**, 100 cycles each, appended to `arms_beamfix.txt`:
1. BFXC20x2_hg (NR=8) and BFXC20x2_rsg (NR=4): dg_form=20.
2. BFXC19x2_hg and BFXC19x2_rsg: dg_form=19.
They take the **next free nodes** in both pools, ahead of everything queued. Do not preempt running arms, so if G_on_hg / G_off_hg have started, let them run.
Behind them, the HG queue is G_on_hg / G_off_hg (if not started), BFX16_hg, BFX15_hg, then the rest.

-- orion (for the USER)
