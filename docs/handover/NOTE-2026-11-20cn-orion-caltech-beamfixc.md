# NOTE 2026-11-20cn orion -> caltech (USER): BEAMFIX combined candidate, HEAD OF BOTH QUEUES (100 cycles); and G_rsg partial log

## 1. The combined candidate (USER priority, as announced in NOTE cm)
New build **B_bfxc** = caltech-runs/beamfixc, commit 7d3d3c03608081843ca27b46dee9916d3737ba3d. It is the xcorner fix plus the bracket forms.
Four arms are appended to `arms_beamfix.txt`, 100 cycles each:
1. BFXC19_hg (NR=8) and BFXC19_rsg (NR=4): dg_form=19, xcorner=3.
2. Then BFXC13_hg and BFXC13_rsg: dg_form=13, xcorner=3.
They take the **next free nodes** in both pools, ahead of everything else, G_on_hg included. Do not preempt anything that is running.
After them, BFX16_hg and BFX15_hg (NOTE cl), then the NOTE ck order.

## 2. G_rsg (LZTEL3): please push the partial run.log now
Please push G_rsg's **current partial run.log.gz and both hst files** to caltech_results/G_rsg_partial/ now, while it runs, plus its info.txt line. Then push the full result when it ends.

-- orion (for the USER)
