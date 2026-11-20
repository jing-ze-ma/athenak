# NOTE 2026-11-20cp orion -> caltech (USER): BFXC20 replaces BFXC13

New build **B_bfxc20** = caltech-runs/beamfixc20, commit 29f989379a6530883fd46e95c1838aa3b921f36f. It is B_bfxc plus one bracket-only commit, adding form 20.
- **Drop BFXC13_hg and BFXC13_rsg** if not started; form 13 already fails on HG.
- **New arms** (100 cycles, appended to `arms_beamfix.txt`): BFXC20_hg (NR=8) and BFXC20_rsg (NR=4), dg_form=20, xcorner=3.
- **HG queue:** BFXC19_hg, G_on_hg (NOTE co), BFXC20_hg, BFX16_hg, BFX15_hg, then the rest.
- **Single-node queue:** BFXC19_rsg, BFXC20_rsg, then the rest.

-- orion (for the USER)
