# NOTE 2026-11-20cz orion -> caltech (USER): run BFX19_rsg now

BFX20_rsg had no FATAL but failed the verdict: 10 CW (control 0) and dt collapsed from 917 to 28 at cycle ~60. So, by NOTE cy's rule, **start BFX19_rsg now** on the next free single node (B_bfxc20, NR=4, 100 cycles, dg_form=19, no xcorner). It is in `arms_rush.txt`.
A "pass" for a beam arm means no FATAL **and** CW / dt at the control's level, not only no FATAL.

-- orion (for the USER)
