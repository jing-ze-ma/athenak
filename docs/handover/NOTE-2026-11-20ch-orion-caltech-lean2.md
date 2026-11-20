# NOTE 2026-11-20ch orion -> caltech (USER): LEAN2 gate rows (3 short arms)

New build **B_lean2** = caltech-runs/lean2, commit ac4b62c5574438bf734a7cf738820f745c6e2a44. It is the v20.1 source (v201e) with the old solver, cell layout and FV rays removed; use the usual he_star_m1 H200 recipe.
The arms are in `docs/handover/caltech_runs_1121/arms_lean2.txt`: LN2_rsg (rt), LN2P_rsg (CLEAN_POISON_GHOST=2) and LN2_hg (rt, ompio). They use the same recipe as LN_* / V201D_*, from rsg.01295 / hegiant.00139.
**Priority:** after BEAMFIX, LF/LX and LZTEL3, before the L2 seeds.
Please report the following, as for LN_*:
- hst and final-bin verdicts vs V201D_rsg/a and V201D_hg/a (bin headers may differ only in the parameter dump, since keys were removed);
- poison vs plain;
- CW / FATAL counts;
- steady s/cycle and GPU MiB peak, side by side with the base.

-- orion (for the USER)
