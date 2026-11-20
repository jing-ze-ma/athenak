# NOTE 2026-11-21d orion -> caltech (USER): 3 arms for the v20g gate, PRIORITY ahead of batch 2/3

New build **B_v20g** = branch caltech-runs/v20g, commit 6b48f352cbdff12bd2c9d9fd708767ebf1a8aa5f (v20f plus the cs face tie). Same cmake line as batch 1.
Arms: `docs/handover/caltech_runs_1121/arms_v20g.txt`. They have the same form as V20F_rsg / V20FP_rsg / V20F_hg (batch 1 P3):
- V20G_rsg: RSG restart check, 4 cycles. Leg a runs straight; leg b runs 1 cycle, writes a rst, restarts and finishes. They must be bitwise equal.
- V20GP_rsg: RSG with CLEAN_POISON_GHOST=1, 4 cycles. It must equal V20G_rsg leg a bitwise; please send both hst files.
- V20G_hg: HG restart check, 4 cycles.

Please start them as soon as nodes free up, **before any batch-2/3 arm**. This is the GPU gate of the v20 candidate. Results go to
`docs/handover/caltech_results/<TAG>/`, as before (res.txt shows the rt comparison).

-- orion (for the USER)
