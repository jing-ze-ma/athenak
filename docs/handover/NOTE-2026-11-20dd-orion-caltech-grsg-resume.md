# NOTE 2026-11-20dd orion -> caltech (USER): fill idle nodes, with the G RSG pair first, then resume the paused arms

No beam arms are queued now. New beam builds (THINF v7) will come soon and go to the **head** of the queue the moment they arrive.
Meanwhile, fill the idle nodes in this order:
1. **G_on_rsg / G_off_rsg** (`arms_rush.txt`): plain B_lztel3, no env, NR=4, 600 cycles, concurrent on the same node type. This is the clean RSG speed and CW pair for the L_z decision. G_off_rsg's hst should equal L1_rsg_off rows 0-599; please state whether it does. Push a p300 partial of each, as for G_*_hg.
2. V202F_hg (v20.2 gate), if not done.
3. **Resume the paused arms:** LX_rsg (from its last dump if one exists, otherwise rerun), the L2_hg pair, LN2_* / V201Dk_*, then LX_hg, D_off_hg and the L2_rsg pair.
Keep the 2-node reservation from 11-21 13:00 PST.

-- orion (for the USER)
