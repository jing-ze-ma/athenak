# NOTE 2026-11-20cv orion -> caltech (USER): RUSH, run all verdict arms at once on separate single nodes

The USER wants the fixes landed within about 30 min, with maximum parallelism. HG fits 4 ranks on one H200 node, as in batch 1.
**Supersedes the layouts and order in NOTEs ct / cu / cr for these arms.** The list is in `docs/handover/caltech_runs_1121/arms_rush.txt`, all NR=4 on ONE node each:
1. **Now, concurrently** (6 single nodes): BFX20_hg, BFX19_hg, BFX20_rsg, BFX19_rsg (100 cycles); G_on_hg and G_off_hg (300 cycles, B_lztel3 plain, no env).
2. Then, on the next free single nodes: BFXC20x2_hg / _rsg, then BFXC19x2_hg / _rsg (100 cycles).
- **USER: pause lower-priority arms (LX_*, D_off_hg, LN2_*, V201Dk_*, L2_*) to free nodes.** They can be restarted later. Let the beam and G arms already running finish.
- G_off_hg at NR=4 is not bitwise-comparable to L1_hg_off (NR=8), so no bitwise verdict is needed. The speed pair is G_on_hg vs G_off_hg on the same layout. If possible, put the two on identical node types.
- Please post each line **as soon as each arm ends** (one commit per arm is fine).

-- orion (for the USER)
