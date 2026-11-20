# NOTE 2026-11-20x orion -> caltech (USER): speed A/B, 8 RSG arms of 10 cycles (no new builds)

v20j was ~3 % slower than v20g on RSG (2.946 vs 2.858 s/cycle). This splits the two changes in it. Builds you already have:
B_v20g (6b48f352c), B_h200 (a774c43b3, the H200 speed commits only), B_lzbook (fb7a8a27f, hst budget columns only),
B_v20j (6b5388315). Arms: `docs/handover/caltech_runs_1121/arms_spg.txt` (SPG_{g,p5,lz,j}_{1,2}; 10 cycles, romio, 4 ranks).
- Round 1: the four `_1` arms start together, 90 s apart, on 4 nodes of the same type (not 02-16).
- Round 2: the four `_2` arms on the same 4 nodes, rotated, so each build runs on two different nodes.
Per arm, report: steady s/cycle, node, GPU peak, and whether user hst cols 1-26 + hydro hst == SPG_g_1 (bitwise expected).
Priority: after the v20j gate and LZCS, ahead of batch 2/3.

-- orion (for the USER)
