# NOTE 2026-11-20bw orion -> caltech (USER): fill every free node until the holders end. Priority list (arms_long.txt)

The USER asks to keep every node busy until 17:46 PST. New arms are appended to `arms_long.txt`, on B_v201e, in this priority:
1. R_fp_rsg / R_tel_rsg / R_xl_rsg: 1000 cycles, NR=4 on single nodes. Start them on free single nodes now (02-04 / 02-16 as RZ3m finishes; r1H1 / testH1 when the L1 RSG pair ends).
2. K_on_hg (NR=8, 2 nodes; after H_xl) and K_on_rsg (NR=4).
3. The L2_* seeds.
REALZ3 / LEAN2 arms, if posted, go ahead of 2-3. There are no new key-off arms: the L1_*_off runs are the same binary and layout, so they serve as
the reference. Do not start an arm that cannot finish before the holder ends; RSG runs at ~2.7 s/cycle on 4 ranks, HG ~3 s/cycle on 8 ranks.

-- orion (for the USER)
