# NOTE 2026-11-21eu orion -> caltech (USER): 2 seam-artifact arms AFTER the F_* / MB smokes (before ST / STD)

This answers a USER question on LQ10 seam artifacts. None of the existing smokes kept bins (armg.sh deletes them), so these two arms provide the data. The files are in 102d22cc9 under docs/handover/fresh_prod_1120/freshgo/ (copy the GL deck into fresh_prod_1120; bsg62m_wb_smoke you already have). Run on **B_v203**, 1 node each, NR=4, with your armg.sh copy, **RT=0 KEEPBIN=1**, 60 cycles, bins at cycle 60 only:
1. `RT=0 KEEPBIN=1 armg.sh SQ_lq $B_v203 decks/bsg62m_wb_smoke.athinput 60 output2/dcycle=60 output3/dcycle=60` (LQ10 fixed)
2. `RT=0 KEEPBIN=1 armg.sh SQ_gl $B_v203 decks/bsg62m_wbgl_smoke.athinput 60 output2/dcycle=60 output3/dcycle=60` (GL 7x18 fixed)
After each arm:
`python3 docs/handover/fresh_prod_1120/freshgo/seam.py <OUT>/<TAG>/a/bin 3.68e12 <athenak tree>/vis/python > caltech_results/<TAG>/seam.txt`
Then commit seam.txt (a few lines) with res / hst as usual, and **delete the bins**.

Order on free nodes: MB_bsg62m + F_* first (NOTE es / et), then SQ_lq / SQ_gl, then the v20.4 gate and ST / STD. Production preempts all of them on the go NOTE.

-- orion (for the USER)
