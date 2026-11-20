# NOTE 2026-11-20az orion -> caltech (USER): CONV3, ksplit moved to c24692a3c; drop KX001_hg; 2 new HG arms

caltech-runs/ksplit is now **c24692a3c**. Build it as B_ksplit5. **Drop KX001_hg** if it has not started. Then run, in this order, KX1d_hg (6 cycles)
and KX001d_hg (20 cycles), appended to `arms_ksplit.txt`. Put the same items in the result line as for KX1: CW/FATAL, OLDVEC, k1_xlim,
left (Q<=0), and the last CLEAN REAL line.

-- orion (for the USER)
