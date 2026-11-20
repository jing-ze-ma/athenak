# NOTE 2026-11-20dg orion -> caltech (USER): thinf v7 smooth xdonor cap (XCP), 3 arms with the TG arms

New build **B_xcap** = caltech-runs/xcap, commit 999d7afaaefb510a55a9b879c7d827687bf0de14 (thinf v7 plus the key rad_m1/clean_stag_xcap_p, read with GetOrAdd; 0 = the old hard cap, bitwise).
Three arms (in `arms_rush.txt`), **1 node each, concurrent with the TG4/TGinf arms**, at the head of the queue:
- XCP8_hg and XCP2_hg: as THF100_hg (hg_thv7.keys, hegiant.00139, NR=4, 100 cycles), with xcap_p = 8 and 2.
- XCP8_rsg: as THF300_rsg (rsg_thv7.keys + clean_conv_kappa=0.1, rsg.01295) but 100 cycles, with xcap_p = 8.
If the key is rejected on the command line, add `clean_stag_xcap_p = P` to <rad_m1> of the deck.
In each line please give, as for THF100_hg:
- CW count, retries, dt at cycle 100, steady s/cycle;
- the histogram of "solves" per stage (CLEAN cycle lines);
- the XDONOR donor and cap counts (the cap count should be 0 with p > 0);
- the full run.log.gz.

-- orion (for the USER)
