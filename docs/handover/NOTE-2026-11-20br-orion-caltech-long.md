# NOTE 2026-11-20br orion -> caltech (USER): LONG robustness pairs on B_v201e (4 arms now + 4 later), fill the nodes

`arms_long.txt` holds the arms, the layout rules, and the expected wall times. The pairs are key-off vs key-on, for 1500 HG / 3000 RSG cycles.
- Start the HG pair (8 ranks on 2 nodes each: hgGH3b and rsgEH2) and the RSG pair (4 ranks, single nodes r1H1 / testH1) now.
- The L2_* second seeds come after them, and after any REALZ3 arms that arrive.
Send back per arm:
- res.txt and both hst files;
- run.log.gz (it may be tens of MB; if it is over 50 MB, grep CLEAN WARNING / CLEAN OLDVEC / retry / FATAL / cycle= lines into a file instead);
- the CW count, the retry count, and the last cycle reached.
Post a progress line every ~30 min with the current cycle of each arm.

-- orion (for the USER)
