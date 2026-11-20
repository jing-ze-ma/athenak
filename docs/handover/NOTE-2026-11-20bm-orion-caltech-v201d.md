# NOTE 2026-11-20bm orion -> caltech (USER): v20.1d gate, 5 arms (new build v201d), HIGH priority (after the V201 rt arms)

New build **B_v201d** = caltech-runs/v201d, commit 774262b563613f59028d289edbf044474628e642. It is v201 plus the kern-fw speed commits, a cubed-sphere seam-geometry fix, and the xlim
key (off). Plain cmake line. Arms are in `arms_v201d.txt`. Xchecks to put in the result lines:
- V201D_rsg/a and V201DP_rsg hst vs V201_rsg/a: equal expected; if not, give the first differing row and max relative difference.
- V201D_hg/a vs V201_hg/a: may differ (seam fix); give the first differing row and max relative difference.
- licmp LI4d vs LI8d into `caltech_results/LI_cmp_v201d.txt`.
- s/cycle next to V201_*.

-- orion (for the USER)
