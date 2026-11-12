# TASK 2026-11-12 viper -> Orion: dt_retry CPU gate, addendum (after NOTE de77c1a92)

## Why
Orion's NOTE found CN vs CF (forced retry, views minimal) DIFFER on CPU, while views=all was bitwise. Viper saw the same
on GPU on the R136 production restart: minimal != straight, views=all == straight. Bisecting with
time/dt_retry_views_add there found **m1.vgd_hr** (the vet_gd half-range ratios; vra_a1 is the same View). Restoring it
alone made the retry bitwise. A cycle that does not rebuild vgd_hr reads the last build's, so after a retry it read
the failed attempt's. Fix: **fork/dtretry-fix-1111 34d198efa** adds m1.vgd_hr to the minimal list.

## Build
As in TASK-2026-11-12-orion-dtretry.md (build_hf recipe), at **34d198efa** (incremental in wt_dtr_fix:
`git -C wt_dtr_fix fetch fork && git -C wt_dtr_fix checkout --detach fork/dtretry-fix-1111`, then make).

## Run (docs/handover/scripts/dtretry_1111/run2.sh, same deck, keys and tables as before)
```
cd <run area>/dtretry_1111 (fresh copy of the scripts dir; RUNS=$PWD/runs2)
RUNS=$PWD/runs2 FIX=<wt_dtr_fix>/build_hf/src/athena TB=<repo>/docs/handover/hgpos-1109/files ./run2.sh > gate2.out 2>&1
```
This runs 10 arms concurrently, 3 ranks each (30 cores), about 6 min. CN/CF/CTN/CTF repeat the failing pairs with
the fix. V1..V6 are CF plus one group of candidate Views added to the minimal list (dt_retry_views_add): if CF still
differs, the group whose V arm equals CN contains the View.

## Return (NOTE-2026-11-12-orion-viper-dtretry-addendum.md on this branch)
Paste gate2.out whole. Commit the NOTE only.
