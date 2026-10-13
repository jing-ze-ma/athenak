# TASK raven -> Caltech: open-ot CVALL, `implicit_cs2_vimp = false` across the FULL IDORT suite (359 run dirs)

From Raven worker OT, via the coordinator. USER: "try the ot fix with all tests".
Bubbles and 2-source are excluded: they are outside the suite.

## 1. Binaries
- OLD = your plain CPU build of `fin3-1012` @ 281dd686 (caltech_ot; md5 9ae9f09a...).
- NEW = 281dd686 + `git am patch/0001-*.patch`, built with the same flags (the CV bundle's binary is fine).
- Report both md5.

## 2. Run
`bash run_cvall.sh <OLD> <NEW> 16 10` on one node; debug-QoS resubmits are fine, finished steps are skipped.

Per dir there are two steps:
1. BIT:
   - OLD and NEW each run 10 cycles with the same threads. NEW names the key in the deck.
   - The final rst data are compared and the result goes to bit.txt (BITWISE / DIFF / FAIL).
   - This is the per-run proof that the key is inactive.
2. FULL:
   - Only for DIFF dirs that are not already in caltech_cv (full_skip.txt).
   - NEW runs to the end into <dir>/full. patm/patmT are two-stage: stage A, then -r for stage B.

Dirs:
- cart: 45 Cartesian tests;
- sp: 194 sp runs (all trees, incl. patm/patmT, nlte2, pinc, shd, core, thin, dif);
- ss: 4 sshell runs;
- p2: 10 P2 runs;
- lad: 106 Cartesian X/T/C/S ladders.

## 3. Report: docs/handover/NOTE-2026-10-13-caltech-open-ot-cvall.md
- md5, job ids;
- bit_all.txt (written at the end; the counts per status);
- every FAIL line;
- `tar czf cvall_runs.tgz runs --exclude='*.rst'`, pushed under docs/handover/caltech_cvall/. Split it if it is > 90 MB.

## Coordinator notes
- Bundle on the fork: docs/handover/open-1013/caltech_cvall/. Priority: after the CV rerun (d57b1ed2) and the od-vetpath shock reruns.
- This is large (359 dirs). If the debug QoS is too slow, you MAY use the L40S node's CPU cores (as you offered for cartlong) for this TASK.
- Post interim counts (bit_all so far) every few hours. METHOD candidate test only; production is NOT affected.

## UPDATE (runner fix): re-pushed run_cvall.sh + shockref/
- Fixes two bundle issues: decks that already carry an rst block ("more than one restart output"), and rs_m* shock refs at Raven paths
  (now in caltech_cvall/shockref/).
- Use the new run_cvall.sh. For any dir whose bit.txt says FAIL, delete bit.txt and rerun; BITWISE / DIFF dirs stay.

## UPDATE: STOP the CVALL FULL step (superseded by NPALL, TASK-2026-10-13-caltech-open-ot-npall.md)
- The key-off variant is rejected. Finish/keep only the BIT step counts (useful cross-check); do not start new FULL runs.
