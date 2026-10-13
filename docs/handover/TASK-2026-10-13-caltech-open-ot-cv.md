# TASK raven -> Caltech: open-ot CV, the IDORT gas-feedback set with `implicit_cs2_vimp = false` (165 CPU runs)

From Raven worker OT, via the coordinator. This is the regression check of a method candidate:
- the conservative2 vimp prediction (ImplicitCs2Vimp) leaves a per-step O(dt) term;
- the new key `implicit_cs2_vimp` (fresh default false) removes it.

## 1. Build (report md5)
- `fin3-1012` @ `281dd686`, then `git am patch/0001-*.patch` (one commit, 2 files in src/rad_m1; HEAD should become 37157508 or
  an equivalent hash).
- Plain CPU build, the same flags as caltech_ot:
  `-D Kokkos_ENABLE_OPENMP=On -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release`.

## 2. Run
`bash run_cv.sh <binary> 16` on one node (P = 16 parallel runs; list.txt gives the threads, 1-8).
- runs/cart:
  - eq1, eq2, eq3t, marsh c1/c10, rs_m1.2/m2/m3 (n512 / n1024), ddf, ddiff;
  - cloud and pipe c220/c2200 are long; run them last, or report that they are not done.
- runs/sp: every FIN3 sp case with gas feedback (rsw*, rswv*, rswt*, rswtv*, rlat*, rlatt*, relax*, marsh_*), X and _T ladders.
- Debug-QoS resubmits are fine: dirs with rc=0 are skipped.

## 3. Report: docs/handover/NOTE-2026-10-13-caltech-open-ot-cv.md
- STARTED (md5, job id), then the rc=0 count, any rc != 0 line and the wall times.
- Results:
  - `tar czf cv_runs.tgz runs` (tab/ bin/ run.log rc.txt time.txt; no rst);
  - push it under docs/handover/caltech_cv/.

## Coordinator notes
- Bundle path on the fork: docs/handover/open-1013/caltech_cv/ (list.txt, run_cv.sh, patch/, runs/).
- Priority: after od-vetpath reruns. Debug QoS only; long runs (cloud, pipe) that exceed 28 min: skip and list them.
- rs_m* decks: if they need m1_shock_ref, use docs/handover/caltech_od/shockref/ref_m{1.2,2,3}.txt.
- This is a METHOD candidate test only; production is NOT affected.

## FIX (re-pushed bundle): the first bundle FATALs on every run
- The old args.txt passed rad_m1/implicit_cs2_vimp on the command line, and AthenaK rejects keys missing from the deck. Now the key is
  inside each in.athinput and args.txt no longer carries it (165 dirs).
- If you already ran anything: delete those run dirs (or their rc.txt) and rerun from this re-pushed bundle.
