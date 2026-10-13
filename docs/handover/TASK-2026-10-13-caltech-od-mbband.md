# TASK raven -> Caltech: OD mbband-1013, the FULL IDORT suite in TWO arms on ONE binary (359 run dirs)

From Raven worker OD, via the coordinator. USER: "fix the vet_sc multi-block band bug and run ALL IDORT tests with the fix",
then "run the IDORT suite in TWO arms on the mbband binary: (a) cells OFF, (b) cells ON". Bubbles and 2-source are excluded.

## 1. Binaries (report both md5)
- OLD = your plain CPU build of `fin3-1012` @ 281dd686 (the caltech_ot / CVALL OLD; md5 9ae9f09a...).
- NEW = branch **`mbband-1013` @ `d98f8868`** (fork; = 9972a56d + the band fix 60882a63 + the vet_path cells commits), built with
  EXACTLY the same flags as OLD (same compiler, OpenMP, MPI, Release, and the same Athena_M1_VGD_NMAX if OLD set one).

## 2. Run
`bash run_mbband.sh <OLD> <NEW> 16 10` on one node (debug-QoS resubmits are fine; finished steps are skipped).
Bundle: docs/handover/open-1013/caltech_mbband/ (decks = the CVALL decks: 45 cart + 194 sp + 4 ss + 10 p2 + 106 lad).

Per dir:
1. BIT: OLD, NEW (a: deck as is, vet_path default none), NEW (b: + rad_m1/vet_path=cells), 10 cycles each, same threads.
   bit.txt = "<OLD vs a> <OLD vs b> <dir> <s OLD> <s a> <s b>". (a) must be BITWISE in every dir (all IDORT runs are
   single-MeshBlock; the band fix only acts in the multi-block sweep). Any (a) DIFF: report it at once.
2. FULL (b): NEW + vet_path=cells to the end in <dir>/full_b, for EVERY dir (cells can switch on after 10 cycles).
   patm/patmT and the other two-stage dirs: stage A, then -r for stage B (two.txt), as in CVALL.

## 3. Report: docs/handover/NOTE-2026-10-13-caltech-od-mbband.md
- md5, job ids; bit_all.txt with the per-arm counts; every FAIL / non-zero rc line; the FULL_b lines (rc, time, vet_gd
  subdivided-step counter where printed).
- `tar czf mbband_runs.tgz runs --exclude='*.rst' --exclude='bitO' --exclude='bitA' --exclude='bitB'`, pushed under
  docs/handover/caltech_mbband/ (split if > 90 MB). Post interim bit_all counts every few hours.

## Coordinator notes
- METHOD candidate test only; production is NOT changed by (a) (he_star runs never use the vet_sc multi-block sweep).
- Priority: after the running od / ot tasks. If the debug QoS is too slow, the L40S node's CPU cores are fine (as for CVALL).
