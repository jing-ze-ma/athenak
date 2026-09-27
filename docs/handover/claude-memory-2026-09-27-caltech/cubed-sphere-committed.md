---
name: cubed-sphere-committed
description: The whole cubed-sphere hydro effort is COMMITTED as 9492a946 -- five sessions of work is no longer at risk in the working tree
metadata:
  type: project
---

Committed 2026-08-27 (session c271a0e2) on branch `polar-average-perf` as **9492a946**,
on top of 1e19d4e7. Everything from [[cubed-sphere-hydro-state]],
[[cubed-sphere-x1-radial]], [[cubed-sphere-panel-frames]], [[cubed-sphere-seam-basis]]
and [[cubed-sphere-seam-interp]] is in that one commit -- it had been uncommitted across
five sessions.

**Re-verified immediately before committing**, all bit-identical to the recorded
references:
* Test A uniform static: mass 24.000000, totE 36.000000, KE ~1e-30, dt 8.78109e-03
* rigid rotation, joint refinement: 1-KE 5.57821e-08 / 2.39364e-09
* linear_wave_hydro: L1 errors bit-identical to the pre-change reference

**Also fixed as part of the commit: the two tracked input files were STALE.**
`inputs/tests/cubed_sphere_{uniform,rigidrot}.athinput` were still written for the OLD
x3-radial convention and carried the retracted "STATUS: FAILING" text. Both are now
x1-radial, match the tested configurations, and were run to confirm they reproduce the
documented numbers. If you re-run from a scratchpad template, **check `ix1_bc`** --
the rigid rotation needs `user` (the exact-state BC), and `reflect` silently gives a
250x worse 1-KE that looks like a regression.

**Style**: the repo-wide `check_athena_cpp_style.sh` is unusable here -- run from
`tst/test_suite/style/` it misses the root `CPPLINT.cfg` and reports 2305 errors on
untouched files. Check instead by diffing violation COUNTS per file against HEAD; every
changed file came out <= its HEAD count, and the two new files are clean.

**What is left, unchanged**: `bvals_fc.cpp` (MHD) needs the basis transform + the seam
resample + seam EMF consistency ([[cubed-sphere-for-hot-jupiter]] item 1); the 2D nx3=1
shell; corner/edge halo buffers are still plain copies.

**NOT PUSHED.** Push to the fork, never upstream ([[use-fork-not-origin]]).
