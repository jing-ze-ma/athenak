---
name: branch-preexisting-breakage
description: "Two pre-existing breakages on the general-eos branch that are NOT EOS related: a non-MPI build failure and stale linwave test scripts"
metadata: 
  node_type: memory
  type: project
  originSessionId: 2cb708dd-0a5f-46c0-b912-dd406f996581
  modified: 2026-08-13T14:08:35.461Z
---

Found 2026-08-13 while running `tst/run_tests.py`. Neither is caused by the general-EOS
work; both predate it and both are trivially fixable if the user wants them fixed.

1. **The code does not build without `-DAthena_ENABLE_MPI=ON`.**
   `src/bvals/physics/bfield_bcs.cpp:451-456` calls `MPI_Allreduce`/`MPI_DOUBLE`/
   `MPI_SUM`/`MPI_COMM_WORLD` outside any `#if MPI_PARALLEL_ENABLED` guard. So the plain
   serial build in CLAUDE.md's quick-start is broken, and the regression harness must be
   given `--cmake=-DAthena_ENABLE_MPI=ON` on top of the usual flag set.

2. **The linwave test scripts passed a stale `problem/vflow`** where the pgen has long read
   `vx0`, so both suites aborted in `run()` (AthenaK exits fatally on a command-line
   override of a parameter the input file does not define). FIXED on the branch in
   `30cb3cbc`; **`main` still has it**, so upstream's CI is not running these suites either.

   `hydro_linwave` passes after the fix. **`mhd_linwave` still fails, and it is a stale
   THRESHOLD, not a code error** — verified 2026-08-13 by running the same test in a
   worktree of `main`, which produces the identical six warnings with bit-identical
   numbers (0.295111 for llf/hlle, 0.296246 for hlld, against `conv_threshold = 0.29`).
   Only rk2 + wenoz fast waves, 6 of 168 checks. The cause is that the test measures
   convergence on its COARSEST pair, nx1 = 16 -> 32, which is not asymptotic for WENO-Z;
   continuing the same run to 64 and 128 gives ratios 0.2963, 0.2553, 0.2503, i.e. clean
   second order once resolved. Raising the threshold to ~0.30, or measuring 32 -> 64,
   would fix it -- but that is an upstream decision, not a branch one.

See also [[fofc-1d-segfault]] and [[freya-build-procedure]] (the harness recipe, which now
needs the MPI flag added).
