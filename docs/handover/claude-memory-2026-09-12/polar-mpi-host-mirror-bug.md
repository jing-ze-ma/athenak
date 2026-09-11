---
name: polar-mpi-host-mirror-bug
description: "FIXED 3882e37f: every multi-rank run with use_polar_boundary aborted on cycle 1 between 2026-08-19 and 2026-08-22, because the polar EMF average's host mirrors were never sized"
metadata:
  node_type: memory
  type: project
---

`mhd_corner_e.cpp` and `mhd_corner_e_uct.cpp` size their polar-average work arrays once,
guarded by `inner_local.extent(0) != ncells1`. That guard NEVER FIRES, because `mhd.cpp`'s
constructor already pre-sizes `inner_local`/`outer_local` when `use_polar_boundary` is on.
So the two host mirrors in the same block, `polar_inner_h`/`polar_outer_h`, stayed at their
default extent 0.

Only the `#if MPI_PARALLEL_ENABLED && nranks > 1` branch touches the mirrors, so serial and
single-rank runs were fine. Anything polar on >1 rank aborted on the FIRST cycle with

    Kokkos::deep_copy extents of views don't match: (0) inner(68)

Introduced 2026-08-19 in **006bae07** ("Stop the polar EMF average running on 68 threads"),
fixed 2026-08-22 in **3882e37f** by testing the mirrors' own extent in the same guard.

**Why it matters beyond the fix:** it was LOUD, so nothing produced wrong answers -- but no
multi-rank dhj run could have worked in that window, and the two dhj inputs are the only
problems in the tree that set `use_polar_boundary`. Nothing in `tst/test_suite/` exercised
a polar boundary at all, let alone under MPI, which is why it survived. There is now
`tst/test_suite/rad/test_rad_dhj_ck_mpicpu.py`, which fails on 006bae07..3882e37f^.

**How to apply:** when a Kokkos View is sized lazily behind an `extent(0) != n` guard,
check that EVERY array in the block shares the condition -- a companion array pre-sized
elsewhere silently disables the whole guard. See [[correlated-k-design]].
