---
name: cs-wb-source-cached
description: The well-balanced gnomonic source (cs_wellbalanced_src) now caches its geometry per (m,k,j) -- 88064f67; its cost fell from 46 % to 4 % of a run, bitwise identical, restart-clean. Also records the dhj restart-test trap (tlim clipping) and the float32 dump ulp
metadata:
  type: project
---

**88064f67 (2026-09-05):** `Coordinates::BuildWBGeometry` stores the cell triad and the
four tangential faces' normals/triads (60 reals per (m,k,j)) in `wb_geom`; the kernel reads
them back with the arithmetic order unchanged. Rebuilt on every call under AMR or if nmb
changes; otherwise once. The old "1.68x wall clock until its face triads are cached" note in
[[cs-mhd-instability-characterized]] / [[cs-dhj-production-retry]] is now obsolete.

    cubed_sphere_mhd_strat, 32^3/panel, serial CPU, 160 cycles:
      old kernel 186.7 s   cached 133.1 s   source OFF 128.1 s   -> 4 % instead of 46 %
    five bin dumps + hst: BITWISE identical old vs cached.

**Restart proof, and TWO traps found on the way:**
1. Ending the first half at `time/tlim` CLIPS its last step to land on tlim, so a restart
   from that file follows a different step sequence than the straight run -- O(1) velocity
   differences after 3 cycles that have NOTHING to do with the code. Restart from a file the
   STRAIGHT run wrote at a real cycle boundary (its own rst output).
2. Done that way, all fields are bitwise identical except a 1.5e-8 in bcc2 (7.5e-9 in velz/
   bcc2 with the source off): the bin dumps are FLOAT32, so that is one ulp at 0.25 -- an
   underlying round-off difference, identical on the pre-cache binary. Not the cache.
3. `cs_test` (build_cs) CANNOT restart: UserProblem returns before enrolling its BCs and its
   file-static parameters are never re-read (src/pgen/cs_test.cpp:282). Pre-existing; use the
   dhj generator for any restart test on the cubed sphere.

Test dirs: scratchpad wb_speed/ and wb_restart/ (session 3ac6262f).
