---
name: localise-by-dilution
description: A flat/uniform error profile is NOT evidence against a boundary origin -- vary the DOMAIN SIZE at fixed dx and dt and see whether the error dilutes as 1/L
metadata:
  type: feedback
---

When an error is suspected of coming from a boundary, the instinct is to profile it
against distance from that boundary and conclude "uniform => not the boundary". **That
inference is wrong**, and it cost most of a session in
[[cubed-sphere-mhd-convergence]]: an error generated at a boundary and then carried
through the volume by waves and advection looks perfectly flat, and I twice reported the
boundary was exonerated when it was the whole cause.

**The test that actually works: change the DOMAIN SIZE, holding dx and dt fixed.** A
boundary-generated error diluted over the volume falls as 1/L; an interior error does not
move. In that session, shell thickness L = 1 / 2 / 3 at fixed dr = 0.0625 and fixed
dt = 6.203e-3 gave a coefficient of **0.375 / 0.186 / 0.124** against 1/L =
0.375 / 0.1875 / 0.125 -- agreement to 1%, which no other hypothesis explains.

**Why:** holding dx fixed keeps the truncation error per cell fixed, and holding dt fixed
removes the temporal scaling, so the ONLY thing varying is how much interior there is to
dilute a fixed boundary flux over. It is the cleanest single discriminator available, and
it is cheap.

Two companions from the same hunt, both worth reusing:
* **Scan the physical parameter the error should be proportional to.** The rate deficit
  was constant to 1% over an 8x range in omega, which is the fingerprint of a ghost PHASE
  error (the ghost error is proportional to omega, so a FRACTIONAL rate deficit has the
  omega divided out). A dimensionless invariant is far more diagnostic than a magnitude.
* **Grep for the offending variable before theorising.** `dt` appears NOWHERE in
  `mhd_fluxes.cpp`, `src/reconstruct/*.hpp` or `hlld_mhd.hpp`, which proved in one command
  that a dt-dependent error could not be built in the interior at all.

See also [[measure-impact-before-claiming]] and [[confounded-tests-rejected]].
