---
name: solar-convection-binary-provenance
description: "The pinned ideal192/table192 binaries CANNOT be rebuilt from git -- they came from uncommitted working state"
metadata:
  type: project
---

Found 2026-09-04 while trying to build a binary consistent with the `ideal192` spin-up.

`soleos/ideal192/athena` (md5 `9546a513`) was built **2026-08-14 16:30**. Every commit that
touches `src/pgen/solar_convection.cpp` is dated **2026-08-15 or later** (`88ace084`,
`6ccac245`, `1f0cbc9a`, `3c1b397b`, `3d1640e8`, `bfa566c2`, `cabf4922`, `d1c289dd`). So that
binary was built from **uncommitted working-tree state that no longer exists anywhere**.
The same is almost certainly true of `table192`.

**Consequence:** you cannot reproduce those runs bit-for-bit, and you cannot build a
"matching" binary to continue them from a restart.

**What to do instead** -- do not assume, measure. Build from HEAD and run the SAME fresh
input for a few hundred cycles, then compare `dt` against the old run's log cycle by
cycle. Done 2026-09-04 (`sunmovie/bincmp` vs `sunmovie/spinup`): identical base-state line,
and **bit-for-bit identical dt through cycle 200**, diverging only at round-off after
(3.500043e-01 vs 3.500046e-01 at cycle 300). That was enough to justify continuing the
spin-up state with the rebuilt binary instead of redoing 3h46 of spin-up.

Going forward the movie chain is pinned at `sunmovie/athena_fixed`, built from
`fix/solar-convection-restart-gravity` (commit `ef9561e2`) -- recorded by commit, not by md5.
Build it out-of-source in `sunmovie/build`; do NOT touch `athenak/build`, which is
configured for `PROBLEM=deep_hot_jupiter_rt`.

See [[solar-convection-restart-gravity-bug]], [[session-state-2026-09-04]].
