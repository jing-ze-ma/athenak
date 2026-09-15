---
name: dhj-general-eos-ck-blowup
description: OPEN - the general/tabulated-EOS blow-up. AUTHORITATIVE SOURCE IS NOW docs/HANDOFF_dhj_ck_eos.md IN THE REPO (commit 397b49b1); read it before proposing anything
metadata:
  type: project
---

**Read `docs/HANDOFF_dhj_ck_eos.md` first** (merged locally 2026-08-24, commit
`397b49b1`, written by the viper session). It lists NINE refuted hypotheses with
evidence, the exact reproducer input deltas, and the one lead left. This file only
records what is NOT in it.

## The two-blow-up split (the thing to internalise)

1. **ideal EOS + correlated-k** -- NaN at t ~ 5.2e3 s. **CLOSED**, fixed by the
   `rt_de_max` limiter (`ea823e3a`).
2. **general/tabulated EOS** -- still open. Dies under ANY RT: general+grey at 2.509
   rotations, general+correlated-k at 0.37-0.67 (1.07 with RKG). 1 rotation = 3.05e5 s.

Failure two is the open one, and the limiter does nothing for it.

## Corrections to what this file used to say

- The viper reproducer is **`bbot = 1.0e1`** and **`use_rkg_sts = false`** with
  `max_eta = 1e13`. Earlier notes here said `bbot = 3.0` and left RKG on -- both wrong.
- The viper meshblock was **`nx2 = nx3 = 16`** (32 blocks on one GPU), not 64.
- Viper's HEAD is `c1d513c2`, the same commit orion has, so **"viper's binary version"
  is not the difference** -- that last open candidate is closed.
- The initial state IS healthy (measured here and independently on viper: `ad_dump_file`
  shows zero genuinely super-adiabatic levels). Still true.

## What the orion session of 2026-08-24 adds

Everything run here used `bbot = 3.0` and `use_rkg_sts = true`, i.e. **not the
reproducer**, so its clean runs prove nothing about failure two:

- `ckrepro` (max_eta 1e14) reached t = 6.89e5 s = 2.26 rot and `ckrepro_eta13` 2.13 rot,
  both clean. Past the 0.37-0.67 rot correlated-k window but short of grey's 2.509.
- The `rt_de_max` on/off A/B was testing **failure one**, which is already fixed -- which
  is why disabling the limiter changed nothing. Trajectories agreed to 7 digits.
- Meshblock 8x8 vs 64x64: identical cycle lines, both clean. Decomposition is not
  physical here.

Genuinely new, and NOT in the handoff doc:
- **[[rt-srclim-warn-rank-local]]** -- the clip count is rank-local and the warning can
  be silent entirely. This matters for the doc's own numbers: any clip count in the
  campaign record is rank 0's share.
- **[[hst-cartesian-volume-on-spherical]]** -- `.hst` weights cells by `dx1*dx2*dx3`, so
  on this spherical mesh its `mass`/`tot-E` are NOT physical integrals. The doc quotes
  "mass to 0.15 %" and refuted hypothesis 9's "loses 67 % of its mass"; if those came
  from `.hst` they are the unweighted integral and should be redone from `.bin`
  ([[dhj-conservation-check]] has the recipe). 67 % is probably real regardless; 0.15 %
  is the kind of number this bug distorts (measured 8x here).

## 2026-08-24 orion: the doc's baseline does NOT reproduce on CPU

Job **190900**, `/orion/ptmp/jinma/Athenak/blowup/`. The exact input is now COMMITTED AND
PUSHED as `inputs/mhd/deep_hot_jupiter_rt_blowup.athinput` (commit `cd2e8815`, md5
`f52295ff8c321338369a91b76a5a1452`) so the viper session can run the identical file --
only `ck_table` / `ck_data_dir` are orion-absolute and need repointing. Built from
the shipped input with EXACTLY the handoff table's deltas -- verified line by line, every
"shipped" value in their column matches `397b49b1`, and `mesh/nx1` is 64). Limiter OFF
(`rt_de_max = -1.0`), so this is the variant the doc says **dies at 0.674 rot**.

**It passed 0.674 rot clean** -- t = 2.0616e5 = 0.6759 rot at cycle 155600, dt steady at
1.068, no NaN, no deceleration. dt fell 10 -> 1.21 by 0.23 rot then FLATTENED and drifted
back up; it does not look like a run heading into a blow-up.

The only remaining difference is the **platform**: viper ran on an APU/GPU, this is the
Freya/orion CPU build (16 ranks x 7 threads, gcc 13 + OpenMPI 4.1). If the failure is
robust across floors, resolution, resistivity, condensation and metallicity but not
across backends, it is likely seeded by numerical detail rather than being a robust
physical instability -- and the Gamma_1 lead would need re-checking on CPU first.

Throughput 14.8 simulated s per wall s, so `tlim = 1.2e6` (3.934 rot) needs a restart
chain past the 12 h slot; checkpoints every 1.0e4.

**Two open asks put to the user, not yet answered:** restart-chain to tlim, and ask the
viper session to re-run its own baseline at `397b49b1` to confirm the failure still
reproduces there (if it does not, the target moved rather than the platform mattering).

## The one lead left, per the doc

Gamma_1: flat 1.4728 for ideal, 1.084-1.666 for the table with jumps of 0.58 between
ADJACENT cells. The clean test is to evaluate Gamma_1 at the interface from the
reconstructed state rather than reconstructing the stored cell-centred value -- same
order, same dissipation. Donor cell is NOT a valid test (changes diffusivity, and
dissipation is already known only to delay). `src/mhd/mhd.hpp:588-596` is the related
code, and **the user said "let's not go that way" on 2026-08-23** -- do not open it
unless they raise it.
