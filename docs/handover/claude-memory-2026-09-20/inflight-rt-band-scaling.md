---
name: inflight-rt-band-scaling
description: "2026-08-21: blocked-band RT scaling COMPLETE. RT_NB=4 is optimal; chain cost is dead linear at 12.3 ms/chain/100cyc without table lookup, 23.3 with it. 88-chain correlated-k would cost ~2100 ms/100cyc (RT ~60% of wall). Decision on whether to adopt correlated-k still open."
metadata:
  node_type: memory
  type: project
  modified: 2026-08-21T17:10:00.000Z
---

**Task:** measure what it costs to expand the deep-hot-Jupiter RT from picket fence (4
chains) toward a SPARC/MITgcm-style correlated-k scheme (11 bands x 8 g-points = 88
chains), before committing to a scheme. Follows [[inflight-rt-kernel-optimization]].

## STATUS: measurement done. Open question is the science decision, not the numbers.

The harness is still UNCOMMITTED in `src/pgen/deep_hot_jupiter_rt.cpp` (also saved as
`bench/polar_ab/rt_bandharness.patch`). Decide whether to commit it, drop it, or keep it
as a patch. `RT_NB` in the pgen is restored to 4 and `build-gpu-bench` holds that build.

## How the harness works

`problem/rt_nchain` (default 4 = production) and `problem/rt_ktab` (default false).
Chains past the first four are synthetic: same lookup, same two recurrences, same flux
accumulation, but their weight comes from a device array `rt_wgt` that holds zero, so the
sum is unchanged and the run stays on the production trajectory. **Verified inert in all
three builds: 88 chains + table lookup is bitwise identical (data) to 4 chains, and the
4-chain harness build is bitwise identical to the committed `fixed200.a` reference.**

`RT_NB` is a compile-time block size (`#define` at the top of the pgen): the sweep steps
`RT_NB` chains together, so private `I_ir_down_c` stays at RT_NB columns no matter how
many chains are requested. A partial block costs a full block (visible as nb16 charging
the same for 4, 8, and 16 chains).

Two harness traps already paid for:
- AthenaK command-line overrides only MODIFY existing parameters, so `rt_nchain`/`rt_ktab`
  must exist in `bench/polar_ab/correctness.athinput` or the run fatals.
- The input dump is embedded in the binary output header, so `cmp` against `fixed200.a` is
  meaningless once a parameter is added. Compare DATA only, with `bench/polar_ab/bindiff.py`.

## MEASURED (jobs 10978050 nb4, 10978077 nb8, 10978078 nb16). RT ms per 100 cycles.

| chains | nb4 | nb4+table | nb8 | nb8+table | nb16 | nb16+table |
|---|---|---|---|---|---|---|
| 4 | 124.8 | 137.1 | 224.7 | 276.5 | 591.2 | 680.2 |
| 8 | 173.5 | 229.0 | 225.2 | 275.9 | 589.9 | 681.6 |
| 16 | 270.7 | 415.3 | 354.1 | 485.5 | 588.2 | 684.5 |
| 32 | 468.8 | 787.4 | 614.3 | 906.5 | 1072.7 | 1285.6 |
| 44 | 614.7 | 1068.1 | 878.3 | 1323.5 | 1552.2 | 1879.5 |
| 88 | 1161.3 | 2096.2 | 1532.7 | 2375.6 | 2994.1 | 3673.4 |

**Marginal cost per chain (ms/100cyc), from the 32->88 slope:**

| | nb4 | nb8 | nb16 |
|---|---|---|---|
| no table | 12.4 | 16.4 | 34.3 |
| table lookup | 23.4 | 26.2 | 42.6 |

**Two conclusions, both firm:**

1. **RT_NB = 4 is the optimum; wider ILP blocks are strictly worse at every chain count**
   (nb8 ~1.3x, nb16 ~2.6x the nb4 cost at 88 chains). Register pressure / occupancy loss
   beats the extra ILP. The block-size axis is closed -- do not revisit it.
2. **Chain cost is dead linear, nothing is free.** This refutes the earlier estimate that
   the first 8-16 chains would be nearly free because the kernel is latency-bound. The ILP
   benefit is already saturated at 4 chains, which is exactly why the 1->4 interleaving in
   `2be5d0d9` won 1.39x and why nothing beyond it is free.

**Table lookup nearly doubles the marginal cost** (12.4 -> 23.4 ms/chain). Any correlated-k
scheme pays this, since per-g-point opacities have to come from a table.

**Consequence for 11 x 8 = 88 chains with table lookup:** RT goes 115 -> ~2100 ms/100cyc,
the span goes ~1530 -> ~3500 ms, so **wall time ~2.3x and RT becomes ~60% of the run**
instead of ~8%. Without the table it is ~1160 ms, ~1.7x wall, ~45% of the run.

**Sub-cycling is NOT the cheap lever it looked like** -- see
[[rt-chain-parallel-split]], which measured t_rad/dt and found 21 % of cells already
radiatively faster than the timestep.

## The table-lookup NaN bug (fixed, verified)

`ktab=true` used to make every cell NaN. `ktab_lookup` clamped the index but computed the
interpolation fraction from the UNCLAMPED coordinate, so a huge extrapolation drove the
opacity negative, flipped the sign of the optical depth increment, and the recurrence
diverged to inf -- which reached the output because `0.0 * inf = NaN`, i.e. the zero-weight
chains are only inert while they stay finite. Fixed by clamping the coordinate before
taking the index; coordinate scaling retuned to (1.2, 0.7) so cells spread over the table.
The `ktab=true` numbers in jobs 10977379 / 10977456 are garbage; the table above supersedes.

## Where things are

`bench/polar_ab/`: `athena_band4/8/16` (all three rebuilt WITH the fix), `bandscale.sh`
(sweep; `BIN=... TAG=... sbatch --export=ALL,BIN=...,TAG=...`), `bandnan.sh`, `bindiff.py`,
`rt_bandharness.patch`. Reference `fixed200.a/bin/dhj.mhd_w_bcc.00001.bin`.
`build-gpu-bench` holds the RT_NB=4 harness build, not the committed production code, and
is shared with another session -- see [[inflight-rt-kernel-optimization]].
Build per [[viper-hip-build-recipe]]. Never write in [[never-write-in-run-dir]].

## Literature the target came from

SPARC/MITgcm (Showman+2009) uses the Marley & McKay two-stream source function with
correlated-k on **11 bands x 8 k-coefficients = 88 column solves**, 0.26-324.86 um, the
binning fixed by Kataria+2013. Lee+2021 (Exo-FMS) found picket fence already reproduces
correlated-k spectra and phase curves closely, so expanding needs a stated science reason
(clouds, disequilibrium chemistry, ultra-hot dissociation, spectral phase curves). Current
frontier is opacity MIXING, not band count: RORR (accurate, too slow in 3D), adaptive
equivalent extinction (Amundsen+2017, >=3x), a DeepSets ML mixer (A&A 2024), a tunable
Monte Carlo k-mixing method (2025).
