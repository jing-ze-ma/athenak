---
name: ck-grav-prod-run
description: ck_grav_prod -- the point-mass-gravity PRODUCTION campaign. GRID REVISED 2026-08-26 to nx1 234 stretched + 2 GPUs. NOT YET LAUNCHED
metadata:
  type: project
---

`/viper/u2/jinma/ATHENAK/bench/ck_grav_prod`, set up 2026-08-26. Same `athena` snapshot as
[[ck-grav-size-run]] (HEAD **c37ebe75**). **This is the science campaign**; it supersedes
[[ck-limb-run]], which is the same setup under the constant-gravity approximation.

| parameter | ck_limb | here |
|---|---|---|
| `problem/grav_point_mass` | false | **true** |
| `mesh/x1max` | 1.433e10 | **2.0556e10** |
| `mesh/nx1`, `meshblock/nx1` | 88 | **200** (dr held at 5.558e7) |
| `time/tlim` | 8.64e7 | 8.64e7 (283 rot) |

**x1max is MEASURED** from ck_grav_size at 21 rotations, not estimated -- both earlier
estimates in [[grav-point-mass-flag]] (1.64e10 and 1.87-1.92e10) were near-IC and are dead.

**The number that changed the answer: "+3 H" was never the delivered standard.** Measuring
[[ck-limb-run]] at 79 rotations the same way, it has **1.76 H** of buffer above its worst
1 % of terminator columns (9.38 H median). Matching that gives 2.0556e10 / nx1 200. A
literal 3 H above the p99 isobar wants 2.166e10 / nx1 220 -- 2.5x ck_limb's cells for a
buffer the accepted run never had. **Calibrate targets against the run the user accepted,
not against a remembered figure.**

Cost ~20 min/rotation, ~95 h = **four 24 h slots**; roughly neutral vs ck_limb because dt
is 2.3x larger with the point mass. `sbatch resubmit_viper.sh` per timeout -- and unlike
ck_grav_size **no `time/tlim` override is needed**, because the input carried the final
tlim from cycle 0 (a restart reads the tlim embedded in the restart file).

NOT yet submitted as of the end of the 2026-08-26 session; smoke test 11058589 on apudev
was still queued. Full rationale in that directory's `NOTES.md`.


## GRID REVISED 2026-08-26 (ninth session) -- nx1 234 stretched, 2 GPUs. NOT LAUNCHED.

The user's requirement: **10 cells per scale height BELOW the 1e-6 bar level.** The uniform
nx1 = 200 grid gives median 7.3 / p05 2.9 there, so it did not meet it. Final design:

| | value |
|---|---|
| `mesh/nx1`, `meshblock/nx1` | **234** (was 200) |
| `mesh/use_grid_stretch_r_poly` | **true** ([[radial-grid-stretch]]) |
| `f_stretch_r_c1..c4` | **+0.080646, -3.645655, +5.697513, -5.725033** |
| submit script | **`--ntasks=2 --gres=gpu:2`** (was gpu:1) |

Everything else is unchanged from the nx1=200 version above; `problem/stellar_tide` is
deliberately absent (OFF) -- see [[stellar-tide-flag]].

Delivers **median 10.1 cells/H below 1e-6 bar** (p05 7.3, p01 5.1), dr 1.73e7 at r/Rp 1.18
to 2.16e8 at the top, dt 22.0 s. Grid costs 1.70x, but the second GPU gives 1.81x
([[meshblock-decomposition-gpu]] correction), so **~90 h = 4 slots -- the same wall time as
the under-resolved single-GPU plan.** 4 GPUs would give 59 h at 31% more node-hours and needs
the `apu` partition (apu1 is MaxNodes=1).

Coefficients were fitted to the `ck_grav_size` H(r) profile at 21 rotations; method and the
metric traps in [[dhj-grid-resolution-design]].

## BUILT AND SMOKE-TESTED 2026-08-27 (tenth session) -- still NOT LAUNCHED, by user request

Everything on the "still to do" list is now done: the input carries nx1 = 234 + the four
stretch coefficients, `athena` is rebuilt from HEAD **1e19d4e7** (BUILD_COMMIT.txt updated),
and both submit scripts are `--ntasks=2 --gres=gpu:2`. Smoke **11105755** on apudev, 2 ranks,
300 cycles: **clean**, 32 meshblocks split 16/16, no fatal, no floors firing.
The user then said **"let's not submit the production run for now"** -- it is ready and held.
`sbatch submit_viper.sh` from that directory is the only step left. src/ is still NOT PUSHED.

**The design's "dt 22.0 s" was WRONG -- measured dt is 10.86 s at cycle 0**, vs 19.15 s for
the uniform nx1 = 200 grid. Root cause MEASURED, not guessed: **the binding CFL direction
switches from azimuthal to radial** once the stretch is on. Full table and the tuning
consequences in [[dt-binding-direction]]. Throughput at cycle 0 is
347 simulated s per wall s on 2 GPUs, against 365 for uniform-200 on ONE GPU -- so the
"1.70x grid cost cancels against 1.81x from the second APU" claim held, but only by luck of
its two errors cancelling. Applying ck_grav_size's measured first-rotation dt settling
(0.64x, then flat) projects **~110 h = five slots, not four**; honest bracket 70-110 h until
the first slot measures the settled dt.
