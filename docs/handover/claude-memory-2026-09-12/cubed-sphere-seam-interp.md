---
name: cubed-sphere-seam-interp
description: The cubed-sphere panel halo now RESAMPLES along the seam -- the hydro scheme is uniformly 2nd order, seams no worse than panel interiors (2026-08-27)
metadata:
  type: project
---

Fixed 2026-08-27 (session e6bb86fb), uncommitted with the rest of the cubed-sphere work.
This closes "BUG 2 / the along-seam offset" left open by [[cubed-sphere-seam-basis]] and
[[cubed-sphere-hydro-fix]]. **The gnomonic hydro scheme is now second order everywhere.**

## The defect and the fix

Across a seam the two charts share the seam-NORMAL coordinate exactly but not the
seam-PARALLEL one. For a source cell with seam-normal angle n and seam-parallel angle a,
its physical point sits at seam-parallel angle `atan(tan(a)/tan|n|)` in the DESTINATION
chart -- not at `a`. The plain index copy therefore handed every ghost the state of a
point up to 0.5 cells (layer 0) or ~1.5 cells (layer 1) away along the seam, an offset
that is resolution-INDEPENDENT in cell units. Ghost value O(dx) wrong -> flux difference
O(dx) -> acceleration O(1).

Inverting the map, the value a ghost needs is the source field at seam-parallel angle
**`atan(tan(a)*tan|n|)`**, which is always closer to the seam midline than the source
cell's own angle, so the stencil never leaves the source block. **The same formula covers
both seam orientations** -- a reversed seam flips the sign of both `a` and the target, and
the existing index reversal already carries that. Implemented as a quadratic (3-point)
Lagrange resample in `bvals_cc.cpp` (`seamval`, wrapping the existing `srcval`), applied
only to same-level FACE buffers on a panel seam; the ng x ng edge/corner buffers stay a
plain copy, which a dimensionally split PLM+HLLC sweep never reads through.

## Measured, joint refinement (nx1 = 8/16/32, nx2 = nx3 = 16/32/64)

One-step residual acceleration and the ghost errors:

| | pre | post | order pre -> post |
|---|---|---|---|
| seam accel | 9.02e-4 -> 3.94e-4 | 1.04e-4 -> 6.22e-6 | 0.6 -> **2.0** |
| interior accel | 9.69e-5 -> 6.01e-6 | identical (halo-only change) | 2.0 -> 2.0 |
| ghost d(rho) | 1.14e-3 -> 2.62e-4 | 9.36e-5 -> 5.91e-6 | 1.0 -> **2.0** |
| ghost d(v_xi) | 5.86e-2 -> 1.41e-2 | 8.55e-5 -> 5.19e-6 | 1.0 -> **2.0** |

**The seam residual is now within 5% of the panel interior at every resolution** -- the
seam is no longer a distinguished error source. t=1 rigid rotation, radial KE (a pure
error field): 2.54e-6/5.17e-7/9.58e-8 (1st order) -> **5.58e-8/2.39e-9/9.44e-11**, ratios
23-25, i.e. 2.3rd order in v_r; a 45x / 216x / 1015x gain. dM/M 6.7e-5 -> 1.2e-5 at nx=16
and below hst precision at nx=64. `iprob=6` scalar ghost max: 7.0e-2/3.6e-2/1.8e-2 ->
5.3e-5/7.8e-6/1.1e-6.

Unchanged: Test A exact (mass 24.000000, totE 36.000000, KE ~1e-30, dt 8.78109e-03
bit-identical); `iprob=5` at round-off (2.2e-16, was 0.0 -- the weights sum to 1 only to
machine epsilon); **linear_wave_hydro BIT-IDENTICAL**, since `cs_seam` can only be set
inside the `use_cubed_sphere` branch.

## RETRACTION: the "interior residual, flat in resolution" was an artifact

[[cubed-sphere-seam-basis]] listed a second open defect -- an interior acceleration flat
under refinement, suspected to be the `src1`/`src2` source terms. **It is not a defect.**
Every earlier convergence scan refined nx2/nx3 only, with **nx1 = 8 held fixed**, and the
interior residual is ordinary RADIAL truncation error, which angular refinement cannot
touch. Refine nx1 as well and it converges at ratio 4.01, 4.01 -- clean second order.
`src1`/`src2` are exonerated. Generalisable lesson: on a grid where one axis is not being
refined, "flat in resolution" says nothing about consistency.

## Also fixed: a stale axis check in build_tree.cpp

`BuildTreeFromScratch` still rejected `nmb_rootx1 != nmb_rootx2` -- a leftover from before
[[cubed-sphere-x1-radial]], when x1/x2 were the tangential pair. The panel-tangential axes
are now x2/x3, so the constraint is `nmb_rootx2 == nmb_rootx3`, and the radial direction
is free. **Multiple MeshBlocks per panel now work and are essentially decomposition-
independent**: 2x2 and 4x4 blocks per panel give 1-KE 2.38785e-9 / 2.39421e-9 vs
2.39364e-9 single-block (0.2% and 0.02%). The resample clamps its stencil to the source
block's ACTIVE range, never its ghosts, and the shift always points toward the seam
midline, which is why decomposition costs almost nothing.

## What is left

* `bvals_fc.cpp` (MHD) still has the plain signed-permutation copy -- it needs BOTH the
  basis transform and this resample, plus seam EMF consistency
  ([[cubed-sphere-for-hot-jupiter]] item 1).
* The 2D (nx3 = 1) shell case, never passing.
* Corner/edge halo buffers are still plain copies. Measured not to matter for PLM+HLLC;
  a higher-order reconstruction would need them.

Harness: session e6bb86fb scratchpad -- `athena.pre` / `athena.post2` binaries,
`jres*`/`Jres*_pre` (one-step residual), `jc_{pre,post}_*` (t=1), `mbt`/`mb4`
(decomposition), `p5_*`/`p6_*`. Build dir is the 8c7040b0 scratchpad's `bcs/`
(PROBLEM=cs_test); `prev/nb` is the built-in-pgen build for linear_wave.
