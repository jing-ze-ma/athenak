---
name: cubed-sphere-mhd-seam
description: The face-centred (bvals_fc) panel seam for MHD -- FIXED and committed 4f19a244; halo is 2nd order and the seam is no longer a distinguished error source. The remaining defect is the panel-CORNER internal energy, which pre-exists
metadata:
  type: project
---

Fixed 2026-08-27/28 (thirteenth session), committed **4f19a244** on top of the WIP
a8148e2d. Supersedes this file's earlier "still wrong" state entirely.
**START HERE for cs+MHD.** Read [[cubed-sphere-mhd]] for the layers below it.

## The one thing to carry forward: the transform was never the problem

a8148e2d guessed the residual was "a residual BASIS error, worst at the panel corners".
It was not. A 40-line standalone harness (strip `#include "athena.hpp"`, `#define
KOKKOS_INLINE_FUNCTION inline`, compile `cubed_sphere.hpp` with plain g++) shows
`TransformFieldToDstNormals` maps exact source projections to exact destination ones to
**<=3.3e-16 on all 24 directed seams, both output components**. Doing that FIRST would
have saved the whole previous session's mis-diagnosis. Generalisable: unit-test the pure
function before blaming it from an end-to-end number.

## The three real defects, all plumbing

1. **The secondary-component average read the SOURCE BLOCK'S OWN GHOSTS at the two END
   faces of every seam.** Only 2 of 33 faces per seam were wrong, by O(1) (0.38 vs
   3e-3..1.7e-2 on the interior faces at nx=32). The halo check maxes over the seam, so
   two bad faces made all 24 seams read O(1) and non-convergent. Now extrapolated from
   the two nearest ACTIVE cells. **This is what a8148e2d misread as a basis error.**
2. **No along-seam resample.** The face twin of [[cubed-sphere-seam-interp]]; same
   quadratic Lagrange, same `atan(tan(a)*tan|n|)` formula. `b.x1f` needs it too even
   though it crosses as a scalar. New wrinkle: the STAGGERING differs per component --
   on a x3-face seam the along-seam axis is xi, on FACES for `b.x2f` and on centres for
   the other two -- handled by an `angles` lambda.
3. **No valid corner halo.** A panel corner is a CUBE VERTEX where only THREE panels
   meet, so the generic corner buffer reaches somewhere meaningless: corner ghosts were
   wrong by ~0.7-1.1 of |B| on EVERY component, **the radial scalar included and flat in
   resolution** -- garbage, not a basis error. New `MeshBoundaryValuesFC::
   FillPanelCornersFC` extrapolates each corner ghost quadratically from the two flanking
   face halos and averages; called at the end of `RecvAndUnpackFC` where every face
   buffer is known unpacked. **This was the entire residual after 1 and 2.**

## Measured (nx2=nx3, nx1=8)

Halo vs exact, max over 24 seams at 16/32/64: `x1f` 1.30e-1/6.65e-2/3.35e-2 ->
7.26e-5/1.08e-5/1.47e-6; `x2f`/`x3f` ~0.43-0.46 FLAT -> 1.10e-3/2.77e-4/6.95e-5, ratios
**3.97, 3.98**. Corner ghosts 0.7-1.08 flat -> 6.33e-3/8.56e-4/1.07e-4.

One-step seam max|v| at nx=16/48/96/128: before 8.63e-3/8.94e-3/8.69e-3/8.62e-3 (FLAT =
inconsistent) -> after 4.46e-4/6.36e-5/1.35e-5/8.37e-6, **at or below the panel interior
at every resolution**.

`linear_wave_mhd` BIT-IDENTICAL; Test A-MHD 1-mom 5.32870e-03 and dt 8.781087e-03 as
documented; cpplint clean.

## RUN THE EVOLVED GATE WITH `user` RADIAL BCs, NOT `reflect`

`mesh/ix1_bc=reflect ox1_bc=reflect` is documented for the **nlim=0** energy check ONLY.
iprob=8 enrols `CSTestRadialBC` via `user_bcs_func`, and the input file's own default is
`user`. Running the evolved test with `reflect` manufactures a radial jump and changed
which resolutions failed -- it cost a confounded scan this session.

## STILL BROKEN, and now the leading defect: the panel-corner ENERGY

After one step the internal energy goes strongly negative in the panel-CORNER cells --
exactly 1.5000 at nlim=0, **-5.2 at nx=48 after one step** -- on all six panels at every
resolution, while |v| there stays at truncation level. **The |v| metric hides it
completely; always print min internal energy and its (i,j,k) too.** It **PRE-EXISTS**
this commit (-21 at the same point with the old bvals_fc; the halo work improves it 4x
and no further), so it is a separate mechanism, presumably the untouched **seam EMF
consistency**: `buffs_fc.cpp` never exchanges the shared seam face, and two panels do not
produce identical EMFs on it ([[cubed-sphere-for-hot-jupiter]] item 1).

Related, and why the gate still cannot be called passed: with the tracked tilt at
**nx=32 and nx=64 specifically**, four of six panels blow up at a corner (max|v| 0.51 and
35, minrho down to 0.46). It is a **knife edge, not a trend** -- 16/48/96/128 are clean,
tilt (0.8,0.1,0.59) is clean at both, `dc`/`hlle` make it far worse, `mb=32` at nx=64 is
clean, and the base fails there too (381 and 128). Deterministic: Debug and Release agree
to all digits, and a Debug build reports **no** bounds violation.

## Harness

Session 28cde4ec scratchpad: `bcs/` (PROBLEM=cs_test, current), `bcs_base/` (HEAD's
bvals_fc for A/B -- rebuild it by `git checkout src/bvals/bvals_fc.cpp`, make, restore),
`bdbg/` (Debug), `bold`/`bnew` (built-in pgen, for the linear_wave bit-identical check),
`t.cpp`+`cs_hdr.hpp` (the standalone transform unit test), `runs/`.
Diagnostics added to `cs_test.cpp`: the halo check now reports a `cnr` row per panel
(corner ghosts) and the localisation reports max|v| in corner cells, the (j,k) of the
seam max, and min rho / min internal energy with its location.
