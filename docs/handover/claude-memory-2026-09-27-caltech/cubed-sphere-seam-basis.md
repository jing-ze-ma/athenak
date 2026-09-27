---
name: cubed-sphere-seam-basis
description: The cubed-sphere panel halo needs a tangent-BASIS transform, not a signed component permutation -- FIXED 2026-08-27, rigid rotation now converges
metadata:
  type: project
---

Fixed 2026-08-27 (session 4c1af17c), uncommitted with the rest of the cubed-sphere work.
**This supersedes the "BUG 2 = along-seam interpolation" diagnosis in
[[cubed-sphere-hydro-fix]]**: the interpolation gap is real but SUBLEADING. The dominant
defect was the vector transform.

## The defect

Across a panel seam the two charts do NOT share a tangent basis. Their eta axes agree
along the shared edge, but the xi axes differ by a shear that vanishes only on the seam
midline: for the 0/1 seam, `e_xi = e_xi' - sqrt(2) sin(eta) e_eta'`, reaching O(1) at the
panel corners. `bvals_cc.cpp` moved the momentum components with a signed axis
permutation, which is exact for the INDICES and exact for the COMPONENTS only at eta = 0.
Measured ghost velocity error: **0.36, flat in resolution** at nx = 16/32/64.

An O(1) ghost velocity makes the flux difference O(1) and the acceleration O(1) -- the
scheme was INCONSISTENT at every seam.

## The fix

New `src/coordinates/cubed_sphere.hpp` (header-only, device-callable): the six panel
frames, `PanelToCart` / `CartToPanel` / `FindPanel` / `PanelTangents`, and
`TransformMomentum`, which raises with the source panel's Gram matrix, rebuilds the vector
in Cartesian, and dots into the destination panel's basis. `m_i = V.e_i`, so the round
trip is exact. IM1 is radial and untouched -- rhat is common to both charts.

`bvals_cc.cpp` calls it per SOURCE CELL in the pack kernel (`srcval` lambda), replacing
`map_vy/map_vz/sign_vy/sign_vz`. The INDEX transpose/reversal is unchanged -- it was
verified correct (below). Off the cubed sphere `cs_xform` is false and the lambda is the
old `a(...)*signvar` textually; **linear_wave_hydro is BIT-IDENTICAL before vs after**,
verified against a reconstructed pre-edit file, not by inspection alone.

## Measured, nx = 16/32/64

| quantity | before | after |
|---|---|---|
| ghost d(v_xi), d(v_eta) | 0.18 / 0.36, FLAT | 5.9e-2 / 3.2e-2 -> 1.4e-2 / 7.0e-3, 1st order |
| ghost d(v_r) | 1.07e-3, 1st order | 9.4e-5 -> 5.2e-6, 2nd order |
| one-step seam max\|v_r\| | 1.07e-3, 1st order | 6.9e-5 -> 4.4e-6, 2nd order |
| t=1 radial KE (pure error) | 5.9e-4 -> 8.4e-4, DIVERGING | 2.54e-6 -> 9.98e-8, ratio ~5 |
| t=1 dM/M | +1.5e-3 | +6.7e-5 -> +3.8e-5 |

Test A (uniform static) is still exact -- mass 24.000000, totE 36.000000, KE ~1e-30, and
dt 8.781087e-03 bit-identical to the pre-existing reference. `iprob=5` radial-halo test is
0.0 on all six panels. The `iprob=6` scalar ghost error is bit-identical (the transform
touches momenta only).

## What is left -- BOTH ITEMS BELOW ARE NOW CLOSED, see [[cubed-sphere-seam-interp]]

The along-seam offset was fixed 2026-08-27 by a quadratic resample in the halo, and the
"interior residual, flat in resolution" was **not a defect at all** -- it is radial
truncation error, and every scan that called it flat had held nx1 fixed at 8. src1/src2
are exonerated. The original text follows.

## What is left (original)

* **The along-seam offset**, i.e. the original BUG 2. Now cleanly isolated and cleanly
  FIRST ORDER: the ghost sits ~0.5 cell (layer 0, panel corner) from the neighbour cell
  that fills it. Fix = interpolate along the seam-parallel index; quadratic or better is
  needed to get the acceleration to 2nd order.
* **The interior residual**, which the seam fix did NOT touch and which is now the LARGER
  term at nx=64 (interior sum v_r^2 3.91e-9 vs seam 1.57e-9). Flat in resolution as an
  acceleration. Not the halo -- it survives an exact-ghost bypass. Suspect the `src1`/
  `src2` horizontal source terms ([[cubed-sphere-hydro-state]] suspect (a)).
* `bvals_fc.cpp` (MHD) has the same signed-permutation shape and is untouched. Face-centred
  B needs its own treatment; see [[cubed-sphere-for-hot-jupiter]] item 1.
