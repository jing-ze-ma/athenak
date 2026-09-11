---
name: cubed-sphere-x1-radial
description: DONE - the cubed sphere now uses x1 = radial, x2/x3 = panel-tangential, matching the spherical-polar grid
metadata:
  type: project
---

Requested by the user 2026-08-27: *"I think for later, we might want to make x1 as radial
for cs to be compatible with spherical grid"*, then *"do the x1 radial refactor first"*.
**DONE and verified the same day** (uncommitted, with [[cubed-sphere-hydro-fix]]).

The cubed sphere was xi = x1, eta = x2, radius = x3; it is now **radius = x1, xi = x2,
eta = x3**, matching `CoordSphericalPolar` (r, theta, phi). The whole change is one CYCLIC
PERMUTATION of the axis labels, x1->x2->x3->x1, which is why it could be verified by
bit-comparison (see below).

Files changed: `mesh.hpp` (face ordering, `NeighborIndexPanel`, `PanelBoundaries` fields),
`meshblock_tree.cpp` (`TransformToPanel`, `FindNeighborGlobal`, `FindNeighborCrossPanel`),
`meshblock.cpp` (see the gotcha below), `mesh.cpp` (the nx2 == nx3 check and the diagnostic
print), `coordinates.cpp` (`CoordGnomonicEquiangle`, `SrcTermsGnomonicEquiangle`,
`GnomonicEquiangleRaiseVel/LowerMom`, the array allocations), `coordinates.hpp` (the six
`GnomonicEquiangle{PrimFace,Flux}X*` routines rotate roles: X1 is now the RADIAL one),
`bvals_cc.cpp`, `hydro_fluxes.cpp` (the PrimFace call sites now pass `k`), `cs_test.cpp`.

**Renames made deliberately, to stop the bug in [[cubed-sphere-hydro-fix]] recurring:**
* `PanelBoundaries::rev_x1/rev_x2` -> **`rev_a`/`rev_b`**, the two TANGENTIAL axes. They
  were never mesh x1/x2, and reading them as such is what reversed the radial index.
* `sin_face1/cos_face1` -> `sin_face_xi/cos_face_xi` (on x2 faces, staggered in j);
  `sin_face2/cos_face2` -> `sin_face_eta/cos_face_eta` (on x3 faces, staggered in k).
  All four, and `sin_cell/cos_cell`, are now indexed **(m,k,j)** -- angle-only, no radial
  extent. `z_ov_rE` stays the RADIAL coefficient, which is also its meaning in
  `SrcTermsSphericalPolar`, so the two grids now agree.

**THE GOTCHA that cost the most time.** `meshblock.cpp` builds neighbours branch by branch
per face type, and the panel transform was applied in six of them -- every branch with a
nonzero x1 or x2 offset -- but **there was no pure x3-face branch call**, because under the
old convention x3 was radial and no seam crossed it. After the swap the x3 face IS a seam.
Symptom: panels 3 and 5 (the ones whose `GetPanelBoundary` entries set `swap_ax`) had
completely EMPTY halos. Fixed by adding the `NeighborIndexPanel(0,0,-l,...)` call there.
The mirror-image vestige -- the transform still sitting in the x1-face branch -- is now
unreachable (its `panel != nt->lloc_.panel` guard can never fire) and was left in place.

**Input files change too**: x1 is now the radial range with physical BCs
(`ix1_bc`/`ox1_bc` = user/reflect) and x2/x3 are [-1,1] with `panel` BCs; `mesh.cpp` now
requires **nx2 == nx3** rather than nx1 == nx2.

**Verification: the refactor is provably a pure relabelling.** Every diagnostic reproduced
its pre-refactor value BIT FOR BIT -- Test A dt 8.78109e-03 and exact to round-off; the
one-step residual seam/interior 1.7443e-2, 1.6616e-2, 1.6242e-2 / 9.6887e-5, 8.1470e-5,
8.0530e-5; the t=1 rigid-rotation arms dM/M +1.5125e-3, +1.3000e-3, +6.2083e-4 and radial
KE 7.10107e-3, 7.64177e-3, 7.89875e-3. The `iprob=5` radial-halo test is 0.0 on all six
panels. Note the history columns MOVE: the radial error is now **1-KE**, not 3-KE.

Still open, unchanged by this: the along-seam interpolation ([[cubed-sphere-hydro-fix]]),
`bvals_fc.cpp` for MHD (its j/k + IVY/IVZ form is now the CORRECT one and only needed the
field renames, but it is untested), and the 2D nx3=1 case.
