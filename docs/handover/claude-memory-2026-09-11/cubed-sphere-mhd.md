---
name: cubed-sphere-mhd
description: MHD on the cubed sphere -- flux path, CT edge lengths, Maxwell source terms and the orthonormal bcc frame are all DONE and committed; the seam halo is the one remaining gap
metadata: 
  node_type: memory
  type: project
  originSessionId: 7d69edbe-adbc-4c35-a568-39cc17d8cf13
  modified: 2026-08-27T16:58:47.178Z
---

Started 2026-08-27 (twelfth session), on top of [[cubed-sphere-committed]] (45f4cb43).
Before this, MHD + `use_cubed_sphere` **ran and silently produced garbage**.

## Fixed and validated this session

1. **`dxedge` was never filled for the gnomonic coordinates** (only for spherical polar),
   so the curvilinear CT branch in `mhd_ct.cpp` multiplied every EMF by an uninitialised
   ZERO edge length: B could never evolve at all. Added in `CoordGnomonicEquiangle`, each
   edge taking its arc from the same `dth_*` the adjoining face area uses. Cross-checked
   `area.x1f` vs `|x2e||x3e|sin(theta)`: 0.0237 -> 0.0121 under 2x refinement, clean O(dx),
   which is exactly the face-arc vs cell-arc difference.
2. **`mhd_fluxes.cpp` never called ANY gnomonic rotation** -- it declared
   `use_cubed_sphere` and did nothing with it -- while `mhd_tasks.cpp` did call
   `SrcTermsGnomonicEquiangle`, which assumes covariant momentum fluxes. Added the six
   `GnomonicEquianglePrimFaceX*` / `GnomonicEquiangleFluxX*` calls mirroring
   `hydro_fluxes.cpp`.
3. **Geometric source terms had no Maxwell stress.** New
   `SrcTermsGnomonicEquiangleMHD`; hydro and MHD share one `...Impl` kernel rather than
   being duplicated the way `SrcTermsSphericalPolar{Hydro,MHD}` are. The rule is the
   general one: `p -> p + B^2/2` wherever p multiplies a metric derivative, and
   `rho v^a v^b -> rho v^a v^b - B^a B^b`.
4. New coordinate helpers `GnomonicEquiangleFaceBX2` and `GnomonicEquiangleEmfX1`.
5. **`bcc` was a non-orthogonal triple** and MHD had no `RaiseVel` counterpart -- see the
   section below, commit 5f635505.

## The frame analysis -- the part worth not re-deriving

Basis `e_xi, e_eta` are unit vectors with `e_xi.e_eta = c = cos_cell`, `s = sin_cell`.
Face normals: `nhat_xi = (e_xi - c e_eta)/s`, and `nhat_xi.nhat_eta = -c`.

**Storage convention, FACES (adopted, and forced by `mhd_ct.cpp`):** `b0.x*f` holds the
FLUX DENSITY `B.nhat` through its own face, so `sum(area*b)` is the physical flux and the
existing `area`/`dxedge` CT update is already right. Do not change this one.

**Storage convention, CELL CENTRES (final, as of 5f635505):** `bcc` holds the ORTHONORMAL
triple `(B.rhat, B.e_xi, B.f2)` with `f2 = (e_eta - c e_xi)/s`. The raw face average is
`(B^r, s*B^xi, s*B^eta)`, whose three directions are NOT mutually orthogonal -- that was
the bug. `w0` holds the CONTRAVARIANT velocity.

Per-sweep orthonormal frames (the Riemann solver permutes v and B identically). The NORMAL
component always comes from `b0.x*f` and is already correct, so only the two face-parallel
slots can be wrong -- and with bcc orthonormal, only ONE sweep needs anything:
* x1: `{rhat, e_xi, f2}` -- slots 1,2 already exactly right, nothing to do
* x2: `{nhat_xi, e_eta, rhat}` -- slot 2 must become `B.e_eta = c*b1 + s*b2`
* x3: `{nhat_eta, rhat, e_xi}` -- slots 0,1 already exactly right, nothing to do

**EMFs: only the x1 sweep needs rotating.** CT integrates E along EDGE directions
(`dxedge`), and the x2/x3 frames' two face-parallel axes ARE edge directions, so those
EMFs pass through untouched. The x1 frame's third axis is `(e_eta - c e_xi)/s`, not
`e_eta`, so `E.e_eta = c*e21 + s*e31`.

## Measured

* `cs_test` gained an MHD path and `problem/b0r`: a **radial monopole** `B_r = b0r*(r0/r)^2`
  with v = 0. It is force-free (curl B = 0) so the exact solution is STATIC, and it is
  exactly div-free on this grid because `area.x1f` carries r^2 while the angular faces
  carry no flux -- independently of the panel map.
* Gate 1, B = 0: MHD reproduces hydro Test A **exactly** -- mass 24.000000, totE 36.000000,
  dt 8.78109e-03, KE ~1e-30. Before the fix it was mass 23.9257 and KE 4e-4.
* Gate 2, monopole: angular momenta and KE at round-off (~1e-31), magnetic energy
  **exactly** constant at 3.56007. Radial residual converges at **2nd order** in nx1:
  1-mom 5.33e-3 / 1.32e-3 / 3.25e-4 (ratios 4.05, 4.05) at nx1 = 8/16/32.
* `cubed_sphere_{uniform,rigidrot}.athinput` re-run: **BIT-IDENTICAL** to the references.

**Watch out: `reflect` is a trap for the radial BC here.** `B_r(r)` is not
reflection-symmetric, so reflect manufactures a jump at the radial boundary -- it gave
mass 24 -> 24.72 and radial KE 2.8e-2, i.e. 3 orders of magnitude worse, and it looks
exactly like an interior bug. `cs_test` now serves `iprob=1` from the exact-state user BC
(`ix1_bc=ox1_bc=user`), including the ghost FACES, since `user` is not a case
`bfield_bcs.cpp` handles. Same class of trap as the `ix1_bc` note in
[[cubed-sphere-committed]].

## Gap 2 is CLOSED -- the orthonormal bcc frame (commit 5f635505)

`bcc` was the raw face average `(B.rhat, B.nhat_xi, B.nhat_eta)`, and those three
directions are NOT mutually orthogonal (`nhat_xi.nhat_eta = -c`), so every sum of their
squares was wrong by an O(1) amount -- the MHD c2p's magnetic energy, `PrimToCons`, and
the fast speed in `mhd_newdt`. **MHD was also missing the `GnomonicEquiangleRaiseVel`
call that `hydro_tasks.cpp` makes after ConsToPrim**, so the primitive velocity held
COVARIANT components where every consumer expects contravariant ones.

New `GnomonicEquiangleRaiseVelMHD` does both: rebuilds `bcc` from the faces (idempotent),
rotates it into `{rhat, e_xi, (e_eta - c e_xi)/s}`, raises the velocity, and subtracts a
magnetic energy that is now genuinely the sum of squares. **Storing bcc that way makes the
flux path SIMPLER, not harder** -- x1 and x3 need no field rotation at all, so
`GnomonicEquiangleFaceBX1/X3` are gone and only `FaceBX2` survives.

**The test that mattered: `cs_test` `iprob=8`, a uniform CARTESIAN field.** Force-free and
exactly static like the monopole, but with tangential components on every panel -- and the
monopole CANNOT see this bug, because a purely radial field is identical under both
conventions. `|B|^2 = b0c^2` analytically, so `u0(IEN)` is seeded with the exact magnetic
energy and the recovered pressure measures what the inversion actually subtracted.
**Seeding that from `bcc` instead would have made the check self-consistent under any
convention rather than correct** -- the same circular-validation trap as
[[eos-inversion-nan-trap]].

| max abs(dp/p), nx2=nx3 = | 16 | 32 | 64 | |
|---|---|---|---|---|
| after | 8.01e-4 | 2.01e-4 | 5.02e-5 | ratios 3.99, 4.00 -- 2nd order |
| before | 9.27e-2 | 1.02e-1 | 1.07e-1 | ~10%, and NOT convergent |

The old error *growing* under refinement is the signature of an inconsistency rather than
an inaccuracy -- exactly what the note at `GnomonicEquiangleRaiseVel` records for the
velocity. No regressions: both hydro inputs bit-identical, MHD at b0r=0 still exactly
reproduces hydro Test A, monopole unchanged to every printed digit.

## The ONE remaining gap

**`bvals_fc.cpp` is still the old signed-permutation seam copy** -- no tangent-BASIS
transform, no along-seam resample, i.e. exactly the two things that took hydro from 1st to
2nd order ([[cubed-sphere-seam-basis]], [[cubed-sphere-seam-interp]]). Neither MHD test
touches it: the monopole's `B_r` is a function of r alone with zero angular face fields, and
`iprob=8` is run at `nlim=0`. **Start there**, and note that the face-centred halo is harder
than the cell-centred one -- the seam faces are SHARED between panels, so the normal
component must agree exactly or div B is broken, and the seam EMFs must be single-valued
([[cubed-sphere-for-hot-jupiter]] item 1). `iprob=8` evolved for a few steps is the natural
gate: it is an exact static state, so any seam error shows up immediately as a force.

Also still open, lower priority: corner/edge halo buffers are plain copies (measured
irrelevant for PLM+HLLC in hydro), the 2D nx3=1 shell, and `ConsToPrimCoarseBndry` has no
gnomonic treatment so SMR/AMR on the cubed sphere is untested.
