---
name: cubed-sphere-hydro-state
description: Cubed-sphere (gnomonic equiangular) hydro in AthenaK -- what was broken, what is FIXED, and the one measured inconsistency that remains
metadata:
  type: project
---

The user wrote the cubed-sphere infrastructure themselves (Feb-Mar 2026: `9e37cfb7`,
`019e6ce3`, `ce3d13fe`, folded into `35503389`, whose own message says the gnomonic grid
was "not successful yet"). It is **far more complete than "should we build it"** implies:
six panel trees, panel-neighbour + orientation tables (`swap_ax`/`rev_x1`/`rev_x2`),
inter-panel exchange for **both** cell- and face-centred fields, the gnomonic metric, a
generic curvilinear CT update, and geometric source terms for hydro and MHD.

## FIXED 2026-08-27 (uncommitted, working tree)

All changes are guarded by `use_cubed_sphere`, so the spherical-polar production path is
untouched; `deep_hot_jupiter_rt` still builds.

1. **`cubed_sphere.cpp` / `solid_body_rot.cpp` do not compile** -- stale `use_wellbalance`
   (split into `_static`/`_dynamic`). They also write kinetic energy into `w0(IEN)`, which
   is the INTERNAL energy density `p/(gamma-1)`, and their velocity transforms are
   single-panel. **Do not diagnose the grid through them.** Wrote `src/pgen/cs_test.cpp`
   instead, plus `inputs/tests/cubed_sphere_{uniform,rigidrot}.athinput`.
2. **2D was instant NaN.** The 2D branch set `r_c = r_l = r_r = 1`, but areas carry
   `(r_r^2-r_l^2)` and the volume `(r_r^3-r_l^3)` -- so every area and volume was
   identically zero and `x_ov_rD`/`y_ov_rC`/`z_ov_rE` divided by it. Now uses the real
   x3 range. **The cubed sphere is a 3D construct**: even fixed, a 2D shell has an
   unbalanced radial pressure source (`z_ov_rE`) with no radial flux to cancel it.
3. **`dt` was exactly 0 in 3D.** `CoordGnomonicEquiangle` set `dx1`,`dx2` but never
   `dx3` (the radial width), so `min_dt3 = dx3/(|v|+cs) = 0`. One line.
4. **The coordinate arrays covered active cells only** (`ks,ke,js,je,is,ie`), leaving all
   geometry zero in ghosts -- `CoordSphericalPolar` had already been changed to loop the
   full array. Now matches. This exposed an **out-of-bounds allocation**: `sin/cos_face1`
   are written at `i+1` and `sin/cos_face2` at `j+1`, but both were sized `ncells1`/
   `ncells2` with no `+1`. Symptom was `malloc(): unaligned tcache chunk`.
5. Hoisted a `pmy_pack->pmesh->three_d` **host-pointer dereference out of a device
   lambda** -- it only worked under Viper's unified memory.

**MILESTONE REACHED: the static uniform state is now exact.** `iprob=1`, 3D, 6 panels:
mass stays 24.000000, total energy 36.000000, momenta at ~1e-17 and KE at ~1e-30 for the
whole run. So the geometry, the geometric pressure source terms, and the inter-panel
exchange are all mutually consistent when v = 0.

## SUPERSEDED 2026-08-27 -- read [[cubed-sphere-seam-basis]] first

The section below is RETRACTED. Its rigid-rotation numbers were measured with a pgen whose
`PanelToCart` had panels 3 and 4 interchanged ([[cubed-sphere-panel-frames]]), so iprob=3
was not a rigid rotation on two of six panels. The real cause was suspect (b), the
panel-seam basis transform -- now fixed; the rigid rotation converges. Suspect (a), the
src1/src2 source terms, is still open and is now the LARGEST remaining term.

## RETRACTED: a non-converging inconsistency with nonzero velocity

`iprob=3` is a rigid-rotation equilibrium (`v = Omega zhat x r`, `p = p0 + rho Om^2 R^2/2`)
-- an EXACT steady solution of the Euler equations. It drifts, and **the drift does not
converge**: spurious radial KE (identically zero in the exact solution) is
**1.03e-2 / 1.06e-2 / 1.08e-2 at nx1 = 16 / 32 / 64** at t=1, where 2nd order would cut it
~16x per doubling. Error is worst on panel seams (|v_r| ~ 0.11) but broadly present in
panel interiors (~0.04).

**What the convention IS, established by experiment not by reading:**
`w0` velocity is CONTRAVARIANT; `u0` momentum is COVARIANT (the tail of
`GnomonicEquiangleFluxX*` lowers the index). I tested the opposite reading by deleting that
lowering -- it **breaks the exact uniform state** (mass 24.07, momentum 2e-2) and makes the
rotation error 5x worse. So the lowering is correct and necessary.

I added `Coordinates::GnomonicEquiangleRaiseVel`/`LowerMom` (metric raise after ConsToPrim,
since `SingleC2P_IdealHyd` does `v = m/rho` and `KE = 0.5 m.m/rho`, both of which assume an
orthonormal basis, while `SrcTermsGnomonicEquiangle` uses the correct
`0.5*rho*(v1*v_1+v2*v_2)`). **That inconsistency is real, but fixing it did NOT change the
drift** -- so it is not the dominant term. Kept, clearly flagged, not claimed as a fix.

**Remaining suspects, untested:** (a) the `src1`/`src2` horizontal source terms, where the
`- rho*v3*v_1/radius` form is commented out and replaced by a radial-flux expression;
(b) the panel-seam basis transform -- the bvals panel map is a signed permutation, but the
two panels at a seam have *different* non-orthogonal bases, which matches the seam-localised
part of the error.

**Method that worked and should be reused:** a uniform static state first (isolates geometry
from velocity), then an exact steady state WITH velocity, then a **convergence scan** -- the
scan is what separates "truncation error" from "inconsistent scheme", and it is the only
reason the remaining bug is known to be real. See [[measure-impact-before-claiming]].

For the MHD outlook see [[cubed-sphere-for-hot-jupiter]].
