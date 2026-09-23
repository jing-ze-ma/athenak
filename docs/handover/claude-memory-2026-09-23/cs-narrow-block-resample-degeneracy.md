---
name: cs-narrow-block-resample-degeneracy
description: The cubed-sphere along-seam resample needs 3 cells along the seam; its stencil bounds INVERT below that and it reads unfilled ghost. Bites on the COARSE array whenever cnx = nx/2 < 3, i.e. any MeshBlock under 6 cells wide, and MeshBlocks are allowed down to 4
metadata:
  type: project
---

Found 2026-08-31 while deriving a least-squares corner stencil, NOT by any test failing --
no input file in the repo uses blocks this narrow.

## The defect

The along-seam resample clamps its 3-point quadratic stencil to
`[blo, bhi] = [start, start + extent - 3]`. When the source's along-seam ACTIVE extent is
under 3 that range **INVERTS** (bhi < blo), the clamp pins the stencil one cell OUTSIDE
the active zone, and it reads unfilled ghost memory.

`extent` is the COARSE extent `cnx = nx/2` on a coarse buffer, so this bites for any
MeshBlock under **6** cells wide -- and `Mesh` allows MeshBlocks down to **4**.

## Measured, refined 4x4 MeshBlocks (iprob=11 halo Linf, |B| ~ 0.2)

|  | tangent | r x EDGE SEAM | r x CORNER SEAM |
|---|---|---|---|
| session start (790d02c0) | 2.4471e-02 | 6.9313e-02 | 2.0057e-01 |
| before the guard | **7.6746e-01** | 3.0473e-01 | **8.0250e-01** |
| with the guard | 2.4471e-02 | 6.2278e-02 | **6.2446e-02** |

**The degeneracy PREDATES the session, but a01ace75 and e4be92c0 widened its blast radius**
by resampling many more buffers -- so those commits made narrow-block meshes 4-31x worse
while making the standard ones much better. Guarded now: fall back to the plain index copy
when the extent is under 3. That costs the O(dx) seam offset the resample exists to remove,
which is bad but BOUNDED; reading unfilled memory is not. Both `bvals_fc.cpp` and
`bvals_cc.cpp` carry the guard.

**Ruled out**: it is NOT the coarse cube-vertex corner fill of 180a9b3e -- identical with
that on and off.

## THE PROCESS LESSON

Every input file in `inputs/tests/` uses wide blocks, so no gate in the repo exercises
`cnx < 3`. **Whenever a change touches stencil bounds, test at the MINIMUM legal MeshBlock
size (4 cells), not just the test-suite meshes.** The bug was found only because deriving a
wider stencil forced me to ask how many cells were actually available.

Reproduce: `inputs/tests/cubed_sphere_resist_smr.athinput` with `mesh/nx2=8 mesh/nx3=8
meshblock/nx2=4 meshblock/nx3=4`, `time/nlim=0`, and read the GHOST SCAN lines.
