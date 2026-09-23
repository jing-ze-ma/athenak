---
name: cubed-sphere-panel-frames
description: The six cubed-sphere panel frames are FORCED by Mesh::panel_neighbors -- panels 3 and 4 are not what they look like, and a duplicate copy had them swapped
metadata:
  type: reference
---

Established 2026-08-27 by propagating `Mesh::panel_neighbors` around the cube. The result
is unique, self-consistent (six distinct normals, all proper rotations, all 24 directed
edges agree), and independently reproduces every entry of `Mesh::GetPanelBoundary`.

Frame `(a, b, n)`, point = `(a*tan(xi) + b*tan(eta) + n)/delta`, n = outward face normal:

| panel | a | b | n |
|---|---|---|---|
| 0 | (0,1,0) | (0,0,1) | (+1,0,0) |
| 1 | (-1,0,0) | (0,0,1) | (0,+1,0) |
| 2 | (0,-1,0) | (0,0,1) | (-1,0,0) |
| 3 | (0,1,0) | (-1,0,0) | (0,0,+1) |
| 4 | (+1,0,0) | (0,0,1) | (0,-1,0) |
| 5 | (0,1,0) | (+1,0,0) | (0,0,-1) |

So **4,0,1,2 walk in longitude and 3 (+z) / 5 (-z) are the two poles**, matching the face
diagram in mesh.hpp. This is counterintuitive: 3 and 4 are not the adjacent pair the
numbering suggests.

**The trap.** `src/pgen/cs_test.cpp` carried an independent copy with **panels 3 and 4
interchanged**. A wrong panel-to-cube assignment is invisible to every single-panel test,
because `CoordGnomonicEquiangle` has no notion of global orientation -- every panel's
metric is the same function of (xi,eta). It only shows up in a GLOBALLY defined field: the
`iprob=3` rigid rotation was a different field on two of six panels, which put an O(1)
mismatch on all their seams. **Every rigid-rotation number recorded before this date is
contaminated**, including the "does not converge" claim in [[cubed-sphere-hydro-state]] and
the 8c7040b0 harness numbers in [[cubed-sphere-hydro-fix]].

The frames now live ONCE, in `src/coordinates/cubed_sphere.hpp`; cs_test.cpp wraps them.
Do not write a second copy.

**The tests that found it**, both new in cs_test.cpp and worth keeping:
* `iprob=6` -- rho = f(cartesian direction), uniform p, v = 0. Exactly steady for ANY rho,
  but NOT invariant under an angular remapping, unlike `iprob=5`. Reports per-panel,
  per-face max ghost error. Flat = wrong panel/orientation; halving = offset only.
  Keep the profile manifestly POSITIVE -- a negative rho is silently floored and reads as
  a fake O(1) active-cell error, which cost a detour.
* `iprob=7` -- rho = 1000*panel + j + k/1000, a TAG. Decoding a ghost names the exact
  (panel, j, k) it was fetched from. This is what proved the halo topology and orientation
  are correct on all 24 faces, so the remaining error had to be geometric, not
  combinatorial.
