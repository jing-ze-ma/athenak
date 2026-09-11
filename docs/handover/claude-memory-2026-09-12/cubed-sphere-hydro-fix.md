---
name: cubed-sphere-hydro-fix
description: Cubed-sphere (gnomonic equiangle) hydro - halo axis bug FIXED; remaining error diagnosed as the panel halo needing along-seam interpolation
metadata:
  type: project
---

Started 2026-08-27 (session 8c7040b0, killed by a network drop; resumed in 5a2885e1).
The user asked how useful a cubed sphere would be for the hot-Jupiter atmosphere and how
hard MHD would be, then: **"try to fix cubed sphere for hydro first. that doesn't work
yet"**. All work is UNCOMMITTED in the tree.

**BUG 1, FIXED: the panel halo exchange used the wrong two axes.**
`src/bvals/bvals_cc.cpp` applied the `PanelBoundaries` swap/reversal to `j`/`k` and
remapped `IVY`/`IVZ`, but the panel-tangential pair is `x1`/`x2` (`Mesh::panel_neighbors`
numbers faces `0:-x1,1:+x1,2:-x2,3:+x2`; `CoordGnomonicEquiangle` puts xi on x1, eta on x2,
RADIUS on x3). It was reversing the RADIAL index across a seam. Now acts on `i`/`j` and
`IVX`/`IVY`, with the buffer transpose written as
`index = si*(i-il) + sj*(j-jl) + ni*nj*((k-kl) + nk*v)`, `(si,sj) = (1,ni)` or `(nj,1)`.
Likely origin: someone part-way between conventions -- **the user wants x1 to be RADIAL
later, for compatibility with the spherical grid** ([[cubed-sphere-x1-radial]]).
`src/bvals/bvals_fc.cpp` has the same pattern and is NOT yet fixed (MHD). **NOTE: after
[[cubed-sphere-x1-radial]] the j/k + IVY/IVZ form is the CORRECT one, and bvals_cc.cpp was
reverted to it -- the original author had written bvals for x1-radial all along.**

**The test that found it (new `iprob = 5` in `src/pgen/cs_test.cpp`):** rho = rho(r) only,
uniform p, v = 0 -- an exact steady state, and invariant under any ANGULAR remapping, so
only radial-index handling can break it. Pre-fix the ghost error was exactly 0.4375 = the
full spread of the profile (a completely reversed radial column) on every panel whose
`GetPanelBoundary` entry has a rev/swap flag; post-fix it is 0.0 everywhere.

**BUG 2 IS SUBLEADING -- see [[cubed-sphere-seam-basis]] (2026-08-27).** The dominant seam
defect was the tangent-BASIS transform (O(1), flat in resolution), not the offset below
(O(dx)). It is now fixed and the rigid rotation converges. Every measured number in this
memory predates [[cubed-sphere-panel-frames]] and is contaminated by the panel 3/4 swap.

**BUG 2, DIAGNOSED NOT FIXED: the halo needs along-seam INTERPOLATION.**
Matching `(1,tan xi,tan eta)` on panel 0 to `(-tan xi',1,tan eta')` on panel 1 gives
  `tan(xi') = -1/tan(xi)`  -> `xi' = -pi/4 + delta` EXACTLY (normal direction is right)
  `tan(eta') = tan(eta)/tan(xi)` -> **`eta' != eta`**, which the copy assumes.
So the ghost's own geometry and the data copied into it are different physical points. The
offset is a pure shift ALONG the seam (an exact match always exists, residual ~1e-9), is a
smooth function of eta, and in units of cells is nearly resolution-INdependent: layer g=0
runs 0 at the seam midline to -0.5 at the panel corner; layer g=1 reaches -1.47. Verified
analytic formula `(atan(tan(eta)/tan(pi/4+(g+0.5)dxi)) - eta)/dxi` matches a numerical
match to all digits. A ghost value wrong by O(dx) makes a flux difference wrong by O(dx),
divided by the cell width -> an **O(1) acceleration error**. Fix = interpolate the ghost
layers along the seam-parallel index after the exchange (>=3-cell stencil for layer 2).

**Measured state after BUG 1's fix** (all values unchanged by the x1-radial refactor;
note the radial error is now the 1-KE history column, not 3-KE) (`bcs/` build in the 8c7040b0 scratchpad):
* Test A, 3D uniform static (`run2`): exact, matches the pre-existing reference.
* `iprob=5` radial halo: exact on all six panels.
* Rigid rotation (`iprob=3`, exact solution, so `3-KE` is pure error), t=1, angular
  refinement 16/32/64: dM/M `+1.51e-3, +1.30e-3, +6.21e-4` (roughly first order), but
  `3-KE` `7.10e-3, 7.64e-3, 7.90e-3` -- still not converging. Refining nx3 jointly
  (8/16/32) changes nothing: `7.10e-3, 7.91e-3, 8.43e-3`.
* **One-step residual is the diagnostic that works**: `|dv3/dt|` = `0.2097, 0.2139,
  0.2164` -> a nonzero constant, and split by region rms`|v3|/dt` is `1.6e-2` at seams vs
  `8.1e-5` in the interior, a factor 200, both flat in resolution. The interior 8e-5 also
  fails to converge -- a second, weaker effect, not yet chased.

**Process lesson worth keeping.** I first localised the error at t=1 and concluded "bulk,
not seam" -- WRONG: by t=1 the seam error has propagated and filled the domain. Localise
the SOURCE with a one-step residual, never with an evolved snapshot. Cf.
[[measure-impact-before-claiming]].

Also: `run1` (the 2D nx3=1 shell with periodic radial BCs) has never passed -- it produced
NaN before any of this work and now gives a spurious radial momentum. Uniform p in a
finite-thickness shell has a real net outward force; the 2D case needs its own thought.
The `conv16/32/64` numbers a second session produced at 13:45 were from a "variant B"
whose own Test A was failing; ignore them.

Harness: `/tmp/claude-28895/-viper-u2-jinma-ATHENAK-athenak/8c7040b0-b255-41e2-ac10-f6dddbad1db0/scratchpad`
-- `bcs/` (build), `run1/2/3`, `conv16/32/64` (angular refinement), `c3d*` (joint),
`res*` (one-step residual), `loc*`, `ghost5/`. Rebuild there, never in `run/`
([[never-write-in-run-dir]]).
