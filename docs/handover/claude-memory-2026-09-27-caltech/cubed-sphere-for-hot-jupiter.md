---
name: cubed-sphere-for-hot-jupiter
description: Whether a cubed sphere is worth it for the dhj campaign, and how hard MHD on it would be
metadata:
  type: project
---

Assessment made 2026-08-27, grounded in the code and in measured dt.

**The timestep argument for a cubed sphere is DEAD.** [[dt-binding-direction]] measured that
the radial stretch made **x1 the binding direction** (dt1 36.2 s vs dt3 63.8 s), so the polar
azimuthal CFL already has 1.76x of slack. Removing the pole buys **nothing** in dt. Do not
re-propose it on those grounds.

**The cell-budget argument is a genuine reallocation.** The current shell is 64 x 128 = 8192
cells, deliberately anisotropic (dtheta 0.84 deg vs dphi 2.81 deg, 3.3:1). A cubed sphere is
quasi-uniform: N=32 gives 6144 cells at 2.8 deg isotropic (0.75x cost); N=64 gives 24576 at
1.4 deg (3.0x cost) -- **2x finer ACROSS the terminator** (the phi direction, where the
day-night gradient is) and 1.7x coarser along it. Worth asking whether the current anisotropy
points the right way for the observable. Estimates from the dt table, not measured.

**The real case is accuracy at the axis, and it is specifically an MHD case.** The polar axis
is a coordinate singularity patched by the phi+pi block-pairing BC; it already forces
`mesh.cpp:155` (polar + MHD requires a single meshblock in r) and is implicated in
[[upper-atm-mottling]]. Near-axis B_phi with CT is where spherical-polar codes classically
fail.

**MHD difficulty, ranked (given what is already in the tree):**
1. **EMF consistency at panel seams -- the crux.** CT holds div B at round-off only if both
   cells sharing an edge use the same EMF x edge length. Across a seam the panels have
   different metrics and rotated reconstruction directions, so they compute *different*
   EMFs. `flux_correct_fc.cpp` is explicitly fine/coarse-only and has **zero** panel
   handling. Needs an explicit seam EMF exchange-and-average.
2. **`mhd_fluxes.cpp` has no gnomonic rotation at all** -- it captures `use_cubed_sphere` at
   line 81 and never uses it, while `hydro_fluxes.cpp` calls `GnomonicEquianglePrimFaceX*`
   and `GnomonicEquiangleFluxX*`. Harder than the hydro version because the rotation mixes
   face-staggered B components.
3. The eight cube corners (three panels meet); every cubed-sphere code special-cases these.
4. AMR across panels: `prolongation.cpp` has no panel handling. `build_tree.cpp` already
   fatals on adaptive + cubed sphere. Not needed here.
5. **The pgen port is the biggest dull cost.** `deep_hot_jupiter_rt.cpp` (~4400 lines) is
   written against **x1 = r**; on the cubed sphere r is **x3**.

**Recommendation given:** not for this campaign. But item 1 is cheap to settle independently
-- advect a field loop across a seam with the existing `cubed_sphere.cpp` MHD branch and
watch div B. Note hydro is NOT yet correct either ([[cubed-sphere-hydro-state]]), so hydro
must be finished first.
