---
name: reflect-wall-mass-leak
description: reflecting radial walls leaked 2-3e-4 mass under a strong blast on sp (always) and cs (since 979edada); FIXED by mirroring the ghost-side wall state in the x1 sweeps
metadata:
  type: project
---

**Bug.** `GridPiecewiseLinearX1` (centroid-position PLM) is not mirror-symmetric about a
reflecting x1 wall: the ghost centroids do not mirror the interior ones, so the two wall
states differ at O(dx^2) and HLLC/HLLD pass a mass flux through a "closed" wall. Strong
blast (blast_p 100, tlim 0.3): sp -2.1e-4 mass; cs -9e-5 (n16) / -3e-4 (n8). It survived
`reconstruct = dc` (x1 Grid-PLM is unconditional on both grids) and HLLE, and scaled like
1/nx1. Invisible in every v_r = 0 test (rigid rotation, toroidal): those were "exact".
Bisected: 13a97399 round-off, 979edada leaks (cs took the sp radial reconstruction).

**Fix (committed 2026-09-06 after 4990eb41).** In `hydro_fluxes.cpp` / `mhd_fluxes.cpp`
x1 sweeps, under `use_spherical_polar || str_r1_`: at an inner/outer `reflect` wall set
wl(is) = mirror of wr(is) (IVX flipped, transverse B copied), and wr(ie+1) likewise. Closed
blasts back to round-off: cs hydro +2.5e-15, cs MHD +1.2e-14 (energy -3e-8 = the
documented O(B^2 h^2) wall leak), sp hydro exactly 0.

**Why:** the user's conservation gates on both grids rest on reflecting walls being closed.
**How to apply:** any "mass drift" on a closed spherical box before this fix is this bug;
re-measure. Two more pre-existing cs_test bugs fixed alongside: iprob 12 total energy
missed B^2/2 (bcc0 unbuilt), and iprob 11 faces were projected, not from a potential.
See [[cs-test-ffdecay-rotaxis]], [[cubed-sphere-seam-conservation]].
