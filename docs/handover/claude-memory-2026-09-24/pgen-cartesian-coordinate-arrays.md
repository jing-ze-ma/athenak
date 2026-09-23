---
name: pgen-cartesian-coordinate-arrays
description: "pcoord's geometry arrays are unallocated 1x1 placeholders on Cartesian meshes; the audit of solar_convection/cooling_convection for unguarded uses, what was fixed and what is still latent"
metadata: 
  node_type: memory
  type: project
  originSessionId: 0955e70d-6014-446b-aa7c-a4f1fadba65d
  modified: 2026-08-15T15:23:08.575Z
---

Written 2026-08-15 after [[sponge-inert-at-192]].

`Coordinates::Coordinates` (`src/coordinates/coordinates.cpp`) allocates x1v/x2v/x3v,
xx1f/xx2f/xx3f, area.*, volume, dx1/dx2/dx3, dxedge, dxface, areaedge, x_ov_rD/y_ov_rC/
z_ov_rE **only** when `use_cubed_sphere || use_spherical_polar` (and sin/cos_cell only for
cubed-sphere). Otherwise they stay the 1x1 Views from the init list. Indexing them on a
Cartesian mesh is an out-of-bounds read that Release builds do not catch. On Cartesian use
`size.d_view(m).x1min/x1max/dx1` + `CellCenterX` instead.

## Audit of solar_convection.cpp / cooling_convection.cpp (both files, same structure)

FIXED:
1. The sponge kernel's `x1v_(m,i)` -- solar_convection only. See [[sponge-inert-at-192]].
2. `Real area_r, area_l, vol;` in `SourceFunc` are assigned ONLY in the
   `if (use_spherical_polar)` branch but used unconditionally at
   `src = bdt*(area_r*(pr-p)+area_l*(p-pl))/vol` when `use_wellbalance_dynamic` is on.
   In Cartesian they were uninitialized STACK locals. Now set to unit faces and
   `vol = size.d_view(m).dx1`, which is the intended `(pr-pl)/dz`. Both files. Not
   triggered by any run so far (the sun runs use `etotgrav`, not wellbalance_dynamic).

STILL LATENT, deliberately not changed (needs a decision):
3. The two MHD user-boundary kernels -- `usrboundaryx1_bfield` and
   `usrboundaryx1_bfieldc` -- are guarded by `pmbp->pmhd != nullptr` ONLY, and index
   `x1f_`, `x1v_`, `x2f_`, `area1/2/3`, `volume`, `z_ov_rE`. They also use spherical
   formulas (`SQR(x1f)` ratios). A Cartesian MHD run would read out of bounds AND get the
   wrong geometry. It is already unusable in that configuration for a separate reason: the
   `pgen_b0` field initialisation is entirely inside `if (use_spherical_polar)`, so `b0` is
   never initialised on a Cartesian mesh. Cleanest fix is a fatal error at setup
   ("MHD requires a spherical-polar mesh in this pgen"), which is a behaviour change.

Clean: the two_stream_RT `dx1` uses ARE properly guarded (`dr = dx1(...)` vs
`size.d_view(m).dx1`).
