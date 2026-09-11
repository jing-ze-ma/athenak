---
name: sponge-inert-at-192
description: "RESOLVED: solar_convection's sponge read pcoord->x1v, which is a 1x1 placeholder View on Cartesian meshes; the out-of-bounds read is why the layer damped at 64^2 and was inert at 192^2"
metadata: 
  node_type: memory
  type: project
  originSessionId: 0955e70d-6014-446b-aa7c-a4f1fadba65d
  modified: 2026-08-15T15:22:53.810Z
---

RESOLVED 2026-08-15, same day it was opened. `problem/sponge = true` damped at 96x64^2 and
produced BIT-IDENTICAL output to no-sponge at 96x192^2.

## Root cause

The sponge kernel took the cell height as `x1v_(m,i)`, i.e. `pcoord->x1v`. That array is
**only allocated and filled for cubed-sphere / spherical-polar meshes**
(`src/coordinates/coordinates.cpp`, the `if (use_cubed_sphere || use_spherical_polar)`
block). On a Cartesian mesh it is still the **1x1 placeholder View** from the constructor's
init list, so `x1v_(m,i)` read past its allocation and returned adjacent heap memory. The
parameter was read, the kernel ran, the writes were fine -- the `x1v > zs` test just
compared against garbage. At 192^2 no garbage exceeded 1.6e8, so not one cell was touched.
At 64^2 some did, so the layer "worked" at whatever cells the neighbouring heap happened to
select -- NOT the top 20%.

Every other use of `x1v_` in the pgen, including the main source loop 30 lines above, is
guarded by `if (use_spherical_polar)` and otherwise uses
`CellCenterX(i-is, indcs.nx1, x1min, x1max)`. The sponge was the one unguarded use.
**Lesson: on a Cartesian mesh in this codebase, ALL of `pcoord`'s geometry arrays
(x1v/x2v/x3v, xx?f, area, volume, dx1/2/3, dxedge, dxface, areaedge, z_ov_rE, sin/cos_cell)
are unallocated 1x1 placeholders. Only `mb_size` (`size.d_view(m).x1min/x1max/dx1`) is
valid.**

## Verified

Same binary, 96x192^2, `sponge` true vs false, 100 s (`soleos/spongefix/{on,off}`): rms v_z
ratio is exactly 1.000 below z/zmax = 0.8, then ramps quadratically to **0.393 at the top
cell** after only 100 s. Density unchanged to 5e-11. That is the intended layer.

## Consequences

- **The 64^2 `sponge08` / `sponge07` results are invalid** and their conclusions must be
  re-derived. See [[solar-convection-general-eos]].
- Two more instances of the same bug class were found in the audit -- see
  [[pgen-cartesian-coordinate-arrays]].
