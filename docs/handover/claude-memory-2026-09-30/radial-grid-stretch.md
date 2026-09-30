---
name: radial-grid-stretch
description: mesh/use_grid_stretch_r_poly -- non-monotonic radial stretch; new coordinates/grid_stretch.hpp; UNCOMMITTED
metadata:
  type: project
---

Written 2026-08-26 (ninth session). **Committed as 1e19d4e7** on `polar-average-perf`, on
top of [[stellar-tide-flag]] (c4aa730f). Not yet pushed. 125 lines across 5 files + the new
`src/coordinates/grid_stretch.hpp`.

**What.** `mesh/use_grid_stretch_r_poly` + `f_stretch_r_c1..c4`, mapping
`r = r0 + (r1-r0) u(xi)`, `u(xi) = xi + sum_k c_k xi^k (1-xi)`. Endpoints fixed for ANY
coefficients; `c = 0` is uniform. The pre-existing `f_stretch_r` (exponential) is monotonic
in cell width and CANNOT refine a mid-domain scale-height minimum while coarsening aloft;
this family can. `Mesh`'s ctor samples `du/dxi` over [0,1] and FATALS if not strictly
positive (a fold gives negative cell widths).

**Three non-obvious things fixed while doing it:**
1. **The pgen computes radii BOTH ways** -- 44 reads of `x1v_` (stretched by Coordinates)
   and 13 direct `CellCenterX/LeftEdgeX` calls (UNSTRETCHED). The old
   `deep_hot_jupiter_rt_old.cpp` called `StretchR` after every one; the current pgen dropped
   them because radial stretching was never used with it. Left alone the coordinates and the
   IC sit on different grids, silently. All 13 now route through `ApplyRStretch`.
2. `StretchR/StretchTheta` were **members of Coordinates**, so calling them in a
   KOKKOS_LAMBDA captured `this` (a host pointer -- works only under unified memory). Moved
   to file scope in a new dependency-free **`src/coordinates/grid_stretch.hpp`**, shared by
   Coordinates and the pgen so they cannot drift. `NSTRETCH_R_POLY` (=4) lives there.
   Headers are not listed in `src/CMakeLists.txt`, so no build change needed.
3. Mesh params are now hoisted into locals before the `par_for` in `CoordSphericalPolar`.

**POST-PROCESSING TRAP.** `bin_convert`/`athena_read` rebuild coordinates as
`np.linspace(xmin,xmax,nx+1)` -- the binary format does not store stretched coordinates.
This is ALREADY true for theta (why every script applies the sinh map to `x2v` by hand); a
stretched radial grid extends it to `x1v`. **Every analysis script must apply the polynomial
map or silently mis-locate every radius.**

**Verified end-to-end** by dumping the actual radii from a run (temporary print, removed):
dr 7.55e7 at the base -> 3.31e7 at r/Rp 1.19 -> 2.41e8 at the top, 7.7:1, non-monotonic as
designed. Style-clean (0 cpplint violations on added lines). CPU + GPU builds green.
A 0.25-rot A/B (uniform 200 vs poly 120) agreed within the run's own temporal variability at
every pressure level, and measured 0.573x cost against a predicted 0.60x.
