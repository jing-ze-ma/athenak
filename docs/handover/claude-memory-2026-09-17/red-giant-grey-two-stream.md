---
name: red-giant-grey-two-stream
description: "The GREY two-stream rewritten on the correlated-k machinery (bd156d77, problem/rt_grey): it removes the thermal runaway and the runs survive past every ck death point, but the top-cell DRAIN and the emergent-flux decay are still there"
metadata:
  type: project
---

`problem/rt_grey` (commit bd156d77) is a NEW grey sweep built on the correlated-k
scaffolding -- the same layer-integrated two-stream coefficients, the column above the
domain, the thermalised intensity at the tau-blend handover, the semi-implicit source,
the column dump -- with ONE band and the opacity taken straight off the conduction
module's table (`rad_kappa_src = table_rho`, the stellar Rosseland table). It is NOT the
old path: that one is the Parmentier picket fence on the Freedman fit and is still there,
still selected by leaving both `rt_ck` and `rt_grey` off. `rt_ck` + `rt_grey` together is
a fatal error. `ck_nquad` picks the angular quadrature for grey too.

The single opacity is the point: emission and absorption then scale together, kappa
cancels out of a thin layer's balance, and its equilibrium temperature is a fixed point.

## What it fixed (controlled, 6x8x8x320, t = 0 -> 6e4, everything else identical)

|                | correlated-k          | grey                      |
|----------------|-----------------------|---------------------------|
| top cell T     | 3311 -> ~1000 K       | 3311 -> 3622 K, settling  |
| next 8 cells   | 2000-3000, scattered  | 3512-3543, coherent       |
| min density    | falls 27x             | 1.52e-11 -> 1.65e-11      |
| max \|v\|      | 25 km/s at the top    | 1.1 km/s, in the interior |
| emergent flux  | 0.74 L, falling       | 0.31 L, rising            |

Same dt (30.66 s at production resolution), same cost.

## What it did NOT fix -- read this before claiming the problem is solved

Run longer and the TOP STILL DRAINS. `greylong` (6x8x8x320, to t = 6e5, no dt collapse,
dt flat at 58.9 s): min density 1.5e-11 -> 5.2e-13 by t = 5e5, |v| 22 km/s at i = 319,
the outer T sliding 3520 (t = 6e4) -> ~2900 (t = 6e5) and oscillating with the top cell
between 1850 and 5200 K. Emergent flux 4.524e9 -> 4.438e9 erg/cm2/s, i.e. still creeping
DOWN, and only 0.31 L to begin with.

So there are TWO defects, not one. Grey removes the THERMAL runaway (T no longer collapses
to 1000 K and the run no longer dies). The MECHANICAL one is still open and is most likely
the outer boundary: `ox1_bc = user` fills the ghost from the INITIAL column with the
radial velocity reflected, so as the cell below evacuates the ghost keeps holding t = 0
density above it -- a dense ghost over a rarefied cell, which is RT-unstable. An
`open_outer` hydrostatic continuation from the top ACTIVE cell, the way `inner_bc = open`
already works, is the obvious next thing to try.

## Runs
- `grey_prod` job 194160 (96 ranks, 6 nodes, full 6x32x32x320, tlim 3e6): past t = 5.9e5
  with dt 30.65 and zero collapses -- **every ck run died between 3.66e5 and 4.44e5**.
  KE is GROWING there (7.2e4 at t = 2e5 -> 1.37e5 at 4.5e5) where ck's was flat; whether
  that is convection or the top oscillation is not yet decided.
- `greytop` / `greylong` in `/orion/ptmp/jinma/Athenak/red_giant/`: the cheap 6-block
  (6x8x8x320) testbed, ~90 s per 1e3 cycles on 6 tasks. Use it for anything like this.

Related: [[red-giant-top-cooling-runaway]] (the diagnosis), [[correlated-k-design]].
