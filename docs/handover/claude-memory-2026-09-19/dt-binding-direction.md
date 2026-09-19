---
name: dt-binding-direction
description: AthenaK dt is a straight MIN over x1/x2/x3; on the dhj grids the binding direction switches from azimuthal to radial once the radial stretch is on
metadata:
  type: project
---

`dt = cfl * min(dt1, dt2, dt3)` in `src/hydro/hydro_newdt.cpp` -- a **straight MIN over
directions**, each `dtN` itself a min of `dxN/(|vN| + cs)` over every cell. The directions
do NOT add harmonically, so the smallest one wins outright and the others are pure slack.

**MEASURED at cycle 0** on the ck_grav_prod setup (job 11105825, three arms, identical
physics, temporary diagnostic print in `hydro_newdt.cpp`, since reverted):

| grid | dt1 (radial) | dt2 (theta) | dt3 (azimuthal) | binding | dt |
|---|---|---|---|---|---|
| uniform nx1 = 200 | 81.20 | 203.62 | **63.84** | x3 | 19.15 s |
| uniform nx1 = 234 | 69.31 | 203.28 | **63.74** | x3 | 19.12 s |
| stretched nx1 = 234 | **36.22** | 203.41 | 63.78 | **x1** | **10.86 s** |

On the uniform grids dt is set by the AZIMUTHAL CFL in the polar row (63.8 s), exactly as
the `f_stretch_theta` comment in the dhj input says. The polynomial radial stretch
([[radial-grid-stretch]]) cuts dt1 to 36.2 s, below dt3, so **x1 takes over** and dt falls
by 63.84/36.22 = 1.76x -- precisely the observed drop.

dt1 falls only 1.91x while the finest cell is refined 2.74x (dr 4.75e7 -> 1.73e7): the
binding radial cell is NOT the finest one, because `dx1/(|v1|+cs)` trades dr against the
local sound speed.

**Why:** I first explained the drop as the fine cells "competing with" the polar CFL, which
is wrong for a min -- and the design had predicted dt would be 22.0 s, i.e. it did not
anticipate the switch at all. Both were fixed by measuring.

**How to apply:**
* The dt penalty is the **unavoidable price of the 10 cells/H requirement**, not a
  misconfiguration. Backing the stretch off until dt1 returns to 63.8 s needs dr_min ~3.3e7,
  which gives only ~5.3 cells/H. No free operating point exists between the two.
* **`f_stretch_theta = 3.0` no longer buys dt** on the stretched grid -- x3 has 1.76x slack.
  Keep it for the equatorial resolution it gives the terminators, not for its dt argument.
* **Angular refinement is now free in dt** until dt3 falls to 36.2 s: `nx3` could roughly
  double from 128 before azimuthal binds again. Free in dt, not in cell count.
* When a dt surprises you on this code, print dt1/dt2/dt3 -- one rebuild and a 20 s apudev
  job answered it. Arms live in `bench/ck_grav_prod/dtdiag/`.

See [[ck-grav-prod-run]], [[radial-grid-stretch]], [[dhj-grid-resolution-design]].

## UPDATE 2026-09-01: in the MHD production run it is NOT the CFL at all

On `cs_dhj_long` (cubed sphere, nx1 = 128 + the refitted stretch, `ohmic_resistivity =
eos`, `max_eta = 1e13`) dt falls from the t=0 advective 28.98 s and then **pins at
1.044661e+01 to seven digits for 15 000+ cycles**.  A flow-limited dt fluctuates; a
constant one is a FIXED limit.  It is the OHMIC term at the eta CAP:

```
  dt = cfl * fac * dr_min^2 / max_eta,   fac = 1/6 in 3D (resistivity.cpp),
       cfl applied at mesh.cpp:894 -- dt = min(dt, cfl_no * presist->dtnew)
     = 0.3 * (1/6) * (4.5709e7)^2 / 1e13 = 10.4466 s      observed 1.044661e+01
```

dr_min is the FINEST RADIAL cell (r/ap = 1.171, the scale-height minimum the stretch
deliberately refines), and eta is pinned at its cap there.

**CONFIRMED by a two-arm apudev A/B (job 11301392), both predictions exact:**
  - `max_eta` 1e13 -> 1e12: dt = **1.044661e+02**, exactly 10x. Scales as 1/eta.
  - the ORIGINAL stretch coefficients (dr_min 3.165e7): dt = **5.007167**, against a
    prediction of 5.0069.  Scales as dr_min^2.

**So the refit's dt gain is 2.09x, not the 1.46x I first measured** -- the 1.46x was the
EARLY ADVECTIVE dt; the asymptotic limit goes as dr_min^2 (see [[radial-stretch-refit]]).
And the knobs on dt here are `max_eta`, `use_rkg_sts`, or the finest cell -- NOT
`cfl_number`, and not the angular grid.

