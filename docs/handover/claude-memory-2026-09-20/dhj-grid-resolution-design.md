---
name: dhj-grid-resolution-design
description: How to size the dhj radial grid on resolution rather than extent -- and the metric mistakes that make it look easy
metadata:
  type: feedback
---

The 2026-08-26 (ninth session) design study behind the N=234 grid in [[ck-grav-prod-run]].

**The user's requirement, stated plainly: 10 cells per scale height BELOW the 1e-6 bar
level.** No constraint above it.

**MEASURE AT THE FEATURE, NOT AGGREGATED.** My first pass reported "terminator p01 2.49
cells/H" for the uniform-200 grid and I recommended a 0.55x cheaper grid on the strength of
it. Wrong: that aggregate mixes two different regions.
* AT the 1e-6 bar observable isobar, uniform-200 gives **p01 9.7, median 17.2** -- it already
  met the standard there.
* The under-resolved place is the H2 dissociation structure much deeper, at p ~ 1e-3..1e-2
  bar.
Equalizing `dr ~ H` therefore ROBS the observable to feed the front. The "40% saving" was an
artifact of the bad metric.

**Two distinct features, ~6 cells apart -- I conflated them at first:**
* **scale-height minimum**, r/Rp = 1.244, p = 2.0e-2 bar, T = 2352 K (coldest), mu = 2.219
  (fully molecular), H = 1.46e8. Genuinely hydrostatic (H/H_hse = 1.01). **This is what sets
  the cost**: 10 cells/H here forces dr = 1.5e7 where the signal speed is high.
* **H2 front** (mu 2.22 -> 1.25), r/Rp = 1.280, p = 2.5e-3 bar. Thickness **2.07e8 cm** (thin
  10%), 3.08e8 median. **It is RESOLVED**: measured thickness changed 1.02x while dr changed
  1.30x between two grids, so it is physical, not grid-limited. Resolving it to 10 cells
  costs only **+11 cells (+0.09x)** -- it is 4.2% of the domain.

**The cost is the TIMESTEP, not the cells.** N 200 -> 234 is 1.07x; dt 32.0 -> 22.0 s is the
rest. Know which CFL limits: on the production grid x3 (azimuthal) limits at 32.0 s, x1 at
41.5 s, x2 at 125.7 s. Radial refinement is free until dt_x1 drops below 32.

**A cruder fit can be CHEAPER and still meet the target.** The degree-4 polynomial hit median
10.1 cells/H at 1.70x, while the idealized profile it was fitted to cost 2.25x -- the
smoothing stops it over-refining the scale-height minimum. Do not assume more coefficients is
better.

**Withdrawn:** I suggested rebalancing `f_stretch_theta` once x1 becomes the limiter. It
cannot help -- theta stretching changes dt, not cell count, and dt is then set by x1 with x3
slack simply unused.
