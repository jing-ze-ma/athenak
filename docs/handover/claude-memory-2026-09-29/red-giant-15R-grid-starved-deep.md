---
name: red-giant-15R-grid-starved-deep
description: "The 1.5 R stretch refit (c1..c4 = +2.301696 -4.096725 +3.607269 -2.841427, nx1 320) leaves 0.92 cells/H_p at the inner wall (1.1 R grid: 5.2) and the star dies at t=3.2e5 from a wall inflow -- NOT the ambient medium (which held rho to 3 % and free-fell as designed). fit_radial_stretch.py's smoothing (~200 pts) hides this: measure cells/H_p on the ACTUAL grid. Need a physics-weighted refit and more radial cells."
metadata:
  type: project
---

R3_full15 (2026-09-09 05:30, 8x8, full production config, no MLT, seed cut, nquad 1):
`dt COLLAPSE cycle=1592 t=3.22e5 cond=4e-10` at cell i=0-6 (r 0.08 R): rho at the wall
x2.9 by 3e5, T(i=6) 1.29e6 -> 2.03e6 then off the EOS table (2.3e10 K, kappa floor,
kappa_rad 1.7e27).  An inflow at the inner wall from cycle 1.

| grid | cells/H_p at i=0 | i=2 | i=5 | i=10 | median |
|---|---|---|---|---|---|
| 1.1 R fit (R2a) | 5.25 | 5.27 | 5.40 | 5.80 | 12.3 |
| 1.5 R fit (R3) | **0.92** | 1.68 | 2.73 | 4.30 | 8.45 |

Delta r = 4.8e10 cm vs H_p = 4.5e10 at the base.  The 320 cells that covered 19 scale
heights now cover ~40 (the star's thin atmosphere from 1.1 to 1.216 R adds ~20, the
ambient ~1.3), and the fit put the loss at the deep end.  The fitter's report ("min 2.20")
was wrong because it smooths H over ~200 fine points.  The input's stretch comment
("5.0 / 10.1 / 18.2") is the OLD grid's numbers.

**Fix (in progress):** a physics-weighted target -- >= 5 cells/H_p everywhere the star
is (r < ~1.1 R), 2-3 in the thin atmosphere above tau ~ 1e-2 (nothing to resolve there),
geometric dr/r = 0.2 in the ambient -- and the smallest nx1 that achieves it (expect
384-448), verified on the actual grid by evaluating the stretch map, not from the fitter.

Also settled the same night: ck_nquad = 2 has NO defect (R2c == R0f to 3 digits; both
Mach 2.5 at r/R 1.074 at 1.4e6 at 8x8 -- a transonic thin atmosphere, an 8x8 pathology
production at 32x32 did not show); the seed cut leaves the deep quiet (v_rms/v_mlt 30-200x
below baseline inside 0.3 R); R2d (MLT chi limiter removed) was worse, MLT stays off.
Related: [[red-giant-seed-across-rcb-proven]], [[red-giant-wb-kills-ambient-medium]].

## R4 (no seed) and the join, 2026-09-09 05:50

- **No seed => no surface convection in 4e5 s**: v_rms at i=140-160 stays 1e-3..3e-2 cm/s
  (noise, no trend) against v_mlt 1.8e5.  From truncation-level amplitude the surface needs
  ~1.5e6 s; the seed (outer envelope only, vpert_rmin) stays.  R3 (seeded) reached
  v_rms/v_mlt 0.3 at t = 2e5.
- **The under-pressured join drives a rarefaction into the star's thin atmosphere.**  The
  background was density-matched at 2e-22 but at 400 K vs the 3364 K skin (8x lower p):
  a rarefaction ran down at c_s, T at 1.08-1.10 R fell to 2200-2400 K, rho at i=175 was
  stripped 17x by 4e5 s (R4), and the 34 km/s free fall shocked against the outflow
  (T 8.8e6 K in shocked thin gas -> dt collapse).  Ram pressure 2.4e-9 dyn/cm^2 is 1e5
  below the thermal pressure at 1.08 R, so it is NOT ram stripping.  In R3 the seeded
  surface flow held the front at ~1.10 R -- a race, not a solution.  Fix: bg_temp = 3364 K
  (pressure-matched join; also the RE temperature of thin gas in the star's light; Bondi
  rho_max 1.1e-20 >> 2e-22) and dfloor 1e-23.  The ambient itself held rho to 3 % and
  free-fell at -3.4e6 cm/s in both runs -- the medium is fine, the JOIN was wrong.

## New grid (agent refit, 2026-09-09 06:10): nx1 = 480
f_stretch_r_c1..c4 = -0.402945 +10.179232 -21.570641 +11.840082 (physics-weighted target:
H_true for rho > 1e-11, 3 H in the thin atmosphere, dr/r = 0.2 in the ambient; verified on
the actual grid).  R5_seed (8x8, bg_temp 3364, dfloor 1e-23, nquad 2, seed cut, no MLT):
dt FLAT at 52.35 s (R3 had 37), wall held, L_rad,out/L 0.78-0.89 -- then an ABRUPT
one-step hydro collapse at t = 5.3720e5 (cycle 10261, rank 3 = panel 3: 52.3 -> 1.3e-7),
i.e. a single cell went bad, not a grind.  Being located (dump 5 at 5e5).  Staged into
prod15/rg.athinput (nx1 480 both blocks + the four c's) -- NOT submitted.
