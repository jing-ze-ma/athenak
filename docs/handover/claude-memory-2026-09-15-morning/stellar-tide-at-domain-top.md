---
name: stellar-tide-at-domain-top
description: The host star's tidal term is missing from the pgen; measured effect is 0.5% at the limb but +6.8 cells at the substellar column on the taller grid
metadata:
  type: project
---

**Computed 2026-08-26.** In the co-rotating frame the full tide+centrifugal acceleration is
`(3 Om^2 x, 0, -Om^2 z)`, x toward the star, z along the spin axis. `deep_hot_jupiter_rt`
ALREADY applies the spin centrifugal part `Om^2 (x, y, 0)` -- that is `oor = SQR(omega)*r*sine`
in the co-rotating source term, ~line 2765. (The "centrifugal term, never enabled" comment at
line ~203 refers to the term in the hydrostatic POTENTIAL, not the momentum source.)
**Missing is the stellar part, `(2 Om^2 x, -Om^2 y, -Om^2 z)`.**

**No stellar mass needed:** the planet is synchronous, so Om is `problem/omega` itself. The
balance point 3 Om^2 r = GM_p/r^2 gives **r_L1 = 4.04e10 cm = 4.28 R_p** from `grav` and `ap`
alone. The taller domain top (1.90e10) is at 0.47 r_L1 -- well inside the Roche lobe, so this
is a correction, not an overflow problem.

| missing term / local g | substellar & antistellar | limb | pole |
|---|---|---|---|
| ck_limb top, 1.433e10 | 1.3 % | 0.6 % | 0.6 % |
| ck_grav_size top, 1.90e10 | 6.9 % | 3.5 % | 3.5 % |

**Effect on r(1e-6 bar)**, integrating hydrostatic balance from 1e-2 bar with the run's own
T(p) and mu (so hydrostatic response only, no dynamical feedback -- fair at 2 %):

| | substellar | leading limb | trailing limb | antistellar |
|---|---|---|---|---|
| ck_limb | +0.23 % (+0.6 cells) | -0.07 % | -0.06 % | +0.05 % |
| ck_grav_size | **+2.07 % (+6.8 cells)** | -0.53 % (-1.5 cells) | -0.38 % | +0.22 % |

**Verdict: not important for the limb** -- half a percent at most, and it moves DOWN (at the
terminator the code applies an unopposed spin-centrifugal push the stellar tide would cancel).
**It shows up at the substellar column**, +6.8 cells on the point-mass grid, which is already
the column with the least headroom ([[grav-point-mass-flag]], [[ck-grav-size-run]]). Size the
production grid against the LIMB as planned and this stays a footnote; size it against the
substellar level with a few cells to spare and the tide is the same size as that margin.

If it is ever implemented, it goes at the same call sites as `grav_point_mass` (momentum
source + the two hydrostatic integrators + the outer-BC ghost extrapolation).
