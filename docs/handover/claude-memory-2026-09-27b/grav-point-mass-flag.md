---
name: grav-point-mass-flag
description: problem/grav_point_mass implements 1/r^2 gravity in deep_hot_jupiter_rt; verified bitwise-off and well-balanced-on, NOT yet committed
metadata:
  type: project
---

**Committed as c37ebe75** (2026-08-26) on `polar-average-perf`. Default is OFF, and the
off path is bitwise identical to f882159f, so nothing changes for existing inputs.
Pushed to the fork 2026-08-26.

**Why:** `<problem>/grav` was applied as a CONSTANT g over a domain spanning r/ap = 1.52,
where a point mass with the same g0 gives 0.43*g0. That compresses the upper atmosphere: at
1e-6 bar the level sits ~5 scale heights too low and H is understated by **2.7x**, which
propagates straight into transmission radii and feature amplitudes -- a first-order
systematic on the limb observable this campaign exists for. See [[ck-limb-run]].

**What:** two helpers, `GravAccAt` (g(r) = g0 (ap/r)^2) and `GravPotAt`
(phi = g0 ap (1 - ap/r)), both reducing to the constant-g forms as r -> ap so the flag is a
no-op at the inner boundary. Applied at all 18 gravity consumers: 10 potential sites
(phicc0 + three phi0 face arrays), the 2 hydrostatic integrators, 2 momentum-source sites,
the Maxwell outer-BC clamp, and **3 RT sites** where tau = kappa p/g closes the unresolved
column above the domain (that column is 2.3x more opaque than the constant-g form gave).

**Verified on apudev (jobs 11025710 off / 11025711 on, production grid, 0.5 rotation):**

* **Flag OFF is BITWISE identical** to a reference binary built from a clean checkout, all
  five arrays. (Raw files differ only in the header, which records the new parameter.)
* **Flag ON gives the right gravity:** -(1/rho)dp/dr across the IC = 927.7 base -> 594.1 mid
  -> **414.5 top vs the analytic 413.6**.
* **Well-balanced property HOLDS.** Deep-interior (r < 1.05e10) rms v_r DECAYS in both arms
  and converges to the same value -- 207 cm/s (off) vs 203 cm/s (on). A missed face array
  would have shown growth. `eos_fail = c2p_it = fofc = 0` in both.

**Two things that fell out, both favourable:**

* **dt is 2.3x LARGER with the point mass** (16.2 s vs 6.96 s at t = 1.5e5), so the correct
  physics is CHEAPER, not more expensive.
* **`eos_dfloor` fires 128x less** (4.75M vs 606M) and `eos_efloor` 42x less. The floor
  activity in the constant-g runs was largely an artefact of artificial compression.

**Watch:** mass drift is HIGHER with the flag on (+0.196 % vs +0.115 % over 0.5 rot) despite
far fewer floor firings, which points at the outer boundary rather than the floors -- and is
consistent with the true 1e-6 bar limb level (~1.55e10) sitting OUTSIDE the current
x1max = 1.433e10. **Turning this flag on requires resizing the domain again**: first-order
x1max ~ 1.64e10 with nx1 ~ 125 to hold dr.

## SIZING: SETTLED 2026-08-26. Both estimates below are DEAD -- see [[ck-grav-size-run]]

**The measured answer is x1max = 2.0556e10, nx1 = 200** ([[ck-grav-prod-run]]). Everything
in the rest of this section is superseded and kept only to show why both estimates failed.

## (superseded) SIZING THE FLAG-ON CAMPAIGN -- do not use these numbers

I estimated x1max ~ 1.64e10 / nx1 ~ 125 by integrating hydrostatic balance with the
CONSTANT-g run's T(p) held fixed. Measuring the actual flag-on arm instead gives a much
larger answer -- terminator top-cell pressure 6.5e-5 bar, isobar outside for 100 % of
columns, implying **x1max ~ 1.87e10 (2H) to 1.92e10 (3H), nx1 = 166-177**.

**But that arm is only 0.49 rotation from the IC, so neither number is trustworthy.**
[[dhj-floors-for-1e-6-bar]] already records this exact trap: measure the EVOLVED state, not
the IC. The old constant-g run's top-cell pressure fell from 7.5e-5 to 3.8e-6 bar as it
relaxed over rotations, so the flag-on requirement will shrink substantially too.

**How to actually do it:** run flag-on on a deliberately OVERSIZED grid (x1max = 1.9e10,
nx1 = 177 holds dr at 5.56e7) for several rotations, measure where the 1e-6 bar isobar
settles at lon = +-90, then trim x1max and nx1 to that plus ~3 scale heights. Cost is
roughly neutral despite ~2x the cells, because **dt is 2.3x larger with the point mass**.
