---
name: red-giant-dust-opacity-kills-lidfree
description: PROVEN 2026-09-09 — the AESOPUS grain-condensation step in the opacity table kills every lid-free red-giant run; problem/opac_tmin=2500 (grains off) survives; fully explicit RT is a negative control
metadata:
  type: project
---

All lid-free 32x32 runs (R9, T1 wall, T2 old semi-implicit, T3 no-Newton, T4 top-ghost RE off)
died at t = 1.9-4.4e5 with the top cooled to 1900-2500 K. Cause: the merged opacity table is
AESOPUS WITH GRAINS; kappa_R at 1e-11 g/cm3 rises 2e-5 -> 1.5e-2 cm2/g between 2600 and 2000 K
(onset 2500-2700 K at photospheric 1e-9), a factor 5-10 per 0.025-dex table row. A cooling,
draining top cell gets radiatively STIFFER instead of more transparent.

- T6_nodust = R9 config + `problem/opac_tmin = 2500` (pgen clamps table rows below it to the
  2500 K row): reached tlim 5e5 with dt 30.65 flat, top cell 1857 K / 6e-13 (colder and thinner
  than R9 ever got) -> the cooling itself is not fatal. Extended to 3e6 as job 194594.
- T5_explicit (`problem/rt_explicit`, de = src*dt): died EARLIER (1.69e5) at tau~18 (i=291,
  T -> 3e10 K in one step). Explicit RT is unstable in the thick layer; the semi-implicit
  update is NOT the problem. Don't bisect the RT update again.
- The clamp does not fix L_out/L (0.35-0.42): that is the spin-up deficit (see
  [[red-giant-flux-deficit-is-spinup]]). Physics: skin T = Teff/2^(1/4) = 3360 K sits within
  ~600 K of condensation; the transient (L_out/L 0.4 -> photosphere 3200 K, skin 2650 K) trips
  it. A hot corona above the star cannot change this ([[red-giant-grey-opaque-lid-bug]]).

**How to apply:** any lid-free run needs opac_tmin ~ 2500 (or a grain-free table) until the
star is at steady state; dust-driven outflow is a separate project. Production candidate
prod12 = lid-free + open top + opac_tmin 2500, pending T6 passing 1.5e6. prod11 (lidded) is
the running control (t 4.1e6, L_out/L 0.62-0.70). Corona agent's ghost/RT-top changes and
opac_tmin are all UNCOMMITTED in build_rg_prod.
