---
name: red-giant-restart-radial-ke-injection
description: FOUND+FIXED 2026-09-09 — the well-balanced cache (wbq0) is all-zero for up to wb_cache_every-1 cycles after a RESTART (gate on ncycle%N, ncycle restored from rst) -> x730 radial-KE kick; fix in hydro/mhd_fluxes (wb_cache_built), COMMITTED+PUSHED 5c0b98e4; prod11 chain protected by hydro/wb_cache_every=1 override
metadata:
  type: project
---

Symptom: T6_nodust restarted at t=5e5 -> hst KE1 jumped 2.06e41 -> 1.5e44 in ~10 cycles, a
spherical radial shell deep in the radiative interior (r/R 0.38: +26 cm/s -> -1.35 km/s); the
star rang, went transonic under the photosphere, NaN'd at 8.6e5 with dt still 30.65 s.

Cause (COMMITTED code, src/hydro/hydro_fluxes.cpp:136 and src/mhd/mhd_fluxes.cpp:160):
`if (wb_cache_every <= 0 || (stage == 1 && ncycle % wb_cache_every == 0)) BuildWBCache()`
-- ncycle is restored from the restart file and wbq0 is freshly zero-allocated, so a restart
at ncycle % 10 != 0 reconstructs against a ZERO background for up to 9 cycles = the
well-balanced scheme silently off in a star balanced to 1e-6. prod_topre's clean restart was
a 1-in-10 fluke. Tests (all from T6's rst): RS1 reproduced (x614 in the first hst row),
RS3 wb_cache_every=1 clean (-0.01%), RS_fix (new flag `wb_cache_built`, forces the build on
the first call after start or restart) clean, identical to RS3 to 4 digits. RS_prod11 with
prod11's own binary: 1-mom x27 in the first row -> prod11 HAS it.

Actions: prod11 chain re-submitted as 194624 -> 194625 with `hydro/wb_cache_every=1` on the
restart line (same binary, physics unchanged, cost unmeasured). T7_nodust_fix (job 194627,
binary md5 ec68df12 = build_rg_prod with the fix) continues the grains-off lid-free star
from T6's t=5e5 rst to 3e6 — the gate for prod12. The open-ghost guard was EXONERATED
(fires only in unread angular-ghost corners) but its fallback to the t=0 column is still a
hazard (p_ghost/p_top 8-130x) — make it a clamped continuation before production.

**How to apply:** the fix is UNCOMMITTED in hydro.hpp/hydro_fluxes.cpp/mhd.hpp/mhd_fluxes.cpp
-- commit it (it affects every WB restart, dhj included). Check hst column 8 across every
restart of any WB run made before it. Related: [[red-giant-dust-opacity-kills-lidfree]],
[[solar-convection-restart-gravity-bug]], [[dhj-restart-loses-gravity-potential]].
