---
name: red-giant-open-top-prod12-recipe-fails
description: 2026-09-10 06:50 — V7 (the full prod12 recipe FROM SCRATCH, open top, opac_tmin=3200, wb_rmax, efloor_from_ekin, random chart-free seed) FAILS at 5.8-6.0e5 by a COLD collapse of the thin atmosphere (energy+density floor cells at 31.6 K blown at 1.4e7 cm/s -> hydro CFL), EARLIER than V4 (9.45e5, clamp 2500); with 3200 the open-top atmosphere is 19x thinner and 250 K cooler by 4e5. The open-top from-scratch problem is NOT solved; the lidded path (B12) is
metadata:
  type: project
---
V7_prod12gate (195374, cancelled 06:30 at 6.08e5, dt 0.2-0.9 s). Clean to 4.4e5 (0 collapse,
0 floor marks, 0 vertex cold cells ever; vertex spread <= 1.0 -> the seed works). Then
5.6e5 -> 5.8e5 -> 6.0e5: cold cells 10 -> 141 -> 717, floor-mark cells 0 -> 3 -> 17, min T
31.6 K (e/rho = 1.423e7 exactly) with rho at dfloor 1e-18, in EDGE/interior cells at i_arr
298-316 (at/above the photosphere, not the top cell). Hydro dt 30.65 -> 4-7 s at 5.918e5
(rank 26) -> 0.2 s; the CFL cell is a floor cell at |v| 1.43e7 (cs 3.95e4): velocity, not
sound speed. A 34055 K conduction-limited cell appears only afterwards. rg.log efloor_de
7.5e6 -> 8.9e9 over 3500 cycles; open-ghost guard fallbacks 10 -> 11259 per report.
V7 vs V4 at 4e5: same min top-cell T (1675 vs 1642 K) but median rho at i=319 3.7e-13 vs
7.1e-12 (19x thinner) and median T 1816 vs 2073 K. Agent's inference (not proven):
opac_tmin=3200 thins/cools the open atmosphere; efloor_from_ekin lets the floor eat KE
silently. => "opac_tmin=3200 fails from scratch" (old note) CONFIRMED and characterised.
OPEN-TOP options left: the hot corona (V8_corona running), or a lid. LIDDED path: B12 (fixed
binary + clamp) ran 3.205e7 -> 3.45e7 with dt flat; census pending.
Files: _analysis_0910/task_V7_cfl.py and the V7 tables in the agent transcript (session
2c3cf987, agent a15f7a87a880eb220).
