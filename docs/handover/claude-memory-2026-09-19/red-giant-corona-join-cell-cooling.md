---
name: red-giant-corona-join-cell-cooling
description: 2026-09-10 07:50 — V8 (6e5 K corona) died at 2.96e5 NOT like R7: the atmosphere held (T -3%), the CORONA evaporated (95% of its mass gone by 2.8e5, accreting at -8.8e6 cm/s onto a cold sink); the trigger is the LAST RADIATIVELY ACTIVE CELL below the kappa cutoff cooling to the T floor (31.6 K) within 2e4 s, coincident with the two-stream "Newton step left e<=0 in 192 cells, rescued to a 99.9% drop" message (also seen in V7 before ITS failure). V9 (2e4 K, 1e-18 floor) instead cools the ATMOSPHERE (3366 -> 2215 K by 1.2e5, join eating inward) = R7's path. Suspect: the RT Newton rescue / the radiative balance of the outermost active cell
metadata:
  type: project
---
V8 join profile (i=437 last active, 438 first corona): t=1e5: 437 at 31.6 K, 3.7e-22, v_r
-8.9e5; 438 (corona) 4.9e5 K, 4.5e-25, -1.4e6. The interface walks inward 3 cells by 2.8e5,
infall reaches -8.8e6 (supersonic -> hydro dt). Event counters 2.2e5: dfloor 73423, efloor
4.6e6. L_out/L started 0.97 (corona sealed the drain, vs V7's 0.60) and decayed to 0.65 as the
corona thinned. R7's signature (atmosphere cooling) is what V9 shows instead.
PHYSICS CHECK: a thin cell absorbing the upward stellar flux should equilibrate at
T ~ T_eff/4^(1/4) ~ 2400 K, never 31 K; reaching the floor in 2e4 s means absorption is
missing or the Newton "rescue" (e -> 1e-3 e) IS the collapse. That message appeared in V7
(2x, at 5.93e5/6.02e5, 1 cell) and in V8 (3x, 192 cells = a shell) right before each failure.
Agent (08:00) investigates two_stream_rt.hpp: the Newton rescue, the J/absorption seen by the
outermost active cell when kappa=0 above (rad_kappa_rmax) or at the open top, and a correct
fallback (explicit/linearised step, or clamp to e_eq, never a fixed 99.9% drop).
Bondi note (V9): at 2e4 K a corona above 1e-18 is 170x over the pgen's 1e-10 Msun/yr bound;
that trade-off is inherent to a cool dense corona.
Files: _analysis_0910/v8.py v9.py v8tab.py v9tab.py. Related: [[red-giant-open-top-prod12-recipe-fails]],
[[red-giant-top-cooling-runaway]], [[red-giant-molecular-opacity-knee-runaway]].

## V9 (2e4 K, dfloor 1e-18, OLD RT) at 3.09e5, 08:15: alive, but the join is ERASED
dt flat (30.63), 0 collapses; corona <T> holds (-6.5%) but drains: M/M0 1 -> 0.78 -> 0.53 ->
0.34 at 1e5/2e5/3e5, <v_r> -2.5e5 -> -1.15e6. The star's top ~10 cells HEAT from 3366 K to
10-17 kK by 2-3e5 and mix with the infalling corona; no density jump left at 3e5 (2.1e-17 ->
1.2e-17 across the old join, was 1.1e-15 / 8.1e-17); velocities no longer 1-D (convection
reached the top). Not a usable configuration as is: the corona accretes onto the over-cooled
(RT bug) cell below the join and is gone within ~5e5. The RT-fixed corona rerun V8f is the
real test; V8f's join cell holds 3353 K at 1.5e4 (was 31 K in V8).
V9 at 4.6e5 (old RT): dt still 30.635, but the ENVELOPE INFLATED THROUGH THE CORONA: below
3.84e12 the 2e4 K gas is replaced by 1700 K stellar gas at 1e-14 (170x the corona density);
the remnant above 3.9e12 heated to 33-36 kK and vents through the open top at +1e6..3e6;
outer mass flux -8.7e-8 (5 orders above V8); L_out/L 0.97 -> 0.58. Not a lid, not a corona.
The 2e4 K corona at 3e-17 has 1e-4 of the pressure the cooled, expanding envelope pushes with.
