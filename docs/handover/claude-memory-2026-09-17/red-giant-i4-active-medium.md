---
name: red-giant-i4-active-medium
description: "I4 (09-11 02:30): FROM SCRATCH with implicit-x1 + cap, clamp OFF, NO opacity cutoff, radiatively active 400 K / 1e-18 free-fall medium, WB on inside (wb_rmax 3.3e12), dfloor 1e-19: PASSED V10's 5.714e5 death (dt 30.634 flat to 6.0e5, one hydro dip to 6 s at 5.896e5, tot-E/E0 = 1.0000000 to 7 digits); medium free-falls at 0.73 v_ff, warms 400 -> 1900 K, thins onto dfloor by 3e5 then refills to 7e-17; star top never floored (min 1382 K), no cold runaway with the clamp off; deep rms v_r 12% QUIETER than V10 (WB on costs nothing); cap burst from 5.72e5 (500-650 cells over cap per 1e4 s at r 3.40-3.48e12, max x_i 1e4), ongoing at 6.0e5; L_out 0.58-0.68 before the burst, 289 during (meaningless)"
metadata:
  type: project
---
Job 195886, dir I4_medium400, binary athena_pin2 (md5 5d9be6a1...). Join at 3.754e12 (V10's
3.720e12). Diffs vs I3: dfloor 1e-19, pfloor 1e-21, rad_kappa_rmax/above REMOVED, bg_rho 1e-18,
bg_temp 400, bg_hydrostatic false, rt_report_r 4e12. Diffs vs V10: + WB on with wb_rmax,
+ rad_implicit_x1 + rad_cap_ang 0.5, clamp off, dfloor/pfloor 10x lower, bg_rho 1e-18.
Medium table (I4 / V10 at 1e5, 3e5, 5e5): rho 1.0e-18/1.0e-17 -> 1e-19/1e-18 (=dfloor!) ->
7.3e-17/9.7e-17; T 605/606 -> 1126/1170 -> 1911/1630 K; v_r -1.2e6 -> -1.8e6 -> +4e5.
Star top band 3.35-3.60e12: T_min 1956 (3e5), 1593 (5e5), 1382 (5.8e5); n(T<1000)=0 always
(V10: 4 floored cells at 3.40e12 at death). L_out/L: 0.94 (9e4), 0.65 (3e5), 0.58 (4.9e5),
0.68 (5.5e5). Deep rms|v_r| r<2e12: 3.18e3 (3e5), 4.75e3 (5e5) vs V10 3.63e3/5.49e3.
OPEN (agent following to the end): does L_out come back to 0.7-1.5 L after the burst without
the cutoff (I2's 2-4 L was the cutoff artefact, [[red-giant-open-top-kappa-cutoff-artefact]])?
Does the cap burst end (I2: 5.98e5)? Medium loading by ejected stellar gas?
Analysis module: _analysis_0910/i4.py (dt_table/cap_table/shell/budget/deep/hst).
Caveat: the cap masks the photospheric stiffness (x_i up to 1e4); dt flat != killer gone.

## POST-MORTEM 09-11 03:30: I4 DIED at 6.313e5 (cancelled), hydro dt 0.03-0.2 s, 18 collapses
Culprit region gid 34 r 3.53-3.57e12 (1.107 R): NOT a floor cell (rho 2.8e-11), NOT cold
(4306 K), NOT hot: a fast OUTFLOWING shell front (+3.9e6 -> +2.6e5 over 6 cells) on a block
edge with 3.9e6 cm/s shear against an infalling angular neighbour. (mesh.cpp prints the cell
only for conduction-limited collapses; hydro ones print nothing -> add the hydro cell to the
report.) ALL open-top runs launch a shell from the photosphere at 5.7e5: M(r>3.6e12) x14 in
I4 (5.6->6.2e5), x140 in I3 by 6.4e5 -- I3 rides it out (4 self-healing collapses), I4 does
not. Discriminator per the agent: I3's pressure-supported 2e4 K corona decelerates the front
(top v +1.5e7) vs I4's unsupported free-fall medium (+2.4e7). Floor and medium opacity ruled
out by the band statistics (I3 has more floor cells and a hotter band).
L_out: I4 never left the burst (289 at 5.82e5, 572 at 6.12e5), so the I2 post-burst 2-4 L
question (cutoff artefact) is NOT tested by I4; the dump photospheric L fell back to 2.3 L by
6.2e5 with Tph50 fixed at 3100 K; the outer-face excess = hot ejected gas radiating in the
medium. Table opacity at the corona/medium conditions = 0.029 (the table's low-rho edge),
so a hydrostatic 2e4 K corona WITH opacity cools in ~2 s to ~700 K and loses its support ->
the agent's I5 (bg_hydrostatic + no cutoff) would likely just reproduce I4. The consistent
design is the user's DENSITY GATE on the opacity (inert below ~1e-14, smoothed), i.e. I3's
corona support without the radius cutoff. Open physics: what launches the shell at 5.7e5
(RCB pile-up ignition overshoot reaching the surface? the downdraft census grows from
4.6e5) -- lidded runs hold it.
