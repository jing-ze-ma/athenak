---
name: red-giant-clamp-off-open-top-ok
description: "OPAC CLAMP UNNECESSARY IN THE OPEN-TOP CASE TOO (I3, 09-10 23:30): I2 vs I3 (opac_tmin=-1) from V9f rst 5e5 to 6.2e5 are indistinguishable -- same dt history, same cap burst 5.7-5.9e5 (I3 peak x_i 200x SMALLER), top of star 50-90 K WARMER without the clamp, 1/3 as many sub-1000 K cells, no 200-500 K parking; open items: 31.6 K floor cells appear in BOTH at 6.0-6.2e5 (12-33 isolated cells i 337-376), and L_out/L ~2.7-3 after the burst (hot evacuated cells radiating?)"
metadata:
  type: project
---
Job 195879, dir I3_implicit_cap_noclamp, binary = athena_pin2. Only diff to I2: opac_tmin -1.
Numbers at 6.0e5 (shell medians): i=395 last active cell 1976 K (I2) vs 2063 (I3); i=347
3623 vs 3532; corona 19476 both. Cold census i=297-395: N(T<1000) 12 vs 3 at 6.0e5, 33 vs 12
at 6.2e5; min T 31.6 K (I2) vs 1513 (I3) at 6.2e5. hst identical to 6 digits.
face-budget L_out/L: 0.61-0.71 before the burst, 5e2-2.6e6 DURING the cap burst (the budget
line is meaningless then), 2.7-3.0 after (6.05e5) -> check whether it returns to ~0.7-0.9 by
7e5 when the hot population is gone (I2 census says the hot cells dissipate by 7e5).
Conclusion: drop opac_tmin from the production config (B13 lidded + I3 open agree).
Script: _analysis_0910/task_I3_cmp.py. See [[red-giant-molecular-opacity-knee-runaway]],
[[red-giant-implicit-radial-diffusion]].

## FOLLOW-UP 09-11 00:30 (I2 to 9e5, I3 to 8.26e5): L_out/L DOES NOT return to 0.7
Both runs settle at L_out/L 1.6-4.1 from 6.2e5 on (I2 5.08 at 8.2e5 decaying to 3.70 at
9.0e5; I3 on the same curve within 10%); L_cut/L grows 0.05 -> 2.1. Carried by a SPARSE hot
population in shells i=380-395: T up to 8e5 K at rho 3e-14..1e-12 at the 5.8e5 peak, ~1e4 K
(vs 2000-2800 median) at 8-9e5. Clamp-independent. The agent's sigma<T^4> check is NOT
independent evidence (optically thin hot cells do not radiate sigma T^4) but the face budget
is the RT's own flux. My reading: v_r 1-4e6 in the thin atmosphere -> shock heating 1e4-1e5 K;
energy source = convective KE (grid-locked cube-edge/vertex downdrafts?). RT semi-implicit
cannot heat above the radiation temperature, and the new operators are conservative, so
these cells are NOT a conduction/RT artefact. OPEN: is the 3 L emergent flux physical
(shock-heated atmosphere) or the vertex-chimney artefact -> check the hot cells' seam class
and the KE flux into i>380. Floor (31.6 K) cells: recycle, 0 -> 32/50 at 6.2-7.0e5 -> 0 from
7.6e5, always at i~384 next to a 2e4-1.3e5 K neighbour with v_r exactly 0 (cold collapse next
to a shock-heated cell). dt 30.634 throughout; I3 had 4 transient hydro collapses 6.24-6.38e5.
Script: _analysis_0910/task_I3_follow.py.
