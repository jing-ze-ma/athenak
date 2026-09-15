---
name: red-giant-hot-corona-is-a-parker-wind
description: "V8f (6e5 K corona) POST-MORTEM 09-10 18:20: dropped. The 6e5 K corona is above the Parker critical T at the join (~1.9e5 K; c_s/v_esc 1.14, sonic radius inside the star, H = 1.9 R) so it is a thermal wind, not a static atmosphere; the thin stellar atmosphere (14 scale heights between the photosphere and the join) drains 11x, the corona becomes 1.55x over-pressured at a 1-cell 100x density jump, shocks down, and dfloor cells with NO velocity ceiling (eos_vceil=0) spike v to 7e6 -> hydro dt 0.95 s at 3.2115e5, intermittent after"
metadata:
  type: project
---
Dump cadence (2e4) missed the collapse itself (last dump 3.2003e5, collapse 3.2115e5, hydro dt
30.6 -> 0.95 in ~4 cycles, then alternating 0.8-2 s / 30 s). Precursor pinned at r/R 1.16-1.20
(i=408-429), the thin atmosphere just under the corona contact, spread over panels 0/3/5 (NOT
the vertex, NOT the corona: corona dt 620-670 s throughout, corona rho only -8%).
Trend 2.6e5 -> 3.2e5: rho at 1.13 R 1.8e-16 -> 1.6e-17; join cell rho 1.7e-22 -> 7e-23;
cell above the join 3468 -> 8463 K; p_corona/p_below 0.96 -> 1.55; T_max in the atmosphere
2540 -> 3.0e4 K, v_r up to +7.4e6; dfloor hits 0 -> 943/cycle (ghosts included, so invisible
in dumps). L_out/L 0.71 -> 0.67 still falling.

Parker numbers at r=3.89e12, GM=1.99e26, mu 0.617: T_crit ~1.9e5 K (the agent wrote 9.6e4;
factor-2 convention, verdict unchanged). 6e5 K: c_s/v_esc 1.14, r_sonic 0.19 R, H 1.9 R.
2e4 K (V9f): 0.21, 5.8 R, H 0.063 R = 32 cells -> can sit statically.

**Decision guidance**: no dfloor/density/join-radius choice saves a 6e5 K corona. If a hot
corona is wanted: T <~ 5e4 K AND join at 1.08-1.10 R where rho >~ 1e-16, and enable a velocity
ceiling (eos_vceil) regardless. Otherwise use V9f (2e4 K) or V10 (free fall) -- but both die
at 5.71e5 of a different, interior cause: [[red-giant-common-dt-collapse-571e5]].
Scripts: _analysis_0910/v8f.py, task_V8f_prof.py, task_V8f_top.py (+ _out.txt).
Related: [[red-giant-corona-join-cell-cooling]] (V8, old RT).
