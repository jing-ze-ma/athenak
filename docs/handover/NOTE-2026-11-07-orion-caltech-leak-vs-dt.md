# NOTE 2026-11-07 Orion -> Caltech: the leak is close to fixed per unit time, not per solve (TASK-2026-11-06t parts 1, 2)

**Result: in the wall band r 20-24, ln(Er_1e-6 / Er_1e-9) grows at a nearly fixed rate per day, about -9e-3 per d, so per step it scales with dt. It is not fixed per solve.** Early values, 0.05 d, from `leak/leak_vs_dt.txt`. Er is volume-weighted over each band, and steps are t/dt(t). The two runs at a given cfl_rad have the same dt.

| cfl_rad | dt [s] | r 20-24 per step | r 20-24 per day | r 24-35 per day | r 35-45 per day | Jacobi mean at 1e-6, 1e-9 |
|---|---|---|---|---|---|---|
| 0.011 | 4.66 | -4.7e-7 | -8.7e-3 | -6.4e-3 | -8e-4 | 1.0, 12 |
| 0.025 | 10.6 | -1.06e-6 | -8.6e-3 | -5.0e-3 | -4e-4 | 1.0, 20 |
| 0.05 | 21.2 | -2.55e-6 | -1.04e-2 | -4.9e-3 | +3e-4 | 1.4, 43 |
| 0.1 | 42.3 | -6.05e-6 | -1.23e-2 | -3.5e-3 | +1.9e-3 | 2.0, 81 |
| 0.25 | 105.7 | -1.97e-5 | -1.61e-2 | +4e-4 | +6.5e-3 | 5.8, 198 |

- **For your sink:** a rate per unit time, not per solve, is the better fit at dt <= ~20 s. It drifts up by ~2x at r 20-24 by dt 106 s.
- **Not simply a single law at all radii:** r 24-35 stays near -5e-3 per d to dt 21 s, then changes sign by 106 s. At large dt the 1e-9 runs are themselves off, so this ratio stops being only the stopping error. dE at 1 d at 1e-9 is -3.8e42 at 10.6 s, -1.26e43 at 42 s and -3.19e43 at 106 s, which is the velocity error at larger dt.
- **Later times aren't constant:** at 1 d the rates are -4e-3 at r 20-24 and -2.8e-3 at r 24-35 per d, at dt 10.6 s, and r 35-45 grows to -9e-3. The table is the early rate only.
- **Viper's cfl_rad 1 runs:** not compared.

Your 7e-7 per solve at dt 10.55 s matches this if taken at r 24-35 (-5.4 to -6.1e-7). In this band-mean measure, r 20-24 is -1.06e-6, higher than the shell-mean -6.6e-7 in my earlier error_limit note. The difference is the volume weighting and the different average. Use this file's numbers with its definition.

## Short 3-D wedge (1/16 lateral, seeded, 112 ranks)
**At the same dt, the 3-D leak matches the column so far, and lateral KE doesn't differ between 1e-6 and 1e-9.** Only 0.025 d is compared, because 1e-9 is slow (it is at ~0.033 d; 43 iterations mean, 95 max).

| t = 0.025 d | r 20-24 per step | r 24-35 per step | r 35-45 per step |
|---|---|---|---|
| wedge, dt 5.52 s | -4.45e-7 | -3.47e-7 | -7e-9 |
| column, dt 4.66 s, at 0.05 d | -4.70e-7 | -3.46e-7 | -4.3e-8 |

Wedge KE at 0.025 d: total 1.250e39 (1e-6) vs 1.234e39 (1e-9); KE_lat 9.00e38 vs 8.86e38. 1e-6 is complete to 0.25 d: KE 4.07e39, KE_lat 1.07e39, x16 = 6.5e40 total. The seed dominates KE this early, so this doesn't yet test whether KE drifts apart later.
Files: `leak/*_ke.txt` (KE, KE_r, KE_lat at hst cadence), `leak_vs_dt.txt` (all bands, dE, Jacobi, cost), `leak.py`.

## Still coming
`prod_seeded` 1 d data in ~25 min. Velocity-fix arms: not reported yet.
