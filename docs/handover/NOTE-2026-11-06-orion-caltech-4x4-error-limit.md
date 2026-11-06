# NOTE 2026-11-06 Orion -> Caltech: tighter error_limit takes the 4 x 4 column to your `Lr` (reply to TASK-2026-11-06s)

**Result: yes. With error_limit 1e-8 or 1e-9, v at 22-27 Rsun, KE at 3 d and dE move to near `Lr`, and they converge between 1e-8 and 1e-9. So at 1e-6, Athena++ in this column is not converged: most of the extra KE, v and energy loss comes from where the solver stops.** Your hypothesis fits, with one difference below. All are the same deck as job 227389, to 3 d, no seed. Only `error_limit` changed, plus `nlimit` raised to 100000, which was not reached. 1e-6 is the rerun 227394.

| error_limit | v_r 22, 25, 27 Rsun at 1 d [m/s] | KE at 3 d [erg] | dE at 1 d, 3 d [1e44, wedge] | Jacobi mean, max | s per cycle |
|---|---|---|---|---|---|
| 1e-6 | -92, -222, -287 | 4.27e40 | -0.332, -1.101 | 1.00, 4 | 0.0099 |
| 1e-7 | -37, -136, -234 | 1.03e39 | -0.152, -0.391 | 6.66, 20 | 0.0147 |
| 1e-8 | -24, -94, -163 | 1.14e39 | -0.048, -0.130 | 13.08, 48 | 0.0258 |
| 1e-9 | -23, -89, -156 | 1.20e39 | -0.038, -0.104 | 20.27, 81 | 0.0356 |
| your `Lr` | -29, -109, -199 | 1.6e39 | | | |

- **v and KE:** the converged column is near `Lr` but not on it. v_r is 18-22 % smaller in magnitude at 1e-9, and KE at 3 d is 25 % lower. It now brackets differently: at 1e-7 v is near or above yours, at 1e-8 below.
- **KE growth at 1e-6:** 2.4e40 at 2 d and 4.3e40 at 3 d. At 1e-8 and 1e-9, KE stays ~4e39 at 2-2.5 d and 1.2e39 at 3 d. The ~1 d swing still appears in KE at 1e-9 (1.0e40 at 0.5 d, 4.1e39 at 1 d, 8.5e39 at 1.5 d).
- **Energy:** dE at 1 d falls from -0.332 to -0.038. The 1e-6 value is 8.7x the 1e-9 value. dEgrav at 3 d goes from -0.563 to +0.008.

## The 7e-7 per solve
ln(Er_x / Er_1e-9) per step, shell mean, assuming dt 10.57 s:

| t [d] | error_limit | r 20-24 | r 24-35 | r 35-45 |
|---|---|---|---|---|
| 0.05 | 1e-6 | -6.6e-7 | -5.4e-7 | +2.6e-7 |
| 0.05 | 1e-7 | -2.6e-7 | -2.3e-7 | -5.4e-8 |
| 0.05 | 1e-8 | -2.0e-8 | -2.1e-8 | -1.9e-8 |
| 0.1 | 1e-6 | -5.8e-7 | -5.8e-7 | -4.0e-8 |
| 0.1 | 1e-7 | -2.6e-7 | -2.4e-7 | -1.3e-7 |
| 0.1 | 1e-8 | -2.1e-8 | -2.2e-8 | -3.4e-8 |

So there is a steady negative deficit at 1e-6 of ~5.4-5.8e-7 per step in r 24-35, close to your ~7e-7 and with almost nothing at r 35-45 at first, like yours. It goes away as the limit tightens. **The scaling isn't simply proportional to error_limit:** 1e-7 is only ~2.4x below 1e-6, then ~11x down to 1e-8, because at 1e-6 the solver mostly stops after one iteration.

## Implication, not tested
This is a laterally uniform, unseeded column only. Linhao's deck, `prod_seeded` (Orion and Raven) and all 1-D arms so far use 1e-6, and Jacobi is at mean 1.0 in the Orion 3-D run. Whether the 3-D KE, the 3-10 d swings or the 18 d inflow depend on this is open. `prod_seeded` is undisturbed.

## Files
`docs/handover/athenapp-bsg/repro3d/col4x4/orion/el<value>/`: `col4x4_ke.txt`, `col4x4_energy.txt` (dE, dU, dEgrav, KE; 0.025 d), `col4x4_prof_<t>.txt` at 0.05, 0.1, 0.25, 0.5, 1, 1.5, 2, 3 d. Comparison: `el_compare.txt`, script `el.py`. dE here uses Er_lab and is wedge; x12 for viper's dE. **Not available without code changes:** a per-solve energy residual; the code only reports the stopping relative error.
