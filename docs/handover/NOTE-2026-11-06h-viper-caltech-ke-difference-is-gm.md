# NOTE 2026-11-06 viper -> Caltech: the 13x in KE is GM (reply to TASK-2026-11-06g)

**Answer: GM carries it, and it is physical.** `af1_cflrad` uses GM = 20.00 Msun; his deck uses 260.795 code = 19.78 Msun. That 1.14 % in gravity alone moves his1d's KE at 0.5 d from 1.17e40 to 1.36e41 (11.7x) and reverses the motion, from outward to inward. **Please check which GM your AthenaK 1-D uses: if 20 Msun, that may be why you follow `af1_cflrad`.**

## Ladder (384 uniform, cfl_rad 0.025, seeded; one change from his1d; files in `athenapp-bsg/cflrad_1106/ladder_1106/`)
KE in pi/3 sr; 1 d values: T/IC cells 0, 1; Lc0 = cell-centred 4 pi r^2 F_lab in cell 0, in L_*; dE full sphere.
| arm | KE 0.1 d | 0.25 d | 0.5 d | 1 d | T/IC | Lc0 | dE [1e45] | 0.5 d max v_r [km/s] @ r, r<45 | KE half-radius |
|---|---|---|---|---|---|---|---|---|---|
| his1d baseline | 6.40e39 | 9.70e39 | 1.17e40 | 5.15e39 | 1.022, 0.991 | -3 | -0.10 | +1.08 @ 33.5 | 33.4 |
| units only (all physics same) | 6.41e39 | 9.71e39 | 1.17e40 | 5.15e39 | same | -3 | -0.10 | +1.08 @ 33.5 | 33.4 |
| **GM = 20 Msun** | **8.21e40** | **1.83e41** | **1.36e41** | **5.23e40** | 1.045, 0.985 | -12 | -0.28 | **-3.57 @ 27.9** | 28.0 |
| nzeta 4 | 6.24e39 | 8.97e39 | 1.06e40 | 4.85e39 | 1.022, 0.991 | -3 | -0.10 | +1.03 @ 33.5 | 34.0 |
| floors = old | 6.40e39 | 9.70e39 | 1.17e40 | 5.13e39 | 1.022, 0.991 | -3 | -0.10 | +1.08 @ 33.5 | 33.4 |
| padded tables | 6.40e39 | 9.70e39 | 1.17e40 | 5.15e39 | 1.022, 0.991 | -3 | -0.10 | +1.08 @ 33.5 | 33.4 |
| seed threshold = old | 6.30e39 | 8.41e39 | 8.31e39 | 6.15e39 | 1.022, 0.991 | -4 | -0.12 | +1.27 @ 41.0 | 41.8 |
| all except GM | 6.15e39 | 7.80e39 | 7.62e39 | 5.80e39 | 1.022, 0.991 | -3 | -0.13 | +1.23 @ 41.0 | 46.6 |
| all five, his units | 8.49e40 | 1.92e41 | **1.58e41** | **5.43e40** | 1.045, 0.985 | -12 | -0.30 | **-3.78 @ 28.2** | 28.2 |
| af1_cflrad | 9.33e40 @ 0.11 d | 1.83e41 @ 0.22 d | **1.58e41** (dump) | **5.41e40** | 1.045, 0.985 | -13 | -0.30 | **-3.78 @ 28.2** | 28.2 |
"All five" reproduces `af1_cflrad` at 0.5 d in KE, max v_r and its r, so nothing else in the two decks matters at this level.

## Item by item
1. **Units:** neutral (6.41e39 vs 6.40e39 at 0.1 d). error_limit, the Jacobi test and floors, once converted, do not depend on units here. This also covers item 5: the old profile file gives the same answer as his.
2. **nzeta:** ~9 % at 0.5 d.
3. **Floors:** nothing near the wall or in KE.
4. **Tables:** no effect at 3 digits, matching your 6e-5 / 2e-4 check. So padded vs original is not the difference. I did not repeat your direct table comparison.
5. **Seed threshold:** 0.7x in KE. The old threshold seeds much further out (1e-11 vs 6e-7 g/cc), and the KE moves outward to r ~ 42-47.
6. **Other deck differences:** Both arms use the same binary and his BC. The remainder (output blocks, comments, `x2` digits, explicit `integrator`, `reduced_factor`, `Compton`) is covered by the row above matching.

## Mechanism
There is no self-gravity, so the point mass must balance the IC. Profile mass is 19.82 Msun at 20 Rsun and 20.00 at 80. A smooth fit of the IC gives (dP_tot/dr)/(rho g) + 1 over r = 21-40 Rsun:
- his GM: median **-0.45 %** (pressure slightly exceeds gravity: net outward)
- 20 Msun: median **+0.68 %** (net inward)

So at 20 Msun the whole interior starts falling, which is the coherent infall front in NOTE-2026-11-06 (r 20 -> 34 Rsun by 1 d). That front was largely this GM choice, as well as cfl_rad 1. The fit is rough (seed scatter, order 6), so use the sign and size only.

## What this changes
- **Two settings differed from his deck in all earlier 1-D and audit arms:** cfl_rad (1 vs 0.025) and GM (20 vs 19.78 Msun). The earlier NOTEs should be read with both in mind.
- **3-D runs:** `prod_seeded` (Raven, Orion copy) uses his deck, so GM = 260.795.
- **A caution:** his Fig. 2 KE is ~7e40 near 0.5-1 d, read by eye. At 1 d the 20 Msun arm gives 5.2e40 and his GM gives 5.2e39. So at this time the 20 Msun arm is closer to his 3-D figure than his own GM in 1-D. This could be 3-D, or this `rusty` deck may not be the published run's. Later on neither fits: `af1_cflrad` is 1-5e40 at 2-10 d, against 1-5e41. Only the 3-D run can separate these. Worth asking Linhao which GM the published run had.

`his1d_30d` is still running; it tests the ~18 d inflow at his GM.
