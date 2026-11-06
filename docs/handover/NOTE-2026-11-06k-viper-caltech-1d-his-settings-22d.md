# NOTE 2026-11-06 viper -> Caltech (USER): 1-D at Linhao's settings, to 22.6 d: no film, no ~18 d inflow, KE ~10-1000x below his Fig. 2

**Result: `his1d_30d` (his units, tables, GM 260.795, cfl_rad 0.025, nzeta 5, seeded, 384 uniform) stays quiet through his inflow time. KE never exceeds 3.3e40 erg, at 3.9 d.** It is still running to 30 d, with no NaN or fatal so far. I will update this note's numbers if 22.6-30 d changes anything.

His Fig. 2 values below are read by eye from the 3-D plot, so treat them as approximate. All this is 1-D against his 3-D, so it does not say his 3-D is wrong.

## KE in the pi/3 sr wedge [erg]
| t [d] | 1-D median | 1-D max | his 3-D, by eye |
|---|---|---|---|
| 0-2 | 9.4e39 | 1.6e40 | 1e40 -> ~7e40 |
| 2-4 | 9.6e39 | 3.3e40 | ~5e40, spiky from ~3 d |
| 4-6 | 1.4e40 | 3.0e40 | 1-5e41, dips ~3e40 |
| 6-8 | 1.0e40 | 2.0e40 | 1-5e41 |
| 8-10 | 3.9e39 | 1.9e40 | up to ~9e41 |
| 10-12 | 3.5e39 | 5.1e39 | ~5e41-1e42 |
| 12-14 | 3.3e39 | 1.0e40 | ~1-2e42 |
| 14-16 | 9.2e39 | 2.0e40 | ~2e42 |
| 16-18 | 4.5e39 | 1.2e40 | rises to ~5e43 |
| 18-20 | 3.3e39 | 8.0e39 | ~4e43, falling |
| 20-22 | 3.7e39 | 9.3e39 | ~1-2e43 |

KE dips at 3.5, 4.4, 5.2, 6.2, 7.2, 8.6, 9.7 d, about 1 d apart and similar by eye to his spikes at 3-10 d. So the ~1 d oscillation appears without 3-D. It fades in 1-D after ~8 d while his grows.

## Wall and energy against time
First-face to third-face N^2 in 1e-10 s^-2; faces and cells as in 11-06b (dr 0.15625). dE is full sphere; dEgrav, dU in pi/3 sr. Lc0 = cell-centred 4 pi r^2 F_lab, cell 0.
| t [d] | T/IC cells 0-3 | rho/IC cell 1 | N^2 faces 0, 1, 2 | rms N^2 cells 0-23 | dEgrav, dU [1e44] | dE [1e45] | Lc0/L_* | R(1e-10 g/cc) |
|---|---|---|---|---|---|---|---|---|
| 1 | 1.022, 0.991, 1.011, 1.009 | 1.02 | -10, +7, +5 | 11 | +0.11, -0.19 | -0.10 | -3 | 54.6 |
| 3 | 1.026, 0.989, 1.012, 1.006 | 1.02 | -26, +12, -5 | 12 | -0.21, -0.30 | -0.61 | -4 | 52.3 |
| 5 | 1.027, 0.989, 1.011, 1.002 | 1.02 | -29, +12, -12 | 14 | -0.21, -0.63 | -1.00 | 0 | 52.0 |
| 10 | 1.031, 0.987, 1.010, 0.993 | 1.03 | -45, +10, -30 | 19 | -0.74, -1.18 | -2.29 | -4 | 52.4 |
| 16 | 1.039, 0.984, 1.007, 0.977 | 1.04 | -70, +6, -60 | 30 | -1.37, -2.06 | -4.11 | -6 | 52.3 |
| 18 | 1.042, 0.982, 1.005, 0.969 | 1.05 | -81, +2, -74 | 36 | -1.74, -2.41 | -4.98 | -9 | 52.4 |
| 22 | 1.049, 0.978, 1.004, 0.954 | 1.06 | -109, +9, -104 | 49 | -2.40, -3.23 | -6.75 | -9 | 52.0 |

## Reading
1. **No wall film to 22 d.** T/IC of cell 1 is 0.978, against 0.26-0.36 at cfl_rad 1 by 5-10 d. There is no violent phase. Jacobi: mean 1.2 iterations, max 18, none capped, over 828000 logged solves, in line with the paper's "less than 10".
2. **But the wall is slowly drifting.** First- and third-face N^2 steepen about steadily to ~-105e-10 by 22 d, alternating with a near-zero second face. T/IC of cell 3 falls 5 %. It is mild but is not stationary, and it could matter on longer runs.
3. **Your 3-cell -15e-10 layer:** at 3 d this is -26, +12, -5 on faces of half your dr. Same scale at the first face, but not the same pattern as your -16, -20, -11.
4. **Energy tracks his until his jump.** At ~16 d he has dEgrav ~ -2.2e44 and dU ~ -2.5e44; this has -1.4e44 and -2.1e44. After his 18 d drop he is near -5.5e44 in both, while this reaches only -2.4e44 and -3.2e44 at 22 d. The loss is ~10 L_* throughout.
5. **Mass flow:** r 25-35 at most ~2 /yr, mostly < 1 /yr, against his ~+-4 /yr at 3-9 d, and no 18 d event. Closer to the wall it rises slowly to ~3.5 /yr at 22 d (r 20.5-25).

## What it implies
- His 18 d inflow does not occur in 1-D at his settings, at least by 22.6 d. Whatever triggers it is likely missing from 1-D, most simply lateral motions, but this run does not prove that. The Orion 3-D run tests it directly.
- If the Orion run is also quiet, his deck does not reproduce his figure and other causes need chasing (including which GM was published; NOTE-2026-11-06h).
- **Use for you:** `his1d_30d` is a cheap Athena++ 1-D reference at his settings, in case you want to test whether your 1-D matches it.

## Files, in `docs/handover/athenapp-bsg/cflrad_1106/`
`final1d_summary.txt`, `final1d.py`, `his1d_30d_wall_profiles.txt` (first 24 cells at the times above, double precision, cgs; columns in its header), `his1d_fig2style.png` (his axes; drawn when the run was at about 4 d, so its panels stop early).
