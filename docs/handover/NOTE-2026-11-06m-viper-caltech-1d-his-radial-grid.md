# NOTE 2026-11-06 viper -> Caltech (USER): his radial grid, unseeded and seeded, to 30 d (reply to TASK-2026-11-06l)

**Result: no inflow at 17-19 d or at any time to 30 d, on his grid with and without the seed. No wall layer at dr 0.3125 at 3 d. But the coarser wall grid raises KE at 2-6 d by ~3-5x, up to 1.5e41, which reaches the bottom of his 1-5e41 range, and energies at 30 d come close to his final values.** All numbers below are 1-D, against his 3-D read by eye.

## Arms
His code, units, tables, GM 260.795, cfl_rad 0.025, nzeta 5; real 1-D SMR worked directly: root 192, level 1 at 35-55 (checked: dr 0.3125 at the wall, 0.15625 at 45, 0.3125 at 70; 256 cells). No FATAL or NaN in any arm, all to 30 d, ~10 min CPU. dt ~5-7 s late, ~10 s at early times in `his1d_g`, so the 1-D dt is not the 3-D one (~4.8 s, set laterally).
- `his1d_g`: seed off. `his1d_gs`: seeded as coded. `his1d_30d`: 384 uniform, seeded, final numbers here (22.6-30 d changed nothing in 11-06k).

## KE in the pi/3 sr wedge [erg]: median (max)
| t [d] | `his1d_g` | `his1d_gs` | `his1d_30d` | his 3-D, by eye |
|---|---|---|---|---|
| 0-2 | 8.4e39 (2.4e40) | 1.2e40 (2.1e40) | 9.4e39 (1.6e40) | 1e40 -> ~7e40 |
| 2-4 | 3.9e40 (1.5e41) | 2.7e40 (1.0e41) | 9.6e39 (3.3e40) | ~5e40, spiky from ~3 d |
| 4-6 | 4.3e40 (1.2e41) | 2.8e40 (9.9e40) | 1.4e40 (3.0e40) | 1-5e41 |
| 6-8 | 1.1e40 (2.4e40) | 8.6e39 (2.3e40) | 1.0e40 (2.0e40) | 1-5e41 |
| 8-10 | 5.2e39 (1.5e40) | 6.4e39 (2.0e40) | 3.9e39 (1.9e40) | up to ~9e41 |
| 14-16 | 6.4e39 (1.7e40) | 1.0e40 (4.5e40) | 9.2e39 (2.0e40) | ~2e42 |
| 16-18 | 5.6e39 (1.3e40) | 7.2e39 (1.8e40) | 4.5e39 (1.2e40) | rises to ~5e43 |
| 18-20 | 6.5e39 (1.2e40) | 7.9e39 (1.6e40) | 3.3e39 (8.0e39) | ~4e43 |
| 28-30 | 8.5e39 (1.2e40) | 6.8e39 (1.1e40) | 4.5e39 (6.2e39) | ~1e43 |
Peak: 1.54e41 at 3.87 d (`g`), 1.01e41 at 3.96 d (`gs`), 3.28e40 at 3.92 d (uniform). Full series at hst cadence: `*_ke.txt`.

## Your `Lg`
KE at 0.1, 0.25, 0.5, 1 d: `g` 4.74e39, 7.88e39, 8.03e39, 6.72e39; `gs` 1.11e40, 2.06e40, 1.31e40, 5.95e39. So your `Lg` 1.36e40 at 0.5 d is 1.7x `his1d_g`, its like-for-like (unseeded), not the 1.17e40 of the seeded uniform arm.

## Wall, on his grid
Faces between cells at r = 20.16, 20.47, 20.78, 21.09 (your cells). N^2 in 1e-10 s^-2.
| t [d] | `g` N^2 faces 0, 1, 2 | `gs` N^2 faces 0, 1, 2 | `g` T/IC cells 0-3 | `g` rho/IC cell 1 |
|---|---|---|---|---|
| 0 | +9, +9, +9 | +22, -3, +8 | 1, 1, 1, 1 | 1 |
| 1 | +4, +10, +8 | -7, +10, +6 | 1.002, 0.998, 0.999, 0.998 | 1.00 |
| 3 | **-5, +11, +6** | -19, +13, +0 | 1.007, 0.995, 0.998, 0.995 | 1.02 |
| 5 | -9, +11, +2 | -24, +14, -6 | 1.008, 0.993, 0.996, 0.989 | 1.02 |
| 10 | -21, +12, -8 | -34, +14, -18 | 1.017, 0.991, 0.997, 0.980 | 1.03 |
| 18 | -39, +14, -25 | -53, +16, -37 | 1.030, 0.989, 0.998, 0.966 | 1.04 |
| 30 | -68, +18, -56 | -77, +28, -69 | 1.050, 0.986, 1.003, 0.941 | 1.06 |

- **Your layer:** at 3 d, your 1-D v* arm has -16, -20, -10.5 at r 20.16, 20.47, 20.78 on cell centres. His unseeded column on your dr, on faces between about those cells, has -5, +11, +6, so only the first face is negative and it is 1/3 yours. Different stencils (face vs centre), so compare your profile from `his1d_g_wall_profiles.txt` with your own script for a clean number. T/IC at 3 d is also far milder here: cell 0 +0.7 % vs your +2.8 %.
- **Slow drift persists on this grid,** roughly half the uniform grid's rate. First face at 30 d: -68 (`g`), -157 (`his1d_30d`). rms N^2 cells 0-23 is 24 vs 68.
- **The seed** moves the wall early: at 0.5-3 d N^2 first face is ~6-14 more negative and T cell 0 is ~2.5 % hotter. It adds ~25-100 % KE at 8-16 d and lowers it ~30 % at 2-6 d.

## Energy [1e44 erg, pi/3 sr]: dEgrav, dU
| t [d] | `his1d_g` | `his1d_gs` | `his1d_30d` | his 3-D, by eye |
|---|---|---|---|---|
| 10 | -1.65, -1.95 | -1.72, -1.71 | -0.74, -1.18 | ~-1.5, ~-1.5 |
| 16 | -2.66, -2.99 | -2.76, -2.86 | -1.37, -2.06 | ~-2.2, ~-2.5 |
| 20 | -3.33, -3.74 | -3.50, -3.73 | -2.06, -2.83 | ~-5.5, ~-5.5 |
| 30 | -5.19, -5.76 | -5.08, -5.24 | -3.72, -4.82 | ~-6.2, ~-5.2 |
On his grid the energy loss follows his to ~16 d (even a little more), misses his sudden step at 18 d, and by 30 d has drifted to about his final values anyway. The total loss is steady at ~9-10 L_*. Whether that is a coincidence, I can't say.

## Against `his1d_30d`, as you asked
- **Wall layer at dr 0.3125:** none like yours, above.
- **Oscillation:** same timing, larger on his grid. Smoothed KE peaks (d): `g` 2.70, 3.21, 3.87, 4.70, 5.64; `gs` 2.70, 3.29, 3.96, 4.75, 5.63, 7.80; uniform 3.22, 3.92, 4.79, 5.70, 6.69, 7.91. Spacing lengthens from ~0.5 to ~1 d. Peak height is ~4.7x higher on his grid (`g` vs uniform), and in all arms it fades after ~6-8 d, while his keeps growing.
- **Inflow:** none. Mdot/Menv at r 25-35 is at most ~1.5 /yr, mostly < 1 /yr, against his ~+-4 /yr at 3-9 d. Near the wall, r 20.5-25, it creeps up to ~3.5 /yr by 30 d, without an event.
- **Jacobi:** none capped; mean 1.05 (`g`), 1.32 (`gs`), 1.57 (uniform).

So the ~18 d inflow is missing in 1-D across uniform and his grid, seeded and not. The Orion 3-D is the test.

## Files, in `docs/handover/athenapp-bsg/cflrad_1106/`
`his1d_g_wall_profiles.txt`, `his1d_gs_wall_profiles.txt` (first 24 cells at every 0.5 d, double precision, cgs), `his1d_g_ke.txt`, `his1d_gs_ke.txt`, `his1d_30d_ke.txt`, `final_g_summary.txt`, `final1d_summary.txt`, `his1d_30d_wall_profiles.txt` (updated to 30 d), `final_g.py`.
