# NOTE 2026-11-06 Orion -> Caltech: 4 x 4 Athena++ column, first results to 4.9 d (reply to TASK-2026-11-06p)

**Result: the column agrees with `his1d_g`. The 1-D line vs 3-D mesh and angles does not explain your v_r and KE differences.** Interim: job 227389 is still running to 25 d, and I'll update these files when it ends. Files are in `docs/handover/athenapp-bsg/repro3d/col4x4/orion/`: `col4x4_ke.txt`, `col4x4_wall_profiles.txt`, `col4x4_dumpstats.txt`, `compare.txt`, `report.py`.

| | 4 x 4 column | `his1d_g` | your `Lr` (from your TASK) |
|---|---|---|---|
| v_r at 22.66 Rsun, 1 d | -126.5 m/s | -121.7 m/s | -42 m/s |
| KE 2-4 d, median (max) | 4.47e40 (1.76e41) | 3.93e40 (1.54e41) | |
| KE peak | 1.76e41 at 3.863 d | 1.54e41 at 3.869 d | |

KE medians in other bands: 0-1 d 6.91e39 vs 6.56e39; 1-2 d 1.72e40 vs 1.61e40; 4-6 d 5.61e40 vs 4.29e40 (this band is incomplete). So the column is ~5-15 % above `his1d_g` to 4 d, with the KE peak at the same time. In the first cell, shell-mean T and rho match to ~0.03 % or better through 2 d. Fr1 differs by 0.1 % at 1 d, 1 % at 2 d and 15 % at 3 d (-1.71e14 vs -1.49e14), and N^2 on the first face at 3 d is -4.2 vs -4.6 x 1e-10.

## Your other questions
- **dt:** 10.60 s at t = 0, matching yours. It stays ~10.6 s to 1.5 d, then falls: 8.4 s at 2 d, 5.7 s at 3 d, 4.5 s at 4.5 d. Per dump in `col4x4_dumpstats.txt`.
- **Jacobi:** mean 1.0-1.26, max 4 since t = 0, none at the cap.
- **Lateral motion: yes, well above round-off.** The t = 0 lateral rms of drho/rho is 1.3e-16 and it is 1.2-2.0e-4 by 0.5-1 d in cells 0-23. rms v_theta is 2-3e3 cm/s at 0.5-1.5 d and v_phi 1-1.4e2 cm/s. v_phi jumps at ~2 d to rms 6e3 cm/s, max 4.7e4. KE_lat [erg]: 1 d 4.1e+37 of 6.6e+39; 2 d 2.9e+37 of 2.4e+40; 3 d 3.4e+37 of 4.3e+40; 3.86 d 4.3e+37 of 1.8e+41; 4 d 4.2e+37 of 1.3e+41. Cause not established; periodic x2 BCs on a sin(theta) wedge are one candidate. It stays a small share of KE.

## Differences from your setup and the request
- 64 ranks with `meshblock/nx1=4` (USER allowed more cores), not 32 blocks: bsg.hst identical over a 400-cycle test, to 6 printed digits.
- Radial faces are explicit, with no level boundaries; `his1d_g` used SMR with prolongation at 35 and 55. The wall (r < 35) is uniform in both.
- N^2 is on faces, as in `his1d_g`, from shell means; its t = 0 first-face value matches to 4 digits (8.784e-10). Fr01 is added.
