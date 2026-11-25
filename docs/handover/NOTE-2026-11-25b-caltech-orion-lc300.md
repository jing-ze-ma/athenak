# NOTE 2026-11-25b caltech -> orion: 300-cycle continuity of rel-f on real states (your 25b item 3)

## t9l62 (11-25 00:42): rst t9l62_t9l62.00003.rst (hard link, c6332), rel-f 80d16d39bd1f47bd3f22b5711c660519 against athena_c4speed_d979d5268, 300 cycles, 02-16; own keys: rad_vet/time2_gtime_slope=true rad_vet/clean_stag_thinf=true rad_vet/clean_stag_lat_sums=true rad_vet/clean_stag_tie_gid=true rad_vet/clean_stag_opac_drag=true rad_vet/vet_gd_allow_clamp=false hydro/vceil=3.0e8 rad_vet/clean_pc_sweep=rb rad_vet/clean_stag_mg_levels=2 time/cfl_number=0.3 rad_vet/clean_conv_max=1e-5; R_ph by rph.py (page agent's definitions, 4-byte dumps); dumps at global cycles that are multiples of 100, so about c 68 / 168 / 268 after the rst (the NOTE's c100 / 200 / 300), and the end
- F300: 2026-11-25 00:27:36 rc=0 (srun) ended: Terminating on cycle limit wall=853 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- C300: 2026-11-25 00:41:47 rc=0 (srun) ended: Terminating on cycle limit wall=850 s CW 0 FATAL 0 inf 0 DT_RETRY 0
| arm | rows | s/cyc (elapsed prints) | dt min/mean/max | its (sum/cycle), solves | CW/FATAL/inf/DT_RETRY | fallback/BE |
|---|---|---|---|---|---|---|
| F300 | 300 | 2.6134 (c6400-6600) | 54.2218/58.9475/66.0362 | 17.9, 6.79 | 0/0/0/0 | 0.000000e+00/0.000000e+00 |
| C300 | 300 | 2.6098 (c6400-6600) | 54.2218/58.9475/66.0362 | 17.9, 6.78 | 0/0/0/0 | 0.000000e+00/0.000000e+00 |

| row | t (rel-f) | L_top/L_in | rel L_top, L_bot, L_mid, L_int, KE, M_tot, Etot, dt (rel-f / running - 1) |
|---|---|---|---|
| 0 | 5.184690e+05 | 2.141651 | 2.22e-16, -1.03e-10, -3.00e-15, 3.23e-13, -4.44e-16, 0.00e+00, -2.22e-16, 0.00e+00 |
| 9 | 5.189596e+05 | 2.138847 | 1.69e-14, -9.61e-10, -7.85e-12, 3.88e-13, -1.45e-13, -4.44e-16, 6.66e-16, -1.63e-14 |
| 49 | 5.211484e+05 | 2.127281 | -1.41e-07, -4.79e-09, -3.66e-04, 4.02e-09, -4.55e-12, -2.00e-15, 2.00e-13, 1.05e-08 |
| 99 | 5.240140e+05 | 2.113626 | -6.62e-06, -1.33e-08, 2.46e-04, 2.13e-08, 4.06e-09, -1.11e-16, -2.55e-11, -1.18e-07 |
| 199 | 5.298664e+05 | 2.090431 | 3.51e-06, 3.68e-09, 1.66e-04, -7.97e-09, -5.97e-09, -1.11e-15, 1.58e-11, -3.46e-08 |
| 299 | 5.360989e+05 | 2.073188 | -3.59e-06, -6.91e-09, 4.15e-04, 1.06e-08, 7.34e-09, -6.66e-16, -2.21e-11, 6.28e-08 |
| 299 | 5.360989e+05 | 2.073188 | -3.59e-06, -6.91e-09, 4.15e-04, 1.06e-08, 7.34e-09, -6.66e-16, -2.21e-11, 6.28e-08 |
- L_top: max |rel-f / running - 1| over 300 rows 2.325e-05
- L_bot: max |rel-f / running - 1| over 300 rows 2.797e-05
- L_mid: max |rel-f / running - 1| over 300 rows 1.128e-03
- L_int: max |rel-f / running - 1| over 300 rows 9.723e-07
- KE_int: max |rel-f / running - 1| over 300 rows 2.724e-08

| dump | cycle | t | R_A (mean tau 2/3) | R_A / first production dump before - 1 | T_A [K] | R_M (median column) | T_M [K] |
|---|---|---|---|---|---|---|---|
| production before | 5943 | 5.076050e+05 | 9.520301e+12 | -2.67e-02 | 16961.93 | 9.185541e+12 | 17071.23 |
| production before | 6332 | 5.184150e+05 | 9.781226e+12 | 0.00e+00 | 16582.88 | 9.493198e+12 | 16736.69 |
| F300 | 6400 | 5.221580e+05 | 9.871817e+12 | 9.26e-03 | 16475.01 | 9.595464e+12 | 16630.06 |
| F300 | 6500 | 5.280440e+05 | 1.001634e+13 | 2.40e-02 | 16312.77 | 9.758380e+12 | 16470.90 |
| F300 | 6600 | 5.340280e+05 | 1.016493e+13 | 3.92e-02 | 16166.17 | 9.918859e+12 | 16329.88 |
| F300 | 6632 | 5.360990e+05 | 1.021645e+13 | 4.45e-02 | 16117.25 | 9.974319e+12 | 16285.82 |
| C300 | 6400 | 5.221580e+05 | 9.871816e+12 | 9.26e-03 | 16474.79 | 9.595456e+12 | 16629.30 |
| C300 | 6500 | 5.280440e+05 | 1.001634e+13 | 2.40e-02 | 16312.67 | 9.758372e+12 | 16472.25 |
| C300 | 6600 | 5.340280e+05 | 1.016493e+13 | 3.92e-02 | 16164.52 | 9.918840e+12 | 16330.11 |
| C300 | 6632 | 5.360990e+05 | 1.021645e+13 | 4.45e-02 | 16119.77 | 9.974328e+12 | 16284.33 |
| production after | 6521 | 5.292430e+05 | 1.004729e+13 | 2.72e-02 | 16280.56 | 9.789622e+12 | 16444.71 |
| production after | 6695 | 5.400260e+05 | 1.031637e+13 | 5.47e-02 | 16034.17 | 1.008700e+13 | 16201.88 |
| production after | 6884 | 5.508170e+05 | 1.059192e+13 | 8.29e-02 | 15833.38 | 1.037768e+13 | 15996.14 |
