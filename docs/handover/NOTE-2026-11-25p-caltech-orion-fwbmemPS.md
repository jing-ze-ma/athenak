# NOTE 2026-11-25p caltech -> orion: star arms of the PROTOTYPE caltech-fwbmem-4a7885d73 (form P + S_pos, 25r) against next 79c1cabef (your 25p; result-changing, not a release)

## bsg (11-25 03:17): rst bsg_bsg3d.00008.rst (the 25b pair's hard link, c29127), prototype 6aa3ab0f against next 9bebf6a3, 300 cycles, 02-03; own keys: none; R_ph by rph.py; dumps at global cycles that are multiples of 100 and the end
- FW300: 2026-11-25 03:17:40 rc=0 (srun) ended: Terminating on cycle limit wall=237 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- NX300: 2026-11-25 02:33:17 rc=0 (srun) ended: Terminating on cycle limit wall=232 s CW 0 FATAL 0 inf 0 DT_RETRY 0
| arm | rows | s/cyc (elapsed prints) | dt min/mean/max | its (sum/cycle), solves | CW/FATAL/inf/DT_RETRY | fallback/BE |
|---|---|---|---|---|---|---|
| FW300 | 300 | 0.6875 (c29200-29400) | 160.3641/164.6545/169.4520 | 36.1, 6.48 | 0/0/0/0 | 0.000000e+00/0.000000e+00 |
| NX300 | 300 | 0.6725 (c29200-29400) | 160.3591/164.6542/169.4520 | 36.2, 6.37 | 0/0/0/0 | 0.000000e+00/0.000000e+00 |

| row | t (prototype) | L_top/L_in | rel L_top, L_bot, L_mid, L_int, KE, M_tot, Etot, dt (prototype / next - 1) |
|---|---|---|---|
| 0 | 4.838736e+06 | 1.019585 | -3.31e-08, -3.98e-09, -2.16e-12, 1.61e-03, 6.34e-10, 1.78e-15, 1.40e-14, 0.00e+00 |
| 9 | 4.840260e+06 | 1.019423 | -2.54e-06, -1.49e-08, -2.15e-09, 1.26e-03, 5.52e-09, 3.33e-15, 1.61e-12, 3.09e-09 |
| 49 | 4.846913e+06 | 1.018745 | -1.25e-05, -7.63e-10, -1.32e-07, 1.89e-03, -1.67e-10, -8.88e-16, 1.80e-11, 3.78e-07 |
| 99 | 4.854984e+06 | 1.017629 | 2.21e-05, 2.10e-07, 6.60e-07, 2.01e-03, -3.36e-08, -1.22e-15, 7.27e-11, 2.73e-06 |
| 199 | 4.871682e+06 | 1.014420 | 2.67e-05, -6.73e-06, 1.50e-06, -1.94e-03, 3.74e-08, -1.89e-15, 2.02e-10, 1.68e-05 |
| 299 | 4.887963e+06 | 1.011192 | 5.18e-06, -8.14e-09, -2.36e-06, 6.80e-03, 1.80e-07, -1.89e-15, 2.52e-10, -1.30e-07 |
| 299 | 4.887963e+06 | 1.011192 | 5.18e-06, -8.14e-09, -2.36e-06, 6.80e-03, 1.80e-07, -1.89e-15, 2.52e-10, -1.30e-07 |
- L_top: max |prototype / next - 1| over 300 rows 8.089e-05
- L_bot: max |prototype / next - 1| over 300 rows 6.732e-06
- L_mid: max |prototype / next - 1| over 300 rows 4.540e-06
- L_int: max |prototype / next - 1| over 300 rows 6.832e-03
- KE_int: max |prototype / next - 1| over 300 rows 1.798e-07

| dump | cycle | t | R_A (mean tau 2/3) | R_A / first production dump before - 1 | T_A [K] | R_M (median column) | T_M [K] |
|---|---|---|---|---|---|---|---|
| production before | 28863 | 4.795360e+06 | 3.649381e+12 | 3.69e-03 | 15710.70 | 3.633354e+12 | 15752.34 |
| production before | 29127 | 4.838570e+06 | 3.635978e+12 | 0.00e+00 | 15710.83 | 3.621075e+12 | 15747.99 |
| FW300 | 29200 | 4.850640e+06 | 3.638528e+12 | 7.01e-04 | 15641.10 | 3.619494e+12 | 15665.84 |
| FW300 | 29300 | 4.867140e+06 | 3.649411e+12 | 3.69e-03 | 15617.39 | 3.624196e+12 | 15611.61 |
| FW300 | 29400 | 4.883610e+06 | 3.660551e+12 | 6.76e-03 | 15582.33 | 3.633482e+12 | 15572.37 |
| FW300 | 29427 | 4.887960e+06 | 3.662740e+12 | 7.36e-03 | 15574.64 | 3.635178e+12 | 15569.43 |
| NX300 | 29200 | 4.850640e+06 | 3.638542e+12 | 7.05e-04 | 15703.03 | 3.619526e+12 | 15739.21 |
| NX300 | 29300 | 4.867140e+06 | 3.649476e+12 | 3.71e-03 | 15666.07 | 3.624299e+12 | 15679.19 |
| NX300 | 29400 | 4.883610e+06 | 3.660672e+12 | 6.79e-03 | 15569.61 | 3.633680e+12 | 15632.31 |
| NX300 | 29427 | 4.887960e+06 | 3.662876e+12 | 7.40e-03 | 15552.48 | 3.635371e+12 | 15627.31 |
| production after | 29388 | 4.881660e+06 | 3.659595e+12 | 6.50e-03 | 15577.71 | 3.633038e+12 | 15636.83 |
| production after | 29657 | 4.924800e+06 | 3.671902e+12 | 9.88e-03 | 15547.73 | 3.650319e+12 | 15576.49 |

- FW300: CLEAN FCONE lines 0 (faces sum 0, max 0 per line, cycles with a cut 0); FCONE_LZ lines 0; CLEAN FCS/bound lines: ['CLEAN FALLBACK', 'CLEAN FLOORMOM', 'CLEAN FWB', 'CLEAN LZB', 'CLEAN OLDVEC', 'CLEAN WBFCT', 'CLEAN XT']; clamped lateral reads (all ranks, end) 0.000000e+00; CLEAN WARNING 0
- FW300: top mass loss Mdot_top mean 1.0272e+08 (first 6.6432e+08, last 2.0386e+07, min 5.3531e+06, max 6.6432e+08); cumulative M_tot change -2.331e-15; L_mid / L_midE - 1: mean 1.699e-09 max abs 7.118e-09; L_intE/L_int - 1 mean 2.683e-03
- NX300: CLEAN FCONE lines 0 (faces sum 0, max 0 per line, cycles with a cut 0); FCONE_LZ lines 0; CLEAN FCS/bound lines: ['CLEAN FALLBACK', 'CLEAN FLOORMOM']; clamped lateral reads (all ranks, end) 0.000000e+00; CLEAN WARNING 0
- NX300: top mass loss Mdot_top mean 1.0354e+08 (first 6.6107e+08, last 1.6998e+07, min 5.5008e+06, max 6.6107e+08); cumulative M_tot change 1.332e-15; L_mid / L_midE - 1: mean 1.873e-09 max abs 6.180e-09; L_intE/L_int - 1 mean 3.796e-03
- L_midE prototype - next: max |diff| / max |next| 4.521e-06, at the last row -2.346e-06
- Mdot_top prototype - next: max |diff| / max |next| 1.142e-02, at the last row 5.126e-03
- E_rad prototype - next: max |diff| / max |next| 1.313e-09, at the last row 1.301e-09
- e_gas prototype - next: max |diff| / max |next| 9.387e-11, at the last row 2.446e-11
- Lz_tot prototype - next: max |diff| / max |next| 3.505e-02, at the last row 3.505e-02
- FW300: L_top / L_in first 1.019585, last 1.011192, mean of the last 100 rows 1.012748
- NX300: L_top / L_in first 1.019585, last 1.011187, mean of the last 100 rows 1.012749
- FW300 CLEAN FCONE lines 0; first: - | last: -
- FW300 CLEAN FWB lines 600; first: CLEAN FWB cycle 29127 t2s 2 faces 19661297 radiation momentum sum |dp| 1.922572e+28 x1 faces sum dp 1.450066e+27 | last: CLEAN FWB cycle 29426 t2s 3 faces 19397924 radiation momentum sum |dp| 1.678538e+28 x1 faces sum dp -1.213240e+26
- NX300 CLEAN FCONE lines 0; first: - | last: -
- NX300 CLEAN FWB lines 0; first: - | last: -
- mass-weighted shell-mean v_r [cm/s] of the last dump (c29427 FW300 / c29427 NX300): r: FW300 | NX300
  1.4023e+12: 3.163e+00 | 3.163e+00; 1.8154e+12: 2.942e+03 | 2.942e+03; 2.2502e+12: 7.126e+03 | 7.126e+03; 2.5491e+12: 8.976e+03 | 8.976e+03; 2.7665e+12: 4.006e+03 | 4.006e+03; 2.9839e+12: 7.228e+03 | 7.228e+03; 3.1904e+12: 3.962e+04 | 3.962e+04; 3.4078e+12: 1.391e+05 | 1.392e+05; 3.6144e+12: 3.905e+05 | 3.921e+05; 3.8372e+12: 1.318e+05 | 1.312e+05; 4.2721e+12: 3.544e+03 | 3.467e+03; 4.6851e+12: 2.872e+03 | 2.814e+03; 5.1199e+12: 2.368e+03 | 2.322e+03; 5.5547e+12: 2.002e+03 | 1.979e+03
- L_int / L_intE and L_mid / L_midE, Mdot_top at rows (0-based) 0, 1, 5, 20, 100, 299:
  - FW300: r0 Li 0.998067 Lm 1.000000 Mdot 6.643e+08; r1 Li 0.998074 Lm 1.000000 Mdot 6.378e+08; r5 Li 0.998111 Lm 1.000000 Mdot 5.473e+08; r20 Li 0.998237 Lm 1.000000 Mdot 3.259e+08; r100 Li 0.998280 Lm 1.000000 Mdot 7.879e+07; r299 Li 0.996238 Lm 1.000000 Mdot 2.039e+07
  - NX300: r0 Li 0.996467 Lm 1.000000 Mdot 6.611e+08; r1 Li 0.996520 Lm 1.000000 Mdot 6.366e+08; r5 Li 0.996701 Lm 1.000000 Mdot 5.452e+08; r20 Li 0.996987 Lm 1.000000 Mdot 3.266e+08; r100 Li 0.996336 Lm 1.000000 Mdot 7.636e+07; r299 Li 0.989542 Lm 1.000000 Mdot 1.700e+07
- top cell (shell mean) FW300 c29200: rho 6.1525e-14, v_r 2.0102e+03, tau (cell centre) 1.9977e-04
- top cell (shell mean) FW300 c29300: rho 6.1525e-14, v_r 2.0937e+03, tau (cell centre) 1.9986e-04
- top cell (shell mean) FW300 c29400: rho 6.1525e-14, v_r 1.9891e+03, tau (cell centre) 1.9981e-04
- top cell (shell mean) FW300 c29427: rho 6.1525e-14, v_r 2.0017e+03, tau (cell centre) 1.9981e-04
- top cell (shell mean) NX300 c29200: rho 6.1525e-14, v_r 1.9986e+03, tau (cell centre) 1.9977e-04
- top cell (shell mean) NX300 c29300: rho 6.1525e-14, v_r 2.0800e+03, tau (cell centre) 1.9985e-04
- top cell (shell mean) NX300 c29400: rho 6.1525e-14, v_r 1.9669e+03, tau (cell centre) 1.9981e-04
- top cell (shell mean) NX300 c29427: rho 6.1525e-14, v_r 1.9786e+03, tau (cell centre) 1.9981e-04

## m35mms (11-25 03:19): rst m35mms_m35mms.00004.rst (the 25b pair's hard link, c22144), prototype 6aa3ab0f against next 9bebf6a3, 300 cycles, 02-09; own keys: rad_vet/time2_gtime_slope=true rad_vet/clean_stag_thinf=true rad_vet/clean_stag_lat_sums=true rad_vet/clean_stag_tie_gid=true rad_vet/clean_stag_opac_drag=true rad_vet/vet_gd_allow_clamp=false hydro/vceil=3.0e8 rad_vet/clean_pc_sweep=rb rad_vet/clean_stag_mg_levels=4 time/cfl_number=0.3 rad_vet/clean_conv_max=1e-5; R_ph by rph.py; dumps at global cycles that are multiples of 100 and the end
- FW300: 2026-11-25 03:19:35 rc=0 (srun) ended: Terminating on cycle limit wall=352 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- NX300: 2026-11-25 02:36:59 rc=0 (srun) ended: Terminating on cycle limit wall=347 s CW 0 FATAL 0 inf 0 DT_RETRY 0
| arm | rows | s/cyc (elapsed prints) | dt min/mean/max | its (sum/cycle), solves | CW/FATAL/inf/DT_RETRY | fallback/BE |
|---|---|---|---|---|---|---|
| FW300 | 300 | 1.0435 (c22200-22400) | 12.0066/12.0068/12.0071 | 54.7, 9.11 | 0/0/0/0 | 0.000000e+00/0.000000e+00 |
| NX300 | 300 | 1.0262 (c22200-22400) | 12.0066/12.0068/12.0071 | 55.2, 9.10 | 0/0/0/0 | 0.000000e+00/0.000000e+00 |

| row | t (prototype) | L_top/L_in | rel L_top, L_bot, L_mid, L_int, KE, M_tot, Etot, dt (prototype / next - 1) |
|---|---|---|---|
| 0 | 2.592125e+05 | 1.013163 | -5.43e-07, 3.81e-10, -6.14e-13, -3.41e-02, 1.25e-07, 6.66e-16, 2.28e-13, 0.00e+00 |
| 9 | 2.593205e+05 | 1.013736 | 7.00e-07, 3.68e-09, -6.39e-12, -3.23e-02, 1.40e-06, -6.22e-15, 6.58e-12, 4.12e-11 |
| 49 | 2.598008e+05 | 1.014394 | -4.48e-05, -1.45e-09, 1.78e-11, -2.79e-02, 7.76e-06, 1.33e-15, 4.75e-11, -3.02e-11 |
| 99 | 2.604012e+05 | 1.014555 | -6.44e-06, 4.51e-09, -2.50e-10, -3.10e-02, 1.28e-05, -5.55e-16, 1.28e-10, 3.95e-11 |
| 199 | 2.616019e+05 | 1.016087 | -4.03e-05, 4.34e-09, -5.83e-08, -2.47e-02, 2.34e-05, -3.33e-16, 3.71e-10, 2.96e-11 |
| 299 | 2.628025e+05 | 1.009928 | -6.65e-06, -9.34e-09, -6.97e-07, -3.26e-03, 3.53e-05, 2.22e-15, 7.39e-10, 1.91e-11 |
| 299 | 2.628025e+05 | 1.009928 | -6.65e-06, -9.34e-09, -6.97e-07, -3.26e-03, 3.53e-05, 2.22e-15, 7.39e-10, 1.91e-11 |
- L_top: max |prototype / next - 1| over 300 rows 1.140e-04
- L_bot: max |prototype / next - 1| over 300 rows 3.848e-08
- L_mid: max |prototype / next - 1| over 300 rows 6.969e-07
- L_int: max |prototype / next - 1| over 300 rows 3.414e-02
- KE_int: max |prototype / next - 1| over 300 rows 3.526e-05

| dump | cycle | t | R_A (mean tau 2/3) | R_A / first production dump before - 1 | T_A [K] | R_M (median column) | T_M [K] |
|---|---|---|---|---|---|---|---|
| production before | 21989 | 2.574120e+05 | 9.583064e+11 | 3.61e-03 | 35908.10 | 9.488446e+11 | 35095.35 |
| production before | 22144 | 2.592000e+05 | 9.548562e+11 | 0.00e+00 | 36130.73 | 9.461009e+11 | 34870.73 |
| FW300 | 22200 | 2.598730e+05 | 9.533899e+11 | -1.54e-03 | 35883.25 | 9.448261e+11 | 34805.55 |
| FW300 | 22300 | 2.610740e+05 | 9.511163e+11 | -3.92e-03 | 35943.85 | 9.438555e+11 | 35023.04 |
| FW300 | 22400 | 2.622740e+05 | 9.516203e+11 | -3.39e-03 | 35925.43 | 9.443215e+11 | 35028.03 |
| FW300 | 22444 | 2.628030e+05 | 9.516405e+11 | -3.37e-03 | 35897.38 | 9.448090e+11 | 35062.34 |
| NX300 | 22200 | 2.598730e+05 | 9.533919e+11 | -1.53e-03 | 36135.07 | 9.448271e+11 | 35056.57 |
| NX300 | 22300 | 2.610740e+05 | 9.511234e+11 | -3.91e-03 | 36264.79 | 9.438694e+11 | 35420.74 |
| NX300 | 22400 | 2.622740e+05 | 9.516319e+11 | -3.38e-03 | 36086.03 | 9.443365e+11 | 35371.07 |
| NX300 | 22444 | 2.628030e+05 | 9.516561e+11 | -3.35e-03 | 35985.03 | 9.448206e+11 | 35356.92 |
| production after | 22294 | 2.610020e+05 | 9.511148e+11 | -3.92e-03 | 36275.66 | 9.438859e+11 | 35417.45 |
| production after | 22444 | 2.628030e+05 | 9.516561e+11 | -3.35e-03 | 35985.65 | 9.448203e+11 | 35353.79 |
| production after | 22594 | 2.646040e+05 | 9.511638e+11 | -3.87e-03 | 35972.37 | 9.462867e+11 | 35181.01 |

- FW300: CLEAN FCONE lines 184 (faces sum 1510, max 29 per line, cycles with a cut 121); FCONE_LZ lines 0; CLEAN FCS/bound lines: ['CLEAN FALLBACK', 'CLEAN FCONE', 'CLEAN FLOORMOM', 'CLEAN FWB', 'CLEAN LZB', 'CLEAN MCAP', 'CLEAN OLDVEC', 'CLEAN XT']; clamped lateral reads (all ranks, end) 0.000000e+00; CLEAN WARNING 0
- FW300: top mass loss Mdot_top mean 5.2557e-04 (first 1.0115e-04, last 3.1362e-04, min 8.0050e-05, max 4.6713e-03); cumulative M_tot change -1.554e-15; L_mid / L_midE - 1: mean -1.018e-12 max abs 1.029e-10; L_intE/L_int - 1 mean 3.167e-03
- NX300: CLEAN FCONE lines 0 (faces sum 0, max 0 per line, cycles with a cut 0); FCONE_LZ lines 0; CLEAN FCS/bound lines: ['CLEAN FALLBACK', 'CLEAN FLOORMOM']; clamped lateral reads (all ranks, end) 0.000000e+00; CLEAN WARNING 0
- NX300: top mass loss Mdot_top mean 6.3687e-03 (first 1.0228e-04, last 4.3068e-03, min 8.7175e-05, max 6.6191e-02); cumulative M_tot change -3.220e-15; L_mid / L_midE - 1: mean -4.577e-13 max abs 6.759e-11; L_intE/L_int - 1 mean -2.173e-02
- L_midE prototype - next: max |diff| / max |next| 6.969e-07, at the last row -6.969e-07
- Mdot_top prototype - next: max |diff| / max |next| 9.649e-01, at the last row -6.033e-02
- E_rad prototype - next: max |diff| / max |next| 2.948e-09, at the last row 2.948e-09
- e_gas prototype - next: max |diff| / max |next| 2.318e-10, at the last row 2.318e-10
- Lz_tot prototype - next: max |diff| / max |next| 3.939e-10, at the last row -2.352e-10
- FW300: L_top / L_in first 1.013163, last 1.009928, mean of the last 100 rows 1.013347
- NX300: L_top / L_in first 1.013164, last 1.009934, mean of the last 100 rows 1.013356
- FW300 CLEAN FCONE lines 184; first: CLEAN FCONE cycle 22144 t2s 3 faces 6 momentum cut sum |dF| V/c^2 1.945872e+20 | last: CLEAN FCONE cycle 22440 t2s 3 faces 6 momentum cut sum |dF| V/c^2 6.185857e+18
- FW300 CLEAN FWB lines 600; first: CLEAN FWB cycle 22144 t2s 2 faces 17132717 radiation momentum sum |dp| 7.699544e+26 x1 faces sum dp -1.958855e+25 | last: CLEAN FWB cycle 22443 t2s 3 faces 17277873 radiation momentum sum |dp| 6.802293e+26 x1 faces sum dp -3.231758e+24
- NX300 CLEAN FCONE lines 0; first: - | last: -
- NX300 CLEAN FWB lines 0; first: - | last: -
- mass-weighted shell-mean v_r [cm/s] of the last dump (c22444 FW300 / c22444 NX300): r: FW300 | NX300
  6.7509e+11: -1.646e+00 | -1.646e+00; 6.9838e+11: -9.478e+01 | -9.478e+01; 7.2344e+11: -1.438e+02 | -1.438e+02; 7.4941e+11: -7.538e+00 | -7.538e+00; 7.7631e+11: 1.738e+03 | 1.738e+03; 8.0309e+11: 8.125e+03 | 8.125e+03; 8.3191e+11: 1.220e+04 | 1.220e+04; 8.6177e+11: 5.289e+04 | 5.290e+04; 8.9270e+11: 1.317e+05 | 1.317e+05; 9.2349e+11: 1.345e+05 | 1.346e+05; 9.5664e+11: 6.150e+05 | 6.277e+05; 9.9098e+11: -8.101e+06 | -8.307e+06; 1.0265e+12: -4.279e+06 | -4.282e+06; 1.0634e+12: 1.190e+01 | -1.180e+02
- L_int / L_intE and L_mid / L_midE, Mdot_top at rows (0-based) 0, 1, 5, 20, 100, 299:
  - FW300: r0 Li 0.993262 Lm 1.000000 Mdot 1.011e-04; r1 Li 0.993642 Lm 1.000000 Mdot 1.021e-04; r5 Li 0.993817 Lm 1.000000 Mdot 9.010e-05; r20 Li 0.994822 Lm 1.000000 Mdot 1.014e-04; r100 Li 0.996834 Lm 1.000000 Mdot 9.506e-05; r299 Li 0.998153 Lm 1.000000 Mdot 3.136e-04
  - NX300: r0 Li 1.028366 Lm 1.000000 Mdot 1.023e-04; r1 Li 1.027791 Lm 1.000000 Mdot 1.017e-04; r5 Li 1.027503 Lm 1.000000 Mdot 9.490e-05; r20 Li 1.026659 Lm 1.000000 Mdot 1.048e-04; r100 Li 1.029151 Lm 1.000000 Mdot 9.422e-05; r299 Li 1.001500 Lm 1.000000 Mdot 4.307e-03
- FW300 CLEAN FCONE per 50 cycles (from c22144; mean faces per stage, mean cut): 0-49: 7, 5.337e+19; 50-99: 8, 2.771e+19; 100-149: 7, 4.728e+20; 150-199: 10, 1.009e+21; 200-249: 6, 6.206e+19; 250-299: 8, 1.583e+20
- top cell (shell mean) FW300 c22200: rho 1.0000e-17, v_r 5.2483e+01, tau (cell centre) 2.5031e-09
- top cell (shell mean) FW300 c22300: rho 1.0000e-17, v_r 6.5563e+01, tau (cell centre) 2.5030e-09
- top cell (shell mean) FW300 c22400: rho 1.0000e-17, v_r 3.7714e+01, tau (cell centre) 2.5028e-09
- top cell (shell mean) FW300 c22444: rho 1.0000e-17, v_r 1.1899e+01, tau (cell centre) 2.5027e-09
- top cell (shell mean) NX300 c22200: rho 1.0000e-17, v_r -1.5944e+01, tau (cell centre) 2.5031e-09
- top cell (shell mean) NX300 c22300: rho 1.0000e-17, v_r -5.6259e+00, tau (cell centre) 2.5030e-09
- top cell (shell mean) NX300 c22400: rho 1.0000e-17, v_r -1.0577e+02, tau (cell centre) 2.5028e-09
- top cell (shell mean) NX300 c22444: rho 1.0000e-17, v_r -1.1795e+02, tau (cell centre) 2.5027e-09
