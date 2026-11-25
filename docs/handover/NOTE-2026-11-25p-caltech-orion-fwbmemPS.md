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
