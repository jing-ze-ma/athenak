# NOTE 2026-11-25m caltech -> orion: star arms of the PROTOTYPE caltech-fwbmem-77d0c0ab7 against next 79c1cabef (your 25m; result-changing, not a release)

## bsg (11-25 02:33): rst bsg_bsg3d.00008.rst (the 25b pair's hard link, c29127), prototype adcfad3f against next 9bebf6a3, 300 cycles, 02-03; own keys: none; R_ph by rph.py; dumps at global cycles that are multiples of 100 and the end
- FW300: 2026-11-25 02:29:25 rc=0 (srun) ended: Terminating on cycle limit wall=240 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- NX300: 2026-11-25 02:33:17 rc=0 (srun) ended: Terminating on cycle limit wall=232 s CW 0 FATAL 0 inf 0 DT_RETRY 0
| arm | rows | s/cyc (elapsed prints) | dt min/mean/max | its (sum/cycle), solves | CW/FATAL/inf/DT_RETRY | fallback/BE |
|---|---|---|---|---|---|---|
| FW300 | 300 | 0.6871 (c29200-29400) | 160.3591/164.6545/169.4520 | 36.5, 6.44 | 0/0/0/0 | 0.000000e+00/0.000000e+00 |
| NX300 | 300 | 0.6725 (c29200-29400) | 160.3591/164.6542/169.4520 | 36.2, 6.37 | 0/0/0/0 | 0.000000e+00/0.000000e+00 |

| row | t (prototype) | L_top/L_in | rel L_top, L_bot, L_mid, L_int, KE, M_tot, Etot, dt (prototype / next - 1) |
|---|---|---|---|
| 0 | 4.838736e+06 | 1.019585 | -3.99e-08, -5.54e-09, -2.70e-12, 2.19e-02, 6.47e-10, -1.55e-15, 2.24e-14, 0.00e+00 |
| 9 | 4.840260e+06 | 1.019425 | -1.30e-06, 1.27e-10, 1.31e-09, 2.11e-02, 4.69e-09, 8.88e-16, 9.37e-13, 1.64e-10 |
| 49 | 4.846913e+06 | 1.018759 | 7.79e-07, -2.00e-09, 5.50e-08, 2.09e-02, -2.31e-08, -3.44e-15, -5.52e-12, 8.08e-08 |
| 99 | 4.854984e+06 | 1.017652 | 4.52e-05, -2.06e-08, 1.28e-07, 2.09e-02, -7.16e-08, 2.22e-16, -4.02e-11, 6.92e-08 |
| 199 | 4.871682e+06 | 1.014496 | 1.01e-04, -6.68e-06, 4.26e-06, 2.20e-02, -1.13e-07, 1.33e-15, -1.79e-10, 1.69e-05 |
| 299 | 4.887963e+06 | 1.011193 | 6.56e-06, -1.02e-08, 1.94e-06, 2.54e-02, -2.74e-07, 1.78e-15, -1.98e-10, -1.54e-07 |
| 299 | 4.887963e+06 | 1.011193 | 6.56e-06, -1.02e-08, 1.94e-06, 2.54e-02, -2.74e-07, 1.78e-15, -1.98e-10, -1.54e-07 |
- L_top: max |prototype / next - 1| over 300 rows 1.014e-04
- L_bot: max |prototype / next - 1| over 300 rows 6.754e-06
- L_mid: max |prototype / next - 1| over 300 rows 4.559e-06
- L_int: max |prototype / next - 1| over 300 rows 2.560e-02
- KE_int: max |prototype / next - 1| over 300 rows 2.736e-07

| dump | cycle | t | R_A (mean tau 2/3) | R_A / first production dump before - 1 | T_A [K] | R_M (median column) | T_M [K] |
|---|---|---|---|---|---|---|---|
| production before | 28863 | 4.795360e+06 | 3.649381e+12 | 3.69e-03 | 15710.70 | 3.633354e+12 | 15752.34 |
| production before | 29127 | 4.838570e+06 | 3.635978e+12 | 0.00e+00 | 15710.83 | 3.621075e+12 | 15747.99 |
| FW300 | 29200 | 4.850640e+06 | 3.638549e+12 | 7.07e-04 | 15735.41 | 3.619520e+12 | 15787.62 |
| FW300 | 29300 | 4.867140e+06 | 3.649508e+12 | 3.72e-03 | 15690.40 | 3.624350e+12 | 15731.66 |
| FW300 | 29400 | 4.883610e+06 | 3.660783e+12 | 6.82e-03 | 15645.43 | 3.633735e+12 | 15685.60 |
| FW300 | 29427 | 4.887960e+06 | 3.663015e+12 | 7.44e-03 | 15636.74 | 3.635485e+12 | 15674.17 |
| NX300 | 29200 | 4.850640e+06 | 3.638542e+12 | 7.05e-04 | 15703.03 | 3.619526e+12 | 15739.21 |
| NX300 | 29300 | 4.867140e+06 | 3.649476e+12 | 3.71e-03 | 15666.07 | 3.624299e+12 | 15679.19 |
| NX300 | 29400 | 4.883610e+06 | 3.660672e+12 | 6.79e-03 | 15569.61 | 3.633680e+12 | 15632.31 |
| NX300 | 29427 | 4.887960e+06 | 3.662876e+12 | 7.40e-03 | 15552.48 | 3.635371e+12 | 15627.31 |
| production after | 29388 | 4.881660e+06 | 3.659595e+12 | 6.50e-03 | 15577.71 | 3.633038e+12 | 15636.83 |
| production after | 29657 | 4.924800e+06 | 3.671902e+12 | 9.88e-03 | 15547.73 | 3.650319e+12 | 15576.49 |

- FW300: CLEAN FCONE lines 0 (faces sum 0, max 0 per line, cycles with a cut 0); FCONE_LZ lines 0; CLEAN FCS/bound lines: ['CLEAN FALLBACK', 'CLEAN FLOORMOM']; clamped lateral reads (all ranks, end) 0.000000e+00; CLEAN WARNING 0
- FW300: top mass loss Mdot_top mean 1.0377e+08 (first 6.6480e+08, last 2.0453e+07, min 5.3682e+06, max 6.6480e+08); cumulative M_tot change 4.441e-15; L_mid / L_midE - 1: mean 1.647e-09 max abs 5.743e-09; L_intE/L_int - 1 mean -1.764e-02
- NX300: CLEAN FCONE lines 0 (faces sum 0, max 0 per line, cycles with a cut 0); FCONE_LZ lines 0; CLEAN FCS/bound lines: ['CLEAN FALLBACK', 'CLEAN FLOORMOM']; clamped lateral reads (all ranks, end) 0.000000e+00; CLEAN WARNING 0
- NX300: top mass loss Mdot_top mean 1.0354e+08 (first 6.6107e+08, last 1.6998e+07, min 5.5008e+06, max 6.6107e+08); cumulative M_tot change 1.332e-15; L_mid / L_midE - 1: mean 1.873e-09 max abs 6.180e-09; L_intE/L_int - 1 mean 3.796e-03
- L_midE prototype - next: max |diff| / max |next| 4.535e-06, at the last row 1.930e-06
- Mdot_top prototype - next: max |diff| / max |next| 9.763e-03, at the last row 5.227e-03
- E_rad prototype - next: max |diff| / max |next| 1.857e-09, at the last row -1.857e-09
- e_gas prototype - next: max |diff| / max |next| 4.774e-10, at the last row -4.774e-10
- Lz_tot prototype - next: max |diff| / max |next| 3.508e-02, at the last row 3.508e-02
- FW300: L_top / L_in first 1.019585, last 1.011193, mean of the last 100 rows 1.012750
- NX300: L_top / L_in first 1.019585, last 1.011187, mean of the last 100 rows 1.012749
- FW300 CLEAN FCONE lines 0; first: - | last: -
- FW300 CLEAN FWB lines 0; first: - | last: -
- NX300 CLEAN FCONE lines 0; first: - | last: -
- NX300 CLEAN FWB lines 0; first: - | last: -
- mass-weighted shell-mean v_r [cm/s] of the last dump (c29427 FW300 / c29427 NX300): r: FW300 | NX300
  1.4023e+12: 3.163e+00 | 3.163e+00; 1.8154e+12: 2.942e+03 | 2.942e+03; 2.2502e+12: 7.126e+03 | 7.126e+03; 2.5491e+12: 8.976e+03 | 8.976e+03; 2.7665e+12: 4.006e+03 | 4.006e+03; 2.9839e+12: 7.228e+03 | 7.228e+03; 3.1904e+12: 3.963e+04 | 3.962e+04; 3.4078e+12: 1.392e+05 | 1.392e+05; 3.6144e+12: 3.942e+05 | 3.921e+05; 3.8372e+12: 1.393e+05 | 1.312e+05; 4.2721e+12: 3.562e+03 | 3.467e+03; 4.6851e+12: 2.886e+03 | 2.814e+03; 5.1199e+12: 2.378e+03 | 2.322e+03; 5.5547e+12: 2.006e+03 | 1.979e+03
