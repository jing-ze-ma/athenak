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

## m35mms (11-25 02:37): rst m35mms_m35mms.00004.rst (the 25b pair's hard link, c22144), prototype adcfad3f against next 9bebf6a3, 300 cycles, 02-09; own keys: rad_vet/time2_gtime_slope=true rad_vet/clean_stag_thinf=true rad_vet/clean_stag_lat_sums=true rad_vet/clean_stag_tie_gid=true rad_vet/clean_stag_opac_drag=true rad_vet/vet_gd_allow_clamp=false hydro/vceil=3.0e8 rad_vet/clean_pc_sweep=rb rad_vet/clean_stag_mg_levels=4 time/cfl_number=0.3 rad_vet/clean_conv_max=1e-5; R_ph by rph.py; dumps at global cycles that are multiples of 100 and the end
- FW300: 2026-11-25 02:31:12 rc=0 (srun) ended: Terminating on cycle limit wall=347 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- NX300: 2026-11-25 02:36:59 rc=0 (srun) ended: Terminating on cycle limit wall=347 s CW 0 FATAL 0 inf 0 DT_RETRY 0
| arm | rows | s/cyc (elapsed prints) | dt min/mean/max | its (sum/cycle), solves | CW/FATAL/inf/DT_RETRY | fallback/BE |
|---|---|---|---|---|---|---|
| FW300 | 300 | 1.0440 (c22200-22400) | 12.0066/12.0068/12.0071 | 54.9, 9.11 | 0/0/0/0 | 0.000000e+00/0.000000e+00 |
| NX300 | 300 | 1.0262 (c22200-22400) | 12.0066/12.0068/12.0071 | 55.2, 9.10 | 0/0/0/0 | 0.000000e+00/0.000000e+00 |

| row | t (prototype) | L_top/L_in | rel L_top, L_bot, L_mid, L_int, KE, M_tot, Etot, dt (prototype / next - 1) |
|---|---|---|---|
| 0 | 2.592125e+05 | 1.013164 | -1.49e-07, 8.31e-10, -3.08e-12, -2.57e-02, 7.82e-08, -5.55e-16, 6.06e-14, 0.00e+00 |
| 9 | 2.593205e+05 | 1.013736 | 2.69e-07, -4.31e-09, -1.35e-11, -2.41e-02, 9.21e-07, -2.66e-15, 3.01e-12, 4.01e-10 |
| 49 | 2.598008e+05 | 1.014437 | -2.91e-06, -1.28e-09, 1.08e-10, -2.02e-02, 5.06e-06, -2.66e-15, 1.64e-11, -1.45e-10 |
| 99 | 2.604012e+05 | 1.014577 | 1.55e-05, 2.09e-09, -2.17e-10, -2.29e-02, 6.84e-06, -1.11e-15, 6.29e-11, 5.56e-10 |
| 199 | 2.616019e+05 | 1.016109 | -1.82e-05, 3.89e-09, -1.93e-08, -1.48e-02, 9.62e-06, 1.33e-15, 1.68e-10, 1.90e-10 |
| 299 | 2.628025e+05 | 1.009937 | 2.97e-06, -5.38e-10, -1.69e-07, 6.09e-03, 2.33e-05, -2.22e-16, 2.68e-10, -8.12e-11 |
| 299 | 2.628025e+05 | 1.009937 | 2.97e-06, -5.38e-10, -1.69e-07, 6.09e-03, 2.33e-05, -2.22e-16, 2.68e-10, -8.12e-11 |
- L_top: max |prototype / next - 1| over 300 rows 1.215e-04
- L_bot: max |prototype / next - 1| over 300 rows 1.889e-08
- L_mid: max |prototype / next - 1| over 300 rows 1.688e-07
- L_int: max |prototype / next - 1| over 300 rows 2.569e-02
- KE_int: max |prototype / next - 1| over 300 rows 2.334e-05

| dump | cycle | t | R_A (mean tau 2/3) | R_A / first production dump before - 1 | T_A [K] | R_M (median column) | T_M [K] |
|---|---|---|---|---|---|---|---|
| production before | 21989 | 2.574120e+05 | 9.583064e+11 | 3.61e-03 | 35908.10 | 9.488446e+11 | 35095.35 |
| production before | 22144 | 2.592000e+05 | 9.548562e+11 | 0.00e+00 | 36130.73 | 9.461009e+11 | 34870.73 |
| FW300 | 22200 | 2.598730e+05 | 9.533904e+11 | -1.54e-03 | 35955.23 | 9.448265e+11 | 34905.79 |
| FW300 | 22300 | 2.610740e+05 | 9.511235e+11 | -3.91e-03 | 36039.63 | 9.438684e+11 | 35186.79 |
| FW300 | 22400 | 2.622740e+05 | 9.516321e+11 | -3.38e-03 | 36011.05 | 9.443353e+11 | 35224.78 |
| FW300 | 22444 | 2.628030e+05 | 9.516559e+11 | -3.35e-03 | 35990.88 | 9.448226e+11 | 35235.15 |
| NX300 | 22200 | 2.598730e+05 | 9.533919e+11 | -1.53e-03 | 36135.07 | 9.448271e+11 | 35056.57 |
| NX300 | 22300 | 2.610740e+05 | 9.511234e+11 | -3.91e-03 | 36264.79 | 9.438694e+11 | 35420.74 |
| NX300 | 22400 | 2.622740e+05 | 9.516319e+11 | -3.38e-03 | 36086.03 | 9.443365e+11 | 35371.07 |
| NX300 | 22444 | 2.628030e+05 | 9.516561e+11 | -3.35e-03 | 35985.03 | 9.448206e+11 | 35356.92 |
| production after | 22294 | 2.610020e+05 | 9.511148e+11 | -3.92e-03 | 36275.66 | 9.438859e+11 | 35417.45 |
| production after | 22444 | 2.628030e+05 | 9.516561e+11 | -3.35e-03 | 35985.65 | 9.448203e+11 | 35353.79 |
| production after | 22594 | 2.646040e+05 | 9.511638e+11 | -3.87e-03 | 35972.37 | 9.462867e+11 | 35181.01 |

- FW300: CLEAN FCONE lines 0 (faces sum 0, max 0 per line, cycles with a cut 0); FCONE_LZ lines 0; CLEAN FCS/bound lines: ['CLEAN FALLBACK', 'CLEAN FLOORMOM']; clamped lateral reads (all ranks, end) 0.000000e+00; CLEAN WARNING 0
- FW300: top mass loss Mdot_top mean 1.3447e-03 (first 9.9015e-05, last 1.4349e-03, min 8.7070e-05, max 1.2614e-02); cumulative M_tot change -2.887e-15; L_mid / L_midE - 1: mean -6.833e-13 max abs 2.306e-12; L_intE/L_int - 1 mean -5.684e-03
- NX300: CLEAN FCONE lines 0 (faces sum 0, max 0 per line, cycles with a cut 0); FCONE_LZ lines 0; CLEAN FCS/bound lines: ['CLEAN FALLBACK', 'CLEAN FLOORMOM']; clamped lateral reads (all ranks, end) 0.000000e+00; CLEAN WARNING 0
- NX300: top mass loss Mdot_top mean 6.3687e-03 (first 1.0228e-04, last 4.3068e-03, min 8.7175e-05, max 6.6191e-02); cumulative M_tot change -3.220e-15; L_mid / L_midE - 1: mean -4.577e-13 max abs 6.759e-11; L_intE/L_int - 1 mean -2.173e-02
- L_midE prototype - next: max |diff| / max |next| 1.688e-07, at the last row -1.688e-07
- Mdot_top prototype - next: max |diff| / max |next| 9.791e-01, at the last row -4.339e-02
- E_rad prototype - next: max |diff| / max |next| 4.657e-10, at the last row 2.815e-10
- e_gas prototype - next: max |diff| / max |next| 2.321e-11, at the last row -7.162e-12
- Lz_tot prototype - next: max |diff| / max |next| 5.412e-10, at the last row 2.735e-10
- FW300: L_top / L_in first 1.013164, last 1.009937, mean of the last 100 rows 1.013354
- NX300: L_top / L_in first 1.013164, last 1.009934, mean of the last 100 rows 1.013356
- FW300 CLEAN FCONE lines 0; first: - | last: -
- FW300 CLEAN FWB lines 0; first: - | last: -
- NX300 CLEAN FCONE lines 0; first: - | last: -
- NX300 CLEAN FWB lines 0; first: - | last: -
- mass-weighted shell-mean v_r [cm/s] of the last dump (c22444 FW300 / c22444 NX300): r: FW300 | NX300
  6.7509e+11: -1.646e+00 | -1.646e+00; 6.9838e+11: -9.478e+01 | -9.478e+01; 7.2344e+11: -1.438e+02 | -1.438e+02; 7.4941e+11: -7.538e+00 | -7.538e+00; 7.7631e+11: 1.738e+03 | 1.738e+03; 8.0309e+11: 8.125e+03 | 8.125e+03; 8.3191e+11: 1.220e+04 | 1.220e+04; 8.6177e+11: 5.290e+04 | 5.290e+04; 8.9270e+11: 1.317e+05 | 1.317e+05; 9.2349e+11: 1.347e+05 | 1.346e+05; 9.5664e+11: 6.531e+05 | 6.277e+05; 9.9098e+11: -8.043e+06 | -8.307e+06; 1.0265e+12: -4.278e+06 | -4.282e+06; 1.0634e+12: 1.822e+01 | -1.180e+02

## t9l62 (11-25 02:46): rst t9l62_t9l62.00003.rst (the 25b pair's hard link, c6332), prototype adcfad3f against next 9bebf6a3, 300 cycles, 02-16; own keys: rad_vet/time2_gtime_slope=true rad_vet/clean_stag_thinf=true rad_vet/clean_stag_lat_sums=true rad_vet/clean_stag_tie_gid=true rad_vet/clean_stag_opac_drag=true rad_vet/vet_gd_allow_clamp=false hydro/vceil=3.0e8 rad_vet/clean_pc_sweep=rb rad_vet/clean_stag_mg_levels=2 time/cfl_number=0.3 rad_vet/clean_conv_max=1e-5; R_ph by rph.py; dumps at global cycles that are multiples of 100 and the end
- FW300: 2026-11-25 02:31:21 rc=0 (srun) ended: Terminating on cycle limit wall=866 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- NX300: 2026-11-25 02:46:26 rc=0 (srun) ended: Terminating on cycle limit wall=890 s CW 0 FATAL 0 inf 0 DT_RETRY 0
| arm | rows | s/cyc (elapsed prints) | dt min/mean/max | its (sum/cycle), solves | CW/FATAL/inf/DT_RETRY | fallback/BE |
|---|---|---|---|---|---|---|
| FW300 | 300 | 2.6443 (c6400-6600) | 54.2218/58.9475/66.0362 | 17.9, 6.77 | 0/0/0/0 | 0.000000e+00/0.000000e+00 |
| NX300 | 300 | 2.6098 (c6400-6600) | 54.2218/58.9475/66.0362 | 17.9, 6.79 | 0/0/0/0 | 0.000000e+00/0.000000e+00 |

| row | t (prototype) | L_top/L_in | rel L_top, L_bot, L_mid, L_int, KE, M_tot, Etot, dt (prototype / next - 1) |
|---|---|---|---|
| 0 | 5.184690e+05 | 2.141651 | 1.49e-13, -1.04e-08, 5.23e-02, 4.14e-13, -3.44e-15, 8.88e-16, 0.00e+00, 0.00e+00 |
| 9 | 5.189596e+05 | 2.138847 | 1.04e-07, 6.44e-09, 4.64e-02, 6.84e-09, 3.66e-12, 1.33e-15, 1.08e-13, 2.72e-09 |
| 49 | 5.211484e+05 | 2.127227 | -2.57e-05, 6.62e-09, 4.21e-02, -1.25e-09, 1.01e-10, 2.66e-15, 1.33e-10, -3.49e-08 |
| 99 | 5.240140e+05 | 2.113628 | 1.04e-06, 9.90e-10, 3.48e-02, 1.59e-11, -4.28e-10, 2.22e-16, 4.58e-10, 2.42e-08 |
| 199 | 5.298664e+05 | 2.090395 | -1.73e-05, -3.92e-09, 1.95e-02, 7.82e-09, 1.35e-09, 1.33e-15, 1.10e-09, 3.02e-07 |
| 299 | 5.360989e+05 | 2.073137 | -2.46e-05, -1.51e-08, 8.59e-03, 1.96e-08, 4.56e-10, 2.66e-15, 1.81e-09, -4.81e-07 |
| 299 | 5.360989e+05 | 2.073137 | -2.46e-05, -1.51e-08, 8.59e-03, 1.96e-08, 4.56e-10, 2.66e-15, 1.81e-09, -4.81e-07 |
- L_top: max |prototype / next - 1| over 300 rows 4.001e-05
- L_bot: max |prototype / next - 1| over 300 rows 2.741e-05
- L_mid: max |prototype / next - 1| over 300 rows 5.228e-02
- L_int: max |prototype / next - 1| over 300 rows 9.944e-07
- KE_int: max |prototype / next - 1| over 300 rows 1.372e-08

| dump | cycle | t | R_A (mean tau 2/3) | R_A / first production dump before - 1 | T_A [K] | R_M (median column) | T_M [K] |
|---|---|---|---|---|---|---|---|
| production before | 5943 | 5.076050e+05 | 9.520301e+12 | -2.67e-02 | 16961.93 | 9.185541e+12 | 17071.23 |
| production before | 6332 | 5.184150e+05 | 9.781226e+12 | 0.00e+00 | 16582.88 | 9.493198e+12 | 16736.69 |
| FW300 | 6400 | 5.221580e+05 | 9.871818e+12 | 9.26e-03 | 16671.86 | 9.595472e+12 | 16726.00 |
| FW300 | 6500 | 5.280440e+05 | 1.001635e+13 | 2.40e-02 | 16514.53 | 9.758387e+12 | 16569.36 |
| FW300 | 6600 | 5.340280e+05 | 1.016495e+13 | 3.92e-02 | 16367.87 | 9.918913e+12 | 16427.53 |
| FW300 | 6632 | 5.360990e+05 | 1.021647e+13 | 4.45e-02 | 16320.59 | 9.974351e+12 | 16378.00 |
| NX300 | 6400 | 5.221580e+05 | 9.871817e+12 | 9.26e-03 | 16475.01 | 9.595464e+12 | 16630.06 |
| NX300 | 6500 | 5.280440e+05 | 1.001634e+13 | 2.40e-02 | 16312.77 | 9.758380e+12 | 16470.90 |
| NX300 | 6600 | 5.340280e+05 | 1.016493e+13 | 3.92e-02 | 16166.17 | 9.918859e+12 | 16329.88 |
| NX300 | 6632 | 5.360990e+05 | 1.021645e+13 | 4.45e-02 | 16117.25 | 9.974319e+12 | 16285.82 |
| production after | 6521 | 5.292430e+05 | 1.004729e+13 | 2.72e-02 | 16280.56 | 9.789622e+12 | 16444.71 |
| production after | 6695 | 5.400260e+05 | 1.031637e+13 | 5.47e-02 | 16034.17 | 1.008700e+13 | 16201.88 |
| production after | 6884 | 5.508170e+05 | 1.059192e+13 | 8.29e-02 | 15833.38 | 1.037768e+13 | 15996.14 |

- FW300: CLEAN FCONE lines 0 (faces sum 0, max 0 per line, cycles with a cut 0); FCONE_LZ lines 0; CLEAN FCS/bound lines: ['CLEAN FALLBACK', 'CLEAN FLOORMOM']; clamped lateral reads (all ranks, end) 0.000000e+00; CLEAN WARNING 0
- FW300: top mass loss Mdot_top mean 3.8629e+05 (first 2.2054e+05, last 3.4629e+05, min 1.5881e+05, max 1.0265e+06); cumulative M_tot change 0.000e+00; L_mid / L_midE - 1: mean 1.888e-02 max abs 3.060e-02; L_intE/L_int - 1 mean 6.753e-08
- NX300: CLEAN FCONE lines 0 (faces sum 0, max 0 per line, cycles with a cut 0); FCONE_LZ lines 0; CLEAN FCS/bound lines: ['CLEAN FALLBACK', 'CLEAN FLOORMOM']; clamped lateral reads (all ranks, end) 0.000000e+00; CLEAN WARNING 0
- NX300: top mass loss Mdot_top mean 2.5157e+05 (first 2.1177e+05, last 2.5731e+05, min 1.2669e+05, max 6.2375e+05); cumulative M_tot change -1.665e-15; L_mid / L_midE - 1: mean -8.393e-03 max abs 2.092e-02; L_intE/L_int - 1 mean 6.879e-08
- L_midE prototype - next: max |diff| / max |next| 6.949e-04, at the last row 3.154e-04
- Mdot_top prototype - next: max |diff| / max |next| 1.182e+00, at the last row 1.427e-01
- E_rad prototype - next: max |diff| / max |next| 9.911e-10, at the last row -9.289e-10
- e_gas prototype - next: max |diff| / max |next| 2.942e-11, at the last row -2.832e-11
- Lz_tot prototype - next: max |diff| / max |next| 6.908e-06, at the last row 3.132e-06
- FW300: L_top / L_in first 2.141651, last 2.073137, mean of the last 100 rows 2.081226
- NX300: L_top / L_in first 2.141651, last 2.073188, mean of the last 100 rows 2.081254
- FW300 CLEAN FCONE lines 0; first: - | last: -
- FW300 CLEAN FWB lines 0; first: - | last: -
- NX300 CLEAN FCONE lines 0; first: - | last: -
- NX300 CLEAN FWB lines 0; first: - | last: -
- mass-weighted shell-mean v_r [cm/s] of the last dump (c6632 FW300 / c6632 NX300): r: FW300 | NX300
  1.6909e+12: 2.463e+02 | 2.463e+02; 2.2091e+12: 2.348e+04 | 2.348e+04; 2.8860e+12: -6.631e+05 | -6.631e+05; 3.7704e+12: -2.393e+06 | -2.393e+06; 4.9597e+12: 5.514e+05 | 5.514e+05; 6.4795e+12: 7.871e+06 | 7.871e+06; 8.4651e+12: 2.118e+07 | 2.118e+07; 1.1135e+13: 8.824e+02 | 7.924e+02; 1.4548e+13: 5.094e+02 | 4.622e+02; 1.9006e+13: 2.922e+02 | 2.662e+02; 2.5000e+13: 1.474e+02 | 1.339e+02; 3.2662e+13: 7.431e+01 | 6.902e+01; 4.2670e+13: 3.848e+01 | 3.555e+01; 5.6130e+13: 1.182e+01 | 1.273e+01

## FCONE / FWB lines (correction: clean_verbose defaults to 0, so the 300-cycle logs print none in either arm and say nothing about cuts): 20 cycles with rad_vet/clean_verbose=1, same rst and keys, 02-16
- t9l62 FW: t9l62 FW rc=0 Terminating on cycle limit CW 0 FATAL 0; FCONE lines 40, FWB lines 40
  - first: CLEAN FCONE cycle 6332 t2s 2 faces 52675329 momentum cut sum |dF| V/c^2 3.412355e+30 | last: CLEAN FCONE cycle 6351 t2s 3 faces 62941253 momentum cut sum |dF| V/c^2 4.180600e+30
  - first: CLEAN FWB cycle 6332 t2s 2 faces 102417981 radiation momentum sum |dp| 9.804047e+30 x1 faces sum dp 5.969097e+30 | last: CLEAN FWB cycle 6351 t2s 3 faces 102247128 radiation momentum sum |dp| 4.753694e+30 x1 faces sum dp 3.508625e+30
  - user.hst against the 300-cycle arm, 20 rows: printed-bitwise
- t9l62 NX: t9l62 NX rc=0 Terminating on cycle limit CW 0 FATAL 0; FCONE lines 40, FWB lines 40
  - first: CLEAN FCONE cycle 6332 t2s 2 faces 28147517 momentum cut sum |dF| V/c^2 1.840912e+29 | last: CLEAN FCONE cycle 6351 t2s 3 faces 28057357 momentum cut sum |dF| V/c^2 1.843953e+29
  - first: CLEAN FWB cycle 6332 t2s 2 faces 135651843 radiation momentum sum |dp| 1.185496e+30 x1 faces sum dp -9.381112e+28 | last: CLEAN FWB cycle 6351 t2s 3 faces 135537483 radiation momentum sum |dp| 1.184113e+30 x1 faces sum dp -9.544726e+28
  - user.hst against the 300-cycle arm, 20 rows: printed-bitwise
- m35mms FW: m35mms FW rc=0 Terminating on cycle limit CW 0 FATAL 0; FCONE lines 40, FWB lines 40
  - first: CLEAN FCONE cycle 22144 t2s 2 faces 163 momentum cut sum |dF| V/c^2 1.965189e+22 | last: CLEAN FCONE cycle 22163 t2s 3 faces 112 momentum cut sum |dF| V/c^2 4.477228e+21
  - first: CLEAN FWB cycle 22144 t2s 2 faces 17132717 radiation momentum sum |dp| 7.855995e+26 x1 faces sum dp -2.253156e+25 | last: CLEAN FWB cycle 22163 t2s 3 faces 17201938 radiation momentum sum |dp| 6.538017e+26 x1 faces sum dp -5.823584e+24
  - user.hst against the 300-cycle arm, 20 rows: printed-bitwise
- m35mms NX: m35mms NX rc=0 Terminating on cycle limit CW 0 FATAL 0; FCONE lines 40, FWB lines 40
  - first: CLEAN FCONE cycle 22144 t2s 2 faces 1439831 momentum cut sum |dF| V/c^2 1.318558e+25 | last: CLEAN FCONE cycle 22163 t2s 3 faces 1376004 momentum cut sum |dF| V/c^2 1.246290e+25
  - first: CLEAN FWB cycle 22144 t2s 2 faces 33481724 radiation momentum sum |dp| 6.263259e+26 x1 faces sum dp 7.979448e+24 | last: CLEAN FWB cycle 22163 t2s 3 faces 33528246 radiation momentum sum |dp| 6.289242e+26 x1 faces sum dp 7.800571e+24
  - user.hst against the 300-cycle arm, 20 rows: printed-bitwise
- bsg FW: bsg FW rc=0 Terminating on cycle limit CW 0 FATAL 0; FCONE lines 0, FWB lines 40
  - first:  | last: 
  - first: CLEAN FWB cycle 29127 t2s 2 faces 19661297 radiation momentum sum |dp| 2.356722e+28 x1 faces sum dp 1.357040e+27 | last: CLEAN FWB cycle 29146 t2s 3 faces 19663000 radiation momentum sum |dp| 2.027297e+28 x1 faces sum dp 2.259849e+25
  - user.hst against the 300-cycle arm, 20 rows: printed-bitwise
- bsg NX: bsg NX rc=0 Terminating on cycle limit CW 0 FATAL 0; FCONE lines 40, FWB lines 40
  - first: CLEAN FCONE cycle 29127 t2s 2 faces 271216 momentum cut sum |dF| V/c^2 6.721382e+25 | last: CLEAN FCONE cycle 29146 t2s 3 faces 261107 momentum cut sum |dF| V/c^2 6.695199e+25
  - first: CLEAN FWB cycle 29127 t2s 2 faces 29933903 radiation momentum sum |dp| 2.058573e+28 x1 faces sum dp 9.326646e+25 | last: CLEAN FWB cycle 29146 t2s 3 faces 29935074 radiation momentum sum |dp| 2.020473e+28 x1 faces sum dp 1.075372e+26
  - user.hst against the 300-cycle arm, 20 rows: printed-bitwise
