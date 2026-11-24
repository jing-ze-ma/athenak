# NOTE 2026-11-24zzzj caltech -> orion: stop of RSG + 4 stars and the test batch (your 24zzzh / 24zzzi)

## U2 (11-24 14:30)
| arm | rows | s/cyc | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries/skipped | clamped reads line | L_top/L_in end | L_bot/L_in end | KE end | dM/M | dE/E | dLz/Lz |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| U2A | 300 | 0.2974 (c22160-22450) | 3.481 / 3.481 / 3.481 | 18.8, 2.00 | 17.6, 2.34 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 1.00009 | 1.00001 | 4.96109e+35 | -3.10e-15 | -2.29e-08 | 2.41e-01 |
| U2B | 300 | 0.3017 (c22160-22450) | 3.481 / 3.481 / 3.481 | 18.9, 2.00 | 17.4, 2.31 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 1.00009 | 1.00001 | 4.96108e+35 | -9.92e-16 | -2.29e-08 | 2.38e-01 |

- matched time (300 rows): L_top rel end -3.435e-07 max 6.369e-05; L_bot end -3.097e-09 max 6.372e-08; KE_int end -1.063e-06 max 1.063e-06
- dt mean ratio-1 1.126e-12; s per simulated s: A 0.08545 B 0.08669 ratio 1.015
- **HARMLESS** (bounds of 24zzzzl)
- U2A: 2026-11-24 14:28:41 rc=0 wall=119 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 8.735251e+01 | faile`; keys in the log: 52:<hydro> floor_audit cycle 22155 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de |53:<hydro> floor_audit cycle 22155 [refresh] cell-c2p 6291456 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) fl|73:<hydro> floor_audit cycle 22156 [update] cell-c2p 12582912 e_c<=0 0 kedom 0
- U2B: 2026-11-24 14:30:39 rc=0 wall=117 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 8.846771e+01 | faile`; keys in the log: 52:<hydro> floor_audit cycle 22155 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de |53:<hydro> floor_audit cycle 22155 [refresh] cell-c2p 6291456 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) fl|73:<hydro> floor_audit cycle 22156 [update] cell-c2p 12582912 e_c<=0 0 kedom 0

## U2: UA vs U2C (11-24 14:33)
| arm | rows | s/cyc | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries/skipped | clamped reads line | L_top/L_in end | L_bot/L_in end | KE end | dM/M | dE/E | dLz/Lz |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| U2A | 300 | 0.2974 (c22160-22450) | 3.481 / 3.481 / 3.481 | 18.8, 2.00 | 17.6, 2.34 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 1.00009 | 1.00001 | 4.96109e+35 | -3.10e-15 | -2.29e-08 | 2.41e-01 |
| U2C | 300 | 0.3024 (c22160-22450) | 3.481 / 3.481 / 3.481 | 18.8, 2.00 | 16.7, 2.25 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 1.00009 | 1.00001 | 4.96108e+35 | 1.61e-15 | -2.29e-08 | 2.36e-01 |

- matched time (300 rows): L_top rel end -1.355e-06 max 2.976e-06; L_bot end -2.546e-08 max 6.191e-08; KE_int end -2.725e-07 max 2.879e-07
- dt mean ratio-1 4.260e-12; s per simulated s: A 0.08545 B 0.08687 ratio 1.017
- **HARMLESS** (bounds of 24zzzzl)
- U2A: 2026-11-24 14:28:41 rc=0 wall=119 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 8.735251e+01 | faile`; keys in the log: 52:<hydro> floor_audit cycle 22155 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de |53:<hydro> floor_audit cycle 22155 [refresh] cell-c2p 6291456 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) fl|73:<hydro> floor_audit cycle 22156 [update] cell-c2p 12582912 e_c<=0 0 kedom 0
- U2C: 2026-11-24 14:33:10 rc=0 wall=121 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 8.838767e+01 | faile`; keys in the log: 52:<hydro> floor_audit cycle 22155 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de |53:<hydro> floor_audit cycle 22155 [refresh] cell-c2p 6291456 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) fl|73:<hydro> floor_audit cycle 22156 [update] cell-c2p 12582912 e_c<=0 0 kedom 0

## U1 (11-24 14:34)
| arm | rows | s/cyc | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries/skipped | clamped reads line | L_top/L_in end | L_bot/L_in end | KE end | dM/M | dE/E | dLz/Lz |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| U1A | 300 | 0.6348 (c6950-7200) | 146.7 / 152.8 / 164.7 | 15.0, 3.39 | 14.3, 3.92 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 1.15812 | 1.01394 | 1.31656e+42 | -2.95e-15 | -1.53e-05 | 5.10e-02 |
| U1B | 300 | 0.6529 (c6950-7200) | 146.7 / 152.8 / 164.7 | 15.9, 3.47 | 14.2, 3.92 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 1.15812 | 1.01394 | 1.31656e+42 | -2.33e-15 | -1.53e-05 | 5.29e-02 |

- matched time (300 rows): L_top rel end -1.239e-06 max 4.930e-04; L_bot end 7.085e-10 max 1.258e-08; KE_int end -7.949e-08 max 7.949e-08
- dt mean ratio-1 1.462e-09; s per simulated s: A 0.00415 B 0.00427 ratio 1.028
- **HARMLESS** (bounds of 24zzzzl)
- U1A: 2026-11-24 14:30:30 rc=0 wall=228 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 1.872151e+02 | faile`; keys in the log: 53:<hydro> floor_audit cycle 6902 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de 0|54:<hydro> floor_audit cycle 6902 [refresh] cell-c2p 16777216 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) fl|85:<hydro> floor_audit cycle 6903 [update] cell-c2p 33554432 e_c<=0 0 kedom 0 
- U1B: 2026-11-24 14:34:17 rc=0 wall=227 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 1.924788e+02 | faile`; keys in the log: 53:<hydro> floor_audit cycle 6902 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de 0|54:<hydro> floor_audit cycle 6902 [refresh] cell-c2p 16777216 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) fl|85:<hydro> floor_audit cycle 6903 [update] cell-c2p 33554432 e_c<=0 0 kedom 0 

## U3: UA vs U3C (11-24 14:38)
| arm | rows | s/cyc | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries/skipped | clamped reads line | L_top/L_in end | L_bot/L_in end | KE end | dM/M | dE/E | dLz/Lz |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| U3A | 300 | 1.0168 (c7000-7250) | 9.862 / 11.5 / 12.01 | 34.1, 4.95 | 16.1, 4.01 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 0.99958 | 1.00016 | 1.57167e+39 | -3.78e-16 | 9.32e-07 | 7.96e-07 |
| U3C | 300 | 1.0451 (c7000-7250) | 9.862 / 11.5 / 12.01 | 34.4, 4.98 | 16.1, 4.01 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 0.99959 | 1.00016 | 1.57167e+39 | 3.59e-15 | 9.32e-07 | 2.01e-05 |

- matched time (300 rows): L_top rel end 1.353e-05 max 1.550e-04; L_bot end -8.920e-10 max 6.483e-08; KE_int end 1.499e-08 max 1.862e-08
- dt mean ratio-1 3.629e-08; s per simulated s: A 0.08844 B 0.09090 ratio 1.028
- **HARMLESS** (bounds of 24zzzzl)
- U3A: 2026-11-24 14:32:31 end: t_U3A|COMPLETED|0:0|00:05:42; log: Terminating on cycle limit CW 0 FATAL 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 3.024931e+02 | faile`; keys in the log: 52:<hydro> floor_audit cycle 6959 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de 0|53:<hydro> floor_audit cycle 6959 [refresh] cell-c2p 22020096 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) fl|71:<hydro> floor_audit cycle 6960 [update] cell-c2p 44040192 e_c<=0 34 kedom 5
- U3C: 2026-11-24 14:38:21 rc=0 wall=350 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 3.088078e+02 | faile`; keys in the log: 52:<hydro> floor_audit cycle 6959 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de 0|53:<hydro> floor_audit cycle 6959 [refresh] cell-c2p 22020096 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) fl|71:<hydro> floor_audit cycle 6960 [update] cell-c2p 44040192 e_c<=0 34 kedom 5

## U1: UA vs U1C (11-24 14:38)
| arm | rows | s/cyc | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries/skipped | clamped reads line | L_top/L_in end | L_bot/L_in end | KE end | dM/M | dE/E | dLz/Lz |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| U1A | 300 | 0.6348 (c6950-7200) | 146.7 / 152.8 / 164.7 | 15.0, 3.39 | 14.3, 3.92 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 1.15812 | 1.01394 | 1.31656e+42 | -2.95e-15 | -1.53e-05 | 5.10e-02 |
| U1C | 300 | 0.6483 (c6950-7200) | 146.7 / 152.8 / 164.7 | 14.9, 3.39 | 14.2, 3.92 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 1.15812 | 1.01394 | 1.31656e+42 | -1.40e-15 | -1.53e-05 | 5.10e-02 |

- matched time (300 rows): L_top rel end 2.756e-08 max 5.346e-06; L_bot end -6.384e-09 max 1.137e-08; KE_int end -4.811e-10 max 4.811e-10
- dt mean ratio-1 7.715e-10; s per simulated s: A 0.00415 B 0.00424 ratio 1.021
- **HARMLESS** (bounds of 24zzzzl)
- U1A: 2026-11-24 14:30:30 rc=0 wall=228 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 1.872151e+02 | faile`; keys in the log: 53:<hydro> floor_audit cycle 6902 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de 0|54:<hydro> floor_audit cycle 6902 [refresh] cell-c2p 16777216 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) fl|85:<hydro> floor_audit cycle 6903 [update] cell-c2p 33554432 e_c<=0 0 kedom 0 
- U1C: 2026-11-24 14:38:37 rc=0 wall=228 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 1.912128e+02 | faile`; keys in the log: 53:<hydro> floor_audit cycle 6902 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de 0|54:<hydro> floor_audit cycle 6902 [refresh] cell-c2p 16777216 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) fl|85:<hydro> floor_audit cycle 6903 [update] cell-c2p 33554432 e_c<=0 0 kedom 0 

## U4 (11-24 14:42)
| arm | rows | s/cyc | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries/skipped | clamped reads line | L_top/L_in end | L_bot/L_in end | KE end | dM/M | dE/E | dLz/Lz |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| U4A | 300 | 2.9935 (c5300-5550) | 10.51 / 32.88 / 46.91 | 14.2, 4.84 | 8.4, 3.61 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 2.44083 | 0.98459 | 1.78357e+44 | 3.09e-16 | -6.01e-05 | -4.71e-06 |
| U4B | 300 | 3.0259 (c5300-5550) | 11.87 / 33.55 / 46.88 | 14.3, 4.82 | 8.4, 3.54 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 2.43851 | 0.98498 | 1.78425e+44 | -2.01e-15 | -6.13e-05 | -6.78e-06 |

- matched time (300 rows): L_top rel end 6.267e-05 max 1.146e-03; L_bot end 9.262e-07 max 2.987e-06; KE_int end -5.683e-07 max 9.815e-07
- dt mean ratio-1 2.062e-02; s per simulated s: A 0.09105 B 0.09018 ratio 0.990
- **NOT HARMLESS by the bounds: dt mean 2.06 % off** (bounds of 24zzzzl)
- U4A: 2026-11-24 14:42:42 rc=0 wall=960 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 8.886084e+02 | faile`; keys in the log: 28:### WARNING (scheme = clean): no effect: time/dt_retry_max rad_vet/clean_stag_thinf_held rad_vet/time2_start rad_vet/vet_col_lat_taucut rad_vet/time2_stage_l|58:<hydro> floor_audit cycle 5279 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de 0|59:<hydro> floor_audit cycle 5279 [refresh] cell-c2p 67108864 e_c<=0 0 kedom 8
- U4B: 2026-11-24 14:42:56 rc=0 wall=974 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 8.981109e+02 | faile`; keys in the log: 28:### WARNING (scheme = clean): no effect: time/dt_retry_max rad_vet/clean_stag_thinf_held rad_vet/time2_start rad_vet/vet_col_lat_taucut rad_vet/time2_stage_l|58:<hydro> floor_audit cycle 5279 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de 0|59:<hydro> floor_audit cycle 5279 [refresh] cell-c2p 67108864 e_c<=0 0 kedom 8
