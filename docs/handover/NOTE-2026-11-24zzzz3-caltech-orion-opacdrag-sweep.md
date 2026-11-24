# NOTE 2026-11-24zzzj caltech -> orion: stop of RSG + 4 stars and the test batch (your 24zzzh / 24zzzi)

## PDz: UA vs PDzB (11-24 14:53)
| arm | rows | s/cyc | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries/skipped | clamped reads line | L_top/L_in end | L_bot/L_in end | KE end | dM/M | dE/E | dLz/Lz |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PDzA | 300 | 0.3024 (c22160-22450) | 3.481 / 3.481 / 3.481 | 18.8, 2.00 | 16.7, 2.25 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 1.00009 | 1.00001 | 4.96108e+35 | 1.61e-15 | -2.29e-08 | 2.36e-01 |
| PDzB | 300 | 0.3031 (c22160-22450) | 3.481 / 3.481 / 3.481 | 18.7, 2.00 | 17.7, 2.32 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 1.00009 | 1.00001 | 4.96108e+35 | 1.74e-15 | -2.29e-08 | 2.38e-01 |

- matched time (299 rows): L_top rel end -1.305e-06 max 2.686e-06; L_bot end -5.995e-10 max 5.977e-08; KE_int end 5.237e-08 max 1.439e-07
- dt mean ratio-1 -9.887e-12; s per simulated s: A 0.08687 B 0.08707 ratio 1.002
- **HARMLESS** (bounds of 24zzzzl)
- PDzA: 2026-11-24 14:33:10 rc=0 wall=121 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 8.838767e+01 | faile`; keys in the log: 52:<hydro> floor_audit cycle 22155 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de |53:<hydro> floor_audit cycle 22155 [refresh] cell-c2p 6291456 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) fl|73:<hydro> floor_audit cycle 22156 [update] cell-c2p 12582912 e_c<=0 0 kedom 0
- PDzB: 2026-11-24 14:53:41 rc=0 wall=115 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 9.021099e+01 | faile`; keys in the log: 52:<hydro> floor_audit cycle 22155 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de |53:<hydro> floor_audit cycle 22155 [refresh] cell-c2p 6291456 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) fl|73:<hydro> floor_audit cycle 22156 [update] cell-c2p 12582912 e_c<=0 0 kedom 0

## PDb: UA vs PDbB (11-24 14:55)
| arm | rows | s/cyc | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries/skipped | clamped reads line | L_top/L_in end | L_bot/L_in end | KE end | dM/M | dE/E | dLz/Lz |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PDbA | 300 | 0.6483 (c6950-7200) | 146.7 / 152.8 / 164.7 | 14.9, 3.39 | 14.2, 3.92 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 1.15812 | 1.01394 | 1.31656e+42 | -1.40e-15 | -1.53e-05 | 5.10e-02 |
| PDbB | 300 | 0.6617 (c6950-7200) | 146.7 / 152.8 / 164.7 | 15.0, 3.39 | 14.2, 3.92 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 1.15812 | 1.01394 | 1.31656e+42 | -3.42e-15 | -1.53e-05 | 5.11e-02 |

- matched time (299 rows): L_top rel end -3.974e-08 max 1.254e-05; L_bot end 1.550e-09 max 1.393e-08; KE_int end -4.154e-11 max 9.155e-11
- dt mean ratio-1 -6.304e-10; s per simulated s: A 0.00424 B 0.00433 ratio 1.021
- **HARMLESS** (bounds of 24zzzzl)
- PDbA: 2026-11-24 14:38:37 rc=0 wall=228 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 1.912128e+02 | faile`; keys in the log: 53:<hydro> floor_audit cycle 6902 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de 0|54:<hydro> floor_audit cycle 6902 [refresh] cell-c2p 16777216 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) fl|85:<hydro> floor_audit cycle 6903 [update] cell-c2p 33554432 e_c<=0 0 kedom 0 
- PDbB: 2026-11-24 14:55:33 rc=0 wall=227 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 1.958953e+02 | faile`; keys in the log: 53:<hydro> floor_audit cycle 6902 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de 0|54:<hydro> floor_audit cycle 6902 [refresh] cell-c2p 16777216 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) fl|85:<hydro> floor_audit cycle 6903 [update] cell-c2p 33554432 e_c<=0 0 kedom 0 

## PDm: UA vs PDmB (11-24 14:57)
| arm | rows | s/cyc | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries/skipped | clamped reads line | L_top/L_in end | L_bot/L_in end | KE end | dM/M | dE/E | dLz/Lz |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PDmA | 300 | 1.0451 (c7000-7250) | 9.862 / 11.5 / 12.01 | 34.4, 4.98 | 16.1, 4.01 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 0.99959 | 1.00016 | 1.57167e+39 | 3.59e-15 | 9.32e-07 | 2.01e-05 |
| PDmB | 300 | 1.0552 (c7000-7250) | 9.862 / 11.5 / 12.01 | 34.2, 5.00 | 16.0, 4.01 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 0.99960 | 1.00016 | 1.57167e+39 | -2.46e-15 | 9.32e-07 | 1.28e-05 |

- matched time (300 rows): L_top rel end 1.073e-05 max 2.157e-04; L_bot end -3.033e-09 max 6.681e-08; KE_int end -1.796e-08 max 2.113e-08
- dt mean ratio-1 4.234e-09; s per simulated s: A 0.09090 B 0.09178 ratio 1.010
- **HARMLESS** (bounds of 24zzzzl)
- PDmA: 2026-11-24 14:38:21 rc=0 wall=350 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 3.088078e+02 | faile`; keys in the log: 52:<hydro> floor_audit cycle 6959 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de 0|53:<hydro> floor_audit cycle 6959 [refresh] cell-c2p 22020096 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) fl|71:<hydro> floor_audit cycle 6960 [update] cell-c2p 44040192 e_c<=0 34 kedom 5
- PDmB: 2026-11-24 14:57:37 rc=0 wall=351 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 3.143855e+02 | faile`; keys in the log: 52:<hydro> floor_audit cycle 6959 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de 0|53:<hydro> floor_audit cycle 6959 [refresh] cell-c2p 22020096 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) fl|71:<hydro> floor_audit cycle 6960 [update] cell-c2p 44040192 e_c<=0 34 kedom 5

## PSm: UA vs PSmB (11-24 15:02)
| arm | rows | s/cyc | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries/skipped | clamped reads line | L_top/L_in end | L_bot/L_in end | KE end | dM/M | dE/E | dLz/Lz |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PSmA | 100 | 1.0782 (c7000-7050) | 12.01 / 12.01 / 12.01 | 37.4, 5.19 | 17.6, 3.98 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 0.98785 | 1.00015 | 1.57604e+39 | 3.21e-15 | 4.68e-07 | -1.55e-07 |
| PSmB | 100 | 1.2015 (c7000-7050) | 12.01 / 12.01 / 12.01 | 39.8, 5.21 | 16.7, 4.03 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 0.98791 | 1.00015 | 1.57604e+39 | 5.67e-16 | 4.68e-07 | 3.22e-07 |

- matched time (100 rows): L_top rel end 6.609e-05 max 1.741e-04; L_bot end -1.230e-08 max 5.768e-08; KE_int end 1.039e-08 max 2.434e-08
- dt mean ratio-1 1.232e-10; s per simulated s: A 0.08979 B 0.10006 ratio 1.114
- **NOT HARMLESS by the bounds: cost ratio 1.114 > 1.05** (bounds of 24zzzzl)
- PSmA: 2026-11-24 15:00:26 rc=0 wall=140 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 1.073944e+02 | faile`; keys in the log: 52:<hydro> floor_audit cycle 6959 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de 0|53:<hydro> floor_audit cycle 6959 [refresh] cell-c2p 22020096 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) fl|71:<hydro> floor_audit cycle 6960 [update] cell-c2p 44040192 e_c<=0 34 kedom 5
- PSmB: 2026-11-24 15:02:56 rc=0 wall=150 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 1.203290e+02 | faile`; keys in the log: 52:<hydro> floor_audit cycle 6959 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de 0|53:<hydro> floor_audit cycle 6959 [refresh] cell-c2p 22020096 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) fl|71:<hydro> floor_audit cycle 6960 [update] cell-c2p 44040192 e_c<=0 34 kedom 5

## PDt: UA vs PDtB (11-24 15:16)
| arm | rows | s/cyc | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries/skipped | clamped reads line | L_top/L_in end | L_bot/L_in end | KE end | dM/M | dE/E | dLz/Lz |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PDtA | 300 | 3.0297 (c5300-5550) | 13.15 / 34.13 / 46.87 | 14.5, 4.80 | 8.5, 3.53 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 2.43645 | 0.98533 | 1.78481e+44 | -1.55e-16 | -6.23e-05 | -5.14e-06 |
| PDtB | 300 | 3.0717 (c5300-5550) | 13.94 / 36.7 / 49.37 | 14.6, 4.76 | 8.4, 3.37 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 2.42708 | 0.98664 | 1.78719e+44 | -9.27e-16 | -6.68e-05 | -5.30e-06 |

- matched time (300 rows): L_top rel end -2.974e-05 max 2.050e-04; L_bot end 1.541e-06 max 3.059e-06; KE_int end -8.994e-07 max 1.353e-06
- dt mean ratio-1 7.524e-02; s per simulated s: A 0.08876 B 0.08369 ratio 0.943
- **NOT HARMLESS by the bounds: dt mean 7.52 % off** (bounds of 24zzzzl)
- PDtA: 2026-11-24 14:59:30 rc=0 wall=971 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 9.008839e+02 | faile`; keys in the log: 28:### WARNING (scheme = clean): no effect: time/dt_retry_max rad_vet/clean_stag_thinf_held rad_vet/time2_start rad_vet/vet_col_lat_taucut rad_vet/time2_stage_l|58:<hydro> floor_audit cycle 5279 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de 0|59:<hydro> floor_audit cycle 5279 [refresh] cell-c2p 67108864 e_c<=0 0 kedom 8
- PDtB: 2026-11-24 15:16:27 rc=0 wall=991 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 9.153046e+02 | faile`; keys in the log: 28:### WARNING (scheme = clean): no effect: time/dt_retry_max rad_vet/clean_stag_thinf_held rad_vet/time2_start rad_vet/vet_col_lat_taucut rad_vet/time2_stage_l|58:<hydro> floor_audit cycle 5279 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de 0|59:<hydro> floor_audit cycle 5279 [refresh] cell-c2p 67108864 e_c<=0 0 kedom 8

## PSt: UA vs PStB (11-24 15:29)
| arm | rows | s/cyc | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries/skipped | clamped reads line | L_top/L_in end | L_bot/L_in end | KE end | dM/M | dE/E | dLz/Lz |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PStA | 100 | 3.0685 (c5300-5350) | 13.15 / 27.02 / 46.39 | 14.7, 5.00 | 8.7, 3.81 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 2.52636 | 1.01453 | 1.75213e+44 | 9.27e-16 | -1.68e-05 | -1.76e-06 |
| PStB | 100 | 3.1325 (c5300-5350) | 12.51 / 25.32 / 41.96 | 13.1, 5.03 | 8.2, 4.02 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0.000000e+00 | 2.52803 | 1.01427 | 1.75120e+44 | 1.24e-15 | -1.58e-05 | -1.78e-06 |

- matched time (94 rows): L_top rel end -1.007e-04 max 1.007e-04; L_bot end -8.171e-07 max 8.171e-07; KE_int end 3.341e-07 max 6.661e-07
- dt mean ratio-1 -6.297e-02; s per simulated s: A 0.11357 B 0.12373 ratio 1.089
- **NOT HARMLESS by the bounds: dt mean -6.30 % off; cost ratio 1.089 > 1.05** (bounds of 24zzzzl)
- PStA: 2026-11-24 15:23:00 rc=0 wall=373 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 3.079916e+02 | faile`; keys in the log: 28:### WARNING (scheme = clean): no effect: time/dt_retry_max rad_vet/clean_stag_thinf_held rad_vet/time2_start rad_vet/vet_col_lat_taucut rad_vet/time2_stage_l|58:<hydro> floor_audit cycle 5279 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de 0|59:<hydro> floor_audit cycle 5279 [refresh] cell-c2p 67108864 e_c<=0 0 kedom 8
- PStB: 2026-11-24 15:29:28 rc=0 wall=388 s CW 0 FATAL 0 inf 0 DT_RETRY 0; last: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 3.157907e+02 | faile`; keys in the log: 28:### WARNING (scheme = clean): no effect: time/dt_retry_max rad_vet/clean_stag_thinf_held rad_vet/time2_start rad_vet/vet_col_lat_taucut rad_vet/time2_stage_l|58:<hydro> floor_audit cycle 5279 [update] cell-c2p 0 e_c<=0 0 kedom 0 sub_tfloor 0 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 0 (kedom 0) floor_de 0|59:<hydro> floor_audit cycle 5279 [refresh] cell-c2p 67108864 e_c<=0 0 kedom 8
