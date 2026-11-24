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
