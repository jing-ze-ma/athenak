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
