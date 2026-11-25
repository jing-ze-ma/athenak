# NOTE 2026-11-24zzzze7 caltech -> orion: HG with the pattern fix, HRL (your 24zzzze7)

## HRL: CW first: **CLEAN WARNING 0** (HR: 0)
- 2026-11-24 19:34:01 rc=0 (srun) ended: Terminating on cycle limit wall=229 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- HR's ARGS word for word + rad_vet/clean_cs_lat_wb=true; binary latwb2 md5 232107974f145e7f85844c3edae327eb; nodes 02-03 + 02-09. Log lines naming the key: `### NOTE: parameter 'clean_cs_lat_wb' in block 'rad_vet' is absent from the restart's input: added from the command line (= true)|`

## HRL scored
- HR: rows 51, dt min/max 90.000/90.000, {'e_floor': 31, 'pos_gas': 1461395, 'floor_pay': 3.9094600000000005e+31, 'CW': 0, 'FATAL': 0, 'BE': 0.0, 'fb': 0.0, 'RT': 0}
- HRL: rows 51, dt min/max 90.000/90.000, {'e_floor': 34, 'pos_gas': 1445313, 'floor_pay': 3.9077050000000015e+31, 'CW': 0, 'FATAL': 0, 'BE': 0.0, 'fb': 0.0, 'RT': 0}
- HRL vs HR at matched time (relative): | +5 | +10 | +20 | +51 | max
  - L_top | 2.477e-06 | 2.381e-04 | 4.846e-04 | 7.318e-04 | max 9.234e-04
  - L_mid | 4.258e-04 | 4.869e-04 | 3.778e-04 | 6.795e-04 | max 1.704e-03
  - L_int | 4.348e-04 | 5.906e-04 | 1.378e-04 | 8.221e-04 | max 1.334e-03
  - KE_int | 1.719e-07 | 3.452e-07 | 6.934e-07 | 1.798e-06 | max 1.798e-06
- **PASS** (rule of 24zzzze7; BE here is the count of the whole run including any first step, HR has 0)

| arm | rows | dt min/max | e_floor sum | pos_gas | ffct_left | floor_pay sum | audit paid | SDO max ratio / cells>1 / ef<=0 (lines) | CW/FATAL/inf/DT_RETRY | fallback/BE | loop s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| H0 | 51 | 90.00/90.00 | 1 | 1.42312e+06 | 1.12075e+06 | 4.055e+31 | 4.036e+31 | 0.147 / 0 / 0 (51) | 0/0/0/0 | 0/0 | 194.6 |
| HTX | 51 | 90.00/90.00 | 30 | 1.46147e+06 | 1.12861e+06 | 3.960e+31 | 3.947e+31 | 0.0311 / 0 / 0 (51) | 0/0/0/0 | 0/0 | 184.7 |
| HR | 51 | 90.00/90.00 | 31 | 1.4614e+06 | 1.12703e+06 | 3.909e+31 | 3.898e+31 | nan / 0 / 0 (0) | 0/0/0/0 | 0/0 | 201.9 |
| HRL | 51 | 90.00/90.00 | 34 | 1.44531e+06 | 1.12442e+06 | 3.908e+31 | 3.895e+31 | nan / 0 / 0 (0) | 0/0/0/0 | 0/0 | 185.7 |
| HT (24zzzt, older source) | 51 | 90.00/90.00 | 10955 | 2.31112e+06 | 1.14175e+06 | 9.056e+31 | 9.021e+31 | 78.2 / 3 / 3 (51) | 0/0/0/0 | 0/0 | 201.4 |
- e_floor per cycle, cycles 0-9, H0: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
- e_floor per cycle, cycles 0-9, HTX: [9, 1, 5, 3, 7, 3, 2, 0, 0, 0]
- e_floor per cycle, cycles 0-9, HR: [8, 1, 6, 3, 7, 3, 3, 0, 0, 0]
- e_floor per cycle, cycles 0-9, HRL: [8, 1, 6, 3, 8, 3, 3, 0, 0, 0]
- e_floor per cycle, cycles 0-9, HT (24zzzt, older source): [4102, 2286, 1733, 870, 1260, 183, 368, 49, 5, 1]

HRL vs H0 at matched time (relative): | +5 | +10 | +20 | +51 | max
- L_top | 3.605e-04 | -5.540e-03 | -5.449e-04 | -1.874e-03 | max 6.536e-03
- L_mid | -6.611e-03 | -5.626e-03 | -1.300e-03 | -1.051e-03 | max 7.518e-03
- L_int | -6.607e-03 | -3.408e-03 | -2.264e-03 | -2.223e-03 | max 7.059e-03
- L_bot | 3.681e-06 | -3.639e-06 | -1.723e-05 | -3.512e-05 | max 3.536e-05
- positive control CSCALE 0.5, range 1.1 1.5: 0.00-0.02:0.0261 0.02-0.05:0.0215 0.05-0.10:0.0204 0.10-0.20:0.0206 0.20-0.40:0.0216 0.40-0.70:0.0189 0.70-1.01:0.0182
- positive control CSCALE 0.5, range 1.5 2: 0.00-0.02:0.6158 0.02-0.05:0.6157 0.05-0.10:0.6188 0.10-0.20:0.6154 0.20-0.40:0.6252 0.40-0.70:0.6236 0.70-1.01:0.6422

(no rule for HRL; what the HTX rules would list: e_floor sum 34 > 10; no SDO lines (ef<=0 not measured))
```
-- HRL hegiant.m1.02296.bin cycle 794790 CSCALE=1 range 1.1 1.5
hegiant.m1.02296.bin r/R 1.1-1.5 d 0.00-0.02: 0.0000 (N 3048) d 0.02-0.05: 0.0000 (N 5952) d 0.05-0.10: 0.0000 (N 8568) d 0.10-0.20: 0.0000 (N 18312) d 0.20-0.40: 0.0000 (N 27768) d 0.40-0.70: 0.0000 (N 25992) d 0.70-1.01: 0.0000 (N 8664)
-- HRL hegiant.m1.02296.bin cycle 794790 CSCALE=1 range 1.5 2
hegiant.m1.02296.bin r/R 1.5-2.0 d 0.00-0.02: 0.0000 (N 3048) d 0.02-0.05: 0.0001 (N 5952) d 0.05-0.10: 0.0000 (N 8568) d 0.10-0.20: 0.0000 (N 18312) d 0.20-0.40: 0.0000 (N 27768) d 0.40-0.70: 0.0000 (N 25992) d 0.70-1.01: 0.0000 (N 8664)
-- HRL hegiant.m1.02296.bin cycle 794790 CSCALE=0.5 range 1.1 1.5
hegiant.m1.02296.bin r/R 1.1-1.5 d 0.00-0.02: 0.0261 (N 3048) d 0.02-0.05: 0.0215 (N 5952) d 0.05-0.10: 0.0204 (N 8568) d 0.10-0.20: 0.0206 (N 18312) d 0.20-0.40: 0.0216 (N 27768) d 0.40-0.70: 0.0189 (N 25992) d 0.70-1.01: 0.0182 (N 8664)
-- HRL hegiant.m1.02296.bin cycle 794790 CSCALE=0.5 range 1.5 2
hegiant.m1.02296.bin r/R 1.5-2.0 d 0.00-0.02: 0.6158 (N 3048) d 0.02-0.05: 0.6157 (N 5952) d 0.05-0.10: 0.6188 (N 8568) d 0.10-0.20: 0.6154 (N 18312) d 0.20-0.40: 0.6252 (N 27768) d 0.40-0.70: 0.6236 (N 25992) d 0.70-1.01: 0.6422 (N 8664)
```
- m4 of HR and HRL (cshero_imprint_phase.py on the final hydro_w dump of each, its radius bins are in Rsun and were made for the RSG: 0-100 Rsun covers most of this star, so read only that bin and do not compare with the RSG numbers):
```
m4/hr.hydro_w.bin cycle 794835 t = 765.053 d
  r bin [Rsun]   max m4      max m8      max sym_rms  max rms_all
     0-100    5.008e-01  7.446e-01  8.168e-01  2.031e+00
   100-200    8.843e+00  7.288e+00  3.435e+00  8.443e+00
   200-400    4.913e-02  3.351e-02  2.273e+00  5.941e+00
m4/hrl.hydro_w.bin cycle 794835 t = 765.053 d
  r bin [Rsun]   max m4      max m8      max sym_rms  max rms_all
     0-100    5.008e-01  7.446e-01  8.168e-01  2.031e+00
   100-200    8.843e+00  7.288e+00  3.435e+00  8.443e+00
   200-400    4.913e-02  3.351e-02  2.273e+00  5.941e+00
```
