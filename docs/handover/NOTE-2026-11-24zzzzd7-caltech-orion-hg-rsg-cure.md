# NOTE 2026-11-24zzzzd7 caltech -> orion: HG and RSG restarts on rel-a c880e2703 (your 24zzzzd7)

## HR: CW first: **CLEAN WARNING 0** (HTX: 0)
- 2026-11-24 19:10:34 rc=0 (srun) ended: Terminating on cycle limit wall=238 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- HTX's ARGS word for word (from its run.log); binary rel-a he_star_vet md5 5cfc3be547ad477553ec5db0e9c07fca; nodes as HTX, 8 ranks. Builds (rel-b on 02-09, latwb2 on 02-03; CPU steps) ran at the same time and share only CPU cores: the loop time is not a timing.

## HR scored (rows vs H0; HTX row for reference; the script's rule list is the HTX rule of zzzzb1, kept as a listing)
- peak GPU memory 75887 MiB
| arm | rows | dt min/max | e_floor sum | pos_gas | ffct_left | floor_pay sum | audit paid | SDO max ratio / cells>1 / ef<=0 (lines) | CW/FATAL/inf/DT_RETRY | fallback/BE | loop s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| H0 | 51 | 90.00/90.00 | 1 | 1.42312e+06 | 1.12075e+06 | 4.055e+31 | 4.036e+31 | 0.147 / 0 / 0 (51) | 0/0/0/0 | 0/0 | 194.6 |
| HTX | 51 | 90.00/90.00 | 30 | 1.46147e+06 | 1.12861e+06 | 3.960e+31 | 3.947e+31 | 0.0311 / 0 / 0 (51) | 0/0/0/0 | 0/0 | 184.7 |
| HR | 51 | 90.00/90.00 | 31 | 1.4614e+06 | 1.12703e+06 | 3.909e+31 | 3.898e+31 | nan / 0 / 0 (0) | 0/0/0/0 | 0/0 | 201.9 |
| HT (24zzzt, older source) | 51 | 90.00/90.00 | 10955 | 2.31112e+06 | 1.14175e+06 | 9.056e+31 | 9.021e+31 | 78.2 / 3 / 3 (51) | 0/0/0/0 | 0/0 | 201.4 |
- e_floor per cycle, cycles 0-9, H0: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
- e_floor per cycle, cycles 0-9, HTX: [9, 1, 5, 3, 7, 3, 2, 0, 0, 0]
- e_floor per cycle, cycles 0-9, HR: [8, 1, 6, 3, 7, 3, 3, 0, 0, 0]
- e_floor per cycle, cycles 0-9, HT (24zzzt, older source): [4102, 2286, 1733, 870, 1260, 183, 368, 49, 5, 1]

HR vs H0 at matched time (relative): | +5 | +10 | +20 | +51 | max
- L_top | 3.580e-04 | -5.777e-03 | -1.029e-03 | -2.604e-03 | max 6.904e-03
- L_mid | -7.034e-03 | -6.110e-03 | -1.677e-03 | -1.729e-03 | max 7.842e-03
- L_int | -7.039e-03 | -3.996e-03 | -2.402e-03 | -3.043e-03 | max 7.441e-03
- L_bot | -1.385e-11 | 2.018e-10 | -6.311e-10 | -1.320e-10 | max 1.128e-09
- positive control CSCALE 0.5, range 1.1 1.5: 0.00-0.02:0.0259 0.02-0.05:0.0215 0.05-0.10:0.0204 0.10-0.20:0.0206 0.20-0.40:0.0216 0.40-0.70:0.0188 0.70-1.01:0.0182
- positive control CSCALE 0.5, range 1.5 2: 0.00-0.02:0.6158 0.02-0.05:0.6159 0.05-0.10:0.6189 0.10-0.20:0.6155 0.20-0.40:0.6252 0.40-0.70:0.6236 0.70-1.01:0.6422

(no rule for HR; what the HTX rules would list: e_floor sum 31 > 10; no SDO lines (ef<=0 not measured))
```
-- HR hegiant.m1.02296.bin cycle 794790 CSCALE=1 range 1.1 1.5
hegiant.m1.02296.bin r/R 1.1-1.5 d 0.00-0.02: 0.0000 (N 3048) d 0.02-0.05: 0.0000 (N 5952) d 0.05-0.10: 0.0000 (N 8568) d 0.10-0.20: 0.0000 (N 18312) d 0.20-0.40: 0.0000 (N 27768) d 0.40-0.70: 0.0000 (N 25992) d 0.70-1.01: 0.0000 (N 8664)
-- HR hegiant.m1.02296.bin cycle 794790 CSCALE=1 range 1.5 2
hegiant.m1.02296.bin r/R 1.5-2.0 d 0.00-0.02: 0.0000 (N 3048) d 0.02-0.05: 0.0001 (N 5952) d 0.05-0.10: 0.0000 (N 8568) d 0.10-0.20: 0.0000 (N 18312) d 0.20-0.40: 0.0000 (N 27768) d 0.40-0.70: 0.0000 (N 25992) d 0.70-1.01: 0.0000 (N 8664)
-- HR hegiant.m1.02296.bin cycle 794790 CSCALE=0.5 range 1.1 1.5
hegiant.m1.02296.bin r/R 1.1-1.5 d 0.00-0.02: 0.0259 (N 3048) d 0.02-0.05: 0.0215 (N 5952) d 0.05-0.10: 0.0204 (N 8568) d 0.10-0.20: 0.0206 (N 18312) d 0.20-0.40: 0.0216 (N 27768) d 0.40-0.70: 0.0188 (N 25992) d 0.70-1.01: 0.0182 (N 8664)
-- HR hegiant.m1.02296.bin cycle 794790 CSCALE=0.5 range 1.5 2
hegiant.m1.02296.bin r/R 1.5-2.0 d 0.00-0.02: 0.6158 (N 3048) d 0.02-0.05: 0.6159 (N 5952) d 0.05-0.10: 0.6189 (N 8568) d 0.10-0.20: 0.6155 (N 18312) d 0.20-0.40: 0.6252 (N 27768) d 0.40-0.70: 0.6236 (N 25992) d 0.70-1.01: 0.6422 (N 8664)
```
