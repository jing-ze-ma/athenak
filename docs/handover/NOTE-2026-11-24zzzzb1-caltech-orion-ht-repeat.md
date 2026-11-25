# NOTE 2026-11-24zzzzb1 caltech -> orion: HG HT repeat with the cs thinf fix (your 24zzzzb1)

- **H0 and HTX cannot run at the same time:** 4 ranks do not fit, and only holder 7333072 (3 nodes) has free nodes, so there is one 2-node pair; the arms run one after the other on 02-03 + 02-09. Production is untouched.
- Args = HT's (24zzzt) word for word plus clean_stag_thinf named in every arm; time2_gtime_slope is not named (as in HT: absent = off). SDO_PRINT=1 in every arm.
- BUILD_OK sdosdo md5=573465bddd359381dd270fd0428e8683 prov=53290b5af83bfd487f5ed12d73532c9407501257
- GATE PASS sdosdo: athena_sdosdo_53290b5af md5 573465bd node hpc-sm-01-15 pass1=3.686e-11 log date 11-24 17:16:37 end: Terminating on cycle limit
- BUILD_OK thinfcs md5=fb0549d6630c36c5e7405ef6e529fec3 prov=87af1d65c6b72a88664ee009042291c7854f26f5
- GATE PASS thinfcs: athena_thinfcs_87af1d65c md5 fb0549d6 node hpc-sm-01-15 pass1=3.686e-11 log date 11-24 17:16:52 end: Terminating on cycle limit

## +1 step check (nlim 794785)
- P1base (ok): 2026-11-24 17:17:37 rc=0 (srun) ended: Terminating on cycle limit wall=43 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps ; peak 69545 MiB per GPU
- P1fix (ok): 2026-11-24 17:18:20 rc=0 (srun) ended: Terminating on cycle limit wall=43 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps ; peak 69913 MiB per GPU

## H0 (11-24 17:22)
- 2026-11-24 17:22:29 rc=0 (srun) ended: Terminating on cycle limit wall=231 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps ; peak GPU memory 75935 MiB; stage time clean stage steps 1.847012e+02; bins kept: 11 files in /resnick/groups/carnegie_poc/jingze/b1_1124/hg/H0/bin
| arm | rows | dt min/max | e_floor sum | pos_gas | ffct_left | floor_pay sum | audit paid | SDO max ratio / cells>1 / ef<=0 (lines) | CW/FATAL/inf/DT_RETRY | fallback/BE | loop s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| H0 | 51 | 90.00/90.00 | 1 | 1.42312e+06 | 1.12075e+06 | 4.055e+31 | 4.036e+31 | 0.147 / 0 / 0 (51) | 0/0/0/0 | 0/0 | 194.6 |
| HT (24zzzt, older source) | 51 | 90.00/90.00 | 10955 | 2.31112e+06 | 1.14175e+06 | 9.056e+31 | 9.021e+31 | 78.2 / 3 / 3 (51) | 0/0/0/0 | 0/0 | 201.4 |
- e_floor per cycle, cycles 0-9, H0: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
- e_floor per cycle, cycles 0-9, HT (24zzzt, older source): [4102, 2286, 1733, 870, 1260, 183, 368, 49, 5, 1]
```
-- H0 hegiant.m1.02296.bin cycle 794790 CSCALE=1 range 1.1 1.5
hegiant.m1.02296.bin r/R 1.1-1.5 d 0.00-0.02: 0.0000 (N 3048) d 0.02-0.05: 0.0000 (N 5952) d 0.05-0.10: 0.0000 (N 8568) d 0.10-0.20: 0.0000 (N 18312) d 0.20-0.40: 0.0000 (N 27768) d 0.40-0.70: 0.0000 (N 25992) d 0.70-1.01: 0.0000 (N 8664)
-- H0 hegiant.m1.02296.bin cycle 794790 CSCALE=1 range 1.5 2
hegiant.m1.02296.bin r/R 1.5-2.0 d 0.00-0.02: 0.0011 (N 3048) d 0.02-0.05: 0.0000 (N 5952) d 0.05-0.10: 0.0000 (N 8568) d 0.10-0.20: 0.0000 (N 18312) d 0.20-0.40: 0.0000 (N 27768) d 0.40-0.70: 0.0000 (N 25992) d 0.70-1.01: 0.0000 (N 8664)
-- H0 hegiant.m1.02296.bin cycle 794790 CSCALE=0.5 range 1.1 1.5
hegiant.m1.02296.bin r/R 1.1-1.5 d 0.00-0.02: 0.0255 (N 3048) d 0.02-0.05: 0.0212 (N 5952) d 0.05-0.10: 0.0200 (N 8568) d 0.10-0.20: 0.0202 (N 18312) d 0.20-0.40: 0.0212 (N 27768) d 0.40-0.70: 0.0189 (N 25992) d 0.70-1.01: 0.0177 (N 8664)
-- H0 hegiant.m1.02296.bin cycle 794790 CSCALE=0.5 range 1.5 2
hegiant.m1.02296.bin r/R 1.5-2.0 d 0.00-0.02: 0.6139 (N 3048) d 0.02-0.05: 0.6156 (N 5952) d 0.05-0.10: 0.6188 (N 8568) d 0.10-0.20: 0.6146 (N 18312) d 0.20-0.40: 0.6245 (N 27768) d 0.40-0.70: 0.6231 (N 25992) d 0.70-1.01: 0.6414 (N 8664)
```

## HTX (11-24 17:26)
- 2026-11-24 17:26:15 rc=0 (srun) ended: Terminating on cycle limit wall=226 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps ; peak GPU memory 75879 MiB; stage time clean stage steps 1.747711e+02; bins kept: 11 files in /resnick/groups/carnegie_poc/jingze/b1_1124/hg/HTX/bin
| arm | rows | dt min/max | e_floor sum | pos_gas | ffct_left | floor_pay sum | audit paid | SDO max ratio / cells>1 / ef<=0 (lines) | CW/FATAL/inf/DT_RETRY | fallback/BE | loop s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| H0 | 51 | 90.00/90.00 | 1 | 1.42312e+06 | 1.12075e+06 | 4.055e+31 | 4.036e+31 | 0.147 / 0 / 0 (51) | 0/0/0/0 | 0/0 | 194.6 |
| HTX | 51 | 90.00/90.00 | 30 | 1.46147e+06 | 1.12861e+06 | 3.960e+31 | 3.947e+31 | 0.0311 / 0 / 0 (51) | 0/0/0/0 | 0/0 | 184.7 |
| HT (24zzzt, older source) | 51 | 90.00/90.00 | 10955 | 2.31112e+06 | 1.14175e+06 | 9.056e+31 | 9.021e+31 | 78.2 / 3 / 3 (51) | 0/0/0/0 | 0/0 | 201.4 |
- e_floor per cycle, cycles 0-9, H0: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
- e_floor per cycle, cycles 0-9, HTX: [9, 1, 5, 3, 7, 3, 2, 0, 0, 0]
- e_floor per cycle, cycles 0-9, HT (24zzzt, older source): [4102, 2286, 1733, 870, 1260, 183, 368, 49, 5, 1]

HTX vs H0 at matched time (relative): | +5 | +10 | +20 | +51 | max
- L_top | 3.593e-04 | -5.752e-03 | -8.685e-04 | -2.119e-03 | max 6.915e-03
- L_mid | -7.056e-03 | -6.184e-03 | -1.470e-03 | -1.967e-03 | max 7.796e-03
- L_int | -7.045e-03 | -4.021e-03 | -2.121e-03 | -2.316e-03 | max 7.437e-03
- L_bot | -1.234e-10 | -3.535e-10 | -2.289e-11 | -1.387e-10 | max 8.633e-10
- positive control CSCALE 0.5, range 1.1 1.5: 0.00-0.02:0.0260 0.02-0.05:0.0215 0.05-0.10:0.0204 0.10-0.20:0.0207 0.20-0.40:0.0216 0.40-0.70:0.0188 0.70-1.01:0.0182
- positive control CSCALE 0.5, range 1.5 2: 0.00-0.02:0.6160 0.02-0.05:0.6159 0.05-0.10:0.6188 0.10-0.20:0.6154 0.20-0.40:0.6252 0.40-0.70:0.6236 0.70-1.01:0.6422

**NOT PASS: e_floor sum 30 > 10** (rules of 24zzzzb1; seam values are compared as printed, 4 decimals)
```
-- HTX hegiant.m1.02296.bin cycle 794790 CSCALE=1 range 1.1 1.5
hegiant.m1.02296.bin r/R 1.1-1.5 d 0.00-0.02: 0.0000 (N 3048) d 0.02-0.05: 0.0000 (N 5952) d 0.05-0.10: 0.0000 (N 8568) d 0.10-0.20: 0.0000 (N 18312) d 0.20-0.40: 0.0000 (N 27768) d 0.40-0.70: 0.0000 (N 25992) d 0.70-1.01: 0.0000 (N 8664)
-- HTX hegiant.m1.02296.bin cycle 794790 CSCALE=1 range 1.5 2
hegiant.m1.02296.bin r/R 1.5-2.0 d 0.00-0.02: 0.0000 (N 3048) d 0.02-0.05: 0.0001 (N 5952) d 0.05-0.10: 0.0000 (N 8568) d 0.10-0.20: 0.0000 (N 18312) d 0.20-0.40: 0.0000 (N 27768) d 0.40-0.70: 0.0000 (N 25992) d 0.70-1.01: 0.0000 (N 8664)
-- HTX hegiant.m1.02296.bin cycle 794790 CSCALE=0.5 range 1.1 1.5
hegiant.m1.02296.bin r/R 1.1-1.5 d 0.00-0.02: 0.0260 (N 3048) d 0.02-0.05: 0.0215 (N 5952) d 0.05-0.10: 0.0204 (N 8568) d 0.10-0.20: 0.0207 (N 18312) d 0.20-0.40: 0.0216 (N 27768) d 0.40-0.70: 0.0188 (N 25992) d 0.70-1.01: 0.0182 (N 8664)
-- HTX hegiant.m1.02296.bin cycle 794790 CSCALE=0.5 range 1.5 2
hegiant.m1.02296.bin r/R 1.5-2.0 d 0.00-0.02: 0.6160 (N 3048) d 0.02-0.05: 0.6159 (N 5952) d 0.05-0.10: 0.6188 (N 8568) d 0.10-0.20: 0.6154 (N 18312) d 0.20-0.40: 0.6252 (N 27768) d 0.40-0.70: 0.6236 (N 25992) d 0.70-1.01: 0.6422 (N 8664)
```

## HX (11-24 17:30)
- 2026-11-24 17:30:00 rc=0 (srun) ended: Terminating on cycle limit wall=225 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps ; peak GPU memory 76011 MiB; stage time clean stage steps 1.734811e+02; bins kept: 11 files in /resnick/groups/carnegie_poc/jingze/b1_1124/hg/HX/bin
| arm | rows | dt min/max | e_floor sum | pos_gas | ffct_left | floor_pay sum | audit paid | SDO max ratio / cells>1 / ef<=0 (lines) | CW/FATAL/inf/DT_RETRY | fallback/BE | loop s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| H0 | 51 | 90.00/90.00 | 1 | 1.42312e+06 | 1.12075e+06 | 4.055e+31 | 4.036e+31 | 0.147 / 0 / 0 (51) | 0/0/0/0 | 0/0 | 194.6 |
| HX | 51 | 90.00/90.00 | 1 | 1.42312e+06 | 1.12075e+06 | 4.055e+31 | 4.036e+31 | 0.147 / 0 / 0 (51) | 0/0/0/0 | 0/0 | 183.6 |
| HT (24zzzt, older source) | 51 | 90.00/90.00 | 10955 | 2.31112e+06 | 1.14175e+06 | 9.056e+31 | 9.021e+31 | 78.2 / 3 / 3 (51) | 0/0/0/0 | 0/0 | 201.4 |
- e_floor per cycle, cycles 0-9, H0: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
- e_floor per cycle, cycles 0-9, HX: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
- e_floor per cycle, cycles 0-9, HT (24zzzt, older source): [4102, 2286, 1733, 870, 1260, 183, 368, 49, 5, 1]

HX vs H0 at matched time (relative): | +5 | +10 | +20 | +51 | max
- L_top | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | max 0.000e+00
- L_mid | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | max 0.000e+00
- L_int | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | max 0.000e+00
- L_bot | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | max 0.000e+00
- positive control CSCALE 0.5, range 1.1 1.5: 0.00-0.02:0.0255 0.02-0.05:0.0212 0.05-0.10:0.0200 0.10-0.20:0.0202 0.20-0.40:0.0212 0.40-0.70:0.0189 0.70-1.01:0.0177
- positive control CSCALE 0.5, range 1.5 2: 0.00-0.02:0.6139 0.02-0.05:0.6156 0.05-0.10:0.6188 0.10-0.20:0.6146 0.20-0.40:0.6245 0.40-0.70:0.6231 0.70-1.01:0.6414

(no rule for HX; what the HTX rules would list: nothing)
```
-- HX hegiant.m1.02296.bin cycle 794790 CSCALE=1 range 1.1 1.5
hegiant.m1.02296.bin r/R 1.1-1.5 d 0.00-0.02: 0.0000 (N 3048) d 0.02-0.05: 0.0000 (N 5952) d 0.05-0.10: 0.0000 (N 8568) d 0.10-0.20: 0.0000 (N 18312) d 0.20-0.40: 0.0000 (N 27768) d 0.40-0.70: 0.0000 (N 25992) d 0.70-1.01: 0.0000 (N 8664)
-- HX hegiant.m1.02296.bin cycle 794790 CSCALE=1 range 1.5 2
hegiant.m1.02296.bin r/R 1.5-2.0 d 0.00-0.02: 0.0011 (N 3048) d 0.02-0.05: 0.0000 (N 5952) d 0.05-0.10: 0.0000 (N 8568) d 0.10-0.20: 0.0000 (N 18312) d 0.20-0.40: 0.0000 (N 27768) d 0.40-0.70: 0.0000 (N 25992) d 0.70-1.01: 0.0000 (N 8664)
-- HX hegiant.m1.02296.bin cycle 794790 CSCALE=0.5 range 1.1 1.5
hegiant.m1.02296.bin r/R 1.1-1.5 d 0.00-0.02: 0.0255 (N 3048) d 0.02-0.05: 0.0212 (N 5952) d 0.05-0.10: 0.0200 (N 8568) d 0.10-0.20: 0.0202 (N 18312) d 0.20-0.40: 0.0212 (N 27768) d 0.40-0.70: 0.0189 (N 25992) d 0.70-1.01: 0.0177 (N 8664)
-- HX hegiant.m1.02296.bin cycle 794790 CSCALE=0.5 range 1.5 2
hegiant.m1.02296.bin r/R 1.5-2.0 d 0.00-0.02: 0.6139 (N 3048) d 0.02-0.05: 0.6156 (N 5952) d 0.05-0.10: 0.6188 (N 8568) d 0.10-0.20: 0.6146 (N 18312) d 0.20-0.40: 0.6245 (N 27768) d 0.40-0.70: 0.6231 (N 25992) d 0.70-1.01: 0.6414 (N 8664)
```

## HT2 (11-24 17:34)
- 2026-11-24 17:33:59 rc=0 (srun) ended: Terminating on cycle limit wall=239 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps ; peak GPU memory 76007 MiB; stage time clean stage steps 1.871276e+02; bins kept: 11 files in /resnick/groups/carnegie_poc/jingze/b1_1124/hg/HT2/bin
| arm | rows | dt min/max | e_floor sum | pos_gas | ffct_left | floor_pay sum | audit paid | SDO max ratio / cells>1 / ef<=0 (lines) | CW/FATAL/inf/DT_RETRY | fallback/BE | loop s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| H0 | 51 | 90.00/90.00 | 1 | 1.42312e+06 | 1.12075e+06 | 4.055e+31 | 4.036e+31 | 0.147 / 0 / 0 (51) | 0/0/0/0 | 0/0 | 194.6 |
| HT2 | 51 | 90.00/90.00 | 10955 | 2.31112e+06 | 1.14175e+06 | 9.056e+31 | 9.021e+31 | 0.25 / 0 / 0 (51) | 0/0/0/0 | 0/0 | 197.3 |
| HT (24zzzt, older source) | 51 | 90.00/90.00 | 10955 | 2.31112e+06 | 1.14175e+06 | 9.056e+31 | 9.021e+31 | 78.2 / 3 / 3 (51) | 0/0/0/0 | 0/0 | 201.4 |
- e_floor per cycle, cycles 0-9, H0: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
- e_floor per cycle, cycles 0-9, HT2: [4102, 2286, 1733, 870, 1260, 183, 368, 49, 5, 1]
- e_floor per cycle, cycles 0-9, HT (24zzzt, older source): [4102, 2286, 1733, 870, 1260, 183, 368, 49, 5, 1]

HT2 vs H0 at matched time (relative): | +5 | +10 | +20 | +51 | max
- L_top | -1.035e-04 | -5.297e-03 | 2.175e-03 | 1.700e-03 | max 6.186e-03
- L_mid | -4.563e-03 | -6.328e-03 | 3.135e-03 | -1.420e-03 | max 6.433e-03
- L_int | -5.102e-03 | -3.666e-03 | 3.153e-03 | -1.370e-03 | max 6.532e-03
- L_bot | -3.882e-10 | -5.298e-11 | -9.217e-10 | -7.234e-10 | max 9.217e-10
- positive control CSCALE 0.5, range 1.1 1.5: 0.00-0.02:0.0272 0.02-0.05:0.0222 0.05-0.10:0.0211 0.10-0.20:0.0210 0.20-0.40:0.0225 0.40-0.70:0.0190 0.70-1.01:0.0177
- positive control CSCALE 0.5, range 1.5 2: 0.00-0.02:0.6107 0.02-0.05:0.6124 0.05-0.10:0.6132 0.10-0.20:0.6110 0.20-0.40:0.6241 0.40-0.70:0.6226 0.70-1.01:0.6407

(no rule for HT2; what the HTX rules would list: seam 1.1 1.5 R d 0.00-0.02: 0.0027 > 0.0000 (H0 0.0000); seam 1.1 1.5 R d 0.02-0.05: 0.0012 > 0.0002 (H0 0.0000); seam 1.1 1.5 R d 0.05-0.10: 0.0005 > 0.0002 (H0 0.0000); seam 1.1 1.5 R d 0.70-1.01: 0.0003 > 0.0002 (H0 0.0000); seam 1.5 2 R d 0.00-0.02: 0.0079 > 0.0011 (H0 0.0011); seam 1.5 2 R d 0.02-0.05: 0.0025 > 0.0002 (H0 0.0000); seam 1.5 2 R d 0.05-0.10: 0.0003 > 0.0002 (H0 0.0000); seam 1.5 2 R d 0.70-1.01: 0.0008 > 0.0002 (H0 0.0000); e_floor sum 10955 > 10; pos_gas 2.31112e+06 > 1.2 x 1.42312e+06; floor_pay 9.056e+31 > 1.2 x 4.055e+31)
```
-- HT2 hegiant.m1.02296.bin cycle 794790 CSCALE=1 range 1.1 1.5
hegiant.m1.02296.bin r/R 1.1-1.5 d 0.00-0.02: 0.0027 (N 3048) d 0.02-0.05: 0.0012 (N 5952) d 0.05-0.10: 0.0005 (N 8568) d 0.10-0.20: 0.0002 (N 18312) d 0.20-0.40: 0.0002 (N 27768) d 0.40-0.70: 0.0001 (N 25992) d 0.70-1.01: 0.0003 (N 8664)
-- HT2 hegiant.m1.02296.bin cycle 794790 CSCALE=1 range 1.5 2
hegiant.m1.02296.bin r/R 1.5-2.0 d 0.00-0.02: 0.0079 (N 3048) d 0.02-0.05: 0.0025 (N 5952) d 0.05-0.10: 0.0003 (N 8568) d 0.10-0.20: 0.0002 (N 18312) d 0.20-0.40: 0.0001 (N 27768) d 0.40-0.70: 0.0002 (N 25992) d 0.70-1.01: 0.0008 (N 8664)
-- HT2 hegiant.m1.02296.bin cycle 794790 CSCALE=0.5 range 1.1 1.5
hegiant.m1.02296.bin r/R 1.1-1.5 d 0.00-0.02: 0.0272 (N 3048) d 0.02-0.05: 0.0222 (N 5952) d 0.05-0.10: 0.0211 (N 8568) d 0.10-0.20: 0.0210 (N 18312) d 0.20-0.40: 0.0225 (N 27768) d 0.40-0.70: 0.0190 (N 25992) d 0.70-1.01: 0.0177 (N 8664)
-- HT2 hegiant.m1.02296.bin cycle 794790 CSCALE=0.5 range 1.5 2
hegiant.m1.02296.bin r/R 1.5-2.0 d 0.00-0.02: 0.6107 (N 3048) d 0.02-0.05: 0.6124 (N 5952) d 0.05-0.10: 0.6132 (N 8568) d 0.10-0.20: 0.6110 (N 18312) d 0.20-0.40: 0.6241 (N 27768) d 0.40-0.70: 0.6226 (N 25992) d 0.70-1.01: 0.6407 (N 8664)
```
