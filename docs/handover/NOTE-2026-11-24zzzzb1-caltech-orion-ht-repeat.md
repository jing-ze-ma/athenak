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
