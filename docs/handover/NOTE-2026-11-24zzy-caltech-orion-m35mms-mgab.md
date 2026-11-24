# NOTE 2026-11-24zzy caltech -> orion: m35mms mg 2 vs 4 A/B/A2 (your 24zzu) and production resume

- rst: m35mms.00005.rst = rst_keep/m35mms.pre-mgab.rst (stopped cleanly at cycle 6278); T2 binary 7d7c4a71, production keys, hst + log only, nothing else on 02-17 (one arm at a time).

| arm | mg | s/cyc c50-200 (mean) | min 50-cycle window | all 50-cycle windows | dt min/mean/max c50-200 | its per t2s2 (sum, solves, max sum) | its per t2s3 | Picard (hst col 22) mean | CW / FATAL / DT_RETRY / fallback |
|---|---|---|---|---|---|---|---|---|---|
| A | 2 | 1.3113 (c6350-6450) | 1.2825 | 1.356 1.340 1.283 | 9.078 / 10.050 / 11.300 | 66.7, 5.82, 83 | 32.1, 4.61, 54 | 0.00 | 0 / 0 / 0 / 0 |
| B | 4 | 1.1740 (c6350-6450) | 1.1529 | 1.198 1.195 1.153 | 9.078 / 10.050 / 11.300 | 46.5, 5.66, 60 | 22.3, 4.39, 40 | 0.00 | 0 / 0 / 0 / 0 |
| A2 | 2 | 1.3108 (c6350-6450) | 1.2818 | 1.361 1.340 1.282 | 9.078 / 10.050 / 11.300 | 66.7, 5.82, 83 | 32.1, 4.61, 54 | 0.00 | 0 / 0 / 0 / 0 |

Final hst row (cycle 6478): 
- rows A/B/A2: 200/200/200; A vs A2 bitwise (user.hst rows): True
- time: A 7.4226496114e+04  B 7.4226495666e+04  rel(B-A) -6.04e-09
- dt: A 9.5012028557e+00  B 9.5012022751e+00  rel(B-A) -6.11e-08
- L_top: A 1.0498471895e+37  B 1.0497838655e+37  rel(B-A) -6.03e-05
- Etot: A 3.0100333776e+44  B 3.0100333779e+44  rel(B-A) 9.31e-11
- M_tot: A 3.7255257204e+29  B 3.7255257204e+29  rel(B-A) -7.56e-16
- Lz_tot: A 5.1886342655e+36  B 5.1886318627e+36  rel(B-A) -4.63e-07

nvidia-smi during arms (5-s samples with util > 50 %): 
- A: 212 samples, SM clock min/mean 1980/1980 MHz, mem clock 3201, power mean/max 451/480 W, T max 49 C, memory 36387 MB
- B: 188 samples, SM clock min/mean 1980/1980 MHz, mem clock 3201, power mean/max 451/472 W, T max 49 C, memory 36325 MB
- A2: 211 samples, SM clock min/mean 1980/1980 MHz, mem clock 3201, power mean/max 454/480 W, T max 49 C, memory 36315 MB

Production m35mms history (300-cycle bins): cycles, s/cyc, its sum/solves t2s2, t2s3, dt
- log.0 c0-300: 0.645 s/cyc | t2s2 9.5/2.00 | t2s3 11.1/2.24 | dt 12.01
- log.0 c300-600: 0.627 s/cyc | t2s2 8.5/2.00 | t2s3 10.2/2.09 | dt 12.01
- log.0 c600-900: 0.650 s/cyc | t2s2 8.7/2.00 | t2s3 13.0/2.32 | dt 12.01
- log.0 c900-1200: 0.628 s/cyc | t2s2 8.9/2.00 | t2s3 9.6/2.16 | dt 12.01
- log.0 c1200-1500: 0.626 s/cyc | t2s2 8.7/2.00 | t2s3 10.1/2.07 | dt 12.01
- log.0 c1500-1800: 0.744 s/cyc | t2s2 8.7/2.00 | t2s3 10.2/2.16 | dt 12.01
- log.1 c2050-2350: 0.626 s/cyc | t2s2 10.1/2.00 | t2s3 9.1/2.04 | dt 12.01
- log.1 c2350-2650: 0.661 s/cyc | t2s2 9.7/2.00 | t2s3 12.9/2.53 | dt 12.01
- log.1 c2650-2950: 0.640 s/cyc | t2s2 11.4/2.00 | t2s3 9.9/2.14 | dt 12.01
- log.1 c2950-3250: 0.688 s/cyc | t2s2 12.0/2.00 | t2s3 15.8/2.57 | dt 12.01
- log.1 c3250-3550: 0.736 s/cyc | t2s2 12.3/2.00 | t2s3 21.4/2.94 | dt 12.01
- log.1 c3550-3850: 0.832 s/cyc | t2s2 12.9/2.00 | t2s3 16.1/2.71 | dt 12.01
- log.1 c3850-4150: 0.735 s/cyc | t2s2 13.6/2.00 | t2s3 17.1/3.17 | dt 12.01
- log.1 c4150-4450: 0.746 s/cyc | t2s2 15.6/2.41 | t2s3 13.1/3.10 | dt 11.80
- log.1 c4450-4750: 0.783 s/cyc | t2s2 21.5/3.28 | t2s3 11.6/2.44 | dt 10.30
- log.1 c4750-5050: 1.018 s/cyc | t2s2 39.9/4.81 | t2s3 19.9/3.53 | dt 9.95
- log.1 c5050-5350: 1.142 s/cyc | t2s2 52.5/4.88 | t2s3 28.2/3.70 | dt 8.76
- log.1 c5350-5650: 1.170 s/cyc | t2s2 48.6/3.99 | t2s3 23.3/2.97 | dt 12.01
- log.1 c5650-5950: 1.129 s/cyc | t2s2 53.9/4.53 | t2s3 25.4/3.64 | dt 11.04
- log.1 c5950-6250: 1.252 s/cyc | t2s2 63.8/5.28 | t2s3 30.4/4.61 | dt 10.38

## Resume
- m35mms production resumed from the pre-stop rst, mg 2: `2026-11-24 12:52:38 chunk 2 STAR=m35mms BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f md5=7d7c4a71a07d3dc3878e7dea1403df23 DECK=/resnick/groups/carnegie_poc/jingze/prod_rc4v_1124/decks/IA123_m` EX hydro/dt_report=0.
