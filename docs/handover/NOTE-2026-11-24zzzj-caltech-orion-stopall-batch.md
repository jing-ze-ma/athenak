# NOTE 2026-11-24zzzj caltech -> orion: stop of RSG + 4 stars and the test batch (your 24zzzh / 24zzzi)
- (11-24 13:07) **m35mms stopped** (clean: yes): final cycle 6959, t 7.9786351968706047e+04 s, last chunk wall elapsed=8.461931e+02 s (log.2); final rst m35mms.00006.rst = rst_keep/m35mms.stop1124b.rst; STOP removed. CW 0, FATAL 0, DT_RETRY 0.
- (11-24 13:07) **m35zams stopped** (clean: yes): final cycle 22155, t 7.7110640711687898e+04 s, last chunk wall elapsed=4.844671e+03 s (log.1); final rst m35zams.00005.rst = rst_keep/m35zams.stop1124b.rst; STOP removed. CW 0, FATAL 0, DT_RETRY 0.
- (11-24 13:09) **rsg stopped** (clean: yes): final cycle 797110, t 7.4606400000000000e+08 s, last chunk wall elapsed=4.145019e+04 s (log.1); final rst rsg.01763.rst = rst_keep/rsg.stop1124b.rst; STOP removed. CW 0, FATAL 0, DT_RETRY 0.
- (11-24 13:09) **t9l62 stopped** (clean: yes): final cycle 2553, t 2.4724774513010253e+05 s, last chunk wall elapsed=1.738359e+03 s (log.2); final rst t9l62.00007.rst = rst_keep/t9l62.stop1124b.rst; STOP removed. CW 0, FATAL 0, DT_RETRY 0.

## CSHERO pre-check time/nlim=0: passed: A B C.

## m35zams-Z1 (11-24 13:12)
| arm | hst rows | s/cyc (cycles) | dt min/mean/max | its t2s2 (sum, solves) | its t2s3 | passes/cyc | L_top/L_in end (mean) | L_bot/L_in | KE_int end (drift) | dM/M | dE/E | CW/FATAL/inf/DT_RETRY | FWB / FCONE lines |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Z1f | 300 | 0.2963 (c22160-22450) | 3.481 / 3.481 / 3.481 | 18.8, 2.00 | 17.6, 2.34 | - | 1.00009 (1.00013) | 1.00001 | 4.96109e+35 (5.29e-02) | -3.10e-15 | -2.29e-08 | 0 / 0 / 0 / 0 | 0 / 0 |
| Z1t | 300 | 0.2920 (c22160-22450) | 3.481 / 3.481 / 3.481 | 19.1, 2.00 | 16.8, 2.27 | - | 1.00009 (1.00013) | 1.00001 | 4.96108e+35 (5.29e-02) | 7.44e-16 | -2.29e-08 | 0 / 0 / 0 / 0 | 0 / 0 |

Pairs on the common first N rows: last-row relative difference (second - first)/first; time, dt, L_top, L_in, KE_int, Etot, M_tot, Lz_tot; and max |L_top/L_in difference| over the rows:
- Z1t vs Z1f (N 300): t -7.65e-14 dt 2.06e-11 L_top -2.97e-07 L_in 0.00e+00 KE -1.19e-06 Etot -7.12e-11 M 1.49e-15 Lz -2.06e-03 | max d(L_top/L_in) 6.37e-05 | dt ratio mean 1.0000

## m35mms-M1-M2-M3 (11-24 13:18)
| arm | hst rows | s/cyc (cycles) | dt min/mean/max | its t2s2 (sum, solves) | its t2s3 | passes/cyc | L_top/L_in end (mean) | L_bot/L_in | KE_int end (drift) | dM/M | dE/E | CW/FATAL/inf/DT_RETRY | FWB / FCONE lines |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| M1 | 100 | 1.1840 (c7000-7050) | 12.01 / 12.01 / 12.01 | 56.6, 5.31 | 25.0, 4.07 | 9.38 | 0.98793 (0.98845) | 1.00015 | 1.57604e+39 (-9.79e-03) | 3.78e-15 | 4.68e-07 | 0 / 0 / 0 / 0 | 200 / 200 |
| M2 | 100 | 1.0604 (c7000-7050) | 12.01 / 12.01 / 12.01 | 37.4, 5.21 | 17.5, 4.02 | 9.23 | 0.98794 (0.98845) | 1.00015 | 1.57604e+39 (-9.79e-03) | 2.83e-15 | 4.68e-07 | 0 / 0 / 0 / 0 | 200 / 200 |
| M3 | 300 | 1.0143 (c7000-7250) | 9.862 / 11.5 / 12.01 | 33.8, 4.98 | 16.0, 3.97 | 8.95 | 0.99957 (0.99228) | 1.00016 | 1.57167e+39 (-1.25e-02) | 1.89e-15 | 9.32e-07 | 0 / 0 / 0 / 0 | 600 / 600 |

Pairs on the common first N rows: last-row relative difference (second - first)/first; time, dt, L_top, L_in, KE_int, Etot, M_tot, Lz_tot; and max |L_top/L_in difference| over the rows:
- M2 vs M1 (N 100): t 3.07e-13 dt -3.09e-10 L_top 4.22e-06 L_in -6.91e-16 KE -1.19e-07 Etot 3.73e-11 M -1.13e-15 Lz -1.07e-06 | max d(L_top/L_in) 1.68e-04 | dt ratio mean 1.0000
- M3 vs M1 (N 100): t 5.47e-13 dt -4.59e-11 L_top 2.48e-05 L_in -8.06e-16 KE 3.32e-07 Etot -1.36e-10 M -1.51e-15 Lz 5.45e-08 | max d(L_top/L_in) 2.16e-04 | dt ratio mean 1.0000
- M3 vs M2 (N 100): t 2.40e-13 dt 2.63e-10 L_top 2.05e-05 L_in -1.15e-16 KE 4.51e-07 Etot -1.73e-10 M -3.78e-16 Lz 1.13e-06 | max d(L_top/L_in) 2.35e-04 | dt ratio mean 1.0000

## Compton C0 vs C1: dt over the first 100 cycles after the rst (your 24zzzm); CW / FATAL / DT_RETRY 0 in both over 300 cycles

| cycles from the rst | C0 dt | C1 dt | C0 L_top/L_in | C1 L_top/L_in |
|---|---|---|---|---|
| 0 | 13.94 | 13.94 | 2.5554 | 2.5554 |
| 10 | 24.97 | 26.45 | 2.5528 | 2.5522 |
| 20 | 14.78 | 15.74 | 2.5507 | 2.5500 |
| 30 | 17.84 | 16.75 | 2.5489 | 2.5483 |
| 40 | 13.96 | 17.13 | 2.5476 | 2.5463 |
| 50 | 13.87 | 22.61 | 2.5460 | 2.5442 |
| 60 | 13.67 | 28.92 | 2.5444 | 2.5411 |
| 70 | 12.99 | 42.66 | 2.5428 | 2.5368 |
| 80 | 29.44 | 23.91 | 2.5405 | 2.5326 |
| 90 | 32.66 | 31.23 | 2.5371 | 2.5299 |
| 99 | 29.75 | 25.32 | 2.5336 | 2.5275 |
- first 100 rows: C0 dt min/mean 10.51 / 19.89, C1 13.94 / 25.67 (mean ratio 1.290); model time advanced C0 1975.5 s, C1 2552.8 s
- first 300 rows: C0 dt min/mean 10.51 / 32.88, C1 13.94 / 35.12 (mean ratio 1.068); model time advanced C0 9848.7 s, C1 10521.3 s

## bsg-thinf-fresh (BSGf = the T2 fresh production c0-300, 02-11, outputs on) (11-24 13:23)
| arm | hst rows | s/cyc (cycles) | dt min/mean/max | its t2s2 (sum, solves) | its t2s3 | passes/cyc | L_top/L_in end (mean) | L_bot/L_in | KE_int end (drift) | dM/M | dE/E | CW/FATAL/inf/DT_RETRY | FWB / FCONE lines |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|

## CSHERO (24zzn / 24zzzg / 24zzzh) on RSG's nodes (11-24 13:27)
- Pre-check passed: A B C. B on 02-03, A / C2 run: B at c600 (t 6.94 d) still has the cube pattern: max m4 of drho 5.7e-4 (800-850 Rsun), 6.0e-4 (850-900), growing since c150; the hero's X was ~3e-4 at c1200. So A and C2 run (they are needed as the comparison; absolute values are not an imprint verdict).. Binary md5 7d7c4a71a07d3dc3878e7dea1403df23, fresh, decks as given, no command-line keys.

### B
- 2026-11-24 13:18:12 rc=0 wall=498 s; CW 0, FATAL 0, created inf 0, DT_RETRY 0; s/cyc 0.4011 (c10-1200) (4 H200)
- imprint (last dump):
```
   850-900    7.332e-04  4.165e-04  1.758e-03  3.689e-03
   900-1000   1.416e-04  9.891e-05  2.416e-04  4.807e-04
  1000-2000   1.836e-06  2.048e-06  7.242e-06  1.176e-05
run_cs/bin/rsg.hydro_w.00008.bin cycle 1200 t = 13.889 d
  r bin [Rsun]   max m4      max m8      max sym_rms  max rms_all
     0-100    1.482e-04  5.417e-05  8.635e-05  8.664e-05
   100-200    1.252e-04  6.288e-05  7.150e-05  7.162e-05
   200-400    3.084e-04  7.056e-05  1.657e-04  1.658e-04
   400-600    7.333e-05  4.535e-05  5.167e-05  5.171e-05
   600-800    1.905e-05  1.533e-05  2.898e-05  6.953e-05
   800-850    1.532e-03  1.437e-03  2.836e-03  6.415e-03
   850-900    7.848e-04  6.836e-04  1.764e-03  4.257e-03
   900-1000   2.228e-04  1.401e-04  3.978e-04  5.466e-04
  1000-2000   2.179e-06  2.369e-06  8.295e-06  1.348e-05
run_cs/bin/rsg.hydro_w.00009.bin cycle 1200 t = 13.889 d
  r bin [Rsun]   max m4      max m8      max sym_rms  max rms_all
     0-100    1.482e-04  5.417e-05  8.635e-05  8.664e-05
   100-200    1.252e-04  6.288e-05  7.150e-05  7.162e-05
   200-400    3.084e-04  7.056e-05  1.657e-04  1.658e-04
   400-600    7.333e-05  4.535e-05  5.167e-05  5.171e-05
   600-800    1.905e-05  1.533e-05  2.898e-05  6.953e-05
   800-850    1.532e-03  1.437e-03  2.836e-03  6.415e-03
   850-900    7.848e-04  6.836e-04  1.764e-03  4.257e-03
   900-1000   2.228e-04  1.401e-04  3.978e-04  5.466e-04
  1000-2000   2.179e-06  2.369e-06  8.295e-06  1.348e-05
```

### A
- 2026-11-24 13:27:22 rc=0 wall=471 s; CW 0, FATAL 0, created inf 0, DT_RETRY 0; s/cyc 0.3746 (c10-1200) (4 H200)
- imprint (last dump):
```
   850-900    8.276e-04  1.591e-04  1.137e-03  1.940e-03
   900-1000   1.316e-04  4.748e-05  2.198e-04  2.475e-04
  1000-2000   1.216e-06  5.119e-07  6.258e-06  7.386e-06
run_cs/bin/rsg.hydro_w.00008.bin cycle 1200 t = 13.889 d
  r bin [Rsun]   max m4      max m8      max sym_rms  max rms_all
     0-100    1.482e-04  5.417e-05  8.635e-05  8.664e-05
   100-200    1.252e-04  6.288e-05  7.150e-05  7.162e-05
   200-400    3.084e-04  7.056e-05  1.657e-04  1.658e-04
   400-600    7.333e-05  4.535e-05  5.167e-05  5.171e-05
   600-800    1.905e-05  1.533e-05  2.898e-05  6.953e-05
   800-850    1.508e-03  1.511e-03  2.784e-03  6.235e-03
   850-900    8.186e-04  1.933e-04  9.770e-04  2.123e-03
   900-1000   1.721e-04  6.609e-05  3.497e-04  3.786e-04
  1000-2000   1.534e-06  5.998e-07  7.191e-06  8.405e-06
run_cs/bin/rsg.hydro_w.00009.bin cycle 1200 t = 13.889 d
  r bin [Rsun]   max m4      max m8      max sym_rms  max rms_all
     0-100    1.482e-04  5.417e-05  8.635e-05  8.664e-05
   100-200    1.252e-04  6.288e-05  7.150e-05  7.162e-05
   200-400    3.084e-04  7.056e-05  1.657e-04  1.658e-04
   400-600    7.333e-05  4.535e-05  5.167e-05  5.171e-05
   600-800    1.905e-05  1.533e-05  2.898e-05  6.953e-05
   800-850    1.508e-03  1.511e-03  2.784e-03  6.235e-03
   850-900    8.186e-04  1.933e-04  9.770e-04  2.123e-03
   900-1000   1.721e-04  6.609e-05  3.497e-04  3.786e-04
  1000-2000   1.534e-06  5.998e-07  7.191e-06  8.405e-06
```

### C
- 2026-11-24 13:25:45 rc=0 wall=374 s; CW 0, FATAL 0, created inf 0, DT_RETRY 0; s/cyc 0.2988 (c10-1200) (4 H200)
- imprint (last dump):
```
   850-900    9.671e-04  1.422e-03  6.708e-03  9.195e-03
   900-1000   8.348e-05  2.118e-04  1.094e-03  1.501e-03
  1000-2000   4.755e-06  6.153e-06  2.711e-05  3.939e-05
run_cs/bin/rsg.hydro_w.00008.bin cycle 1200 t = 14.047 d
  r bin [Rsun]   max m4      max m8      max sym_rms  max rms_all
     0-100    1.501e-04  5.562e-05  8.686e-05  8.717e-05
   100-200    1.253e-04  5.261e-05  7.179e-05  7.192e-05
   200-400    3.084e-04  6.952e-05  1.661e-04  1.661e-04
   400-600    7.500e-05  4.612e-05  5.259e-05  5.263e-05
   600-800    1.950e-05  1.533e-05  2.883e-05  6.911e-05
   800-850    1.295e-03  1.590e-03  8.659e-03  1.273e-02
   850-900    1.415e-03  1.186e-03  7.625e-03  1.069e-02
   900-1000   1.500e-04  2.780e-04  1.315e-03  1.733e-03
  1000-2000   5.589e-06  7.315e-06  3.229e-05  4.688e-05
run_cs/bin/rsg.hydro_w.00009.bin cycle 1200 t = 14.047 d
  r bin [Rsun]   max m4      max m8      max sym_rms  max rms_all
     0-100    1.501e-04  5.562e-05  8.686e-05  8.717e-05
   100-200    1.253e-04  5.261e-05  7.179e-05  7.192e-05
   200-400    3.084e-04  6.952e-05  1.661e-04  1.661e-04
   400-600    7.500e-05  4.612e-05  5.259e-05  5.263e-05
   600-800    1.950e-05  1.533e-05  2.883e-05  6.911e-05
   800-850    1.295e-03  1.590e-03  8.659e-03  1.273e-02
   850-900    1.415e-03  1.186e-03  7.625e-03  1.069e-02
   900-1000   1.500e-04  2.780e-04  1.315e-03  1.733e-03
  1000-2000   5.589e-06  7.315e-06  3.229e-05  4.688e-05
```
- Files: cshero_results/. B hydro_w c0/600/1200 kept on Caltech (/resnick/groups/carnegie_poc/jingze/cshero_1124/B/run_cs/bin/rsg.hydro_w.{00000,00004,00008}.bin), other dumps deleted.

## bsg-thinf-fresh (BSGf = T2 fresh production c0-300 on 02-11 with outputs; BSGt on 01-15, s/cyc not comparable) (11-24 13:29)
| arm | hst rows | s/cyc (cycles) | dt min/mean/max | its t2s2 (sum, solves) | its t2s3 | passes/cyc | L_top/L_in end (mean) | L_bot/L_in | KE_int end (drift) | dM/M | dE/E | CW/FATAL/inf/DT_RETRY | FWB / FCONE lines |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| BSGf | 300 | 0.6963 (c0-300) | 164.6 / 168.2 / 168.8 | 15.5, 3.65 | 19.9, 3.89 | - | 1.06351 (8.52386) | 0.99543 | 1.04088e+41 (nan) | 3.10e-16 | -3.12e-04 | 0 / 0 / 0 / 0 | 0 / 0 |
| BSGt | 301 | 0.7833 (c0-300) | 164.6 / 168.2 / 168.8 | 15.9, 3.61 | 20.6, 3.92 | - | 1.06356 (8.51445) | 0.99551 | 1.03847e+41 (nan) | -4.66e-16 | -3.12e-04 | 0 / 0 / 0 / 0 | 0 / 0 |

Pairs on the common first N rows: last-row relative difference (second - first)/first; time, dt, L_top, L_in, KE_int, Etot, M_tot, Lz_tot; and max |L_top/L_in difference| over the rows:
- BSGt vs BSGf (N 300): t 3.46e-10 dt 1.84e-10 L_top 6.70e-05 L_in 0.00e+00 KE -3.99e-07 Etot -2.10e-10 M 2.02e-15 Lz -1.14e-02 | max d(L_top/L_in) 4.87e+00 | dt ratio mean 1.0000

## t9l62-T1-T2 (11-24 13:33)
| arm | hst rows | s/cyc (cycles) | dt min/mean/max | its t2s2 (sum, solves) | its t2s3 | passes/cyc | L_top/L_in end (mean) | L_bot/L_in | KE_int end (drift) | dM/M | dE/E | CW/FATAL/inf/DT_RETRY | FWB / FCONE lines |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T1 | 300 | 2.3968 (c2600-2850) | 96.84 / 96.84 / 96.84 | 6.1, 3.00 | 8.1, 2.21 | 5.21 | 1.02447 (1.01754) | 0.99594 | 1.05974e+42 (2.25e+00) | -7.73e-16 | -2.07e-06 | 0 / 0 / 0 / 0 | 600 / 600 |
| T2 | 300 | 2.4024 (c2600-2850) | 96.84 / 96.84 / 96.84 | 5.7, 2.97 | 8.0, 2.24 | 5.21 | 1.02445 (1.01764) | 0.99594 | 1.05974e+42 (2.25e+00) | 3.09e-15 | -2.08e-06 | 0 / 0 / 0 / 0 | 600 / 600 |

Pairs on the common first N rows: last-row relative difference (second - first)/first; time, dt, L_top, L_in, KE_int, Etot, M_tot, Lz_tot; and max |L_top/L_in difference| over the rows:
- T2 vs T1 (N 300): t 9.59e-11 dt 8.62e-10 L_top -1.57e-05 L_in 1.72e-16 KE 2.93e-06 Etot -1.13e-08 M 7.73e-16 Lz 6.19e-07 | max d(L_top/L_in) 3.05e-03 | dt ratio mean 1.0000

## RSG-thinf-from-01809 (rc3) (11-24 13:45)
| arm | hst rows | s/cyc (cycles) | dt min/mean/max | its t2s2 (sum, solves) | its t2s3 | passes/cyc | L_top/L_in end (mean) | L_bot/L_in | KE_int end (drift) | dM/M | dE/E | CW/FATAL/inf/DT_RETRY | FWB / FCONE lines |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| RSGf | 100 | 1.4913 (c817400-817450) | 956.5 / 956.5 / 956.5 | 6.9, 4.54 | 6.3, 5.55 | - | 0.94473 (0.94668) | 1.07306 | 2.07574e+46 (-2.11e-04) | -5.79e-14 | 2.09e-07 | 0 / 0 / 0 / 0 | 0 / 0 |
| RSGt | 100 | 1.5208 (c817400-817450) | 956.5 / 956.5 / 956.5 | 7.9, 5.00 | 6.7, 5.62 | - | 0.94573 (0.94811) | 1.07306 | 2.07574e+46 (-2.11e-04) | -5.91e-14 | 2.04e-07 | 0 / 0 / 0 / 0 | 0 / 0 |

Pairs on the common first N rows: last-row relative difference (second - first)/first; time, dt, L_top, L_in, KE_int, Etot, M_tot, Lz_tot; and max |L_top/L_in difference| over the rows:
- RSGt vs RSGf (N 100): t 0.00e+00 dt 0.00e+00 L_top 1.06e-03 L_in 0.00e+00 KE 2.53e-11 Etot -5.68e-09 M -1.24e-15 Lz -4.81e-05 | max d(L_top/L_in) 3.11e-03 | dt ratio mean 1.0000

## t9l62-C0t (C0 = Compton-arm baseline thinf false, C0t = same + thinf=true, small-dt phase, t9l62.00015) (11-24 13:52)
| arm | hst rows | s/cyc (cycles) | dt min/mean/max | its t2s2 (sum, solves) | its t2s3 | passes/cyc | L_top/L_in end (mean) | L_bot/L_in | KE_int end (drift) | dM/M | dE/E | CW/FATAL/inf/DT_RETRY | FWB / FCONE lines |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C0 | 300 | 3.0014 (c5300-5550) | 10.51 / 32.88 / 46.91 | 14.2, 4.84 | 8.4, 3.61 | - | 2.44083 (2.50902) | 0.98459 | 1.78357e+44 (2.72e-02) | 3.09e-16 | -6.01e-05 | 0 / 0 / 0 / 0 | 0 / 0 |
| C0t | 300 | 2.9804 (c5300-5550) | 11.07 / 33.69 / 46.96 | 14.5, 4.83 | 8.5, 3.54 | - | 2.43799 (2.50604) | 0.98506 | 1.78438e+44 (2.76e-02) | 1.55e-16 | -6.15e-05 | 0 / 0 / 0 / 0 | 0 / 0 |

Pairs on the common first N rows: last-row relative difference (second - first)/first; time, dt, L_top, L_in, KE_int, Etot, M_tot, Lz_tot; and max |L_top/L_in difference| over the rows:
- C0t vs C0 (N 300): t 5.02e-04 dt 1.63e-02 L_top -1.16e-03 L_in 0.00e+00 KE 4.51e-04 Etot -1.45e-06 M 3.09e-16 Lz -2.12e-06 | max d(L_top/L_in) 8.64e-03 | dt ratio mean 1.0247
