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
