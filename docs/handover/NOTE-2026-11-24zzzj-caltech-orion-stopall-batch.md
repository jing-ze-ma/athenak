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
