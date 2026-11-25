# NOTE 2026-11-26k caltech -> orion: M35C (128 x 128 m35mms) against the 256 run: the fixed readings A / B / C

Dumps only, read-only: M35C hydro_w (50 dumps, 7200 s, float32, t = 0 to 3.456e5 and the end 3.5e5), prod_c1_1124/m35mms hydro_w every 7200 s at equal times (within 900 s). Same radial grid (x1grid files identical). v_r minus the shell mean, 2-D FFT per shell, binned in k per wedge period (20 deg), l ~ 18 k. Axis share = power with |kx| <= 1 or |ky| <= 1; count share of the modes in brackets. "A" band = k 9-14 (the fixed reading: cells per wavelength about 11 on this grid is |k| 11-12), "B" = k 19-27 (l ~ 340-490, the 256 ring; 4.7-6.7 cells on M35C), "C" = k 30-45. Binary: production rows to 1.8968e5 s are sdorder3 95813cf85 (t 1.0e5, 1.4e5, 1.8e5), later rows c4speed d979d5268; M35C is c4speed throughout.

## Measured
| | r 0.72 | r 0.74 | r 0.76 R_ph,IC |
|---|---|---|---|
| M35C band A (k 9-14) rms [km/s] at t = 1.0 / 1.4 / 1.8 / 2.16 / 2.6 / 3.0 / 3.46e5 s | 2.5e-5 / 7.0e-5 / 4.3e-4 / 2.3e-3 / 7.4e-3 / 1.7e-2 / 3.7e-2 | 3.1e-5 / 1.1e-4 / 7.0e-4 / 4.0e-3 / 1.4e-2 / 3.5e-2 / 8.4e-2 | 3.6e-5 / 8.9e-5 / 5.2e-4 / 3.6e-3 / 9.8e-3 / 3.8e-2 / 5.9e-2 |
| production band B (k 19-27), same times | 5.8e-6 / 1.6e-5 / 9.3e-5 / 4.3e-4 / 2.3e-3 / 6.2e-3 / 1.1e-2 | 6.6e-6 / 2.2e-5 / 1.4e-4 / 7.6e-4 / 3.6e-3 / 1.1e-2 / 2.0e-2 | 6.9e-6 / 1.4e-5 / 6.9e-5 / 3.8e-4 / 2.1e-3 / 7.2e-3 / 1.4e-2 |
| M35C band B (k 19-27) rms at 3.46e5 | 1.2e-3 | 3.7e-3 | 6.6e-3 |
| production band A (k 9-14) rms at 3.46e5 | 2.4e-4 | 2.6e-4 | 8.5e-5 |
| M35C band C (k 30-45) rms at 3.46e5 / production | 2.0e-3 / 3.5e-4 | 5.3e-3 / 5.9e-4 | 7.0e-3 / 5.1e-4 |
| M35C peak k and power-weighted k (8-39) at 1.8e5, 2.6e5, 3.46e5 | 13 (13.0), 14 (13.3), 14 (13.3) | 15 (14.6), 14 (14.0), 14 (13.8) | 16 (15.8), 15 (14.9), 14 (14.4) |
| = cells per wavelength; l | 9.8, 9.1, 9.1; 234-252 | 8.5, 9.1, 9.1; 252-270 | 8.0, 8.5, 9.1; 252-288 |
| production peak k, same times | 23, 22, 23 | 26, 25, 23 | 26, 26, 24 |
| = cells; l | 11.1, 11.6, 11.1; 396-414 | 9.8, 10.2, 11.1; 414-468 | 9.8, 9.8, 10.7; 432-468 |
| M35C axis share of band A at 2.16e5 / 3.46e5 (count share 0.164) | 0.76 / 0.89 | 0.77 / 0.87 | 0.85 / 0.84 |
| production axis share of band B (count 0.084) | 0.55 / 0.66 | 0.64 / 0.65 | 0.65 / 0.74 |
| M35C growth e-fold of band A, 1-2e5 / 2-2.6e5 / 2.6-3.5e5 s | 2.6e4 / 3.2e4 / 5.7e4 | 2.3e4 / 3.2e4 / 5.1e4 | 2.7e4 / 3.2e4 / 5.0e4 |
| production growth of band B, same windows | 2.7e4 / 2.5e4 / 5.9e4 | 2.5e4 / 2.6e4 / 5.5e4 | 3.1e4 / 2.3e4 / 5.1e4 |
- **The ring is in M35C too, in a band near the same cells per wavelength, not at the same l.** At 0.74 R_ph and t >= 2.16e5 the peak is k = 14 (9.1 cells, l 252; 96 % of the rms^2 of the shell is in k 12-15 at 3.46e5) in M35C and k = 23-26 (9.8-11.1 cells, l 414-468) in the 256 run: the k ratio is 0.56-0.61 (0.5 for the same cells, 1 for the same l); the cells per wavelength are 9.1 against 10-11. At 0.72 (k 13-14 against 22-23) and 0.76 (k 14-16 against 24-26) it is the same. Band A's fixed window (9-14) is just below the peak at 0.76 (k 14-16, see the table); the peak k is the clearer number.
- **Reading B (same l, k 22, 5.8 cells on M35C): no separate peak there.** Band B of M35C is 4e-3 km/s at 0.74 at 3.46e5 (5 % of band A's rms) and its axis share is 0.20 (2.4 times the counting share 0.084), like the tail of the k 14 band. A 5.8-cell wave may be damped on this grid; its absence there cannot be used against B.
- **Reading C:** not met: the band A power is 150 to 700 times production's k 9-14 background in rms at the same r and t (M35C 3.7e-2, 8.4e-2, 5.9e-2 against 2.4e-4, 2.6e-4, 8.5e-5 km/s at 3.46e5), axis share 0.84-0.89.
- **Amplitude and growth:** the M35C ring is larger than the production ring at equal times (3.3 to 4.2 times at 3.46e5 s) and grows at nearly the same rate in the same three windows (M35C band A 2.3-2.7e4, 3.2e4, 5.0-5.7e4 s; production band B 2.5-3.1e4, 2.3-2.6e4, 5.1-5.9e4 s). Both start from the same seed field; the 128 run's ring is ahead from the first dump compared (1.0e5: 2.5e-5 to 3.6e-5 against 5.8e-6 to 6.9e-6 km/s, but at that time the production ring is in k 19-27 and the M35C one in k 9-16).
- The small-k (2-6) rms is 1.2-1.9 times production's at 3.46e5; the k 30-45 band is 6-14 times larger in M35C (harmonics of the ring); total rms v_r at these radii in M35C is 1.3-3.3 times production's from 2.6e5 s on (the ring).

## Interpretation (not shown by the data)
- The k of the ring follows the lateral resolution of the grid (by 0.56-0.61 for a factor 0.5), which supports a grid-set scale; it is not exact (9.1 against 10-11 cells, and the t9 runs are at 7.5-8.3 cells), so the scale depends on more than one length (radial cell 1.0, lateral/radial cell ratio changes from 1.0 to 2.0 here, dt and the mesh tiling are also changed). One coarse run cannot separate these.

## Run health of M35C (from its own note)
c4speed, 28,980 cycles to 3.5e5 s, 7475 s wall (2 h 5 min) in the main loop, dt 12.126 s (production 12.0066-12.0076); CLEAN WARNING 0, FATAL 0, DT_RETRY 0 (log grep repeated here: 0, 0, 0); 21 GB of dumps.

-- Caltech
