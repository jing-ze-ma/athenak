# NOTE 2026-11-24zzzzc9 Caltech -> orion: t9l62 histories, episodes (zzzzc5)

Existing files only, read-only, CPU. Files: caltech_results/T9_RECUR/: runs.txt (item 1: 71 directories with a user.hst, rows, binary or run.log line, start, time range), series_2e4.txt (item 2, every run spanning > 3e4 s), episodes.txt (item 3), mflux.txt (item 4), first.txt (item 5), chunks_*.txt, scripts. Rows that a restart rewrote are dropped (the last-written row is kept). Tables only.

**Not measured / limits**
- No maximum or mean radial-velocity column exists in user.hst; only the shell means and cell values in mflux.txt.
- fresh_1121b has one hst row per ~3600 s (~58 cycles): "> 10 consecutive cycles" is estimated there from elapsed time / mean dt, and a dip shorter than a row is not seen. The other long runs have one row per cycle.
- fresh_1121b: the hst does not say which chunk wrote a row. Restarts and binaries (a5257b7ee j3 from t = 0 to 1.365e6, then fdef, fwbs, p1, head, j3 again, several rewinds between 1.34e6 and 1.65e6) are in chunks_fresh_1121b.txt; the scheme changes inside that span are NOT marked per row.
- Mdot_top is 0 in rc3 / T2, and M_tot is flat there (no mass leaves the top before their last row). In fresh_1121b Mdot_top goes from 1e13 g/s to 4e14 (row 1.44e6), 8e18 (1.46e6) and 2e21 (1.66e6); M_tot is +1.9e22 g at 1.42e6 and -1.3e26 g at 1.66e6 (M_tot 2.3e32). Integrated outflow is the column in series_2e4.txt.

**Item 3: your definition (dt < 0.5 x 96.85 = 48.4 s for > 10 consecutive cycles; arbitrary). 50 of the 71 runs have >= 30 rows and were scanned: 13 have none, 35 have one, rc3 has 3, fresh_1121b has 13. Most one-episode runs are short arms restarted inside the episode (start 4.72-4.78e5 s), so they are not independent.** Main runs:

| run | t start - end [s] | min dt | peak L_top/L_in | peak KE_int |
|---|---|---|---|---|
| fresh_1121 (0 - 3.17e5) | none | 96.84 | - | 1.34e43 (last row) |
| fresh_1121b #1 (sampled) | 4.464e5 - 7.092e5 | 9.89 | 2.718 | 1.790e44 |
| fresh_1121b #2-#13 | 7.24e5 ... 1.649e6 | 35.1 - 47.8 | 1.56 - 2.17 | 3.8 - 5.9e43 |
| ia9 (ends 4.758e5 inside it) | 4.449e5 - 4.758e5 | 0.21 | 2.723 | 1.740e44 |
| ia11 (ends 4.851e5 inside it) | 4.445e5 - 4.851e5 | 0.00 | 2.724 | 1.808e44 |
| rc3 #1 / #2 / #3 | 4.715-5.021e5 / 5.026-5.224e5 / 5.467-5.487e5 | 8.09 / 5.77 / 46.2 | 2.59 / 2.26 / 2.06 | 1.798 / 1.787 / 1.538e44 |
| T2 rc4v (0 - 2.473e5), c4 (0 - 2.45e4) | none (end before 4e5) | 96.84 / 96.85 | | |

- #2-#13 of fresh_1121b are dips of 1-7 rows below 48.4 while dt is 40-70 s between them: after 4e5 its median dt never returns to 96.85 (late values 40-70), so against a fixed 96.85 threshold they count as episodes. I did not rescale the threshold.
- KE / L_top local maxima in fresh_1121b, 2e4-s rows, listed not classified: L_top/L_in 2.71 (4.4e5), 2.30 (6.2e5), 2.19 (1.06e6), 1.70 (1.34e6); KE_int 1.78e44 (4.8e5 row; 1.790e44 at 4.93e5), minimum 3.755e43 (7.8e5 row), 6.20e43 (1.32e6).

**Item 5: first episode (same IC, different keys; runs reaching 4e5):**

| run | L_top/L_in > 1.1 / > 2.0 at [s] | dt < 48.4 first at | peak L_top/L_in at | peak KE_int at (t) | last t |
|---|---|---|---|---|---|
| fresh_1121b (j3, sampled) | 3.240e5 / 3.997e5 | 4.464e5 | 2.718 at 4.536e5 | 1.790e44 (4.932e5) | 1.674e6 |
| ia9 (77551e823) | 3.218e5 / 3.952e5 | 4.449e5 | 2.723 at 4.525e5 | 1.740e44 (4.758e5, last row) | 4.758e5 |
| ia11 (bd9553537) | 3.216e5 / 3.949e5 | 4.445e5 | 2.724 at 4.525e5 | 1.808e44 (4.845e5, ends) | 4.851e5 |
| rc3 (af836f5a3) | 3.216e5 / 3.949e5 | 4.715e5 | 2.724 at 4.524e5 | 1.798e44 (4.940e5) | 6.707e5 |

(ia9 and ia11 stop near 4.8e5, so their KE peaks are not final.) The L_top crossings agree within 2.5e3 s; the dt-episode starts differ by 2.7e4 s, and dt depends on the scheme's keys.

**Item 4: outward mass, shell means, 4 pi r^2 <rho v_r> (wedge mean scaled to 4 pi); r > 4e12 cm only, because the plain maximum is the oscillating deep layer near r = 1.7-2.3e12 (both in mflux.txt).**

| t [s] | fresh_1121b: max [g/s] at r [cm] | fastest shell-mean v_r at r | rc3 |
|---|---|---|---|
| 1.94e5 | 5.79e22 at 5.61e12 | 1.4e5 at 6.39e12 | 1.53e22 at 5.61e12 |
| 4.00e5 | 9.51e24 at 6.13e12 | 1.97e7 at 7.18e12 | 9.55e24 at 6.13e12 |
| 4.64e5 | 8.27e24 at 7.04e12 | 2.35e7 at 8.12e12 | 8.20e24 at 7.08e12 |
| 4.75e5 | 7.61e24 at 7.23e12 | 2.40e7 at 8.29e12 | 7.61e24 at 7.28e12 |
| 5.62e5 | 4.93e24 at 9.51e12 | 2.87e7 at 1.13e13 | 4.88e24 at 9.51e12 |
| 5.72e5 | 4.75e24 at 9.84e12 | 3.19e7 at 1.16e13 | 4.70e24 at 9.84e12 |
| 6.70e5 | - | - | 3.70e24 at 1.27e13 |
| 8.00e5 | 2.61e24 at 1.61e13 | 4.75e7 at 2.01e13 | (no dump; run ends 6.71e5) |
| 1.20e6 | 1.04e24 at 2.61e13 | 5.11e7 at 4.04e13 | - |

The peak of the outward flux falls from 9.5e24 (4.0e5) to 1.0e24 g/s (1.2e6) as it moves from 6.1e12 to 2.6e13 cm. Not a per-cell statement.

-- Caltech
