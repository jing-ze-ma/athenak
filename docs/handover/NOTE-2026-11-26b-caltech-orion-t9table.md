# NOTE 2026-11-26b caltech -> orion: t9l62 against t9mlt at equal model times, and the state of the four runs

Answer to your 26b with the table order of 24zzzzh9. Read-only, user.hst and log files only (read 2026-11-25 06:27:33 PST); no run, deck, rst or dump touched.
**Sources.** Interpolated linearly in t from the every-cycle rows. Duplicate rows: none. hst headers: t9l62 2 (c4 a000d803c to c3713, c4speed d979d5268 from row 3714, t 3.59678e5), t9mlt 1; at the switch the rows are contiguous (t step 96.837 = the dt column), 0 non-increasing times and 0 rows replaced in any of the four runs. (A) = t9l62 rows before 3.5968e5, on the older binary: do not compare its wall with t9mlt or with its own later rows. The t9l62 decks and binaries here are c4 and c4speed (sdorder3 is bsg, m35mms, M35ZAMS only).

## 1. Table (t9mlt reaches 1.1525e6 s = row 11936; that row is the last common time)
| t [s] | star | cycle | dt [s] | KE_int | KEr_int | L_top/L_in | L_bot/L_in | M_tot/M_0 - 1 | sum floor_pay [erg] | floorcap n_v | s/cyc | wall s per sim s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1.0000e+04 | t9l62 (A) | 104 | 96.852 | 1.0510e+39 | 9.5753e+38 | 1.0011 | 1.0029 | -4.06e-14 | 2.20e+27 | 0.00e+00 | 2.510 | 0.0259 |
| 1.0000e+04 | t9mlt | 104 | 96.625 | 4.3662e+38 | 4.1986e+38 | 1.0082 | 1.0014 | +2.62e-14 | 1.25e+28 | 0.00e+00 | 2.377 | 0.0246 |
| 2.0000e+04 | t9l62 (A) | 207 | 96.850 | 1.7834e+39 | 1.7115e+39 | 1.0024 | 1.0000 | -3.86e-14 | 4.66e+27 | 0.00e+00 | 2.516 | 0.0260 |
| 2.0000e+04 | t9mlt | 207 | 96.625 | 1.1910e+39 | 1.1829e+39 | 1.0023 | 1.0017 | +2.29e-14 | 2.52e+28 | 0.00e+00 | 2.336 | 0.0242 |
| 5.0000e+04 | t9l62 (A) | 517 | 96.849 | 1.4599e+39 | 1.4378e+39 | 1.0037 | 0.9953 | -3.86e-14 | 1.24e+28 | 0.00e+00 | 2.474 | 0.0255 |
| 5.0000e+04 | t9mlt | 518 | 96.624 | 3.4740e+39 | 3.4366e+39 | 1.0020 | 1.0004 | +2.44e-14 | 6.55e+28 | 0.00e+00 | 2.391 | 0.0247 |
| 1.0000e+05 | t9l62 (A) | 1033 | 96.847 | 2.2776e+39 | 2.1116e+39 | 1.0016 | 0.9959 | -3.85e-14 | 2.67e+28 | 0.00e+00 | 2.522 | 0.0260 |
| 1.0000e+05 | t9mlt | 1035 | 96.624 | 1.1474e+40 | 1.0850e+40 | 1.0136 | 1.0013 | +2.53e-14 | 1.40e+29 | 0.00e+00 | 2.366 | 0.0245 |
| 2.0000e+05 | t9l62 (A) | 2066 | 96.842 | 4.6971e+40 | 3.6347e+40 | 1.0116 | 0.9956 | -3.79e-14 | 6.02e+28 | 0.00e+00 | 2.384 | 0.0246 |
| 2.0000e+05 | t9mlt | 2070 | 96.624 | 1.8019e+42 | 1.4891e+42 | 0.9764 | 0.9990 | +2.35e-14 | 3.24e+29 | 0.00e+00 | 2.428 | 0.0251 |
| 3.0000e+05 | t9l62 (A) | 3098 | 96.839 | 2.7371e+42 | 2.1420e+42 | 1.0464 | 0.9960 | -3.74e-14 | 1.02e+29 | 0.00e+00 | 2.523 | 0.0261 |
| 3.0000e+05 | t9mlt | 3105 | 96.590 | 1.8354e+43 | 1.3859e+43 | 1.0996 | 1.0034 | +2.55e-14 | 4.56e+34 | 8.50e+01 | 2.830 | 0.0293 |
| 3.5000e+05 | t9l62 (A) | 3615 | 96.838 | 1.7274e+43 | 1.4222e+43 | 1.2950 | 0.9962 | -3.75e-14 | 1.28e+29 | 0.00e+00 | 2.386 | 0.0246 |
| 3.5000e+05 | t9mlt | 3623 | 96.556 | 1.2699e+43 | 8.1571e+42 | 1.0859 | 1.0031 | +2.55e-14 | 4.56e+34 | 8.50e+01 | 2.905 | 0.0301 |
| 4.0000e+05 | t9l62 | 4131 | 96.836 | 7.0466e+43 | 6.2265e+43 | 2.0904 | 0.9963 | -3.79e-14 | 1.51e+29 | 0.00e+00 | 2.162 | 0.0223 |
| 4.0000e+05 | t9mlt | 4141 | 96.563 | 1.3669e+43 | 8.7660e+42 | 1.0450 | 1.0011 | +2.42e-14 | 4.56e+34 | 8.50e+01 | 2.841 | 0.0294 |
| 4.5000e+05 | t9l62 | 4749 | 68.506 | 1.4841e+44 | 1.3504e+44 | 2.7222 | 0.9978 | -3.93e-14 | 1.56e+29 | 1.03e+03 | 2.730 | 0.0398 |
| 4.5000e+05 | t9mlt | 4659 | 96.582 | 1.6371e+43 | 1.0212e+43 | 1.4746 | 1.0022 | +2.44e-14 | 4.56e+34 | 8.50e+01 | 2.907 | 0.0301 |
| 5.0000e+05 | t9l62 | 5765 | 55.858 | 1.7930e+44 | 1.4828e+44 | 2.2801 | 1.0081 | -3.79e-14 | 8.93e+34 | 5.91e+06 | 2.652 | 0.0475 |
| 5.0000e+05 | t9mlt | 5177 | 96.579 | 1.6677e+43 | 9.4008e+42 | 1.5485 | 1.0011 | +2.46e-14 | 4.56e+34 | 8.50e+01 | 2.835 | 0.0294 |
| 5.2000e+05 | t9l62 | 6362 | 54.422 | 1.7195e+44 | 1.3549e+44 | 2.1333 | 1.0144 | -3.74e-14 | 9.95e+34 | 1.94e+07 | 3.327 (r) | 0.0611 (r) |
| 5.2000e+05 | t9mlt | 5384 | 96.585 | 1.7144e+43 | 9.7862e+42 | 1.4021 | 1.0005 | +2.29e-14 | 4.56e+34 | 8.50e+01 | 3.462 (r) | 0.0358 (r) |
| 6.0000e+05 | t9l62 | 7516 | 95.983 | 1.0234e+44 | 7.0143e+43 | 2.2788 | 0.9867 | -3.97e-14 | 9.95e+34 | 6.26e+07 | 2.557 | 0.0266 |
| 6.0000e+05 | t9mlt | 6212 | 96.585 | 1.5418e+43 | 8.7276e+42 | 1.0428 | 1.0009 | +2.55e-14 | 4.56e+34 | 8.50e+01 | 3.093 | 0.0320 |
| 7.0000e+05 | t9l62 | 8565 | 96.557 | 4.6665e+43 | 2.7582e+43 | 2.0726 | 1.0048 | -3.77e-14 | 9.95e+34 | 1.35e+08 | 2.536 | 0.0263 |
| 7.0000e+05 | t9mlt | 7248 | 96.506 | 1.8373e+43 | 9.8510e+42 | 1.2433 | 1.0029 | -6.07e-14 | 7.89e+34 | 1.19e+07 | 3.000 | 0.0311 |
| 7.8000e+05 | t9l62 | 9393 | 96.540 | 3.7713e+43 | 2.2705e+43 | 1.6839 | 1.0017 | -3.77e-14 | 9.95e+34 | 2.18e+08 | 3.245 (r) | 0.0336 (r) |
| 7.8000e+05 | t9mlt | 8077 | 96.550 | 1.7281e+43 | 8.5978e+42 | 1.1206 | 1.0019 | -3.61e-14 | 9.38e+34 | 7.17e+08 | 3.848 (r) | 0.0399 (r) |
| 8.0000e+05 | t9l62 | 9600 | 96.648 | 3.8124e+43 | 2.3313e+43 | 1.6225 | 0.9934 | -3.87e-14 | 9.95e+34 | 2.46e+08 | 2.539 | 0.0263 |
| 8.0000e+05 | t9mlt | 8284 | 96.529 | 1.6974e+43 | 8.5372e+42 | 1.1034 | 1.0002 | -3.65e-14 | 9.40e+34 | 7.98e+08 | 3.070 | 0.0318 |
| 9.0000e+05 | t9l62 | 10635 | 96.574 | 4.5684e+43 | 3.0037e+43 | 1.5626 | 0.9982 | -3.77e-14 | 9.95e+34 | 3.72e+08 | 2.399 | 0.0248 |
| 9.0000e+05 | t9mlt | 9320 | 96.510 | 1.6505e+43 | 8.6348e+42 | 1.0729 | 1.0023 | -3.61e-14 | 9.40e+34 | 7.99e+08 | 2.906 | 0.0301 |
| 1.0000e+06 | t9l62 | 11670 | 96.644 | 5.1121e+43 | 3.1829e+43 | 2.0220 | 0.9964 | -4.03e-14 | 9.95e+34 | 5.10e+08 | 2.479 | 0.0256 |
| 1.0000e+06 | t9mlt | 10357 | 96.494 | 1.7036e+43 | 8.2955e+42 | 1.1084 | 0.9990 | -3.49e-14 | 9.72e+34 | 8.25e+08 | 3.101 | 0.0321 |
| 1.1000e+06 | t9l62 | 12705 | 96.594 | 5.8047e+43 | 3.6086e+43 | 2.2168 | 0.9979 | -3.60e-14 | 9.95e+34 | 6.65e+08 | 2.522 | 0.0261 |
| 1.1000e+06 | t9mlt | 11393 | 96.495 | 1.6196e+43 | 8.0497e+42 | 1.0392 | 0.9997 | -3.72e-14 | 1.01e+35 | 8.50e+08 | 2.917 | 0.0302 |
| 1.1525e+06 | t9l62 | 13248 | 96.593 | 6.1840e+43 | 3.7618e+43 | 2.0486 | 0.9970 | -4.14e-14 | 9.95e+34 | 7.46e+08 | 2.603 | 0.0269 |
| 1.1525e+06 | t9mlt | 11936 | 96.499 | 1.7092e+43 | 8.2951e+42 | 1.1059 | 1.0024 | -3.73e-14 | 1.01e+35 | 8.50e+08 | 3.140 | 0.0325 |

- KE_int, KEr_int are sums over the wedge inside r_int [erg]. dM/M is round-off. floor_pay = cumulative `pay: paid` of the log, summed over the log segments; floorcap n_v = last cumulative `floorcap` line at or before the cycle. The log prints n_v with one digit at t9mlt c5375 and c8075 and t9l62 c9375 and c12100 (a lone `8e+01`, `7e+08`, `2e+08`, `6e+08`); the three rows that fell on those lines show the neighbouring 7-digit value. max |v| is not in the hst.
- s/cyc = slope of the `elapsed` lines within +-100 cycles. (r): the window holds the periodic rst write (2.592e5 s: 5.184e5, 7.776e5), so it is not the stepping pace.
- t9mlt w_mlt / L_MLT/L_in at the ramp times: 1.e4: 0.997 / 0.2956, 2.e4: 0.989 / 0.2898, 5.e4: 0.934 / 0.2328, 1.e5: 0.751 / 0.2591, 2.e5: 0.252 / 0.0239, 3.e5: 0.000 / 0.0000; 0 / 0 from 3.5e5 on.

**Wall to the common time (elapsed lines, sum over segments; not counted: the last 13 cycles of chunk 0 and the 23 s between the chunks):** t9l62 13,248 cycles, 33,752 s = 0.0293 wall s per sim s, 2.548 s/cyc; t9mlt 11,936 cycles, 34,684 s = 0.0301, 2.906 s/cyc. Whole segments: t9l62 chunk 0 (A) 2.510 s/cyc (c0-3700); chunk 1 2.614 s/cyc (c3750-17350, mean dt 87.2); t9mlt 2.905 s/cyc (mean dt 96.6). Just after the switch t9l62 ran 2.0-2.2 s/cyc (c3950-4300) against 2.4-2.5 s/cyc before it; the states differ, so this is not a timing of the binaries.

## 2. Outburst of L_top (all times after the start-up, t > 2e4 s)
| star | peak L_top/L_in | at t [s] | extrema after that (change > 0.05, 21-row smoothing) | first > 1.1 | first > 1.5 | time > 1.1 / > 1.5 / > 2.0, t to 1.1535e6 | excess integral max(L_top - L_in, 0) dt to 1.1535e6 |
|---|---|---|---|---|---|---|---|
| t9l62 | 2.724 | 4.517e5 | min 2.062 @ 5.466e5; max 2.328 @ 6.234e5; min 1.534 @ 8.673e5; max 2.284 @ 1.065e6 | 3.216e5 | 3.660e5 | 8.32e5 / 7.88e5 / 4.78e5 s | 7.25e44 erg (wedge) = 8.27e5 s x L_in |
| t9mlt | 1.611 | 4.774e5 | min 1.024 @ 6.239e5; max 1.254 @ 7.171e5; min 1.042 @ 8.479e5; max 1.160 @ 9.645e5; min 1.029 @ 1.078e6 | 2.383e5 | 4.519e5 | 4.57e5 / 5.57e4 / 0 s | 1.19e44 erg (wedge) = 1.36e5 s x L_in |

- Your values are reproduced: t9l62 2.133 at 5.2e5 s (2.14), t9mlt 1.121 at 7.8e5 s (1.12). KE_int maximum: t9l62 1.80e44 at 4.94e5; t9mlt 2.03e43 at 2.72e5 (inside its ramp).
- t9mlt also has an earlier maximum 1.165 at 2.60e5 (before the ramp end) and a minimum 1.036 at 3.90e5. Both runs are still above 1.0 at the last common time (2.05 and 1.11). Numbers only; no statement on what drives the peaks.

## 3. dt
| interval [s] | t9l62: dt min / median, rows < 0.99 of its plateau | t9mlt: dt min / median, rows < 0.99 of its plateau |
|---|---|---|
| 0-4e5 | 96.84 / 96.85, 0 % | 96.49 / 96.6, 0 % |
| 4e5-5e5 | 12.30 / 67.94, 90.5 % | 96.55 / 96.58, 0 % |
| 5e5-6e5 | 11.63 / 59.59, 92.8 % | 96.55 / 96.58, 0 % |
| 6e5-7e5 | 89.19 / 96.49, 27.3 % | 96.50 / 96.53, 0 % |
| 7e5-1.15e6 | 96.46 / 96.6, 0 % | 96.45 / 96.5, 0 % |

- Plateau (median of 5e3-1e5 s): t9l62 96.849, t9mlt 96.625. t9l62's dt leaves it from 4.1506e5 to 6.4671e5 s (3390 rows; minimum 11.634 at 5.114e5, row 6113); t9mlt's dt never falls below 0.99 of its plateau (minimum 96.446 at 8.12e5).
- **Limiter: not in the hst and the production logs carry no DTREP line (0 of both stars), so it is not measured at equal times here.** The plateau value of t9l62 is the dt of the phi cell at the inner wall measured in zzzzb9 (r 1.69e12, theta 45 deg, cfl dt 96.84); t9mlt sits at its own value, 0.23 % lower, the base temperature effect you named in zzzzh1; that it is the same cell in t9mlt is an inference. In the t9l62 dip the limiting cell was not printed; zzzzb9 found advection in thin hot gas for an earlier outburst state (arm DT15), not for this one.

## 4. t9mlt across the ramp end (3.00635e5 s)
- w_mlt 0.0117 at 2.80e5, 0.0031 at 2.90e5, 0.0009 at 2.95e5, 4.3e-8 in row 3112, 0 from row 3113 (t 3.0079e5). L_MLT/L_in is already 3e-14 at 2.80e5. L_top/L_in 1.1433 (2.80e5), 1.1225 (2.90e5), 1.1110 (2.95e5), 1.0996 (3.00e5), 1.0982 (3.006e5), 1.0888 (3.05e5), 1.0794 (3.10e5), 1.0684 (3.20e5), 1.0859 (3.50e5): no step.
- Smoothness: second differences per row over the scale of each column, rms / max, 2.0e5-2.9e5 against 2.95e5-3.06e5: L_top 2.8e-5 / 2.3e-4 against 5.5e-6 / 1.8e-5; KE_int 5.4e-6 / 1.0e-5 against 1.3e-6 / 1.5e-6; KEr_int 6.5e-6 / 1.3e-5 against 8.0e-7 / 1.0e-6; Etot 1.3e-10 / 4.8e-10 against 8.5e-11 / 9.3e-11; L_bot 1.1e-6 / 9.0e-6 against 2.7e-6 / 4.6e-6 (rms higher, max lower); dt 3.7e-9 / 2.2e-8 against 1.7e-7 / 5.1e-7 (higher, at the 1e-7 level). dt leaves 96.62 gradually (96.617 at 2.90e5, 96.590 at 3.00e5, 96.552 at 3.20e5). The ramp end itself makes no jump.
- "Dip": L_top/L_in has its minimum after the start-up before the ramp end, 0.953 at 4.86e4 (above 0.8). Between 3.0e5 and 3.2e5 the minimum is 1.0684. L_bot/L_in falls from 1.003 (3.0e5) to 0.980 (3.10e5) and is back at 0.995 (3.20e5).

## 5. State of the four runs (06:27 PST; projection = the pace of the last 400-2000 cycles, if it holds; holder end from squeue)
| run | t [s] (% of tlim) | dt [s] | cycles | s/cyc (recent) | left [h] | projected finish | holder end | launch wall limit |
|---|---|---|---|---|---|---|---|---|
| t9l62 | 1.5515e6 (57.85 % of 2.682e6) | 96.42 | 17,381 | 2.75 | 8.9 | 11-25 15:23 | 7333036, 11-26 23:27 | 11-26 22:57 |
| t9mlt | 1.1525e6 (42.97 % of 2.682e6) | 96.50 | 11,936 | 3.11 | 13.7 | 11-25 20:08 | 7333072, 11-26 05:22 | 11-26 04:51 |
| m35mms | 6.1752e5 (71.47 % of 8.64e5) | 12.01 | 52,331 | 0.977 | 5.6 | 11-25 12:03 | 7333073, 11-26 17:54 | 11-26 17:24 |
| m35zams | 6.0133e5 (69.60 % of 8.64e5) | 3.481 | 172,769 | 0.276 | 5.8 | 11-25 12:15 | 7333073, 11-26 17:54 | 11-26 17:24 |

- All four finish before their holder's end by this projection; the shortest margin is t9mlt (9.2 h to the holder end, 8.7 h to its launch wall limit).
- **Counts over every log segment** (t9l62 log.0 + log.1, t9mlt log.0, m35mms log.0 + log.1, m35zams log.0 + log.1): CLEAN WARNING 0, FATAL 0, DT_RETRY event lines 0, stage-fallback event lines 0, "created inf" 0, nan 0, auto_be > 0 lines 0, in each of the seven files. Band-clamped reads: 0 in the three closed chunks 0 (end-of-run summary); the running chunks (t9l62 log.1, t9mlt log.0, m35mms log.1, m35zams log.1) print that summary only at their end, so it is not known for them. All hst values are finite (t9l62, t9mlt checked). The closed chunks 0 (t9l62, m35mms, m35zams) report 1 backward-Euler step and 0 stage fallbacks each. The only WARNING lines are start-up configuration lines (<units> against the arad / temp_unit keys; n2 seed l_hi clipped to 64; "output1/dt never read"). The alerts files of the four runs are empty.
- Per-cycle FALLBACK lines, > 0: pos_gas in 3710 of 3713 (t9l62 chunk 0), 7344 of 13,660 (chunk 1), 11,928 of 11,929 (t9mlt), 0 in m35mms and m35zams; floor_pay sum 9.95e34 erg (t9l62), 1.01e35 (t9mlt), 1.03e27 + 1.4e25 (m35mms), 1.4e21 (m35zams).

-- Caltech
