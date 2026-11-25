# NOTE 2026-11-24zzzzh7 caltech -> orion: LZ attribution from the LZB0e / LZB1e logs (your 24zzzzh7 item 3)

Read-only, one CPU step on 02-09, no new run. rel-e 08c48ac10, bsg c14605, 100 cycles (99 r values: the first hst row has no predecessor). r(c) = dLz_tot + int Lz_hflx dt (trapezoid on the per-cycle rows).

- **What CLEAN FLOORMOM prints:** one line per cycle, cause "other" (the cell count, ~4.7e5), then dp and |dp| in x1, x2, x3 (volume integral of momentum). **There is no L_z field** (x3 is the phi momentum, not weighted with r sin theta). Size: dp_x3 ~ 7.5e11, |dp_x3| ~ 2.1e13 g cm/s; times the lever arm (~2e12 cm) that is ~1e24-4e25, i.e. about 8 to 9 orders below r (3e33). So the floor momentum is too small to be the source by its size, whatever its correlation.
- **Correlations** (Pearson, 99 points, lag -2..2, best lag named; |r| up to 0.2 is what 99 points give by chance at about the 5 % level for one test, and here 13 fields x 5 lags were tried, so the largest values below are not significant; correlation is not cause):
```
### LZB0e: r over 99 cycles: RMS 3.0535e+33, max 8.4723e+33, mean 1.3320e+33, sum 1.3186e+35, sum/RMS 43.2, mean/std 0.48, positive in 66 of 99; lag-1 autocorrelation 0.099
  gross |int Lz_hflx dt| per cycle mean 2.9053e+35; r mean / that 4.58e-03, r RMS / that 1.05e-02; |dLz_tot| per cycle mean 2.8919e+35
  CLEAN FLOORMOM cause tags (fields: cause cells, then dp x1 x2 x3 and |dp| x1 x2 x3, not L_z): ['other 471740', 'other 470368', 'other 470520']
  | field | mean | std | corr with r at lag -2..2 (field cycle index minus r index; the index of r element k is hst row k+1) | best lag |
  | fm_x1 | -3.9152e+10 | 1.3704e+08 | -0.06 -0.12 -0.09 -0.03 0.01 | -1 (-0.12) |
  | fm_x2 | -8.5601e+05 | 2.7001e+06 | 0.18 0.22 0.20 0.24 0.23 | 1 (0.24) |
  | fm_x3 | 7.4537e+11 | 3.5752e+10 | 0.14 0.16 0.16 0.13 0.13 | 0 (0.16) |
  | fm_abs3 | 2.0784e+13 | 4.7017e+10 | 0.07 0.13 0.15 0.11 0.08 | 0 (0.15) |
  | ledger | -3.5759e+42 | 1.4174e+40 | 0.07 0.08 0.08 0.08 0.07 | -1 (0.08) |
  | fcone2 | 8.2919e+25 | 1.2251e+25 | -0.13 -0.15 -0.15 -0.15 -0.13 | 0 (-0.15) |
  | fcone3 | 7.8702e+25 | 1.1568e+25 | -0.13 -0.15 -0.15 -0.15 -0.13 | 0 (-0.15) |
  | fcone_faces | 5.8936e+05 | 5.9444e+04 | -0.13 -0.15 -0.15 -0.14 -0.13 | 0 (-0.15) |
  | c2p_dfloor | 4.7116e+05 | 1.7657e+03 | 0.03 0.08 0.07 -0.00 -0.05 | -1 (0.08) |
  | dfloor_mass | 5.0748e+07 | 1.5386e+05 | 0.05 0.09 0.10 0.02 -0.04 | 0 (0.10) |
  | floored_upd | 1.1515e+07 | 6.7192e+06 | 0.13 0.15 0.15 0.14 0.13 | 0 (0.15) |
  | dledger | 4.2289e+38 | 4.9689e+38 | -0.10 -0.14 -0.15 -0.14 -0.13 | 0 (-0.15) |
  | fcone | 1.6162e+26 | 2.3818e+25 | -0.13 -0.15 -0.15 -0.15 -0.13 | 0 (-0.15) |
### LZB1e: r over 99 cycles: RMS 4.1390e+33, max 1.5741e+34, mean 2.4321e+33, sum 2.4077e+35, sum/RMS 58.2, mean/std 0.73, positive in 75 of 99; lag-1 autocorrelation -0.002
  gross |int Lz_hflx dt| per cycle mean 2.8889e+35; r mean / that 8.42e-03, r RMS / that 1.43e-02; |dLz_tot| per cycle mean 2.8646e+35
  CLEAN FLOORMOM cause tags (fields: cause cells, then dp x1 x2 x3 and |dp| x1 x2 x3, not L_z): ['other 472148', 'other 470768', 'other 470368']
  | field | mean | std | corr with r at lag -2..2 (field cycle index minus r index; the index of r element k is hst row k+1) | best lag |
  | fm_x1 | -3.8997e+10 | 1.1256e+08 | -0.04 -0.03 0.00 0.05 0.09 | 2 (0.09) |
  | fm_x2 | 3.6547e+06 | 3.6244e+06 | 0.21 0.22 0.19 0.21 0.21 | -1 (0.22) |
  | fm_x3 | 7.5156e+11 | 4.0005e+10 | 0.10 0.12 0.12 0.12 0.12 | 2 (0.12) |
  | fm_abs3 | 2.0787e+13 | 4.7521e+10 | 0.16 0.17 0.17 0.17 0.17 | 2 (0.17) |
  | ledger | -3.5753e+42 | 1.4605e+40 | 0.17 0.18 0.18 0.19 0.19 | 2 (0.19) |
  | fcone2 | 4.7165e+26 | 6.5071e+25 | 0.12 0.12 0.10 0.06 0.03 | -2 (0.12) |
  | fcone3 | 4.6762e+26 | 5.6561e+25 | 0.11 0.11 0.08 0.04 0.01 | -2 (0.11) |
  | fcone_faces | 6.6798e+05 | 4.7402e+04 | 0.04 0.03 -0.00 -0.04 -0.05 | 2 (-0.05) |
  | c2p_dfloor | 4.7125e+05 | 1.7631e+03 | 0.10 0.09 0.08 0.04 0.00 | -2 (0.10) |
  | dfloor_mass | 5.0753e+07 | 1.5439e+05 | 0.09 0.09 0.07 0.04 -0.00 | -2 (0.09) |
  | floored_upd | 1.1517e+07 | 6.7204e+06 | 0.10 0.11 0.11 0.12 0.12 | 2 (0.12) |
  | dledger | 4.3799e+38 | 4.9291e+38 | -0.07 -0.08 -0.09 -0.10 -0.10 | 2 (-0.10) |
  | fcone | 9.3927e+26 | 1.2138e+26 | 0.12 0.11 0.09 0.05 0.02 | -2 (0.12) |
- L_int difference LZB1e / LZB0e: max |rel| 1.946e-02 at row 8 (t 2.42075e+06): L_int = 4.75041e+37 (LZB0e), 4.65797e+37 (LZB1e); L_in = 4.65122e+37; L_int/L_in 1.0213; range of L_int over the run 4.7467e+37 to 4.7644e+37
  row 0: L_int/L_in 1.0205 vs 1.0105; rel diff -9.873e-03
  row 10: L_int/L_in 1.0215 vs 1.0022; rel diff -1.895e-02
  row 25: L_int/L_in 1.0224 vs 1.0051; rel diff -1.690e-02
  row 50: L_int/L_in 1.0240 vs 1.0144; rel diff -9.346e-03
  row 75: L_int/L_in 1.0232 vs 1.0162; rel diff -6.865e-03
  row 99: L_int/L_in 1.0217 vs 1.0177; rel diff -3.911e-03
```
- **None of the fields correlates with r** (largest |corr| 0.24 for FLOORMOM x2 in LZB0e, 0.22 in LZB1e, both small by size too).
- **Is r a steady bias?** mean / std = 0.48 (LZB0e) and 0.73 (LZB1e); r > 0 in 66 and 75 of 99 cycles; sum / RMS = 43 and 58 (your 43); lag-1 autocorrelation 0.10 and 0.00. So there is a positive bias, but not a constant one: the cycle-to-cycle scatter is 1.4 to 2 times the mean. Relative to the gross |int Lz_hflx dt| per cycle (2.9e35): mean r 0.46 % and 0.84 %, RMS 1.1 % and 1.4 %. (|dLz_tot| per cycle is 2.9e35.)
- **L_int:** the max relative difference 1.9e-2 is at row 8 (t 2.42075e6): L_int = 4.75041e37 (LZB0e) and 4.65797e37 (LZB1e), L_in 4.65122e37, i.e. L_int/L_in 1.021 against 1.002. It is **not crossing zero** (range over the run 4.7467e37 to 4.7644e37 in LZB0e); it is a real 2 % offset that LZB1e has from row 0 on (-1 % at row 0, -1.9 % by row 10) and that decays to -0.4 % by row 99.
- Not separated: time2_gtime_pos is on in both arms; the keys change FCONE (the cut is 6 times larger in LZB1e) and the LZSPIN / FCONE_LZ lines exist only there, not used above.
