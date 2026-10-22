# NOTE 2026-10-22 caltech -> viper: your interior cube-pattern metric on our N=128 runs: SAME pattern, SAME amplitude as N=32 (re e4dda86f)

Your scr/interior.py (masks unchanged) on Caltech prodF dumps, N=128, f55be371 (arith face): RSG E 0-102 d (every 3 d), HG 0-3 d.
Full tables: docs/coord/interior_n128/{rsgE,hg}_prodF_arith.txt.

## RSG E N=128 (r/R_ph0 0.7 and 0.9; columns: v_r z(FC-ED), rms_FC/rms_ED - 1, v_r rms [cm/s] | Trad z)
| t (d) | 0.7 R_ph | 0.9 R_ph |
|---|---|---|
| 0 | IC: Trad/e z +0.22 (geometric) | +0.22 |
| 3 | +0.22 +0.31 6.9e2 | +0.22 +0.30 1.7e2 |
| 21 | +0.22 +0.33 4.4e3 | +0.22 +0.31 1.1e3 |
| 39 | +0.21 +0.34 5.6e3 | +0.23 +0.33 3.2e3 |
| 51-66 | z falls to +0.06 as resolved convection develops | +0.20, ratio +0.11-0.17 |
| 75-84 | re-emerges +0.27/+0.24 | +0.21 +0.18-0.24 |
| 90-102 | +0.22 +0.27 1.4e4 | **+0.58 -> +0.79, ratio +1.0, Trad z -0.38** (strong cube-locked mean flow at 0.9 R_ph) |
- From 3 d: v_r z = +0.22, rms ratio +0.30 on every shell = your N=32 numbers (+0.21, +0.30). v_r rms at 0.9 R_ph 1.7e2 (3 d) and
  3.2e3 cm/s (39 d) = your N=32 values (1.7e2, 3.2e3). => amplitude independent of N (4x lateral resolution): set by cell SHAPE,
  consistent with a non-well-balanced HSE on cs cells.
- After ~87 d (N=128, cliff runaway active since ~49 d) the 0.9 R_ph flow becomes strongly cube-locked (z 0.8, FC 2x faster).
## HG N=128 0-3 d: same signature from 1 d (v_r z +0.22, ratio +0.30-0.35 at 0.3-0.95 R_ph).

## Caltech plan (USER go): both globals restart FRESH now with implicit_face_opac = harmonic (fixes cliff runaway + dt; HG N=128
Kirchhoff restart cut its imprint -2.4 % -> -1.5 %). The interior cube-locked flow will remain until a well-balanced cs HSE fix;
your static HSE test (N=32/64, with/without well-balancing) is the key next result. If it reproduces, a second fresh start follows.
