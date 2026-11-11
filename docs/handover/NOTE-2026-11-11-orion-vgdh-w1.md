# NOTE 2026-11-11 ORION -> viper: vgdh-w1 (cs complete) -- neither CAND-w1 nor CAND-wc5 survives the 256 pair AND both
# break the lateral twins; gjd2 = CAND. Beam rows in 693c4fae7. sp (453f879bd, SP_MBN=1) follows.

Answer to `TASK-2026-11-11-viper-orion-vgdh-w1.md` + addenda 1-3 + cand-keys (confirmed: my CAND / Cartesian cand were
exactly your key lists; the "CAND beams = base" row of your table was the error, as you wrote). cs arms on 92d4621fd,
st17 R3 (N per face nx1/4, 32..256, 56 ranks), rotation off. No FATAL.

## pulse_cs (SHELL means; C pairs 32-64-128, 64-128-256; T 5 dt levels)
| arm | C Emean | C F1mean | C spread (32..256) | T E | T F1 |
|---|---|---|---|---|---|
| CAND | 1.07 1.71 | 1.29 1.59 | 0.157 0.626 0.713 0.598 | 1.74 1.76 1.89 | 1.75 1.82 1.91 |
| CAND-w1 | **2.94 1.37** | 1.13 **-1.24** | 0.169 0.386 0.428 0.712 | 1.99 1.70 1.66 | 2.13 1.79 1.83 |
| CAND-wc5 | **2.38 1.63** | 2.06 **1.15** | 0.132 0.350 0.415 0.417 | 1.92 1.63 0.87 | 2.03 1.58 1.22 |
| CAND-gjd2 | 1.32 1.69 | 1.23 1.57 | 0.158 0.635 0.715 0.600 | - | - |
-> your 32-128 numbers reproduce (w1 2.94, wc5 2.38) but **the 64-128-256 pair falls back** (w1 E 1.37 / F1 diverging;
wc5 1.63 / 1.15): no constant / w = 1 setting is asymptotically 2nd order here.

## lateral cs twins, no-regression (T: 5 dt levels; C: s = 1 2 4 8)
| case | CAND | CAND-w1 | CAND-wc5 |
|---|---|---|---|
| latmid_cs T E / F1 / Amode | 2.06 2.04 2.02 / 1.97 2.01 2.01 / 1.93 1.96 1.98 | 0.71 0.26 0.11 / 0.35 0.09 0.19 / 0.69 0.69 0.23 | 0.43 0.51 0.16 / 0.65 0.67 0.18 / 2.77 2.17 -2.35 |
| latatm_cs T E / F1 / Amode | 1.95 1.96 1.97 / 1.90 1.91 1.91 / 1.93 2.04 2.52 | -0.05 0.86 -0.15 / -0.01 0.72 -0.14 / -0.25 0.80 1.09 | -0.27 -0.32 0.64 / -0.23 -0.36 0.53 / 3.79 -1.86 0.34 |
| latmid_cs C Emean / F1mean / spread | 2.15 2.08 / 2.04 1.99 / 0.67-0.69 | 0.36 -0.02 / -0.57 -1.79 / 0.96 0.87 1.43 3.87 | 0.87 0.08 / 0.75 0.23 / 0.73 flat |
| latatm_cs C Emean / F1mean / spread | 1.78 1.67 / 1.79 2.03 (col) / 0.63-0.65 | 0.53 0.66 / 0.34 1.11 / 4.2 6.8 5.0 4.7 | 1.84 1.41 / 1.72 1.64 / 1.00 0.87 0.73 0.68 |
**Both w1 and wc5 FAIL the no-regression** (latmid / latatm T collapse to 0-0.9); wc5 keeps the P_2 spread bounded
(unlike w1) but the order is gone.

## Beams (from 693c4fae7)
cand 0.5703 / 1.1846 / 0.0543 / 0.1656 (rF 0.77); w1 0.5703 / 1.2458 / 0.0723 / 0.0465 (rF 1.04); wc5 0.6261 (lum 0.81) /
1.2278 / 0.4116 / 0.0879 (rF 0.92).

## sp (pulse_u SP_FULL, one angular block per level, 453f879bd): CAND, geo3, geo4, col, lf0, lf2, wsx0, w1, wc5
running (the 512 level on one rank is slow); follow-up NOTE.
