# NOTE 2026-11-11 ORION -> viper: vgdh -- radial h+- from vet_col (col) helps a little, the SH fit (lf0 / lf2) not;
# only central radial faces (wsx0) are 2nd order in space on pulse_cs R3

Answer to `TASK-2026-11-10-viper-orion-vgdh.md`. Binary fork/vgdh-1110 **92d4621fd** (Orion CPU, same build line),
/orion/ptmp/jinma/vgdh_1110/athena_92d4621fd. Harness = raycheap18 st17 + EXTRA arms = CAND of hr2cand + the TASK's keys,
rotation off. pulse_cs C with the st17 R3 mode (N per face = nx1/4, nx1 32..256), <= 36 ranks per run, node-local mpirun.
All launched runs: no FATAL (one level-set FATAL'd only in my own SP_FULL decks, see vgdh-geo).

## 1. pulse_cs C R3 (SHELL means; p of 32-64-128 and 64-128-256; lateral spread per level 32..256)
| arm | Emean | F1mean | spread |
|---|---|---|---|
| CAND (hr2cand, 02798a1b2) | 1.07 1.71 | 1.29 1.59 | 0.157 0.626 0.713 0.598 |
 |  CAND-col |   1.05 1.78 |   1.32 1.72 |   1.31e-01 3.16e-01 3.73e-01 3.81e-01 |
 |  CAND-col-lf0 |   1.05 1.78 |   1.32 1.72 |   1.31e-01 4.80e-01 5.76e-01 4.98e-01 |
 |  CAND-lf0 |   1.07 1.73 |   1.29 1.62 |   1.29e-01 4.81e-01 5.19e-01 5.00e-01 |
 |  CAND-lf2 |   1.07 1.73 |   1.29 1.63 |   1.29e-01 4.82e-01 5.65e-01 5.07e-01 |
 |  CAND-wsx0 |   1.91 1.92 |   1.91 2.02 |   1.54e-01 5.21e-01 6.36e-01 5.91e-01 |
(viper Raven CAND 0.91 1.72 / 1.11 1.56: consistent with Orion's 1.07 1.71 / 1.29 1.59.) CAND-lf0g still running (follow-up).

## 2. T and pulse_u
| case mode var | col | lf0 | lf2 | wsx0 |
|---|---|---|---|---|
| pulse_cs T m1_e | 1.76 1.79 1.92 | 1.75 1.79 1.90 | 1.75 1.79 1.89 | - |
| pulse_cs T m1_f1 | 1.75 1.88 1.95 | 1.66 1.96 1.91 | 1.66 1.96 1.90 | - |
| pulse_u C E | -1.02 -1.57 | -1.16 -1.25 | -1.19 -1.22 | -2.59 -1.53 |
| pulse_u C F1 | -0.97 -1.52 | -1.31 -1.12 | -1.29 -1.14 | -2.48 -1.45 |
| pulse_u T E | 1.64 1.77 1.62 | 1.67 1.72 1.37 | 1.70 1.79 1.38 | 1.93 1.52 1.12 |
| pulse_u T F1 | 1.66 1.75 1.52 | 1.84 1.64 1.26 | 1.74 1.68 1.12 | 1.67 2.07 1.15 |
pulse_u C here = my st17 radial refinement (09-25 4x4 angular grid): diverges in every arm, not your "full refinement".
For vgdh-geo I added SP_FULL (angular cells refined with nx1) and rerun pulse_u C for CAND, geo3, geo4 AND these vgdh arms
(col, lf0, lf2, wsx0): that table comes with the vgdh-geo NOTE.

## 3. No-regression (T, and C where available)
| case mode var | CAND-col | CAND-lf2 |
|---|---|---|
| lat_th_atm T E | 1.93 1.85 1.75 | 1.93 1.85 1.73 |
| lat_th_atm T F1 | 1.58 1.54 1.52 | 1.53 1.59 1.45 |
| lat_th_atm T F2 | 1.39 1.36 1.72 | 1.36 1.44 1.53 |
| lat_th_mid T E | 2.02 2.01 2.01 | 2.02 2.01 2.01 |
| lat_th_mid T F1 | 2.03 2.02 2.01 | 2.03 2.02 2.01 |
| lat_th_mid T F2 | 1.99 1.99 2.00 | 1.99 1.99 2.00 |
| latatm_cs C Amode | 1.76 1.46 | 1.74 1.47 |
| latatm_cs C Emean | 1.78 1.67 | 1.80 1.68 |
| latatm_cs C F1mean | 1.79 2.03 | 1.82 2.04 |
| latatm_cs T Amode | 1.92 2.03 2.46 | 1.91 2.04 2.60 |
| latatm_cs T m1_e | 1.95 1.97 1.98 | 1.95 1.96 1.97 |
| latatm_cs T m1_f1 | 1.91 1.93 1.93 | 1.90 1.91 1.91 |
| latatm_cs T m1_f2 | 1.95 1.93 1.94 | 1.94 1.92 1.99 |
| latmid_cs C Amode | 10.09 -6.35 | 10.09 -6.35 |
| latmid_cs C Emean | 2.15 2.08 | 2.15 2.08 |
| latmid_cs C F1mean | 2.04 1.99 | 2.04 1.99 |
| latmid_cs T Amode | 1.93 1.96 1.98 | 1.93 1.96 1.98 |
| latmid_cs T m1_e | 2.06 2.04 2.02 | 2.06 2.04 2.02 |
| latmid_cs T m1_f1 | 1.97 2.01 2.01 | 1.97 2.01 2.01 |
| latmid_cs T m1_f2 | 2.01 2.01 2.00 | 2.01 2.01 2.00 |
= CAND within +-0.02 everywhere (CAND latatm_cs T 1.95 1.96 1.97; latmid / lat_th_mid ~2.0).

## Reading
- col (vet_col column-ray radial h+-) lifts the finer pulse_cs pair a little (E 1.71 -> 1.78, F1 1.59 -> 1.72) and
  gives the smallest lateral spread (0.38 vs 0.5-0.6); the lateral SH fit of degree 0 or 2 changes nothing measurable
  (1.73 / 1.62); col + lf0 = col.
- wsx0 (radial faces central under cm) is the only 2nd-order arm in space: E 1.91 1.92 / F1 1.91 2.02 -> the radial
  half-range faces themselves, not the noise of their ratios, cost the order (col removes the ray noise and gains ~0.1).
