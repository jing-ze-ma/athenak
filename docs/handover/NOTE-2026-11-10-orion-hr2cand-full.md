# NOTE 2026-11-10 ORION -> viper: hr2cand -- full (all but CAND pulse_cs R3 C, still running)

Follows the partial ee7227cd2; same binary 02798a1b2 and arms. All launched runs: no FATAL. sp table = the partial's plus
the retried CAND / CAND-nogj / CENN arms; cs table new.

## sp (p of all pairs)
| case mode var | PRODN | CAND | CAND-nogj | CENN |
|---|---|---|---|---|
| lat_ph_atm T E | 0.87 0.16 0.06 | 1.93 1.74 2.15 | 1.93 1.83 2.11 | 1.65 1.52 1.38 |
| lat_ph_atm T F1 | 1.05 -0.20 -0.54 | 1.83 1.11 2.53 | 1.78 1.32 2.45 | 1.89 1.68 1.50 |
| lat_ph_atm T F3 | 1.48 0.88 0.40 | 1.88 1.46 2.47 | 1.86 1.69 2.13 | 1.76 1.59 1.49 |
| lat_th_atm T E | 0.09 -0.18 -0.02 | 1.93 1.84 1.76 | 1.96 1.93 1.83 | 1.58 1.40 1.29 |
| lat_th_atm T F1 | -0.23 -0.25 -0.05 | 1.53 1.55 1.54 | 1.55 1.61 1.56 | 1.24 1.13 1.06 |
| lat_th_atm T F2 | 0.10 -0.10 -0.03 | 1.39 1.42 1.63 | 1.39 1.42 1.60 | 1.28 1.21 1.08 |
| lat_th_mid T E | 2.02 2.01 2.01 | 2.02 2.01 2.01 | 2.02 2.01 2.01 | 2.02 2.01 2.01 |
| lat_th_mid T F1 | 2.03 2.02 2.01 | 2.03 2.02 2.01 | 2.03 2.02 2.01 | 2.03 2.02 2.01 |
| lat_th_mid T F2 | 1.99 1.99 1.99 | 1.99 1.99 2.00 | 1.99 1.99 2.00 | 1.99 1.99 2.00 |
| pulse_u C E | 0.03 -0.77 | -1.25 -1.21 | -1.19 -1.15 | -2.56 -1.49 |
| pulse_u C F1 | 0.02 -0.81 | -1.26 -1.21 | -1.29 -1.04 | -2.48 -1.42 |
| pulse_u T E | 0.23 0.04 0.29 | 1.66 1.92 1.26 | 1.74 1.70 1.72 | 1.78 1.44 1.31 |
| pulse_u T F1 | 0.28 0.29 0.26 | 1.77 1.51 1.28 | 1.55 1.84 1.42 | 1.61 1.99 0.97 |

## cs (T: bin fields; C: SHELL means, R3 = N per face nx1/4, nx1 32..256; spread = max(Emax-Emin)/Emean per level)
| case mode var | PRODN | CAND | CAND-nogj | CENN |
|---|---|---|---|---|
| latatm_cs T Amode | - | 1.93 2.04 2.52 | 1.99 2.06 2.28 | 1.53 1.25 1.23 |
| latatm_cs T m1_e | - | 1.95 1.96 1.97 | 1.98 1.99 1.93 | 1.91 1.85 1.77 |
| latatm_cs T m1_f1 | - | 1.90 1.91 1.91 | 1.96 1.97 1.85 | 1.84 1.69 1.57 |
| latmid_cs T Amode | - | 1.93 1.96 1.98 | 1.93 1.96 1.98 | 1.93 1.96 1.98 |
| latmid_cs T m1_e | - | 2.06 2.04 2.02 | 2.06 2.04 2.02 | 2.06 2.04 2.02 |
| latmid_cs T m1_f1 | - | 1.97 2.01 2.01 | 1.97 2.01 2.01 | 1.97 2.01 2.01 |
| pulse_cs C Emean | 0.59 0.32 | - | 1.13 1.29 | 1.94 1.98 |
| pulse_cs C F1mean | 0.34 0.30 | - | 0.98 1.13 | 1.89 2.02 |
| pulse_cs C spread | 2.84e-01 4.18e-01 5.61e-01 6.86e-01 | - | 1.58e-01 6.23e-01 7.03e-01 5.91e-01 | 1.84e-01 4.12e-01 4.78e-01 4.99e-01 |
| pulse_cs T m1_e | - | 1.74 1.76 1.89 | 1.73 1.78 1.85 | 1.99 1.96 1.88 |
| pulse_cs T m1_f1 | - | 1.75 1.82 1.91 | 1.73 1.87 1.92 | 1.94 1.95 1.86 |
(PRODN cs T of the lateral twins and pulse_cs from the first pass: latatm_cs T Amode 0.02 -0.26 -0.07 / E 0.49 -0.06 -0.01;
latmid_cs 1.93 1.96 1.98 / 2.06 2.04 2.02; pulse_cs T E 0.25 0.09 0.22 / F1 0.30 0.25 0.26.)

## Reading
- **TIME: CAND restores it on every case where production is flat**: latatm_cs T E 0.49 -0.06 -0.01 -> **1.95 1.96 1.97**
  (F1 1.90 1.91 1.91, Amode 1.93 2.04 2.52); pulse_cs T 0.25 0.09 0.22 -> 1.74 1.76 1.89; lat_th_atm 0.09 -0.18 -0.02 ->
  1.93 1.84 1.76; lat_ph_atm 0.87 0.16 0.06 -> 1.93 1.74 2.15; pulse_u 0.23 0.04 0.29 -> 1.66 1.92 1.26. CAND-nogj ~ CAND
  in T. No regression on latmid_cs / lat_th_mid (~2.0 everywhere).
- **SPACE: pulse_cs R3 C: CENN (central radial faces) 1.94 1.98 / F1 1.89 2.02; CAND-nogj 1.13 1.29 / 0.98 1.13; PRODN
  0.59 0.32** (CAND R3 C pending) -> as on Raven: the half-range radial faces cost the space order; vgdh is aimed there.
- **pulse_u C (64..512, my harness = radial refinement only, 09-25 4x4 angular grid) diverges in every arm** (PRODN 0.03
  -0.77, CAND -1.25 -1.21, CENN -2.56 -1.49): not your "full refinement" -- your spfull decks are not on any branch, so for
  vgdh-geo I added an opt-in SP_FULL (nx2 = nx3 = max(4, nx1/16), 4x4 angular blocks) and rerun pulse_u that way there.
- The pulse_cs lateral spread grows with n in every arm incl. CENN (0.18 -> 0.50; CAND-nogj 0.16 -> 0.59-0.70).
