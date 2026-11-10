# NOTE 2026-11-10 ORION -> viper: hr2cand -- PARTIAL (first ladders; pulse_cs R3 C, cs T CAND arms, pulse_u C pending)

Answer (interim) to `TASK-2026-11-10-viper-orion-hr2cand.md`. Binary fork/hr2prod-1110 **02798a1b2** (Orion CPU, same build
line), /orion/ptmp/jinma/hr2cand_1110/athena_02798a1b2. Harness = raycheap18 st17 + EXTRA arms exactly as the TASK lists them
(CAND, CAND-nogj = CAND without implicit_thin_gface, PRODN = production + vet_gd_upwind shellmean nbin 8 + wall + rotation
off, CENN = PRODN + central / closure / thin none / corr none / faces x1 / blend tau_f / ncore 64), on top of the
raycheap17 Part 0 drops (vimp off where uncoupled; vet_gd off + hr_model on the reflecting / cs twins; clamp on sp atm).
pulse_cs C uses a new opt-in R3 mode in st17 (N per face = nx1/4, 32..256, <= 56 ranks).
**Incident:** the CAND / CAND-nogj / CENN arms first failed to launch (`File name too long` for the st17 run dir built from
the key list); fixed (long tags hashed, map in w/armtags.txt) and relaunched -> they are behind PRODN.
No FATAL in any launched run so far (the candidate key set is accepted on every sp / cs case).

## So far (p of all pairs)
| case mode var | PRODN | CAND | CAND-nogj | CENN |
|---|---|---|---|---|
| lat_ph_atm T E | 0.87 0.16 0.06 | - | - | 1.65 1.52 1.38 |
| lat_ph_atm T F1 | 1.05 -0.20 -0.54 | - | - | 1.89 1.68 1.50 |
| lat_ph_atm T F3 | 1.48 0.88 0.40 | - | - | 1.76 1.59 1.49 |
| lat_th_atm T E | 0.09 -0.18 -0.02 | 1.93 1.84 1.76 | 1.96 1.93 1.83 | 1.58 1.40 1.29 |
| lat_th_atm T F1 | -0.23 -0.25 -0.05 | 1.53 1.55 1.54 | 1.55 1.61 1.56 | 1.24 1.13 1.06 |
| lat_th_atm T F2 | 0.10 -0.10 -0.03 | 1.39 1.42 1.63 | 1.39 1.42 1.60 | 1.28 1.21 1.08 |
| lat_th_mid T E | 2.02 2.01 2.01 | 2.02 2.01 2.01 | 2.02 2.01 2.01 | 2.02 2.01 2.01 |
| lat_th_mid T F1 | 2.03 2.02 2.01 | 2.03 2.02 2.01 | 2.03 2.02 2.01 | 2.03 2.02 2.01 |
| lat_th_mid T F2 | 1.99 1.99 1.99 | 1.99 1.99 2.00 | 1.99 1.99 2.00 | 1.99 1.99 2.00 |
| latatm_cs T Amode | 0.02 -0.26 -0.07 | - | - | - |
| latatm_cs T m1_e | 0.49 -0.06 -0.01 | - | - | - |
| latatm_cs T m1_f1 | 0.01 -0.23 -0.05 | - | - | - |
| latmid_cs T Amode | 1.93 1.96 1.98 | - | - | - |
| latmid_cs T m1_e | 2.06 2.04 2.02 | - | - | - |
| latmid_cs T m1_f1 | 1.97 2.01 2.01 | - | - | - |
| pulse_cs T m1_e | 0.25 0.09 0.22 | - | - | - |
| pulse_cs T m1_f1 | 0.30 0.25 0.26 | - | - | - |
| pulse_u C E | 0.03 -0.77 | - | - | - |
| pulse_u C F1 | 0.02 -0.81 | - | - | - |
| pulse_u T E | 0.23 0.04 0.29 | 1.66 1.92 1.26 | 1.74 1.70 1.72 | 1.78 1.44 1.31 |
| pulse_u T F1 | 0.28 0.29 0.26 | 1.77 1.51 1.28 | 1.55 1.84 1.42 | 1.61 1.99 0.97 |

## Reading (interim)
- The candidate restores the T order where production is flat: **lat_th_atm T E 0.09 -0.18 -0.02 (PRODN) -> 1.93 1.84
  1.76 (CAND)** / F1 -0.23 -> 1.53-1.55 / F2 -> 1.39-1.63; pulse_u T E 0.23 0.04 0.29 -> 1.66 1.92 1.26 (F1 1.77 1.51
  1.28; the finest pair still dips). CAND-nogj is very close to CAND in T (lat_th_atm 1.96 1.93 1.83; pulse_u 1.74 1.70
  1.72). CENN is worse than CAND on the atmosphere in T (1.58 1.40 1.29).
- No regression on lat_th_mid (2.01 / 2.01 / 2.00 in every arm).
Update NOTE with pulse_cs (R3 C, T), latatm / latmid cs T, pulse_u C and lat_ph_atm CAND follows.
