# NOTE DeltaAI -> Raven (2026-10-18 16:25 UTC): hg145 deep-convection saturation check, DeltaAI 0-28 d (TASK 2c703c2a)

## Verdict
1. **KE growth is slowing but NOT saturated at 28 d.** The e-fold time rises 1.7 d (2-4 d) -> 3.0 d (10-18 d) -> 4.7 d (24-26 d) -> 5.6 d (26-28 d).
   KE_int is 1.38e42 erg at 28 d (wedge): inside your 1-2e42 MLT estimate and still rising, so it should reach the low 2e42s, not level at 1e42.
2. **Resolved convection carries only part of the non-diffusive flux.** (L_conv,res + L_adv)/(L_in - L_diff) at 28 d: 20-30 Rsun 0.08,
   **30-60 Rsun 0.51** (0.10 at 24 d, 0.25 at 26 d, rising fast), 60-100 Rsun 0.00. The adaptive sub-grid closure still carries the rest
   (L_MLT/L_in from the hst falls 1.00 -> 0.78 over 16-28 d as resolved transport grows). Not at 1 yet anywhere.
3. **Flags:** KE < 3e42 OK. Band rms v_r max 13.2 km/s (19-30 Rsun) < 20 OK. **But max|v_r| ~40 km/s in EVERY band at 28 d**, incl. 8-19 Rsun
   (rms only 1.1 there) and 60-100 Rsun (rms 1.1; jumped 3.6 -> 26 km/s between 20 and 22 d). These are narrow plumes and overshoot well beyond
   the Fe-bump zone. Worth a look in the slices: the 28 d slices on the page show them rising from the 8-30 Rsun region.

## Method (scripts + full output: docs/handover/hg145-sat-1018/)
- Shell means are solid-angle weighted as in reduce.py; Omega = pi/3; L_in from hst col 7 (wedge), last occurrence per time (restart rows dropped).
- L_conv,res = Omega r^2 (Fenth + Fke) from red_*.npz (Fenth uses the ideal-gas h = 5/3 e/rho, as reduce.py).
- **Important:** m1_f1 is the LAB-frame flux. It already contains the radiation enthalpy the flow carries, 4/3 E v_r, and that is the main
  convective channel in this radiation-dominated zone. So "L_rad from the m1 flux" in the TASK would count most convective transport as
  radiative. sat2.py splits it: L_adv = Omega r^2 4/3 (<E v_r> - <E><v_r>) (from the bins), L_diff = Omega r^2 <F1> - L_adv.
  Plain TASK ratio L_conv,res/(L_in - L_rad,lab) is in sat.txt: ~0 until 16 d, then -0.18 (20-30), +0.09 (30-60), 0 (60-100) at 28 d.
- Ideal-gas Fenth is NEGATIVE in 20-30 Rsun (-0.055 L_in at 28 d), so with L_adv +0.087 the net there is +0.032 L_in. Ionisation enthalpy
  is missing from h, so treat the 20-30 Rsun gas-enthalpy number as uncertain.

## Table (28 d last 3 dumps and every 2 d; L in units of L_in)
| t [d] | L_adv 20-30/30-60 | L_conv,res+L_adv 20-30/30-60 | L_diff 20-30/30-60/60-100 | ratio 20-30/30-60 | rms v_r / max|v_r| km/s 8-19, 19-30, 30-60, 60-100 | KE_int |
|---|---|---|---|---|---|---|
| 0 | 0 / 0 | 0 / 0 | 0.53 / 0.42 / 0.70 | 0 / 0 | 0/0, 0/0, 0/0, 0/0 | 0 |
| 8 | 0.000 / 0.000 | 0.001 / 0.000 | 0.53 / 0.42 / 0.70 | 0.001 / 0.000 | 0.02/1.2, 0.81/2.3, 0.58/2.0, 0.11/0.5 | 3.8e39 |
| 14 | 0.004 / 0.002 | 0.004 / 0.002 | 0.53 / 0.42 / 0.70 | 0.008 / 0.003 | 0.05/4.1, 2.28/7.3, 1.62/5.9, 0.20/0.9 | 2.9e40 |
| 18 | 0.015 / 0.006 | 0.014 / 0.006 | 0.48 / 0.47 / 0.71 | 0.027 / 0.011 | 0.15/13.3, 4.54/15.6, 3.10/10.9, 0.27/1.3 | 1.1e41 |
| 20 | 0.027 / 0.010 | 0.024 / 0.011 | 0.62 / 0.52 / 0.71 | 0.063 / 0.023 | 0.28/21.7, 6.24/22.1, 4.22/14.0, 0.33/3.6 | 2.1e41 |
| 22 | 0.043 / 0.018 | 0.034 / 0.019 | 0.58 / 0.59 / 0.72 | 0.080 / 0.048 | 0.49/29.3, 8.16/29.0, 5.63/17.6, 0.42/26.5 | 3.8e41 |
| 24 | 0.060 / 0.029 | 0.041 / 0.032 | 0.63 / 0.69 / 0.75 | 0.112 / 0.104 | 0.75/35.1, 10.0/33.7, 7.32/21.2, 0.56/33.4 | 6.3e41 |
| 26 | 0.076 / 0.044 | 0.043 / 0.049 | 0.60 / 0.81 / 0.77 | 0.107 / 0.253 | 0.99/40.9, 11.7/36.4, 9.18/24.8, 0.79/37.5 | 9.7e41 |
| 27.5 | 0.086 / 0.057 | 0.037 / 0.063 | 0.61 / 0.86 / 0.81 | 0.093 / 0.444 | 1.06/39.1, 12.9/38.0, 10.6/27.5, 1.02/39.7 | 1.27e42 |
| 28.0 | 0.087 / 0.061 | 0.032 / 0.068 | 0.60 / 0.87 / 0.83 | 0.081 / 0.507 | 1.06/39.8, 13.2/39.5, 11.1/28.4, 1.10/40.0 | 1.38e42 |
Note L_diff in 30-60 Rsun rose 0.42 -> 0.87 L_in from 16 to 28 d: the diffusive flux itself took over most of the load there as the
structure adjusted, so the gap left for convection shrank to 0.13 L_in. 60-100 Rsun has no resolved transport. No GPU used.
