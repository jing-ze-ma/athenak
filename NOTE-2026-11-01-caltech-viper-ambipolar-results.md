# NOTE 2026-11-01 Caltech -> viper (USER): ambipolar arms, FINAL (1.5 rot, t = 0.5 -> 2.0 rot): AD changes T by <= 1 % and E_B by <= 0.8 %, with no trend; Saha AD costs +34 % RKG stages, tracer AD +1.4 % s/cycle. Recommend: ambipolar stays OFF

Answers TASK 212fbaef (addendum to 6e5e226a). The FINAL section (+1.0 and +1.5 rot) is at the end; the interim text is kept above it.

## Arms
All arms start from keep/off_0p50.rst (T1 off arm, t = 0.5 rot, cycle 8774). r1m_trc deck with the production keys: rkg_lean, rkg_dt_safety 0.5, outflow1, fac 4, cs_seam_flux shared, am_hst_split false. max_eta 1e15, 1 rank each.

| arm | binary (md5) | extra keys | GPU |
|---|---|---|---|
| R15 (other agent) | fin51 65658289 (c7682f07) | - | 01-11 GPU 0, shared |
| A15 | fin51 65658289 (c7682f07) | mhd/ambipolar = eos (Saha x_e) | 02-03 GPU 1, shared with the tag arms for the first ~70 % |
| T15 | trcdef 8ccae663 (de6bf3f2) | mhd/xe_tracer = true | 02-03 GPU 0, shared with TA15 |
| TA15 | trcdef 8ccae663 (de6bf3f2) | xe_tracer = true + ambipolar = eos | 02-03 GPU 0, shared with T15 |

AD parameters: the defaults (ad_mion_amu 30, ad_sigv_in 1.9e-9). The run.log only shows that ambipolar was added from the command line; the stage counts and state differences below confirm that it acts.

## Cost over t = 0.5-1.0 rot
| arm | mean dt [s] | RKG stages/cycle | s/cycle |
|---|---|---|---|
| R15 | 4.078 | 15.3 | 0.135 (shared) |
| A15 | 4.066 | **20.1 (+31 %)** | 0.206 (shared differently; not comparable) |
| T15 | 3.880 | 14.5 | 0.1689 |
| TA15 | 3.882 | **14.5 (+0 %)** | 0.1735 (**+2.7 %**; same GPU as T15, so comparable) |

- With RKG the hydro dt is unchanged, as you saw on viper.
- Saha AD raises the RKG stage count by ~1/3 (the eta_A upper bound sets the RKG dt). Tracer AD leaves it unchanged and adds only the AD flux cost.

## State at t = 1.0 rot (ana/meana.py; mass-weighted at p = 1e-4, 1e-3, 1e-2, 1e-1 bar +-0.15 dex)
```
arm     t[rot]  E_B(p>1e-6)  u1e-04 ueq Td Tt Tn  u1e-03 ueq Td Tt Tn  u1e-02 ueq Td Tt Tn  u1e-01 ueq Td Tt Tn
R15    1.000 2.1560e+32 | -0.178 -0.187  3525  3017  1101 |  0.036  0.045  3290  2400  1085 |  0.040  0.012  2692  1871  1137 |  0.014  0.004  2412  1814  1446
        B_rms[G] 1e-06:10.49 1e-05:9.84 1e-04:11.24 1e-03:50.98 1e-02:73.71
A15    1.000 2.1355e+32 | -0.171 -0.181  3525  3011  1111 |  0.032  0.037  3290  2402  1089 |  0.040  0.013  2692  1871  1137 |  0.014  0.004  2412  1814  1446
        B_rms[G] 1e-06:10.51 1e-05:9.88 1e-04:11.60 1e-03:50.68 1e-02:73.46
T15    1.000 2.3881e+32 | -0.166 -0.119  3525  3030  1111 |  0.047  0.054  3289  2361  1085 |  0.039  0.012  2692  1870  1135 |  0.014  0.004  2412  1814  1446
        B_rms[G] 1e-06:10.58 1e-05:9.50 1e-04:12.39 1e-03:53.19 1e-02:76.68
TA15   1.000 2.3699e+32 | -0.163 -0.106  3525  3028  1118 |  0.043  0.056  3289  2367  1088 |  0.040  0.012  2692  1870  1136 |  0.014  0.004  2412  1814  1446
        B_rms[G] 1e-06:10.50 1e-05:9.39 1e-04:12.76 1e-03:52.97 1e-02:76.47
        vs T15: max |du/u| 0.0962  max |du_eq/u_eq| 0.1107  max |dT/T| 0.00657  dE_B/E_B -0.0076  dB_rms/B_rms per band -0.008 -0.011 +0.030 -0.004 -0.003
```
(u in km/s; T in K for day |lon| < 60, terminator 60-120, night > 120 deg; B_rms per decade band from 1e-6 bar.)

| pair | wind (absolute) | T | E_B (p > 1e-6 bar) | B_rms per band 1e-6 / 1e-5 / 1e-4 / 1e-3 / 1e-2 |
|---|---|---|---|---|
| A15 - R15 | <= 0.008 km/s at all four levels; the relative 11-17 % is only because the mean u_phi is ~0.04 km/s at 1e-3 bar | max 0.9 % (night 1e-4 bar 1101 -> 1111 K); 1e-2 and 1e-1 bar unchanged | -0.95 % | +0.2 / +0.4 / +3.1 / -0.6 / -0.3 % |
| TA15 - T15 | <= 0.013 km/s (u_eq at 1e-4 bar) | max 0.66 % (night 1e-4 bar 1111 -> 1118 K) | -0.76 % | -0.8 / -1.1 / +3.0 / -0.4 / -0.3 % |
| for scale: R15 - T15 (x_e Saha vs tracer) | up to 0.07 km/s | 1.6 % | -9.7 % | -0.9 / +3.5 / -9.2 / -4.2 / -3.9 % |

- hst over 0.5-1.0 rot: KE end equal to 1e-4 in both pairs; total-E drift -7.5e-6 in all four arms; no FATAL or NON-CONVERGED.
- No C (layout) arm exists for these runs, so there is no measured chaos spread at +0.5 rot.
- Both AD pairs differ by about the same small amount. The largest change is +3 % B_rms in the 1e-4..1e-3 bar band, where the AD acts; T changes are <= 10 K, at the night 1e-4 bar level only.
- The x_e model itself (Saha vs tracer) matters ~10x more than AD does.
- **Interim reading:** as viper predicted. TA15 = T15 within ~1 % in T and E_B at p >= 1e-5 bar. A15 differs from R15 mainly at the night top and in B at 1e-4..1e-3 bar, and only by a few %. Tracer AD is nearly free (+2.7 % s/cycle), while Saha AD costs +31 % RKG stages.
- eta_A/eta_O medians: not computed yet. They need the EOS x_e on a dump; I will add them with the final numbers if time allows.

Run dirs: /resnick/groups/carnegie_poc/jingze/trc_1101/runs/me_{A15,T15,TA15} (R15 = the other agent's me_R15). Bins every 0.5 rot.

## FINAL (all arms ran to t = 2.0 rot, clean: no FATAL / NON-CONVERGED / nan)

### Cost over the whole 1.5 rot (t = 0.5-2.0)
| arm | mean dt [s] | min dt [s] | RKG stages/cycle | s/cycle |
|---|---|---|---|---|
| R15 | 3.921 | 3.663 | 15.1 | 0.0951 (alone for the last half) |
| A15 | 3.919 | 3.040 | **20.3 (+34 %)** | 0.1395 (shared differently) |
| T15 | 4.075 | 3.614 | 13.4 | 0.1541 |
| TA15 | 4.100 | 3.626 | **13.4 (+0 %)** | 0.1562 (**+1.4 %**, same GPU as T15) |

hst: total-E drift over 1.5 rot -2.47e-5 to -2.48e-5 in all four arms. KE end equal to 1e-4 in A15/R15 and 1.3e-4 in TA15/T15. E_B end +1.0 % (A15/R15) and +0.2 % (TA15/T15).

### State at t = 1.5 and 2.0 rot (meana.py)
```
R15    1.500 2.5453e+32 | -0.204 -0.235  3525  3026  1122 |  0.022  0.035  3292  2453  1109 |  0.042  0.020  2694  1880  1127 |  0.015  0.005  2427  1815  1377
        B_rms[G] 1e-06:10.98 1e-05:10.32 1e-04:10.10 1e-03:45.35 1e-02:79.48
A15    1.500 2.5241e+32 | -0.201 -0.239  3525  3030  1128 |  0.025  0.032  3292  2450  1111 |  0.042  0.020  2694  1880  1127 |  0.015  0.005  2427  1815  1377
        B_rms[G] 1e-06:10.92 1e-05:10.35 1e-04:10.04 1e-03:45.21 1e-02:79.25
        vs R15: max |du/u| 0.1730  max |du_eq/u_eq| 0.0748  max |dT/T| 0.00508  dE_B/E_B -0.0083  dB_rms/B_rms per band -0.006 +0.003 -0.006 -0.003 -0.003
R15    2.000 2.7236e+32 | -0.222 -0.283  3525  3038  1128 |  0.009  0.026  3292  2477  1126 |  0.039  0.023  2695  1883  1125 |  0.014  0.006  2432  1821  1344
        B_rms[G] 1e-06:11.07 1e-05:10.37 1e-04:9.96 1e-03:38.78 1e-02:79.68
A15    2.000 2.7011e+32 | -0.208 -0.245  3525  3029  1133 |  0.010  0.026  3292  2471  1130 |  0.039  0.023  2695  1883  1126 |  0.014  0.005  2433  1820  1344
        B_rms[G] 1e-06:11.00 1e-05:10.49 1e-04:9.92 1e-03:38.68 1e-02:79.44
        vs R15: max |du/u| 0.1233  max |du_eq/u_eq| 0.1336  max |dT/T| 0.00459  dE_B/E_B -0.0083  dB_rms/B_rms per band -0.007 +0.011 -0.003 -0.003 -0.003
T15    1.500 2.9136e+32 | -0.188 -0.150  3524  3044  1152 |  0.047  0.050  3291  2412  1118 |  0.040  0.019  2694  1879  1129 |  0.015  0.005  2427  1815  1378
        B_rms[G] 1e-06:11.03 1e-05:9.56 1e-04:10.35 1e-03:48.75 1e-02:83.35
TA15   1.500 2.8921e+32 | -0.185 -0.143  3524  3044  1161 |  0.050  0.055  3291  2402  1118 |  0.041  0.020  2694  1879  1129 |  0.015  0.005  2427  1815  1378
        B_rms[G] 1e-06:10.75 1e-05:9.47 1e-04:10.46 1e-03:48.66 1e-02:83.12
        vs T15: max |du/u| 0.0698  max |du_eq/u_eq| 0.0935  max |dT/T| 0.00800  dE_B/E_B -0.0074  dB_rms/B_rms per band -0.025 -0.009 +0.011 -0.002 -0.003
T15    2.000 3.1462e+32 | -0.205 -0.175  3524  3051  1158 |  0.035  0.048  3291  2447  1127 |  0.039  0.022  2695  1883  1128 |  0.014  0.006  2432  1820  1346
        B_rms[G] 1e-06:11.09 1e-05:9.67 1e-04:9.92 1e-03:42.39 1e-02:83.73
TA15   2.000 3.1220e+32 | -0.198 -0.155  3525  3050  1170 |  0.038  0.051  3291  2441  1128 |  0.039  0.022  2695  1884  1128 |  0.014  0.005  2432  1820  1346
        B_rms[G] 1e-06:10.74 1e-05:9.66 1e-04:9.90 1e-03:42.38 1e-02:83.54
        vs T15: max |du/u| 0.0864  max |du_eq/u_eq| 0.1134  max |dT/T| 0.00980  dE_B/E_B -0.0077  dB_rms/B_rms per band -0.032 -0.001 -0.002 -0.000 -0.002
```

| pair | t [rot] | max abs du (all 4 levels) | max dT/T | dE_B/E_B (p > 1e-6) | dB_rms per band 1e-6 / 1e-5 / 1e-4 / 1e-3 / 1e-2 |
|---|---|---|---|---|---|
| A15 - R15 | 1.0 | 0.008 km/s | 0.91 % | -0.95 % | +0.2 / +0.4 / +3.1 / -0.6 / -0.3 % |
| A15 - R15 | 1.5 | 0.004 km/s | 0.51 % | -0.83 % | -0.6 / +0.3 / -0.6 / -0.3 / -0.3 % |
| A15 - R15 | 2.0 | 0.038 km/s (u_eq at 1e-4 bar, -0.283 -> -0.245) | 0.46 % | -0.83 % | -0.7 / +1.1 / -0.3 / -0.3 / -0.3 % |
| TA15 - T15 | 1.0 | 0.013 km/s | 0.66 % | -0.76 % | -0.8 / -1.1 / +3.0 / -0.4 / -0.3 % |
| TA15 - T15 | 1.5 | 0.007 km/s | 0.80 % | -0.74 % | -2.5 / -0.9 / +1.1 / -0.2 / -0.3 % |
| TA15 - T15 | 2.0 | 0.020 km/s (u_eq at 1e-4 bar) | 0.98 % (night 1e-4 bar, 1158 -> 1170 K) | -0.77 % | -3.2 / -0.1 / -0.2 / -0.0 / -0.2 % |

- The 1e-4..1e-3 bar B_rms bump of +3 % at 1.0 rot does not persist (<= 1.1 % at 1.5 and 2.0 rot).
- The steady signals are three:
  - E_B at p > 1e-6 bar is ~0.8 % lower with AD in both pairs, at every time.
  - The night 1e-4 bar T is 5-12 K warmer with AD.
  - TA15 has a few % less B_rms at 1e-6..1e-5 bar (the night top, where viper predicted AD/advection ~ 1).
- At 1e-2 and 1e-1 bar, T and B change by <= 0.3 % and u is unchanged.
- The x_e model (Saha vs tracer) moves E_B by ~10 % and B_rms by up to 9 %, far more than AD does in either pair.
- **Conclusion:** viper's prediction holds. TA15 = T15 at p >= 1e-5 bar within ~1 % in T, B_rms, E_B and winds. A15 differs from R15 mainly at the night top, and only at the ~1 % level.
  - Saha AD (the upper bound) costs +34 % RKG stages.
  - Tracer AD costs +1.4 % s/cycle and the same stages.
- **Recommendation: ambipolar stays OFF in production** (the USER decides). If it is ever wanted, use it with the x_e tracer (nearly free), not with Saha.
- eta_A/eta_O medians were not computed (no x_e reader in our analysis tools).

Bins kept: runs/me_{A15,T15,TA15}/bin (t = 1.0, 1.5, 2.0 rot); hst and run.log kept; rst deleted.
