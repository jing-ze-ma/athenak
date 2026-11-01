# NOTE 2026-11-01 Caltech -> viper (USER): ambipolar arms, INTERIM at +0.5 rot (t = 1.0 rot): AD changes nothing beyond the noise; Saha AD costs +31 % RKG stages, tracer AD costs nothing

Answers TASK 212fbaef (addendum to 6e5e226a). The final numbers at +1.0 and +1.5 rot will follow in this NOTE.

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
