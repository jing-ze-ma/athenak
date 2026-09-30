# synth2_emis: angle-dependent, R1000 emission post-processor for WASP-121b (1x / 10x, rot 300)

State: `w1x/bin/dhj.hydro_w.00150.bin` and `w10x/bin/dhj.hydro_w.00150.bin`, t = 3.3046e7 s =
**rot 300.00** (read in place). T and p come from the runs' own EOS dumps (`synth_rot300*/diag/eos_table_*.txt`).
The old-tool comparison uses `olr_00600.txt` (rot 300.00 + 380 s, 22 cycles later).
Numbers are in `out/products.txt`, `out/hotspot.txt`, `out/validate_1d.txt` and `out/gcm2s_*.txt`.
Scripts are here and copied to `docs/handover/scripts/w121/synth2_emis/`.

## Method
- **Radiative transfer (`rtcore.py`).**
  - Geometry: one column per cell column, along the local vertical (plane-parallel slant path). The emission angle is mu = n.o, with the sub-observer longitude 180 - 360 phi, as in synth.py. There are 72 phases.
  - Physics: absorption only, LTE source B(T). The source is linear in tau inside each cell. Each layer's emission is weighted by r^2 of its own radius. This is the column form of "intensity is conserved along a ray", and it keeps the planet's luminosity L = 4 pi r_ph^2 F_ph.
  - Disk flux: F_p d^2 = sum over columns of J_top(mu) mu dOmega.
  - Not included: no scattering, so no reflected light and no Rayleigh scattering.
- **Opacity (`prt_table.py`).**
  - Code and data: petitRADTRANS 3.4.0, correlated-k at R = 1000, 0.3-28 um, 16 g-points, random-overlap mixing (pRT's own). Species: H2O POKAZATEL, CO HITEMP, CO2 UCL-4000, CH4 HITEMP, OH MoLLIST, TiO McKemmish, VO VOMYT, FeH MoLLIST, Fe Kurucz, Na and K Allard, SiO SiOUVenIR, CIA H2-H2 and H2-He, and H- bound-free/free-free.
  - Installation: pRT is in a venv, `/viper/ptmp2/jinma/prt_venv` (python-waterboa 3.13). The system python 3.9 has no headers, so pRT could not be built there.
  - Data: under `/viper/ptmp2/jinma/prt_data`, fetched from the Keeper Seafile API with `keeper.py`, because pRT's own fetcher needs Chrome.
  - Chemistry: FastChem 3 (pyfastchem 4.0.3) with ions (e-, H-), thermal dissociation and local equilibrium condensation (no rainout). Abundances are the runs' own inputs: Asplund09 for 1x, and metals +1 dex with H and He fixed for 10x (`ckhitemp_0925/fastchem_input`, `wasp121_0925/met10/fc_in`).
  - Premixed (T, p) table: 48 x 46 points, 150-6000 K, 1e-12 to 1e3 bar.
  - Range limits: pRT clamps T above each line list's range (2995-4000 K). Below 1e-3 bar the dayside reaches 3300-5000 K, so line opacity there is evaluated at the clamp.
- **Star.** PHOENIX-ACES-AGSS-COND-2011 HiRes, trilinear interpolation at T_eff 6628 K, log g 4.24, [Fe/H] +0.13. Its integral is 0.985 of sigma T_eff^4. Beyond 5.5 um a scaled blackbody is used. R* = 1.461 Rsun, Rp = 1.742 RJ, a/R* = 3.7844 (Sing+2024).
- **Throughputs.**

  | Band | Throughput source | Wavelength window |
  |---|---|---|
  | NRS1, NRS2 | SVO JWST/NIRSpec.G395H_F290LP (starts at 2.87 um) | 2.70-3.72 and 3.82-5.15 um (ME23) |
  | NIRISS SOSS o1 / o2 | CRDS jwst_niriss_spectrace_0023.fits | 0.85-2.85 um / 0.6-0.85 um |
  | WFC3 | SVO HST/WFC3_IR.G141 | 1.12-1.64 um |
  | IRAC1, IRAC2 | SVO Spitzer/IRAC.I1 and .I2 | full response |

  All bands are photon-weighted: Fp/Fs = sum(P R lambda dlambda) / (R*^2 sum(F* R lambda dlambda)).
- **Products.**
  - Day and night flux: phi = 0.5 and phi = 0.
  - Amplitude: max - min of a 2-harmonic fit.
  - Offsets: off2 is the peak of the 2-harmonic fit; off1 is the peak of the first harmonic alone (+ = east).
  - T_b: from pi <B>_w(T_b) = y <F*>_w (R*/Rp)^2, using the same star as the curve.
  - T_day and T_night: from the bolometric disk flux (0.3-28 um plus a Planck tail, 0.05-0.23 %) as (P/(Rp^2 sigma))^(1/4).
  - epsilon = 8/(3 (T_d/T_n)^4 + 5), and A_B = 1 - 1.5 (T_d^4 + 5/3 T_n^4)/T0^4 with T0 = T* (R*/a)^(1/2) (Cowan & Agol 2011).

## Validation
**(a) Consistency with the GCM (`emis11.py`, `gcm2s.py`).** I ran the same columns through the run's own 11-band hiT2 Exo-FMS opacity (ck11.py = rr_lib_met transcription). Results against the olr rt_Fb:

| Run | Ray RT with r^2 (this tool) | Plane-parallel | GCM's own spherical two-stream closure, re-implemented |
|---|---|---|---|
| 1x | L = 2.233e30 vs 2.794e30 erg/s (**0.799**; the same 0.78-0.86 in every band, day 0.80, night 0.79) | 1.38 | L = 0.960 (per column median 0.946, 5-95 % 0.91-1.00) |
| 10x | L = 2.492e30 vs 2.858e30 (**0.872**) | 1.31 | L = 1.023 (bands 0.98-1.06) |

- The GCM closure is ck_spherical (two_stream_rt.hpp:1002): S = (I+ + I-)/2 is passed continuously at faces, so the mean intensity does not dilute through the transparent upper shell. The down ray is carried negative instead.
- For a thin shell this gives L_GCM / L_true = 2x/(1+x), with x = A_top/A_ph = (r_top/r_ph)^2. For the 1x state, x = 1.72 and the formula gives 1.26; the measured ratio is 1/0.80 = 1.25.
- The gap is therefore a GCM thermal-RT bias, not a post-processing error: the run's rt_Fb, hst Lir and old synth numbers are ~25 % (1x) / ~15 % (10x) above what its T field actually radiates. The remaining 4 % / 2 % comes from cell discretisation (half-cell steps in the GCM vs linear source here).

**(b) 1-D checks (`validate_1d.py`).**
- An isothermal column (1500 K and 2500 K) gives I/B = 1.000000 in all 4536 bins at mu = 0.05-1.
- Against pRT's own `calculate_flux` on the same column (hydrostatic, same angle grid):

  | Column | 0.3-28 um | 0.85-5.15 um bands | 0.6-0.85 um |
  |---|---|---|---|
  | Dayside | 1.033 | 1.003-1.015 | 1.137 |
  | Nightside | 0.996 | 0.991-0.997 | 1.003 |

  The dayside 0.6-0.85 um excess comes from interpolating the premixed table in T across the steep TiO/VO/H- region of the 3000-5000 K upper atmosphere.

**(c) Old to new on 1x, one change at a time.** Figures are NRS1 day/night and NRS2 day/night in ppm; T_day/T_night in K.

| Step | NRS1 d/n | NRS2 d/n | T_day/T_night |
|---|---|---|---|
| OLD (GCM rt_Fb, Lambertian, 11 bands, BB star) | 3713/845 | 4257/1291 | 2654/1694 |
| A1: ray RT on the GCM 11-band opacity (**GCM closure bias**) | 2985/667 | 3347/996 | 2512/1593 |
| A2: + angle-dependent intensity (vs Lambertian) | 2986/700 | 3341/1024 | 2499/1617 |
| A3b: + pRT opacity and FastChem, binned to the 11 bands | 3028/720 | 3374/1077 | 2492/1622 |
| A3: + full R1000 resolution (top-hat windows) | 2963/693 | 3457/1129 | same |
| A3t: + real throughputs, photon-weighted | 3012/753 | 3432/1126 | same |
| **NEW**: + PHOENIX star (F*_PHOENIX/F*_BB = 0.936 NRS1, 0.896 NRS2, 1.02 SOSS1, 1.05 WFC3) | 3219/805 | 3830/1257 | 2488/1622 |
| NEWc: chemistry at min(T, 2000 K), i.e. no dissociation or ionisation | 4711/900 | 4761/1317 | 2872/1681 |

- Offsets change by less than 1 deg at every step.
- At full resolution, Lambertian (NEWl) vs angle-dependent (NEW) differs by 34 ppm at NRS1 night and 29 ppm at NRS2 night; the day flux changes by less than 10 ppm.
- The NRS2 dayside deficit (new 3830 vs 4924 ppm observed) does not come from the coarse bands. It is -1094 ppm, larger than the old -667 ppm once the closure bias is removed, even though the PHOENIX star adds +400 ppm.
- Opacity inconsistency: pRT line lists and FastChem differ from the GCM's Exo-FMS premixed table. This changes 11-band luminosity by -1.1 % on 1x (A2 vs A3b) and NRS night flux by +3 to +5 %. It is a post-processing inconsistency, stated rather than hidden.

## Results table
Photon-weighted, PHOENIX star, 72 phases. The obs column gives the observed value and its source (LITERATURE.md).

| Quantity | old 1x | new 1x | old 10x | new 10x | obs (source) |
|---|---|---|---|---|---|
| NRS1 day / night [ppm] | 3713 / 845 | 3219 / 805 | 3716 / 278 | 3479 / 326 | 3924 / 136 (ME23) |
| NRS2 day / night [ppm] | 4257 / 1291 | 3830 / 1257 | 4283 / 539 | 4261 / 632 | 4924 / 630 (ME23) |
| NRS1 / NRS2 offset (2-harm; 1st harm) [deg E] | 17.2 / 18.1 | 17.8; 20.3 / 18.6; 20.7 | 5.9 / 6.7 | 6.4; 7.2 / 7.0; 7.6 | 3.36 / 2.66 (ME23) |
| NRS1 / NRS2 amplitude [ppm] | 3052 / 3179 | 2557 / 2732 | 3462 / 3784 | 3140 / 3624 | -- |
| NRS1 T_b day / night [K] | 2851 / 1574 | 2528 / 1471 | 2852 / 1150 | 2624 / 1136 | 2762 / 926 (ME23, PHOENIX) |
| NRS2 T_b day / night [K] | 2785 / 1544 | 2469 / 1465 | 2795 / 1122 | 2620 / 1144 | 2768 / 1122 (ME23) |
| SOSS o1 day / night [ppm], offset [deg] | 850 / 143, 20.0 | 681 / 127, 20.9 | 963 / 56, 7.9 | 829 / 51, 8.9 | 1156, 5.1 +- 1.4 E (Splinter+25) |
| SOSS o2 day [ppm] | 253 | 187 | 407 | 331 | 363 (Splinter+25) |
| WFC3 day / night [ppm], offset [deg] | 859 / 129, 20.4 | 642 / 111, 21.2 | 998 / 52, 8.1 | 782 / 43, 9.2 | c1 ~1160-1195; ~6 E (ME22) |
| IRAC1 day / night [ppm] | 3945 / 1051 | 3208 / 943 | 3969 / 392 | 3447 / 460 | 4077 / 553 (Davenport+25) |
| IRAC2 day / night [ppm] | 4294 / 1299 | 4047 / 1245 | 4317 / 542 | 4534 / 558 | 5121 / 1045 (Davenport+25) |
| T_day / T_night (bolometric) [K] | 2654 / 1694 | 2488 / 1622 | 2789 / 1337 | 2667 / 1300 | 2717 +- 17 / 1562 +- 19 (Splinter+25) |
| epsilon | 0.347 | 0.370 | 0.129 | 0.138 | 0.246 +- 0.014 (Splinter+25) |
| A_B implied by T_d, T_n | 0.295 | 0.445 | 0.267 | 0.384 | 0.277 (Splinter+25; the run imposes 0.277) |
| disk-averaged 4 pi <P> [erg/s] | 2.96e30 | 2.33e30 | 3.05e30 | 2.60e30 | run absorbs L_sw ~2.77e30 |

- "Old" is synth.py recomputed on olr_00600. It matches the old RESULTS.md, for example NRS1 3713/845 and 10x 3716/278.
- Scattering and reflected light are absent, so SOSS o2 and WFC3 are lower bounds.

## Masked hot knot (1x, coordinator 09-30)
- **Mask** (`hotspot.py`): columns with |lat| < 20 and 95 < lon < 135 deg (lat = asin z, lon = atan2(-y, -x), dhjcs). That is 222 columns.
- **Replacement:** in those columns, at p < 1e-3 bar, T is replaced by the mean of the reference columns (140 < lon < 170, |dlat| < 3 deg) interpolated in log p. p is kept; rho is rescaled at fixed p by mu_new T_old/(mu_old T_new); chemistry and opacity follow from the new (T, p) through the table.
- **The mask removes the knot** (`out/hotspot_T1e-4.png`): T at 1e-4 bar, maximum over the box in 4-deg bins, goes from 3171 K to 2196 K. The warmer equatorial tongue east of lon 135 deg is outside the box and stays.
- **Method:** only the masked columns are re-solved. The rest of the disk sum is identical.

| Band | max dFp/Fs [ppm] at phi | ampl [ppm] / d | off2 / d [deg] | off1 / d [deg] | knot-column contribution p16/p50/p84 [bar] | fraction from p < 1e-3 bar |
|---|---|---|---|---|---|---|
| NRS1 | -36.0 at 0.194 | 2557 / +1.3 | 17.80 / -0.70 | 20.28 / -0.92 | 8e-4 / 2.6e-3 / 3.6e-2 | 0.21 |
| NRS2 | -28.2 at 0.194 | 2732 / +0.8 | 18.55 / -0.50 | 20.71 / -0.66 | 8e-4 / 6.6e-3 / 7.0e-2 | 0.17 |
| SOSS o1 | -4.6 at 0.194 | 598 / -0.2 | 20.85 / -0.35 | 22.63 / -0.52 | 1.1e-3 / 1.6e-2 / 0.16 | 0.11 |
| SOSS o2 | -5.4 at 0.181 | 184 / +0.5 | 18.00 / -1.20 | 18.23 / -1.88 | 2.8e-4 / 8e-4 / 2.6e-3 | 0.56 |
| WFC3 | -3.2 at 0.194 | 574 / -0.2 | 21.20 / -0.25 | 22.99 / -0.38 | 1.6e-3 / 2.5e-2 / 0.16 | 0.08 |
| IRAC1 | -21.0 at 0.194 | 2415 / -0.0 | 18.80 / -0.45 | 21.45 / -0.58 | 1.1e-3 / 1.1e-2 / 7.0e-2 | 0.13 |
| IRAC2 | -35.6 at 0.194 | 2965 / +1.4 | 18.00 / -0.55 | 20.11 / -0.76 | 8e-4 / 6.6e-3 / 3.6e-2 | 0.20 |
| CO 4.5-4.8 um (top-hat, photon) | -47.0 at 0.194 | 3275 / +2.7 | 17.35 / -0.65 | 19.30 / -0.89 | 4.6e-4 / 4.1e-3 / 2.5e-2 | 0.25 |

- phi = 0.194 puts the sub-observer point at lon 110 deg, where the knot faces the observer.
- The knot, at 1e-5 to 1e-4 bar, lies above the 16th percentile of every band's contribution function except SOSS o2. It reaches the IR bands only through their low-p wing.
- Effect: amplitude changes by at most 3 ppm, and offsets move by 0.3-0.9 deg (up to 1.9 deg on SOSS o2) toward the west. The knot explains less than 1 deg of the 15-deg offset excess over ME23.
