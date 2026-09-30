# WASP-121b: published observables for comparison with the 3-D RHD runs

Web survey done 2026-09-29 (no code, no jobs). The source is given for every number. Where a
number was read from the full text (PDF) rather than the abstract, that is stated. Items marked
**[unverified]** came from a search snippet or secondary quotation and were not checked in the
primary paper.

Sign conventions: a phase-curve peak **before** mid-eclipse means the hot spot is **east** of the
substellar point (the direction of rotation, as for a super-rotating jet). A peak **after**
mid-eclipse means the hot spot is west. In high-resolution transit work, "morning" is the leading
limb (seen mostly early in the transit) and "evening" is the trailing limb (seen mostly late in
the transit). A negative velocity is a blueshift, i.e. motion towards the observer.

---------------------------------------------------------------------------------------------------

## 5. System parameters (listed first because the other sections use them)

| Quantity | Value | Source |
|---|---|---|
| P | 1.2749255 (+0.0000020/-0.0000025) d | Delrez+2016, MNRAS 458, 4025, arXiv:1506.02471 |
| P (update) | 1.27492504 +- 1.5e-7 d | Seidel+2025 Ext. Data Table 2, via Bourrier+2020 |
| M_p | 1.183 (+0.064/-0.062) M_J (Delrez+2016); 1.157 +- 0.070 M_J (Bourrier+2020 HEARTS III, A&A 635, A205, arXiv:2001.06836); **1.170 +- 0.043 M_J** (Sing+2024, AJ 168, 231, arXiv:2501.03844, from the JWST planetary RV) | |
| R_p | 1.865 +- 0.044 R_J (Delrez+2016); 1.7402 +- 0.0393 R_J (Seidel+2025); 1.7420 +- 0.0060 R_J (Sing+2024 Table 2) | |
| M_* | 1.353 (+0.080/-0.079) Msun (Delrez); 1.38 +- 0.02 (Bourrier+2020); **1.330 +- 0.019** (Sing+2024, independent of stellar models) | |
| R_* , T_eff | 1.458 +- 0.030 Rsun, 6460 +- 140 K (Delrez); 1.44 +- 0.03 Rsun (Bourrier+2020); T_eff 6628 +- 66 K (Sing+2024 Table 2) | |
| a/R_*, a | 3.7844 +- 0.0069 (Sing+2024); 3.81 +- 0.02, a = 0.02582 +- 0.00055 au (Seidel+2025); 3.8002 +- 0.005 (Splinter+2025 fit) | |
| T_eq | 2358 +- 52 K (Delrez+2016 as tabulated by Seidel+2025); 2409 +- 24 K (Sing+2024) | |
| T_irr = T_* sqrt(R_*/a) | 3320 +- 72 K | Mikal-Evans+2022 (Methods) |
| Irradiation | ~7.1e9 erg s^-1 cm^-2 | Delrez+2016 |
| K_p (planet RV semi-amplitude) | **215.7 +- 1.1 km/s** (NRS1 212.7 +- 1.8, NRS2 218.2 +- 1.3) | Sing+2024 |
| v_sys | 38.198 +- 0.002 km/s (Bourrier+2020); 38.64 +- 0.06 km/s (Seidel/Prinoth+2025 fits) | |
| Obliquity | nearly polar: 88.1 +- 0.25 deg or 91.11 +- 0.20 deg; projected lambda = -87.08 +- 0.28 deg | Bourrier+2020 |
| Roche | a is ~1.15 x the Roche limit | Delrez+2016 |
| Age | 1.11 +- 0.14 Gyr | Sing+2024 |

What to use: Sing+2024 for M_p, M_*, a/R_* and K_p (these values are model independent). Note that
Vaulato+2025 used M_* = 1.38 and found Delta K_p = -15 +- 3 km/s, and they state that this result
depends on M_* (see section 3).

---------------------------------------------------------------------------------------------------

## 1. Phase curves

### 1a. Broadband and white-light phase-curve parameters

| Band / instrument (epoch) | Eclipse (day) depth Fp/F* | Nightside Fp/F* | Amplitude | Peak offset | T_b,day | T_b,night | Source |
|---|---|---|---|---|---|---|---|
| TESS 0.6-1.0 um (2018-19, Sectors 7/8) | 482 (+41/-39) ppm **[unverified]** | 98 (+22/-18) ppm **[unverified]** | ~400-500 ppm | no significant offset | 3012 (+40/-42) K | 2022 (+254/-602) K | Daylan+2021, AJ 161, 131, arXiv:1909.03000 |
| TESS (same data, independent analysis) | -- | -- | -- | hot spot at the substellar point | 2870 +- 50 K (pure thermal) | < 2200 K (3 sigma) | Bourrier+2020b, A&A 637, A36, arXiv:1909.03010; the alternative pure-reflection fit gives A_g ~ 0.37 |
| HST/WFC3 G141 1.12-1.64 um, 2 epochs (Mar 2018, Feb 2019) | semi-amplitude c1 = 1195 ppm (2018), 1158 ppm (2019) (read from Ext. Data Table 1; errors garbled in the PDF extraction) | -- | -- | broadband peak ~6 deg **before** eclipse = east (Morello+2023 quoting ME22); the l=2 map peaks ~9 deg east at ~3200 K, coldest nightside ~1200 K | 2760 +- 100 K (blackbody, computed by Morello+2023, not by ME22) | 1665 +- 65 K (same) | Mikal-Evans+2022, Nat. Astron. 6, 471, doi:10.1038/s41550-021-01592-w, arXiv:2202.09884 |
| Spitzer 3.6 um (Jan 2018, Program 13242) | 4230 +- 80 ppm (F_day max) | 50 +- 240 ppm | -- | **5.9 +- 1.6 deg after eclipse = WEST** | 2670 (+55/-40) K | 710 (+270/-710) K | Morello+2023, A&A 676, A54, arXiv:2307.00669 (Table 5) |
| Spitzer 4.5 um (same program) | 5090 +- 90 ppm | 710 +- 270 ppm | -- | **5.0 (+3.4/-3.1) deg WEST** | 2700 (+70/-50) K | 1130 (+130/-160) K | Morello+2023 |
| Spitzer 3.6 um, reanalysis (apparently the same Program 13242 data) | 4077 +- 59 ppm | min 553 +- 96 ppm | 1771 +- 95 ppm | -0.78 +- 1.87 deg (consistent with 0) | 2779 +- 40 K (abstract) / 2543 +- 42 K (Sec. 3 text) | 1259 +- 67 K (abstract) / 1082 +- 61 K (text) | Davenport+2025, AJ, doi:10.3847/1538-3881/adc0a4, arXiv:2503.12521 |
| Spitzer 4.5 um, reanalysis | 5121 +- 76 ppm | min 1045 +- 107 ppm | 2048 +- 109 ppm | 0.42 +- 1.74 deg | 2905 +- 51 K | 1349 +- 54 K | Davenport+2025 |
| Spitzer eclipses only (2017, Program 13038) | 3685 +- 114 ppm (3.6), 4684 +- 121 ppm (4.5) | -- | -- | -- | 2490 +- 77 K, 2562 +- 66 K | -- | Garhart+2020 as quoted by Morello+2023 and ME23. Morello attributes the ~400-550 ppm deficit to the neglected phase blend |
| **JWST NIRSpec G395H NRS1 2.70-3.72 um** (2022 Oct 14-15, 37.8 h, GO-1729) | **3924 +- 7 ppm** | **136 +- 8 ppm** | -- | **3.36 +- 0.11 deg before eclipse = EAST** | 2762 (+30/-32) K | **926 +- 12 K** | Mikal-Evans+2023, ApJL 943, L17, doi:10.3847/2041-8213/acb049, arXiv:2301.03209 |
| **JWST NIRSpec G395H NRS2 3.82-5.15 um** | **4924 +- 9 ppm** | **630 +- 10 ppm** | -- | **2.66 +- 0.12 deg EAST** | 2768 +- 39 K | **1122 +- 10 K** | Mikal-Evans+2023 |
| **JWST NIRISS/SOSS Order 1 (~0.85-2.85 um)** | 1156 +- 14 ppm | -- | -- | **5.1 +- 1.4 deg EAST** | -- | -- | Splinter+2025, AJ 170, 323, arXiv:2509.09760 (Table 1, Sec. 3.2) |
| NIRISS/SOSS Order 2 (~0.6-0.85 um) | 363 (+21/-22) ppm | -- | -- | 10.5 +- 9.9 deg (not significant) | -- | -- | Splinter+2025 |

Offsets per wavelength: Splinter+2025 find near-zero offsets above 1.5 um and eastward offsets that
grow towards shorter wavelengths. Frazier+2026 (arXiv:2605.01589) confirm the trend of increasing
eastward offset with decreasing wavelength below ~1.4 um and say that none of their GCMs explains it.
Yang+2026 (arXiv:2607.29057) invert the NIRSpec phase curve and find a hot-spot offset that grows
with altitude, from ~4 deg to ~9 deg.

Splinter+2025 also fit a "substellar radius" R_sub/R* = 0.1414 +- 0.0028 against Rp/R* = 0.1222
(Order 1): the dayside appears ~16 % larger than the transit radius (1.157 +- 0.023 R_p). This
tests the dayside scale height and inflation directly.

### 1b. Energy budget (bolometric)

| Quantity | Value | Source |
|---|---|---|
| **T_day (effective)** | **2717 +- 17 K** | Splinter+2025 (NIRISS/SOSS; the band holds 50-83 % of the bolometric flux) |
| **T_night (effective)** | **1562 +- 19 K** | Splinter+2025 |
| **Bond albedo A_B** | **0.277 +- 0.016**; ~0.3 from an energy-balance fit | Splinter+2025 |
| **Heat redistribution epsilon** | **0.246 +- 0.014** | Splinter+2025 |
| Energy-balance fit | mixed-layer pressure ~1 bar; "unexpectedly slow" effective winds of 0.2 km/s | Splinter+2025 |
| Geometric albedo (0.6-0.85 um) | 0.093 +- 0.029 (3 sigma < 0.175) | Splinter+2025 |
| A_B, epsilon (HST) | 0.14 +- 0.08, 0.29 +- 0.02 | Mikal-Evans+2022 |
| A_B, epsilon (Spitzer) | 0.37 (+0.07/-0.09), 0.013 (+0.034/-0.013) at 3.6 um; 0.32 (+0.08/..), 0.077 (+0.040/..) at 4.5 um; combined blackbody epsilon = 0.07 (+0.05/..) | Morello+2023 (lower error bars cut off in the extraction) |
| 3-D GCMs vs NIRISS | GCM emission is ~12 % above the observed emission; strong drag preferred; nightside clouds suggested | Frazier+2026, arXiv:2605.01589 |

### 1c. Variability between epochs

- HST WFC3, 2018 vs 2019 phase curves: the semi-amplitudes 1195 and 1158 ppm are similar
  (Mikal-Evans+2022). Changeat+2024 (ApJS 270, 34, arXiv:2401.01465) reanalysed all the HST data
  and claim (i) a shift of the "putative" hot-spot offset between the two phase curves and
  (ii) varying transit and eclipse spectra. The five eclipses span an average temperature range of
  311 K (1e5-1e3 Pa). Their GCM predicts ~5-10 % variability in disk-averaged flux with a period of
  ~5 planet days. **No per-epoch offset values are given.** Their retrieval uses a hot spot of
  size 50 deg with offset 30 deg, so the "~50 deg offset" some papers quote is a hot-spot size.
- Spitzer vs JWST: Morello's westward Spitzer offsets disagreed with the eastward JWST offsets at
  5.7 sigma (3.6 um vs NRS1) and 2.5 sigma (4.5 um vs NRS2). Morello+2023 offered weather as one
  possible explanation. The Davenport+2025 reanalysis gives offsets consistent with zero and
  reduces the tension to 2.2 sigma and 1.3 sigma. **The "westward Spitzer hot spot" is therefore
  probably a reduction artefact, not weather.**
- TESS dayside T_b 3012 K (Daylan) vs 2870 K (Bourrier) comes from the same data analysed two ways.
- Ground-based ESPRESSO abundances are consistent across epochs (Maguire+2023, below).
  Mikal-Evans+2020 (MNRAS 496, 1638, arXiv:2005.09631) did not reproduce the 1.25 um "VO" feature
  of the single 2016 eclipse.

### How to compute the same from a simulation
Phase curve: for each orbital phase phi, compute
F_p(phi) = sum over the visible hemisphere of I_nu(mu, cell) * mu * dA / d^2.
Use the outgoing top-of-atmosphere specific intensity in each correlated-k band, or at least the
flux under an Eddington/cos-mu weighting, with the sub-observer longitude = 180 deg - 360 phi.
Then:
- band-integrate with the instrument response (TESS, IRAC 1/2, NRS1/NRS2, SOSS);
- divide by the stellar flux in the same band (PHOENIX, T_eff 6460-6630 K, R_* 1.44-1.46 Rsun);
  this gives Fp/F* in ppm;
- the day and night values are phi = 0.5 and phi = 0;
- the offset is the phase of maximum. Before eclipse (phi < 0.5) means east.
Brightness temperature: invert Planck on the band-averaged disk flux, as the papers do. T_eff
day/night: the bolometric OLR averaged over each hemisphere with mu weighting.
epsilon, A_B: Cowan & Agol (2011) formulas applied to the two T_eff.
Our run lacks scattering (A_B = 0 by construction), so expect T_day to be too hot compared with
A_B ~ 0.28 (roughly (1-0.28)^0.25 = 0.92, i.e. ~8 % in T).

---------------------------------------------------------------------------------------------------

## 2. Day-night contrast (summary)

| Probe | T_day | T_night | Delta T | Source |
|---|---|---|---|---|
| Bolometric (NIRISS) | 2717 +- 17 | 1562 +- 19 | ~1155 K | Splinter+2025 |
| NRS1 / NRS2 | 2762 / 2768 | 926 / 1122 | ~1840 / ~1650 K | ME23 |
| WFC3 G141 | 2760 +- 100 | 1665 +- 65 | ~1100 K | ME22 via Morello+2023 |
| IRAC 3.6 / 4.5 (Davenport) | 2779 / 2905 | 1259 / 1349 | ~1520 / ~1560 K | Davenport+2025 |
| TESS | 3012 | 2022 (+254/-602) | poorly constrained | Daylan+2021 |

ME23 note that the nightside temperatures lie below the silicate condensation temperatures, and that
the 12-sigma difference between NRS1 and NRS2 T_night points to a nightside cloud deck and/or a
cooling-with-height nightside profile. The Parmentier+2018 cloud-free GCMs overpredict T_night.
A cloud-free, drag-free GCM should therefore be expected to be **too hot on the nightside**.

---------------------------------------------------------------------------------------------------

## 3. Winds from high-resolution spectroscopy

### 3a. Transmission (terminator) line-of-sight velocities

| Tracer (probed layer) | Morning / leading limb | Evening / trailing limb | Full-transit / other | Source |
|---|---|---|---|---|
| **Fe I CCF, ESPRESSO 4-UT** (deepest optical tracer, around 1e-3..1e-5 bar by their Fig. 1 contribution function; they state that only relative pressures are secure) | **-4.12 +- 0.15 km/s** | **-6.90 +- 0.11 km/s** | de-projected day-to-night flow **~6-10 km/s** | Seidel+2025, Nature 639, 902, arXiv:2502.12261 (Ext. Data Table 3) |
| All-species CCF, 4-UT (planet mask) | -2.68 +- 0.14 | -7.51 +- 0.16 | evening LOS -> actual wind 9.8 +- 1.6 km/s | Seidel+2025 |
| All-species CCF, 1-UT (F6 mask) | -2.80 +- 0.28 | -7.66 +- 0.16 | -- | Borsa+2021, A&A 645, A24, arXiv:2011.01245 (numbers as tabulated by Seidel+2025) |
| **Na I doublet**, MERC retrieval (mid layers, above Fe) | jet **13.7 +- 6.1 km/s** (net REDshift, i.e. away from the observer: super-rotation); T_iso 2404 +- 662 K; v_vertical 16.6 +- 10.1 | jet **26.8 +- 7.3 km/s** (net blueshift); T_iso 3350 +- 470 K; v_vertical 26.7 +- 9.3 | morning-to-evening heating **950 +- 560 K**; jet confined to less than half the latitude range (~+-30 deg); a secondary Na feature sits at ~20 km/s | Seidel+2025 (Ext. Data Table 1); evening segment from Seidel+2023, A&A 673, A125 |
| **H-alpha** (upper atmosphere, < ~1 microbar, Parker-type outflow) | +5.2 +- 1.4 km/s (red: zonal super-rotation) | -19.2 +- 1.4 km/s | extends to 1.44 R_p (Borsa) / 1.54 +- 0.04 R_p (Maguire) | Seidel+2025; Borsa+2021; Maguire+2023 |
| Fe (HARPS "reloaded RM") | -- | -- | blueshift **-5.2 +- 0.5 km/s**, width 14.3 +- 1.2 km/s | Bourrier+2020 HEARTS III |
| Metals, HARPS CCF | -- | -- | mean blueshift ~ **-5 km/s** in the low layers (day-to-night wind) | Hoeijmakers+2020 HEARTS IV, A&A 641, A123, arXiv:2006.11308 |
| Exospheric H-alpha, Fe II, Ca II | -- | -- | Delta v_sys = -1.8, -7.7, -16.2 km/s; heights 1.54, 1.17, 2.52 R_p | Maguire+2023, MNRAS 519, 1030, arXiv:2211.09621 **[velocities unverified: from a search snippet]** |
| Fe, three ESPRESSO transits | -- | -- | ~ **-6 km/s**, consistent across epochs | Maguire+2023, as quoted by Gandhi+2023 (AJ 165, 242, arXiv:2305.17228). Gandhi find similar first- and last-half winds for WASP-121b |
| Ti I, V I (4-UT) | Ti I at rest in the planet frame in the first half | Ti I increasingly blueshifted late in the transit; V I shows a smaller shift | Ti must sit at the equator of the evening terminator | Prinoth+2025, A&A 694, A284, arXiv:2502.12262 (numbers only in the figures) |
| **CO (IGRINS K band)** | Q1: +1.8 +- 1.0 | Q2 -3.2 +- 0.6; Q3 -1.6 +- 0.6; Q4 -3.1 +- 1.5 | blueshift grows through the transit | Wardenier+2024, PASP, doi:10.1088/1538-3873/ad5c9f, arXiv:2406.09641 |
| **H2O (IGRINS)** | Q1 not detected (evening clouds?); Q2 -3.3 +- 1.0 | Q3 +1.9 +- 1.2; Q4 +2.3 +- 1.2 (redshift) | H2O fits the **strong-drag** GCM (tau_drag = 1e4 s); Fe fits the weak or no-drag GCMs | Wardenier+2024 |
| Fe + CO + V (HARPS + NIRPS + CRIRES+) | -- | -- | Delta K_p = K_p,retr - K_p,orb = **-15 +- 3 km/s** (with M_* = 1.38); K_p,retr = 202.99 (+2.84/-2.92); at 1e-4..1e-3 bar | Vaulato+2025, A&A (Nov 2025), arXiv:2509.00151 |

### 3b. Emission (dayside) velocities

| Measurement | Value | Source |
|---|---|---|
| K_p from the JWST NIRSpec phase curve | NRS1 212.7 +- 1.8, NRS2 218.2 +- 1.3 km/s; the difference between detectors (pressure levels) is 5.5 +- 2.2, of which rotation explains ~1 km/s, leaving **4.5 +- 2.2 km/s** of wind difference; de-projected **~5.2 km/s** equatorial wind | Sing+2024 |
| IGRINS dayside (pre + post eclipse) | K_p = 215.28 (+0.35/-0.34), dv_sys = 1.20 (+0.13/-0.11) km/s; H2O effective K_p ~4 km/s larger than CO/OH; OH net blueshift ~1.5 km/s | Smith+2024, AJ 168, 293, arXiv:2410.19017 |
| CRIRES+ + ESPRESSO dayside | H2O blueshifted by 2-3 km/s relative to CO, Fe, Ni (asymmetric H2O from dissociation) | Pelletier+2025, AJ, arXiv:2410.18183 **[number from the search summary]** |
| Na in emission | tentative strong velocity offsets for Na I, absent for species probing deeper layers | Hoeijmakers+2024 (Mantis IV, A&A 685, A?, arXiv:2210.12847), as cited by Seidel+2025 |

### How to compute the same from a simulation
Transmission: ideally, ray-trace through the 3-D T, rho, v field (Wardenier-style) with a
Fe/Na/CO/H2O opacity, and cross-correlate time-resolved spectra.
Cheap proxy:
- at each limb (morning = west terminator, longitude -90 deg; evening = east, +90 deg) and each
  pressure level, compute the line-of-sight velocity v_LOS = -(u_phi * sin(angle) + ...);
- at the terminator of a tidally locked planet, v_LOS = -u_r(day-to-night component along the
  star-planet line) + rotation;
- concretely: v_LOS = -v_{x} (the component along the star-planet axis pointing to the observer)
  plus the solid-body rotation term +-Omega R cos(lat) ~ +-7 km/s at the equator
  (2 pi R_p / P ~ 7.1 km/s);
- average over the terminator annulus with a tau~1 slant-path weighting at ~1e-3..1e-5 bar (Fe),
  ~1e-4..1e-6 bar (Na; our top is 1e-8 bar), and 1e-3 bar (CO/H2O at IGRINS).
Include the rotation angle: the planet turns by ~2 pi * T14/P = 2 pi * 2.9 h / 30.6 h ~ 34 deg
during transit, so early exposures see the west limb plus the nightside-facing part, and late
exposures see the east limb.
Compare with -4 / -7 km/s (Fe, morning/evening) and a Na jet of +14 / -27 km/s.
Emission: disk-integrated, intensity-weighted v_LOS over the visible dayside at the NRS1/NRS2
photosphere, compared as a K_p shift.

---------------------------------------------------------------------------------------------------

## 4. Temperature-pressure structure and limb asymmetry

| Observable | Value | Pressure | Source |
|---|---|---|---|
| Dayside inversion (H2O in emission at 1.3-1.6 um, 5 sigma) | T rises from ~2500 K to ~2800 K | ~30 -> 5 mbar | Evans+2017, Nature 548, 58, arXiv:1708.01076 (T-P numbers **[unverified: search snippet]**) |
| Confirmation (5 WFC3 eclipses) | blackbody ruled out at > 6 sigma; inversion confirmed; 1.25 um VO feature not reproduced | -- | Mikal-Evans+2020 |
| G102 0.8-1.1 um | H- emission shortward of 1.1 um; [M/H] = 1.09, C/O = 0.49 | -- | Mikal-Evans+2019, MNRAS 488, 2222, arXiv:1906.06326 |
| Dayside vs nightside profile | H2O feature in emission on the dayside (warming with height), in absorption on the nightside (cooling with height); H2O dissociated at low p on the dayside; nightside cold enough for CaTiO3 (and Mg, Fe, V condensates) | ~0.01-0.1 bar (their GCM comparison at 10 mbar) | Mikal-Evans+2022 |
| Inversion base (all HST) | strong inversion starting at **~0.1 bar**; hottest region ~300 K hotter than the rest of the dayside (1e5-1e3 Pa); -0.77 < log Z < 0.05; 0.59 < C/O < 0.87 | 1e5-1e3 Pa | Changeat+2024 |
| TiO / Ti | TiO not detected; Ti cold-trapped (Hoeijmakers+2024 Mantis IV); Ti I later detected at ~5 sigma per spectrum, ~19 sigma stacked, only on the evening equatorial limb (Prinoth+2025); NIRISS dayside Ti deficient (arXiv:2508.18341) | -- | as listed |
| Nightside chemistry | CH4 on the nightside (3.1-5.1 sigma) -> strong vertical mixing / disequilibrium; dayside H2O, CO, SiO | NRS phase curve | Evans-Soma+2025, Nat. Astron., arXiv:2506.01771 |
| Terminator heterogeneity (transmission) | H2O and H2 dissociated on the dayside; more H2O on the nightside than the dayside; SiO at 5.2 sigma | G395H | Gapp+2025, AJ 169, 341, arXiv:2506.02199 |
| **Rotational transit (limb asymmetry vs time)** | CO absorption grows by ~200 ppm as the planet rotates; H2O roughly constant. Radius-rate coefficient R1/R* = 409 +- 53 ppm/h (NIRISS white light); CO band 650 +- 188 (NIRISS), 580 +- 90 ppm/h (NIRSpec). **Stronger longitudinal T gradient across the evening than the morning terminator**, i.e. the east half of the dayside is hotter. SPARC/MITgcm underestimates the rate (day-night contrast too small; morning clouds missing) | -- | Gapp+2026, Nat. Astron., arXiv:2606.19487 |
| Static morning/evening limb split | no evidence for asymmetric limbs in the white or spectroscopic transit light curves (catwoman fit) | -- | Gapp+2025 (as summarised by search) **[unverified]** |
| Na-derived limb temperatures | morning 2404 +- 662 K, evening 3350 +- 470 K | Na layer | Seidel+2025 |
| 2.5-D inversion of the NIRSpec phase curve | dayside inversion extends to both limbs; limb temperatures differ; offset grows with altitude from 4 to 9 deg; more CH4 on the nightside; clouds on the nightside and morning side | -- | Yang+2026, arXiv:2607.29057 |
| Dayside composition (JWST 0.6-12 um incl. MIRI LRS) | Si/O = 3.54 (+0.86/-0.69) x stellar, Si/C = 3.05 (+1.12/-0.80) x stellar | -- | Kahle+2026, arXiv:2606.21855 (submitted to A&A) |
| Dayside C/O (IGRINS) | 0.70 (+0.07/-0.10); refractory/volatile 3.83 (+3.62/-1.67) x stellar | -- | Smith+2024 |
| Volatile/refractory (CRIRES+ + ESPRESSO) | 1.75 (+0.57/-0.41) x stellar | -- | Pelletier+2025 |

Disagreement: Smith+2024 (refractory/volatile super-stellar) and Pelletier+2025 (volatile-enriched
relative to refractory) are in tension, with large error bars. The metallicity estimates also
differ: ME22 retrieve [M/H] ~ 0.7 (5x solar) for the dayside, Changeat+2024 find subsolar to solar,
and Evans-Soma+2025 find super-stellar C/H and O/H. **Our 1x and 10x solar runs bracket the range.**

How to extract: take the dayside-mean and nightside-mean T(p), weighted by mu over each visible
hemisphere, and compare with an inversion base near 0.1 bar and ~2500 -> 2800 K between 30 and
5 mbar. For the limbs, average T and the H2O/CO abundance (if the chemistry is tracked) over
+-30 deg of longitude around each terminator at 1-10 mbar, and check that the evening limb is
hotter and has a steeper east-west gradient than the morning limb.

---------------------------------------------------------------------------------------------------

## Time variability ("weather") claims, collected

1. Changeat+2024: hot-spot shift between the 2018 and 2019 HST phase curves; spread of five
   eclipses 311 K; GCM variability 5-10 % on ~5 planet days. **Disputed** by the consistency of the
   ME22 amplitudes; no per-epoch offsets are published.
2. Spitzer-westward vs JWST-eastward offsets: largely removed by the Davenport+2025 reanalysis
   (see 1c).
3. Seidel+2025: the hot-spot location "is likely variable" (citing Changeat); not measured.
4. Against weather: consistent ESPRESSO abundances over months to years (Maguire+2023), and
   Gandhi+2023 winds agree with Maguire's three epochs (~-6 km/s).
5. NIRSpec pre-eclipse spectra have slightly higher flux than post-eclipse spectra (~1 sigma and
   0.3 sigma for the two eclipses) (Yang+2026 summary) **[unverified detail]**.

Simulation test: the standard deviation over time of the disk-integrated band flux at phase 0.5
and of the phase of peak, with a window of ~5 orbits. Compare with a ~5-10 % flux spread and a few
degrees of offset.

---------------------------------------------------------------------------------------------------

## References (arXiv id / DOI)

- Bourrier V. et al. 2020a, A&A 635, A205, HEARTS III, arXiv:2001.06836
- Bourrier V. et al. 2020b, A&A 637, A36, TESS phase curve, arXiv:1909.03010
- Borsa F. et al. 2021, A&A 645, A24, arXiv:2011.01245
- Changeat Q. et al. 2024, ApJS 270, 34, arXiv:2401.01465
- Davenport B. et al. 2025, AJ, doi:10.3847/1538-3881/adc0a4, arXiv:2503.12521
- Daylan T. et al. 2021, AJ 161, 131, arXiv:1909.03000
- Delrez L. et al. 2016, MNRAS 458, 4025, arXiv:1506.02471
- Evans T.M. et al. 2017, Nature 548, 58, doi:10.1038/nature23266, arXiv:1708.01076
- Evans-Soma T.M. et al. 2025, Nat. Astron., arXiv:2506.01771
- Frazier R.C. et al. 2026, arXiv:2605.01589
- Gandhi S. et al. 2023, AJ 165, 242, arXiv:2305.17228
- Gapp C. et al. 2025, AJ 169, 341, arXiv:2506.02199
- Gapp C. et al. 2026, Nat. Astron., doi:10.1038/s41550-026-02887-6, arXiv:2606.19487
- Hoeijmakers H.J. et al. 2020, A&A 641, A123, HEARTS IV, arXiv:2006.11308
- Hoeijmakers H.J. et al. 2024, A&A, Mantis IV (Ti cold trap), arXiv:2210.12847
- Kahle K.A. et al. 2026, arXiv:2606.21855
- Maguire C. et al. 2023, MNRAS 519, 1030, arXiv:2211.09621
- Mikal-Evans T. et al. 2019, MNRAS 488, 2222, arXiv:1906.06326
- Mikal-Evans T. et al. 2020, MNRAS 496, 1638, doi:10.1093/mnras/staa1628, arXiv:2005.09631
- Mikal-Evans T. et al. 2022, Nat. Astron. 6, 471, doi:10.1038/s41550-021-01592-w, arXiv:2202.09884
- Mikal-Evans T. et al. 2023, ApJL 943, L17, doi:10.3847/2041-8213/acb049, arXiv:2301.03209
- Morello G. et al. 2023, A&A 676, A54, arXiv:2307.00669
- Pelletier S. et al. 2025, AJ, doi:10.3847/1538-3881/ad8b28, arXiv:2410.18183
- Prinoth B. et al. 2025, A&A 694, A284, arXiv:2502.12262
- Seidel J.V. et al. 2023, A&A 673, A125 (high-velocity Na feature)
- Seidel J.V. et al. 2025, Nature 639, 902, doi:10.1038/s41586-025-08664-1, arXiv:2502.12261
- Sing D.K. et al. 2024, AJ 168, 231, arXiv:2501.03844
- Smith P.C.B. et al. 2024, AJ 168, 293, doi:10.3847/1538-3881/ad8574, arXiv:2410.19017
- Splinter J. et al. 2025, AJ 170, 323, arXiv:2509.09760
- Vaulato M. et al. 2025, A&A, arXiv:2509.00151
- Wardenier J.P. et al. 2024, PASP, doi:10.1088/1538-3873/ad5c9f, arXiv:2406.09641
- Yang Y. et al. 2026, arXiv:2607.29057

Not found or not covered: a JWST MIRI phase curve (only a MIRI/LRS dayside spectrum, in Kahle+2026);
a second JWST NIRSpec phase-curve epoch; Sicilia+2025; per-epoch HST offsets. The He I 10833
outflow (CRIRES+ 2024; arXiv:2510.09809, helium over more than half the orbit) lies above our domain.
