# Seidel et al. 2026: ultra-hot Jupiter wind survey, a note for picking dhj targets

Written 2026-11-01. Source PDF (ESO science-release copy of the accepted paper):
`/viper/ptmp2/jinma/dhj_planets_1101/src/eso2606a.pdf` (text: `src/eso2606a.txt`).
Downloaded from https://www.eso.org/public/archives/releases/sciencepapers/eso2606/eso2606a.pdf

## 1. Bibliographic data and summary

- **Title:** "Magnetic field strengths of hot giant exoplanets consistent with Solar System values"
- **Authors:** J. V. Seidel, V. Parmentier, B. Prinoth, T. Hood, N. Mehta, V. De Lia, K. Batygin (Caltech GPS),
  T. Guillot, R. Van den Broeck, H. Beltz, B. Thorsbro, F. Debras, D. D. B. Koll, T. D. Komacek, E. Rauscher,
  L. Pino, M. Brogi, J. P. Wardenier, J. L. Bean, B. Benneke, J.-M. L. B. Desert, P. Drake, S. Gandhi,
  M. Hammond, D. Kasper, M. R. Line, E. K. H. Lee, S. Pelletier, A. Seifahrt, A. Simonnin, P. C. B. Smith,
  K. B. Stevenson (32 authors).
- **Journal:** Nature Astronomy, published online 2 June 2026, **DOI 10.1038/s41550-026-02870-1**
  (ESO press release eso2606).
- **arXiv:** 2606.22455, submitted 21 June 2026 (verified 11-01 on arxiv.org/abs/2606.22455; the first search missed it).
- **Other Seidel-led work** that matters here, but is not the 2026 sample paper: Seidel et al. 2025,
  "Vertical structure of an exoplanet's atmospheric jet stream" (WASP-121b, ESPRESSO 4-UT, arXiv:2505.08926),
  and Seidel et al. 2021, "Into the storm" (WASP-76b, A&A 653, A73). Seidel is also a co-author on a 2026
  WASP-31b ESPRESSO paper (arXiv:2602.15501), which is not a wind sample. The sample paper is the
  Nature Astronomy paper above.

**Summary.** The paper uses archival and PI transit time series from ESPRESSO (VLT, 1-UT HR mode and 4-UT MR
for WASP-121b) and MAROON-X (Gemini-North), and reduces all of them the same way. The pipeline is molecfit
telluric correction, a StarRotator Doppler-shadow (RM) fit, systemic velocity re-measured per night (a
rotationally broadened Voigt fit to the out-of-transit CCF), and cross-correlation with a single 3000 K
isothermal **Fe I** template (tayph). Fe I is detected in 7 ultra-hot Jupiters. For each planet the paper
reports **one number**: the net line-of-sight (blueshift) "day-to-night wind" at the terminator, read off a
stacked Kp-vsys map marginalised over Kp. The wind falls from about 7 km/s (WASP-76b, Teq 2229 K) to about
1 km/s (WASP-189b, 2638 K). The authors attribute this to Ohmic drag and use a Koll & Komacek (2018)
heat-engine scaling, calibrated with k0 = 0.25 against the SPARC/MITgcm grid, to infer B_max of 2 G or less for
the five planets between 2200 and 2500 K, and up to about 6 G for the two hottest. A cross-correlation
contribution function puts the Fe signal at **1-100 mbar, peaking at 30 mbar**. The authors say this is an
upper bound and the true level is probably higher. The drag times needed run from about 7 h at 2200 K to under 15 min for
WASP-189b (ED Fig. 5).

## 2. Planet table: what the paper itself measures

The paper prints **no table of wind values**: Fig. 1 is the only source. I digitised Figs. 1 and 3 from the
vector PDF paths, so the conversion is exact to about 0.01 km/s. The Teq values are the paper's own
zero-albedo recomputation (ED Table 1). The digitised Teq of WASP-121b is 2339 K, but ED Table 1 lists 2358 K.
Errors are the plotted 1-sigma bars, which are symmetric. B_max is the midpoint and range of the Fig. 3
"original from data" bar (digitised). The table's species, asymmetry, dKp, width, rotation, altitude and
phase entries are blank because the paper does not report them: it measures one Fe I offset per planet.

| Planet | Teq [K] (paper) | Instr., nights | Species used/detected | Net LOS blueshift ("day-night wind") [km/s] | B_max [G] (Fig. 3) | Limb asym. / dKp / width / rotation / altitude / phase-resolved |
|---|---|---|---|---|---|---|
| WASP-76b  | 2229 +/- 36 | E x2, MX x3 | Fe I | **6.95 +/- 0.36** | 1.6-1.9 | not in this paper |
| KELT-20b (MASCARA-2b) | 2263 +/- 68 | MX x1 | Fe I | **5.52 +/- 1.30** | 1.7-2.8 | not in this paper |
| WASP-121b | 2358 (ED T1) +/- 51 | E-4UT x2 used (1-UT nights not used) | Fe I | **4.75 +/- 0.30** | 2.0-2.4 | not in this paper |
| WASP-178b (KELT-26b) | 2398 +/- 77 | E x2 | Fe I | **4.79 +/- 1.19** | 1.6-2.7 | not in this paper |
| TOI-1518b | 2491 +/- 38 | MX x3 | Fe I | **3.82 +/- 0.50** | 1.9-2.5 | not in this paper |
| HAT-P-70b | 2559 +/- 77 | MX x2 | Fe I | **2.55 +/- 0.74** | 2.3-6.8 | not in this paper |
| WASP-189b | 2638 +/- 28 | E x1, MX x2 | Fe I | **1.05 +/- 1.59** (consistent with 0) | 2.9-7.5 | not in this paper |

Rejected targets (ED Table 1, with the reasons in Sec. 5.1.1):
- pulsating hosts: KELT-4Ab (1827 K), KELT-7b (2048 K), MASCARA-4b (2405 K);
- no Fe detected: KELT-11b (1712 K), WASP-172b (1740 K), WASP-173Ab (1880 K, high mass);
- other: KELT-17b (Am host), WASP-19b (low S/N);
- outside the survey: WASP-33b (pulsator) and KELT-9b (in a class of its own).

What else the paper says:
- The systematics, all of which bias B lower, are a contribution function 1-2 dex lower in pressure, a limb
  temperature up to 200 K above Teq, and an LOS projection that underestimates the true wind by up to 30 %.
  With all of them combined, B_max falls to about 0.5-1 G (Fig. 3).
- The authors find no correlation of the wind with orbital period, mass, radius or gravity.
- Methods 5.1.6 warns that partial transits bias the mean velocity, citing the TOI-1518b disagreement with
  PEPSI, and that different species probe different layers.

**Limb-resolved and time-resolved information on these planets is in the earlier single-planet papers that the
sample paper cites**, not in the sample paper itself. These are what to compare against for the extra
observables:
- WASP-76b: Ehrenreich et al. 2020 (Nature; Fe asymmetry), Kesseli & Snellen 2021, Seidel et al. 2021,
  Pelletier et al. 2023 (MAROON-X).
- WASP-121b: Borsa et al. 2021, Maguire et al. 2023, Seidel et al. 2025 (vertical jet structure, ESPRESSO 4-UT).
- WASP-178b: Damasceno et al. 2024 (arXiv:2406.08348).
- TOI-1518b: Simonnin et al. 2025 (A&A 698, A314; six species, time-resolved, "strong drag");
  Basinger et al. 2025 (PEPSI).
- WASP-189b: Prinoth et al. 2023 (time-resolved, arXiv:2308.04523), Prinoth et al. 2024 (MAROON-X atlas).

I did not extract numbers from those papers for this note.

## 3. Suitability for our code

What our code assumes:
- synchronous rotation, circular orbit;
- equilibrium chemistry with ck tables for an H2 atmosphere;
- no clouds;
- a cubed-sphere grid down to about hundreds of bar;
- 11 bands, with UV below 0.26 um folded into band 10 (KELT-20b setup);
- existing setups for WASP-121b and KELT-20b.

The host and orbit facts below are **literature values, not from Seidel et al.** Where the archive table is
quoted in Sec. 6, it supersedes the approximate values here.

| Planet | Host (lit., approx.) | Orbit / rotation (lit.) | Gravity (lit.) | Flags for our code |
|---|---|---|---|---|
| WASP-76b | F7, Teff ~6300 K, slow rotator | P 1.81 d, circular, aligned | low: Rp ~1.85 RJ, Mp ~0.9 RJ, g ~6-7 m/s^2 | **Large H (puffy), but not extreme.** F host, so UV is within what the bands handle. Coolest member, so the nightside is the coolest in the sample: check that the ck tables reach the nightside T at 1-100 mbar (Fe/TiO-VO condensation is not modelled; no clouds). Best wind lever arm (fastest wind, smallest error). |
| KELT-20b | A2, Teff ~8700-9000 K, vsini ~115 km/s | P 3.47 d, circular, aligned | mass only an upper limit (~2 MJ assumed in our setup) | **Fast-rotating A host: strong short UV** (about 14 % below 0.26 um, folded into band 10). Gravity uncertain. Setup already exists. The single MAROON-X night gives a large error (1.3 km/s). |
| WASP-121b | F6, Teff 6628 K (archive), vsini 13.5 | P 1.27 d, circular, near-polar orbit | g ~9 m/s^2 | Baseline; the current production target. Most precise wind (0.30 km/s), plus limb- and altitude-resolved data from Seidel 2025. |
| WASP-178b | A1, Teff ~9300 K, small vsini (probably seen near pole-on) | P 3.34 d, circular | Mp ~1.7 MJ, Rp ~1.8 RJ, g ~13 m/s^2 | **Hottest host in the sample: strongest short UV** (more than KELT-20). The tables and band set need checking for UV and Mg/SiO-type absorbers. Error 1.2 km/s. |
| TOI-1518b | A, Teff ~7300 K, vsini ~85 km/s | P 1.90 d, circular, strongly misaligned | Mp ~2 MJ (MAROON-X), Rp ~1.9 RJ | **Fast-rotating A host (UV)**, though cooler than KELT-20/WASP-178. Short period, so close to WASP-121b in rotation. Time-resolved six-species data exist (Simonnin 2025): **richest comparison set after WASP-121b.** |
| HAT-P-70b | A2, Teff ~8450 K, vsini ~100 km/s | P 2.74 d, retrograde | mass only an upper limit (several MJ) | **Gravity unknown** (mass upper limit), so the setup is ill-posed. Strong UV. Wind error 0.74 km/s. |
| WASP-189b | A6, Teff ~8000 K, vsini ~93 km/s | P 2.72 d, polar orbit | Mp ~2 MJ, Rp ~1.6 RJ, g ~18 m/s^2 | **Hottest planet (Teq 2640 K, dayside > 3000 K): beyond the high-T end** of our tables; check the H- / H2-dissociation / ionisation coverage. Wind consistent with zero. MHD-dominated regime (drag time < 15 min), so hydro will fail by construction. Strong UV. |
| WASP-127b (**not in the paper**; Caltech-named) | G5, Teff ~5600 K, slow rotator | P 4.18 d, circular | **very low: Mp ~0.16-0.18 MJ, Rp ~1.3 RJ, g ~2.5 m/s^2** | **Not an ultra-hot Jupiter (Teq ~1400 K), so it is outside this paper.** Flags: **huge H** (about 4x WASP-121b), so the grid needs to be deeper in height and the dt is affected; **cool: clouds expected and the nightside is near or below the low-T end of the ck tables**; slower rotation (4.2 d). Comparison data: Nortmann et al. 2025 (CRIRES+, supersonic equatorial jet with a limb-resolved signature, cool poles). |

Eccentric orbits: none in the sample. All are close-in and probably circular; none is flagged eccentric in
the literature as far as I know.

## 4. Ranking next to WASP-121b

The criteria are, in order:
1. a precise Seidel wind;
2. extra resolved data;
3. physics within our code (UV, T range, known gravity);
4. lever arm in Teq along the paper's trend.

1. **TOI-1518b.** Precise (3.82 +/- 0.50 km/s) and three MAROON-X nights. Published time-resolved data in
   six species. Teq is 130 K hotter than WASP-121b, which tests the drag trend. Period is close to
   WASP-121b's. The A host's UV is moderate (Teff ~7300 K). **Best Caltech-named choice.**
2. **WASP-76b** (not Caltech-named). **Better than KELT-20b or WASP-127b.** It has the most precise and
   fastest wind (6.95 +/- 0.36 km/s), an F host (no UV problem), the best-studied limb asymmetry in Fe, and
   it anchors the cool end of the trend. Our hydro model should reach its speed if drag is weak there; this
   is the cleanest test of whether the WASP-121b wind is set by drag.
3. **KELT-20b.** The setup is ready and Teq is close to WASP-76b's. But the wind error is large
   (1.3 km/s, one night), the mass is only an upper limit, and the A-star UV is strong. Useful, but it
   discriminates weakly against Seidel.
4. **WASP-127b.** Not in Seidel 2026, so there is no comparison number from this paper. It is cool and very
   low gravity, and clouds matter. Its value is a different test (CRIRES+ jet and cool poles), so it is
   off-axis for this comparison.

Others in the paper:
- **WASP-178b** is a reasonable 5th choice. Its Teq lies between WASP-121b and TOI-1518b, but the wind is
  imprecise and the host UV is the hardest in the sample.
- **HAT-P-70b** is poor: the gravity is unknown.
- **WASP-189b** is poor: it lies beyond the tables' hot end, its wind is consistent with zero, and only MHD
  with strong drag is meaningful there. Choose it later, as an MHD test only.

Recommended set (superseded by the parameter filter in Sec. 6(e)): **WASP-121b, TOI-1518b, WASP-76b,
KELT-20b**. Add WASP-127b only if Caltech wants the cooler, low-g regime, accepting that there is no
cloud model. After the filter, KELT-20b drops out (no measured mass) and WASP-189b and WASP-178b come in;
see 6(e).

## 5. Observables

**Direct from model output** (velocity, T and density fields; no line RT):
- **Terminator-mean LOS velocity**, weighted over 1-100 mbar or at 30 mbar, as a single number to compare
  with Seidel's net blueshift per planet. Use an opacity- or density-weighted mean over the terminator annulus,
  projected on the line of sight and including synchronous rotation.
- The same quantity **separately for the leading and trailing limbs**, and versus pressure (1, 10, 30,
  100 mbar), to compare with the literature limb asymmetries and Seidel 2025's altitude structure. This is
  model-only; Seidel 2026 has no such number.
- **Day-night flow speed**: the substellar-to-antistellar mass-flux-weighted velocity and the equatorial jet
  speed versus pressure. Compare these with the paper's quoted factor-of-1.3 LOS-to-true conversion and with
  the drag times of ED Fig. 5. Our hydro run sets the zero-drag reference; the MHD runs test the B < 2 G claim
  directly, which is the paper's own call for "self-consistent MHD follow-up".
- **Limb temperature** at 1-100 mbar versus Teq: this checks the paper's assumption of T_limb ~ Teq, which
  it says can be off by up to 200 K for WASP-121b.
- **Electron fraction or conductivity** at 30 mbar on the limb from our EOS. Compare it with the paper's
  Saha tables (ED Fig. 3).

**Needs post-processing** (line RT, which we do not have yet):
- **Fe I CCF velocity offsets and Kp-vsys maps.** These are the paper's actual observable. They need
  high-resolution Fe I opacities, ray tracing through the 3-D terminator, the paper's 3000 K template and
  instrument broadening (2.14 km/s ESPRESSO 1-UT, 4.28 km/s 4-UT MR, 3.53 km/s MAROON-X), and the same
  marginalisation over Kp. Any **dKp**, **line widths** and **phase-resolved (ingress/egress) shift curves**
  need this too.
- **The real contribution function** of Fe lines. The paper admits that its 30 mbar figure is probably too
  deep, so a model-based contribution function would itself be a result.
- Species other than Fe (Na, H-alpha, Ca+, Mg, which probe higher layers) need line RT plus chemistry and
  ionisation beyond what the ck tables carry.

Practical shortcut until line RT exists: weight the LOS velocity on terminator rays by an approximate Fe I
line opacity (equilibrium Fe abundance times a Boltzmann factor at T) and take the absorption-weighted
centroid. This is cruder than full RT but captures the altitude weighting.

## 6. Sample along the trend (for reproducing the paper's points)

### (a) The paper's main relation

- **x = Teq** (zero albedo, recomputed by the authors). **y = net line-of-sight Fe I velocity**, called the
  "day-to-night wind speed at the terminator", in km/s.
- It is shown in **Fig. 1**, and the same points appear in Fig. 2 with model curves. Teq values are in
  ED Table 1. The paper prints **no table of velocities**; the values in Sec. 2 are digitised from the
  vector Fig. 1.
- Range covered: Teq **2229-2638 K**, v **6.95 to 1.05 km/s**, monotonically decreasing apart from
  WASP-121b and WASP-178b, which are equal within errors.
- **The paper fits no line.** It overlays Ohmic-drag heat-engine curves for B = 1, 2, 3, 10, 30 and 100 G,
  plus a shear curve (v proportional to Teq^(4/3)) and a shock curve (v proportional to sqrt(Teq)).
- The data follow about 2 G between 2200 and 2500 K and rise to 3-6 G at the hot end. The derived
  relation, B_max vs Teq, is Fig. 3; the drag time vs Teq is ED Fig. 5.
- My own weighted linear fit to the digitised points (not in the paper) is
  **v = 4.44 +/- 0.21 km/s - (1.24 +/- 0.18) km/s per 100 K x (Teq - 2400 K)**, with chi^2 = 4.9 for
  5 degrees of freedom.
- The two quantities to reproduce per planet are therefore **v_LOS(Teq)**, in hydro and MHD, and, from the
  MHD runs, whether B of about 2 G reproduces the slope.

### (b) Picks

Sorted by Teq. The four core picks are evenly spaced at about 135 K; two more are optional. **The final
selection, after the system-parameter filter, is in 6(e) and supersedes this table's roles.**

| Role | Planet | Teq [K] (paper) | v_LOS [km/s] (paper Fig. 1, digitised) | B_max [G] (Fig. 3 bar) | Host Teff [K] | P [d] | M_p [MJ] | R_p [RJ] | Why |
|---|---|---|---|---|---|---|---|---|---|
| **low end (core)** | WASP-76b | 2229 +/- 36 | 6.95 +/- 0.36 | 1.6-1.9 | 6329 | 1.810 | 0.894 +/- 0.014 | 1.854 | Smallest error in the sample, F host, circular, well-measured mass. |
| optional (low) | KELT-20b | 2263 +/- 68 | 5.52 +/- 1.30 | 1.7-2.8 | 8720 | 3.474 | <3.38 (upper limit) | 1.741 | Setup already exists. Large error and mass only an upper limit, so it duplicates WASP-76b's Teq. |
| **mid-low (core, baseline)** | WASP-121b | 2358 +/- 51 | 4.75 +/- 0.30 | 2.0-2.4 | 6628 | 1.275 | 1.17 +/- 0.04 | 1.742 | Production target; most precise wind. |
| optional (mid) | WASP-178b | 2398 +/- 77 | 4.79 +/- 1.19 | 1.6-2.7 | 9360 | 3.345 | 1.66 +/- 0.12 | 1.81 | Fills 2358-2491. Hottest host (strong UV), error 1.2 km/s. |
| **mid-high (core)** | TOI-1518b | 2491 +/- 38 | 3.82 +/- 0.50 | 1.9-2.5 | 7299 | 1.903 | 1.83 +/- 0.47 | 1.878 | Second-smallest relative error; time-resolved six-species data (Simonnin 2025). A host, vsini 85 km/s. |
| **high end (core)** | WASP-189b | 2638 +/- 28 | 1.05 +/- 1.59 | 2.9-7.5 | 8000 | 2.724 | 1.99 +0.16/-0.14 | 1.619 | The only hot-end planet with a measured mass. Polar orbit (irrelevant for synchronous rotation). Hot dayside: check the ck/EOS high-T coverage first. v is consistent with 0, so treat it as an upper limit of about 2.6 km/s. |

The host Teff, P, M_p and R_p columns come from the NASA Exoplanet Archive pscomppars table, queried on
2026-11-01 (`/viper/ptmp2/jinma/dhj_planets_1101/src/pscomp.csv`). The archive's "mass" for KELT-20b
(3.38 MJ) and HAT-P-70b (6.78 MJ) has no error bars and is an upper limit. The paper instead takes those two
masses from transit scaling (Gandhi et al. 2023), and its Supplementary Table 1 (not in the ESO PDF) holds the
values it used. Archive Teq values (for example WASP-189b at 3353 K, which includes a heat-redistribution
choice) differ from the paper's. **Use the paper's Teq for the x axis.**

Not picked: **HAT-P-70b** (2559 K, 2.55 +/- 0.74 km/s). Its mass is only an upper limit (6.78 MJ), so the
gravity, and with it the setup, is undetermined. Use it only if a mass appears.

### (c) Where the Caltech-named planets fall, and the gaps

- **WASP-121b**: in the middle of the trend (2358 K, 4.75 km/s), right on the fitted line (4.44 + 0.52 =
  4.96 predicted).
- **KELT-20b**: at the low end (2263 K, 5.52 +/- 1.30), next to WASP-76b. It adds little leverage, and its
  error bar is large.
- **TOI-1518b**: upper-middle (2491 K, 3.82 +/- 0.50), on the line (3.31 predicted, within 1 sigma).
- **WASP-127b**: **off the relation**. At Teq about 1400 K it lies below the paper's Teq > 1600 K
  Ohmic-drag window, and Fe I is not detectable there. The paper's trend cannot be compared with it.
  Its own test is the CRIRES+ jet.

So the four Caltech planets cover only **2263-2491 K**, three of the seven Teq points. They miss both ends:
**WASP-76b** is needed for the fast, low-drag end, and **WASP-189b** (or HAT-P-70b) for the slow, strong-drag
end that sets the slope.

Gaps that remain after the core four:
1. **2491-2638 K**, filled only by HAT-P-70b, whose mass is unknown;
2. nothing below 2229 K. Fe I is not detected in the paper's cooler targets, so there is no data point to
   reproduce there;
3. the hot end (WASP-189b) has the largest error and is consistent with 0 km/s. That end of the slope is
   pinned more by HAT-P-70b and TOI-1518b than by WASP-189b.

Hydro runs alone should give v rising with Teq (shear/shock scaling). Reproducing the falling trend needs the
MHD runs (or a drag parametrisation), which is the paper's own conclusion.

### (e) System-parameter filter (hard filter / strong ranking factor)

**Source:** NASA Exoplanet Archive pscomppars, queried 2026-11-01
(`/viper/ptmp2/jinma/dhj_planets_1101/src/pscomp_err.csv`). Relative errors are the mean of the + and -
1-sigma errors divided by the value. The reference after each mass is the source the archive adopts.
The archive Teq uses its own assumptions; the Seidel Teq error (zero albedo, propagated Teff and a/R*) is in
the Sec. 2 table.

**Filter used:**
- **pass**: an RV (or equivalent) mass with sigma_M/M <= 15 %, R_p <= 5 %, Teq and a/R* <= 5 %, Teff <= 5 %,
  and e consistent with 0;
- **marginal**: one quantity outside these limits but measured;
- **fail**: no measured mass (upper limit only) or disputed basics.

| Planet | Teq (paper) | M_p [MJ], rel. err, source | R_p [RJ], rel. err | a/R*, rel. err | Teq (archive), rel. err | Host Teff, rel. err | [Fe/H] | P [d] | e | Flags | Verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| WASP-76b | 2229 | 0.894, **1.5 %** (Ehrenreich+2020, RV) | 1.854, 4.1 % | 4.08, 1.0 % | 2228, 5.5 % | 6329, 1.0 % | +0.37 +/- 0.05 | 1.80988 | 0 (fixed, circular) | Inflated radius; Teq error the largest of the passes but < 6 % | **PASS** |
| KELT-20b | 2263 | **< 3.38 (upper limit only**, Lund+2017; archive lim flag = 1) | 1.741, 4.1 % | 7.42, 1.8 % | 2262, 3.2 % | 8720, 2.9 % | -0.29 +/- 0.22 (loose) | 3.47411 | 0 (fixed) | **No measured mass**, so g unknown (fast-rotating A host prevents RV). Albedo disputed (0.32 adopted in our setup); [Fe/H] uncertain | **FAIL** |
| WASP-121b | 2358 | 1.17, **3.7 %** (Sing+2024) | 1.742, 0.3 % | 3.784, 0.2 % | 2409, 1.0 % | 6628, 1.0 % | +0.17 +/- 0.05 | 1.274925 | 0.0085 +/- 0.0013 (Pelletier+2026; measured, small) | Best constrained of all. The small measured e has negligible dynamical effect but is formally nonzero | **PASS** |
| WASP-178b | 2398 | 1.66, 7.2 % (Hellier+2019, RV) | 1.81, 5.0 % | 7.17, 2.9 % | 2470, 2.4 % | 9360, 1.6 % | +0.21 +/- 0.16 | 3.34483 | 0 (fixed) | Hot A1 host (UV); R at the 5 % limit | **PASS** (borderline R) |
| TOI-1518b | 2491 | 1.83, **25.7 %** (Simonnin+2025, MAROON-X RV) | 1.878, 2.2 % | 4.109, 0.5 % | 2546, 1.4 % | 7299, 1.4 % | -0.10 +/- 0.12 | 1.902611 | 0.01 **(upper limit**, Cabot+2021) | Mass measured but loose (g +/- 26 %). Earlier work gave only an upper limit (< 2.3 MJ, Cabot+2021) | **MARGINAL** |
| HAT-P-70b | 2559 | 6.78, **no error: an upper limit** (Zhou+2019; archive lim flag 0 but the paper gives a limit) | 1.87, 6.7 % | 5.45, 7.2 % | 2562, 1.9 % | 8450, 7.3 % | -0.06 +/- 0.08 | 2.744325 | 0 (fixed) | **No measured mass**; R, a/R* and Teff all > 5 % | **FAIL** |
| WASP-189b | 2638 | 1.99, 7.5 % (Lendl+2020, RV) | 1.619, 1.3 % | 4.60, 0.6 % | 3353 (*), 0.9 % | 8000, 1.0 % | +0.29 +/- 0.13 | 2.724033 | 0 (fixed) | (*) The archive Teq uses a different redistribution; the paper's zero-albedo value is 2638 +/- 28 K. Polar orbit, which does not matter for synchronous rotation | **PASS** |
| WASP-127b (not in paper) | 1400 (archive) | 0.165, 11.7 % (Seidel+2020) | 1.311, 2.1 % | 7.81, 1.3 % | 1400, 1.7 % | 5620, 1.5 % | -0.19 +/- 0.01 | 4.17806 | 0 (fixed) | Parameters fine, but off the paper's relation (Teq < 1600 K, no Fe I) | **PASS on parameters**, not usable for this trend |

"e = 0 (fixed)" means the discovery fits assumed a circular orbit. No candidate has a measured significant
eccentricity, so none is flagged eccentric.

**Caltech-named planets:**
- WASP-121b: **pass**.
- TOI-1518b: **marginal** (mass at 26 %).
- KELT-20b: **fail** (no mass; disputed albedo).
- WASP-127b: passes on parameters but is **off the relation**.

**Final selection along the trend** (parameter filter applied):

| Teq slot | Pick | v_LOS (paper) [km/s] | Constraint status |
|---|---|---|---|
| low, 2229 K | **WASP-76b** | 6.95 +/- 0.36 | pass |
| mid, 2358 K | **WASP-121b** | 4.75 +/- 0.30 | pass |
| mid, 2398 K | **WASP-178b** | 4.79 +/- 1.19 | pass (borderline R; strong host UV) |
| mid-high, 2491 K | **TOI-1518b** | 3.82 +/- 0.50 | marginal (mass 26 %). Kept because it is the only planet with a well-measured wind between 2400 and 2640 K; run it at the nominal mass and treat g as a +/- 26 % systematic |
| high, 2638 K | **WASP-189b** | 1.05 +/- 1.59 | pass (check the hot-end ck/EOS coverage) |

Excluded: KELT-20b and HAT-P-70b (masses are upper limits only).

**Gap after filtering:** 2491-2638 K. HAT-P-70b would fill it but fails the filter, and WASP-178b sits next
to WASP-121b, so the upper half of the trend rests on TOI-1518b (marginal mass) and WASP-189b (a wind
consistent with zero).

If KELT-20b is still wanted because its setup exists, run it only as a mass-bracketing pair (for example
1 and 3.4 MJ) and do not use it as a trend point.
