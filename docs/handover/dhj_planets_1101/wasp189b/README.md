# WASP-189b: hydro deep-hot-Jupiter setup package (hot end of the Seidel+2026 trend)

Package: `/viper/ptmp2/jinma/dhj_planets_1101/wasp189b/`. It was built with the same pipeline as
`kelt20_0925` and `dhj_planets_1101/toi1518b`: the face-0 route-B ck RCE, then the face-5 exact RCE, then the transit
radius, the IC, `planet_setup.py`, the W121-class radial grid, the deck and the CPU gates. Caltech runs it;
viper ran only CPU checks. Nothing was committed or pushed.

Seidel et al. 2026 (Nat. Astron., DOI 10.1038/s41550-026-02870-1) put WASP-189b at Teq 2638 +- 28 K with a wind of
v = 1.05 +- 1.59 km/s, which is consistent with zero. Its target in this trend is the slow, strong-drag end.

---------------------------------------------------------------------------------------------------
## 0. TABLE COVERAGE VERDICT (read this first)

**Verdict: OK. No new tables are needed, at any of the three available metallicities (1x, 3x, 10x).**
- **T:** there is no clamping anywhere in the domain.
- **p:** the ck/CE tables clamp at p < 1e-8 bar. That is the top 23 of 124 radial cells (18.5 %), in the region the accuracy rule does not count
  (p < 1e-6 bar), and it is the same clamp every W121 run already has.

Evidence (scripts/table_coverage.py; outputs in gate/table_coverage_*.txt):

| item | table range | WASP-189b need (1-D estimates) | status |
|---|---|---|---|
| ck k-tables `Premixed_{1x,3x,10x}_g8_11_hiT2` | T 100-10100 K (native Exo-FMS data to 6100 K; hiT2 adds 6100-10100 K), p 1e-8 - 1e3 bar; hiT2 equals the native table at every native node (max dlog k = 0) | global-mean RCE (3x): 2340-4534 K in 1e-9..283 bar. **Dayside-mean (2^1/4 Teq): 2600-4611 K. Substellar-flux bound (sqrt2 Teq, 4x forcing): 2948-5225 K** | inside the NATIVE data (< 6100 K) everywhere |
| CE `FastChem_ck_{1x,3x,10x}_int_hiT2` (mu, H2, He, H, e-, H-) | same T and p as the ck tables | same | inside the native range |
| CIA H2-H2 / H2-H / H2-He / He-H | to 3000 / 2500 / 9900 / 10000 K, held flat beyond (Exo-FMS practice) | H2-H2 is held flat above 3000 K, but H2 is already mostly dissociated there (x_H2 3e-2 at 1e-3 bar, 3000 K) | harmless |
| H- bf/ff (John 1988, from FastChem n(H-), n(e-)) | analytic | the dominant continuum above ~3000 K | included |
| general EOS (`general_eos = table`, built at start-up: Saha for H, He, He+; H2 dissociation with rot/vib partition functions; Na/K/Ca/Al/Mg/Fe donors with condensation; ideal, non-degenerate) | deck: log rho -19..0, log T 1.5..6 (381 x 451 nodes); off-table T is clamped, rho is extrapolated linearly | every IC row lies inside it (log rho -15.5..-2.1). The dump used by the python tools stops at log rho -14: 30 IC rows above ~1e-8 bar fall below it, but they are on the deck's own table | OK |
| radiation p in the EOS | off (`eos_radiation = false`) | aT^4/3 / p_gas ~ 1e-6 at 4500 K, 1e-6 bar | OK |

Physics validity at the dayside T and p:
- **H2 dissociation** is THE process in this regime. The table's mu agrees with the FastChem mu at the
  same (p, T), so the EOS and the opacity chemistry share one composition. At 1 bar mu goes 2.38 (2000 K) -> 2.21 (3000 K) -> 1.53
  (4000 K); at 1e-3 bar it is 1.32 at 3000 K, i.e. H2 is gone above the photosphere on the dayside. FastChem gives
  2.225 / 1.538 at the same points (gate/table_coverage_tables.txt).
- **H ionization** sets in only above ~5000 K at 1e-8..1e-6 bar (mu 1.25 at 1e-8 bar, 5000 K). This is the domain top
  under the strongest forcing (up to 5225 K in the 4x bound), and the EOS includes it. The electron fraction at 4000-5000 K
  (3e-4..4e-3) comes from metals and H, and it feeds H-.
- **Not included** (and not needed here): non-ideal and degeneracy terms (rho < 1e-2 g/cc), second ionization of the metals,
  and H- in the EOS (abundance 1e-7, energetically negligible).
- **Band-grid limit (code warning):** above 5572 K more than 1 % of the local Planck function lies blueward of 0.26 um, and
  the 11-band grid folds it into the bluest band. The 1-D estimates stay below this (max 5225 K). A local
  mu0 = 1 substellar column could reach ~5800 K at p < 1e-6 bar, so only the uncounted top could touch it. This is
  monitoring advice, not a blocker. Check the dayside T(p) at the first output.
- **Stellar UV:** the host is an A star, so 6.2 % of its flux lies below 0.26 um (ATLAS9; a blackbody would give 8.0 %). That flux is folded into
  the bluest band (0.26-0.42 um; band 10 in the code's descending order). The bluest band then holds 27.7 % of the stellar flux
  (code start-up line: "band fractions 0.277014 (bluest)").
- **Deep T beyond 6100 K** occurs only in the IC table rows below 1000 bar (2 rows, max 6230 K at 2240 bar). That is
  below the inner wall (deepest ghost 413 bar, 4789 K).

Caveat: these are 1-D estimates. The global and dayside RCE columns are mu0-sphere averages. The 4x run is a bound for
the substellar flux, but it deposits that flux with the averaged mu0 profile, not as a true mu0 = 1 column.

---------------------------------------------------------------------------------------------------
## 1. System parameters (Lendl et al. 2020, A&A 643, A94; the NASA archive default set)

| quantity | value | source |
|---|---|---|
| P_orb = P_rot | 2.7240330 d (235356.5 s) | Lendl+2020 |
| R_p/R_* | 0.07045 (CHEOPS band) | Lendl+2020 |
| R_p | 1.619 +- 0.021 R_J (1.1575e10 cm) | Lendl+2020 |
| M_p | 1.99 +0.16/-0.14 M_J (7.5 %, RV) | Lendl+2020 |
| a/R_* | 4.600 +0.031/-0.025 | Lendl+2020 |
| R_*, M_* | 2.362 +- 0.030 R_sun, 2.030 +- 0.066 M_sun | Lendl+2020 |
| a | 4.600 R_* = 0.05053 AU (7.559e11 cm) | from a/R_* |
| Teff, log g_* | 8000 +- 80 K, 4.00 (from M_*, R_*) | Lendl+2020 |
| host [Fe/H] | +0.29 +- 0.13 | Lendl+2020 (archive) |
| v sin i_* | 93.1 km/s | Anderson+2018 |
| e | 0 (fixed) | Lendl+2020 |
| orbit | polar (irrelevant for a synchronous planet) | Lendl+2020 |
| **Teq** | **2637.5 K = Teff sqrt(R_*/2a)** (zero albedo, full redistribution; the code convention, T_irr = sqrt2 Teq) | = Seidel+2026 2638 +- 28 K. The archive's 3353 K uses a different redistribution |
| g(R_p) | 1881.8 cm/s^2 (log g 3.27) | GM_p/R_p^2 |
| Omega | 2.669647e-5 s^-1 | 2 pi/P |

**Inconsistency, flagged:** Kepler's law with P and a = 4.600 R_* gives M_* + M_p = 2.32 M_sun, against Lendl's 2.03 M_sun
(+14 %); with M_* = 2.03 instead, a = 0.04835 AU (scripts/system.py). Teq and the tide use a = (a/R_*) R_* (the
light-curve value, consistent with Seidel's Teq). The pgen's Roche option derives G M_* from omega and roche_a, so it
would use the 2.32 M_sun value, i.e. a tide 14 % stronger than with the literature M_*. The tide is small either way (sec. 6).

## 2. Metallicity (own section: the user's ranking of evidence)

**(1) WASP-189b's own retrieved abundances (the only planet-specific constraint).** Gandhi et al. 2023, AJ 165, 242
("Retrieval survey of metals in six ultra-hot Jupiters"; ESPRESSO/HARPS high-resolution transmission) gives these log VMRs for WASP-189b
(their Table 2): Fe -4.94 +0.24/-0.26, Mg -4.67 +0.35/-0.36, Ni -5.65, Cr -7.21, Mn -7.72, Ca -8.31, Ti -8.97, V -10.20;
Na < -8.45 (2 sigma). They state that WASP-189b is the only planet in the sample whose Fe lies more than 1 sigma below the stellar
value, by ~0.7 dex. They attribute this to Fe ionisation (it is the hottest planet in the sample) or to a mass-metallicity trend.
- My conversion with Asplund+2009 (Fe 7.50, Mg 7.60) and H/particle ~0.92 (atomic H at ~3900 K) gives
  **[Fe/H]_atm ~ -0.40 +0.24/-0.26**, [Mg/H] ~ -0.23 +- 0.35. With an H2-dominated normalisation the values would be ~0.2 dex higher.
- Ionisation bias: a Saha estimate of the Fe+/Fe ratio at 3900 K, 1e-3 bar, n_e from FastChem gives ~0.3, i.e. +0.1 dex.
  So the corrected Fe is about solar to 0.3 dex subsolar. Ca and Ti are far below solar, which points to ionisation or cold-trapping,
  not to bulk metallicity.
- Vaulato et al. 2025 (A&A, arXiv:2507.21229, NIRPS/HARPS) detect Fe only in the optical and find hydride/Fe ~0.5 dex above
  equilibrium. They give no bulk metallicity. No JWST spectrum of WASP-189b was found.
- So the retrieval favours **~1x solar (or slightly sub-solar)**, roughly 0.3-0.7 dex below the host.

**(2) Mass-metallicity expectation (M_p = 1.99 M_J).**
- Welbanks+2019 (ApJL 887, L20) report that gas-giant Na/K follow the solar-system C trend, log(X/H)/solar = -1.11 log(M/M_J) + 0.38,
  which gives **+0.05 dex (1.1x solar)** at 1.99 M_J. Their H2O fit (-1.09 log M - 0.95) is strongly sub-solar
  (-1.3 dex) and is usually read as oxygen or cloud depletion, not bulk metallicity. The scatter is ~0.3-0.5 dex.
- Thorngren+2016 BULK enrichment, Z_p/Z_* = 9.7 M^-0.45, gives 7x stellar for the interior. That is an upper bound for the envelope,
  not an atmospheric value.

**(3) Host [Fe/H] = +0.29 +- 0.13 (context only).**

**Choice:**

| deck | [M/H] | vs host | vs Gandhi Fe | tables | status |
|---|---|---|---|---|---|
| **nominal `deck/w189_1x.athinput`** | 0 (1x) | -0.29 dex | +0.40 dex (1.5 sigma; ~+0.3 after the ionisation correction) | repo `data/exo_fms_ck` 1x hiT2 (ck + CE), EOS X 0.7381 Y 0.2485 | built and gated |
| **sensitivity `deck/w189_3x.athinput`** | +0.477 (3x) | +0.19 dex | +0.88 dex | CK3X hiT2 = the WASP-121b production set | built and gated |

- 1x is the closest existing table to the planet's own retrieval and to the mass-metallicity expectation.
- 3x is the upper end of the bracket. It is also the metallicity of the WASP-121b, TOI-1518b and WASP-178b runs on the trend, so a
  3x WASP-189b run is the like-for-like trend point.
- **Missing set:** a sub-solar set (0.3x, [M/H] = -0.5) would bracket the retrieval's central value from below. It needs new
  Premixed ck + FastChem CE tables with the hiT2 extension. The EOS needs only keys (eos_xh/yhe, met). A 10x set exists
  (`kelt20_0925/ckdata10`, hiT2, plus a 10x EOS dump `wasp121_0925/met10/dump`) but it lies above every constraint, so I do not recommend it.

Table sets found on viper (all with the hiT2 T extension to 10100 K, p 1e-8..1e3 bar; table_coverage_tables_{1x,3x,10x}.txt):
- 1x: `/viper/u2/jinma/ATHENAK/athenak/data/exo_fms_ck` (= `kelt20_0925/ckdata1`).
- 3x: `/viper/ptmp2/jinma/w121_s10fresh_1027/pkg/ck/CK3X` (= `ck3x_0930/ckdata3`, identical md5s).
- 10x: `/viper/ptmp2/jinma/kelt20_0925/ckdata10`.

## 3. Host SED -> 11-band shortwave flux (`sed/`)

- `sed/make_sed.py 8000 4.00 sw_band_flux_W189_11.txt 0.29`: Castelli & Kurucz 2004 ATLAS9 (STScI CDBS ck04models
  ckp02/ckp05 at 8000/8250 K, downloaded into `sed/atlas9/`), linear in [M/H] between +0.2 and +0.5. The band edges
  are those of sw_band_flux_W121_11.txt; the file is ascending in wavelength, in erg/s/cm^2 at the stellar surface (only its shape is used).
- Integrated F / sigma Teff^4 = 0.9973.
- **Fraction below 0.26 um, folded into the 0.26-0.42 um band (band 10 in the code's order): 6.16 %** (blackbody 7.97 %;
  [M/H] = 0 model 7.13 %, file `sw_band_flux_W189_11_p00.txt`).
- Band shares (model / BB): 0.26-0.42 um 27.7 % / 35.3 %; 0.42-0.61 um 34.1 / 26.5; 0.61-0.85 um 18.2 / 17.8; 0.85-1.32 um 12.8 / 12.9;
  longer bands < 5 % (`sed/make_sed.out`).
- The file is used for both metallicity decks (`ck/CK3X/sw_flux/` and `ck/CK1X/sw_flux/`).
- The strong short UV (< 0.26 um, 6 %) is absorbed high up in the real atmosphere. Here it goes into the 0.26-0.42 um band, whose
  opacity deposits it at ~3.5e-4 bar (transit per-band p_ref; rce/rceA/transit.txt).

## 4. Model choices

- **Bond albedo 0.10 (adopted, not measured).** Deline+2022 (A&A 659, A74, CHEOPS phase curve) give only A_g < 0.48.
  Their 96.5 ppm eclipse is consistent with thermal emission and inefficient redistribution, i.e. a dark dayside.
  0.10 is also what TOI-1518b uses. The RCE and the deck use the same value (`problem/albedo`).
- **T_int = 490.5 K** from the Thorngren+2019 fit inside rce_w121.py, with F = 4 sigma Teq^4. The code applies it as sigma T_int^4 at the ck cut
  (route B).
- **R_p at the computed transit pressure:** CHEOPS passband (0.35-1.10 um, top hat), chord optical depth on the
  final RCE column (`scripts/transit_w189.py` = wasp121_0925/transit.py with a CHEOPS band added).
  - 3x: **p_ref = 1.195e-3 bar** (range 9.1e-4 .. 1.4e-3 bar for T -/+ 300 K).
  - 1x: see sec. 9.
  - planet_setup.py then puts R_p = 1.619 R_J on the terminator area-equivalent radius of that isobar (rot_potential).
- **Rsolver / keys:** the W121 hfresh template unchanged: lhllc (cb_sensor default 1), `<hydro>/am_vector = true`,
  `problem/am_transport = false` (required), cs_seam_flux positive, sponge_top on (1e-7..1e-6 bar, radial), sponge_bottom
  off, ck T4/C2/E4 keys, ck_nquad 2, ck_sph_face 6, rk2, cfl 0.3. Changed keys are marked `# W189:` in each deck.
  `deck/w189_*.json` lists them.

## 5. IC (`ic/ic_w189_{1x,3x}.txt`) and hydrostatic check

- Route-B ck RCE (rce_w121.py, global mean), then the face-5 exact RCE (w121fresh_1001/ic/rce_f5.py, F5_INLINE=1), then make_ic.py:
  isothermal rows to 1e-10 bar at the top and the exact EOS adiabat below the RCE bottom (283 bar) to 2240 bar.
- The RCE and the IC use the simulation's own tables (same ck/CE/CIA files, the 3x or 1x EOS dump, the same band fluxes).
- f5 tolerance RCE_RESTOL = 2e-5 (`finish_f5.sh`). With 1e-7 the residual plateaus at ~4e-6 in the top cells and the
  run does not end; toi1518b made the same change. 3x f5: OLR balance 0.99981 (rce/rceB/f5/log.txt).
- 3x column: 1e-6 bar 4182 K, 1e-4 bar 3658 K, 1e-3 bar 3373 K, 1e-2 bar 2397 K, 0.1-1 bar 2350 K (isothermal deep
  radiative zone), 10 bar 2726 K, 100 bar 3908 K. There is a thermal inversion above 1e-2 bar.
- **Gate (CPU, gate/gate_ic_3x.txt; fresh t = 0 start, C16 angular, the deck's radial grid, ghost zones written):**
  - T vs the IC table: |dT/T| <= 7e-4 (1e-6..1e-3 bar), 3.7e-3 (1e-3..1 bar, the inversion kink), 6.6e-4 (1..300 bar).
  - r(p_ref) over all columns within 2.3e-6 / 5.0e-6 of the tool's pole / equator radii.
  - Hydrostatic residual |dp/dr + rho g|/(rho g): median 3.0e-3 / 4.7e-3 / 3.9e-3 (pole g; same three bands), 95 % <= 2.1e-2.
    This is the centred-difference truncation at ~4 cells/H.
  - Deepest active cell 227 bar / 4394 K (pole), 291 bar / 4552 K (equator).
  - Min density 1.3e-18, top cells 2.2-3.0e-18 g/cc, near dfloor 1e-18 (as W121's raised top).

## 6. Radial grid, tide, a/R_*

- grid112.py carries the W121 hfresh cells/H, n_H(p) = H/dr of every W121 cell, to this planet. x1max = r(p_top) with the same
  p_top = 6.12e-18 bar as W121 (isothermal continuation). It then minimaxes the 8 coefficients + 2 bumps.
- **nx1 = 124 is the production choice**: every pressure band reaches the W121 minimum cells/H at p >= 1e-6 bar.

| band [bar] | W121 nx1 112: min cells/H | W189 nx1 112 | **W189 nx1 124** |
|---|---|---|---|
| 1e-6..1e-4 | 2.98 | 3.24* | **3.26** |
| 1e-4..1e-2 | 4.04 | | **4.29** |
| 1e-2..1 | 4.10 | | **4.34** |
| 1..10 | 3.80 | | **3.96** |
| 10..300 | 4.01 | | **4.06** |

  (*) nx1 112 has max dr/target 1.12, i.e. 12 % coarser than W121 in places. nx1 124 has max 1.005 (grid/g3x/grid.log).
  - 124 + 4 ghosts = 128 cells per x1 line; the single radial MeshBlock (ck one-block constraint) is kept.
  - 3x keys: x1min = ap = 1.090810e10, x1max = 1.567297e10 cm (x1max/ap 1.437; W121 1.868).
- **Tide** (point-mass g; scripts/system.py): 3 (r/a)^3/q (stellar tide + orbital centrifugal, substellar):
  - W189 3x grid: 9.6e-3 at ap, 1.15e-2 at R_p, 2.86e-2 at x1max (2 (r/a)^3/q, the stellar tide alone: 6.4e-3 / 7.7e-3 / 1.9e-2).
  - W189 1x grid: 9.8e-3 at ap, 3.01e-2 at x1max (gate/tide.txt).
  - Spin centrifugal om^2 r/g: 3.7e-3 / 4.4e-3 / 1.1e-2 (already in rot_potential).
  - W121, for comparison: 9.1e-2 / 0.12 / 0.59 (ap / R_p / x1max 2.1e10).
  - **WASP-189b is ~20x less tidally distorted than WASP-121b.** a/R_* = 4.600; ap/r_L1 = 0.228, x1max/r_L1 = 0.327
    (roche_anchor.py).
- **Roche switch for Caltech** (fork/fin51-r-fix c8aac483c, inputs/production/tools/roche_anchor.py, copied to scripts/;
  output gate/roche_anchor_3x.txt; the decks keep roche off):
  ```
  roche_anchor.py --mp 1.99 --rp 1.619 --pref 1.195e-3 --a_au 0.050529 --porb 2.7240330 --teff 8000 \
     --rstar 2.362 --mstar 2.030 --profile ic/ic_w189_3x.txt \
     --eos /viper/ptmp2/jinma/w121prod_0929/met3/dump/eos_table.txt --x1max 1.567297e10
  -> <mesh> x1min = 1.093066e+10  <problem> ap = 1.093066e+10  grav = 2110.0416  omega = 2.669647e-05
     roche = true  roche_a = 7.559031e+11
  ```
  - The p_ref surface is 1.61723 R_J at the pole, 1.62077 at the terminator equator, 1.63179 substellar and 1.63156 antistellar.
  - Inner wall p ranges 194 bar (pole) .. 507 bar (substellar), all on the adiabat extension (T 4294..4935 K).
  - Max Phi on x1max is 0.46 Phi(L1).
  - Switching it on also changes x1min (by +0.21 %). The radial grid keys then need a rerun of grid112.py with the new ap.
  - 1x: sec. 9.

## 7. Cost per rotation and rotations needed (scripts/cost.py; estimate scaled from the measured W121 hfresh)

Scaling: W121 hfresh measured 0.0277 s/cycle, dt 13.79 s, 3.7 min/rot on 1 viper node. This planet's dt is scaled by the ratio of
the radial-CFL estimates on the IC columns, and the cost per cycle by nx1.

| deck | dt (expected) | cycles/rot | cost/rot vs W121 | min/rot, 1 viper node | 300 rot |
|---|---|---|---|---|---|
| **1x (nominal), nx1 124** | 7.85 s | 29964 | 4.15x | 15.4 | 76.8 h |
| 3x, nx1 124 | 7.66 s | 30728 | 4.26x | 15.8 | 78.8 h |
| 3x, nx1 112 (does not meet the cells/H) | 8.72 s | 26994 | 3.38x | 12.5 | 62.5 h |

- dt is smaller than W121's because g is higher (smaller H, so smaller dr) and c_s is higher. The rotation is 2.14x longer.
- Per unit simulated time the cost is 1.99x W121.
- **Rotations:** t_rad = c_p p/(4 sigma g T^3) is shorter than W121's at every level: 0.55x at 100 bar (1.3e6 s = 5.7 rot, against W121's
  22 rot) and 0.47x at 1 bar.
  - 140 W189 rotations = W121's 300-rotation physical time (3.3e7 s).
  - **Recommend 150 rot** (40 h on one viper node; it also covers ~25 deep t_rad at 100 bar). 300 rot is for parity in
    rotation count (79 h).
  - These are estimates from the IC columns, not GPU measurements. The first Caltech segment should measure dt.

## 8. Dynamo field (scripts/dynamo.py, gate/dynamo.txt)

- Method: energy-flux scaling (Christensen+2009; Reiners & Christensen 2010), B_dyn = 4.8e3 G (M L^2/R^7)^(1/6) (solar units),
  with L = 4 pi R_p^2 sigma T_int^4. The surface polar dipole is scaled from Jupiter (8.3 G):
  - **WASP-189b, T_int 490.5 K: B_dyn = 155 G, B_dip,pole ~ 62 G** (47-81 G for T_int 400-600 K).
  - WASP-121b (T_int 567 K) on the same scaling: 166 G / 66 G. The W121 MHD runs used 10 G and 100 G.
- With Teq 2638 K and a dayside > 3500 K, the ionisation fraction is far higher than in W121. Seidel's ~zero wind and their
  Ohmic-drag interpretation make this the planet where **hydro is expected to overestimate the wind**. Read the hydro run as
  the no-drag reference, and plan an MHD (or drag) arm.

## 9. 1x nominal deck (`deck/w189_1x.athinput`): the same chain with W189_MET=1x

- Tables: `ck/CK1X` (symlinks to the repo `data/exo_fms_ck` 1x hiT2 ck + CE, cia, ray; identical cia/ray files) and the 1x EOS
  dump `deepconv_0925/dump` (X 0.7381, Y 0.2485, [M/H] 0).
- Deck keys: `rad_met 0`, `problem/met 0`, `eos_xh 0.7381`, `eos_yhe 0.2485`, `ck_data_dir ck/CK1X`, 1x ck/CE table names.
  Everything else is as in the 3x deck.
- RCE: face 0 (rce/rce1A, rce1B), then f5 (rce/rce1B/f5) seeded from the face-0 checkpoint at it 120 (res 3.0e-6).
  OLR balance 0.99963.
  Column: 1e-6 bar 4192 K, 1e-4 bar 3726 K, 1e-3 bar 3458 K, 1e-2 bar 2574 K, 0.1-1 bar 2340-2380 K, 10 bar 2472 K,
  100 bar 3590 K.
- **p_ref (CHEOPS) = 2.369e-3 bar** (1.90e-3 .. 2.64e-3 bar for T -/+ 300 K). That is 2x deeper than at 3x, because the opacity is lower.
- planet_setup: ap = x1min = 1.095459e10 cm, grav = 2100.8301.
- Grid nx1 124 (grid/g1x, planet column on the 1x EOS):
  - x1max = 1.595318e10 cm.
  - Band minima 3.33 / 4.39 / 4.42 / 4.06 / 4.05 cells/H at p >= 1e-6 bar, all >= W121.
  - Max dr/target 1.017.
- Coverage (gate/table_coverage_1x.txt):
  - IC 2342-5814 K, global-mean column 2342-4214 K.
  - 0 rows above 6100 K; 41 IC rows (20 in the RCE column) at p < 1e-8 bar are p-clamped.
  - Every row is on the deck EOS table.
- CPU gates: parse rc 0; fresh gz start rc 0 (gate/gate_ic_1x.txt):
  - |dT/T| vs the table <= 4.1e-4 / 1.6e-3 / 7.4e-4 (the three bands).
  - r(p_ref) within 1.1e-6 / 5.9e-6 of the tool.
  - Hydrostatic median 3.1e-3 / 4.8e-3 / 4.4e-3 (pole g).
  - Deepest active cell 226 bar / 4073 K.
- Roche (gate/roche_anchor_1x.txt):
  ```
  roche_anchor.py --mp 1.99 --rp 1.619 --pref 2.369e-3 --a_au 0.050529 --porb 2.7240330 --teff 8000 \
     --rstar 2.362 --mstar 2.030 --profile ic/ic_w189_1x.txt \
     --eos /viper/ptmp2/jinma/deepconv_0925/dump/eos_table.txt --x1max 1.595318e10
  -> x1min = ap = 1.097734e+10  grav = 2092.1317  roche = true  roche_a = 7.559031e+11
  ```
  Inner wall 189 bar (pole) .. 541 bar (substellar), 3962..4647 K; x1max/r_L1 0.333.
- Cost: 15.4 min/rot on 1 viper node (sec. 7). t_rad at 100 bar is 5.8 rot (0.56x W121 in seconds).

## 11. GPU smoke (apudev job 12527967, vipa1326, 1 node x 2 MI300A, HSA_XNACK=1, HSA_NO_SCRATCH_RECLAIM=1)

- `smoke/smoke.sh`: the hfresh binary (md5 38c3fa6c...) runs each deck for 40 cycles. Logs are in `smoke/run_{1x,3x}/run.log`.
- **Both decks: rc 0, "Terminating on cycle limit", no FATAL.**
- dt: cycle-0 dt 8.49 s (1x) / 8.32 s (3x); 40-cycle mean dt 8.47 s (1x, t = 338.7 s) / 7.84 s (3x, t = 313.4 s). The
  sec. 7 estimate was 7.85 / 7.66 s.
- This was a smoke test, not a timing measurement. The min/rot in sec. 7 remain scaled estimates.

## 10. Files, commands, md5

- **Decks:** `deck/w189_1x.athinput` (nominal), `deck/w189_3x.athinput` (sensitivity), with the override lists `deck/*.json`.
  Template: `/viper/ptmp2/jinma/w121hfresh_1101/deck/hfresh.athinput` (fork/fin51-r-fix 182477553), via `scripts/make_deck.py` and
  `scripts/build_deck.py`.
- **Run dir:** `ck -> <package>/ck` and `ic -> <package>/ic`. All table paths are relative.
  - `ck/CK3X/{ck,CE_tables,cia,ray}` symlink to `w121_s10fresh_1027/pkg/ck/CK3X`.
  - `ck/CK1X/*` symlink to the repo `data/exo_fms_ck`.
  - `sw_flux/sw_band_flux_W189_11.txt` is a real file in both.
  - When copying to Caltech, dereference the links (`tar -h`) or point them at Caltech's CK3X / 1x copies; the md5s are below.
- **IC:** `ic/ic_w189_1x.txt`, `ic/ic_w189_3x.txt`. RCE runs are in `rce/` (rceA, rceB [3x final], rce1A, rce1B [1x final],
  day_mean, day_sub). Grids are in `grid/g3x`, `grid/g1x` (prelim_rceA is the first pass).
- **Scripts:** `scripts/`, copied from toi1518b and adapted. New ones: table_coverage.py, transit_w189.py, run_rce_day.sh, build_deck.py, dynamo.py,
  roche_anchor.py (fork c8aac483c). Also `chain.sh` and `finish_f5.sh`. **Every command line is in `COMMANDS.sh`.**
- **Gate outputs:** `gate/`: parse_{1x,3x}/parse.log, run_gz_{1x,3x}/, gate_ic_{1x,3x}.txt, table_coverage*.txt, cost_*.txt,
  roche_anchor_{1x,3x}.txt, dynamo.txt, tide.txt.
- **Binaries:**
  - CPU gates: `/viper/ptmp2/jinma/builds/bin/athena_dhj_cpu_8bcd50c8` (fork/fin51-r-fix 8bcd50c89, a descendant of 182477553;
    the extra commits are M1 halo, cs LaunchBounds and the MHD remap, none of which is on this hydro path).
  - GPU smoke: `w121hfresh_1101/bin/athena_dhj_gpu72_18247755_hfresh`.
- **Process note:** another worker's `pkill -f "rce_w121.py B 2 20"` killed my first 1x face-0 launch (rce1B, "Terminated" at
  ~05:24). It was relaunched. Every RCE used for an IC ended with a converged f5 "final OLR" line (3x 0.99981, 1x 0.99963);
  rceA, rceB, rce1A, day_mean and day_sub ended with their TINT line.
- md5 (`md5sums.txt`):
```
77b98d895b9fca20a9c844c552430d57  deck/w189_1x.athinput
915a25bc21320304e640eb89fa1d69e8  deck/w189_3x.athinput
2d8b4572294bc63f46fd1584b6810e5f  deck/w189_1x.json
e280ea25936c666d3c9dd10985631120  deck/w189_3x.json
44a3dbdd978949c2ed61caf7081f90d2  ic/ic_w189_1x.txt
40bf42d360558137c210abd30f6d80fa  ic/ic_w189_3x.txt
15548436fb3c550ad472020642121415  ck/CK3X/sw_flux/sw_band_flux_W189_11.txt
2873c2776b57ee8b42d1ea4285f471af  grid/g3x/gridkeys_124.txt
a24b521554418f744ceff91cfc31f6fb  grid/g1x/gridkeys_124.txt
545397fbccd3865486db3f387d38ccd8  ck/CK3X/ck/Premixed_3x_g8_11_hiT2.txt
a36a7a73c135a4002380f97cc1bb66bc  ck/CK3X/CE_tables/FastChem_ck_3x_int_hiT2.txt
2df0d027b8f42e937091b66402a311a1  ck/CK1X/ck/Premixed_1x_g8_11_hiT2.txt
b6643409062f9ef3b070667c67b62da8  ck/CK1X/CE_tables/FastChem_ck_1x_int_hiT2.txt
6823797cd33e5263741fc6dcf497b633  /viper/ptmp2/jinma/builds/bin/athena_dhj_cpu_8bcd50c8
38c3fa6cc0c27720b5643ee7c7b02bfe  /viper/ptmp2/jinma/w121hfresh_1101/bin/athena_dhj_gpu72_18247755_hfresh
```
