# NOTE 2026-10-08 (orion): RSG wind -- radiation cannot hold up or drive the extended atmosphere (Gamma_gas <= 0.25); Betelgeuse needs ~10-20 km/s turbulent/wave support; observed Mdot = dust-free shock-launched tail; dust only beyond ~12-20 R*

Answer to TASK-2026-10-08-orion-rsg-wind.md, steps 1-3 (step 3 as a two-phase gas + dust test, section 4; low-P opacity + non-LTE, section 5; grain growth, section 6; observed densities, section 7; rescaled grains, section 8; Betelgeuse wind budget, section 9; CNO, step 4, not done).
Orion owns this thread (user 10-08). Python on CPU, no AthenaK changes.

Files: `docs/handover/rsg-wind-orion-1008/` -- `scripts/` (rt.py, rsglib.py with MESA reader + ft/mft cases,
sobolev.py, run_sobolev.py, mesa_sob.py), `tables/` (ft_table.md, mesa_table.md, sobolev_tables.md, mesa_sob.md),
`plots/`. Work dirs (npz, logs) on orion: `/orion/ptmp/jinma/rsg_wind_1008/{repro,ft,mesa,sobolev}`.
MESA profiles: `/orion/ptmp/jinma/Arepo/ICs_golden/ic/mesa-profile/{RSGg,golden}/` (orion only; `MESA_DIR` in
rsglib.py). Environment: `module load python-waterboa/2025.06` (system python3 has numpy 1.19, too old),
`CK_DATA=<repo>/docs/handover/rsg-ck-1008/data/`.

## 0. Reproduction of viper's numbers: exact
gmap_table.md, band_thin.md and all 28 non-MARCS rows of col_table.md (hse1/3/10, w5/w6 x v10/v30, both stars,
clamp + extrap, grey-T bracket, F_top/F*, band fluxes) are identical to viper's (`diff` empty). MARCS columns not
run (models not committed, not downloaded).

## 1. Fuller & Tsuna (2024) chromosphere
Fuller & Tsuna 2024, OJAp, arXiv:2405.21049 ("Boil-off of red supergiants"): convective shocks (v_con ~ 10 km/s,
Maxwellian in v_s) support a cool dense chromosphere; time average (eq. 7, prefactor exactly 1, checked):
rho(r) = rho_ph (R/r)^2 exp[-(v_esc/v_con) sqrt(1 - R/r)] for R <= r <= R_d = 5 R; Mdot = 4 pi R^2 rho_ph v_con
exp[-(v_esc/v_con) sqrt(1 - R/R_d)] (eq. 15); dust wind beyond R_d (here v -> 30 km/s, cut at 10 R). FT's 15 Msun
model: v_esc/v_con = 12.5. This puts ~1e3-1e4 x more gas at 1-3 R than viper's hse10 / 1e-5 Msun/yr wind.

On viper's grey hse photospheres (ft8-ft15 = v_esc/v_con 8-15, `tables/ft_table.md`): max Gamma_F 0.013-0.10
(golden16), 0.009-0.035 (Betelgeuse-like); the DENSER the chromosphere the LOWER Gamma_F (bands saturate at 1-2 R,
F_top/F* drops to 0.67-0.70, the 2-4 R layer cools). The chromosphere stays optically thin (tau_R = 1 within 1.02 R*).

## 2. Nine MESA RSG profiles + FT chromosphere (task step 1, `tables/mesa_table.md`)
MESA rho, T up to the surface (ck tau_R ~ 100 at the inner boundary), FT chromosphere on top from the MESA surface
density; v_con "FT way" = MESA conv_vel where tau = P_rad c/(P v_con) (0.98-0.99 R), and the fixed v_esc/v_con = 12.5.
m15lgl5.1 is exactly golden16 (M 15.26, L 1.149e5, R 669.9, Teff 4106).

| model | M | log L | Teff | v_esc/v_con FT | Mdot_FT [Msun/yr] | static GF FT-way / 12.5 | Sobolev A, v/H: FT-way / 12.5 |
|---|---|---|---|---|---|---|---|
| m10lgl4.4a2 | 9.85 | 4.40 | 3748 | 19.2 | 5e-9 | 0.028 / 0.016 | 0.068 / 0.053 |
| m10lgl4.6 | 10.3 | 4.60 | 4237 | 13.5 | 9e-7 | 0.051 / 0.048 | 0.12 / 0.11 |
| m10lgl4.8 | 9.90 | 4.75 | 3993 | 11.7 | 6e-6 | 0.062 / 0.069 | 0.13 / 0.15 |
| m10lgl5.0 | 9.73 | 4.97 | 3757 | 9.5 | 7e-5 | 0.026 / 0.087 | 0.10 / 0.21 |
| m15lgl5.1 | 15.26 | 5.06 | 4106 | 11.4 | 1.4e-5 | 0.078 / 0.090 | 0.18 / 0.20 |
| m20lgl5.3 | 20.5 | 5.25 | 4088 | 11.6 | 1.7e-5 | 0.087 / 0.101 | 0.20 / 0.23 |
| m20lgl5.4 | 19.5 | 5.41 | 3867 | 9.9 | 1.1e-4 | 0.057 / 0.141 | 0.16 / 0.30 |
| m20lgl5.5 | 19.4 | 5.50 | 3776 | 9.1 | 2.8e-4 | 0.030 / 0.151 | 0.13 / **0.34** |
| m26_p455 | 25.0 | 5.36 | 4187 | 11.9 | 1.6e-5 | 0.097 / 0.107 | 0.23 / 0.25 |

(clamp; extrap within ~0.01.) Peak always at 3-4 R, T 1500-1700 K, rho ~1e-15 g/cm^3, kF/kR 60-110; at fixed
v_esc/v_con Gamma grows with L. Mdot_FT is extremely sensitive to v_con: FT-way gives 5e-9 to 3e-4 Msun/yr (log
L >= 5.4 far above observed rates; golden16 ~10 x FT's own estimate).
Caveats: with MESA rho, T the ck kappa_R is below MESA's, ck tau_R(R) = 0.05-0.29, F/F* up to 1.3-1.5 just below R
and F_top/F* = 1.0-1.3 (static Gamma_F biased high by ~10-30 %); Lucy T 360-550 K below MESA T(R) at the join;
lambda iteration last dT/T 0.01-0.07.

## 3. Velocity-gradient (Sobolev) desaturation (task step 2, `tables/sobolev_tables.md`, `tables/mesa_sob.md`)
Each (band, g-point) treated as a line ensemble: tau_S = k_bg rho v_D/|dv/dr|, v_D = sqrt(2kT/(mu m_H) + xi^2)
(mu 20 molecules / 56 Fe; xi 0, 2, 5 km/s), Gamma_Sob = sum_b F_b/F sum_g w_g k_bg (1 - e^-tau_S)/tau_S / kappa_Edd;
flux A = unattenuated f_b(Teff) (upper bound; Gamma_Sob <= Gamma_thin verified exactly), B = attenuated by each
band's window (smallest-g) opacity from the photosphere. Checks: dv/dr -> inf gives Gamma_thin to 4 digits,
dv/dr -> 0 gives 0, thin value matches gmap_table.md.

- beta-law winds (golden16 / Betelgeuse-like, beta 0.5-2, v_inf 10-40 km/s, Mdot 1e-7..1e-4, 1728 cases): max
  0.31 (mu 20) / 0.33 (mu 56) at Mdot 1e-7, v_inf 40, beta 0.5, 1.45 R, 2500 K (97 % the 0.26-0.42 um Fe band).
  Gamma_Sob FALLS with Mdot: at 1e-5 0.08-0.20, at 1e-4 0.03-0.12. At moderate Mdot the peak is on the 1700 K
  plateau at 2.3-3.8 R and already 50-95 % of the thin value: desaturation is nearly complete there, the thin limit
  (0.29 golden16, 0.20 Betelgeuse-like) is the cap.
- viper's winds, FT and MESA columns: Sobolev = 1.5-4 x static. MESA overall max **0.34** (m20lgl5.5, v_esc/v_con
  12.5, |dv/dr| = v_con/H, flux A; 0.26 with flux B), at 3.19 R, 1692 K, 2e-15 g/cm^3; carried by 0.61-0.85 um
  (TiO/VO/K) 0.58, 0.42-0.61 um 0.24. A beta = 1 acceleration to 30 km/s gives a SHALLOWER gradient at 3-4 R than
  v_con/r (6-20 % lower Gamma_Sob), so v_con/H is the upper bracket.
- Approximations (mostly overestimates): g-points are not lines, line overlap/blanketing ignored, radial streaming,
  LTE, equilibrium chemistry, 1x solar (not CNO-processed); the FT gradients v_con/r, v_con/H are proxies.

## 4. Two-phase gas + dust (task step 3, `tables/twophase_tables.md`, `scripts/twophase.py`)
Added 10-08 (user: "is it possible to have a multi-phase gas"). Optically thin radiative equilibrium of gas in the
diluted stellar field with the ck opacities, all roots T(rho or P, r); then a warm + cool phase in pressure balance
on the MESA columns (m15lgl5.1, m20lgl5.3, m20lgl5.5), cool volume filling factor f = 0.01/0.1/0.3, Sobolev gas
force per phase, dust in the cool phase, clump porosity (1 - e^-tau)/tau for clump size l = 0.01 r, 0.1 r.

**(a) Thermal bistability is real at 3-5 R, absent at <= 2 R** (isobaric/Field, P = 1e-4 dyn/cm^2, clamp):

| star | r/R | stable warm / cool T [K] | rho_c/rho_w |
|---|---|---|---|
| golden16 | 1.5 / 2 | single phase (3358 / 2940) | - |
| golden16 | 3 / 4 / 5 | 2744/939, 2579/795, 1928/708 | 5.4, 6.0, 5.0 |
| m20lgl5.5 | 3 / 4 | 2496/894, 1902/763 | 5.2, 4.6 |

Mechanism: the warm phase is heated in the 0.26-0.85 um bands (TiO/VO-type absorbers) and cools at 4.4-8.7 um;
below ~1300 K the visible absorbers leave the equilibrium chemistry, heating collapses and the gas runs to an IR-balanced
cool root at 700-940 K. This is the likely origin of the 1800-2100 K lambda-iteration oscillation in the columns (but
the thin-RE warm root is 2500-2750 K, the columns' plateau ~1700 K). **Caveat (binding):** all wind states lie below
the ck/FastChem table edge (1e-8 bar = 1e-2 dyn/cm^2; wind P ~ 1e-7..1e-3 dyn/cm^2), so the chemistry is taken at
1e-8 bar (extrap changes only the band-mean k: warm roots move <= 200 K, cool roots <= 1 K). The cool-root T and the
instability window are clamp-dependent; at the real pressure these transitions should sit at lower T.

**(b) Two-phase gas without dust lowers the force:** the cool phase holds 34-40 % of the mass (f = 0.1) with gas
Gamma_c ~0.01; mass-weighted gas Gamma <= 0.27 (single phase 0.34).

**(c) Dust in the cool phase -- the grain temperature is the binding criterion.** T_d = Teff W^(1/(4+p)); golden16 at
3 R: 2017 / 1689 / 1256 K for p = 1 (absorbing) / 0 / -1 (nearly transparent iron-free silicate, Al2O3), 4 R: 1793 /
1458 / 1032 K. The gas criterion (T_c < T_cond) is met wherever a cool phase exists but the grains would evaporate.
- Al2O3 (delta <= 1e-4): <= 0.3 cm^2/g of gas vs kappa_Edd 0.8-1.7 -- ~6x short, never drives.
- absorbing silicates (p = 1) and p = 0: do not condense inside 5 R (p = 0 only m20lgl5.5 at 5 R, optimistic T_cond).
- iron-free silicate, p = -1 (delta 4e-3, max supply 12 cm^2/g of gas at kappa_d 3000): the only working case.
  Minimum delta f_cond kappa_d for Gamma_c > 1: 1.7-1.8 cm^2/g gas (golden16, 3-4 R), 1.5 (m20lgl5.3), 0.8-1.0
  (m20lgl5.5); large clumps (l = 0.1 r) cannot reach Gamma_c > 1 at 2 R (golden16, m20lgl5.3) or anywhere inside 4 R
  (m20lgl5.5).

| model | silicate T_cond | Gamma_c > 1 (needs f_cond x kappa_d [cm^2/g dust]) | mass-weighted Gamma > 0.5 |
|---|---|---|---|
| m15lgl5.1, m20lgl5.3 | 1200 K (optimistic) | 4 R only, >= ~900 | 0.50-0.58 at 4 R with f = 0.3 |
| m15lgl5.1, m20lgl5.3 | nominal | no | no |
| m20lgl5.5 | 1200 K | 3 R, >= ~300 | 0.88 at 3 R |
| m20lgl5.5 | nominal | 4 R, >= ~300 | 0.56 at 4 R (f = 0.1) |

kappa_d >= 300-900 cm^2/g(dust) means ~0.1-1 um grains with Q_pr ~ 1 (scattering), the Hoefner (2008) AGB mechanism;
small iron-free grains absorb far too little. The kappa_d scan (300/1000/3000) is an ASSUMED range (literature values
not verified here; geometric limit 2300-7800 x Q_pr cm^2/g for a = 1-0.3 um). With Asplund 2009 Mg+Si the silicate
delta is ~2.5e-3, raising the requirement 1.6x.
Approximations: optically thin RE, single-temperature grains, no drag/drift, no grain-growth kinetics (can grains reach
~0.3 um in the time available?), no gas heating by dust, both phases see the unattenuated stellar flux.

## 5. Below the 1e-8 bar table edge: scaled ck extension (fallback), non-LTE thermal balance, DACE rebuild, physical non-LTE
Added 10-08. Work dirs `/orion/ptmp/jinma/rsg_wind_1008/{lowP,nlte}` (tables and npz there, not committed).

**5a. Low-P ck extension (fallback, `lowP/data/ck/Premixed_1x_g8_11_hiT2_lowP.txt`, P 1e-14..1000 bar, 52 nodes).**
Lines are Doppler-limited at <= 1e-8 bar (Lorentz = Doppler only near 0.3 bar), so below the edge k per gram changes
only through chemistry. Viper's summed carrier-group taper (data/exo_fms_ck/tools/carriers.py) FAILS in pressure:
inside the table (1e-8..1e-5 bar vs 1e-8) median error 0.36 dex, p90 1.95 dex (carriers change: H2O/TiO/CO
dissociate, CN/Fe/Fe+ take over above ~2500 K). Used instead: a non-negative fit k_bg = sum_s c_s(T,b,g) X_s(T,P) over
the 31 premixed species with pyfastchem 4.0.3 local equilibrium condensation (in-table median 0.004 dex), then
k(P < 1e-8) = k(1e-8) x clip(S(P)/S(1e-8), 0, 1). Holdouts: 1 decade p90 0.24 dex, 2 decades p90 0.82 dex; predicted
INCREASES (TiO/Fe evaporation, 700-2500 K) failed validation and are cut, so the table is a LOWER BOUND on k below
1e-8 bar; > 2-3 decades of extrapolation unvalidated. Rows at P >= 1e-8 bitwise identical to the original. Effect:
thin Gamma_F at 1500-2250 K drops 2-10x vs clamp (TiO/VO dissociation), max still 0.38 (2500 K). A proper rebuild from
DACE per-species 1e-8 bar cross sections (Lee et al. 2021's own source) is section 5c.

**5b. Non-LTE thermal balance, simple two-level eps model (`tables/nlte_tables.md`, `scripts/nlte.py`) -- SUPERSEDED by 5d (its "warm phase lost, gas 700-930 K" conclusion is wrong: Fe II metastable heating and rotational cooling were missing).** At the wind densities (rho ~1e-15,
n ~ 3e8 cm^-3, P ~ 1e-10 bar) collisions (C = n q, q 1e-12..1e-10 cm^3/s) are far slower than radiative decay, so the
thermalisation fraction eps = C/(C + A) is ~1e-3..1e-4 for IR vib-rot (A 10-100 /s) and <= 1e-6 for electronic/atomic
lines (A 1e6-1e8 /s). Net heating H = sum_b (eps_b kl_b + kc_b)(W B_b(Teff) - B_b(T)) (H-/CIA continuum eps = 1;
Rayleigh scattering removed). LTE reproduces section 4 exactly (golden16 2745/939, 2580/795, 1928/708 K).
- Only R = eps_el/eps_IR matters (identical results for eps_IR 1e-2 and 1e-4). Pure two-level scattering
  (R <= 1e-4): the warm phase DISAPPEARS; the gas sits on the IR-balanced cool root, within 10 K of LTE: golden16
  1220/927/790/707 K at 2/3/4/5 R, m20lgl5.5 ~1170/888/761/684 K. Bistability needs R >~ 0.5-1 (most absorbed visible
  energy degraded into vibration and quenched, e.g. 2536/931 K at 3 R for R = 0.5); R >= 3 gives a single 2750-3000 K
  phase. The continuum supplies < 0.5 % of the absorbed power at P <= 1e-9 bar and cannot hold the warm phase.
- Physical R (orion's reading): the LTE warm phase is heated mainly in 0.26-0.42 um = atomic Fe lines, which resonance-
  scatter with no vibrational cascade (R ~ 1e-4); only TiO electronic bands could reach R ~ 0.5. So the warm phase is
  most likely lost and the wind gas is cool (700-930 K at 3-5 R), where the gas force is ~0.01.
- Is RE reached? LTE: t_th 40-3000 s << t_dyn (flow r/v 3e7-5e8 s, FT shock interval R/v_con 6e7-1e8 s). Non-LTE cool
  root: t_th/t_dyn ~ 0.01-0.2 at 1e-10 bar, 0.1-20 at 1e-11 bar -> at P <~ 1e-11 bar (and 1e-10 bar for slow
  collisions) T is set by shocks/adiabatic expansion, not RE; at P >= 1e-9 bar RE always holds.
- Force: absorption from the stellar beam with isotropic re-emission carries no net momentum, so eps changes Gamma only
  through T.
- Flags: equilibrium chemistry assumed (H2 formation between the phases ~1e10 s > t_dyn); single two-level eps per band
  (no line-specific A, q, no radiative pumping of IR levels); 1-2 um bands treated as vib-rot although TiO/CN/FeH have
  electronic bands there; lowP table is a lower bound (LTE warm roots within 150 K of clamp).

**5c. DACE per-species rebuild of the ck table below 1e-8 bar (`scripts/premix_lib.py`, `premix_lowP.py`,
`build_lowP.py`, `gate.py`; tables and the 16 GB raw DACE slab in `/orion/ptmp/jinma/rsg_wind_1008/dace/`, not
committed).** Lee et al. 2021 (arXiv:2106.11664) built the Exo-FMS table from HELIOS-K/DACE cross sections, premixed at
high resolution with equilibrium-chemistry VMRs, then sorted into k; DACE's grid stops at 1e-8 bar, hence the table
edge. Rebuild: the 1e-8 bar slab of the 32 premixed species from DACE (dace-query API; line lists as in Lee 2021 Table
1, v1.0; atoms Kurucz, available only >= 2500 K), premixed with pyfastchem 4.0.3 equilibrium-condensation VMRs, sorted
per band, value at the table's 4+4 Gauss g-points (this convention matches; the sub-interval mean does not).
- Gate at 1e-8 bar vs the table, 500-3000 K, b2-b10: median 0.063 dex (method OK), p90 0.88 dex -> FAILED, for three
  specific reasons: (1) SiO: FastChem has SiO ~5e-5 at 1200-2300 K carrying 98-99 % of b0-b3, the table behaves as if
  SiO were absent at 1200-1700 K and ~1e-4 x FastChem at 1900-2900 K (cause unresolved: GGchem Si chemistry, SiO cross
  section, or SiO missing); (2) Fe+: every DACE atom file has a ~2-3e-3 cm^2/g smooth floor (untruncated wings); the
  table keeps it for Fe I but not Fe+ (wing cutoff in Lee's Fe+) -> floor removed; (3) a 500-900 K low-g IR deficit
  (chemistry CH4/HCl or wing treatment; unresolved, irrelevant for the wind). b6-b9 (0.42-1.3 um, TiO/VO) match at p90
  0.12-0.26 dex. With the Fe+ floor removed and SiO fitted per T: median 0.044, p90 0.29 (circular for SiO).
- Two tables (P 1e-14..1000 bar, 52 nodes; rows >= 1e-8 bar = the original, bitwise; new rows = table(1e-8) x
  premix(P)/premix(1e-8), cross sections held at 1e-8 bar = Doppler limit): **A** `Premixed_1x_g8_11_hiT2_dace_lowP_SiOfc.txt`
  (FastChem SiO, physical) and **B** `..._SiOtab.txt` (SiO suppressed as the table implies). A vs B differ only in
  b0-b6 at 1000-2900 K (p90 0.24 dex below 1700 K, 0.83 dex at 1900-2900 K).
- Unlike the scaled fallback (5a, decreases only), the premix shows INCREASES at low P in cool gas: at 1000 K and 1e-10
  bar band k x1.1-28 (Fe gas x36, SiO x90 as they evaporate from condensates); also at 600-1300 K and x1.5 in b7 at
  2500 K. At 2500 K A drops harder than 5a (SiO dissociates), B does not.
- Thin Gamma_F, golden16, rho 1e-17 -> 1e-14:

  | T [K] | clamp at 1e-8 bar | 5a scaled | A | B |
  |---|---|---|---|---|
  | 1000 | 0.018 | 0.018 | 0.13 -> 0.056 | 0.069 -> 0.048 |
  | 1500 | 0.25 | 0.14 -> 0.23 | 0.26 -> 0.24 | 0.25 -> 0.24 |
  | 2000 | 0.19 | 0.02 -> 0.14 | ~0.16 | ~0.16 |
  | 2500 | 0.38 | 0.20 -> 0.38 | 0.014 -> 0.32 | 0.018 -> 0.32 |

  Grid maximum 0.32 (A, B) vs 0.38: the gas-force ceiling drops slightly; cool ~1000 K gas gets 3-7x more (0.05-0.13).
- Caveats: atoms below 2500 K are held at their 2500 K cross sections (the 1000 K Fe increase is uncertain -- exactly
  the T range of the non-LTE cool root); rows above 6100 K are distorted by the carried-over ratio + g-monotonicity
  (32 % of those entries changed, median 0.26 dex; irrelevant for the wind, a per-band ratio would fix it); B's SiO
  factor is circular.

**5d. Non-LTE thermal balance with physical rates (`tables/nlte2_tables.txt`, `scripts/nlte2/`; work dir
`/orion/ptmp/jinma/rsg_wind_1008/nlte2`).** Replaces 5b.
- Model: statistical-equilibrium model atoms in the diluted stellar field J = W B(Teff) with Sobolev escape (dv/dr =
  v/r, v 10/30 km/s): Fe I (248 levels, 7134 lines, Barklem 2018 Fe+H rates), Fe II (121 levels), Ti, Cr, Ca, Mn, Ni,
  Mg, Na, K, Al, Si, O, C (Kurucz gfemq E1+M1+E2); photoionisation/recombination (rough, x3-10); metastable quench by
  H/H2 q_forb 1e-12..1e-10 cm^3/s (no data; controls the atomic heating; S_H of Drawin irrelevant). Molecular cooling:
  Neufeld & Kaufman 1993 non-LTE cooling functions (H2O, CO rot/vib, H2; transcribed by eye from the scan), coupled to
  the stellar field (radiative pumping). TiO/VO electronic heating bracketed (fluorescence fraction x vib quench).
  H-/CIA eps = 1.
- Physics: Fe is photoionised (Fe I 1e-6..1e-2 of Fe), so **Fe II** carries the atomic heating via UV-pumped
  metastables quenched by collisions; non-LTE molecular cooling is ~1e-3 of LTE band cooling, dominated by H2O/CO
  ROTATIONAL lines (they thermalise far better than the vibrational bands).
- Stable roots, nominal (golden16; m20lgl5.5 similar but colder):

  | r/R | P [bar] | stable T [K] | warm-phase heating / cooling |
  |---|---|---|---|
  | 2 | 1e-11..1e-9 | 330-510 cold / 1870-2260 warm / 3390-3490 hot | Fe II quench 64-91 % / CO rot or vib |
  | 3 | 1e-11, 1e-10 | < 150-160 cold / ~1800 warm / 3000-3200 hot | same |
  | 3 | 1e-9 | 247 | - |
  | 4 | 1e-11 | < 150 / 1628 | - |
  | 4-5 | others | < 150 (no root above the grid edge) | cold: H2 / IR-band pumping heats, H2O rot cools |

  m20lgl5.5: warm ~1730-1790 K only at 2 R; 3-5 R cold (~210 K or < 150 K).
- The warm phase is NOT robust: it disappears (except at 2 R) if the stellar UV below 300 nm is 1 % of a blackbody
  (real RSG photospheric UV is far below a blackbody; chromospheric UV, e.g. Mg II h&k, is not modelled), moves by
  > 1000 K within the q_forb range, and appears everywhere at the high TiO-heating bracket.
- RE reached? t_th at the warm roots 3e5-7e7 s vs flow 3e7-5e8 s and shock interval 6e7-1e8 s: RE at 1e-9 bar,
  marginal at 1e-10 bar, not at 1e-11 bar. The cold "< 150 K" roots have t_th 3e6-1e8 s, i.e. they are approached
  but not reached; there the gas T follows shocks and adiabatic expansion.
- Flags: blackbody UV (biggest risk to the warm phase); q_forb unknown; TiO fluorescence fraction unknown; equilibrium
  chemistry/condensation (Fe condenses below ~1300 K, removing the atomic heating from the cool phase -- a reason for the
  bistability); NK93 hand-transcribed; OH omitted, SiO crude; ice formation and adiabatic cooling below 150 K not
  modelled.

## 6. Silicate grain growth in the Fuller & Tsuna chromosphere (`tables/grains_*`, `scripts/grains/`, `plots/grains_*`)
Added 10-08. Work dir `/orion/ptmp/jinma/rsg_wind_1008/grains` (venv with miepython 3.3.0, pyfastchem 4.0.3).
- Model: Mie optics of amorphous Mg2SiO4 (Jaeger et al. 2003 via the Kitzmann & Heng 2018 LX-MIE compilation; variants
  lowk = k x0.1 below 5 um, and k floors 3e-4 / 1e-3 below 8 um standing for ~0.3 % / ~1 % Fe), size-dependent grain
  temperature from radiative equilibrium in the diluted Planck field; growth by Mg accretion (Mg is the key species,
  caps dust/gas at 2.05e-3), da/dt = V_mon alpha (n_Mg vbar/4)/2 [1 - sqrt(T_d/T_g)/S], S from the JANAF ln K of
  forsterite (FastChem logK_condensates) at the GRAIN temperature, key-species depletion; seeds 1 nm, n_seed/n_H
  1e-16..1e-13, switched on where S > 1; ballistic parcels g(1 - Gamma_gas - Gamma_d), Gamma_gas 0/0.1/0.3, T_gas
  400-1700 K; density = FT eq. 7 (or mass-conserving); launch flux Mdot(>v0) = 4 pi R^2 rho_ph v_con exp(-v0/v_con)
  (the exponential tail implied by FT eqs. 7 and 15: gas reaching r was launched at v_esc sqrt(1 - R/r)). 20,736
  trajectories. Checks: Mie vs Bohren & Huffman, blackbody grain T, analytic growth, forsterite S = 1 at 1350 K / 1e-4
  bar (FastChem 1345 K), 1036-1108 K at the wind densities; gas-grain heating negligible.
- Optics: kappa_pr ~ 11 cm^2/g(dust) for a <= 0.01 um, 3040 at 0.1 um, 6410 at 0.3 um (Q_pr 0.82), 2400 at 1 um --
  scattering-dominated, the same for all variants; only 0.1-1 um grains reach the 300-900 cm^2/g section 4 needs.
- Grain temperature (Mie) is much colder than the p = -1 estimate of section 4 for pure forsterite (T_d ~ W^0.49):
  seeds first grow at 1.4 R (lowk) / 2.0 R (pure) / 2.35 R (3e-4) / 3.7 R (1e-3) in golden16, 1.3 / 1.7 / 2.0 / 2.9 R
  in m20lgl5.5. Large grains run hotter (pure 0.3 um grains at 2 R exceed T_cond: growth stalls at 0.05-0.1 um until
  the parcel moves out).
- Growth vs dynamics (alpha = 1): t(0.1 um) = 1e6-1e7 s at 1.5-2 R (FT n_H 4e10-3e11) < free fall 1-5e7 s;
  1e8 s at 3 R, 3e8 s at 4 R, 7e8 s at 5 R (golden16) >> free fall. So grains must start at ~1.5-2 R in the dense
  chromosphere; growth that starts at 3-5 R (Fe-bearing grains) is too slow: golden16 then has no dust-driven wind at
  any alpha, m20lgl5.5 only with alpha = 1. alpha = 0.01 essentially never works; seed number and T_gas matter weakly.
- Dynamics: parcels launched at 20-50 km/s turn round at 1.05-1.7 R; reaching 3-5 R needs 76-83 km/s (golden16) or
  61-67 km/s (m20). Dust-free ballistic escape (v0 > v_esc sqrt(1 - Gamma_gas)) already gives Mdot = 4e-6 / 7e-6 /
  2.6e-5 Msun/yr (golden16, Gamma_gas 0 / 0.1 / 0.3) and 1.1-4.7e-4 (m20). Dust-driven winds (dust-free orbit bound,
  1309 of 20,736) need nearly iron-free grains, alpha 0.1-1, v0 >= 50-70 km/s: golden16 Mdot 6e-6..8e-4 Msun/yr,
  v_inf median 92 km/s; m20lgl5.5 2e-4..7e-3, median 95 km/s. Observed RSGs: 1e-7..1e-5 Msun/yr, 10-40 km/s.
- **Reading:** if pure forsterite formed efficiently at ~2 R in the FT chromosphere, RSG winds would be 1-3 orders of
  magnitude too strong and too fast (dust lowers the launch threshold from ~93 to 60-70 km/s, x30-100 in the
  exponential tail). Observations therefore require inefficient grain formation (Fe content, low sticking, Mg+
  photoionisation: Mg+/Mg ~ 1e3-6e4 in an undiluted-blackbody-UV estimate, shock destruction) and/or a chromosphere
  LESS DENSE than the MESA-based FT profile. Note that FT's own eq. 15 with the MESA v_con already gives 1.4e-5
  Msun/yr for golden16, ~10x FT's quoted ~1e-6 for 15-20 Msun: the density normalisation is exponentially sensitive
  to v_esc/v_con (section 2).
- Flags: drift 10-60 km/s (Gamma_d and v_inf upper limits); no attenuation by the dust (inconsistent at Mdot >= 1e-4);
  the parcel density is a prescription (growth happens at 1.5-2 R where it is highest); the result hinges on k at
  0.4-2 um (a 1 % Fe floor removes the golden16 winds); fall-back parcels regrow grains (no shock sputtering); no
  nucleation calculation.

## 7. Observed RSG atmosphere densities vs the FT chromosphere (`density_survey/`)
Added 10-08. Page with plots (rho [g/cm^3] vs r [cm]): https://claude.ai/artifact/JYF6Z6B9btSFQ6s5a6jrqW (private;
same content as `density_survey/density_check.html`). Literature from web; digitized values in
`density_survey/dent24_fig3_digitized.txt`; working copies of the papers in `/orion/ptmp/jinma/rsg_wind_1008/dens_survey`.
- Best current mean density for Betelgeuse: the Harper semi-empirical radio model as updated by Dent et al. 2024
  (arXiv:2404.06501, Fig. 3, read off; R* = 1014 Rsun at 222 pc): n_H = 8.3e11 / 9.3e10 / 7.6e9 / 1.6e9 / 6.8e8 / 3.7e8
  cm^-3 at 1.2 / 1.5 / 2 / 3 / 4 / 5 R* (rho = n_H x 2.27e-24 g), T peaking ~3800 K at 2 R*. Harper, Brown & Lim 2001
  (ApJ 551, 1073) use the 11.15 um diameter (56 mas = 1.32 x the 42.5 mas photospheric diameter); rescaled to the
  photospheric R* and 222 pc (n_H ~ d^-1/2) it agrees with Dent to within 20 % at 2-5 R*. Other points: Ohnaka et al.
  2009 CO layer ~2e10 cm^-3 at 1.4-1.5 R* (derived, +-5x); Antares MOLsphere (O'Gorman et al. 2020 App. B) n_H ~3e10
  at ~1.3 R* (+-10x). Observed mean outflow inside 5 R* < 5 km/s (quasi-static / turbulent).
- Ratio observed (Dent 2024) / FT: vs golden16 as built 0.34 / 0.22 / 0.36 at 1.5 / 2 / 3 R* (0.47, 0.56 at 4, 5 R*);
  vs FT scaled to Betelgeuse (18 Msun, 1014 Rsun, 222 pc) 0.16 / 0.083 / 0.12; at 168 pc (764 Rsun) 0.43 / 0.28 /
  0.48; Antares MOLsphere ~0.02-0.03. The observed profile falls as ~r^-2.85 beyond ~2 R*, steeper than FT's r^-2.
- Mass loss: Betelgeuse (1-4)e-6 Msun/yr, v_inf ~10-15 km/s; Antares 2e-6; Beasor et al. 2020/2023 prescription at
  log L 5.06, 15 Msun: 2.8-2.9e-6. FT eq. 15 with MESA inputs gives 1.4e-5: 5-10x too high, consistent with the
  density excess.
- **Bottom line:** the MESA-based FT chromosphere is too dense at 1.5-3 R* by ~3-5x (golden16 as built) to ~2-12x
  (scaled to Betelgeuse), largest near 2 R*; uncertainty ~x3 (distance/radius, metal ionisation fraction n_H ~
  x_e^-1/2, epoch variability; radio measures <n_H^2>, so clumping makes the true mean even lower). A 10-25 % lower
  v_con (~6-7 km/s instead of MESA's 8.2) puts FT onto the observations. Rescaled grain runs: section 8; the
  wind budget built on this observed profile: section 9.

## 8. Grain growth with the chromospheric density lowered (`tables/grains_rescale_*`, `scripts/grains_rescale/`)
Added 10-08. 42,240 trajectories: golden16 and m20lgl5.5; density x f_rho in {1, 0.3, 0.1, 0.03, 0.01} (profile and
launch flux) or v_esc/v_con in {12.5, 15, 20} at fixed rho_ph; optics {lowk, pure, Fe3e-4, Fe1e-3}, alpha {0.1, 1},
seeds {1e-15, 1e-13}, Gamma_gas {0, 0.1}, T_gas {400, 1200} K; launch bins of 2.5 km/s weighted by the exponential tail.
- **No dust-driven wind is observed-like.** Dust-driven Mdot is either >= 7.5e-5 (golden16) / 1.3e-5 (m20) Msun/yr with
  v_inf ~60-190 km/s, or exactly 0 once the density is low enough for observed rates (f_rho <= 0.1 or v_esc/v_con >= 15):
  grains still form (1.7-2.7 R pure, 2.9-5.8 R Fe1e-3) but never reach 0.1 um. Grains grow only where the launch
  density is high, and that density sets a high Mdot -- there is no gentle 1e-6 Msun/yr dust wind in this model.
- Every case in the observed box (1e-7..1e-5 Msun/yr, 10-40 km/s) comes from the DUST-FREE ballistic tail (gas launched
  above v_esc sqrt(1 - Gamma_gas)), v_inf 24-36 km/s. golden16 ballistic Mdot (Gamma_gas 0 / 0.1): f = 1: 4.1e-6 / 7.3e-6;
  0.3: 1.2e-6 / 2.2e-6; 0.1: 4.1e-7 / 7.3e-7; q = 12.5: 1.2e-6 / 2.4e-6; q = 15: 8.5e-8 / 1.8e-7.
- Consistency: the density survey (section 7) puts Betelgeuse at 0.22-0.36 of golden16; at f ~0.3 the ballistic tail
  gives 1.2-2.4e-6 Msun/yr = Betelgeuse's observed rate. The mismatch is v_inf: ~31 km/s vs 10-15 (24-27 km/s at best,
  only at v_esc/v_con = 20 with Mdot ~1e-9).
- Low dust-driven v_inf needs a peak Gamma_d just above 1 (median v_inf 25 km/s for Gamma_d < 1.5; 74-133 km/s for 1.5-6);
  those runs are near-escape launches with Gamma = 1 reached late (~4-5 R) and carry ~15 % of the dust-driven mass.
- Flags: drift and dust attenuation still ignored (dust-driven Mdot, v_inf upper limits; the ballistic tail is unaffected);
  the ballistic tail extrapolates FT's exp(-v0/v_con) beyond v_esc; 68 of 42,240 runs timed out near thresholds.

## 9. Betelgeuse wind budget from the observed density (`tables/betelgeuse_tables.md`, `scripts/betelgeuse/`)
Added 10-08. Observed profile = Dent et al. 2024 (n_H digitized, spline over 1.08-5 R*, r^-2.85 beyond; T from six
read-off points, Harper's shape beyond 5 R*), Harper 2001 as a second profile. Set A: 222 pc, R 1014 Rsun, Teff 3650 K,
L 1.64e5 Lsun, M 18/20; set B: 168 pc, R 764 Rsun, L 1e5, M 16.5/19 (n_H x 1.15).
- Required force, steady wind with the observed rho and Mdot 2e-6 (set A, M 18): Gamma_req = 0.97 / 0.90 / 0.92 / 0.93 /
  0.92 / 0.95 / 1.00 / 1.11 at 1.2 / 1.5 / 2 / 3 / 5 / 10 / 20 / 30 R* with no turbulent pressure. Thermal pressure gives
  only 3-10 % (gas scale height 0.007 R* vs observed density scale height 0.13-0.95 R*); inertia <= 0.008. Steady-wind
  v = 0.3 / 1.0 / 1.7 / 3.2 / 4.5 km/s at 2 / 5 / 10 / 20 / 30 R* (observed < 5 km/s inside 5 R*). Gamma_req > 1 only
  beyond 12-20 R* (set and Mdot dependent).
- Available radiation force on the gas (LTE, DACE table A; old table within 30 %): static RE column <= 0.03 inside 5 R*,
  max 0.23 at 9 R*; Sobolev <= 0.002-0.009 inside 5 R*, max 0.10-0.19 at 8 R*; cool-phase brackets (300-2000 K) <= 0.25.
  Non-LTE (5d) would lower it.
- **Missing support:** even with the most generous gas force at each radius, 66-87 % of gravity is unaccounted for at
  1.2-30 R*: an equivalent turbulent/wave pressure v_t = 17.6 / 11.2 / 14.6 / 15.9 / 12.9 / 9.3 / 6.7 / 5.9 km/s at 1.2 /
  1.5 / 2 / 3 / 5 / 10 / 20 / 30 R* (set B 13-22 km/s), i.e. 6-28 x the gas pressure -- the same order as the observed
  MOLsphere / turbulent velocities (+-10-30 km/s, Ohnaka). The Dent profile fits FT eq. 7 with v_esc/v_con = 10.35
  (v_con 8-9.4 km/s, rms 0.11 dex), whose escaping tail gives Mdot 2.9-4.2e-6 (observed) at v_inf ~32-37 km/s.
- Dust: in the slow observed flow, seeds grow from 1.7 R* (pure) / 2.8 R* (Fe1e-3) and reach 0.1 um by 1.8-5 R* (100-1000
  yr transit) -- consistent with the observed inner dust (~1.5 R*, Haubois et al. 2019; Great Dimming dust, Montarges et
  al. 2021). Efficient growth gives Gamma_d 3-6 (would overdrive the wind); Fe-bearing, low-sticking grains <= 0.1-1.6.
  Ballistic parcels: dust-driven winds only for pure/Fe3e-4, alpha 1, seeds 1e-13, Mdot 5e-5..1.7e-4 (excluded).
- **Bottom line:** Betelgeuse's extended atmosphere is held up by non-radiative (turbulent / pulsation-shock / wave)
  pressure of ~10-20 km/s, not by radiation; shock-lifted gas sets the mass budget (the FT tail fitted to the observed
  density reproduces the observed Mdot); dust is needed only beyond ~12-20 R* at Gamma_d ~1, so the grain efficiency must
  be moderate or episodic.
- Flags: T from six points; density beyond 5 R* extrapolated (Harper agrees within 30 %); the steady dust run is
  kinematic (no feedback, drift, attenuation); dust-shell radii from a quick search (papers not opened).

## 10. Verdict
- Line + molecular force on real MESA RSG structures, with a Fuller & Tsuna chromosphere and Sobolev desaturation:
  **Gamma = 0.1-0.34** in LTE (no case of 1,800+ above 0.35; the DACE low-P opacity lowers the ceiling to 0.32),
  peaking at 3-4 R, T 1500-1700 K, mainly 0.61-0.85 um (TiO/VO). It cannot drive RSG mass loss. With physical
  non-LTE rates (5d) the gas beyond ~3-4 R is cold (a few hundred K and falling, not in RE at <= 1e-10 bar), where the
  molecular force is ~0.01-0.1; a warm (1600-2300 K) Fe II-heated phase exists mainly at 2-3 R and depends on the
  stellar UV, so the molecular contribution to driving is small.
- Molecules set up cold, dense, SiO-rich gas, but the grain temperature (radiative) decides: with Mie optics pure
  iron-free forsterite seeds survive from ~2 R (1.7 R for m20lgl5.5), where the FT chromosphere is dense enough for
  growth to 0.1-0.5 um within a free-fall time (section 6); Fe-bearing grains form only at 3-4 R, too late to grow.
- In the FT chromosphere efficient dust formation OVERPRODUCES the wind (Mdot 1e-5..1e-2 Msun/yr, v_inf ~90 km/s vs
  observed 1e-7..1e-5 and 10-40 km/s), and dust-free ballistic escape alone already gives ~4e-6 Msun/yr for golden16.
  So dust plausibly REGULATES rather than enables RSG mass loss, and observed rates constrain the chromospheric density
  and/or the grain formation efficiency. Observed Betelgeuse/Antares densities (section 7) confirm the MESA-based FT
  chromosphere is ~3-12x too dense at 1.5-3 R* (v_con ~6-7 km/s would fit).
- Biggest uncertainties: the chromospheric density (v_con normalisation, FT time-averaging vs observed RSG atmospheres),
  the stellar/chromospheric UV (Fe II heating, Mg photoionisation), grain composition (Fe) and sticking, drift and
  dust attenuation, the Fe II metastable quench rate, non-equilibrium chemistry, the low-P opacity (5c).
- With the observed Betelgeuse density (section 9) the picture is: the extended atmosphere is held up by non-radiative
  pressure of ~10-20 km/s (radiation on the gas gives <= 0.03 inside 5 R*, <= 0.25 anywhere); the shock-lifted
  ballistic tail sets the observed Mdot (section 8: with the density lowered to the observed level, every observed-like
  wind is dust-free, and dust either overdrives or does nothing); dust is needed only beyond ~12-20 R* at Gamma_d ~1. The
  open mismatch is v_inf (~31-37 km/s from the tail vs observed 10-15 km/s).
- Next GPU-side work (non-grey M1/VET binning) is not justified by the gas force. Next analysis: a self-consistent dusty
  wind with drift and attenuation; the origin of the 10-20 km/s support (pulsation/convective shocks vs Alfven waves);
  the supernova-CSM view of the observed profile (M(<r), effective Mdot(r)).
