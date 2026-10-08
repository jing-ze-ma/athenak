# Can molecules drive RSG mass loss? 1-D ck feasibility test (10-08)

Python only, no AthenaK changes. Work dir `/viper/ptmp2/jinma/rsg_ck_1008/`.
Scripts: `rsglib.py` (opacity lookups, stars), `gmap.py` (test a), `rt.py` (test b), `plots.py`.

## Stars

| model | M [Msun] | L [Lsun] | log L | R [Rsun] | Teff [K] | log g | kappa_Edd [cm^2/g] | source |
|---|---|---|---|---|---|---|---|---|
| golden16 (MESA) | 15.26 | 1.149e5 | 5.06 | 669.9 | 4106 | -0.03 | 1.734 | Golden RSG Structure page data (`golden_page_data.json`: M, R*, kE_glob; Teff = MESA surface T 4106.3 K = (L/4 pi sigma R^2)^1/4) |
| Betelgeuse-like | 18 | 1.259e5 | 5.10 | 912 | 3600 | -0.23 | 1.868 | set by hand per the brief |

golden16 is a genuine RSG (15 Msun, log L 5.06, R 670 Rsun), at the warm edge (Teff 4106 K).
No MESA rho profile of golden16 exists on viper (the run lives on Fritz/Arepo; the page holds only
T, kappa_R, H_p), so the "MESA atmosphere" is rebuilt as a grey hydrostatic Lucy atmosphere
with the same M, R, L (MESA's own atmosphere is grey too).

MARCS spherical 15 Msun, solar, xi = 5 km/s models (downloaded from marcs.astro.uu.se, `dl_marcs/`):
`s4000_g+0.0_m15` for golden16, `s3600_g-0.5_m15` for Betelgeuse (also `s3600_g+0.0`). These are
the reference: proper non-grey, line-blanketed radiative-equilibrium T(tau), down to tau_R = 1e-5.

## Opacities

Exo-FMS ck `ck/Premixed_1x_g8_11_hiT2.txt` + `CE_tables/FastChem_ck_1x_int_hiT2.txt` (1x solar,
equilibrium condensation, NOT CNO-processed), continuum CIA + Rayleigh (H2, He, H, e-) + H- (John 1988)
from `ck_lib.py` (copy of `data/exo_fms_ck/tools/ck_lib.py`). P(rho, T) iterated with the FastChem mu.
Pressures below the table edge (1e-8 bar, i.e. rho < ~1e-14 at 2000-3000 K: the whole wind region):
* `clamp`: k and composition at 1e-8 bar (upper bound on molecular abundance);
* `extrap`: band-mean k extrapolated in log P from the two lowest nodes, slope clipped >= 0 (lower bound).

**Check against MARCS** (same T, rho): ck kappa_R is 20-35 % below MARCS's own kappa_R between
tau_R 1e-3 and 1, within 5 % at tau_R = 1e-5, and 13-19 % above it at tau_R = 60 (T 7800-8100 K); the OPLIB+AESOPUS table that the
AthenaK red_giant runs use (`data/stellar_opac/rosseland_gs98_x0.7_z0.014.txt`) agrees with ck
to within 20 % at the same points. So the ck table is a valid RSG photospheric opacity at the Rosseland level.

## (a) Thin-limit Gamma map (upper bound)

Gamma_F = sum_b f_b(Teff) sum_g w_g k_bg / kappa_Edd, with f_b the B_nu(Teff) band fractions. Full
tables in `gmap_table.md`, figure `gmap.png`. golden16, clamp:

| T [K] \ log rho | -17 | -16 | -15 | -14 | -13 | -12 | -11 | -10 | -9 |
|---|---|---|---|---|---|---|---|---|---|
| 1500 | 0.25 | 0.25 | 0.25 | 0.25 | 0.25 | 0.25 | 0.25 | 0.26 | 0.099 |
| 2000 | 0.19 | 0.19 | 0.19 | 0.19 | 0.19 | 0.28 | 0.35 | 0.38 | 0.39 |
| 2500 | 0.38 | 0.38 | 0.38 | 0.38 | 0.39 | 0.38 | 0.37 | 0.37 | 0.47 |
| 3000 | 0.059 | 0.059 | 0.059 | 0.059 | 0.086 | 0.29 | 0.57 | 0.66 | 0.66 |
| 3500 | 0.01 | 0.01 | 0.01 | 0.01 | 0.013 | 0.032 | 0.13 | 0.47 | 0.9 |
| 4000 | 0.011 | 0.011 | 0.011 | 0.011 | 0.012 | 0.015 | 0.027 | 0.084 | 0.36 |
| 5000 | 0.029 | 0.029 | 0.029 | 0.029 | 0.026 | 0.025 | 0.026 | 0.035 | 0.059 |
| 6000 | 0.14 | 0.14 | 0.14 | 0.14 | 0.094 | 0.065 | 0.056 | 0.059 | 0.083 |

Same grid, other means (golden16 kE, clamp):

| T [K] | Gamma_P (gas T) rho 1e-15 / 1e-11 / 1e-9 | Gamma_R ck | Gamma_R OPLIB+AESOPUS | Gamma_F thin Betelgeuse |
|---|---|---|---|---|
| 1500 | 0.22 / 0.22 / 0.22 | 9e-4 / 9e-4 / 1e-3 | 5e-3 / 0.028 / 0.037 | 0.18 / 0.19 / 0.074 |
| 2500 | 0.035 / 0.037 / 0.082 | 9e-6 / 1.4e-5 / 1.9e-3 | 2e-5 / 1.7e-5 / 3.7e-3 | 0.16 / 0.16 / 0.24 |
| 3500 | 0.004 / 0.05 / 0.34 | 2e-5 / 4.5e-5 / 4e-4 | 7e-5 / 6e-5 / 4e-4 | 0.004 / 0.055 / 0.38 |

* Thin-limit Gamma_F never reaches 1 in the molecular layer: max 0.9 at 3500 K, 1e-9 g/cm^3 (the
  photosphere itself), 0.2-0.4 at 1500-2500 K, < 0.06 at 3000-4000 K below 1e-13 g/cm^3.
  Betelgeuse set (cooler spectrum, kE 1.87): max 0.38. `extrap` lowers the low-rho 2000-3500 K
  values by up to 100x (3000 K, 1e-17: 0.059 -> 7e-4).
* kappa_F(thin)/kappa_R = 1e3-1e4: grey Rosseland M1 underestimates the thin-limit force by 3-4 dex.
* Band decomposition (`band_thin.md`): the thin Gamma_F is NOT molecular for T >= 2000 K. The
  0.26-0.42 um band (atomic Fe/Fe II lines, 3 % of the flux) carries 0.14-0.88 of it; H2O/CO bands
  (2.9 and 5.8 um) give only 0.01-0.02. At 1500 K the 0.61-0.85 um band (TiO/VO/K) gives 0.12 of 0.25.

## (b) Saturated columns (static spherical RT, 88 ck problems)

Method: short-characteristic formal solution along p-z rays (16 core rays + one tangent ray per
shell), LTE pure absorption S = B_band(T) (Rayleigh/e- scattering counted as absorption), diffusion
inner boundary at tau_R ~ 100 (T < 9800 K). Structures: grey Lucy T^4 = Teff^4 (W(r) + 3/4 tau'),
hydrostatic gas pressure with scale height x1, x3, x10 (`hse1/3/10`); winds rho = max(rho_hse,
Mdot/4 pi r^2 v) to 10 R* for Mdot 1e-6/1e-5 Msun/yr at v 10/30 km/s (`w6v10` ...); MARCS T, rho as
given (no iteration). Then 14 damped (geometric-mean) lambda iterations to ck radiative equilibrium
in shells with tau_R < 3, rho(r) held fixed. The lambda iteration does not move T at tau_R > 1e-5
(T at tau_R 1e-3: MARCS 3005 K vs our 3265 K) and oscillates in the 1800-2100 K bistable cooling
regime (last-iteration max dT/T 0.1-0.3 in hse, 0.01-0.07 in winds); a grey-T run (no iteration,
`*_grey.npz`) is given as the bracket. Emergent flux F_top/F* = 0.80-0.83 (grey structures; ck
kappa_R gives a 15-20 % flux deficit), 1.08/1.20 for MARCS models placed at our R; Gamma_F below
uses kappa_F (flux-normalised) / kappa_Edd.

Max Gamma_F above the photosphere (tau_R < 1) (`col_table.md`):

| star | case | Gamma_F max (clamp / extrap) | where (clamp): r/R*, T, rho | kF/kR there | thin Gamma_F there | grey-T Gamma_F max |
|---|---|---|---|---|---|---|
| golden16 | MARCS | 0.005 / 0.005 | top 1.031, 2496 K, 3.8e-12 | 494 | 0.37 | (RE model) |
| golden16 | hse1 | 0.037 / 0.039 | 1.09, 1954 K, 1e-16 | 1.6e3 | 0.21 | 0.009 |
| golden16 | hse3 | 0.080 / 0.073 | 1.33, 1780 K, 8e-18 | 274 | 0.28 | 0.015 |
| golden16 | hse10 | 0.118 / 0.110 | 1.80, 1715 K, 2.4e-17 | 125 | 0.29 | 0.050 |
| golden16 | w6v10 | 0.113 / 0.053 | 3.08, 1509 K, 2.4e-16 | 122 | 0.25 | 0.110 |
| golden16 | w6v30 | 0.085 / 0.082 | 2.89, 1895 K, 9e-17 | 2.2e3 | 0.28 | 0.119 |
| golden16 | w5v10 | 0.101 / 0.105 | 3.58, 1670 K, 1.8e-15 | 86 | 0.29 | 0.072 |
| golden16 | w5v30 | 0.109 / 0.108 | 3.30, 1662 K, 7e-16 | 93 | 0.28 | 0.091 |
| Betelgeuse | MARCS | 0.006 / 0.006 | top 1.073, 2355 K, 1.5e-12 | 655 | 0.13 | (RE model) |
| Betelgeuse | hse1 | 0.046 / 0.044 | 1.10, 1820 K, 1.4e-17 | 349 | 0.18 | 0.016 |
| Betelgeuse | hse3 | 0.041 / 0.063 | 1.32, 1821 K, 2.3e-17 | 320 | 0.18 | 0.012 |
| Betelgeuse | hse10 | 0.109 / 0.066 | 1.93, 1501 K, 1.2e-17 | 128 | 0.19 | 0.104 |
| Betelgeuse | w6v10 | 0.076 / 0.038 | 2.21, 1436 K, 2.5e-16 | 118 | 0.17 | 0.082 |
| Betelgeuse | w6v30 | 0.078 / 0.036 | 2.00, 1490 K, 1.0e-16 | 97 | 0.18 | 0.090 |
| Betelgeuse | w5v10 | 0.064 / 0.068 | 2.37, 1699 K, 2.2e-15 | 56 | 0.20 | 0.060 |
| Betelgeuse | w5v30 | 0.072 / 0.074 | 2.41, 1502 K, 7.2e-16 | 85 | 0.19 | 0.066 |

In the hydrostatic photosphere (tau_R 1e-5..0.1) Gamma_F = 2e-4 to 1e-3 in BOTH our grey HSE and the
MARCS RE models (at equal tau_R within 2x). Saturation cuts the thin-limit value by 100-1000x there:
the strong-line g-points are optically thick, their source function is nearly isothermal, so they
carry almost no net flux; flux leaks out in windows (kappa_F/kappa_R only 1-100 at tau_R > 1e-5).

Emergent band flux fractions at the top (`col_*.png` lower right):

| case | 38 | 12 | 5.8 | 3.9 | 2.9 | 2.2 | 1.6 | 1.0 | 0.71 | 0.50 | 0.32 um |
|---|---|---|---|---|---|---|---|---|---|---|---|
| B_nu(4106 K) | 0.000 | 0.003 | 0.016 | 0.016 | 0.046 | 0.052 | 0.185 | 0.302 | 0.218 | 0.131 | 0.031 |
| golden16 MARCS | 0.000 | 0.002 | 0.014 | 0.016 | 0.050 | 0.063 | 0.268 | 0.254 | 0.177 | 0.128 | 0.027 |
| golden16 hse1 | 0.000 | 0.003 | 0.017 | 0.019 | 0.060 | 0.075 | 0.304 | 0.254 | 0.158 | 0.096 | 0.015 |
| golden16 w5v10 | 0.000 | 0.003 | 0.021 | 0.020 | 0.062 | 0.076 | 0.302 | 0.255 | 0.153 | 0.094 | 0.014 |
| B_nu(3600 K) | 0.000 | 0.004 | 0.023 | 0.022 | 0.061 | 0.067 | 0.222 | 0.314 | 0.187 | 0.087 | 0.014 |
| Betelgeuse MARCS | 0.000 | 0.003 | 0.019 | 0.021 | 0.066 | 0.083 | 0.335 | 0.257 | 0.131 | 0.075 | 0.008 |
| Betelgeuse w5v10 | 0.001 | 0.005 | 0.028 | 0.027 | 0.083 | 0.106 | 0.398 | 0.235 | 0.084 | 0.030 | 0.003 |

Flux is pushed out of the optical/UV (where the strong lines are) into the 1.3-2 um H- opacity minimum.

## (c) Verdict

* **Gamma_F peaks at 0.04-0.12** in every column (both stars, both low-P treatments, grey or
  iterated T), at T = 1450-1950 K, rho 1e-17..2e-15 g/cm^3, r = 1.1-1.3 R* (hydrostatic, x1-x3),
  1.8-1.9 R* (x10) and 2-4 R* (winds). It sits on the edge of the molecular cooling cliff
  (T drops from ~2000 K to ~1300 K), carried by the 0.61-0.85 um and 0.26-0.42 um bands.
* In the actual photosphere and the MARCS RE atmospheres Gamma_F is 1e-3 or less.
* **Molecules alone cannot drive an RSG wind** (Gamma_F > 1 nowhere; even the unsaturated
  thin-limit upper bound stays below 1 outside the photosphere). They help only marginally:
  ~0.1 is in the "not" band (< 0.1) to the low edge of "help". A 10 % cut in effective gravity
  raises the scale height by 1.1x, which matters only together with another levitation mechanism
  (pulsation/convection shocks, turbulence) that already lifts gas to 2-4 R*.
* For comparison, Mdot = 1e-5 Msun/yr at 10 km/s needs only Mdot v c / L = 0.04 of the photon
  momentum, and G M Mdot / R* / L = 0.06, so energy and momentum are not limiting: Gamma < 1 is.
* For M1/VET: kappa_R underestimates the molecular-layer force by 10-1000x (kF/kR 56-2e3 at the peaks).
  The thin-limit Planck/flux mean overestimates it by 2-1000x (photosphere 100-1000x, peak 2-5x).
  A grey run with kappa_R gives Gamma ~ 1e-5..1e-3 there; the saturated truth is 1e-3..0.1.

Caveats: 1x solar premixed (CNO-processed C/O and N would change CO/H2O/CN; not tested); no dust
(Al2O3/silicates form below ~1500 K, exactly where Gamma_F peaks: the dust question is separate and
not answered here); LTE, equilibrium chemistry, pure absorption (no scattering of the line
radiation, no velocity/Sobolev desaturation: in a wind with dv/dr the strong g-points desaturate,
so the static result is a lower bound for the wind case, between this and the thin-limit bound
0.2-0.4); table edge 1e-8 bar (clamp vs extrap bracket shown; extrap is up to 2x lower); ck
optical band has little TiO at T >= 3000 K (0.61-0.85 um mean line k drops from ~1 at 2000 K to
3e-3 at 3000 K, 1e-4 bar), while TiO defines M-type photospheres: the optical contribution near
the photosphere may be underestimated (thin-limit effect at most ~0.1 in Gamma); lambda iteration
not converged in thick bands (bracketed by grey T and by the MARCS RE models, which agree to 2x).

Files: `RSG_CK.md`, `gmap.png`, `gmap_table.md`, `band_thin.md`, `col_golden16.png`,
`col_betelgeuse.png`, `col_table.md`, `col_*.npz`, `runs.log`, `runs_grey.log`,
`dl_marcs/*.mod` (plus `dl_marcs/marcs_st_mod.tar`, 150 MB, deletable).
