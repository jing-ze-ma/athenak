# NOTE 2026-10-08 (orion): RSG wind -- line + molecular force with real MESA structures, the Fuller & Tsuna chromosphere and Sobolev desaturation: Gamma <= 0.34, cannot drive

Answer to TASK-2026-10-08-orion-rsg-wind.md, steps 1 and 2 (dust, step 3, and CNO, step 4, not done).
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

## 4. Verdict
Line + molecular radiation force on real MESA RSG structures, with a Fuller & Tsuna chromosphere and Sobolev
desaturation, reaches **Gamma = 0.1-0.34** (no case of 1,800+ above 0.35), peaking at 3-4 R where T passes
1500-1700 K, carried mainly by the 0.61-0.85 um molecular band (TiO/VO). It cannot drive RSG mass loss (Gamma > ~0.5
nowhere), but can reduce the effective gravity by 20-30 % in exactly the layer where dust condenses, which is the
open lever: the dust estimate (step 3) is next. The 1700 K plateau is the bistable molecular cooling regime viper
saw; a two-phase (clumpy, thermally unstable) medium there is a candidate for dust formation in the cool phase.
Next GPU-side work (non-grey M1/VET binning) is not justified by the gas force alone.
