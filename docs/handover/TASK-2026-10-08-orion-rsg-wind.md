# TASK for Orion: RSG wind -- does molecular / line / dust radiation force matter? (CPU, Python, no AthenaK)

## 1. Context (viper step 1, done 10-08)
Science goal (user): RSG wind mass loss, potentially driven by molecules. Viper's 1-D test is in
`docs/handover/rsg-ck-1008/RSG_CK.md` (+ tables, plots, scripts). Verdict: **molecules alone cannot drive.** The
saturated flux-weighted Gamma_F = kappa_F/kappa_Edd peaks at 0.04-0.12 (hydrostatic x1-x10 scale heights, r^-2 winds
1e-6/1e-5 Msun/yr), at T 1450-1950 K, rho 1e-17..2e-15, r 1.1-1.3 R* (hse) or 2-4 R* (winds); thin-limit upper bound
0.2-0.4, mostly from the 0.26-0.42 um band (atomic Fe lines), not H2O/CO. Stars: golden16 (15.26 Msun, log L 5.06,
Teff 4106 K, R 670 Rsun; atmosphere rebuilt grey/Lucy -- no density profile on viper) and a Betelgeuse-like set
(18 Msun, log L 5.10, 3600 K); MARCS s4000_g+0.0 / s3600_g-0.5 as non-grey references.
Open caveats this task should close: (i) real MESA RSG structures instead of Lucy rebuilds, (ii) static RT has no
velocity-gradient desaturation (winds are a lower bound), (iii) dust forms below ~1500 K exactly where Gamma_F
peaks, (iv) 1x solar not CNO-processed, (v) the ck table edge at 1e-8 bar.

## 2. Files on this branch (`rsg-ck-1008`)
`docs/handover/rsg-ck-1008/`: `scripts/` (rsglib.py, ck_lib.py, gmap.py, rt.py, plots.py, golden_page_data.json),
`data/` (ck/Premixed_1x_g8_11_hiT2.txt, CE_tables/FastChem_ck_1x_int_hiT2.txt, cia/, ray/, wavelengths_GCM_11.txt,
rosseland_gs98_x0.7_z0.014.txt). The scripts find data/ relative to themselves (CK_DATA overrides). MARCS models:
download from marcs.astro.uu.se yourself (s4000_g+0.0, s3600_g-0.5; spherical 15 Msun solar) -- not committed.
First reproduce viper's gmap and golden16/Betelgeuse column numbers (RSG_CK.md tables) to ~1 %.

## 3. Steps
1. **Your MESA RSG profiles.** Run the column test (hse x1/x3/x10 + winds) on 3-6 real RSG profiles spanning
   12-25 Msun and early-late RSG (log L 4.7-5.5, Teff 3400-4200 K), using each profile's own rho(r), T(r) up to its
   photosphere and the scripts' atmosphere above. Table: star, Gamma_F peak, where (r, T, rho), kappa_F/kappa_R.
2. **Velocity-gradient desaturation (most important).** For the wind columns, replace the static per-g-point
   transfer by a Sobolev / CAK-type estimate: treat each band's k-distribution as a line-strength distribution and
   compute the force multiplier M(t) with t = kappa_e rho v_th / |dv/dr| (or the per-g Sobolev optical depth
   k_g rho v_th/|dv/dr|), for v(r) beta laws (beta 0.5-2, v_inf 10-40 km/s) and Mdot 1e-7..1e-4 Msun/yr. Report
   Gamma_F with desaturation vs the static value. State the approximation (ck g-points are not individual lines).
3. **Dust, order of magnitude.** At the Gamma_F peak region (T < 1500 K): condensation of Al2O3 / silicates
   (equilibrium or a simple condensation-temperature switch), grey-ish dust opacity per gram of gas
   (kappa_dust ~ 1-10 cm^2/g x condensed fraction; state the source), Gamma_dust vs Gamma_F.
4. **Optional, only if cheap:** sensitivity to CNO processing (C/O 0.3 vs solar) via FastChem abundances of CO, H2O,
   TiO in the bands that matter; and the ck table extrapolation below 1e-8 bar (clamp vs extrap bracket as viper).

## 4. Rules
CPU only, Python; no AthenaK changes, no jobs needed (Orion compute node fine for long scans). Do not push large
data. Report in `docs/handover/NOTE-2026-10-08-orion-rsg-wind.md` on this branch: tables for steps 1-3, plots, and
a verdict: can line + molecular + dust force drive or substantially help (Gamma > ~0.5) RSG mass loss, and where.
