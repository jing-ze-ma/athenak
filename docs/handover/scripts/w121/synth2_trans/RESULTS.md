# synth2_trans: WASP-121b transmission Fe I CCF and Na D velocities, rot-300 (2026-09-30)

State: `w1x/bin/dhj.hydro_w.00150.bin` and `w10x/bin/dhj.hydro_w.00150.bin`, t = 3.3046e7 s =
rotation 300.0 (P = 1.27492504 d). Both are read in place. Grid: C32 (24 MeshBlocks of 16x16), nx1 76 (1x) / 74 (10x).

## Method (trans.py, chem.py, ana.py)
- The geometry is the same as synth.py `winds`: dhjcs panel frames, nearest column plus linear interpolation in r,
  180 position angles x 120 impact parameters (rc[0]..rc[-1]), and 801 chord samples.
  Velocities are the rotating-frame wind plus solid-body rotation. Morning = leading limb = +e1.
- Chemistry: FastChem 3 (pyfastchem 3.1.3) with Asplund+2009 abundances. Metals are x10 for the 10x run, with He/H unchanged as in the input.
  Local equilibrium condensation without rainout. Tables of Fe I, Na I, e-, H-, H, H2 and He on a (log T, log p) grid.
- Lines: Kurucz gfall (gf2600.all, gf1100.all; kurucz.harvard.edu/atoms, air wavelengths).
  - Fe mask: the 300 strongest Fe I lines in 380-790 nm (ranked by gf lambda e^{-E/kT} at 2500 K; lines less than 20 km/s apart are excluded).
    The mask actually spans 380-668 nm, because the strong Fe lines are in the blue.
  - Blends: 1242 Fe I lines inside the windows.
  - Na D2/D1.
  - Partition functions: Barklem & Collet 2016 (CDS J/A+A/588/A96).
  - Voigt profile: thermal plus Gamma_rad plus van der Waals (Kurucz gamma_W x (n_H+n_H2+0.4 n_He)(T/1e4)^0.3).
    Profiles use a first-order expansion in a (1 % error at a=0.05) and the exact Faddeeva function for a>0.02.
  - Fe damping uses the medians over the mask lines (log Gr 7.88, log Gw -7.79).
- Continuum: H- bound-free and free-free (John 1988) plus Rayleigh scattering by H2, H and He (own code).
  CIA, TiO/VO and other species' lines are not included.
- Spectral grid: 0.5 km/s pixels (R = 600 000). Fe windows are +-40 km/s around each mask line; Na windows are +-100 km/s.
- Transit: 15 epochs evenly spaced between 1st and 4th contact (planet rotation angle -17.2..+17.2 deg).
  Each chord is weighted by the local stellar intensity (a/R* 3.7844, b 0.10, quadratic limb darkening 0.35/0.28).
- Measurement:
  - The spectra are convolved with an R = 70 000 Gaussian line-spread function (ESPRESSO 4-UT MR mode).
  - Fe: CCF with a weighted binary mask. The weights are the line depths of a windless, non-rotating run of the same state (`*_template`).
    The RV is the centre of a Gaussian-plus-constant fit over +-30 km/s, minus the template's own centre (which is 0.00).
  - Na: the same fit applied to D2 and D1, then averaged.
  - Positive RV means redshift.
- Line-forming pressures: a contribution function per sample, defined as the marginal line absorption dtau_line x e^{-tau}
  (summed over window pixels and mask-weighted for Fe). The table gives its 16/50/84 % points in log10 p [bar].

## Results (out/table_main.txt; old = synth_rot300*/winds.txt, rho ds, 'ray' weighting, alpha 0)

| quantity | old 1x | new 1x | old 10x | new 10x | observed (source) |
|---|---|---|---|---|---|
| Fe RV morning [km/s] | +2.22 | **-0.66** | +3.37 | **+0.48** | -4.12 +- 0.15 (Seidel+2025 ESPRESSO 4-UT) |
| Fe RV evening | -13.21 | **-11.86** | -12.21 | **-11.65** | -6.90 +- 0.11 (Seidel+2025) |
| Fe RV whole transit | -7.32 | **-6.75** | -5.56 | **-6.48** | -5.2 +- 0.5 (Bourrier+2020); ~-6 (Maguire+2023 via Gandhi+2023) |
| Fe RV first / second half | -- | -4.65 / -8.36 | -- | -4.06 / -8.45 | Gandhi+2023: first and second half similar |
| Na D RV morning | +0.83 | **+1.67** | +1.33 | **+2.75** | MERC jet 13.7 +- 6.1 red (Seidel+2025); a jet speed, not a centroid |
| Na D RV evening | -13.04 | **-8.95** | -14.26 | **-7.65** | MERC jet 26.8 +- 7.3 blue (Seidel+2025/2023) |
| Na D RV whole transit | -8.49 | -4.84 | -8.57 | -2.37 | -- |
| Fe line-forming log p [bar] 16/50/84, morning | window -5..-3 | -5.63 / -4.60 / -3.50 | window | -6.27 / -5.29 / -4.20 | Seidel: about -5..-3 |
| Fe, evening | | -5.59 / -4.38 / -3.26 | | -6.19 / -5.14 / -4.02 | |
| Na D, morning | window -6..-4 | -6.87 / -4.51 / -3.47 | window | -7.86 / -6.94 / -4.79 | |
| Na D, evening | | -6.27 / -5.55 / -3.47 | | -7.13 / -6.23 / -4.74 | |
| fraction of the contribution inside the old window (Fe / Na, both limbs) | 1 | 0.56 / 0.43 | 1 | 0.35 / 0.32 | |
| Fe CCF mean line depth / Na D depth (R 70k, per epoch, both limbs) | -- | 1489 / 3042 ppm | -- | 1474 / 2995 ppm | |

Across all transit sets, the fitted RVs change by 0.25 km/s or less when the line-spread function is R = 0 or 140 000 instead of 70 000.

## Domain top (out/top.txt; step 6)
- Top-cell pressure at the terminators (|lat|<30, lon +-90 +-10):
  - 1x: morning median 10^-8.55 bar, evening 10^-5.90 bar (range 10^-7.1..10^-5.05).
  - 10x: morning 10^-8.63 bar, evening 10^-6.42 bar.
  - T at the top is about 3400-4000 K.
- **With an isothermal hydrostatic extension above x1max up to r = 2.2e10 cm (extrapolation: T and v of the top cell):**
  - Contribution from above x1max: 1x Fe 1.9 % and Na 0.7 % (evening), 0.0-0.1 % (morning); 10x Fe 2.9 % evening.
  - Absorption lost by the cut, as a fraction of the evening excess equivalent width: 1x Fe 1.7 %, Na 0.9 %; 10x Fe 2.3 %, Na 0.5 %.
  - RV with the extension (`*_ext`):
    - 1x: morning -0.66, evening -12.01, whole transit -6.97 (Fe); +1.71 / -9.11 (Na).
    - 10x: +0.49 / -11.49 / -6.68 (Fe); +2.75 / -7.74 (Na).
  - The domain is therefore high enough: Na D and Fe cores form below the top at the evening limb.

## Validation
- **7a** (out/valid_7a.txt): static isothermal 2500 K atmosphere, no wind, no rotation.
  - Fe CCF shift -0.002 km/s; Na 0.000 km/s.
  - With solid-body rotation only, the code gives Fe morning/evening +5.222/-5.198 km/s. The analytic value is +5.237/-5.208: the static spectrum
    shifted by Omega b_eff cos(theta) and averaged over position angle, with b_eff weighted by line-centre absorption.
  - Full-disk width sigma: code 8.21, analytic 8.23 km/s.
  - Na D2: code +-5.30, analytic +-5.15. The b_eff proxy is saturated at the Na core, which accounts for this difference.
- **7b** (out/{1x,10x}_rhotest.txt): the new chord sampler with rho ds in fixed windows reproduces the old numbers to 0.01 km/s.
  - 1x, alpha 0: +2.22 / -13.21 / -7.32 (ray) and +4.13 / -11.38 / -4.51 (win).
  - 10x: +3.36 / -12.20 / -5.56.
- **7c** (out/table_sens.txt; 5 epochs; baseline 1x_s5base: morning -0.68, evening -12.09, whole transit -6.99):
  - Fe x0.1: -1.30 / -12.38 / -7.39.
  - Fe x10: +0.25 / -11.73 / -6.55.
  - Extension pressure x0.1: -0.68 / -12.14 / -7.07 (Fe).
  - Extension pressure x1: -0.68 / -12.24 / -7.24.
  - Extension pressure x10: -0.67 / -12.54 / -7.78 (Fe); Na -9.66 evening.
  - Gas-only chemistry (no condensation): -0.68 / -12.07 / -7.01.
- Decomposition (5 epochs):
  - No rotation: Fe -4.89 / -5.97 / -5.46, i.e. winds only.
  - No wind: Fe +4.97 / -5.78 / -1.74, i.e. rotation only.

## Files
Scripts are in `docs/handover/scripts/w121/synth2_trans/`. Outputs (`out/*.npz`, `out/*.txt`) and data (`data/`) are in
`/viper/ptmp2/jinma/w121prod_0929/synth2_trans/`. The runs are listed in `driver.sh`.
