# 3-D RHD massive-star envelope simulations on the Marchant & Bodensteiner (2024) HRD

Directory: `bench/analysis_0916/hrd_3d/`

## Products

| file | what |
| --- | --- |
| `hrd_3d.png` / `hrd_3d.pdf` | the figure: published 3-D RHD envelope simulations over the MB24 Fig. 2 HRD |
| `hrd_3d_nobg.png` / `.pdf` | same points on plain axes (no background image) |
| `hrd_3d_dark.png` / `.pdf` | dark version for black slides: background inverted to dark greys, axes/text white |
| `hrd_3d_dark_transparent.png` | same, with transparent figure margins |
| `sims_3d.csv` | the model table (machine readable) |
| `HR_fig2_axes.png` | background, raw: MB24 Fig. 2 rasterised and cropped to the axes frame |
| `HR_fig2_axes_grey.png` | background of the light version: legend text masked out, desaturated to grey |
| `HR_fig2_axes_dark.png` | background of the dark version: the grey image inverted and contrast-raised, `out = min(255, 2*(255 - g))` |
| `page/index.html` + `page/img/` | the stand-alone page (generated) |
| `make_hrd.py` | the plotting script (calibration constants in its docstring) |
| `prep_background.py` | builds the two background rasters from `HR_fig2_axes.png` (mask, desaturate, fade tracks) |

## Background image and attribution

Background reproduced from **Marchant & Bodensteiner (2024), "The Evolution of Massive
Binary Stars", ARA&A 62, Figure 2** (arXiv:2311.01865). The image is
`figures/HR.pdf` from the arXiv source tarball (`https://arxiv.org/e-print/2311.01865`),
the first `\includegraphics` of the review (`2_binarity_constraints.tex` line 18),
rasterised with ghostscript at 300 dpi (`-sDEVICE=png16m -r300`) and cropped to the
axes frame. Original caption, verbatim from the source:

> Binary fractions and the index of the period distribution $\pi$ for different classes
> of stars and evolutionary stages across the HRD. Stellar evolution tracks of different
> masses are shown, covering evolution from the zero-age MS until core-carbon depletion
> (or until just before the hydrogen envelope is removed in the AGB phase for the
> $5\;\mathrm{M_\odot}$ and $7\;\mathrm{M_\odot}$ models). The dashed line indicates the
> Humphreys-Davidson (HD) limit, above which there is a lack of observed stars
> (Humphreys & Davidson 1979). See text for references.

The review's tracks were computed with MESA at Z = 0.0142; the figure data are archived
at https://doi.org/10.5281/zenodo.10009234.

## Pixel <-> data calibration

Measured on the full-resolution raster (6250 x 4063 px) by locating the major tick
marks inside the frame (columns/rows that are dark over >= 80 % of a 30 px band next to
the spine). Frame spines: left x = 340.0, right x = 6149.5, top y = 406.5,
bottom y = 3656.5 px.

x axis (`log10 Teff`, decreasing to the right), ticks 5.0 / 4.5 / 4.0 / 3.5 at
px 1710.5 / 3080.5 / 4450.5 / 5820.5 (spacing 1370.0 px per 0.5 dex):

    log10 Teff = 5.624270070 - 3.649635036e-4 * px      residuals < 1e-15 dex (exact)

y axis (`log10 L/Lsun`), ticks 6 / 5 / 4 / 3 at py 1056.5 / 1778.5 / 2500.5 / 3223.0
(spacing 722.17 px per dex):

    log10 L = 7.462853890 - 1.38475382e-3 * py          max residual 2.8e-4 dex (0.2 px)

Frame corners therefore map to log Teff = 5.5002 (left) / 3.3799 (right) and
log L = 6.9000 (top) / 2.3995 (bottom).
The cropped background spans full-res pixels x = [343.5, 6145.5], y = [410.5, 3652.5],
i.e. the `extent` used by `make_hrd.py`.


## Background processing (`prep_background.py`)

From `HR_fig2_axes.png` (the 300 dpi crop of `figures/HR.pdf`, axes frame only):

1. **mask** the f_bin legend text block with white, pixels `x = [8, 975]`, `y = [535, 1570]`
   (the mask stops short of the "40 Msun" track label; the axes frame is untouched).
2. **desaturate**: `g = 255 - (255 - luminance)*0.80`, `luminance = 0.299R + 0.587G + 0.114B`.
3. **fade the MESA evolutionary tracks.** Thin dark features are the pixels lying more than 25
   levels below a 41x41 grey closing of `g`; a connected component is called a *track* when its
   bounding box spans more than 250 px, which separates the long curves from the region names, the
   mass labels and the dashed HD-limit segments. Then
   * light version (`HR_fig2_axes_grey.png`): `g -> g + 0.50*(closing - g)` on track pixels, so
     the track median goes **158 -> 186** while region shading (183-228) and the region-name
     glyphs (62-169) are unchanged;
   * dark version (`HR_fig2_axes_dark.png`): `d = min(255, 2*(255 - g))`, then `d -> 0.62*d` on
     track pixels, so the track median goes **192 -> 119 (#777)** while the region shading stays
     at ~55-145 and the region names / HD limit / mass labels stay near white.

## Added region: "stripped He" (ours, not MB24's)

`make_hrd.py` draws one band that is **not** part of Marchant & Bodensteiner (2024): a shaded
**parallelogram tilted along the helium main sequence**, the way MB24's O / B / OBe bands are
tilted along the hydrogen MS. Corners, in (log Teff, log L):

    (4.75, 2.8) - (5.11, 5.3) - (4.81, 5.3) - (4.45, 2.8)

The **hot edge** (4.75, 2.8) -> (5.11, 5.3) follows the **Woosley (2019)** He-star ZAMS
(3 Msun: log L 3.91, log Teff 4.91; 5 Msun: 4.47, 4.98; 8 Msun: 4.90, 5.04 - a slope of ~7 dex in
L per dex in Teff). The **cool edge** is that line offset **0.30 dex cooler**, because stripped
stars that keep a residual hydrogen layer (**Goetberg et al. 2018**) are cooler at fixed L.
Style: facecolor #cdcdcd (light) / #5c5c5c (dark), alpha 0.72, zorder 1 so the tracks and mass
labels still show through; the label "stripped He" is 25 pt serif in #8a8a8a (light) / #d0d0d0
(dark), rotated to run along the band midline (4.60, 2.8) -> (4.96, 5.3) - the angle is computed
from `ax.transData` after `tight_layout`, so it follows the band exactly. The band is drawn only
on the versions that carry the MB24 background.

Why the region exists at all: MB24 has no stripped-helium box (its cWR box only covers
log L > 5.5), so this is the sdO / stripped-He strip below the classical WR region, from the
Goetberg et al. (2018) grid (2-7 Msun, log L 3.0-4.6, Teff 50-100 kK) and the Woosley (2019)
He stars up to ~10 Msun (log L ~5.0, Teff ~125 kK). Our own 3 Msun He star (log L 3.91,
Teff 81 kK) lies inside the band.

The three code names sit in the lower-left corner at log Teff 5.47 and log L 3.40 / 3.00 / 2.60 at
18 pt, clear of the band, and each carries the year of that code's first massive-star-envelope
paper: **Athena++ (2015-)** from Jiang et al. 2015, **MPI-AMRVAC (2022-)** from Moens et al. 2022,
**AREPO (2025-)** from Ma et al. 2025 and **AthenaK (2026-)** for this work. The three published
years were checked against the `paper` column of `sims_3d.csv` (earliest year per code: 2015,
2022, 2025); the list now has four lines at log L 3.55 / 3.20 / 2.85 / 2.50 (step 0.35).

## Build-up versions (dark, for slides)

Three dark figures with **identical framing, axes, ticks and fonts**, so they overlay exactly in a
talk (verified: the axes spines land on the same pixels, rows 29/1345 and columns 205/2270 of the
2300 x 1520 raster, in all three and in `hrd_3d_dark.png`):

| file | content |
| --- | --- |
| `hrd_3d_dark_1_background.png` / `.pdf` | background only: grey inverted HRD, the stripped-He band, axes; no markers, no legend text |
| `hrd_3d_dark_2_codes.png` / `.pdf` | + the Athena++, MPI-AMRVAC and AREPO points and their three legend lines (no AthenaK) |
| `hrd_3d_dark_3_athenak.png` / `.pdf` | + the AthenaK stars and its legend line - identical to `hrd_3d_dark.png` |

Each also exists with transparent figure margins:
`hrd_3d_dark_1_background_transparent.png`, `hrd_3d_dark_2_codes_transparent.png`,
`hrd_3d_dark_3_athenak_transparent.png`. They are produced by the same `draw()` call with a
`codes=[...]` filter (`[]`, the three published codes, then `None` for all), so nothing but the
plotted content differs.

## Conventions used for the plotted (Teff, L)

One point per simulated *star* (a box/patch simulation is plotted at the global Teff and L of the
star it was cut from). **No text labels are drawn on the figure** - the identities are in
`sims_3d.csv` and in the page table. Where a paper gives two effective temperatures, the plotted
value is the simulated **photospheric / tau = 1** temperature (the one the Jiang model names
T<Teff/kK>L<logL> encode); the flux-defined `Teff,F` is kept in the `Teff_alt_K` column.

**Markers.** All markers are filled; there are no hollow markers. Two choices are recorded here
and in the table rather than flagged on the figure: for **M35ZAMS** and **M35MMS** the
**Schultz et al. 2023** values (log L = 5.2 and 5.4) are plotted where Schultz et al. 2022 give 5.0
and 5.2, and for **YSG1L4.7 / YSG2L5.1** the plotted Teff is the **midpoint** of the range those
papers quote (4800-9000 K and 5000-10000 K). Teff for RSG1L4.5 / RSG2L4.9 was computed from the
papers' own published L and Rphot.

**Legend and labels.** The figure has **no legend box** and no per-point labels: the legend is the
three code names set in large bold colour (22 pt) in the freed lower-left corner, and the model
identities live in `sims_3d.csv` / the page table. Colours: Athena++ golden yellow (#f2c400,
circles), MPI-AMRVAC blue (#1f5fd6, squares), AREPO purple (#7d3ac1, diamonds); each marker has a
thin dark edge (#33302a) in the light version and a thin light edge (#f2efe8) in the dark one.
In the **dark** version the blue, purple and red-orange are swapped for more luminous variants
that read on black, **#5b9bff** (MPI-AMRVAC), **#b57bff** (AREPO) and **#ff7a5c** (AthenaK), for
the markers and the code names alike;
the yellow #f2c400 is unchanged.
Marker size 24 pt, tick labels 18 pt, axis labels 22 pt, code names 18 pt, no title.

**Codes.** `AthenaK` (red-orange stars, #e8452c light / #ff7a5c dark, ms 30) carries **our own
three runs**, which are work in progress and not published: the 15 Msun B-star iron-peak box
(prod_w7, running), the 3 Msun stripped-He-star box (box_w8, running) and the planned global
4 Msun He-star presupernova model (in development); their parameters come from the run inputs.
`Athena++` (one colour and marker; it includes the three Jiang et al. 2015 and four
Jiang et al. 2017 local boxes, which used the original Athena with the Jiang+2012 VET radiation
module, the Athena++ predecessor), `MPI-AMRVAC`, `AREPO`. **No *published* AthenaK massive-star
envelope simulation exists; the AthenaK points here are ours, unpublished.**

## Not plotted (in the CSV, `plot = no`)

* `T9L6.2Z0.1`, `T19L6.4Z0.1`, `T19L6.4Z2` (Schultz, Bildsten & Jiang 2020, ApJ 902, 67, Table 1):
  genuine new runs, but metallicity variants that land on top of the corresponding Jiang et al.
  (2018) points.
* `Ma, Pakmor, Justham & de Mink 2025` (AREPO-IDORT, arXiv:2503.16627): L, Teff, R could not be read
  from the accessible text (probably the same 10 Msun model as AREPO-RSG, unconfirmed).
* The four magnetised local boxes `StarB1-StarB4` of Jiang et al. 2017 (ApJ 843, 68) share the
  StarTop stellar parameters and are folded into that point (marked in the tag).

## Considered and excluded

| paper | reason |
| --- | --- |
| **CO5BOLD RSG models** - Chiavassa, Freytag, Masseron & Plez 2011 (A&A 535, A22; 4 models) and Ahmad, Freytag & Hoefner 2023 (A&A 669, A49, Table A.1; 23 models with M >= 5 Msun); Freytag, Steffen & Dorch 2002 (AN 323, 213) | **excluded by choice** (user request, 2026-09-17): global star-in-a-box RSG models from a different code family; they were compiled and are available in the git history of `sims_3d.csv` if wanted |
| **Chen, Jiang, Goodman & Lin 2024, ApJ 974, 106** (arXiv:2408.12017), 50 Msun TAMS star in an AGN disc, `Iso` + 7 accreting runs | **excluded by choice** (user request): a 3-D Athena++ envelope simulation, but its subject is accretion from an AGN disc, and its Teff is not published |
| Goldberg, Jiang & Bildsten 2022, ApJ 933, 164 (arXiv:2206.04134) | supernova explosion hydrodynamics on the existing 3-D RSG envelopes; its Table 1 is the source of Rphot = 796 / 902 Rsun |
| Schultz, Tsang, Bildsten & Jiang 2023, ApJS (arXiv:2209.14772) | SEDONA Monte-Carlo post-processing of existing Athena++ models, no new hydrodynamics |
| Jiang 2023, Galaxies 11, 105 (arXiv:2310.07829) | review; no new simulations |
| Antoni, Jiang & Quataert 2025 (arXiv:2509.16308) | 1-D explosion calculations using 3-D envelopes as initial conditions |
| Sundqvist, Owocki & Puls 2018; Driessen et al. 2022; Debnath et al. 2024, 2025; ud-Doula et al. 2025 | 2-D (the Debnath 2024 "O4" model is the 2-D ancestor of the 3-D Delbroek et al. 2025 run) |
| Poniatowski et al. 2021; Moens et al. 2022 A&A 657 A81 (FLD method paper); Debnath et al. 2026; Schillemans et al. 2026 (LBV) | 1-D, or method/benchmark papers with no 3-D stellar model |
| Delbroek et al. 2026 (arXiv:2604.25391); Gonzalez-Tora et al. 2025 A&A 699 A244 | re-use existing 3-D AMRVAC models, no new simulation |
| Chiavassa et al. 2009, 2010 | CO5BOLD (same model st35gm03n07), excluded with the rest of the CO5BOLD set |
| Ohlmann / Roepke common-envelope AREPO work; GR-Athena++ / AthenaK NR papers; AGB and RGB models | not massive-star envelope RHD |

## Caveats and open points

1. **log L conflict for the two 35 Msun Schultz models.** Schultz, Bildsten & Jiang 2022 (and the
   model names T42L5.0 / T32L5.2) give log L = 5.0 and 5.2; the Schultz, Bildsten & Jiang 2023
   ApJL 951 L42 Table 1 gives 5.2 and 5.4 for what its own footnotes call the same runs. The 2023
   values are plotted (filled markers, like all others) and both are recorded in the CSV.
2. **Teff definition** differs between papers (flux-based vs tau = 1); both are given where published.
   Moens et al. (2022) Teff and Rph are emergent, time- and angle-averaged, with FWHM ~ 20 kK scatter
   in Teff; the Ma et al. (2025) AREPO and Goldberg et al. (2025) YSG values are time averages over
   pulsations.
3. Teff computed by us from a paper's own L and R: RSG1L4.5 (2721 K) and RSG2L4.9 (3230 K). These are
   filled markers - the inputs are published.
4. ADS bibcodes marked in the CSV were partly constructed from verified journal volume/page in the
   standard ADS format rather than looked up individually.

## Full table

| tag | code | paper | ref | star type / phase | M (Msun) | log L/Lsun | Teff (K) | R (Rsun) | plotted | provenance |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| StarDeep | Athena++ | Jiang, Cantiello, Bildsten, Quataert & Blaes 2015, ApJ 813, 74 | 2015ApJ...813...74J / arXiv:1509.05417 | post-MS 40 Msun star; local plane-parallel box at the Fe opacity peak | 40 | 5.952 | 7130 | 431 | yes | Sec. 2 text: 'effective temperature 7.13e3 K, radius 431 Rsun and luminosity 8.96e5 Lsun'; logL = log10(8.96e5). Run with the original Athena plus the Jiang+2012 VET radiation module (the Athena++ predecessor); grouped with Athena++ here. |
| StarMid | Athena++ | Jiang, Cantiello, Bildsten, Quataert & Blaes 2015, ApJ 813, 74 | 2015ApJ...813...74J / arXiv:1509.05417 | evolved 80 Msun star; local box at the Fe opacity peak | 80 | 6.100 | 11700 | 274 | yes | Sec. 2 text: 'effective temperature 1.17e4 K, radius 274 Rsun and luminosity 1.26e6 Lsun'. |
| StarTop (+ StarB1-4) | Athena++ | Jiang et al. 2015, ApJ 813, 74; magnetic runs StarB1-B4 in Jiang et al. 2017, ApJ 843, 68 | 2015ApJ...813...74J / arXiv:1509.05417; 2017ApJ...843...68J / arXiv:1612.06434 | near-ZAMS 80 Msun O star; local box at the Fe opacity peak | 80 | 5.952 | 47500 | 14 | yes | Jiang+2015 Sec. 2: 'effective temperature 4.75e4 K, radius 14 Rsun and luminosity 8.96e5 Lsun'. The four magnetised runs StarB1-StarB4 of Jiang+2017 (Table 2; B_z0/B_y0 = 60/121, 191/382, 382/764, 382/3819 G) use the same StarTop stellar conditions, so they fall on this point. |
| T9L6.2 (M56HG) | Athena++ | Jiang, Cantiello, Bildsten, Quataert, Blaes & Stone 2018, Nature 561, 498 | 2018Natur.561..498J / arXiv:1809.10187 | LBV / Hertzsprung gap (56 Msun core) | 56 (core) + 0.13 (env) | 6.2 | 9060 (alt 15000) | 342.8 | yes | logL and name from Sec. 2 ('log(L/Lsun)=6.2 and Teff=9000K (T9L6.2)'); plotted Teff is the simulated photospheric value, 'photosphere ... at 342.8 Rsun ... averaged radiation temperature 9.06e3 K'. Masses from Schultz+2020 ApJ 902 67 Table 1. Teff_alt = flux-defined Teff,F = 15 kK from Schultz+2023 ApJL 951 L42 Table 1 (model M56HG). |
| T19L6.4 (M80HG) | Athena++ | Jiang et al. 2018, Nature 561, 498 | 2018Natur.561..498J / arXiv:1809.10187 | LBV / Hertzsprung gap (80 Msun core) | 80 (core) + 0.011 | 6.4 | 18700 (alt 26000) | 102 | yes | Sec. 2 for logL/name; 'photosphere at a smaller radius 102 Rsun with a higher effective temperature 1.87e4 K'. Masses from Schultz+2020 Table 1. Teff_alt = Teff,F = 26 kK (Schultz+2023 Table 1, M80HG). |
| T19L6 | Athena++ | Jiang et al. 2018, Nature 561, 498 | 2018Natur.561..498J / arXiv:1809.10187 | Hertzsprung gap (35 Msun core) | 35 (core) | 6.0 | 18900 | 63.7 | yes | Sec. 2: 'one with log(L/Lsun)=6.0 and Teff=19000K (T19L6)'; later 'Teff = 1.89e4 K with time averaged photosphere radius of 63.7 Rsun'. Core mass from Schultz+2020 Table 1. |
| T9L6.2Z0.1 | Athena++ | Schultz, Bildsten & Jiang 2020, ApJ 902, 67 | 2020ApJ...902...67S / arXiv:2009.01238 | as T9L6.2 but Z = 0.1 Zsun | 56 + 0.13 | 6.20 | 9000 | 322.5 (r_ph) | no | Table 1 (M_core, M_env, Teff, logL, r_ph, r_Fe, tau_Fe, Z). Not plotted: coincides with T9L6.2. |
| T19L6.4Z0.1 | Athena++ | Schultz, Bildsten & Jiang 2020, ApJ 902, 67 | 2020ApJ...902...67S / arXiv:2009.01238 | as T19L6.4 but Z = 0.1 Zsun | 80 + 0.011 | 6.40 | 19000 | 102.2 (r_ph) | no | Table 1. Not plotted: coincides with T19L6.4. |
| T19L6.4Z2 | Athena++ | Schultz, Bildsten & Jiang 2020, ApJ 902, 67 | 2020ApJ...902...67S / arXiv:2009.01238 | as T19L6.4 but Z = 2 Zsun | 80 + 0.011 | 6.40 | 19000 | 104.7 (r_ph) | no | Table 1. Not plotted: coincides with T19L6.4. |
| M35ZAMS (T42L5.0) | Athena++ | Schultz, Bildsten & Jiang 2022, ApJL 924, L11; also 2023, ApJL 951, L42 | 2022ApJ...924L..11S / arXiv:2110.13944; 2023ApJ...951L..42S / arXiv:2306.08034 | ZAMS 35 Msun O star | 35 (core) | 5.2 | 42000 (alt 42000) | 7.86 (r_tau=1) | yes | Schultz+2023 ApJL 951 L42 Table 1: Teff,F = Teff(tau=1) = 42 kK, logL = 5.2, Mcore 35. CONFLICT: Schultz+2022 (and the model name T42L5.0) give logL = 5.0. Plotted with the 2023 value; both are recorded here. |
| M35MMS (T32/T35L5.2) | Athena++ | Schultz, Bildsten & Jiang 2022, ApJL 924, L11; also 2023, ApJL 951, L42 | 2022ApJ...924L..11S / arXiv:2110.13944; 2023ApJ...951L..42S | mid-main-sequence 35 Msun O star | 35 (core) | 5.4 | 37000 (alt 34000) | 13.15 (r_tau=1) | yes | Schultz+2023 Table 1: Teff,F = 34 kK, Teff(tau=1) = 37 kK, logL = 5.4. CONFLICT: Schultz+2022 give logL = 5.2 and Teff 32 kK for the same run (footnote e: 'T32L5.2 in Schultz et al. 2022 and T35L5.2 in Schultz et al. 2023'). Plotted with the 2023 value. |
| M13TAMS | Athena++ | Schultz, Bildsten & Jiang 2023, ApJL 951, L42 | 2023ApJ...951L..42S / arXiv:2306.08034 | terminal-age MS 13 Msun B star | 13 (core) + 0.019 | 4.4 | 27000 (alt 25000) | - | yes | Table 1 of that paper (new model): Mcore 13, Menv 0.019, Teff,F = 25 kK, Teff(tau=1) = 27 kK, logL = 4.4. |
| RSG1L4.5 | Athena++ | Goldberg, Jiang & Bildsten 2022, ApJ 929, 156 | 2022ApJ...929..156G / arXiv:2110.03261 | red supergiant (16.4 Msun) | 16.4 (12.8 inside the inner boundary) | 4.5 | 2721 | 800 | yes | Table 1 (m_c = 12.8, M_final = 16.4) and Sec. 4: 'RSG1L4.5 has log(L/Lsun) ~ 4.5 and Rphot ~ 800 Rsun'. Teff computed by us: 5772 K (L/Lsun)^(1/4) (R/Rsun)^(-1/2) = 2721 K. |
| RSG2L4.9 | Athena++ | Goldberg, Jiang & Bildsten 2022, ApJ 929, 156 | 2022ApJ...929..156G / arXiv:2110.03261 | red supergiant (12.9 Msun) | 12.9 (10.79 inside the inner boundary) | 4.9 | 3230 | 900 | yes | Table 1 and Sec. 4: 'RSG2L4.9 has log(L/Lsun) ~ 4.9 and Rphot ~ 900 Rsun'. Teff computed by us as above = 3230 K. |
| YSG1L4.7 | Athena++ | Goldberg, Jiang, Bildsten & Cantiello 2025 (Type IIb progenitors in 3D) | arXiv:2508.12486 / 2025arXiv250812486G | partially stripped yellow supergiant (SN IIb progenitor) | 8.2 (core 8.154, env 0.051) | 4.7 | 6900 | 135-185 (pulsating) | yes | Sec. 2-3 text: log<L>/Lsun = 4.7; 'The Teff of YSG1L4.7 varies from ~4800 K to 9000 K'; R_ph ranges 135-185 Rsun. Plotted Teff = 6900 K, the midpoint of the quoted range (our choice; no single Teff is given). |
| YSG2L5.1 | Athena++ | Goldberg, Jiang, Bildsten & Cantiello 2025 (Type IIb progenitors in 3D) | arXiv:2508.12486 / 2025arXiv250812486G | partially stripped yellow supergiant (SN IIb progenitor) | 11.98 (core 11.77, env 0.21) | 5.1 | 7500 | 120-250 (pulsating) | yes | Sec. 2-3 text: log<L>/Lsun = 5.1; 'YSG2L5.1 ranges from ~5000 K to 10000 K'. Plotted Teff = 7500 K, midpoint of the quoted range (our choice). |
| Moens+22 G1 | MPI-AMRVAC | Moens, Poniatowski, Hennicker, Sundqvist, El Mellah & Kee 2022, A&A 665, A42 | 2022A&A...665A..42M / arXiv:2203.01108 | classical Wolf-Rayet / stripped He star (wind-launching envelope) | 10.0 | 5.33 | 98000 | 1.61 | yes | Table 2 of the paper (3-D grid): M*=10 Msun and log L are inputs (Table 1); <Teff> and <R_ph> are emergent, time- and angle-averaged quantities at <tau_Ross>=2/3, with large scatter (FWHM ~20 kK in Teff). Teff converted from log10<Teff>. |
| Moens+22 G2 | MPI-AMRVAC | Moens, Poniatowski, Hennicker, Sundqvist, El Mellah & Kee 2022, A&A 665, A42 | 2022A&A...665A..42M / arXiv:2203.01108 | classical Wolf-Rayet / stripped He star (wind-launching envelope) | 10.0 | 5.47 | 68000 | 3.97 | yes | Table 2 of the paper (3-D grid): M*=10 Msun and log L are inputs (Table 1); <Teff> and <R_ph> are emergent, time- and angle-averaged quantities at <tau_Ross>=2/3, with large scatter (FWHM ~20 kK in Teff). Teff converted from log10<Teff>. |
| Moens+22 G3 | MPI-AMRVAC | Moens, Poniatowski, Hennicker, Sundqvist, El Mellah & Kee 2022, A&A 665, A42 | 2022A&A...665A..42M / arXiv:2203.01108 | classical Wolf-Rayet / stripped He star (wind-launching envelope) | 10.0 | 5.64 | 78000 | 3.64 | yes | Table 2 of the paper (3-D grid): M*=10 Msun and log L are inputs (Table 1); <Teff> and <R_ph> are emergent, time- and angle-averaged quantities at <tau_Ross>=2/3, with large scatter (FWHM ~20 kK in Teff). Teff converted from log10<Teff>. |
| Moens+22 G4 | MPI-AMRVAC | Moens, Poniatowski, Hennicker, Sundqvist, El Mellah & Kee 2022, A&A 665, A42 | 2022A&A...665A..42M / arXiv:2203.01108 | classical Wolf-Rayet / stripped He star (wind-launching envelope) | 10.0 | 5.74 | 76000 | 4.27 | yes | Table 2 of the paper (3-D grid): M*=10 Msun and log L are inputs (Table 1); <Teff> and <R_ph> are emergent, time- and angle-averaged quantities at <tau_Ross>=2/3, with large scatter (FWHM ~20 kK in Teff). Teff converted from log10<Teff>. |
| Delbroek+25 O4 | MPI-AMRVAC | Delbroek, Sundqvist, Debnath, Moens, Backs, Van der Sijpt, Verhamme & Schillemans 2025, A&A 704, A104 | 2025A&A...704A.104D / arXiv:2507.12057 | early O star, unified 3-D atmosphere + line-driven wind ("box-in-a-star") | 58.3 | 5.79 | 40200 | 16.2 | yes | Table 1: M=58.3 Msun, log L=5.79, <Teff>=40.2 kK, <R*>=16.2 Rsun, log g=3.79. Sect. 2.1: new 3-D runs using the set-up and parameters of the 2-D O4 model of Debnath et al. 2024. Re-used without new simulations by Delbroek et al. 2026 (arXiv:2604.25391). |
| Ma+25 RSG 10 Msun | AREPO | Ma, Justham, Pakmor, Chiavassa, Ryu & de Mink 2025 (AREPO-RSG), arXiv:2510.14875 | 2025arXiv251014875M / arXiv:2510.14875 | red supergiant, core C burning; global full-sphere RHD with dusty wind | 9.78 | 5.037 | 3331 | 990 | yes | Table 1 / Appendix A: L = (1.09 +- 0.24)e5 Lsun (log by the agent), Teff = 3331 +- 376 K, R_Ross = 990 +- 132 Rsun, M = 9.78 Msun. 3-D time averages of a pulsating star. |
| Ma+25 RSG 20 Msun | AREPO | Ma, Justham, Pakmor, Chiavassa, Ryu & de Mink 2025 (AREPO-RSG), arXiv:2510.14875 | 2025arXiv251014875M / arXiv:2510.14875 | red supergiant, core C burning | 19.60 | 5.417 | 3594 | 1317 | yes | Table 1: L = (2.61 +- 0.58)e5 Lsun, Teff = 3594 +- 220 K, R_Ross = 1317 +- 111 Rsun, M = 19.60 Msun. |
| Ma+25 AREPO-IDORT RSG | AREPO | Ma, Pakmor, Justham & de Mink 2025 (AREPO-IDORT), arXiv:2503.16627 | 2025arXiv250316627M / arXiv:2503.16627 | red supergiant, 10 Msun; first global 3-D RHD envelope simulation with AREPO+IDORT | 10 |  |  | - | no | Abstract only; L, Teff, R could not be read from the accessible text. Very likely the same 10 Msun model reported in AREPO-RSG (Ma+2025), but not confirmed. NOT PLOTTED (no L/Teff). |
| AREPO (Ma, in prep.) A | AREPO | Ma et al., in preparation (unpublished) | - | RSG (in preparation) | - | 4.6 | 3981 | - | yes | user-supplied, unpublished; coordinates given directly as (log Teff, log L) = (3.60, 4.6) |
| AREPO (Ma, in prep.) B | AREPO | Ma et al., in preparation (unpublished) | - | RSG (in preparation) | - | 5.2 | 3715 | - | yes | user-supplied, unpublished; coordinates given directly as (log Teff, log L) = (3.57, 5.2) |
| B star FeCZ box (prod_w7) | AthenaK | this work (Ma et al., in preparation) | - | 15 Msun ZAMS B star, iron-opacity-peak box (running) | 15 | 4.28 | 29360 | 5.34 | yes | this work, in progress: AthenaK box run prod_w7; parameters from the run inputs (L = 1.9e4 Lsun -> log L 4.28, Teff 29360 K, R 5.34 Rsun). |
| He star FeCZ box (box_w8) | AthenaK | this work (Ma et al., in preparation) | - | 3 Msun stripped He star, FeCZ box (running) | 3 | 3.91 | 81280 | 0.455 | yes | this work, in progress: AthenaK box run box_w8, initialised on the Woosley (2019) He/2 3 Msun He-star model; parameters from the run inputs. |
| He4 presupernova global (planned) | AthenaK | this work (Ma et al., in preparation) | - | 4 Msun He star at core collapse, global model (in development) | 3.15 | 4.78 | 48978 | 3.41 | yes | this work, in development: global AthenaK model of the Woosley (2019) 4 Msun He star at collapse (M = 3.15 Msun); parameters from the run inputs. |
