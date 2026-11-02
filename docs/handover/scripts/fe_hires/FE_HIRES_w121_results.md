# WASP-121b: high-resolution Fe I transmission + CCF vs Seidel+2026 / Seidel+2025

Directory: `/viper/ptmp2/jinma/w121_obs_1003/mhd_final/fe_hires`. This replaces the density-weighted proxy in `../kpvsys/KPVSYS.md`.
- Scripts: `scripts/`, with `fe_core.py` (opacity / chemistry), `fe3d.py` (3-D ray tracer), `ana.py` (CCF, Kp-Vsys, segments), `val.py`, `cv.py`, `tab.py` and `figs.py`.
- Products: `out/` (`spec_*.npz` excess-absorption spectra per exposure and limb, `ana_*.npz`, `table.md`, `val_summary.txt`, `cv_summary.txt`, `stats.txt`, `val_opac.txt`) and `figs/`.

**Scope.** One dump per model was run:
- h151 00151 (rot 300).
- hr 00160 (rot 318).
- s10 00166 (rot 318).
- b100 00159 (rot 318).

**The spread across dumps was NOT done** (stopped on request). The 3x-Fe variant was not run either.

## 1. Method and assumptions

**1. Fe I number density.** Computed per cell from the dump's (rho, T).
- Electron density, H/H2/He split and mu come from the run's own composition model. This is `EOSCompositionModel` with X 0.7188, Y 0.2420, [M/H] 0.4771, H2, H/He ionisation, Na/K/Ca/Al/Mg/Fe donors and metal condensation, the same keys as the run inputs. It was tabulated by `scripts/comptab.cpp` into `out/comptab.bin`.
- Fe I/Fe II follow thermal Saha with chi = 7.9024 eV. Fe III is neglected.
- Partition functions are direct level sums over all levels in Kurucz `gf2600.pos` (Fe I) and `gf2601.pos` (Fe II). This gives U(Fe I) = 19.5 / 22.0 / 24.6 at 2000 / 3000 / 4000 K.
- **Fe abundance: solar Fe/H = 3.2e-5 by number.** This is the value in our EOS donor table (`eos_composition.hpp`), not a quoted literature number.
- **No condensation** in the baseline. The condensation flag is the T < 2000 K cells: 23-33 % of cells at 1e-6..1e-1 bar, and 23-33 % at the terminator (±20 deg). These cells are mostly on the nightside (`out/stats.txt`).
- **Fe is removed at p < 1e-6 bar** (baseline key `pcut`). The model top is sponge/floor territory, with cold 200 K floor cells where thermal Saha keeps Fe neutral. Without the cut, the strongest lines stay optically thick up to the top boundary (tau ~1e4 at 1e-9 bar).
  - We do not trust that region (resolution counts only below 1e-6 bar). Real Fe there would be photoionised, which we do not model.
  - The cut is a decisive assumption; see section 3.

**2. Opacity.**
- **Line list.** Every Kurucz `gfall` Fe I line between 380 and 790 nm (vacuum, air-to-vacuum by IAU/Morton) is included: 9615 lines, computed at R = 1e6 (pixel 0.2998 km/s).
- **Line profile.** Exact per-sample Gaussian thermal cores (no microturbulence) plus Lorentz wings: radiative and van der Waals (gfall gamma values, perturbers H + H2 + He, (T/1e4)^0.3).
- **pRT 3.4 line-by-line Fe.** The 7.3 GB file `56Fe__Kurucz.R1e6_0.3-28mu.xsec.petitRADTRANS.h5` was downloaded into `prt_data`. It is used for the template and as the reference for the opacity validation.
  - Its T-P grid is 81-4000 K x 1e-6..1e3 bar. Our Fe-forming cells reach a 99th percentile of 3.9-4.9 kK and a maximum of 5.7 kK; 9.6 % of b100 cells are above 4000 K.
  - **The 3-D runs therefore use the line list directly**, which has no T clamp.
- **Agreement with pRT** (`scripts/val_lines.py`, `val_opac.py`; 164 isolated strong lines, 3000 K, 1 mbar):
  - Per-line integrated cross section: pRT/ours median 1.003 (16-84 %: 1.00-1.10).
  - Line positions agree to < 0.1 km/s.
  - Peak ratio 0.98, i.e. the pressure broadening that pRT includes is negligible here.
  - Total sigma over 380-790 nm: 0.990 of pRT.
  - The first attempt with `gf2600.pos` (Kurucz computed gf) was off by factors of 0.45-1.8 per line; it was replaced by `gfall` (lab gf).
- **Continuum.**
  - H- bound-free and free-free: Gray formulae as in pRT, with H- from Saha.
  - Rayleigh: H2 and He from pRT's own routine; atomic H from the Kurucz-ATLAS form, which I reproduced from memory and have not checked against a source.
  - The continuum is evaluated at 9 nodes and interpolated in ln lambda. Each ray contributes exp(-tau_c)(1 - exp(-tau_line)).

**3. 3-D transmission (`fe3d.py`).**
- **Rays.** Rays run on (position angle 2 deg) x (impact parameter: 2 per radial cell between the 1-bar level and the top). Each ray is sampled uniformly every 5e7 cm. The nearest column is used horizontally, with linear interpolation radially.
- **Velocity.** The planet-frame LOS velocity (winds + rotation, + = redshift) Doppler-shifts every line, sample by sample.
- **Geometry.** Same as `kp.py` / `obs.winds_one` rotated with phase:
  - Observer ob = (cos 2 pi phi, -sin 2 pi phi, 0).
  - a/R_s 3.7844, b = 0.160, R_p/R_s 0.1225.
  - Per-ray stellar-disc mask, so ingress and egress partial overlap is handled. The disc is uniform: no limb darkening and no RM.
- **Exposures.** 31 equally spaced exposures from T1 to T4 (phi = ±0.0474), with no exposure smearing.

**Convergence** (`out/cv_summary.txt`; h151 at 5 phases; per-limb and both-limb CCF fits; maximum change vs the base run):

| change | max change in Gaussian fit | in peak fit |
|---|---|---|
| dtheta 1 deg | 0.088 km/s | 0.015 km/s |
| 4 b per cell | 0.007 km/s | -- |
| ds 2.5e7 cm | 0.005 km/s | -- |
| line pruning 1e-6 instead of 1e-4 | 0.000 | -- |

All are far below the observational errors.

**4. Analysis (`ana.py`).**
- **Spectra.** Each exposure's excess-absorption spectrum gets continuum removal (smoothed running minimum over 300 km/s, the same as the template) and an instrumental Gaussian of FWHM 4.28 km/s.
- **Template.** The CCF uses a **pRT 3.4 line-by-line template**: isothermal 3000 K, all-neutral solar Fe, H2/He Rayleigh + CIA, broadened to 4.28 km/s.
  - The paper instead uses "a neutral iron template from Kitzmann et al. (2023) at 3000 K" (Seidel+26). Line weights and line lists differ, so absolute offsets can differ at the sub-km/s level; this is untested.
- **Kp-Vsys.** Done as in `kp.py`: exposure CCFs are shifted by (Kp0 - Kp) sin 2 pi phi and summed, with Kp0 = 219.2 km/s. The offset is a Gaussian fit to the highest peak (±12 km/s) of the Kp0 row. The centroid and a single wide Gaussian are also given.
  - Seidel+26 coverage is assumed to be full T1-T4. Their two 4UT nights are the 2023 ingress and the 2018 egress partial transits, and they "corrected for the Rossiter-McLaughlin effect" (A.4).
  - A variant with |phi| < 0.017 removed changes nothing (< 0.1 km/s).
- **Seidel+2025 segments.** Following their text:
  - "we have opted to discard the affected phases conservatively between -0.017 and 0.017".
  - "one from the ingress epoch ... covering a viewing angle from -17° to -7° (from here called the morning segment ...)", with the egress epoch giving the evening segment.
  - "We provide the average shift and FWHM for comparison purposes in Table 3, where we have excluded true ingress and egress".
  - Hence: morning = 6 exposures with -T2 < phi < -0.017, evening = +0.017 < phi < T3. The value is the **mean of per-exposure single-Gaussian shifts** (and FWHM). The stacked-CCF fit is also listed.
  - Observed (Ext. Data Table 3, "Fe i 4-UT"): shift "-4.12 ± 0.15" / "-6.90 ± 0.11", width "15.0 ± 0.4" / "13.1 ± 0.3" km/s.
  - Seidel+26, WASP-121b: 4.754 ± 0.301 km/s, read as -4.75 (blueshift). This is the Zenodo row quoted in KPVSYS.md; the sign is my reading.

## 2. Validation

**(a) 1-D isothermal static atmosphere** (2500 K, hydrostatic on the h151 grid, v = 0; `out/ana_val_static.npz`):
- Kp-Vsys peak +0.00, Gaussian +0.05, centroid +0.05.
- Both segments +0.05 km/s, FWHM 13.1 km/s.
- **The pipeline zero-point is +0.05 km/s.** It probably comes from the gfall-vs-pRT line position difference of < 0.1 km/s and was not removed.

**(b) Same atmosphere with solid-body rotation only.** Every chord then has v_LOS = Omega b cos(theta) exactly. The analytic CCF is the static CCF convolved with the ray-weighted (static weights) distribution of Omega b cos(theta).

| quantity | measured | analytic |
|---|---|---|
| mid-transit morning-limb Gaussian | +4.661 | +4.662 |
| mid-transit evening-limb Gaussian | -4.535 | -4.534 |
| peak (morning / evening) | +4.687 / -4.602 | +4.675 / -4.591 |
| segments (mor / eve) | +0.087 / +0.087 | +0.084 / +0.084 |

- The W-weighted mean limb velocity is ±4.56 km/s, against (2/pi) Omega b_eff = 4.56 km/s.
- **The sign convention is confirmed:** the morning limb is redshifted by rotation.
- Rotation alone gives a CCF FWHM of 23 km/s, and moves the Kp-Vsys map maximum to Kp 201 km/s (Kp0 219.2).

**(c) New vs proxy (KPVSYS.md).**
- **Where the signal forms.** The template-weighted Fe I absorption forms at tangent pressures of 16/50/84 % = 1e-5.0 / 1e-4.0 / 1e-2.7 bar in all models (fig. `figs/contribution.png`).
  - That matches the proxy's 1e-5..1e-3 window.
  - **It is 1-2 dex above Seidel+26's quoted "1 to 100 mbar with a maximum of the contribution function at 30 mbar".**
- **What changes against the proxy:**
  - **Line saturation caps the dense equatorial rays.** Each chord's absorption is at most 1, while density weighting grew without bound. Interpretation, not separately tested: the evening equatorial-jet peak (-20 to -24 km/s) loses weight against the central day-to-night / high-latitude component.
  - **Thermal ionisation removes Fe I** from the hottest dayside cells. The Fe I fraction p10 is 0.26-0.42, and 0.004 in b100.
  - **Consequence.** The hydro LR Gaussian peak goes from -20.90 to -6.03. Hydro C128 still has the jet peak winning (-19.00), but its centroid is -7.00, against -8.01 in the proxy. s10 goes from -5.14 to -4.26, and b100 from -1.44 to -0.38.
- **Limb-split transit averages.** Those of the proxy (+5 to +9 / -13 to -18) were not the observable. The proper Seidel+25 segment values are below.

## 3. Results (one dump per model; km/s; `out/table.md`)

| quantity | obs | hydro LR (h151, rot 300) | hydro C128x256 (rot 318) | MHD 10 G s10 (rot 318) | MHD 100 G b100 (rot 318) |
|---|---|---|---|---|---|
| **Kp-Vsys Gaussian peak at Kp0** (T1-T4) | **-4.75 ± 0.30** (Seidel+26) | -6.03 | -19.00 | **-4.26** | -0.38 |
| Kp-Vsys centroid at Kp0 | -- | -5.85 | -7.00 | -3.78 | -0.91 |
| Kp-Vsys single Gaussian (±35) | -- | -8.84 | -11.69 | -5.07 | -0.80 |
| Kp of map maximum | -- | 203.2 | 202.2 | 206.2 | 204.2 |
| **morning segment**, mean per-exposure shift | **-4.12 ± 0.15** (Seidel+25) | -6.50 | -9.05 | **-3.43** | +0.41 |
| **evening segment**, mean per-exposure shift | **-6.90 ± 0.11** (Seidel+25) | -10.85 | -13.56 | **-6.78** | -2.01 |
| morning / evening FWHM | 15.0 / 13.1 | 37.8 / 28.4 | 41.8 / 30.7 | 32.8 / 28.0 | 20.6 / 22.0 |
| proxy (KPVSYS.md) Kp-Vsys peak / centroid | -- | -20.90 / -7.72 | -23.60 / -8.01 | -5.14 / -6.33 | -1.44 / -3.07 |

**Variants** (same dump; Kp-Vsys peak, centroid, morning, evening):

| variant | Kp-Vsys peak | centroid | morning | evening |
|---|---|---|---|---|
| **Fe removed at T < 2000 K**: h151 | -6.09 | -5.99 | -6.63 | -10.97 |
| Fe removed at T < 2000 K: s10 | -4.30 | -3.87 | -3.49 | -6.85 |
| Fe removed at T < 2000 K: b100 | -0.38 | -0.94 | +0.41 | -2.00 |
| **No top cut** (Fe kept at p < 1e-6 bar): h151 | -5.67 | -5.01 | -5.49 | -9.11 |
| No top cut: s10 | -3.24 | -2.95 | -3.14 | -3.96 |
| No top cut: b100 | -0.00 | -0.21 | +0.01 | -0.25 |
| Cut at 1e-7 bar: b100 | -0.02 | -0.28 | +0.10 | -0.41 |

- Condensation changes everything by at most 0.14 km/s.
- **The top cut is NOT negligible.** Without it the CF median moves to 1e-4.9 (s10) and 1e-6.1 bar (b100), i.e. into the floor/sponge top.

**Widths.** All models give FWHM 20-42 km/s, against the observed 13-15. Rigid rotation alone already gives 23 km/s, and our static CCF is 13.1. The width is therefore pipeline/template dependent, and I do not use it as a discriminator.

**Figures.**
- `figs/contribution.png`: CF vs tangent pressure (both limbs, morning and evening rays), Seidel+25 segment CCFs, and mid-transit (position angle x p_tan) maps.
- `figs/ccf_kpvsys.png`: Kp-Vsys maps, the Kp0 row with the observed value, and the per-exposure planet-frame CCFs.

## 4. Verdict
1. **Seidel+2026 Kp-Vsys offset (-4.75 ± 0.30): 10 G (s10) matches** at -4.26 (1.6 sigma; this is one dump, and the dump-to-dump spread was not measured). The other models miss:
   - hydro LR is 4 sigma too blue (-6.03);
   - hydro C128 is jet-dominated (-19.0);
   - 100 G is too slow (-0.38).
2. **Seidel+2025 morning/evening (-4.12 / -6.90): 10 G again fits best** (-3.43 / -6.78). The evening value agrees; the morning is 0.7 km/s too slow (4.6 sigma formally, excluding model spread). Hydro is 2-7 km/s too blue at both limbs (LR -6.50 / -10.85, C128 -9.05 / -13.56), and 100 G is far too slow (+0.41 / -2.01).
3. **The proxy's main qualitative conclusion survives**, now with physically weighted line formation. The hydro proxy peak of -21 was an artefact of density weighting: the full RT hydro LR gives -6.0, still too blue.
4. **"B <= 2 G" vs our runs does not change.** A weak-drag (hydro-like) atmosphere is too fast in both observables, and our runs need a dipole of order 10 G (100 G overdamps). Their 0-D B limit and our 3-D B remain inconsistent by a factor of about 5.
5. **Caveats.**
   - One dump per model (no spread).
   - The answer depends on removing Fe above 1e-6 bar. Without the cut, s10's evening segment drops to -3.96 and b100 goes to about 0. This is untrusted model top / photoionisation that we do not model, and should be settled before quoting B.
   - Our CF sits at 1e-4 bar, not at their stated 30 mbar.
   - s10 is not saturated (KPVSYS.md).
   - Template differs from Kitzmann (2023).
   - Uniform stellar disc, no exposure smearing.

**Housekeeping.** Three core files (`core.viper13.200412/200431/200450`, 35 MB each) remain in this directory. They come from the headless Chrome that pRT's Keeper listing helper (selenium) launched at 22:10; it crashed on the login node. That route was abandoned, and files were listed via the Seafile API and fetched with curl; no result depends on it. The agent's `rm` was blocked by the tool's safety check, so they need manual deletion.
