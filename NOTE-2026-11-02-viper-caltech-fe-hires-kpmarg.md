# NOTE 2026-11-02 viper -> Caltech (USER): fe_hires Kp marginalisation

Follow-up to NOTE-2026-11-02-viper-caltech-fe-hires-tool.md. Changed:
`docs/handover/scripts/fe_hires/{ana.py, planet.py, README.md, seidel2026_ed_table1.md}`.

## 1. What the paper says (Seidel+2026, ESO PDF eso2606a.pdf; arXiv:2606.22455 v1 text identical)
- Sec. 2, p. 6: "Given the uncertainty on the orbital velocity, we marginalise over the orbital
  velocity as well as over the planetary wind velocity. The resulting detection maps for each target
  can be found in Extended Data Figure 2."
- Methods 5.1.6, p. 17: "The offset relative to the measured system velocity, which encodes the
  atmospheric movement, is determined from the stacked Kp-Vsys diagrams. To generate these diagrams,
  we scan over a range of possible projected orbital velocities Kp = v_orb sin i, avoiding
  assumptions about the true orbital velocity or inclination and preventing bias from their
  uncertainties. Additionally, the map is corrected for the system velocity as derived for each
  individual night prior to stacking in Kp-Vsys space based on S/N per night of observation."
- Fig. 1 caption, p. 7: "The uncertainty on the wind speed is the fit uncertainty of the observed
  cross-correlation peak, taking into account the photon noise of the data, as well as the
  propagated uncertainty of the system velocity measurement."
- ED Fig. 2 caption, p. 15: significance = "dividing the diagram by its standard deviation outside
  the velocity ranges affected by the star, planet, Rossiter-McLaughlin effect, and tellurics."
- **Not given anywhere** (main text, Methods, ED captions, arXiv): the Kp range or prior, the fit
  function for the planet peak, or how the Kp marginal is formed. vsys per night is a pymultinest
  Voigt fit to the out-of-transit stellar CCF (5.1.5, p. 16); the paper does not say median or mode.
  No Supplementary text on the Kp step (Suppl. Tables 1-2 = planet parameters, vsys).

## 2. Implemented (ana.py, new lines labelled `Kp-marg`; all old numbers unchanged)
- **(a) = the reading of the Sec. 2 sentence** (paper ambiguous, so flat prior Kp0 +- 50 km/s):
  - each Kp row: Gaussian+offset fit (+-12 km/s around the row maximum, the existing `peak`);
  - noise sigma = map std at |v| > 30 (ED Fig. 2 normalisation); row ln L = S^2/2, S = row peak
    above the |v| > 30 median in sigma units;
  - offset posterior = L-weighted mixture of N(mu(Kp), sigma_mu(Kp)) -> mean, median, 16/84 %.
- **(b)** rows averaged over the same Kp window, then the usual peak/wide/centroid fits.
- Planet JSON keys (optional): `kp_marg_range` [lo, hi] (null = Kp0 +- 50), `kp_marg_snr`
  (null = model map noise; else the map peak is scaled to this observational S/N).
- Done for `full`, `rmcut`, and the mor/eve segment maps. The segments use one side of the transit
  only, so Kp and v are degenerate: their Kp-marg numbers are printed but not meaningful.

## 3. WASP-121b s10 dump 00166 (rot 318), ana only on run_w121/out (log run_w121/ana_kpmarg.log)

| km/s | full | rmcut |
|---|---|---|
| Kp0 row Gaussian peak (old, unchanged) | -4.26 | -4.28 |
| centroid / single Gaussian (old) | -3.78 / -5.07 | -3.73 / -4.99 |
| map max v (Kp) (old) | -3.60 (206.2) | -3.40 (205.2) |
| **Kp-marg (a)** median [16, 84] | **-4.16 [-4.29, -4.04]** (<Kp> 205.7, S_max 13.4) | -4.11 [-4.26, -3.95] |
| (a) with kp_marg_snr = 10 / 5 (mean) | -4.17 / -4.24 | |
| **Kp-marg (b)** peak / wide | -4.82 / -5.05 | -5.03 / -4.96 |
| Seidel+2026 Fig. 1 (digitised) | 4.75 +- 0.30 blueshift | |

- Segments unchanged: mor -3.43, eve -6.78 (mean per-exposure).
- Kp-marg moves the s10 number by +0.10 (a) or -0.56 (b) km/s from the Kp0 row value.

## 4. How to run
```
python ana.py --planet planets/<p>.json --out OUT <tag>     # same command as before
```
ana only (about 10 s per dump); no fe3d rerun. Output: `Kp-marg ...` lines in the summary and
`kpm_<set>` dicts in `OUT/ana_<tag>.npz`.

## 5. Where the paper contradicts our notes (only the Kp step is implemented)
- **WASP-121b nights.** 5.1.1, p. 12: "For WASP-121 b, we only used the data taken in 4-UT mode".
  The tool NOTE sec. 3 "correction" (3 x 1-UT + 2 x 4-UT) and seidel2026_ed_table1.md were
  wrong; seidel2026.md ("E-4UT x2, 1-UT not used") was right. ed_table1 note now fixed.
- **WASP-76b nights.** 5.1.1: "we discarded the night of 18 Oct 2019 for WASP-76 b" -> 2 E + 3 MX
  (seidel2026.md right; the tool NOTE's "3 x E + 3 x MX" is the table, not what was used).
- **"Gaussian at Kp0, as in the paper"** (kpvsys/KPVSYS.md; the tool NOTE item-2 Kp row): the paper
  states no Kp0 row and no Gaussian for the planet peak. Gaussians appear only in the RM model
  (5.1.4); vsys uses a rotationally broadened Voigt (5.1.5).
- **"marginalised over Kp"** (tool NOTE item 2, seidel2026.md): the paper marginalises over Kp
  **and** over the wind velocity (Sec. 2); no range or prior is given.
- **Template.** Item-2 row matches (Kitzmann+2023, 3000 K, isothermal, hydrostatic + chemical
  equilibrium). Not in our notes: the mean out-of-transit **CCF** is divided after cross-correlation,
  and a 3rd-order colour correction is applied (5.1.3, p. 14).
- **Phases.** Not in our notes: exposures with S/N <~ 20 or deviant S/N are rejected per night
  (5.1.3). 5.1.6 warns that partial transits distort the mean velocity. No phase lists published.
- **RM.** Item-2 row matches (StarRotator correction, quadratic LD). The paper additionally excludes
  the RM velocity range from the significance std (ED Fig. 2); no RM phase mask. Our `rmcut` mask
  is Seidel+2025 only.
- **LSF.** Matches (2.14 / 4.28 / 3.53 km/s). With 4-UT-only data, WASP-121b is one LSF (4.28).
- **Fe pressure level.** Matches (1-100 mbar, peak 30 mbar, 3.2 / 5.2.2). The paper adds that the
  CF method "is unfortunately flawed" and the true levels "likely ... sit higher" (p. 22), at most
  two dex.
