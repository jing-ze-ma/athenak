# NOTE 2026-11-01 viper -> Caltech (USER): WASP-121b Fe wind from our post-processed 3x resistive MHD run vs Seidel 2026, and the plan this feeds

FYI. Nothing to run. This is what our model already gives for the WASP-121b point of the Seidel et al. 2026 trend.

## Result (fe_hires, 10-07; viper /viper/ptmp2/jinma/w121_obs_1003/mhd_final/fe_hires/FE_HIRES.md)
- **Method:**
  - full high-resolution Fe I transmission: Kurucz gfall, 9615 lines at R = 1e6, validated against petitRADTRANS 3.4 to 1 %;
  - the run's own composition model;
  - CCF with a pRT template, a Kp-Vsys map, and the Seidel+25 morning/evening segments.
- **Runs:** WASP-121b 3x metals, C32 x nx1 112, one dump each at rot 318:
  - s10 = s10_e15: 10 G polar dipole at R_p, ohmic_resistivity eos, max_eta 1e15, RKG2 STS;
  - b100: 100 G;
  - hydro C32 (h151) and C128 (hr).

| run | Kp-Vsys peak [km/s] | morning / evening [km/s] |
|---|---|---|
| **observed** (Seidel+25; WASP-121b point of Seidel+26) | **-4.75 +- 0.30** | **-4.12 / -6.90** |
| s10, 10 G resistive MHD | **-4.26** | **-3.43 / -6.78** |
| hydro C32 | -6.0 | |
| hydro C128 | -19 | |
| b100, 100 G MHD | -0.4 | |

- **10 G matches** the net shift (1.6 sigma) and both limbs. Hydro is too fast; 100 G nearly stops the wind.
- **Tension with the paper:** Seidel+26 infer B <~ 2 G, and our match needs a 10 G dipole. The field geometry, the eta cap (max_eta 1e15) and the ionization treatment differ from their drag model, so the B values are not directly comparable.
- The 10 G run was not saturated (field still growing at rot 318-320).

## Caveats
1. **Pre-fix solver.** s10_e15 covered rot 300-320 with the OLD LHLLD, which was anti-dissipative (see NOTE-2026-11-01-viper-caltech-w121-status, item 2). The damage we found is deep, far below the Fe-forming layer, but the run has not been repeated on fin51-r-fix.
2. **One dump per run**, so there is no dump spread.
3. **Top-cut dependence.**
   - Fe is cut at p < 1e-6 bar. Without the cut, the s10 evening limb is -3.96 and b100 is about 0.
   - The MHD runs put 27-50 % of the contribution within 1 dex of the top cell, so the top boundary matters.
4. **Line widths are 1.5-3x too broad** in every run (incl. 100 G). Unresolved.
5. **Formation height:** Fe forms at ~0.1 mbar in our model, vs the 1-100 mbar the paper assumes.
6. **A cruder pipeline gives the same ranking but redder morning limbs.** synth2_trans (mhd_final/synth2_trans/RESULTS.md) puts the morning limb 4-8 km/s redder than fe_hires; no arm gives a blueshifted morning limb there.

## The plan this feeds (written 10-08, on hold until now)
- **Stage 1, WASP-121b B calibration.** Extend 10 G (unsaturated); new 3 G and 30 G arms; a 10 G arm with max_eta 1e13 to isolate the eta-cap effect.
- **Stage 2, the trend.** The five-planet sample (NOTE-2026-11-01-viper-caltech-dhj-sample-seidel2026): hydro, then MHD at the best B from stage 1 (and a second B if the trend is steep). Same geometry and eta cap. Per-planet metallicity, with a bracket where it is poorly constrained.
- **Analysis per run with the same fe_hires pipeline:**
  - Kp-Vsys offset, morning/evening, FWHM, formation pressure;
  - wind vs Teq plotted over the Seidel+26 points.
- **Before stage 2:** resolve the too-broad Fe widths and the top-cut dependence. Rerun stage 1 on the fixed solver (fin51-r-fix >= 182477553).
- **Correction to seidel2026.md:** our 10-08 notes have the paper as arXiv:2606.22455, submitted 2026-06-23; the package note says no arXiv version was found. This will be fixed in the next push.
