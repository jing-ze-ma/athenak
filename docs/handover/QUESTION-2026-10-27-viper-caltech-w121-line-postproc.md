# QUESTION viper -> caltech (USER): have you post-processed any WASP-121b run with real line RT (Fe I CCF / Na D, or R~1000 emission)? Reply first, do not start anything

The user remembers that we did real line post-processing for the WASP-121b limb velocities. viper finds only this:
- **`synth2_trans`** (09-30): Fe I CCF and Na D through the 3-D state.
  - FastChem 3 chemistry, Kurucz lines, Barklem & Collet partition functions, Voigt profiles, H- and Rayleigh continuum.
  - R 600k grid convolved to ESPRESSO R 70k, 15 transit epochs, a CCF with a weighted mask.
- **`synth2_emis`** (09-30): angle-dependent emission at R1000.
- Both were applied ONLY to the w121prod_0929 hydro 1x and 10x states at rot 300.
- Scripts: `docs/handover/scripts/w121/synth2_trans/` and `synth2_emis/`. Results: viper `/viper/ptmp2/jinma/w121prod_0929/synth2_{trans,emis}/RESULTS.md`.
- The later comparisons used only the old rho-ds wind average in fixed pressure windows: the Raven 1x/3x/10x arms at rot 270-300 and the
  MHD arms b100_e13 / s10_e15 at rot 314-320. That method put the Fe morning limb at about +5 km/s,
  against -0.7 with the real CCF on the same 1x state.

**Please answer:**
1. Has Caltech run synth2_trans, synth2_emis, or any other line-RT post-processor (pRT, gCMCRT, your own) on a WASP-121b state?
   - If yes: which run, rotation, field strength or metallicity, the result path, and the Fe/Na morning/evening RVs.
2. In particular: has anything run on an MHD state (Caltech or Raven arms)?
3. Do you hold WASP-121b MHD restarts or bins that viper does not have? (viper/Raven have w121_mhd_1003 b100_e13 and s10_e15 at rot 320.)

Nothing to start. A separate TASK may follow: a magnetized WASP-121b test in the accretor allocation (accretor paused at an rst, then resumed).
