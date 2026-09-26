---
name: cs-seam-flux-positivity-0926
description: ROOT CAUSE of the WASP-121b crashes - cs seam flux averaging is not positivity-preserving; fix <mesh>/cs_seam_flux = positive|upwind (merged 2d07119a, default average)
metadata:
  type: project
---
Found 09-26 (/viper/ptmp2/jinma/seam_0926; docs/dev/cs_seam_defect_0926.md).
- **Mechanism.** src/bvals/flux_seam_cc.cpp (~l.473/475) averages the two panels' seam fluxes. Each panel's flux comes from along-seam resampled ghosts (bvals_cc.cpp ~l.441). An emptied seam cell looks denser to the other panel, which drains it (38x its mass per step at the first collapse). The legacy dfloor keeps the momentum, so the cell runs away and dt collapses.
- **Scope.** All 4 sponge arms died this way (base at rot 27). The ck cadence/tolerance is cleared: every=1 at tol 1e-9 still crashes.
- **Healthy base shows no seam anomaly** except mild extra roughness within 3 cells of night-side cube vertices.
- **Fix.** <mesh>/cs_seam_flux = upwind | positive (default average, bitwise). Conservative. Crash reruns are clean; rigid rotation is within 2 % of average.
- **Recommendation.** positive for dhj, after a multi-rotation check from a base restart. A new key needs an -i overlay on restart.
**How to apply:** every future cs dhj run sets cs_seam_flux = positive once checked. The floor switches (dfloor_keep_velocity/temperature, vceil) are a second net. Related: [[sparc-sponge-campaign-0925]], [[index-cubed-sphere]].
- **Audit 09-26 (seamaudit_0926):** MHD uses the same cc seam flux, so cs_seam_flux covers it (use positive for MHD too). EMF/fc path: no bug (div B ~4e-5 at seams/vertex). OPEN: a vertex defect. An isothermal hole on a cube vertex goes NaN in 2 cycles under all options; mid-seam is clean. Cavity reproducer: cs_test cav_*, branch cs-seam-mhd-0926.
- **Multi-rotation check DONE 09-26** (rot 20->32): positive and upwind both stable (0 warnings; base crashed at 27.04). Mean T within 0.6 K at 1e-2..10 bar; deep eint rms vs base grows to 6.4e-3 as scatter, no bias; the 1e-6..1e-4 bar layer is chaotic. VERDICT: use positive in all cs dhj runs.
