---
name: bstar-w7-linear-saturation
description: B star prod_w7 at 10 turnovers (09-17): single tilted stripe mode (L/4 lanes, coherent tau 1-1000); KE1 grew +1.2/turnover through turnover 8, saturated at 9-10; developed convection expected 12-15
metadata:
  type: project
---

B star prod_w7 (bench/bstar_fecz/prod_w7, entropy seed k 2-16, cache-every-stage), 2026-09-17:
- history KE1 ln-ratio per turnover: 2 +0.94, 3 +1.07, 4 +1.16, 5 +1.21, 6 +1.24, 7 +1.23,
  8 +1.13, 9 +0.73, 10 +0.03 (KE1 2.6e34 erg, horizontal 1.05e34 at turn 10).
- dump 10 maps (analysis_0916/w7_dump10/bstar/b_maps.png, I looked): ONE coherent tilted
  stripe pattern, four downflow lanes, wavelength ~L/4 (~95 Mm), identical from tau 1 to
  tau 1000; rms v_z 0.14 km/s (tau 1) -> 1.39 km/s (tau 1000), max 6.1 km/s at tau ~2900.
  This is the end of the LINEAR phase with a single dominant mode, not developed convection.

**How to apply:** redo maps/cuts/iso3d at ~15 turnovers (dump 15) to see the nonlinear
breakup; do not call the 10-turnover state "developed". Products: analysis_0916/w7_dump10
(slices_w7.py + iso3d.py --vthr-pct 70 --vthr-tau 300, commands in NOTES.txt).
See [[he-w8-onset-slow]], [[fmode-wb-cache-culprit]].

f-MODE at 10.4 turnovers (analysis_0916/fmode_w7_t10/RESULTS.txt): the OSCILLATORY (2,2)
component (|f-f_f| <= 2 df) is flat, +0.006/turnover (2.29e-7 -> 2.39e-7), same as the
+0.004+-0.004 at 5.5 turnovers -> the cache fix holds through saturation. The naive band fit
(+0.42/turn) is ALIASED: the (2,2) bin family carrier fell to 9e-6 Hz (10 % in-band) because a
near-DC (4,2) convective pattern now carries 69 % of the non-mean surface power (= the tilted
stripe mode in the maps). Always split bands into f_f and slow parts from now on.
