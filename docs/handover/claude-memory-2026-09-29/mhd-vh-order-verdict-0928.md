---
name: mhd-vh-order-verdict-0928
description: 09-28 MHD He box v_h time-order loss (1.48 at "77 G" = 77 code = 273 G) is an additive ~0.1 cm/s (1e-8 c_s) error, flow/field/solver-independent; negligible; no code change
metadata:
  type: project
---
/viper/ptmp2/jinma/mhdvh_0928/RESULTS.md. Excluded: rst conversion/caches, lhlld (hlle same), M1 coupling (coupling
off still degrades), wall guard, GS05 corner emf in the box, FOFC, PLM, tolerances. Order drops with B (1.81/1.66/1.48/
1.15 at B 0/7.7/77/243.5 code) because v_h there is Mach 1e-5; with a Mach-0.1 shear the absolute error only grows
~2-3x (0.04 -> 0.11 cm/s) while v grows 1000x -> additive ~0.1 cm/s, also at B=0. Negligible for real flows.
Side findings: (1) GS05 corner-emf upwind switch drops fixed-grid TIME order on 2-D oblique Alfven/slow waves
(0.40/1.48 -> 2.0 with averaging); (2) lhlld low-Mach correction erratic (2e-3 of amplitude) at Mach 1e-3, beta 5000
(phi ~0.03); hlle clean; (3) "77 G" label wrong: code units are Heaviside-Lorentz, 77 code = 273 G.
Revisit (1) on a production-like MHD case at the WASP-121b MHD phase ([[rmdx-analysis-0928]]).
