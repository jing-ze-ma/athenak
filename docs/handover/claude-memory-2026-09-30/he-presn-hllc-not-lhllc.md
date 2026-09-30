---
name: he-presn-hllc-not-lhllc
description: User 09-30 - He presupernova runs should use plain hllc (hydro) / hlld (MHD), not the low-Mach lhllc/lhlld: its convection is near-sonic (MLT Mach ~0.4-0.5)
metadata:
  type: feedback
---
User 09-30: "for he psn we probably just want hllc/hlld. the convection is near sonic, no?" Yes: 4 Msun presn He star
FeCZ v_MLT 1.45e7 cm/s, Mach ~0.42 (sizing: presn He stars Mach 0.35-0.55); 3-D wedge plumes v_r ~90 km/s.
**Why:** lhllc/lhlld only help at low Mach (they scale the velocity-jump dissipation by phi ~ Mach); at Mach ~0.4 the
benefit is small and the low-Mach fix brings its pressure-velocity-decoupling side effects ([[lhllc-radial-oddeven-0930]]).
**How to apply:** He presn inputs (he_presn_m1_wedge.athinput on he-presn-m1, future He-star global runs): rsolver hllc
(hlld for MHD). Keep lhllc for dhj (low Mach deep, with the dhj-only radial phi fix) and ask before changing the He box
(FeCZ box convection is much less efficient/lower Mach; the lhllc x1 check there is running). Related: [[dhj-rsolver-lhllc-lhlld]].
