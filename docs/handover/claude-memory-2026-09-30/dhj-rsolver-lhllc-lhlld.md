---
name: dhj-rsolver-lhllc-lhlld
description: user 09-29: all dhj runs use rsolver lhllc (hydro) and lhlld (MHD) from now on (was hllc / hlld); tell Caltech and DeltaAI
metadata:
  type: feedback
---
dhj (deep_hot_jupiter_rt, WASP-121b etc.): <hydro>/rsolver = lhllc and <mhd>/rsolver = lhlld from now on.
**Why:** user 09-29 after the hlld-vs-lhlld dissipation test (mhdvh_0928: plain HLL-type solvers damp low-Mach flow
3-60x more) and the finding that the WASP-121b production input used plain hllc for low-Mach deep flows.
**How to apply:** every new dhj input / worker brief; the note NOTE-2026-09-29-dhj-rsolver.md tells Caltech and
DeltaAI. The Caltech WASP-121b fresh start (NOTE-2026-09-28-caltech-w121prod.md, started 09-28/29 with hllc) needs
the user's decision to restart with lhllc. WASP-121b production runs on CALTECH AND (user 09-29) is also queued on viper: /viper/ptmp2/jinma/w121prod_0929 (1x 12018343, 10x 12018345, + afterany links; binary cda4da33).
