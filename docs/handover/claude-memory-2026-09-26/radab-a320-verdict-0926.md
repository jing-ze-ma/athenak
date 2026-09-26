---
name: radab-a320-verdict-0926
description: a320 vs a128 radial A/B verdict 09-26 - spin up on a128 (14.8x cheaper), then remap to 256/320 for the science stretch; a128 not converged at 1e-4..1e-2 bar (45 K, kinks)
metadata:
  type: project
---
/viper/ptmp2/jinma/radab_0924/ana_0926 (cmp.txt, ana.txt). Old planet: start hyd4 rot 91, compared over 0.5-3 rot, p > 1e-6 bar.
- **Below 1e-2 bar:** a128 is within 1-13 K (<= 0.4 %) of a320, not growing. Winds, jet and deep retrograde flow agree within scatter.
- **1e-4..1e-2 bar:** a128 has a 45 K mean offset (~20x noise), night side +25 K, and 650-800 kinked night columns that a320 does not have. The offset sets up within 0.5 rot.
- **Cost:** a320 is 14.8x a128 per rotation (4.1x ms/cycle, 3.6x smaller dt).
- **Recommendation:** spin up on a128, then remap to 256 8-coef (or 320) for the science stretch. This fits [[two-phase-relax-then-accurate]]. 256 itself was not tested.
- These runs predate the cs seam fix. No collapses. a320 gained 8.6e-4 mass right after the remap, then stayed flat.
