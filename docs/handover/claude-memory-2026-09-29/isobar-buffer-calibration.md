---
name: isobar-buffer-calibration
description: How to size x1max from the 1e-6 bar isobar, and the two validations any such measurement must pass
metadata:
  type: feedback
---

The recipe that produced the [[ck-grav-prod-run]] grid, worth reusing whenever a dhj domain
is resized.

**Measure on the ISOBAR** ([[dhj-isobar-vs-shell]]), over terminator columns
lon = +-90 +-20, |lat| < 30, **log-interpolating between cells** -- snapping to cell centres
quantizes r99 to dr = 5.558e7 and hides whether it has converged.

**Express the buffer in SCALE HEIGHTS, and calibrate the target against the accepted run,
not a remembered number.** ck_limb was accepted as meeting the science goal with 1.76 H
above its worst 1 % of terminator columns; the "~3 scale heights" in [[ck-limb-run]] refers
to the median-ish case and taking it literally over-builds the grid by 2.5x. Cells and
scale heights disagree violently: under point-mass gravity H at the isobar is 22 cells, so
40 cells of headroom is only 0.85 H.

**Two validations the inversion must pass before any of it is quoted**
([[eos-inversion-nan-trap]]):
1. `rho k T / (p m_H)` on the p = 1 barye surface must give **mu ~ 1.25** (neutral H+He),
   tight. The broken inversion gives ~1.95.
2. hydrostatic mu from the measured dp/dr with the analytic g(r) should agree -- it does to
   2 % in ck_limb. Where it does not (0.76 under point-mass gravity, i.e. H measured 1.31x
   hydrostatic) that is **wind support, not an inversion error**: the same 1.3 is already
   recorded for the ck_limb dayside top.

**Why: sizing decisions cost multi-day campaigns, and the two traps here (measure the
evolved state not the IC; validate the inversion) have each already cost a retraction.**

**How to apply:** oversize deliberately so the isobar never leaves the domain -- then the
answer is a measurement. Check the SIGN of the residual drift: if it is inward and H has
converged, the buffer can only grow and the grid can be chosen before the run finishes.
Scripts: `sizetrend2.py` / `isocheck.py` in the 2026-08-26 session scratchpad.
