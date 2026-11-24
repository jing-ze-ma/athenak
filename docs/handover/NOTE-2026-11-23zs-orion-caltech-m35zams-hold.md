# NOTE 2026-11-23zs orion -> Caltech: do NOT start an M35ZAMS 128² production run yet (ray-set symmetry defect found)

Our Raven M35ZAMS run (the same deck as your smoke) shows a spurious net +φ radiation flux in the optically thin layer just above the photosphere, present from t = 0:
- mean F_φ/(cE) +0.13-0.26 right after the steady start;
- the thin layer spins up coherently in +φ (v_φ ~ 2-3e7 cm/s by t 2.3e4 s), and that sets dt.

Suspects: the staggered GL 7x18 set, or the interpolated φ-wall maps on this narrow 3.67° wedge. We're testing both now.

- Your M35ZAMS smoke: finish as planned and post the timing / health report. It's still useful for s/cycle and memory. But don't continue it as a production start.
- Please also say whether the smoke's log shows v_φ growth in the hst (if you have a v_φ or KE_φ column) or bins.
- t9l62 and m35mms continue as they are. We're checking whether they're affected and will tell you.

-- orion (coordinator)
