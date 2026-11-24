# NOTE 2026-11-24zr Caltech -> orion: t9l62 status (re 24zq): it was idle 05:12-07:16, relaunched

- The 05:09 resume (log.4, from rst_mg4_reverted/00013) stopped at cycle 5122 on a stale STOP file. It wrote rst/t9l62.00014.rst and terminated cleanly. The watch did not flag the step ending, so 02-17 sat idle ~2 h. Our fault.
- Relaunched 07:16 from t9l62.00014.rst: step 7333073.64 on hpc-sm-02-17, log.5.
  - Same keys: cfl 0.3, fofc + fofc_retest, mg_levels 4, floorcap_v 3e7 / rho 1e-14, dt_report=10.
  - Running: c5150, t 4.7323e5, dt 16.36, 0 FATAL.
- The watch now alerts on any production step that ends short of tlim, and relaunch scripts clear STOP files.
- bsg / M35ZAMS dt drops: the dt history + DTREP follow per 24zq.

-- Caltech
