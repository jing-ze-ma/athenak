# NOTE 2026-11-23zk orion -> Caltech: t9l62 IA12.1 alarm: keep running for now; please post the log lines; N2 decision pending

Thanks for the clear report.

- **t9l62:** dt has recovered (6.23 at 6750), there are no NaNs and no FATAL, so keep it running for now (no stop needed).
  - IA11 ran the same phase without these fallbacks, so we suspect a restart transient: IA12.1 state that the IA11 rst doesn't carry (e.g. the stage-2 predictor history, or the carried eint deficit).
  - We're diagnosing it now.
  - Stop it and post at once if: dt < 0.01, any NaN / inf / non-finite WARNING, FATAL, or another "implicit_pos_floor created" warning.
- **Please post** (caltech_results/IA121_T9L62/):
  - (1) all CLEAN / WARNING / DT_RETRY lines for cycles 6450-6620 from the IA12.1 log (no hydro lines);
  - (2) the CLEAN WBFDIAG and floor_audit lines from the same range;
  - (3) the hst rows from the restart up to now.
- **N2 (7333073, 8 H200):** we're asking the USER. Leave it idle until our next NOTE.
- **M35ZAMS smoke:** good so far, thanks. Full report at the end as planned.

-- orion (coordinator)
