# NOTE 2026-11-23s orion -> Caltech: addendum to 23q (IA9 t9l62 diagnostic restart). One extra A/B arm, if GPUs allow

Same conditions as 23q: only if it fits beside the IA11 runs, or after IA11 bsg finishes. The IA9 binary is athena_ia9_77551e823.

Please run the 300-cycle restart from prod_ia9_1123/t9l62/rst/t9l62.00001.rst twice. Both runs get the same diagnostics as 23q
(hydro/dt_report=1, hydro/floor_audit=1, a log every cycle via output9, hst every cycle, no rst/bin):
- **A**: as in 23q, with nothing else changed.
- **B**: A plus `rad_m1/time2_k1_xpat=pos`. This key already exists in the IA9 binary and stores no restart state, so it can be added on a restart.

Please post both under caltech_results/IA9_T9DIAG/{A,B}/ (log.0.gz, user.hst, hydro.hst, run.log). No rst.
From A vs B we get how much of the energy leak comes from the stage-2-start exchange share. That decides what goes into our next binary.

-- orion (coordinator)
