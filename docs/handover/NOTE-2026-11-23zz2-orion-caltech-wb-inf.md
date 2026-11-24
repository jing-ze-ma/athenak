# NOTE 2026-11-23zz2 orion -> Caltech: please push arm_WA / arm_WB files, and check WB's cycle 6469

Your 23zz lists caltech_results/IA121_T9L62/arm_{WA,WB}/ (grep.txt.gz, log.0.gz, hst, skip_de_per_cycle.txt), but they aren't on bsg-files-1009 yet. Please push them, including the files so far for still-running arms.

The WB "implicit_pos_floor created inf" at cycle 6469 t2s 1 needs one check (`grep -n "cycle 6469" log.0` in arm_WB):
- If it's followed by "stage solve failed, the step is redone with backward Euler" and "... of the failed attempt is discarded with it", the inf was discarded. WB is valid; keep running.
- If neither line follows, the inf was kept: stop WB (test stop rule) and post the CLEAN cycle / WBFDIAG lines for cycles 6465-6470.

Post WA's and WB's status line (cycle, t, dt, CW, fallbacks, DT_RETRY) too.

-- orion (coordinator)
