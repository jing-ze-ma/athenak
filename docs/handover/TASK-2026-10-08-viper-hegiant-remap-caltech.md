# TASK for viper (from Caltech, user 10-08): N897 grid keys + remap recipe, and the fresh-start reference numbers

User plan on Caltech: run the N445 fresh start (files from 96b09994) on H200 to the remap point t = 9.004208e5 s (10.42 d),
remap onto the N897 grid with `docs/handover/scripts/he_remap_rst.py`, smoke, then continue to 30 d. Caltech smoke of the
fresh start (1 H200, 60 cycles) was clean: 0.955 s/cycle, Picard mean 4.95 max 10, NON-CONV 0, he_ic_balance 0.00369707 at
2.15098e11, last hst t 1202.8153622134355, mass 1.5435376857928201e32. Leave hegiant897 / DeltaAI copies alone; the user decides.

Please push to `docs/handover/hegiant-files-1007/` (fast-forward, never force) and reply as `NOTE-2026-10-0x-viper-hegiant-remap-caltech.md`:
1. The **reference smoke 12123795** table promised in NOTE-2026-10-08 (dt at cycles 0/1/10/20/59/60, last hst rows of hydro and
   user hst, Picard solves/mean/max/NON-CONV, NEWTON-FALLBACK count, he_ic_balance line) + `ref_smoke/` (run.log, hst).
2. The **N897 input**: `hegiant_n897_rst_embedded.athinput` (the input embedded in hegiant.00022.rst), with the 3 file paths as @BUNDLE@.
3. The **exact remap command** that made hegiant.00022.rst from the scout128 restart (all he_remap_rst.py arguments, the source rst
   name/time, which helper scripts were needed), and its checks (mass/energy conservation numbers you got).
4. The **restart keys** for the first N897 link (restart_keys.txt, output last_time handling, restart_refill_ghosts, opac_newton_slope_max)
   and anything changed between scout128 and N897 beyond the grid (e.g. the scaffold ramp, `implicit_opac_newton_slope_max = 3`).
5. Whether the scout128 run used any key changes mid-run before 10.42 d (run.cfg edits between links).
