# R136a1 open top (r136open_1106) + E_tot-drop timeline

Unchanged copies from `$J/r136open_1106/`: REPORT.md, r136open.png, analysis scripts an.py, wbforce.py (predicted WB gravity
source above wb_rmax vs measured dv/dt), force.py (dv/dt force split per radius), tcmp.py, plot.py; `r136noise_1106_an1.py`
(= `$J/r136noise_1106/an1.py`, the loader both wbforce.py and force.py exec); `arms/<arm>.arm.info` for B, BL, C, M, ML, T0-T6
(keys of each test arm).

Gravity-source analysis behind `problem/he_wb_src_linear_rmin`: the numbers are in REPORT.md ("Force imbalance above 40 Rsun");
the wbforce.py / force.py stdout was not saved to files on Caltech (only printed in the agent session), so no table beyond REPORT.md exists.
`he_wb_src_linear_rmin.patch` = `diff -u` of the test source `$J/r136open_1106/build/exp/src/pgen/he_star_m1.cpp` against
`git show dc2dc20c:src/pgen/he_star_m1.cpp` (binary `builds_rstfixx/bin/athena_gpu_he_star_m1_dc2dc20c_srclin_fmadoff`). Not committed to any branch.

`r136a1_runB_Etot_drop_timeline.txt`: grep/sed/awk excerpts of `$J/r136a1_1104/runB/prod/run.log` (+ r136.user.hst rows) with
every command printed above its block: all dt COLLAPSE lines (+3 lines of cell info), posfix rounds >= 1e4 cells, counts of
NOT ADMISSIBLE / BE redo / NON-CONVERGED, the 6.641 d window (cycle 131896) and the 6.997 d window (cycles 137263-138002 incl. 137913),
the counters_rollback line (end of the segment at 11-07 08:18: BE redo 157, helper solves 858, E-floor energy dropped 0 / inf / inf),
and hst E_rad / e_gas / Etot rows around both times.
