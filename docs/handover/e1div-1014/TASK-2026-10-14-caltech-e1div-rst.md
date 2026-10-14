# TASK Caltech (DR, e1div-1014 part 2): two short E1 restart arms from rst 00022 -- diagnosis only, NOT production

Finding so far (from your profiles): the E1 envelope undergoes a growing, near-spherical radial pulsation (atmosphere inflates to ~1.2 R_ph,
falls back at -35 km/s near 1e7 s, rebounds +33 km/s at 1.05e7, falls again at -60..-80 km/s by 1.24e7). Picard non-convergence starts at
1.12 R_ph at cycle 15198 (t ~1.105e7) and follows the infall inward to 1.03 R_ph, where it diverges. dt stays 727.27 s (hydro CFL, cfl_number 0.3).
Question for these arms: does a halved time step carry the run through the infall (time-step problem), or does it fail the same way?

## Setup (do NOT touch e1_1013/prod, its rst or its run.cfg; same binary dad00933 ff7b7f11, same launch e1_1013/e1launch.sh, 1 node x 2 H200)
```
P=<abs>/e1_1013/prod; W=<abs>/e1div_rst; mkdir -p $W/R0/rst $W/R1/rst
ln $P/rst/agcar3d.00022.rst $W/R0/rst/ ; ln $P/rst/agcar3d.00022.rst $W/R1/rst/     # hard links (no copy); read-only use
K="output1/last_time=1.1e7 output2/last_time=1.1e7 output3/last_time=1.1e7 output4/last_time=1.1e7 output5/last_time=1.1e7 output5/dt=1e30"
# R0: baseline, must reproduce link 1 (expect DIVERGED near t 1.2425e7, cycle 17085)
cd $W/R0; e1launch.sh <BIN> -r $W/R0/rst/agcar3d.00022.rst $K time/tlim=1.26e7 output2/dt=2.5e4 output3/dt=2.5e4 -d $W/R0 -t 05:40:00 > run.log 2>&1
# R1: one change, time/cfl_number = 0.15 (dt -> ~364 s)
cd $W/R1; e1launch.sh <BIN> -r $W/R1/rst/agcar3d.00022.rst $K time/cfl_number=0.15 time/tlim=1.30e7 output2/dt=2.5e4 output3/dt=2.5e4 -d $W/R1 -t 05:40:00 > run.log 2>&1
```
- One 6 h job per arm (or both back to back in one 12 h job). output5/dt=1e30 = no new rst (disk). If the last_time keys are rejected, drop them
  (the bins then simply start at the next multiple of 2.5e4).
- If R0 is not identical to link 1 up to its divergence, say so (cmp of the first hst rows after 1.1e7 vs prod/agcar3d.user.hst).

## Report (push to docs/handover/e1div-1014/caltech_out_rst/, NOTE-2026-10-14-caltech-e1div-rst.md)
- For each arm: rc, end t/cycle, FATAL/DIVERGED line, counts of NOT ADMISSIBLE and "NON-CONVERGED after", the NC-local lines (all), Picard mean/max,
  dt at the end, s/cycle; the user.hst + hydro.hst.
- Run e1div_prof.py on each arm (`python3 e1div_prof.py $W/R0 $W/R0_out 6`, same for R1) and push the prof_* (skip the slices unless < 40 MB).
