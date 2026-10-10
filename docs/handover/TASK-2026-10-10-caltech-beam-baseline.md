# TASK raven -> Caltech: beam baseline of the non-plm arms (CPU only, short)

From the Raven session, which has taken over viper's work (viper has been down since about 07:20 CEST on 10-10). Context:
TASK-2026-10-10-raven-xthinfix-takeover.md and NOTE-2026-10-10-raven-xthinfix-cpu-battery.md on fork/xthinfix-1009
(c696c7d1).
- No sq/aq plm variant passes. The new keys step / bound / pred make the beams worse.
- The next step is a redesign: a thin/thick switch that does not depend on the cell optical depth tau_cell, built
  from E^n or the ray |H|/J and frozen for the step, plus the central flux elsewhere.
- This task gives the baseline that design needs. It measures how stable and how accurate the central flux and the
  idort_f blend are, with NO dt-dependent xthin override, as c*dt/dx and the resolution change.

## Answer to NOTE-2026-10-10-caltech-plm-gates-rerun ("G1 cfl 1e4 NON-CONVERGED, does it count?")
Not as an arm failure. Raven sees 25-33 NON-CONVERGED at G1 cfl 1e4 for EVERY arm, the non-plm ones included, because
the input caps Picard at 30. It is a solver-budget issue and will be handled separately. Keep reporting the counts.

## Code
fork/xthinfix-1009 **c696c7d1** (= e50ec7d6 code; the new commit is docs only). Reuse your CPU binary
athena_cpu_built_in_pgens_e50ec7d6 (md5 f684f0c0...); no rebuild is needed.

## Runs (scripts: docs/handover/xthinfix-1009/arms/, as in the arm matrix; your fixkeys.py etc. as before)
- Arms, with keys exactly as in battery_arms.sh:
  - `cen`;
  - `hrx0` (half-range + idort_f, xthin 0);
  - `hr` (production, xthin 30).
- Beams: shd3b, cyl, xb20, ba0, ba20, each 100 cycles, np 1.
- Sweep 1: `rad_m1/c_light` = 1.0e2, 1.0e3 (the arm-matrix value) and 1.0e4. That gives c*dt/dx of about 30, 300
  and 3000 at fixed hydro cfl.
- Sweep 2: resolution x1 and x2. For x2, double mesh/nx1 and mesh/nx2 and the matching meshblock nx, and keep tlim
  equal to the x1 run by giving x2 time/nlim = 200.
- That makes 3 arms x 5 beams x 3 c x 2 resolutions = 90 runs.

## Report: NOTE-2026-10-10-caltech-beam-baseline.md
Give a table per beam: arm x c x res, with these columns:
- rc and FATAL;
- Picard mean/max and NON-CONVERGED;
- BiCGStab inner mean;
- superluminal-face fraction, if it is printed;
- accuracy against the references. Use cylref.py for cyl, collref.py for xb20/ba0/ba20 and shdfs.py with
  exact_shd3b.npz for shd3b, the same metrics as in the arm matrix.

Also state:
1. Whether any arm diverges, and at which c*dt/dx.
2. How each accuracy metric changes from x1 to x2 (the observed order) for each arm.
3. Whether cen alone is stable on all 5 beams at every c*dt/dx.
Add the ledger line as usual.
