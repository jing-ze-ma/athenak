# NOTE 2026-11-17 Orion -> viper: steady start with clean_closure = predict converges to a stale closure (affects the cell layout)

Found by Orion's WB check of the staggered study. It is a shared-code bug in the steady start (clean_init = steady); the patch is attached.

- **Symptom:** the reg wb Cartesian deck (wb/rhse_cart.athinput) names no clean_closure, so the default predict is used. With force_reference = clean
  and hydro wellbalance_static + reconst, max v/c_s reaches 3.8e-4 in the CELL layout when forced to predict (4.2e-4 staggered),
  against 5.7e-9 with extrap1.
- **Cause:**
  - The predict closure is cached per (time, dt), and the steady start keeps time and dt fixed within each phase. So the start computed only
    2 formal solutions in 7 passes and converged to the steady state of a STALE closure. It stopped at 1e-5 and reported "converged".
  - force_reference = clean then captured that state's force, and the run relaxed towards the scheme's real steady state with a constant force error.
  - The first step after the start reused the stale closure too.
- **Fix:** 8b6979103 (+2 lines, rad_m1_clean.cpp): the steady start rebuilds the predict closure in every pass and for the first step.
  - Cell with predict: 3.8e-4 -> 3.9e-8. Staggered: 4.2e-4 -> 5.7e-9.
  - Patch: handover_1117_steadypredict/0001-*.patch. It is on Orion's local base 788e3391c (push 14 line); please re-apply it on your head.
  - Reg harness keys-off bitwise. The cell layout's default closure (extrap1) is unaffected.
- **Also, for force_reference (staggered-only so far):** the reference was taken at Newton iterate 0 and marked set even if that solve failed.
  dt_retry does not checkpoint clean members (see NOTE 42eb96e99), so a bad reference survived. Check whether the cell layout's capture has
  the same exposure.
