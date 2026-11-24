# NOTE 2026-11-24zzw orion -> Caltech: question: which preconditioner made bsg much faster, before the new staggered scheme? And idort names for gate A

1. **USER question:** before the new staggered scheme (so roughly before 2026-11-16: the old prod forms / collocated clean scheme, the early bsg runs), bsg ran much faster with a preconditioner other than pure multigrid. Please search your bsg notes, logs, decks and results (docs/handover NOTEs, caltech_results, prod decks of Oct to mid-Nov, caltech-runs/ branches) and tell us:
   - the name and key (e.g. line / tridiagonal radial solve, block Jacobi, DSA / HYPRE, PCR, ADI, two-level, `implicit_precond` or similar), and where it lives (branch + sha);
   - the measurement: bsg s/cyc or Krylov its for it vs pure mg, the deck and the date;
   - whether it was on GPU;
   - why it was dropped or not carried into the staggered scheme, if you know.
   No new runs. Quote the source line / file for each number; say "not found" where you can't find it.
2. **Gate A fix (24zzq):** `--only idort` selected 0 tests because --only takes test names. Please rerun in the kept CPU build (wt_cpu, 29046c3f2): `--only idort_beams,idort_beams_mb,idort_shadow,idort_shadow_gl,idort_cylinder,idort_radshock_m2,idort_radshock_m5,idort_thermal_eq2`, and post the PASS / FAIL lines. Gate B "key on" for wb_sp / wb_cs fails at start-up: those decks need `rad_vet/vet_scatter=true` and `rad_vet/temp_unit_kelvin`. That's a bad request of mine; no rerun is needed.

-- orion (coordinator)
