# NOTE 2026-11-01 viper -> Caltech (USER): tracer T2 with BOTH tracers is not bitwise on 89213fcb (ck restart format bug, fixed in ac98652f2); defaults branch pushed

## 1. Restart bug that the tracers exposed (fixed on fork/tracers-defaults-1101 @ ac98652f2)
- **Cause:** the implicit-ck restart section snapshots u0 per variable and counted ALL variables, passive scalars included. The VEF snapshot has 6 slots.
  - With both tracers (7 variables), CkRstCollect silently dropped the whole VEF section. The restarted run then re-forms the VEF factors and says so in its log: "...the factors are re-formed... (this restart is not bitwise)". That gives a ~1e-3 dt offset after the restart.
  - When a tracer is switched on at a restart and the run stops before the first ck refresh, the final rst write dereferences an unstaged slab and segfaults (rc 139).
- **Fix:** the XS/VEF u0 slabs carry only the gas variables (the ck pass never reads scalars), plus a null guard.
  - Key off (no scalars): the slab list is unchanged.
  - Side effect: hydro runs with passive scalars (acctrc) no longer write the unused scalar slabs. Cross-version restart compatibility is being checked.
- **Effect on your tracers-test-1101 (89213fcb) tests:**
  - T1 (no restart): unaffected.
  - T2 continuation with BOTH tracers: NOT bitwise there (expected, from this bug).
  - Single-tracer continuation: fine.
  - A restart that switches a tracer on and stops at cycle 0 can segfault at the final write.
  - **Please run T2 (and any T1 restart) with ac98652f2 or later, or read a both-tracer T2 mismatch as this bug.**

## 2. Defaults (USER decision): fork/tracers-defaults-1101
- These apply on **fresh starts only**, through the existing recipe::NewDefaults. On a fresh start the resolved value is embedded in the input, so restarts keep it. A restart without the key keeps the old behaviour. An explicit key always wins.
  - `valfven_noinfall_fac = 4` when valfven_max_cgs > 0 (radial x1).
  - `floor_tracer = true` when valfven_max_cgs > 0.
  - `xe_tracer = true` when ohmic_resistivity or ambipolar = eos and the EOS has the electron-fraction table (otherwise false, with a log line).
- When Mdrain/Edrain do not fit in NHISTORY_VARIABLES, they are dropped with a NOTE instead of a FATAL.
- These go to fin51-r-fix only after your T0-T2 pass.
- The same mechanism will make `<mesh>/cs_seam_fc_nolimit` default on (branch csvertex-defaults-1101), after your A/B passes.

## 3. x_e tracer now advects electrons per H nucleus (d6e73a6ba on the defaults branch)
- The scalar is y = n_e/n_H (n_H proportional to rho) instead of z = n_e/n_tot. This conserves electron number per gram when mu changes, ~1.3 dayside -> ~2.3 nightside.
- Conversion at use: z = y n_H/n_tot for eta_O; rho_i = y n_H m_ion for ambipolar.
- Your 89213fcb xe arms use z. dt and robustness results stand; night F would be up to ~1.7x lower than with y.
- Restart files from tracers-test-1101 with x_e are NOT compatible with the defaults branch (test files only).

## 4. FYI from our ambipolar merge (e2cf02799, gates passed)
- On the s10 state, your rkg_dt_safety eigenvalue bound is 16.49 vs a measured eigenvalue of 17.20 (rkg_lmax_diag), i.e. 4 % below the true value.
- Safety 0.5 covers it. Do not raise rkg_dt_safety above ~0.9 without a margin.
