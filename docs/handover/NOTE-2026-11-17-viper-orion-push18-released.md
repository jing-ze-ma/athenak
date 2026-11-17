# NOTE 2026-11-17 viper -> Orion (USER): push 18 released, clean-cc = b53e35565

USER go at 00:15 PST. fork/clean-cc fast-forwarded 5785ee51d -> b53e35565 (= fork/clean-cc-push18). Thanks for the CPU
regression (feba59618): A bitwise, restarts bitwise, timing flat.

Contents: 5eb7adf3f (clean FV rays get the vet_gd source of every shell), 84fe1b09f (dbg_rst_fp diagnostic), 88e5e0486
(CleanFvRows, bitwise), 1e4cfd388 (clean_fv_lag_frac OR rule, opt-in), b53e35565 (style). clean_yres default stays false.

GPU gates: gate A 30/30 bitwise on MI300A and A100; restarts bitwise incl. R136 (defaults and extrap); forced retry bitwise.
R136 moves under defaults (m1 5e-5 / hydro 4e-4 over 10 cycles), from 5eb7adf3f alone (the restart-bug fix).

Your notes taken for push 18b (in progress on clean-cc-physfix-viper): 40b0b48c5 (no-cut decided in the constructor),
ParameterDump current widths (9db95d9ea), FATAL for clean_rays = fv1 without vet_gd keys, the extrap tol on linear waves
(F1 35 % on rw: clean_fv_lag stays opt-in), and the summary ratio floor. The 12-rule check of clean-cc is running; NOTE to
follow.

-- viper
