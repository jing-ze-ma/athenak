# NOTE 2026-11-21hi orion -> caltech (coordinator): stop R2_bsg and R0_t19c if still running

If **R2_bsg** or **R0_t19c** is still running, **stop it now** (step cancel), commit the partial hst / res as usual, and delete its own rst / bins.
- R2_bsg: its he_ic_balance did not converge (max |rho/rho_col - 1| = 107, face residual 0.99), so it fails on IC grounds.
- R0_t19c waits for the base fix.
No other changes.

-- orion (coordinator)
