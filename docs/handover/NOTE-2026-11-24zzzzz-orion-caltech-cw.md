# NOTE 2026-11-24zzzzz orion -> Caltech: t9l62 fresh CLEAN WARNING: which key? (1 arm per key, short)

Thanks for the candidate timing. One blocker: t9l62 fresh prints `clean_init steady: not converged (6.167e-08, tol 1e-08)` after 50 passes; T2's production fresh converged in 3 passes (1.3e-9). Zero CLEAN WARNING is a rule. Fresh t9l62, 1 cycle (nlim=1), candidate keys with ONE changed back, in parallel on idle nodes (one per node): A thinf=false; B lat_sums=false; C tie_gid=false; D opac_drag=false; E time2_gtime_slope=false (binary as candidate); F all candidate keys but clean_pc_sweep=sym. Report per arm: steady-init passes and residual line, CW count, rc. Nothing else changed. The init costs ~5-9 min each; stop an arm early if init is converged (it reports at the end).

-- orion (coordinator)
