# TASK 2026-10-22 caltech -> viper (USER): judge fixes by the SEAM/PANEL IMPRINT, not by the cliff disturbance (re 925609bd)

USER (10-22): "it's fine they have disturbance, as long as there is no significant imprint of the seam" -- please check on RSG E.
Acceptance metric for every arm = panel contrast: T_eff (or L at 1.2-1.4 R_ph) face-centre (|x|,|y|<0.3) / edge (max>0.8) - 1, and
corr(T_eff residual, cell dOmega); plus dt / retries. The cliff-shell sd is secondary.

## HG N=32 re-scored by panel contrast (L(1.2-1.4) c-e; base 1.9e-5 at 1.0 d, 1.06e-3 at 1.5 d, 4.6e-3 at 2.5 d)
| arm | contrast | at |
|---|---|---|
| nx1 1408 (2x radial) | -7.9e-6 | 2.5 d (600x below base) |
| nx1 704 photosphere-packed x2 (same cell count) | -5e-7 | 1.5 d |
| nx1 800 packed x2 | +1.6e-7 | 1.5 d |
| mlt_fsub_smooth 2 | +9.6e-6 (corr +0.02) | 1.5 d |
| implicit_face_weight equal | +5e-6 | 1.33 d |
| opac_lag_passes 4 (fresh; seed caveat) | -2e-6 | 1.5 d |
Restart arms from a growing state (lag/convergence keys) = base. Caltech now runs packed grids / sm2 / 2x nx1 fresh to 4 d.

## Asked of viper (RSG E, N=32)
1. Score your base and 40 d restart arms (fsub 4/8, rebuild) by the panel contrast above, vs time to 60 d. In particular: does the
   panel contrast stay small while the cliff disturbance grows (your 48 d: |fc/ed - 1| <= 1.6e-4 with e/rho sd 0.17)? When does it
   become significant (say > 1e-3), and does it correlate with dOmega?
2. If cheap: a photosphere-packed radial grid at the same nx1 (x2 finer over ~0.99-1.0 R_ph) and 2x nx1, fresh to ~60 d, same metric.
