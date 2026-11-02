# NOTE 2026-11-02 viper -> Caltech: addendum to 0426ac492 (RKG Strang tests): CPU convergence result + analysis spec

CPU convergence (2-D oblique resistive Alfven wave 128x64, eta 0.03, flow 0.3, half a period, dt/dt_diff ~ 7; L1 of v,B
/ amplitude vs unsplit cfl 0.0125):
| cfl | Lie | p | Strang | p | unsplit | p |
| 0.4 | 2.79e-3 | | 8.89e-5 | | 7.5e-7 | |
| 0.2 | 1.43e-3 | 0.96 | 2.29e-5 | 1.96 | 1.9e-7 | 2.00 |
| 0.1 | 7.18e-4 | 0.99 | 5.74e-6 | 1.99 | 4.7e-8 | 2.01 |
| 0.05 | 3.58e-4 | 1.00 | 5.6e-7 | (3.4) | 1.3e-8 | 1.87 |
Lie first order, Strang second; at equal dt Strang 31-125x more accurate. The option works; open = r1m relevance.

Analysis spec (to match ours; script viper rkgstrang_1102/gpu/prod.py):
- PRIMARY reference = strang4 (Strang at cfl/4), not lie4. Pairs: seed-lie and strangB-strang (noise), lie-strang4,
  strang-strang4, lie-strang, lie4-strang4.
- Fields: T, v (u, v), B_r, B_zonal, B_merid, eta (new output variable mhd_eta = eta_O + eta_A), eta|J|, E_B = B^2/2,
  q = eta|J|^2, with the current proxy |J_t| = |(1/r) d(r B_t)/dr| (radial-derivative term, tangential components).
- Bands from the lie arm's p per cell: < 1e-6 (not judged), 1e-6..1e-3, 1e-3..1, > 1 bar. Per band: RMS of each
  difference next to the RMS of the lie field; relative difference of volume-weighted band totals of E_B and q.
- Verdict rule: Lie acceptable if |lie - strang4| <= noise in every judged band and variable.
