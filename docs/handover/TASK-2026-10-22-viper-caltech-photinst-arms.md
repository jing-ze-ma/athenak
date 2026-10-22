# TASK 2026-10-22 viper -> caltech (USER): run the photosphere-instability MECHANISM arms (opacity / coupling candidates)

USER (10-22): Caltech runs the candidate arms from NOTE c5f8b9e8 (items 1-4); viper does NOT run them. Viper keeps: the RSG E
base onset run (N=32, fresh to ~60 d, f55be371), the time2_vet_col = rebuild arm, and the physical-vs-numerical criterion (item 5);
the base onset time / e-fold will be posted for comparison.

Arms (one key set each vs your base, fixed grid, equal sim time; metric = e-fold time of the photospheric T sd/mean and the
panel contrast; on whichever deck reproduces the onset fastest — HG at N=32 grows in ~1 d per your TASK):
1. Lagged opacity in the implicit solve: implicit_opac_lag off / step (vs pred); implicit_opac_lag_passes larger;
   implicit_opac_update and implicit_tsolve_opac (+ _mode) variants.
2. Face opacity across the cliff: implicit_face_weight (distance vs harmonic / upwind); implicit_marshak_face; implicit_thin_hcface.
3. Gas-radiation exchange: kappa_P = kappa_R (eth = 1); implicit_gas_newton off / implicit_gas_newton_switch; H (or He) recombination
   in the T solve (check the eint(T) jump across the ionisation front against the cell's dlnT).
4. Opacity-table interpolation: table log T spacing vs dlnT per cell at the cliff (RSG E: 0.12 per cell); a finer-resampled table
   (same physics) as an arm.
Your no-effect list (TASK 9f4e78db) is taken as done. Please NOTE the e-fold per arm on bsg-files-1009.
