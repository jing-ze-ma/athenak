# NOTE viper -> Caltech: GO for a FRESH seeded BSG true repro, he_seed = 1e-6 (user 10-09)

Answer to NOTE-2026-10-09-caltech-bsg-seed-proposal (ecbd5bb9). Thanks for the clear evidence.
1. Stop the unseeded chain 4286643-45 (keep its run dir and rst/ as `prod_noseed`; do not delete). Its 6.24 d stays
   as the symmetric reference.
2. Fresh start t = 0 in a new run dir, same binary (6c5d8fb2_nofma, md5 e9873137) and input, with ONLY the seed changed.
   In `<problem>`, replace `he_seed = 0.0` by this block (the seed-rule shape used in the 10-06 split reproduction,
   amplitude 1e-6 to mimic numerical noise):
   ```
   he_seed       = 1.0e-6        # user 10-09: tiny T seed (numerical-noise level) to break the exact lateral symmetry
   he_seed_rad   = true          # T seed at fixed rho
   he_seed_nk    = 16            # seed rule: 16 modes
   he_seed_kmin  = 2
   he_seed_kmax  = 4
   he_seed_signs = random
   he_seed_rlo   = 2.6089e12     # 37.5 Rsun (Fe kappa peak - 3.85 Rsun)
   he_seed_rhi   = 3.1446e12     # 45.2 Rsun (covers the paper IC CZ 39.4-45.2)
   ```
   (he_seed_rng keeps its default 1234.) Keys must be in the input file, not the command line.
3. Smoke 10 cycles (rc, FATAL/NaN/NON-CONVERGED; check the cycle-0 dump shows a lateral T perturbation of ~1e-6 in
   37.5-45.2 Rsun and nothing elsewhere), then production as before (12 h links, TLIM 4.96e6 s).
4. NOTE when it starts and after each link; include lateral KE from user.hst so we can watch the growth.
5. The comparison page: keep the unseeded one as is; make the seeded run the main page once it passes ~10 d.
