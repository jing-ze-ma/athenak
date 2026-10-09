# NOTE viper -> Caltech: BSG FIRST on hpc-sm-02-11 (user 10-09)

User decision: the seeded BSG true repro gets the 2 H200 now; the He giant fresh run yields.
1. He giant fresh (chain 4288960-65, hegiant_1009/fresh897): stop it NOW (it is ~40 min in). If a restart file
   exists keep it; keep the run dir (do not delete). Cancel its pending links.
2. Start the seeded BSG (prod_seed, he_seed 1e-6 per NOTE-2026-10-09-viper-caltech-bsg-seed-GO) on those 2 H200 at once
   (smoke first if not yet done), 12 h links as before; energy-budget numbers in every link NOTE.
3. He giant afterwards: re-queue the fresh N897 chain on Caltech's shared gpu partition (6 h links, backfill) as a
   SECOND copy; viper re-submits it on Raven now. First start wins: push NOTE-2026-10-09-caltech-hegiant-fresh-2.md the
   moment it starts; before each link check for NOTE-2026-10-09-viper-caltech-hegiant-CANCEL.md (viper pushes it if Raven
   starts first). If it restarts from an existing rst, say so (t, cycle).
