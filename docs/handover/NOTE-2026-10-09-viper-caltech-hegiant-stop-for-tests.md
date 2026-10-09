# NOTE viper -> Caltech: stop the He giant fresh run; its 2 H200 go to radiation-scheme TESTS (user 10-09)

1. Stop the He giant fresh chain 4288960-65 now (keep the run dir and any rst; cancel pending links). Reason: its
   scheme (implicit_blend_xthin = 30) was found to lose convergence (viper order study); it will be restarted on the fixed
   scheme. Push a short NOTE with t / cycle reached.
2. Use the freed GPUs for the radiation-scheme test batteries viper sends next
   (TASK-2026-10-09-caltech-xthinfix-tests.md on this branch, coming shortly): short jobs, 1-2 H200 + CPU.
3. The seeded BSG (prod_seed, link 1 4294781) stays queued as it is; production keeps its place.
