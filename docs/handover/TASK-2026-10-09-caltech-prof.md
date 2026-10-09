# TASK viper -> Caltech: GPU cost breakdown of the production radiation scheme on H200 (short job; user 10-09)

Thanks for NOTE-2026-10-09-caltech-short-jobs. This is the profiling task DeltaAI/Delta have not run yet
(TASK-2026-10-09-deltaai-prof.md + scripts docs/handover/deltaai-prof-1009/ on this branch): first result wins — if
NOTE-2026-10-09-deltaai-prof.md or rad-beam-1008 NOTE-2026-10-09-delta-prof-batch2.md appears before your job starts,
skip it. H200 numbers are the most relevant (BSG/He giant production runs there).
- Binary: athena_gpu_he_star_m1_98835d99_nofma (md5 bc6f4e9e) — the AG Car/He giant production code (6c5d8fb2 is the
  same code plus the opt-in lowmem key; either is fine, say which).
- Layout: 1 node x 2 H200 if available, else 1 H200 (say which); one job, <= 1 h.
- Cases: (a) AG Car B, 60 cycles, fresh (bundle docs/handover/agcar-prod-1009/, SETUP.sh); (b) BSG reduced 4x1 mesh,
  30 cycles (bundle docs/handover/bsg-hrdet-1009/, as the DeltaAI mem tasks), dumps off.
- Arms per case: plain (also the smoke: stop on rc != 0 / FATAL), Kokkos simple-kernel-timer at 10 and at 60 (resp. 30)
  cycles (kokkos-tools master, profiling/simple-kernel-timer; the short run is subtracted), and rad_m1/implicit_timers = 3
  via an input copy (command-line keys must already exist in the input).
- Report NOTE-2026-10-09-caltech-prof.md: s/cycle; top 25 kernels with % of total; grouped: hydro (fluxes/update/c2p),
  rad_m1 implicit (assembly, Krylov vector ops, mg preconditioner, Picard update), vet_gd sweep (+twin), half-range pass,
  vet_col/lat, halo/MPI pack-unpack (+ time outside kernels), outputs; radiation vs hydro fraction; Picard mean/max,
  inner its, NON-CONVERGED. prof_group.py in deltaai-prof-1009/ does the grouping (label file labels_98835d99.tsv).
