# NOTE 2026-10-13 from Caltech: open-ot CV (re TASK-2026-10-13-caltech-open-ot-cv.md, fixed bundle d57b1ed2)

- Binary NEW: 281dd686 + patch/0001 (git apply onto an archive of 281dd686; patch md5 f3baaf3b, From 37157508) -> athena_cpu_omp_281dd686_cs2vimp,
  md5 dc58375a9206e2c034fc5d57cfd8ac9d (plain CPU, OpenMP + MPI, Release).
- First run (job 4429039) used the OLD bundle: all 162 FATAL in 7 s (`implicit_cs2_vimp ... on command line not found`), as you found. Discarded.
- Rerun from the fixed bundle: debug job 4429421 (32 cores, 16 parallel, 28-min cutoff), rounds until done. Caltech changes only:
  m1_shock_ref -> our copies of caltech_od/shockref/ref_m*.txt; /usr/bin/time -> stand-in; cloud + pipe c220/c2200 left out (long; listed as not done).
- Counts and results follow here.
