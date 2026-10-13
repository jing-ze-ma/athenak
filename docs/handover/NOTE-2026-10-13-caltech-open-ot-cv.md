# NOTE 2026-10-13 from Caltech: open-ot CV (re TASK-2026-10-13-caltech-open-ot-cv.md, fixed bundle d57b1ed2)

- Binary NEW: 281dd686 + patch/0001 (git apply onto an archive of 281dd686; patch md5 f3baaf3b, From 37157508) -> athena_cpu_omp_281dd686_cs2vimp,
  md5 dc58375a9206e2c034fc5d57cfd8ac9d (plain CPU, OpenMP + MPI, Release).
- First run (job 4429039) used the OLD bundle: all 162 FATAL in 7 s (`implicit_cs2_vimp ... on command line not found`), as you found. Discarded.
- Rerun from the fixed bundle: debug job 4429421 (32 cores, 16 parallel, 28-min cutoff), rounds until done. Caltech changes only:
  m1_shock_ref -> our copies of caltech_od/shockref/ref_m*.txt; /usr/bin/time -> stand-in; cloud + pipe c220/c2200 left out (long; listed as not done).
- Counts and results follow here.

## INTERIM (10-12 23:35 PDT, after 3 debug rounds)
- **155 / 165 rc=0, 0 failures.**
- Not done: rs_m1.2/CV_n512, rs_m2/CV_n512, rs_m2/CV_n1024, rs_m3/CV_n512, rs_m3/CV_n1024 (cut by the 28-min limit while 16 x 8-thread runs
  shared 32 cores; alone they take 6-18 min) -> extra round 4432264 runs only these 5 (6 threads each).
- Not run on Caltech (long, > 28 min): cart/ddf/CV_d, cart/ddiff/CV_d (cut twice), cart/cloud/CV_r12g2, cart/pipe/CV_c220, cart/pipe/CV_c2200.
- Results tarball after the rs round.

## DONE (10-13 00:05 PDT): 160 / 165 rc=0, 0 failures
- rs round 4432264 (6 threads): rs_m1.2/CV_n512 306 s, rs_m2/CV_n512 333 s, rs_m2/CV_n1024 886 s, rs_m3/CV_n512 411 s, rs_m3/CV_n1024 905 s, all rc=0.
- Not done on Caltech (long, > 28 min in debug): cart/ddf/CV_d, cart/ddiff/CV_d, cart/cloud/CV_r12g2, cart/pipe/CV_c220, cart/pipe/CV_c2200.
- Results under docs/handover/caltech_cv/ (the 160 rc=0 dirs, paths runs/<list dir>/..., no rst):
  - cv_runs_tab.tgz (9.9 MB): tab/ run.log rc.txt time.txt;
  - cv_runs_bin.tgz.part00 + .part01 (70 MB total, split -b 44M): bin/. Rebuild: `cat cv_runs_bin.tgz.part0* > cv_runs_bin.tgz; md5sum -c cv_runs_bin.tgz.md5` (051ecf937e0dfb60c601fdee9cbf0b26).
