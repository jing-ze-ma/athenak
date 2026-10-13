# NOTE 2026-10-13 from Caltech: od-vetpath (re TASK-2026-10-13-caltech-od-vetpath.md)

## QUEUED (10-12 20:32 PDT)
- Binary: open-od-1013 @ 87bdbd1f plain CPU (OpenMP + MPI, Release, default pgens) -> athena_cpu_omp_87bdbd1f, md5 71fd84a69735c592b4969130774f4f59.
- Debug QoS job 4424769 (1 node x 32 cores, 28-min cutoff, 8 runs x 4 threads); repeated until the short runs are done (rc=0 dirs skipped).
- **Two bugs in caltech_od/run_cart.sh (fixed in our copy run_cart_caltech.sh, behaviour otherwise unchanged):**
  1. list.txt entries already start with `runs/`, but one() does `cd $HERE/runs/$1` -> the cd fails for every run (and the run would start in the wrong dir).
     Fix: `cd $HERE/$1 || return`.
  2. one() starts with `rm -rf bin tab rc.txt`, so a rerun never skips finished runs (needed for repeated 30-min debug jobs). Fix: skip dirs whose rc.txt is 0.
  (+ /usr/bin/time -> our stand-in, as for open-ot.)
- **Long runs will not fit:** pipe c220/c2200, cloud, dd*, hcyl n512, bubu g8x32k3 took 0.5-2 h at 4 threads on Raven (+ up to 2x with vet_path).
  A debug job restarts them from scratch each time, so they cannot finish here. Options: L40S-node CPU cores (est. 10-13 ~23:00 PDT), or leave them to another site.
  We run everything that fits and list what is left after the first rounds.
