# TASK 2026-11-15 viper -> Orion: cc_fixedD (156 radiation-wave runs, clean M1 scheme)

Self-contained CPU task: one-rank singletons, no MPI ranks needed.

## Code
Branch `fork/clean-cc`, commit 87e045577. Build (built-in pgens, no PROBLEM):
```
git worktree add --detach <wt> 87e045577
module purge; module load gcc/13 openmpi/4.1 cmake/3.30
cmake -S <wt> -B <wt>/build_cc -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release \
      -D CMAKE_CXX_COMPILER=mpicxx -D CMAKE_CXX_FLAGS=-ffp-contract=off
make -C <wt>/build_cc -j32
```

## Runs
Package `docs/handover/task_1115_cc_fixedD/`: decks in runs/<case>/rw.athinput, runs.txt (one command per run, relative to
the package), cases.txt, scorer rwfarm.py (numpy only; references from rw/radwave_disp.py).
Copy the directory somewhere writable, then:
```
./run_all.sh <wt>/build_cc/src/athena 32     # runs all, skips finished ones, then scores
```
Use sbatch/srun if Orion requires it (one core per run). Runs take seconds to ~30 min each.

## What to send back (NOTE-2026-11-15-orion-viper-cc_fixedD.md on fork/bsg-files-1009)
1. table.md and score.txt (L1 errors, orders, box-total energy / momentum drift per run).
2. clean_lines.txt (last CLEAN line, E-floor energy line and rc per run).
3. Any failed run (rc != 0): the last 20 lines of its log.txt.
Also commit the runs/*/farm.json files (small) under docs/handover/task_1115_cc_fixedD/results/ if convenient.
