# TASK 2026-11-15 viper -> Orion: cc_idort_head (24 one-rank CPU runs)

Self-contained CPU task in `docs/handover/task_1115_cc_idort_head/`; one-rank singletons (use sbatch/srun if Orion
requires it, one core per run). Copy the directory somewhere writable and run from there.

## Build (built-in pgens, no PROBLEM; commits on the fork)
```
module purge; module load gcc/13 openmpi/4.1 cmake/3.30
git -C <athenak> worktree add --detach <wt_dc67dba7> dc67dba7b
cmake -S <wt_dc67dba7> -B <wt_dc67dba7>/build -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release \
      -D CMAKE_CXX_COMPILER=mpicxx -D CMAKE_CXX_FLAGS=-ffp-contract=off
make -C <wt_dc67dba7>/build -j32
```

Run: `./run_all.sh <wt_dc67dba7>/build/src/athena 32`

## Runs
IDORT ladder of the cell-centred clean scheme, arm CC1 (rad_m1 keys: scheme = clean,
clean_layout = cell, clean_closure = extrap, clean_maxpic = 4, clean_pc_sweeps = 2).
24 one-rank runs (runs.txt: deck, resolution factor, tlim): nltek1/nltek4 x0.5/1/2,
hcyl x0.5/1/2, shd/shdg x0.25/0.5/1, cb/cbg/ws66 x0.5/1/2. The x2 beam runs are the slowest;
each run has an 8 h timeout.
```
./run_all.sh <wt>/build_<sha8>/src/athena 24     # skips finished runs, then score.txt + clean_lines.txt
```
## What to send back
1. score.txt (L1 error and order per deck; nltek against the exact transport solution).
2. clean_lines.txt (last CLEAN line: marks, BE steps, neg / nreal; wall time; rc per run).
3. Any failed or timed-out run: the last 20 lines of its run.log.

Reply: NOTE-2026-11-15-orion-viper-cc_idort_head.md on fork/bsg-files-1009 with the items above; record
the commit(s) and binary md5s.
