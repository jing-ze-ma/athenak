# TASK 2026-11-15 viper -> Orion: restyle_m3 (40 bitwise A/B comparisons, restyle)

Self-contained CPU task: one-rank runs, two binaries per problem generator, outputs
compared byte for byte.  Every line must come back IDENTICAL.

## Code
Branch `fork/restyle-1114`: A = c9d1a898e (reference), B = 349870ed6 (candidate).
Problem generators: built_in_pgens cs_test deep_hot_jupiter_rt he_star_m1 red_giant sp_test .  For each <prob> and each X in A B (one tree per commit,
one build dir per prob):
```
git worktree add --detach <wt_X> <sha_X>; git -C <wt_X> submodule update --init
module purge; module load gcc/13 openmpi/4.1 cmake/3.30
cmake -S <wt_X> -B <wt_X>/build_<prob> -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release \
      -D PROBLEM=<prob>
make -C <wt_X>/build_<prob> -j32
cp <wt_X>/build_<prob>/src/athena <bindir>/athena_<prob>_X
```
Both commits must be built with the same compiler and flags.

## Runs
Package `docs/handover/task_1115_restyle_m3/`: suite.txt (name|prob|deck|keys|links, @PKG@ = the package dir), decks/,
data/ (tables named by the decks), tools/ (mkdeck.py, bincmp.py, cmp_ab.sh).
Copy the directory somewhere writable, then
```
./run_all.sh <bindir> 16        # TMO=<s> overrides the 7200 s per-run timeout
```
Use sbatch/srun if Orion requires it (one core per run; the he_star decks need ~20 GB
and up to ~1 h each).

## What to send back (NOTE-2026-11-15-orion-viper-restyle_m3.md on fork/bsg-files-1009)
1. table.txt (one line per comparison: IDENTICAL / DIFFER / FAIL, files, cycles).
2. For every DIFFER or FAIL: the differing file names and the last 20 lines of both
   runs' stdout.txt.
3. The compiler and module versions used.
