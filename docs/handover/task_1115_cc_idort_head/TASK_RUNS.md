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
