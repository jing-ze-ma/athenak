# CPU farm for one-off test batches

## Orion (primary): orion_task.sh
Packages runs as a self-contained TASK on fork/bsg-files-1009; Orion replies with a NOTE.
```
./orion_task.sh rw  <caselist> <name> [rev=clean-cc]    # rw/run_radwave.py case lines
./orion_task.sh pkg <dir> <name> <rev> [<rev2> ...]     # any ready package
DRY=1 ./orion_task.sh ...                               # build the package only, push nothing
```
- Each rev must be on the fork or a local branch (pushed as fork/<branch>; never main).
- `pkg`: `<dir>` holds `run_all.sh` (called on Orion as `./run_all.sh <bin> [<bin2> ...] 32`,
  one binary per rev, in order: two revs = A/B or bitwise comparisons) and `TASK_RUNS.md`
  (what the runs are and what to send back). pkg_idort/ is an example.
- `rw`: decks are written by run_radwave.py itself (byte-identical to local runs); the
  package carries rwfarm.py, rw_run_all.sh (as run_all.sh) and rw_TASK_RUNS.md.
- Output: docs/handover/task_<mmdd>_<name>/ + TASK-<date>-orion-<name>.md, committed from a
  scratch worktree /viper/ptmp2/jinma/orion_tasks/wt_<mmdd>_<name> (remove with
  `git worktree remove` when done). poll_orion.sh <name>... saves the reply NOTEs to notes/.

## rwfarm.py (radiation-wave ladders without running locally)
```
python3 rwfarm.py emit  CASELIST ROOT [--rel] > runlist   # decks + one command per run ($ATHENA)
python3 rwfarm.py score CASELIST ROOT                     # -> <tag>/farm.json (+ conservation)
python3 rwfarm.py table CASELIST ROOT                     # space / time ladder tables, orders
```

## Raven (fallback, currently not to be used: user order 11-15 = CPU tests on Orion only)
`ship.sh <rev> [tag]` bundles a local revision, builds none_cpu in a Raven job and prints the
binary; `submit.sh <runlist> <athena> [walltime]` runs an absolute-path run list on 1-4 general
nodes (72 runs at a time per node); `status.sh <runlist>` counts done / failed / running.

## Run lists prepared 11-15
- cases_fixedD.txt (156 = rw/eddS.txt + rw/eddT2.txt), clean_closure = eddington, ref = edd,
  arm keys ts=hesdirk2 vimp=true rm_scheme=clean, both layouts (cell / staggered):
  space ladders N 64/128/256/512, nt = 2N, amp 1e-6, tol 1e-14 / nl 1e-13, maxnl 60,
  P 1/3 and 3.3 x tau 0.1/10/1000 (48 runs); time-only ladders N256 nt 512/1024/2048/4096,
  amp 1e-4 and 1e-5, tol = tol_nl = 1e-15, maxnl 80 (96 runs); plus 12 amp 1e-6 lines
  (P 3.3, tau 10: nt 1024 tol 1e-15 and 1e-16; N4096 nt 512..4096).
- pkg_idort/runs.txt (24): arm CC1 = scheme clean, clean_layout cell, clean_closure extrap,
  clean_maxpic 4, clean_pc_sweeps 2; nltek1/nltek4 x0.5/1/2, hcyl x0.5/1/2, shd/shdg
  x0.25/0.5/1, cb/cbg/ws66 x0.5/1/2 (decks = blend_1022/p2/decks).
