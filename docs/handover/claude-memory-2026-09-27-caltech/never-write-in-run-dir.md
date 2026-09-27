---
name: never-write-in-run-dir
description: run/ holds 165 GB of simulation output and must never be written to; read-only is fine
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 727372a4-3557-4a71-b294-5104e4fbe37f
  modified: 2026-08-16T16:18:00.915Z
---

`/viper/u2/jinma/ATHENAK/athenak/run/` holds ~165 GB of simulation output and plots across 32 entries. Never create, modify, or delete anything inside it. Reading files there is fine and sometimes required (e.g. `run/2500_test/ohm_128/submit_viper.sh` is the reference SLURM submission script).

**Why:** the user asked for this explicitly at the very start of the session, before anything else. It is irreplaceable output from long GPU runs.

**How to apply:** Put scratch run directories in `/viper/u2/jinma/ATHENAK/bench/` or the session scratchpad instead. Note `run/` is untracked but NOT in `.gitignore`, so never use `git add -A`/`git stash -u` in this repo — they would try to swallow all 165 GB. Stage files by explicit path. A sparse-checkout rule (`/*`, `!/run/`) is configured to keep git from writing the handful of stray PNGs the branch tracks under `run/`; leave it in place.
