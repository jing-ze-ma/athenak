---
name: job-watch-10min
description: user 09-24 - background job watchers poll every 10 min (sleep 600), not faster
metadata:
  type: feedback
---

User (09-24): "I think you can let the monitoring be 10 min probably."
When waiting on Slurm jobs, use ONE background watcher (`while squeue -h -j ... | grep -q .; do sleep 600; done; sacct ...`,
run_in_background) that notifies on completion; poll interval 10 min. No foreground waits, no faster polling.

**Why:** jobs run 5-60 min; faster polling costs nothing useful.
**How to apply:** for agents' GPU/CPU jobs that nothing else will announce (agent stopped after submitting),
arm a 600 s watcher, then resume the agent / analyse. See [[no-blocking-waits]], [[no-heavy-monitoring]],
[[housekeeping-minimal-tokens]].
