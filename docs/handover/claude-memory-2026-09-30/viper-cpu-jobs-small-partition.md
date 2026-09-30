---
name: viper-cpu-jobs-small-partition
description: user 09-28: NO CPU partitions on viper for us (small/general denied); only the login node might serve CPU work (unconfirmed); all batch jobs go to apu/apudev
metadata:
  type: feedback
---
We cannot use viper's CPU partitions (small / general refuse this account, also with -A mpa_cpu; user confirmed
09-28: "we cannot use cpu here"). The only CPU we can use is probably the login node (user "not sure").
**Why:** tried 09-28 to move a CPU-only gate off the full apu queue; permission denied.
**How to apply:** never propose small/general. CPU-only work: short tests on the login node within the RAM cap
(claude-work.slice ~70 GB, keep analysis RSS <= 40 GB, nice) or on apu/apudev GPU nodes. On apu, request only the
time actually needed so backfill can place the job (scontrol update JobId=... TimeLimit=... for pending jobs).
apudev = 2 nodes (often 1 down), 15 min. Related: [[viper-2-gpus-per-node]], [[apudev-for-short-jobs]],
[[login-node-ram-cap-pool-size]].
