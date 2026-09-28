---
name: cpu-checks-on-login-node
description: user 09-28: CPU-only checks (CPU builds, A/B bitwise, tst *_cpu/*_mpicpu, small physics tests) run on the viper LOGIN node, niced, <=16 cores, <=40 GB -- never as apu/apudev jobs
metadata:
  type: feedback
---
CPU-only work runs directly on the login node (viper13: 256 cores, claude-work.slice 70 GB cap):
`nice -n 10`, at most ~16 cores / 16 MPI ranks, RSS <= 40 GB, short: only checks that take <= ~30 min wall on the login node (user: "if they don't take too long"); longer CPU work stays in the GPU queue.
GPU nodes (apu/apudev) are only for GPU runs: timing, GPU bitwise/determinism, production-size runs.
**Why:** 09-28 night CPU-only gates (d2Bcpu, efix_tst) sat hours in the full apu queue while the login node was
idle; CPU partitions small/general are denied to this account ([[viper-cpu-jobs-small-partition]]).
**How to apply:** put this rule in every worker brief that includes CPU gates.
