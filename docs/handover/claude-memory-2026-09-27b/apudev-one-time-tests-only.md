---
name: apudev-one-time-tests-only
description: viper SLURM partition rule - apudev is for one-off short tests only; production and chained (dependency / auto-resubmit) jobs go to apu
metadata:
  type: feedback
---

On viper, submit to the `apudev` partition ONLY for one-time, short test jobs (smoke tests,
A/B checks, bisection probes that finish in one allocation). Never submit chained jobs
(dependency chains, self-resubmitting restart chains, ensembles meant to run to completion)
to `apudev`. Those go to `apu`.

**Why:** the user said so explicitly on 2026-09-10 ("stop using apudev for chained jobs. only
use it for one-time test"). apudev is a shared development queue; occupying it with chains
blocks other people's quick tests.

**How to apply:** any sbatch script with `--dependency`, a resubmit-on-exit loop, or a
production-length run must carry `--partition=apu`. When in doubt, use apu. See
[[viper-hip-build-recipe]] and [[delegate-simple-tasks-to-opus]] (the Opus agent that
writes submit scripts must be told this rule).

**REPEATED by the user 2026-09-10 (second time): "don't run chained long jobs on apudev".** No
`--dependency` chains, no self-resubmitting submit.sh, nothing that needs more than one 15-min link
on apudev, ever. A single short one-off test (a restart-writer, a few-hundred-cycle twin) is fine.
Everything that continues goes to `apu`.
