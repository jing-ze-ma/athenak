---
name: measure-cost-on-gpu
description: User 2026-09-22 - performance/cost verdicts must be measured on the GPU (apudev), not on CPU; production runs are GPU
metadata:
  type: feedback
---

User, 2026-09-22 ~09:00: "you should always try the cost on gpus since the production runs are gpus."

**Why:** the CPU bottleneck is not the GPU bottleneck (the frozen-operator ck solve was judged a miss at
3.3x on CPU where the four B-dependent recurrences dominate; on the APU, table lookups, memory traffic
and the g-point parallelism weigh differently). Every earlier cost number (rt-cost-profile, ck sweep
split) was GPU-measured for the same reason.

**How to apply:** any deliverable whose verdict is a cost ratio (RT solvers, reconstruction, halo work,
sweeps) gets a GPU timing arm on apudev (15 min) from a clean-snapshot HIP build before the verdict is
written; CPU timings are for correctness gates and rough guidance only, and must be labelled as such.
Briefs to agents must say so. See [[apudev-for-short-jobs]], [[rt-cost-profile-0916]].

**Tightened by the user 09-22 ~23:40: "from now on we should mainly do timing comparison on GPUs; all the
productions will be on GPUs".** Every timing / speed-up / cost comparison is measured on the GPU (apudev, same
binary, interleaved arms, repeats). CPU runs are for correctness gates (bitwise, convergence) only; do not
report CPU timings as results, at most as a side note. Put this in every agent brief that has a cost gate.
