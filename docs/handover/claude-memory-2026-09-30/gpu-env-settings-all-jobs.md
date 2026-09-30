---
name: gpu-env-settings-all-jobs
description: User rule 09-23 - every GPU job submission (production, test, agent jobs) exports the validated ROCm env settings; currently HSA_NO_SCRATCH_RECLAIM=1 (+ HSA_XNACK=1); add newly validated ones everywhere
metadata:
  type: feedback
---

User 2026-09-23: "if you find the right combination, remember to add those to all submissions".
Current validated set (bitwise-identical output, measured on GPU):
- `export HSA_XNACK=1` (existing)
- `export HSA_NO_SCRATCH_RECLAIM=1` -- prod4 rot-283 restart 22.5 -> 17.0 s GPU time (1.32x), identical
  output (tests_m1/runs_3k_gpu3d/README_HALO.md). Keeps GPU scratch allocated between kernels; ROCm
  runtime setting, not APU-specific.
ENV SWEEP DONE 09-23 (tests_env_tuning/README.md): NOTHING NEW to adopt. Baseline 17.23 +- 0.15 s (9 runs).
Within noise: HIP_FORCE_DEV_KERNARG, HSA_ENABLE_INTERRUPT=0, GPU_MAX_HW_QUEUES, --cpu-bind, UCX_PROTO_ENABLE.
Slower: HSA_ENABLE_SDMA=0 +3.7 %, UCX_RNDV_THRESH=inf +3.8 %, OMPI_MCA_pml=ob1 +71 %.
REJECTED: HSA_XNACK=0 (3 % faster on prod4 but GPU memory access fault in M1; Kokkos requires XNACK on MI300A).
Slurm default binding already NUMA-correct; binary has no OpenMP. Old pending line was:
binding, GPU_MAX_HW_QUEUES, UCX. Add every setting it validates HERE and in .claude/agents/worker.md.

**Why:** these are free speed-ups with identical output; forgetting them in one submit script silently
costs 25-30 % of GPU time.
**How to apply:** every sbatch script written by me or an agent exports the full set; production submit.sh
files too. Already-queued jobs keep the environment of their snapshot: applying a new setting to queued
production links = edit submit.sh + cancel/resubmit the pending links (done 09-23 for prod4/hyd4 links
2-4 with the user's OK). See [[cs-mhd-prod4-run]], [[cs-hyd4-prod-run]], [[measure-cost-on-gpu]].
