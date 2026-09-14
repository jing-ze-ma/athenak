---
name: apudev-for-short-jobs
description: USER PREFERENCE (2026-09-15) - any job that finishes within 15 minutes should go on the apudev partition (-p apudev, --time 00:15:00, QOS cap 15 min; same --constraint=apu --gres=gpu:N, HSA_XNACK=1) instead of waiting on apu priority; builds, narrow-box 1e4-s arms (~8-10 min), 1-D 300-s arms, cost-box 200-cycle timings qualify. Longer jobs stay on apu.
metadata:
  type: feedback
---
**Why:** apu priority waits of 10-40 min were dominating the turnaround of small gate arms.
**How to apply:** in every agent brief that submits Slurm jobs, state the apudev rule for sub-15-min jobs.
