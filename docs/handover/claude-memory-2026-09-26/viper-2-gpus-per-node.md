---
name: viper-2-gpus-per-node
description: RECURRING MISTAKE - every Viper apu/apudev node has exactly 2 GPUs (MI300A); apudev max 1 node; never gpu:4 or 4 ranks/node
metadata:
  type: feedback
---
Every Viper GPU node, in both the apu and the apudev partition, has gres gpu:2, 96 CPUs and 220 GB (checked with sinfo on 09-26). There are NO 4-GPU nodes.
- apudev JOB_SIZE is 1 node, so at most 2 GPUs on apudev.
- More than 2 GPUs -> apu, N nodes x 2: `-p apu --nodes=N --ntasks-per-node=2 --gres=gpu:2 --constraint=apu --cpus-per-task=24`.
- 8 GPUs = 4 nodes. `--gres=gpu:4` is rejected as "Invalid feature specification" or sits pending forever.

**Why:** the user said on 09-26 "this is a recurring problem for you". Job 11980276 (gpu:4, 2 nodes on apudev) could never run, and I then wrongly told the user that apudev nodes have 4 GPUs.
**How to apply:** before every sbatch, mine or an agent's, check the header against this. Put this rule in every agent brief that submits GPU jobs. Also see [[apudev-for-short-jobs]] and [[gpu-env-settings-all-jobs]].
