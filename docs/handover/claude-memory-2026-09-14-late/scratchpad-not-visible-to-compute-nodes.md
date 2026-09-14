---
name: scratchpad-not-visible-to-compute-nodes
description: "The session scratchpad under /tmp is node-local, so SLURM jobs cannot see it; stage batch work on /orion instead"
metadata: 
  node_type: memory
  type: project
  originSessionId: d0224e66-3587-4359-a545-9608682177c3
  modified: 2026-08-07T14:13:58.477Z
---

On orion, `/tmp` is **node-local**. The session scratchpad (`/tmp/claude-.../scratchpad`) exists only on the login node, so anything submitted with `sbatch` from there fails: the batch script runs but every `srun` step dies immediately (exit 1) and even the `#SBATCH -o` log file cannot be written, leaving no diagnostics.

**How to apply:** stage input files, submit scripts, and run directories for SLURM jobs on shared storage under `/orion/u/jinma/...` (e.g. a scratch dir alongside the repo such as `/orion/u/jinma/ATHENAK/eostest`), not in the scratchpad. Keep the scratchpad for login-node-only work. Do not put scratch runs inside the repo, and never inside [[run-directory-untouchable]].

Symptom to recognize: `sacct` shows the job `COMPLETED` in ~10s while every `.0/.1/.2` step shows `FAILED 1:0`, and no log files appear.
