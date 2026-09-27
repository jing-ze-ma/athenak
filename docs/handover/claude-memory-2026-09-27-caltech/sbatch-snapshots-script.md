---
name: sbatch-snapshots-script
description: Slurm copies the batch script at submission - editing submit.sh (BIN=...) after sbatch does NOT change a queued job; the 09-15 prod_m3 chain ran the OLD binary for 200 cycles because of this. Always cancel + resubmit after editing a submit script; verify with `scontrol write batch_script <id> -`
metadata:
  type: feedback
---

On 09-15 the prod_m3 chain (11709068-70, submitted 09-14) was released after I edited submit.sh to point at the new binary; the stored script still had BIN=athena_pinned_final, so link 1 ran the stale binary (no PCR line in out.txt, 695 ms/cycle, race signature). The user had to cancel; resubmitted as 11710002 -> 11710003 -> 11710004, losing the queue position. The stale outputs (incl. an rst at t=0 that the restart-aware script would have picked up) are in prod_m3/stale_run_oldbinary_0915.
**Why:** sbatch snapshots the script; only files read at RUN time (the athinput via -i, the binary file contents at the stored path) follow later edits - and the binary path in the script is frozen.
**How to apply:** after any edit to a submit script of a queued job: scancel + resubmit, then `scontrol write batch_script <id> - | grep BIN`; on a restart-aware chain also clear stale rst/ before resubmitting. Verify the startup line of the intended feature (e.g. 'solver pcr nseg 32') in the first minute of every production link.
