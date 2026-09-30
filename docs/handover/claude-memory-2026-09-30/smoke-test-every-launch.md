---
name: smoke-test-every-launch
description: User rule 09-30: ALWAYS smoke-test on apudev before any production/chained launch or relaunch, with the exact script, binary, input and command-line KEYS
metadata:
  type: feedback
---
Before submitting (or resubmitting after ANY edit to script, input, binary or KEYS) a production or chained job,
run an apudev smoke test (a few cycles, 2 GPUs, same HSA env) that uses the SAME binary, input file and the
SAME command-line override string as the real job, writing into a separate smoke dir. Only submit/release the
real chain when the smoke reaches cycling without FATAL. Applies to agents too: state it in every launch brief.

**Why:** 09-30 the He presn hllc runs 12037023/12037044 died in 9 s after hours in the queue: r2.sh/r3.sh passed
rad_m1/implicit_opac_newton=true but he3d_M1.athinput lacked the key (parser fatal on unknown command-line keys).
A 1-minute apudev smoke would have caught it. User: "make sure this doesn't happen again and always do smoke test".

**How to apply:** run the exact srun line from the job script with -d <smoke dir> time/nlim=~10 appended;
check run.log for FATAL / cycling. The user granted scontrol hold/release (09-30) -> hold queued links while
smoke-testing a fix, release after it passes. Related: [[apudev-for-short-jobs]], [[sbatch-snapshots-script]].
