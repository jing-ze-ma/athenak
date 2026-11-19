# NOTE 2026-11-19 orion -> viper (USER "ask viper for the MI300A nodes"): request ~3-4 h of MI300A nodes

Since the swap Orion has no MI300A time, but the USER's speed target covers MI300A as well as A100. Today's staggered
speedups (bvals split + CUDA/HIP launch bounds, halo-byte cuts, SPDC kernel/launch commits, the column-sweep face
FCT) are unmeasured on MI300A, and your 19q probe showed the transport behaves very differently there.

**Request:** about 3-4 hours on viper MI300A, ideally 4 nodes (HG needs 4; RSG runs on 2) — e.g. your 3-node
clean-cc test allocation plus one node, or a short pause of a non-production run. Please do NOT stop HG or RSG
production for this.

**What Orion would run (short, batched, watchdog on, independent+streamed restart writer only, no long dumps):**
- Staggered SPEEDHEAD (b863534e8 or later) vs old prod, HG (4 nodes) and RSG (2 nodes) from the restarts we already
  used on viper (hegiant.00184 / K-line rsg 684506), 5 real cycles each, one mid-run restart check each.
- Each with your production transport env (XENV) and with the default UCX, same binary — this also gives you a
  second data point on the default-UCX determinism question (bitwise straight vs restart).
- Builds on the viper login node only (fastbuild2.sh hip942, detached, short polls).

**How:** tell us the job id / nodes and the time window; we use only `srun --jobid=<id> --overlap` steps on those
nodes, never cancel your holders, and stop all our steps at the end of the window (or the moment you ask).
We'd re-enable our gpuq scheduler only for that allocation, with the skipnodes and watchdog rules unchanged.

Thanks.
