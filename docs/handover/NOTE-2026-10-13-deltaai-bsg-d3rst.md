# NOTE DeltaAI: swap the running BSG chain to d3rst-1013 @ 9ceb9cdc (TASK-2026-10-13-deltaai-bsg-d3rst) -- PICKED UP

**Picked up** 19:47 CDT (TASK 63b07f61). The DeltaAI chain is between links (link 2 ended t 9.38 d; link 3 pending, est. ~01:00-02:40 CDT).
- Build: d3rst-1013 @ 9ceb9cdc, same flags as d3bf42d7 (build_inc_deltaai.sh hes_gpu); running now (199 objects). md5 follows.
- Step 2 runs INSIDE the next link job (ghx4-interactive allows one job per user, and a separate batch job would wait ~1 day):
  scratch dir bsgprod_1012/d3rstcheck/<job>, not the run dir, no lock.
  - A: BSG deck fresh, 25 cycles continuous, rst at cycle 12 (`output5/dt=1050`), hst every cycle, dumps at t 0 and 2187.5.
  - B: new srun from A's rst 00001 to cycle 25.
  - PASS = cycle-25 dumps byte-identical (cmp) A vs B, user.hst rows 13-25 identical, and A's cycle-25 dump payload identical to
    the 10-12 d3bf42d7 smoke (same deck; paycmp, header excluded).
  - On PASS the job points BINARY at 9ceb9cdc and the SAME job continues the production link from the newest rst (the one cold
    slope start); on FAIL it keeps d3bf42d7 and I report.
- The pending link 3 (3368862, submitted with the old script) is resubmitted with the hooked script once the build is done (its queue age is lost).
  The 12 h ghx4 job 3367703 reads BINARY at start, so it runs d3bf42d7 before the swap and 9ceb9cdc after it.

- 19:50 CDT: built **athena_hes_gpu_9ceb9cdc0ed4, md5 54c0a97171c33a72a6dd8f11c3a74bfb**. Swap armed (SWAP_PENDING): the next interactive link
  runs the scratch restart check first and swaps only on PASS. That link is job **3369777** (est. start ~07:35 CDT 10-13; the previous
  link-3 job 3368862 was cancelled while pending because it had the old script). The 12 h ghx4 job 3367703 has the old script: if it
  starts first it runs d3bf42d7 (no check), and the swap waits for the next interactive link.
