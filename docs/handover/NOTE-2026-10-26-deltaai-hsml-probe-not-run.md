# NOTE 2026-10-26 DeltaAI -> Delta + viper (USER): Hsml probe of QUESTION 5d602691 NOT run on DeltaAI. It cannot start before Delta's diag job. Probe ready, job3 still HELD

**Why not:**
- ghx4-interactive allows **1 job per user**, and the held job3 (3462816) fills that slot. sbatch failed with QOSMaxSubmitJobPerUserLimit.
- The regular ghx4 partition estimates a start ~5.5 days from now (2026-11-01).
- Freeing the slot would mean cancelling job3. Even then, interactive waits today were 5-7 h, which is no sooner than Delta's diag job (est. ~17:54 CDT).
- USER: don't run. 0 SU used.

**Ready here if wanted later (local only, not in git):**
- Probe binary: g2m d9c0152bce + an init print, md5 31b1e582.
  - The print sits right before and right after the init `density();`, for every gas cell with Hsml > 0.2 box: task, ID, position, Hsml/box, plus a per-task count and a total.
  - It prints both the initial guess and the converged value, the same pair as Delta's V1.
- Job script: run 1 at 576 ranks (job3 layout, 288/node, 72/island) and run 2 at 256 ranks (128/node). Each stops at the first Sync-Point (PASS) or TERMINATE (FAIL). 20-min limit, about 3-6 SU. The dry run passed.
- It needs the interactive slot. Running it means cancelling job3 first, or waiting until job3 has run. USER decides.

**job3 3462816:** HELD (JobHeldUser), 0 SU. It is released only on viper's go. If the fix changes the binary or the deck, I cancel it and resubmit only after the USER says go.
