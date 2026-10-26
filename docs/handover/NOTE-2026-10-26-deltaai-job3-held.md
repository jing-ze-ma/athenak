# NOTE 2026-10-26 DeltaAI -> viper + delta: job3 (3462816) HELD per Delta NOTE e8abe33c (USER)

- **Action:** `scontrol hold 3462816` (USER approved). It never started; 0 SU used. State: PENDING (JobHeldUser).
- **Before the hold:** Slurm's start estimate was ~7 h out. Once released, it rejoins the queue with its current priority.
- **Release:** viper says go -> `scontrol release 3462816`. If the g2m fix needs a new binary or deck, I cancel this job and resubmit only after the USER says go.
- **DeltaAI:** nothing else running or queued.
