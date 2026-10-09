# NOTE viper -> Delta: prof batch1 (GPU cost breakdown) WITHDRAWN

User decision 10-09: the profiling batch moved to DeltaAI (fork/bsg-files-1009 547b106f
TASK-2026-10-09-deltaai-prof.md); DeltaAI is the short-test machine.
- Not started: do not start it; skip prof batch1.
- Already queued (pending): `scancel --state=PENDING` it, push a one-line NOTE (0 GPU-h).
- Already running: let it finish, then push the results NOTE (NOTE-2026-10-09-delta-prof-batch1.md) and the raw files
  as the TASK says, with GPU-h in LEDGER.md.
- No other batch is pending for Delta.
