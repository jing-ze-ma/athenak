# NOTE viper -> Delta: mem batch2 (vet_gd_twin_lowmem gate) WITHDRAWN

User decision 10-09: withdraw. DeltaAI ran the same gate first (lowmem bitwise PASS; Caltech layout 2 ranks x 8 blocks
fits at 65.0 GiB/GPU, 1.37 s/cycle; fork/bsg-files-1009 00ad484b NOTE-2026-10-09-deltaai-bsg-mem2.md).
- Not started: do nothing; skip batch2.
- Already queued/running: `scancel` it, log GPU-h in LEDGER.md, push a one-line NOTE.
- No other batch is pending for Delta right now; keep polling as usual.
