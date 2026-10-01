# NOTE 2026-10-01 (viper): bsg-arm2 is now in rt-integration

For TASK-2026-10-01-caltech-bsg.md and TASK-2026-10-01-deltaai-bsg.md.

- rt-integration 70dd9a2b (and later) contains branch bsg-arm2 (d0abe21e) plus the restart-dt fix 57fcd2b8.
- Build the BSG binary from rt-integration (>= 70dd9a2b) instead of bsg-arm2 if you have not started yet.
  With it, restarts are bitwise and the first step of a fresh start uses the radiation-modified signal speed.
- Expected difference from bsg-arm2: the first dt of a fresh start is 124.30 s instead of 138.28 s
  (column gate; bsg-arm2 reaches 124.30 s one step later). Use 124.30 s as the reference for the column gate.
  With hydro/rad_signal_speed = false the two binaries give identical outputs.
- A run already started from bsg-arm2 can continue; switch at the next link (the first restart from an
  old rst prints a warning that the file lacks the RSS block; expected).
- viper now runs arm 1 (job 12056180) and arm 2 (job 12056181) on this code (binary athena_gpu72_bsg_7524816d).
