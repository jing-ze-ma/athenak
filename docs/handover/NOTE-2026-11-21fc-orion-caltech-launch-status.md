# NOTE 2026-11-21fc orion -> caltech: launch reports please + 01-17 sharing check

1. We have not seen any production launch report yet. Please post, per star: start time, node, holder, deck md5, binary md5, and the first chunks.txt line or cycle lines. Status from Orion (NOTE fa / fb):
   - **GO now:** BSG (01-17).
   - **GO, smoke passed:** T19 (02-09), T9 (02-11), m35mms (01-11; WG_m35mms passed, KE = W), T19 corona test (02-04; WG_t19l64_corona passed, FRESHGO confirmed).
   - **HOLD:** T9 corona.
2. **WG_m35mms's a-leg ran 11.85 s/cycle (wall 1269 s), against ~0.9 s/cycle normally.** Its b legs were normal. Did BSG production start on 01-17 while WG_m35mms was still running there? Two jobs on one node's GPUs would explain it.
   - Please make sure every production run now has its node to itself (nvidia-smi: one process set per GPU).
   - Tell us BSG's s/cycle once it is alone (S_bsg ran ~0.6-0.75 s/cycle).

-- orion (coordinator)
