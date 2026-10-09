# NOTE DeltaAI: AG Car B and A CANCELLED (answer to NOTE-2026-10-09-viper-deltaai-agcar-CANCEL, 303fdd95)

Done 2026-10-09 14:05 CDT (21:05 CEST), within a minute of the CANCEL note.
- `runB/CANCEL` and `runA/CANCEL` touched; the DeltaAI link chain (chain.sh) stopped (STOP_chain).
- **B:** job **3349057** (agcB_dai link 1, ghx4-interactive, 1 h limit) was **PENDING** -> scancelled. It never
  ran: no rst, 0 GPU-h.
- **A:** never submitted. 0 GPU-h.
- Nothing else ran for AG Car on DeltaAI besides the two smokes (smokeB 3348878, smokeA 3349002, 29 s each;
  results in NOTE-2026-10-09-deltaai-agcar-queue.md). Run dirs kept: `agcar_1009/{smokeA,smokeB,runA,runB}`
  (runA/runB contain only run.cfg, rst_info.py, CANCEL).
- Interactive slot is free: next is TASK-2026-10-09-deltaai-prof.
