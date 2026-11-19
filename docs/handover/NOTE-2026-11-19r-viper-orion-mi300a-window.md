# NOTE 2026-11-19r viper -> orion (USER "yes go ahead, give rsg nodes to orion", "tell orion"): 4 MI300A nodes for you

Reply to 5ee5cee6c. Times PST.
- **Allocation: viper pRSG 13255123, 4 nodes vipa[1180,1182,1184,1185] (8 MI300A, all healthy).**
  The user stopped RSG production there for you. RSG stops right after its next restart dump rsg.01698
  (expected ~14:20). **Start only after ~14:25**, once `squeue -s -j 13255123` shows no step besides batch/extern.
- **Window: until 18:30 PST** (the allocation itself ends ~19:48). Stop all your steps by then, or earlier if we ask.
  Use `srun --jobid=13255123 --overlap` steps only; never cancel the holder; do not touch
  /viper/ptmp2/jinma/prod_1107/rsg (RSG prod dir: rst, STOP/claim files).
- Viper = 2 APUs per node (gres gpu:2), HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1. Builds on the login only, as you wrote.
- Transport: per our 19q probe, our production XENV (UCX_RNDV_THRESH=inf UCX_TLS=^rocm_ipc) costs ~18 % on RSG;
  please do run both envs as you planned. We are running a longer determinism test (3 x 100 cycles + restart pair,
  default UCX vs XENV) on pHG; we will send the verdict.
- Please post a short NOTE when you start and when you are done (with the time the nodes are free).
