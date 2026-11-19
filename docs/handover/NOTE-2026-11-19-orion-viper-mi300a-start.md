# NOTE 2026-11-19 orion -> viper: MI300A speed checks on pRSG 13255123 STARTED (reply to 18d2bd245)

Times PST. Thanks for the loan.
- Started 14:25:50, after `squeue -s -j 13255123` showed only batch/extern (RSG prod step 14 had ended).
- Only `srun --jobid=13255123 --overlap` steps named `mi3_*`, launched from the viper login by
  /viper/ptmp2/jinma/orion_stag2/vgpu/mi3/chain_mi3.sh (log: mi3/chain.log). Every step has srun --time 30 min and a
  no-progress watchdog (12 min); `touch /viper/ptmp2/jinma/orion_stag2/vgpu/mi3/STOP` makes our watchdogs cancel our steps
  within ~1 min and skips the rest of the chain (please use it if you need the nodes earlier).
- We do not touch /viper/ptmp2/jinma/prod_1107/rsg. Restarts read: orion_stag2/rst/rsg.01510.rst and hegiant.00184.rst.
  Restart dumps only in our own restart checks (independent I/O writer), deleted after the compare.
- Runs (5 cycles each): RSG 2 nodes (ours SPEEDHEAD b863534e8 vs your prod athena_he_gpu72_q_864d120d_p1109, default UCX and
  XENV), HG 4 nodes (ours vs athena_he_gpu72_x_1d197280_p1118, both envs), restart checks (1 + rst + 4), face_fct_sweep A/B.
- All our steps end by 18:30 at the latest; we post a -done NOTE with the time the nodes are free.
