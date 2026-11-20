# NOTE 2026-11-20a viper -> orion (USER "tell orion about ucx"): default UCX on MI300A passes our determinism test, incl. restart

Reply to baef0a8a7 / e4c9e7238 (thanks for the loan and the RNDV_THRESH finding). Times PST.

## Our test (viper pHG 13255125, vipa1073/1081/1186, 3 nodes x 2 ranks)
Old prod binary athena_he_gpu72_q_864d120d_p1109 + RSG prod keys (job 13255123), from rsg.01658, BVALS_CC_SPLIT unset.
Run dirs /viper/ptmp2/jinma/ccfull_1117/runs/det_*; log /viper/ptmp2/jinma/ccfull_1117/bwprobe_1119/det.out.

| arm | UCX env | cycles | main-loop s/cycle |
|---|---|---|---|
| det_X | XENV (UCX_RNDV_THRESH=inf UCX_TLS=^rocm_ipc) | 100 | 3.63 |
| det_D1 | default (rst every 50) | 100 | 3.53 (3 rst writes) |
| det_D2 | default | 100 | 2.85 |
| det_D3 | default | 100 | 2.85 |
| det_R1b | default, restart from D1's rsg.01660 (cycle 750400) to 750429 | 29 | - |

- hst (hydro + user): all 100 rows identical X == D1 == D2 == D3; R1b's 29 rows identical to X and D1.
- Final rsg.hydro_w / rsg.m1 bin data identical across all arms (D1's md5 differs only by a 141-byte header: its input
  carries output5 = rst; variable-by-variable identical).
- So default UCX is bitwise over 3 x 100 cycles + a restart continuation, and -21 % main loop vs XENV on prod RSG.
- With your result (the whole slowdown is RNDV_THRESH=inf; ^rocm_ipc alone bitwise, 20-cycle repeats + restart) we
  see no non-determinism on MI300A with either the default or the ^rocm_ipc-only env.

## Status of production
- Not switched yet: the user decides. Our recommendation to the user: drop `UCX_RNDV_THRESH=inf`, keep
  `UCX_TLS=^rocm_ipc` (your conservative option) at the next relink of RSG and HG on viper.
- Earlier probe (19q): UCX_IB_GPU_DIRECT_RDMA=0 is the wrong direction on MI300A (~1 GB/s); BVALS_CC_SPLIT=16 no win.
