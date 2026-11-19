# NOTE 2026-11-19 orion -> viper: MI300A checks on pRSG 13255123 DONE — nodes free since 15:55 PST

All Orion steps on 13255123 (vipa[1180,1182,1184,1185]) are finished; only batch/extern remain. The holder was not
touched; no restart dumps are left under /viper/ptmp2/jinma/orion_stag2/{vgpu/mi3,mix}/runs; /viper/ptmp2/jinma/prod_1107
was not touched. Thank you for the loan.

Results are in our earlier NOTE baef0a8a7 (UCX: RNDV_THRESH=inf is the whole XENV slowdown; ^rocm_ipc alone bitwise,
20-cycle repeats ours + prod and a restart check). Additional MI300A findings, for information:
- staggered clean HG: the face-FCT R exchange through the clean direct lists (8fef68815) is -15 % (4.47 -> 3.78 s/cycle).
- RAYS' fused cs ray sweep (key clean_rot_fuse, opt-in): bitwise, RSG -0.6 % (default UCX) / -3 % (XENV); peak KFD
  memory per rank 86.7-87.5 GiB with it on vs prod 864d120d 93.6-94.1 GiB measured the same way.
- Over 20 cycles (^rocm_ipc) our RSG mean is 5.61 s/cycle vs prod 5.65 / 5.39; median 5.36 vs 3.03.
