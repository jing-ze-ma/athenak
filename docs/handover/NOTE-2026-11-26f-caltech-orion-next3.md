# NOTE 2026-11-26f caltech -> orion: CUDA suite of caltech-next3-c3f133b44 (your 26f)

## Builds (worktree wt_orion1121/next3; CUDA 12.9 + g++13; 02-16, 4 cores each)
- BUILD next3 rc=0 wall=2201 s TUs=250
- BUILD_OK next3 md5=ad2f2007f6d5d795de7cd73a6ed76db5 prov=c3f133b44b603f7326cac383b08c182088d545e8
- GATE PASS next3: athena_next3_c3f133b44 md5 ad2f2007 node hpc-sm-02-16 pass1=3.686e-11 log date 11-25 10:31:26 end: Terminating on cycle limit
- BUILD suite rc=0 wall=2169 s
- BUILD_OK suite md5=3069594e5eb3b661083052b6a1ca235a prov=c3f133b44b603f7326cac383b08c182088d545e8
- Tiny gate deck as before (clean_face_hctau = false in <rad_m1>).
