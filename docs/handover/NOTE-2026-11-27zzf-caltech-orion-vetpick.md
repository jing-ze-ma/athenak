# NOTE 2026-11-27zzf caltech -> orion: VET-solve kernel pick, batched 1-GPU timing (your 27zzf)

## Builds (CUDA 12.9 + g++13, star_vet, HOPPER90, -ffp-contract=off; 02-09). C and D are the B worktree (vpb) built in separate build dirs with -D Athena_EXTRA_DEFS=ATHENA_LB_RHS=2 / ATHENA_CUDA_BODY_GLOBAL=512
- BUILD vpa rc=0 wall=1886 s TUs=280
- BUILD_OK vpa md5=1cf9bfdea9a06c5cc568c93a293eefc8 prov=f42d790321f8581c64aa4731dd41309320c3e453
- GATE PASS vpa: athena_vpa_f42d79032 md5 1cf9bfde node hpc-sm-02-09 pass1=3.720e-11 log date 11-26 02:42:46 end: Terminating on cycle limit
- BUILD vpb rc=0 wall=1921 s TUs=280
- BUILD_OK vpb md5=a45c5fd1901f2313b6e47002ec2b9cfa prov=0fbba37d48721da7b6eb8117185fffe55e137035
- GATE PASS vpb: athena_vpb_0fbba37d4 md5 a45c5fd1 node hpc-sm-02-09 pass1=3.720e-11 log date 11-26 02:43:26 end: Terminating on cycle limit
- BUILD vpc rc=0 wall=1907 s TUs=280
- BUILD_OK vpc md5=bbc5be4c03c589662f1b0598c4dbf075 prov=0fbba37d48721da7b6eb8117185fffe55e137035
- GATE PASS vpc: athena_vpc_0fbba37d4 md5 bbc5be4c node hpc-sm-02-09 pass1=3.720e-11 log date 11-26 02:43:10 end: Terminating on cycle limit
- BUILD vpd rc=0 wall=1889 s TUs=280
- BUILD_OK vpd md5=75a656644efe3bca54182fd4c3ff6a20 prov=0fbba37d48721da7b6eb8117185fffe55e137035
- GATE PASS vpd: athena_vpd_0fbba37d4 md5 75a65664 node hpc-sm-02-09 pass1=3.720e-11 log date 11-26 02:42:47 end: Terminating on cycle limit
