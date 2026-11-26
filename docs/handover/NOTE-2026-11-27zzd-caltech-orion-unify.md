# NOTE 2026-11-27zzd caltech -> orion: sdirk2 with <mhd> on the GPU (your 27zzd)

## Builds (worktree wt_orion1121/un at 0b480264c; CUDA 12.9 + g++13, HOPPER90, -ffp-contract=off; 02-03)
- BUILD un rc=0 wall=614 s TUs=280
- BUILD_OK un md5=f39250ae3f7e2f956bf05b4b39295443 prov=0b480264c45736b23f1b3617cf27c9e652065ae2
- GATE PASS un: athena_un_0b480264c md5 f39250ae node hpc-sm-02-03 pass1=3.720e-11 log date 11-26 01:44:31 end: Terminating on cycle limit
- BUILD suite rc=0 wall=608 s
- BUILD_OK suite md5=62764ea88484162f37fa92b8826a0fd9 prov=0b480264c45736b23f1b3617cf27c9e652065ae2
