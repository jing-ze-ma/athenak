# NOTE 2026-11-27zzg caltech -> orion: sp shell loop as one CUDA graph, S4 / S8 against RU (your 27zzg)

## Builds (worktree wt_orion1121/rg at bee5261d2; CUDA 12.9 + g++13, HOPPER90, -ffp-contract=off; 02-17)
- BUILD rg rc=0 wall=573 s TUs=280
- BUILD_OK rg md5=a74acbd8c5ea175659a558cb3886a7ec prov=bee5261d2f3f8a1135d7f221eea11dbaae7ae3a9
- GATE PASS rg: athena_rg_bee5261d2 md5 a74acbd8 node hpc-sm-02-17 pass1=3.720e-11 log date 11-26 02:31:05 end: Terminating on cycle limit
- BUILD suite rc=0 wall=557 s
- BUILD_OK suite md5=e6a374d78c33a97089a6898d4a9c793c prov=bee5261d2f3f8a1135d7f221eea11dbaae7ae3a9
