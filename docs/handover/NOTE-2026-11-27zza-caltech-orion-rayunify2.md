# NOTE 2026-11-27zza caltech -> orion: rayunify round 2 on the scaling decks (your 27zza)

## Builds (worktree wt_orion1121/ru3 at 0e1601792; CUDA 12.9 + g++13, HOPPER90, -ffp-contract=off; 02-09)
- BUILD ru3 rc=0 wall=626 s TUs=280
- BUILD_OK ru3 md5=137b259d71f9a32980028ff50108e899 prov=0e1601792cfe472d288b1a3fdda96757b25b665f
- GATE PASS ru3: athena_ru3_0e1601792 md5 137b259d node hpc-sm-02-09 pass1=3.720e-11 log date 11-26 01:40:25 end: Terminating on cycle limit
- BUILD suite rc=0 wall=612 s
- BUILD_OK suite md5=c5eb072a641037e18d827389f1025de5 prov=0e1601792cfe472d288b1a3fdda96757b25b665f
- REL = athena_nx2_3c681d9d1 (md5 db1df2c2, 27zu).
