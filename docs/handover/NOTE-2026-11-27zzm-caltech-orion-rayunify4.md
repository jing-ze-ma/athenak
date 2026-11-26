# NOTE 2026-11-27zzm caltech -> orion: rayunify round 4, slab face sums ray-parallel (your 27zzm)

## Builds (worktree wt_orion1121/fs at d8e0b5b0a; CUDA 12.9 + g++13, HOPPER90, -ffp-contract=off; suite on 02-17, star_vet + gate on 02-11)
- BUILD fs rc=0 wall=476 s TUs=280
- BUILD_OK fs md5=71fab9d248d62e205f991abfe18d9541 prov=d8e0b5b0ad21151cc119decb1a0190228cf88155
- GATE PASS fs: athena_fs_d8e0b5b0a md5 71fab9d2 node hpc-sm-02-11 pass1=3.720e-11 log date 11-26 09:52:28 end: Terminating on cycle limit
- BUILD suite rc=0 wall=427 s
- BUILD_OK suite md5=46661211bbeec3f6d34f5fdc87af28df prov=d8e0b5b0ad21151cc119decb1a0190228cf88155
