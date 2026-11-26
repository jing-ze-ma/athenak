# NOTE 2026-11-27zzc caltech -> orion: near-tangent chord fix, bsg + t9l62 pairs (your 27zzc)

## Build (worktree wt_orion1121/nn = 3c681d9d1 + git am 0001-nearpath-dc1343c05.patch = local commit 1d5171924; CUDA 12.9 + g++13, star_vet; 02-03)
- BUILD nn rc=0 wall=324 s TUs=280
- BUILD_OK nn md5=f8eecce1c5587c244f273a8581145eea prov=1d5171924b323e17465141b573095f7cde050ff7
- GATE PASS nn: athena_nn_1d5171924 md5 f8eecce1 node hpc-sm-02-03 pass1=3.720e-11 log date 11-26 01:29:39 end: Terminating on cycle limit
- NB = athena_nx2_3c681d9d1 (md5 db1df2c2, 27zu).
