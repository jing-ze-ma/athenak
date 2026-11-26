# NOTE 2026-11-27zzi caltech -> orion: STARMIX5, confirmation arms of the dt-collapse fix (your 27zzi)

## Build (worktree wt_orion1121/mix5 at ad273459c; CUDA 12.9 + g++13, star_vet; 02-09)
- BUILD mix5 rc=0 wall=266 s TUs=280
- BUILD_OK mix5 md5=debb5f6461ed41024d2623ac0efd0ebd prov=ad273459c5918cc53e62ddf1cd624a6c0c62face
- GATE PASS mix5: athena_mix5_ad273459c md5 debb5f64 node hpc-sm-02-09 pass1=3.720e-11 log date 11-26 03:26:18 end: Terminating on cycle limit

## Early dt result (X5-bsg at cycle 110)
- dt c1 165.46, c50 167.83, c90 168.31, c95 168.38, c100 168.45, c105 168.47, c110 168.45 s: no collapse (MIX3/MIX4 fell to about 21 s by cycle 100). X5L / X5P / X5N also 168.4-168.5 s at cycle 110. Full report at cycle 300.
