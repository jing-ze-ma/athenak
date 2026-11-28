# NOTE 2026-11-27zzw caltech -> orion: EXACT3, the cone cut's options on t9l62 / bsg (your 27zzw)

## Build (worktree wt_orion1121/ex3 at 44fb4d8ca; CUDA 12.9 + g++13, star_vet; 02-17)
- BUILD ex3 rc=0 wall=528 s TUs=280
- BUILD_OK ex3 md5=aeadd2cff903dde0ff41504adc079a7f prov=44fb4d8ca0d1becf7c2acffac148f1c661fee7d9
- GATE PASS ex3: athena_ex3_44fb4d8ca md5 aeadd2cf node hpc-sm-02-17 pass1=3.720e-11 log date 11-26 15:05:06 end: Terminating on cycle limit

## Status
- All four arms ran their 300 cycles on 02-17 (11-26, ended 15:51), binary athena_ex3_44fb4d8ca md5 aeadd2cf: E3-t9, E3F-t9 and E3FE-t9 on 4 H200 each, E3FE-bsg on 1 H200. Each ends with 'Terminating on cycle limit'.
- No readings yet. All holders have ended and analysis is not allowed on the login node. Items 1-6 (script lzab_1127/e3an.py) follow once a compute allocation exists. Logs and dumps are kept in lzab_1127/arms/E3*.
