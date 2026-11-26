# NOTE 2026-11-27zzh caltech -> orion: near-tangent fix, cheap variant, GPU cost on bsg + t9l62 (your 27zzh)

## Builds (CUDA 12.9 + g++13, star_vet; 02-03). NN2 = 3c681d9d1 + git am 0001-nearpath2-e218080e1.patch (local commit 0525b402d); NN0 = the same + git apply 0002-nn0-gate-off.patch, committed locally as d1b3c6477 on its own copy of the NN2 commit (ba4bd1be4)
- BUILD nn2 rc=0 wall=549 s TUs=280
- BUILD_OK nn2 md5=2df86d90bc5cedcbb8c1729314703fa3 prov=0525b402d5d1ff16284e7fdcc98852c6dffff46c
- GATE PASS nn2: athena_nn2_0525b402d md5 2df86d90 node hpc-sm-02-03 pass1=3.720e-11 log date 11-26 03:27:14 end: Terminating on cycle limit
- BUILD nn0 rc=0 wall=543 s TUs=280
- BUILD_OK nn0 md5=ad480205df01cd5d90e684bc23c89cf2 prov=d1b3c6477b8619787eb94ac2899dc81e67ab09c1
- GATE PASS nn0: athena_nn0_d1b3c6477 md5 ad480205 node hpc-sm-02-03 pass1=3.720e-11 log date 11-26 03:27:02 end: Terminating on cycle limit
