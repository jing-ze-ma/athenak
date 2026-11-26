# NOTE 2026-11-27zzp caltech -> orion: post-push round rerun on caltech-next-41bbe91d6 (your 27zzp, UPDATE)

## Builds and tiny gate (worktree wt_orion1121/nx5 at 41bbe91d6; CUDA 12.9 + g++13, HOPPER90, -ffp-contract=off; built_in_pgens and star_vet; node 02-09)
- BUILD nx5 rc=0 wall=444 s TUs=280
- BUILD_OK nx5 md5=33d1a7a992a8054cb40af42038e7dd9e prov=41bbe91d645ecb648281bbe4f08ec0c8552c0923
- GATE PASS nx5: athena_nx5_41bbe91d6 md5 33d1a7a9 node hpc-sm-02-17 pass1=3.720e-11 log date 11-26 14:23:07 end: Terminating on cycle limit
- BUILD suite rc=0 wall=395 s
- BUILD_OK suite md5=fbd69dfc4c549406861c2e27cfd8055b prov=41bbe91d645ecb648281bbe4f08ec0c8552c0923
