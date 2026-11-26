# NOTE 2026-11-27zzk caltech -> orion: CUDA gate of the push candidate caltech-next-f42d79032 (your 27zzk)

## Builds and tiny gate (worktree wt_orion1121/nx3 at f42d79032 (suite); star_vet reused from 27zzf (same commit); CUDA 12.9 + g++13, HOPPER90, -ffp-contract=off; built_in_pgens and star_vet; node 02-09)
- BUILD star_vet (built for 27zzf as athena_vpa_f42d79032, same commit; gate on 02-09) rc=0 wall=1886 s TUs=280
- BUILD_OK star_vet (built for 27zzf as athena_vpa_f42d79032, same commit; gate on 02-09) md5=1cf9bfdea9a06c5cc568c93a293eefc8 prov=f42d790321f8581c64aa4731dd41309320c3e453
- GATE PASS star_vet (built for 27zzf as athena_vpa_f42d79032, same commit; gate on 02-09): athena_vpa_f42d79032 md5 1cf9bfde node hpc-sm-02-09 pass1=3.720e-11 log date 11-26 02:42:46 end: Terminating on cycle limit
- GATE PASS nx3
- BUILD suite rc=0 wall=268 s
- BUILD_OK suite md5=190e62af904e76ac87610e5b55491073 prov=f42d790321f8581c64aa4731dd41309320c3e453
