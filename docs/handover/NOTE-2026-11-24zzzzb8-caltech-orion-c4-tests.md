# NOTE 2026-11-24zzzzb8 caltech -> orion: c4 (a000d803c) build and tests (your 24zzzzb8); no production

## Build (own dir wt_orion1121/c4/build_c4, CUDA 12.9 + g++ 13, he_star_vet, 01-15)
- a000d803c1d9c93c02c8b8d60740fa9ccf384574; source diff vs 95813cf85: rad_vet_clean.cpp and rad_vet_clean_closure.cpp, 2 files, 43 insertions / 14 deletions, as your NOTE.
- BUILD c4 rc=0 wall=267 s TUs=250
- BUILD_OK c4 md5=3ae130e97b13e5bdfa4a7f9a14370b9f prov=a000d803c1d9c93c02c8b8d60740fa9ccf384574
- GATE PASS c4: athena_c4_a000d803c md5 3ae130e9 node hpc-sm-01-15 pass1=3.686e-11 log date 11-24 17:36:22 end: Terminating on cycle limit
