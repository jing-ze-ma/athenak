# NOTE 2026-11-24zzzzh8 caltech -> orion: LZB0f, the FOFC-flag fix (your 24zzzzh8)

## Builds (own worktree wt_orion1121/fofcflag at 86f90f362 = aaffc366d + the fix, descends from rel-e 08c48ac10; CUDA 12.9 + g++ 13, HOPPER90; CPU steps on 02-09; test binary)
- BUILD fofcflag rc=0 wall=695 s TUs=250
- BUILD_OK fofcflag md5=6a831ffe6486d4c961a91cf1520b90b3 prov=86f90f3628141be86835f3b10ebd49ab7506aa38
- GATE PASS fofcflag: athena_fofcflag_86f90f362 md5 6a831ffe node hpc-sm-02-09 pass1=3.686e-11 log date 11-24 21:40:53 end: Terminating on cycle limit
- BUILD suite rc=0 wall=687 s
- BUILD_OK suite md5=e4717efc2f0ca1a5ef6b333e2da8f1bf prov=86f90f3628141be86835f3b10ebd49ab7506aa38
- Tiny gate deck copy with clean_face_hctau = false, as for rel-b and rel-e.
