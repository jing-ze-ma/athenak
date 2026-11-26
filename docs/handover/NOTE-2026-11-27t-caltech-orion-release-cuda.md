# NOTE 2026-11-27t caltech -> orion: CUDA suite of the release code caltech-rel-82231a337 (your 27t)

## Builds (worktree wt_orion1121/rel82; CUDA 12.9 + g++13; 02-03 (suite) and 02-17 (star_vet), 28 cores each)
- BUILD rel82 rc=0 wall=299 s TUs=280
- BUILD_OK rel82 md5=e80fcc9c65def95726b5dfd80f9b861a prov=82231a337f4e5f5a5a4c314dfaaeda8e88b2ae14
- GATE FAIL rel82: athena_rel82_82231a337 md5 e80fcc9c node hpc-sm-02-17 pass1= log date 11-25 17:05:55 end: 
- BUILD suite rc=0 wall=291 s
- BUILD_OK suite md5=ee84ed64f900c6eb485736a344307211 prov=82231a337f4e5f5a5a4c314dfaaeda8e88b2ae14
- Tiny gate deck as before (clean_face_hctau = false in <rad_m1>).
