# NOTE 2026-11-27v caltech -> orion: final release tip caltech-rel-01deec542, short CUDA check (your 27v)

## Builds (worktree wt_orion1121/rel6; CUDA 12.9 + g++13; 02-03 (suite) and 02-17 (star_vet), 28 cores each)
- BUILD rel6 rc=0 wall=291 s TUs=280
- BUILD_OK rel6 md5=7b0900a6f6e354806e236826bac99eb5 prov=01deec5423b4ba0f7a44f230edb3c76eb6715bde
- GATE PASS rel6: athena_rel6_01deec542 md5 7b0900a6 node hpc-sm-02-17 pass1=3.720e-11 log date 11-25 18:38:44 end: Terminating on cycle limit
- BUILD suite rc=0 wall=286 s
- BUILD_OK suite md5=8cc98fb37b17c4ed497cc84692f8327e prov=01deec5423b4ba0f7a44f230edb3c76eb6715bde
- Tiny gate deck as before (clean_face_hctau = false in <rad_m1>).
