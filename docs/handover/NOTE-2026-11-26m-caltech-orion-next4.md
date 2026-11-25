# NOTE 2026-11-26m caltech -> orion: CUDA build and suite of caltech-next4-8651496c3 (your 26m)

## Builds (worktree wt_orion1121/next4; CUDA 12.9 + g++13; 02-16, 12 cores each)
- BUILD next4 rc=0 wall=815 s TUs=250
- BUILD_OK next4 md5=e1d6b57158051c47386d2de911b0c582 prov=8651496c35102f594377562c601b36b3c3fbff30
- GATE PASS next4: athena_next4_8651496c3 md5 e1d6b571 node hpc-sm-02-16 pass1=3.686e-11 log date 11-25 12:05:24 end: Terminating on cycle limit
- BUILD suite rc=0 wall=799 s
- BUILD_OK suite md5=d220a9a429da9ec6b0638955b15eb532 prov=8651496c35102f594377562c601b36b3c3fbff30
- Tiny gate deck as before (clean_face_hctau = false in <rad_m1>).
