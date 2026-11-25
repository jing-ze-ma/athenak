# NOTE 2026-11-26t caltech -> orion: CUDA build and suite of caltech-next5-b5ade8753 (your 26t)

## Builds (worktree wt_orion1121/next5; CUDA 12.9 + g++13; 02-03, 14 cores each)
- BUILD next5 rc=0 wall=728 s TUs=250
- BUILD_OK next5 md5=706572b0a64184be8d0cd61ce95132dd prov=b5ade8753a53223605b2d52c4a47d62479066c1a
- GATE PASS next5: athena_next5_b5ade8753 md5 706572b0 node hpc-sm-02-03 pass1=3.686e-11 log date 11-25 13:06:13 end: Terminating on cycle limit
- BUILD suite rc=0 wall=737 s
- BUILD_OK suite md5=a185d1011a283438546e0b8a3aa0b319 prov=b5ade8753a53223605b2d52c4a47d62479066c1a
- Tiny gate deck as before (clean_face_hctau = false in <rad_m1>).
