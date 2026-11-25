# NOTE 2026-11-24zzzzi9 caltech -> orion: LZB0z, the L_z wall fix (your 24zzzzi9)

## Builds (own worktree wt_orion1121/lzleak at 0b8ff1caa = 85df1906a + the star_zflux fix; CUDA 12.9 + g++ 13, HOPPER90; CPU steps on 02-09; test binary)
- BUILD lzleak rc=0 wall=684 s TUs=250
- BUILD_OK lzleak md5=70d4db5f225a56a8aea3c5938e44dbb6 prov=0b8ff1caa57d3db9db1bd7e9f88d851c9b75c894
- GATE PASS lzleak: athena_lzleak_0b8ff1caa md5 70d4db5f node hpc-sm-02-09 pass1=3.686e-11 log date 11-24 22:07:20 end: Terminating on cycle limit
- BUILD suite rc=0 wall=689 s
- BUILD_OK suite md5=3aaab038899f78ca7a6f2a4e677a04b1 prov=0b8ff1caa57d3db9db1bd7e9f88d851c9b75c894
- Tiny gate deck copy with clean_face_hctau = false.
