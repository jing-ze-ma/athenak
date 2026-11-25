# NOTE 2026-11-24zzzzd6 caltech -> orion: GPU gate of rel-a c880e2703 (your 24zzzzd6)

## 1. Builds (own worktree wt_orion1121/rela at c880e2703, CUDA 12.9 + g++ 13, HOPPER90, CPU steps on 02-03)
- BUILD rela rc=0 wall=718 s TUs=250
- BUILD_OK rela md5=5cfc3be547ad477553ec5db0e9c07fca prov=c880e2703605969324c5e1695fae47fbfb1509fc
- GATE PASS rela: athena_rela_c880e2703 md5 5cfc3be5 node hpc-sm-02-03 pass1=3.686e-11 log date 11-24 18:42:42 end: Terminating on cycle limit
- BUILD suite rc=0 wall=713 s
- BUILD_OK suite md5=51bcde3e842c725d5e0397982afda953 prov=c880e2703605969324c5e1695fae47fbfb1509fc
- star_vet build follows; the suite starts on the built_in binary.

## 1. Builds so far (11-24 18:50, own worktree wt_orion1121/rela at c880e2703, CUDA 12.9 + g++ 13, HOPPER90, CPU steps on 02-03)
- BUILD rela rc=0 wall=718 s TUs=250
- BUILD_OK rela md5=5cfc3be547ad477553ec5db0e9c07fca prov=c880e2703605969324c5e1695fae47fbfb1509fc
- GATE PASS rela: athena_rela_c880e2703 md5 5cfc3be5 node hpc-sm-02-03 pass1=3.686e-11 log date 11-24 18:42:42 end: Terminating on cycle limit
- BUILD suite rc=0 wall=713 s
- BUILD_OK suite md5=51bcde3e842c725d5e0397982afda953 prov=c880e2703605969324c5e1695fae47fbfb1509fc
- **star_vet: not built.** -D PROBLEM=star_vet fails at CMake; the branch has src/pgen/he_star_vet.cpp and no star_vet source or alias. Nothing edited. The error:
```
45:CMake Error at src/CMakeLists.txt:263 (target_sources):
46-  Cannot find source file:
47-    pgen/star_vet.cpp
48-  Tried extensions .c .C .c++ .cc .cpp .cxx .cu .mpp .m .M .mm .ixx .cppm
49-  .ccm .cxxm .c++m .h .hh .h++ .hm .hpp .hxx .in .txx .f .F .for .f77 .f90
50-  .f95 .f03 .hip .ispc
51--- Generating done (2.6s)
52-CMake Generate step failed.  Build files cannot be regenerated correctly.
-- Generating done (2.6s)
CMake Generate step failed.  Build files cannot be regenerated correctly.
```

## 3. Old-rst check, M35ZAMS (4 ranks, 01-17)
- wr8: 95813cf85 restarts the T2 stop rst (m35zams.stop1124b, c22155) with time2_gtime_slope = false and clean_stag_thinf = false (my reading of "keys false"; the other candidate keys as in production) for 8 cycles and writes its rst (c22163). base16: 95813cf85 from that rst for 16 cycles. rel8: rel-a from the same rst for 8 cycles, writing its own rst; rel8b: rel-a from that for 8 more. Compared by printed t.
- wr8: 2026-11-24 18:50:22 rc=0 (srun) ended: Terminating on cycle limit wall=35 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- base16: 2026-11-24 18:50:51 rc=0 (srun) ended: Terminating on cycle limit wall=29 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- rel8: 2026-11-24 18:51:30 rc=0 (srun) ended: Terminating on cycle limit wall=39 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- rel8b: 2026-11-24 18:51:57 rc=0 (srun) ended: Terminating on cycle limit wall=27 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- user.hst rel-a first 8: 8 rows, printed-identical to the base 16-cycle run at the same t: 8; first different or missing t: none
- user.hst rel-a second 8: 8 rows, printed-identical to the base 16-cycle run at the same t: 8; first different or missing t: none
- user.hst: base 16 rows; BITWISE (printed)
- hydro.hst rel-a first 8: 8 rows, printed-identical to the base 16-cycle run at the same t: 8; first different or missing t: none
- hydro.hst rel-a second 8: 8 rows, printed-identical to the base 16-cycle run at the same t: 8; first different or missing t: none
- hydro.hst: base 16 rows; BITWISE (printed)
