# NOTE 2026-11-27zp caltech -> orion: rayunify round 1 on H200 (your 27zp)

## Builds / gate not OK (stopped, as the note says)
- BUILD ru rc=2 wall=466 s TUs=221
- BUILD_FAIL ru
- BUILD suite rc=2 wall=453 s
- BUILD_FAIL
```
             ^

/resnick/groups/carnegie_poc/jingze/wt_orion1121/ru/src/rad_vet/rays/sc_cart_sweep_mb.cpp(266): error: A type local to a function ("lambda [](int, int, int, Real)->void") cannot be used in the return or parameter types of the operator() of 
      auto writer = [=] __attribute__((host)) __attribute__((device))(const int m, const int k, const int j) {
                    ^

1 error detected in the compilation of "/resnick/groups/carnegie_poc/jingze/wt_orion1121/ru/src/rad_vet/rays/sc_cart_sweep_mb.cpp".
[ 79%] Building CXX object src/CMakeFiles/athena.dir/rad_vet/rays/sc_cs_sweep.cpp.o
make[2]: *** [src/CMakeFiles/athena.dir/build.make:2781: src/CMakeFiles/athena.dir/rad_vet/rays/sc_cart_sweep_mb.cpp.o] Error 2
make[2]: *** Waiting for unfinished jobs....
make[1]: *** [CMakeFiles/Makefile2:373: src/CMakeFiles/athena.dir/all] Error 2
make: *** [Makefile:136: all] Error 2
srun: error: hpc-sm-02-11: task 0: Exited with exit code 2
--
             ^

/resnick/groups/carnegie_poc/jingze/wt_orion1121/ru/src/rad_vet/rays/sc_cart_sweep_mb.cpp(266): error: A type local to a function ("lambda [](int, int, int, Real)->void") cannot be used in the return or parameter types of the operator() of 
      auto writer = [=] __attribute__((host)) __attribute__((device))(const int m, const int k, const int j) {
                    ^

1 error detected in the compilation of "/resnick/groups/carnegie_poc/jingze/wt_orion1121/ru/src/rad_vet/rays/sc_cart_sweep_mb.cpp".
make[2]: *** [src/CMakeFiles/athena.dir/build.make:2781: src/CMakeFiles/athena.dir/rad_vet/rays/sc_cart_sweep_mb.cpp.o] Error 2
make[2]: *** Waiting for unfinished jobs....
/resnick/groups/carnegie_poc/jingze/wt_orion1121/ru/src/hydro/hydro.hpp(573): warning #177-D: variable "wb_option_" was declared but never referenced
        const auto wb_option_ = wb_option;
                   ^

Remark: The warnings can be suppressed with "-diag-suppress <warning-number>"

/resnick/groups/carnegie_poc/jingze/wt_orion1121/ru/src/hydro/hydro.hpp(575): warning #177-D: variable "igm1" was declared but never referenced
        Real igm1 = 1.0/(eos.gamma-1.0);
--
make[1]: *** [CMakeFiles/Makefile2:373: src/CMakeFiles/athena.dir/all] Error 2
make: *** [Makefile:136: all] Error 2
srun: error: hpc-sm-02-11: task 0: Exited with exit code 2
```

- The full error line (cut off above), the same in both builds (star_vet and built_in_pgens):
```
/resnick/groups/carnegie_poc/jingze/wt_orion1121/ru/src/rad_vet/rays/sc_cart_sweep_mb.cpp(266): error: A type local to a function ("lambda [](int, int, int, Real)->void") cannot be used in the return or parameter types of the operator() of an extended __device__ or __host__ __device__ lambda
      auto writer = [=] __attribute__((host)) __attribute__((device))(const int m, const int k, const int j) {
                    ^

```
- One error only. Nothing was run; the 27zp arms wait for a fixed tip.

## Builds and gate (worktree wt_orion1121/ru2 at 8375031bb (27zs fix of 97d437d68); CUDA 12.9 + g++13, HOPPER90, -ffp-contract=off; node 02-09)
- BUILD ru2 rc=0 wall=926 s TUs=280
- BUILD_OK ru2 md5=1dcc12e9ae4f8a6c8d8d9558775da1a9 prov=8375031bb0b81f20aaa25cb160d21c5af55e587d
- GATE PASS ru2: athena_ru2_8375031bb md5 1dcc12e9 node hpc-sm-02-09 pass1=3.720e-11 log date 11-26 00:00:23 end: Terminating on cycle limit
- BUILD suite rc=0 wall=931 s
- BUILD_OK suite md5=02dfb6940c9005dcb4326deba37adf9f prov=8375031bb0b81f20aaa25cb160d21c5af55e587d

## 1. Bitwise REL vs RU (64^3, hot blob, vet dumps bin_double at cycle 0 and 6; 02-11)
- B1c: 0 / 0 bin files; whole file identical 0, data after the parameter block identical 0 -> DIFFERENT
- B4c: 0 / 0 bin files; whole file identical 0, data after the parameter block identical 0 -> DIFFERENT
- B1s: 0 / 0 bin files; whole file identical 0, data after the parameter block identical 0 -> DIFFERENT
- B4s: 0 / 0 bin files; whole file identical 0, data after the parameter block identical 0 -> DIFFERENT
- B4sm: 0 / 0 bin files; whole file identical 0, data after the parameter block identical 0 -> DIFFERENT
- B1c_REL FATAL: ### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel6/src/parameter_input.cpp at line 469 Block name 'output1' on command line not found 
- B1c_RU FATAL: ### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/ru2/src/parameter_input.cpp at line 469 Block name 'output1' on command line not found 
- B4c_REL FATAL: ### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel6/src/parameter_input.cpp at line 469 Block name 'output1' on command line not found 
- B4c_RU FATAL: ### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/ru2/src/parameter_input.cpp at line 469 Block name 'output1' on command line not found 
- B1s_REL FATAL: ### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel6/src/parameter_input.cpp at line 469 Block name 'output1' on command line not found 
- B1s_RU FATAL: ### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/ru2/src/parameter_input.cpp at line 469 Block name 'output1' on command line not found 
- B4s_REL FATAL: ### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel6/src/parameter_input.cpp at line 469 Block name 'output1' on command line not found 
- B4s_RU FATAL: ### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/ru2/src/parameter_input.cpp at line 469 Block name 'output1' on command line not found 
- B4sm_REL FATAL: ### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel6/src/parameter_input.cpp at line 469 Block name 'output1' on command line not found 
- B4sm_RU FATAL: ### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/ru2/src/parameter_input.cpp at line 469 Block name 'output1' on command line not found 
    B1c_REL rc=1 wall=13 s bin athena_rel6_01deec542_suite args: -i /resnick/groups/carnegie_poc/jingze/scspeed_1127/cart2.athinput mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 problem/vet_test=hot_blob problem/blob_rho=1.0 problem/blob_temp=2.0 problem/bg_rho=0.1 problem/bg_temp=1.0 problem/blob_x1=0.5 problem/blob_x2=0.5 problem/blob_x3=0.5 problem/blob_a1=0.25 problem/blob_a2=0.25 problem/blob_a3=0.25 meshblock/nx1=64 meshblock/nx2=64 meshblock/nx3=64 output1/file_type=bin output1/variable=vet output1/bin_double=true output1/dcycle=6 time/nlim=6
    B1c_RU rc=1 wall=10 s bin athena_ru2_8375031bb_suite args: -i /resnick/groups/carnegie_poc/jingze/scspeed_1127/cart2.athinput mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 problem/vet_test=hot_blob problem/blob_rho=1.0 problem/blob_temp=2.0 problem/bg_rho=0.1 problem/bg_temp=1.0 problem/blob_x1=0.5 problem/blob_x2=0.5 problem/blob_x3=0.5 problem/blob_a1=0.25 problem/blob_a2=0.25 problem/blob_a3=0.25 meshblock/nx1=64 meshblock/nx2=64 meshblock/nx3=64 output1/file_type=bin output1/variable=vet output1/bin_double=true output1/dcycle=6 time/nlim=6
    B4c_REL rc=1 wall=13 s bin athena_rel6_01deec542_suite args: -i /resnick/groups/carnegie_poc/jingze/scspeed_1127/cart2.athinput mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 problem/vet_test=hot_blob problem/blob_rho=1.0 problem/blob_temp=2.0 problem/bg_rho=0.1 problem/bg_temp=1.0 problem/blob_x1=0.5 problem/blob_x2=0.5 problem/blob_x3=0.5 problem/blob_a1=0.25 problem/blob_a2=0.25 problem/blob_a3=0.25 meshblock/nx1=64 meshblock/nx2=32 meshblock/nx3=32 output1/file_type=bin output1/variable=vet output1/bin_double=true output1/dcycle=6 time/nlim=6
    B4c_RU rc=1 wall=8 s bin athena_ru2_8375031bb_suite args: -i /resnick/groups/carnegie_poc/jingze/scspeed_1127/cart2.athinput mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 problem/vet_test=hot_blob problem/blob_rho=1.0 problem/blob_temp=2.0 problem/bg_rho=0.1 problem/bg_temp=1.0 problem/blob_x1=0.5 problem/blob_x2=0.5 problem/blob_x3=0.5 problem/blob_a1=0.25 problem/blob_a2=0.25 problem/blob_a3=0.25 meshblock/nx1=64 meshblock/nx2=32 meshblock/nx3=32 output1/file_type=bin output1/variable=vet output1/bin_double=true output1/dcycle=6 time/nlim=6
    B1s_REL rc=1 wall=4 s bin athena_rel6_01deec542_suite args: -i /resnick/groups/carnegie_poc/jingze/scspeed_1127/sp2.athinput mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 problem/vet_test=hot_blob problem/blob_rho=1.0 problem/blob_temp=2.0 problem/bg_rho=0.1 problem/bg_temp=1.0 problem/blob_x1=2.8 problem/blob_x2=1.5707963 problem/blob_x3=0.5 problem/blob_a1=0.6 problem/blob_a2=0.25 problem/blob_a3=0.25 meshblock/nx1=64 meshblock/nx2=64 meshblock/nx3=64 output1/file_type=bin output1/variable=vet output1/bin_double=true output1/dcycle=6 time/nlim=6
    B1s_RU rc=1 wall=10 s bin athena_ru2_8375031bb_suite args: -i /resnick/groups/carnegie_poc/jingze/scspeed_1127/sp2.athinput mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 problem/vet_test=hot_blob problem/blob_rho=1.0 problem/blob_temp=2.0 problem/bg_rho=0.1 problem/bg_temp=1.0 problem/blob_x1=2.8 problem/blob_x2=1.5707963 problem/blob_x3=0.5 problem/blob_a1=0.6 problem/blob_a2=0.25 problem/blob_a3=0.25 meshblock/nx1=64 meshblock/nx2=64 meshblock/nx3=64 output1/file_type=bin output1/variable=vet output1/bin_double=true output1/dcycle=6 time/nlim=6
    B4s_REL rc=1 wall=14 s bin athena_rel6_01deec542_suite args: -i /resnick/groups/carnegie_poc/jingze/scspeed_1127/sp2.athinput mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 problem/vet_test=hot_blob problem/blob_rho=1.0 problem/blob_temp=2.0 problem/bg_rho=0.1 problem/bg_temp=1.0 problem/blob_x1=2.8 problem/blob_x2=1.5707963 problem/blob_x3=0.5 problem/blob_a1=0.6 problem/blob_a2=0.25 problem/blob_a3=0.25 meshblock/nx1=64 meshblock/nx2=32 meshblock/nx3=32 output1/file_type=bin output1/variable=vet output1/bin_double=true output1/dcycle=6 time/nlim=6
    B4s_RU rc=1 wall=8 s bin athena_ru2_8375031bb_suite args: -i /resnick/groups/carnegie_poc/jingze/scspeed_1127/sp2.athinput mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 problem/vet_test=hot_blob problem/blob_rho=1.0 problem/blob_temp=2.0 problem/bg_rho=0.1 problem/bg_temp=1.0 problem/blob_x1=2.8 problem/blob_x2=1.5707963 problem/blob_x3=0.5 problem/blob_a1=0.6 problem/blob_a2=0.25 problem/blob_a3=0.25 meshblock/nx1=64 meshblock/nx2=32 meshblock/nx3=32 output1/file_type=bin output1/variable=vet output1/bin_double=true output1/dcycle=6 time/nlim=6
    B4sm_REL rc=1 wall=8 s bin athena_rel6_01deec542_suite args: -i /resnick/groups/carnegie_poc/jingze/scspeed_1127/sp2.athinput mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 problem/vet_test=hot_blob problem/blob_rho=1.0 problem/blob_temp=2.0 problem/bg_rho=0.1 problem/bg_temp=1.0 problem/blob_x1=2.8 problem/blob_x2=1.5707963 problem/blob_x3=0.5 problem/blob_a1=0.6 problem/blob_a2=0.25 problem/blob_a3=0.25 meshblock/nx1=64 meshblock/nx2=32 meshblock/nx3=32 rad_vet/vet_gd_halo_ipc=false output1/file_type=bin output1/variable=vet output1/bin_double=true output1/dcycle=6 time/nlim=6
    B4sm_RU rc=1 wall=8 s bin athena_ru2_8375031bb_suite args: -i /resnick/groups/carnegie_poc/jingze/scspeed_1127/sp2.athinput mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 problem/vet_test=hot_blob problem/blob_rho=1.0 problem/blob_temp=2.0 problem/bg_rho=0.1 problem/bg_temp=1.0 problem/blob_x1=2.8 problem/blob_x2=1.5707963 problem/blob_x3=0.5 problem/blob_a1=0.6 problem/blob_a2=0.25 problem/blob_a3=0.25 meshblock/nx1=64 meshblock/nx2=32 meshblock/nx3=32 rad_vet/vet_gd_halo_ipc=false output1/file_type=bin output1/variable=vet output1/bin_double=true output1/dcycle=6 time/nlim=6

## VOID: the "rayunify: bitwise" post above (bed8b3e7e)
- That post is void, and its DIFFERENT verdicts are wrong.
  - All 10 runs, REL and RU alike, stopped at start-up with `Block name 'output1' on command line not found`. cart2/sp2 have no <output1> block, and the command line cannot create one.
  - So no bin file was written. The "0 / 0 files -> DIFFERENT" lines are an artefact of my comparison script.
- The builds and the gate are valid: ru2 8375031bb, star 1dcc12e9 / suite 02dfb694, pass1 3.720e-11.
- Rerun now with copies of cart2/sp2 plus an <output1> block (file_type bin, variable vet, bin_double true, dcycle 6), nothing else changed (rayunify_1127/*_out1.athinput). The real bitwise result follows.

## 1. Bitwise REL vs RU (64^3, hot blob, vet dumps bin_double at cycle 0 and 6 from an <output1> block added to copies of cart2/sp2; 02-09, packed with non-timing arms)
- B1c: 3 / 3 bin files; whole file identical 3, data after the parameter block identical 3 -> IDENTICAL
- B4c: 3 / 3 bin files; whole file identical 3, data after the parameter block identical 3 -> IDENTICAL
- B1s: 3 / 3 bin files; whole file identical 3, data after the parameter block identical 3 -> IDENTICAL
- B4s: 3 / 3 bin files; whole file identical 3, data after the parameter block identical 3 -> IDENTICAL
- B4sm: 3 / 3 bin files; whole file identical 3, data after the parameter block identical 3 -> IDENTICAL
    B1c_REL rc=0 wall=3 s bin athena_rel6_01deec542_suite args: -i /resnick/groups/carnegie_poc/jingze/rayunify_1127/cart2_out1.athinput mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 problem/vet_test=hot_blob problem/blob_rho=1.0 problem/blob_temp=2.0 problem/bg_rho=0.1 problem/bg_temp=1.0 problem/blob_x1=0.5 problem/blob_x2=0.5 problem/blob_x3=0.5 problem/blob_a1=0.25 problem/blob_a2=0.25 problem/blob_a3=0.25 meshblock/nx1=64 meshblock/nx2=64 meshblock/nx3=64 time/nlim=6
    B1c_RU rc=0 wall=4 s bin athena_ru2_8375031bb_suite args: -i /resnick/groups/carnegie_poc/jingze/rayunify_1127/cart2_out1.athinput mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 problem/vet_test=hot_blob problem/blob_rho=1.0 problem/blob_temp=2.0 problem/bg_rho=0.1 problem/bg_temp=1.0 problem/blob_x1=0.5 problem/blob_x2=0.5 problem/blob_x3=0.5 problem/blob_a1=0.25 problem/blob_a2=0.25 problem/blob_a3=0.25 meshblock/nx1=64 meshblock/nx2=64 meshblock/nx3=64 time/nlim=6
    B4c_REL rc=0 wall=12 s bin athena_rel6_01deec542_suite args: -i /resnick/groups/carnegie_poc/jingze/rayunify_1127/cart2_out1.athinput mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 problem/vet_test=hot_blob problem/blob_rho=1.0 problem/blob_temp=2.0 problem/bg_rho=0.1 problem/bg_temp=1.0 problem/blob_x1=0.5 problem/blob_x2=0.5 problem/blob_x3=0.5 problem/blob_a1=0.25 problem/blob_a2=0.25 problem/blob_a3=0.25 meshblock/nx1=64 meshblock/nx2=32 meshblock/nx3=32 time/nlim=6
    B4c_RU rc=0 wall=6 s bin athena_ru2_8375031bb_suite args: -i /resnick/groups/carnegie_poc/jingze/rayunify_1127/cart2_out1.athinput mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 problem/vet_test=hot_blob problem/blob_rho=1.0 problem/blob_temp=2.0 problem/bg_rho=0.1 problem/bg_temp=1.0 problem/blob_x1=0.5 problem/blob_x2=0.5 problem/blob_x3=0.5 problem/blob_a1=0.25 problem/blob_a2=0.25 problem/blob_a3=0.25 meshblock/nx1=64 meshblock/nx2=32 meshblock/nx3=32 time/nlim=6
    B1s_REL rc=0 wall=4 s bin athena_rel6_01deec542_suite args: -i /resnick/groups/carnegie_poc/jingze/rayunify_1127/sp2_out1.athinput mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 problem/vet_test=hot_blob problem/blob_rho=1.0 problem/blob_temp=2.0 problem/bg_rho=0.1 problem/bg_temp=1.0 problem/blob_x1=2.8 problem/blob_x2=1.5707963 problem/blob_x3=0.5 problem/blob_a1=0.6 problem/blob_a2=0.25 problem/blob_a3=0.25 meshblock/nx1=64 meshblock/nx2=64 meshblock/nx3=64 time/nlim=6
    B1s_RU rc=0 wall=8 s bin athena_ru2_8375031bb_suite args: -i /resnick/groups/carnegie_poc/jingze/rayunify_1127/sp2_out1.athinput mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 problem/vet_test=hot_blob problem/blob_rho=1.0 problem/blob_temp=2.0 problem/bg_rho=0.1 problem/bg_temp=1.0 problem/blob_x1=2.8 problem/blob_x2=1.5707963 problem/blob_x3=0.5 problem/blob_a1=0.6 problem/blob_a2=0.25 problem/blob_a3=0.25 meshblock/nx1=64 meshblock/nx2=64 meshblock/nx3=64 time/nlim=6
    B4s_REL rc=0 wall=12 s bin athena_rel6_01deec542_suite args: -i /resnick/groups/carnegie_poc/jingze/rayunify_1127/sp2_out1.athinput mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 problem/vet_test=hot_blob problem/blob_rho=1.0 problem/blob_temp=2.0 problem/bg_rho=0.1 problem/bg_temp=1.0 problem/blob_x1=2.8 problem/blob_x2=1.5707963 problem/blob_x3=0.5 problem/blob_a1=0.6 problem/blob_a2=0.25 problem/blob_a3=0.25 meshblock/nx1=64 meshblock/nx2=32 meshblock/nx3=32 time/nlim=6
    B4s_RU rc=0 wall=7 s bin athena_ru2_8375031bb_suite args: -i /resnick/groups/carnegie_poc/jingze/rayunify_1127/sp2_out1.athinput mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 problem/vet_test=hot_blob problem/blob_rho=1.0 problem/blob_temp=2.0 problem/bg_rho=0.1 problem/bg_temp=1.0 problem/blob_x1=2.8 problem/blob_x2=1.5707963 problem/blob_x3=0.5 problem/blob_a1=0.6 problem/blob_a2=0.25 problem/blob_a3=0.25 meshblock/nx1=64 meshblock/nx2=32 meshblock/nx3=32 time/nlim=6
    B4sm_REL rc=0 wall=9 s bin athena_rel6_01deec542_suite args: -i /resnick/groups/carnegie_poc/jingze/rayunify_1127/sp2_out1.athinput mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 problem/vet_test=hot_blob problem/blob_rho=1.0 problem/blob_temp=2.0 problem/bg_rho=0.1 problem/bg_temp=1.0 problem/blob_x1=2.8 problem/blob_x2=1.5707963 problem/blob_x3=0.5 problem/blob_a1=0.6 problem/blob_a2=0.25 problem/blob_a3=0.25 meshblock/nx1=64 meshblock/nx2=32 meshblock/nx3=32 rad_vet/vet_gd_halo_ipc=false time/nlim=6
    B4sm_RU rc=0 wall=7 s bin athena_ru2_8375031bb_suite args: -i /resnick/groups/carnegie_poc/jingze/rayunify_1127/sp2_out1.athinput mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 problem/vet_test=hot_blob problem/blob_rho=1.0 problem/blob_temp=2.0 problem/bg_rho=0.1 problem/bg_temp=1.0 problem/blob_x1=2.8 problem/blob_x2=1.5707963 problem/blob_x3=0.5 problem/blob_a1=0.6 problem/blob_a2=0.25 problem/blob_a3=0.25 meshblock/nx1=64 meshblock/nx2=32 meshblock/nx3=32 rad_vet/vet_gd_halo_ipc=false time/nlim=6

## 2-3. Speed and production restart (02-11 alone; the 8-GPU arms on 02-11 + 02-17; REL/RU order A B A B)
| arm | REL run 1 / run 2 | RU run 1 / run 2 | RU / REL (means) |
|---|---|---|---|
| H1c | 0.3066 / 0.3069 | 0.2846 / 0.2835 | 0.926 |
| H1s | 0.2864 / 0.2863 | 0.2867 / 0.2870 | 1.002 |
| H4c | 0.2160 / 0.2169 | 0.1376 / 0.1365 | 0.633 |
| H4s | 0.1377 / 0.1392 | 0.1378 / 0.1389 | 0.999 |
| H4m | 0.1447 / 0.1423 | 0.1408 / 0.1406 | 0.981 |
| H4Lc | 0.7064 / 0.6980 | 0.6118 / 0.6135 | 0.872 |
| H4Ls | 0.6728 / 0.6583 | 0.6630 / 0.6672 | 0.999 |
| H8Lc | 0.5573 / 0.5482 | 0.4262 / 0.4174 | 0.763 |
| H8Ls | 0.6471 / 0.6053 | 0.5470 / 0.5586 | 0.883 |

### Timer lines (run 1 of each)
- H1c REL: FATAL 0, CW 0; `<rad_vet> vet_sc: formal solutions=7.000000e+00 SC seconds=8.879761e-01 (1.268537e-01 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=8.879761e+299`
- H1c RU: FATAL 0, CW 0; `<rad_vet> vet_sc: formal solutions=7.000000e+00 SC seconds=7.386554e-01 (1.055222e-01 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=7.386554e+299`
- H1s REL: FATAL 0, CW 0; `<rad_vet> vet_gd: rank 0 fenced seconds: source 6.603040e-03, sweeps 4.415423e-01 (of which per-shell halo+walls 1.585400e-04, of which MPI wait 0.000000e+00), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)` | `<rad_vet> vet_gd halo (rank 0): exchanges 0.000000e+00, bytes sent 0.000000e+00 (0.000000e+00 per exchange), post 0.000000e+00 s, partners 0; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- H1s RU: FATAL 0, CW 0; `<rad_vet> vet_gd: rank 0 fenced seconds: source 1.182758e-02, sweeps 4.665976e-01 (of which per-shell halo+walls 1.847000e-04, of which MPI wait 0.000000e+00), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)` | `<rad_vet> vet_gd halo (rank 0): exchanges 0.000000e+00, bytes sent 0.000000e+00 (0.000000e+00 per exchange), post 0.000000e+00 s, partners 0; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- H4c REL: FATAL 0, CW 0; `<rad_vet> vet_sc: formal solutions=7.000000e+00 SC seconds=8.226337e-01 (1.175191e-01 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=8.226337e+299`
- H4c RU: FATAL 0, CW 0; `<rad_vet> vet_sc: formal solutions=7.000000e+00 SC seconds=2.842150e-01 (4.060214e-02 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=2.842150e+299`
- H4s REL: FATAL 0, CW 0; `<rad_vet> vet_gd: rank 0 fenced seconds: source 7.218820e-02, sweeps 3.262788e-01 (of which per-shell halo+walls 1.984775e-02, of which MPI wait 6.088966e-02), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)` | `<rad_vet> vet_gd halo (rank 0): exchanges 1.813000e+03, bytes sent 6.990299e+08 (3.855653e+05 per exchange), post 1.592800e-04 s, partners 3; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- H4s RU: FATAL 0, CW 0; `<rad_vet> vet_gd: rank 0 fenced seconds: source 7.859976e-02, sweeps 3.750765e-01 (of which per-shell halo+walls 7.255627e-03, of which MPI wait 5.715836e-02), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)` | `<rad_vet> vet_gd halo (rank 0): exchanges 1.813000e+03, bytes sent 6.990299e+08 (3.855653e+05 per exchange), post 3.180380e-04 s, partners 3; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- H4m REL: FATAL 0, CW 0; `<rad_vet> vet_gd: rank 0 fenced seconds: source 7.432265e-02, sweeps 4.422493e-01 (of which per-shell halo+walls 2.793288e-01, of which MPI wait 2.149081e-01), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)` | `<rad_vet> vet_gd halo (rank 0): exchanges 1.813000e+03, bytes sent 6.990299e+08 (3.855653e+05 per exchange), post 1.425242e-02 s, partners 3; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- H4m RU: FATAL 0, CW 0; `<rad_vet> vet_gd: rank 0 fenced seconds: source 7.291705e-02, sweeps 5.251696e-01 (of which per-shell halo+walls 3.141498e-01, of which MPI wait 2.796857e-01), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)` | `<rad_vet> vet_gd halo (rank 0): exchanges 1.813000e+03, bytes sent 6.990299e+08 (3.855653e+05 per exchange), post 1.495764e-02 s, partners 3; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- H4Lc REL: FATAL 0, CW 0; `<rad_vet> vet_sc: formal solutions=7.000000e+00 SC seconds=1.844388e+00 (2.634840e-01 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=1.844388e+300`
- H4Lc RU: FATAL 0, CW 0; `<rad_vet> vet_sc: formal solutions=7.000000e+00 SC seconds=1.234357e+00 (1.763367e-01 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=1.234357e+300`
- H4Ls REL: FATAL 0, CW 0; `<rad_vet> vet_gd: rank 0 fenced seconds: source 2.232548e-01, sweeps 1.337808e+00 (of which per-shell halo+walls 1.788804e-01, of which MPI wait 1.642368e-01), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)` | `<rad_vet> vet_gd halo (rank 0): exchanges 3.619000e+03, bytes sent 3.168571e+09 (8.755377e+05 per exchange), post 4.644980e-04 s, partners 3; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- H4Ls RU: FATAL 0, CW 0; `<rad_vet> vet_gd: rank 0 fenced seconds: source 2.061266e-01, sweeps 1.327876e+00 (of which per-shell halo+walls 1.352039e-02, of which MPI wait 1.559929e-01), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)` | `<rad_vet> vet_gd halo (rank 0): exchanges 3.619000e+03, bytes sent 3.168571e+09 (8.755377e+05 per exchange), post 4.855490e-04 s, partners 3; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- H8Lc REL: FATAL 0, CW 0; `<rad_vet> vet_sc: formal solutions=7.000000e+00 SC seconds=1.834918e+00 (2.621311e-01 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=1.834918e+300`
- H8Lc RU: FATAL 0, CW 0; `<rad_vet> vet_sc: formal solutions=7.000000e+00 SC seconds=9.220981e-01 (1.317283e-01 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=9.220981e+299`
- H8Ls REL: FATAL 0, CW 0; `<rad_vet> vet_gd: rank 0 fenced seconds: source 1.912859e-01, sweeps 4.311391e+00 (of which per-shell halo+walls 4.285531e+00, of which MPI wait 1.584512e+00), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)` | `<rad_vet> vet_gd halo (rank 0): exchanges 3.619000e+03, bytes sent 2.287342e+10 (6.320370e+06 per exchange), post 0.000000e+00 s, partners 5; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- H8Ls RU: FATAL 0, CW 0; `<rad_vet> vet_gd: rank 0 fenced seconds: source 2.563097e-01, sweeps 3.836933e+00 (of which per-shell halo+walls 3.806099e+00, of which MPI wait 1.562247e+00), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)` | `<rad_vet> vet_gd halo (rank 0): exchanges 3.619000e+03, bytes sent 2.287342e+10 (6.320370e+06 per exchange), post 0.000000e+00 s, partners 5; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`

## 3. bsg production restart P4 (4 GPUs, cycle 29657 + 20), order REL RU REL RU; s/cycle over cycles 5-20
- P4_REL1: cpu time used 20.452 s / 20 cycles = 1.0226 s/cycle (no per-cycle elapsed lines: the deck has ndiag 50 and my time/ncycle_out key does not exist); clean stage steps 15.774 s; FATAL 0 CW 0; `<rad_vet> vet_gd: rank 0 fenced seconds: source 0.528824, sweeps 1.60185 (of which per-shell halo+walls 0.0552884, of which MPI wait 0.475601), tables (rotation) 0; in-solve rebuilds 0 (0 s)`
- P4_RU1: cpu time used 21.477 s / 20 cycles = 1.0739 s/cycle (no per-cycle elapsed lines: the deck has ndiag 50 and my time/ncycle_out key does not exist); clean stage steps 16.000 s; FATAL 0 CW 0; `<rad_vet> vet_gd: rank 0 fenced seconds: source 0.542678, sweeps 1.63607 (of which per-shell halo+walls 0.020381, of which MPI wait 0.471434), tables (rotation) 0; in-solve rebuilds 0 (0 s)`
- P4_REL2: cpu time used 21.212 s / 20 cycles = 1.0606 s/cycle (no per-cycle elapsed lines: the deck has ndiag 50 and my time/ncycle_out key does not exist); clean stage steps 15.797 s; FATAL 0 CW 0; `<rad_vet> vet_gd: rank 0 fenced seconds: source 0.520717, sweeps 1.62246 (of which per-shell halo+walls 0.055806, of which MPI wait 0.46749), tables (rotation) 0; in-solve rebuilds 0 (0 s)`
- P4_RU2: cpu time used 20.736 s / 20 cycles = 1.0368 s/cycle (no per-cycle elapsed lines: the deck has ndiag 50 and my time/ncycle_out key does not exist); clean stage steps 15.689 s; FATAL 0 CW 0; `<rad_vet> vet_gd: rank 0 fenced seconds: source 0.49699, sweeps 1.6137 (of which per-shell halo+walls 0.0198146, of which MPI wait 0.444285), tables (rotation) 0; in-solve rebuilds 0 (0 s)`
- kp_P4_REL: cpu time used 23.681 s / 20 cycles = 1.1840 s/cycle (no per-cycle elapsed lines: the deck has ndiag 50 and my time/ncycle_out key does not exist); clean stage steps 17.209 s; FATAL 0 CW 0; `<rad_vet> vet_gd: rank 0 fenced seconds: source 0.517473, sweeps 1.8546 (of which per-shell halo+walls 0.533013, of which MPI wait 0.463575), tables (rotation) 0; in-solve rebuilds 0 (0 s)`
- kp_P4_RU: cpu time used 21.789 s / 20 cycles = 1.0895 s/cycle (no per-cycle elapsed lines: the deck has ndiag 50 and my time/ncycle_out key does not exist); clean stage steps 17.148 s; FATAL 0 CW 0; `<rad_vet> vet_gd: rank 0 fenced seconds: source 0.526021, sweeps 1.91383 (of which per-shell halo+walls 0.178828, of which MPI wait 0.471659), tables (rotation) 0; in-solve rebuilds 0 (0 s)`
- P4 pair 1 bitwise: final hydro_w (bsg3d.hydro_w.00116.bin) identical yes; user.hst byte-identical yes
- P4 pair 2 bitwise: final hydro_w (bsg3d.hydro_w.00116.bin) identical yes; user.hst byte-identical yes
- P8: not possible. The rst holds 4 MeshBlocks of 256 x 128 x 128, and a restart cannot change the MeshBlock size (27zm test FATAL: `meshblock/nx2: restart file 128 -> now 64`). The 8-GPU numbers are H8L above.

### Summary
- Cartesian path, RU against REL:
  - 1 GPU: 0.93x;
  - 4 GPUs, 128^3, banded: 0.63x (SC per call 0.118 -> 0.041 s);
  - 4 GPUs, 256^3: 0.87x;
  - 8 GPUs, 256^3: 0.76x.
- sp path:
  - 1-4 GPUs: unchanged (0.98-1.00x).
  - per-shell halo+walls, rank 0, fenced: H4 0.0198 -> 0.0073 s; H4L 0.179 -> 0.0135 s. The IPC fused exchange fires now.
  - The MPI path (H4m, ipc false) does not gain: halo+walls 0.279 -> 0.314 s.
  - 8 GPUs (2 nodes, MPI between nodes): 0.88x. halo+walls there is 4.29 -> 3.81 s, and it still dominates the sweep (4.31 / 3.84 s).
- bsg P4 (production layout, 4 GPUs, 20 cycles):
  - s/cycle REL 1.023 / 1.061, RU 1.074 / 1.037: equal within the 4 % run-to-run spread.
  - The ray share is small: fenced sweeps 1.60 vs 1.62 s per 20 cycles. halo+walls 0.055 -> 0.020 s.
  - kp (fenced, one extra arm each): REL vet_vgd_shell 1.072 s (5000 launches); RU vet_vgd_shell_int 1.010 s (4960 launches) + the edge items.
  - **Bitwise P4: yes**: the final hydro_w and user.hst of REL and RU are byte-identical in both pairs.
- P8: not possible from the rst (4 MeshBlocks; a restart cannot change the MeshBlock size). The 8-GPU numbers are H8L.
- My error: the P4 arms used time/ncycle_out=1, which is not a key, so no per-cycle elapsed lines were written. Their s/cycle is cpu time / 20 cycles, start-up included.
- Logs (first 200 / last 100 lines) and kp.0 are in docs/handover/caltech_results/RAYUNIFY/<arm>/. No rst/bin is kept.
