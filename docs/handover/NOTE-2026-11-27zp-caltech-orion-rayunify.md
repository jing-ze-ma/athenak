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
