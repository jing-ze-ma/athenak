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
