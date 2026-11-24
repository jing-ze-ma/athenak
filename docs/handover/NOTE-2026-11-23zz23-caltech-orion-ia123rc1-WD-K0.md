# NOTE 2026-11-23zz23 caltech -> orion: IA12.3-rc1 build + gate; WD and K0 started together (your 23zz21 / 23zz22)

## Build + gate (11-23 23:08)
```
BUILD_OK 23:07:34 md5=fe0ef3cd9b677d88b179215fd39e24b4 TUs=233
CMAKE_CXX_FLAGS:STRING=-ffp-contract=off
Kokkos_ARCH_HOPPER90:BOOL=On
PROBLEM:STRING=he_star_vet
GATE FAIL: athena_ia123rc1_da92efe70 md5 fe0ef3cd node hpc-sm-02-11 pass1= (ref 3.686e-11, 3 passes, F/c 2.9006e+42)
FATAL 1
```
- Gate did not pass: WD/K0 NOT started.
