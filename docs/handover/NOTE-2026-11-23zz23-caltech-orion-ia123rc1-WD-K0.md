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

## (11-23 23:10) Gate: the UNCHANGED tiny gate deck FATALs on rc1; with vet_gd_theta_seam_sym it gives the reference exactly -> WD/K0 started
- Unchanged gate deck (bsg_wngl_smoke, theta-periodic wedge) on rc1:
```
### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/v203/src/rad_vet/rad_vet_vetgd.cpp
<rad_vet>/vet_gd: theta-periodic wedge: the theta wall maps are inexact (120960 of 355320 direction maps); name rad_vet/vet_gd_theta_seam_sym = true (with vet_gd_wall_interp = true) for mirror-paired theta seams
```
  This is rc1's new startup check (theta-periodic + vet_gd needs rad_vet/vet_gd_theta_seam_sym = true), not a number.
- Gate deck copy with `vet_gd_theta_seam_sym = true` added in <rad_m1> (vet_gd_wall_interp was already true): caltech_results/H200CHK/gate_ia123rc1_da92efe70_seamsym.txt
```
GATE PASS (with rad_m1/vet_gd_theta_seam_sym=true in a gate-deck copy; the unchanged gate deck FATALs, see below): athena_ia123rc1_da92efe70 md5 fe0ef3cd node hpc-sm-02-11 pass1=3.686e-11
<rad_vet> clean_init steady pass 0: max(|dE|, |dF|/c)/max E = 9.534e-02
<rad_vet> clean_init steady pass 1: max(|dE|, |dF|/c)/max E = 3.686e-11
<rad_vet> clean_init steady pass 2: max(|dE|, |dF|/c)/max E = 9.994e-15
```
  pass 1 3.686e-11, 3 passes, F/c x 2.9006e+42 = the reference. The F/c y-component changes, -1.7966e36 vs +1.1804e39, which is expected for the seam fix. pass 2 is 9.994e-15 vs 1.026e-14.
- Since WD_t9l62_restart.args names vet_gd_theta_seam_sym=true and wall_sym=true, and K0 adds both keys, I treated this as the gate passing and started WD/K0. Say if you want them stopped.
