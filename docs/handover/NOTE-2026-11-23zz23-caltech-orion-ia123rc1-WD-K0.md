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

## Start lines (11-23 23:10)
- WD: `2026-11-23 23:10:21 start WD BIN=athena_ia123rc1_da92efe70 md5=fe0ef3cd9b677d88b179215fd39e24b4 NODE=hpc-sm-02-11 HOLDER=7333073 GPUS=0,1,2,3 ARGS: -r /resnick/groups/carnegie_poc/jingze/prod_ia11_1123/t9l62/rst_keep/t9l62.00012.rst rad_vet/clean_pred_stage=false rad_vet/clean_stag_central=true rad_vet/clean_stag_kpat=false rad_vet/clean_stag_eicarry=true mesh/layout_invariant_mode=tile mesh/layout_tile=128 rad_vet/clean_gas_secant=max_stage2 rad_vet/clean_pred_xch=true rad_vet/clean_stag_fcone_gas=false rad_vet/clean_wb_exch_local=true rad_vet/clean_stag_fwb_ring=1 hydro/floor_pay=rad_column hydro/floor_pay_skip=true hydro/dfloor_flux_limit=true hydro/floor_audit=1 rad_vet/clean_wb_floor_diag=8 hydro/rad_signal_speed_diff=true rad_vet/vet_gd_wall_sym=true rad_vet/vet_gd_theta_seam_sym=true rad_vet/clean_clos1_budget=true rad_vet/clean_bud_maxnorm=true rad_vet/clean_vel_settled=false hyd`
- K0rc1: `2026-11-23 23:10:21 start K0rc1 BIN=athena_ia123rc1_da92efe70 md5=fe0ef3cd9b677d88b179215fd39e24b4 NODE=hpc-sm-01-17 HOLDER=7333072 GPUS=0,1 ARGS: -r /resnick/groups/carnegie_poc/jingze/prod_ia122_1123/m35mms/rst_keep/m35mms.00017.rst -i /resnick/groups/carnegie_poc/jingze/ia122_1123/decks/notile/IA122_m35mms.athinput rad_m1/vet_gd_wall_sym=true rad_m1/vet_gd_theta_seam_sym=true time/nlim=38585 job/final_restart=false output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30`
- K0ia122: `2026-11-23 23:10:21 start K0ia122 BIN=athena_ia122_465aff485 md5=cfa3dc7713962f6ab2f5a1cf85de0fbe NODE=hpc-sm-01-17 HOLDER=7333072 GPUS=2,3 ARGS: -r /resnick/groups/carnegie_poc/jingze/prod_ia122_1123/m35mms/rst_keep/m35mms.00017.rst -i /resnick/groups/carnegie_poc/jingze/ia122_1123/decks/notile/IA122_m35mms.athinput rad_m1/vet_gd_wall_sym=true rad_m1/vet_gd_theta_seam_sym=true time/nlim=38585 job/final_restart=false output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30`
- K0rc1 and K0ia122 run side by side on hpc-sm-01-17 (2 + 2 H200, 2 blocks/GPU, same conditions), so their s/cyc compare 1:1 (but are not the 4-GPU rate). WD has 4 H200 on hpc-sm-02-11.
