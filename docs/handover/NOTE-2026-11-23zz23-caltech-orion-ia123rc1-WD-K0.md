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

## WD finished (11-23 23:11)
- 2026-11-23 23:10:21 start WD BIN=athena_ia123rc1_da92efe70 md5=fe0ef3cd9b677d88b179215fd39e24b4 NODE=hpc-sm-02-11 HOLDER=7333073 GPUS=0,1,2,3 ARGS: -r /resnick/groups/carnegie_poc/jingze/prod_ia11_1123/t9l62/rst_keep/t9l62.00012.rst rad_vet/clean_pred_stage=false rad_vet/clean_stag_central=true rad_
- Echo checks: ### NOTE: parameter 'clean_gas_secant' in block 'rad_vet' is absent from the restart's input: added from the command line (= max_stage2);### NOTE: parameter 'floor_from_ekin' in block 'hydro' is absent from the restart's input: added from the command line (= true);### NOTE: parameter 'floor_audit_worst' in block 'hydro' is absent from the restart's input: added from the command line (= true);
  never-read list: 
- CW 0, stage fallbacks 0, DT_RETRY 0, pos_floor created 0, FATAL 0
- Files: caltech_results/IA121_T9L62/arm_WD/grep.txt (`grep -E "floor_audit|DT_RETRY|CLEAN WARNING|FALLBACK"`, 6 lines), run.log, log.0.gz, hst.

## (11-23 23:12) reply to 23zz24
- (1) FATAL line of the unchanged gate deck: see the Gate section above (`theta-periodic wedge: the theta wall maps are inexact (120960 of 355320 direction maps); name rad_vet/vet_gd_theta_seam_sym = true ...`).
- (2) The gate on a COPY of the gate deck with vet_gd_theta_seam_sym = true + vet_gd_wall_sym = true in <rad_m1> gives no FATAL, so no key had to be dropped (implicit_pos_floor_solve / cs_cross_mono are not in this gate deck).
  **pass 0 9.534e-02 | pass 1 3.686e-11 | pass 2 9.994e-15**, converged after 3 passes, F/c 2.9006e+42 -1.7966e+36 -2.0678e+33 (caltech_results/H200CHK/gate_ia123rc1_da92efe70_seam_wall.txt). Identical to seam_sym alone.
- (3) pass 1 equals 3.686e-11 and pass 2 is 1e-14 <= 1e-13, so your criterion is met. WD and K0 have been running in parallel since 23:21 (start lines above).

## (11-23 23:13) **WD CRASHED in its first cycle (rc1 code error, not a physics stop)**: please fix on your side
- WD (rc1 da92efe70, md5 fe0ef3cd; WD_t9l62_restart.args + time/nlim=6480 job/final_restart=false + rst/bin dt 1e30; from rst_keep/t9l62.00012; hpc-sm-02-11, 4 ranks) started 23:10:21. It aborted at 23:11:18 on all 4 ranks (core dumped).
  It got through startup, the first floor_audit lines for cycle 6457 and the dt_retry checkpoint line, then:
```
  what():  Error: Kokkos::deep_copy with no available copy mechanism: from source view ("cons") to destination view ("cons").
There is no common execution space that can access both source's space
(Cuda) and destination's space (Host), so source and destination
must be contiguous and have the same layout.

```
- So a Kokkos::deep_copy of a Cuda view 'cons' into a Host view 'cons' with no common space. It sits right after `<time> dt_retry checkpoint: 1434 members, 406 Views ...`, so it is likely the floor_audit_worst / dt_retry path.
  The K0rc1 arm (m35mms, same binary, no floor_audit_worst / floor_from_ekin) runs normally, so the trigger is in WD's extra keys.
- The 'CW 0 ... FATAL 0' WD line above is just the counts of a 1-cycle log; the run did NOT complete. The echo checks did pass: floor_from_ekin = true, floor_audit_worst = true and clean_gas_secant = max_stage2 are all 'added from the command line'; the never-read list is empty.
- Files: caltech_results/IA121_T9L62/arm_WD/crash_context.txt (startup tail + backtrace), log.0.gz, run.log, grep.txt. I did not retry and did not patch anything. hpc-sm-02-11 is free again.

## K0 finished (11-23 23:24)
- K0rc1: 2026-11-23 23:22:59 rc=0 last: elapsed=6.100928e+02 cycle=38550 time=3.043792e+05 dt=1.200603e+01; CW 0, stage fallbacks 0, DT_RETRY 0, FATAL 0; s/cyc 3.6004 (c38400-38550)
- K0ia122: 2026-11-23 23:23:34 rc=0 last: elapsed=6.391209e+02 cycle=38550 time=3.043792e+05 dt=1.200603e+01; CW 0, stage fallbacks 0, DT_RETRY 0, FATAL 0; s/cyc 3.7778 (c38400-38550)
- user.hst rc1 vs IA12.2: BITWISE identical (200 rows)
- hydro.hst rc1 vs IA12.2: BITWISE identical (200 rows)
- Both ran side by side on hpc-sm-01-17 (2+2 H200, 2 blocks/GPU), so the s/cyc ratio is the comparison, not the absolute rate. Files: caltech_results/IA123_K0/{K0rc1,K0ia122}/.

## (11-23 23:28) rc1-wd 8197e0c4f (your 24a): build + gate + GPU check
- BUILD_OK 23:27:45 athena_ia123rc1wd_8197e0c4f md5=438f78cba817222ba1929faf88abd70f TUs=53; rc1 binary untouched: md5 fe0ef3cd
- GATE PASS (seam+wall sym deck copy): athena_ia123rc1wd_8197e0c4f md5 438f78cb node hpc-sm-02-11 pass1=3.686e-11 pass2=9.994e-15
- GPU check, 1 H200, `tests/vet_clean/run_tests.py --exe athena_ia123rc1wd_8197e0c4f --only floor_pay,inf_reject --jobs 1 --threads 1`:
```
clean VET test suite: 2 tests, full, 1 at a time x 1 threads, binary /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc1wd_8197e0c4f
FAIL floor_pay             2 s    <- ERROR run: return code 1 (see /resnick/groups/carnegie_poc/jingze/ia123_1123/wd2/vet_clean_runs/floor_pay/run/run.log)
FAIL inf_reject            2 s    <- ERROR run: return code 1 (see /resnick/groups/carnegie_poc/jingze/ia123_1123/wd2/vet_clean_runs/inf_reject/run/run.log)
0 passed, 2 failed, 4 s
srun: error: hpc-sm-02-11: task 0: Exited with exit code 1
vet_clean rc=1
```
