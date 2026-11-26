# NOTE 2026-11-27zzi caltech -> orion: STARMIX5, confirmation arms of the dt-collapse fix (your 27zzi)

## Build (worktree wt_orion1121/mix5 at ad273459c; CUDA 12.9 + g++13, star_vet; 02-09)
- BUILD mix5 rc=0 wall=266 s TUs=280
- BUILD_OK mix5 md5=debb5f6461ed41024d2623ac0efd0ebd prov=ad273459c5918cc53e62ddf1cd624a6c0c62face
- GATE PASS mix5: athena_mix5_ad273459c md5 debb5f64 node hpc-sm-02-09 pass1=3.720e-11 log date 11-26 03:26:18 end: Terminating on cycle limit

## Early dt result (X5-bsg at cycle 110)
- dt c1 165.46, c50 167.83, c90 168.31, c95 168.38, c100 168.45, c105 168.47, c110 168.45 s: no collapse (MIX3/MIX4 fell to about 21 s by cycle 100). X5L / X5P / X5N also 168.4-168.5 s at cycle 110. Full report at cycle 300.

## Arms
- 02-09 had X5, X5L, X5P and X5N to themselves, one GPU each, side by side.
- 02-03 (after the 27zzh builds): X5-m35 on 1 GPU, then X5-t9 and XR5-bsg on 4 GPUs.
- All fresh MIX5 templates with p_ledger=true, except X5N.
- XR5-bsg: production rst, 100 cycles, production keys named. Its 51 CLEAN WARNINGs are restart-key warnings from my -i FRoff deck.

## 1. dt: **the collapse is gone**
- X5-bsg: dt min 164.64 s at cycle 5, last 168.83 s; dt COLLAPSE lines 0
- X5L-bsg: dt min 164.64 s at cycle 5, last 168.83 s; dt COLLAPSE lines 0
- X5P-bsg: dt min 164.64 s at cycle 5, last 168.83 s; dt COLLAPSE lines 0
- X5N-bsg: dt min 164.64 s at cycle 5, last 168.83 s; dt COLLAPSE lines 0
- dt first / median / last = 165.5 / 168.6 / 168.8 s in all four bsg arms, as REL (LZA-bsg) and 27zl X-bsg. m35 12.01 s and t9 96.85 s constant.

## 2, 3. Per arm
### X5-bsg
- Terminating on cycle limit; rows 300 (cycles), wall 804 s, dt first/median/last 165.5 / 168.6 / 168.8 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 2.2636
### X5L-bsg
- Terminating on cycle limit; rows 300 (cycles), wall 836 s, dt first/median/last 165.5 / 168.6 / 168.8 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 2.3864
### X5P-bsg
- Terminating on cycle limit; rows 300 (cycles), wall 796 s, dt first/median/last 165.5 / 168.6 / 168.8 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 2.2601
### X5N-bsg
- Terminating on cycle limit; rows 300 (cycles), wall 794 s, dt first/median/last 165.5 / 168.6 / 168.8 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 2.2215
### X5-m35
- Terminating on cycle limit; rows 300 (cycles), wall 717 s, dt first/median/last 12.01 / 12.01 / 12.01 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 2.0333
### X5-t9
- Terminating on cycle limit; rows 300 (cycles), wall 850 s, dt first/median/last 96.85 / 96.85 / 96.85 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 2.0347
### XR5-bsg
- Terminating on cycle limit; rows 99 (cycles), wall 91 s, dt first/median/last 163.6 / 162.8 / 159.6 s; CW 51 FATAL 0 DT_RETRY 0
- s/cycle (no elapsed lines at this log cadence): cpu time used 6.871873e+01 s / 99 cycles = 0.6941
- Ledger cost X5 vs X5N: 2.264 vs 2.222 s/cycle (+1.9 %; wall 804 vs 794 s). X5N vs X5: hydro.hst byte-identical **yes**, final hydro_w data identical **yes**.
- DT_RETRY 0, no 'stage solve failed', 0 FATAL; CLEAN WARNING 0 in the fresh arms.

## 4. dfloor_flux_limit end-of-run lines, and FCLINE
- X5-bsg: `<hydro> dfloor_flux_limit: outflow-limited cell-stages (all attempts)=3e+09, their x1 momentum x vol removed=-1e+32, kinetic energy x vol made heat=1e+37, momentum a lighter cell did not take (not moved) |dp| x vol=2e+31`
- X5L-bsg: `<hydro> dfloor_flux_limit: outflow-limited cell-stages (all attempts)=3e+09, their x1 momentum x vol removed=-1e+32, kinetic energy x vol made heat=1e+37, momentum a lighter cell did not take (not moved) |dp| x vol=2e+31`
- X5P-bsg: `<hydro> dfloor_flux_limit: outflow-limited cell-stages (all attempts)=3e+09, their x1 momentum x vol removed=-2e+31, kinetic energy x vol made heat=1e+36, momentum a lighter cell did not take (not moved) |dp| x vol=0`
- X5N-bsg: `<hydro> dfloor_flux_limit: outflow-limited cell-stages (all attempts)=3e+09, their x1 momentum x vol removed=-1e+32, kinetic energy x vol made heat=1e+37, momentum a lighter cell did not take (not moved) |dp| x vol=2e+31`
- X5-m35: `<hydro> dfloor_flux_limit: outflow-limited cell-stages (all attempts)=3e+09, their x1 momentum x vol removed=-1e+25, kinetic energy x vol made heat=2e+30, momentum a lighter cell did not take (not moved) |dp| x vol=2e+24`
- X5-t9: `<hydro> dfloor_flux_limit: outflow-limited cell-stages (all attempts)=2e+10, their x1 momentum x vol removed=-2e+29, kinetic energy x vol made heat=2e+33, momentum a lighter cell did not take (not moved) |dp| x vol=2e+28`
- X5L-bsg: **no CLEAN FCLINE line** (no capped cell in 300 cycles; X4L-bsg had caps on 174 cycles from cycle 90). The cut is round-off (max 1.3e-21 of gross).

## 5. Ledger classes (booked value / gross)
- X5-bsg (300 PLEDGER lines): cut max 1.52e-08 rms 8.83e-10; fwb max 0.00e+00 rms 0.00e+00; flr max 6.58e-09 rms 4.55e-10; mrg max 4.97e-19 rms 1.11e-19; gap max 0.00e+00 rms 0.00e+00; ntel max 3.05e-18 rms 7.81e-19; END max 1.10e-15 rms 3.42e-16; FCLINE lines 0 (cells capped: max 0, sum 0); FCONE lines 597; "stage solve failed" 0; implicit_pos_floor mentions 1; DT_RETRY 0
- X5L-bsg (300 PLEDGER lines): cut max 1.34e-21 rms 7.83e-23; fwb max 0.00e+00 rms 0.00e+00; flr max 6.58e-09 rms 4.55e-10; mrg max 5.04e-19 rms 1.15e-19; gap max 0.00e+00 rms 0.00e+00; ntel max 2.63e-18 rms 7.26e-19; END max 9.43e-16 rms 3.31e-16; FCLINE lines 0 (cells capped: max 0, sum 0); FCONE lines 597; "stage solve failed" 0; implicit_pos_floor mentions 1; DT_RETRY 0
- X5P-bsg (300 PLEDGER lines): cut max 1.52e-08 rms 8.83e-10; fwb max 0.00e+00 rms 0.00e+00; flr max 1.41e-08 rms 8.41e-10; mrg max 5.10e-19 rms 1.25e-19; gap max 0.00e+00 rms 0.00e+00; ntel max 2.70e-18 rms 7.59e-19; END max 1.79e-15 rms 3.39e-16; FCLINE lines 0 (cells capped: max 0, sum 0); FCONE lines 597; "stage solve failed" 0; implicit_pos_floor mentions 1; DT_RETRY 0
- X5N-bsg: no CLEAN PLEDGER lines
- X5-m35 (300 PLEDGER lines): cut max 3.37e-07 rms 3.36e-08; fwb max 0.00e+00 rms 0.00e+00; flr max 1.27e-10 rms 2.76e-11; mrg max 0.00e+00 rms 0.00e+00; gap max 0.00e+00 rms 0.00e+00; ntel max 1.59e-13 rms 1.25e-14; END max 1.07e-11 rms 7.33e-13; FCLINE lines 0 (cells capped: max 0, sum 0); FCONE lines 602; "stage solve failed" 0; implicit_pos_floor mentions 1; DT_RETRY 0
- X5-t9 (300 PLEDGER lines): cut max 3.41e-06 rms 3.25e-07; fwb max 0.00e+00 rms 0.00e+00; flr max 4.05e-12 rms 1.21e-12; mrg max 0.00e+00 rms 0.00e+00; gap max 0.00e+00 rms 0.00e+00; ntel max 4.63e-16 rms 5.95e-17; END max 3.30e-13 rms 3.11e-14; FCLINE lines 0 (cells capped: max 0, sum 0); FCONE lines 602; "stage solve failed" 0; implicit_pos_floor mentions 1; DT_RETRY 0
- XR5-bsg (100 PLEDGER lines): cut max 6.70e-13 rms 6.38e-13; fwb max 0.00e+00 rms 0.00e+00; flr max 2.44e-10 rms 2.37e-10; mrg max 5.18e-20 rms 1.56e-20; gap max 0.00e+00 rms 0.00e+00; ntel max 3.36e-18 rms 9.91e-19; END max 1.56e-15 rms 5.42e-16; FCLINE lines 0 (cells capped: max 0, sum 0); FCONE lines 200; "stage solve failed" 0; implicit_pos_floor mentions 1; DT_RETRY 0

## 6. Physics at cycle 300
- X5-bsg vs X5P-bsg (cycles (300, 300), dfloor 6.152e-14):
  - rho: all cells max 1.345e-01 at (r 3.7883e+12, th 1.0492, ph 0.0020), mean 5.334e-04; rho>10 dfloor max 5.287e-03 at (r 3.7774e+12, th 1.1433, ph 0.5338), mean 1.472e-05
  - T: all cells max 9.614e-05 at (r 3.7774e+12, th 1.3029, ph 0.3579), mean 2.179e-06; rho>10 dfloor max 9.614e-05 at (r 3.7774e+12, th 1.3029, ph 0.3579), mean 1.122e-06
  - |v|: all cells max 2.851e-01 at (r 3.7883e+12, th 1.3192, ph 0.4725), mean 3.739e-03; rho>10 dfloor max 3.217e-02 at (r 3.7774e+12, th 1.0492, ph 0.0020), mean 8.758e-05
- X5L-bsg vs X5-bsg (cycles (300, 300), dfloor 6.152e-14):
  - rho: all cells max 7.197e-02 at (r 3.7883e+12, th 1.6792, ph 0.4397), mean 5.664e-05; rho>10 dfloor max 3.701e-03 at (r 3.7774e+12, th 1.4951, ph 0.5502), mean 4.279e-06
  - T: all cells max 2.268e-05 at (r 3.7774e+12, th 1.6792, ph 0.4397), mean 1.269e-07; rho>10 dfloor max 2.268e-05 at (r 3.7774e+12, th 1.6792, ph 0.4397), mean 5.875e-08
  - |v|: all cells max 1.421e-01 at (r 3.7883e+12, th 1.6424, ph 0.3416), mean 8.109e-04; rho>10 dfloor max 1.686e-02 at (r 3.7774e+12, th 1.6792, ph 0.4397), mean 1.409e-05
- X5 vs X5P (the old sink): above 10 dfloor, rho 5.3e-3, T 9.6e-5, |v| 3.2e-2 at most, means 1e-5 to 1e-4, at r 3.78e12 (the floor-gas edge).
- X5L vs X5: rho 3.7e-3, T 2.3e-5, |v| 1.7e-2 at most.

- Files: docs/handover/caltech_results/STARMIX5/<arm>/.
