# NOTE 2026-11-24zzzzd2 caltech -> orion: SDO_MOM, t9l62.00015 (your 24zzzzd2)

## Build (own worktree wt_orion1121/sdomom at 373e2cc1c, CUDA 12.9 + g++ 13, CPU step on 02-09)
- BUILD sdomom rc=0 wall=403 s TUs=250
- BUILD_OK sdomom md5=58aaa0e74c962098c1d5feeed248897b prov=373e2cc1c1c6852ccaf37befe0cd544ceaa5e8c9
- GATE PASS sdomom: athena_sdomom_373e2cc1c md5 58aaa0e7 node hpc-sm-02-09 pass1=3.686e-11 log date 11-24 18:30:47 end: Terminating on cycle limit

## SDO_MOM arm (11-24 19:44)
- 2026-11-24 19:39:01 start mom hpc-sm-01-17 ranks 4 SDO_MOM=1 (SDO_CELL unset) BIN athena_sdomom_373e2cc1c md5=58aaa0e74c962098c1d5feeed248897b ARGS: -r /resnick/groups/carnegie_poc/jingze/compton_1124 ... ARGS = Con's word for word (nlim 5350 was Con's own)
- 2026-11-24 19:43:47 rc=0 (srun) ended: Terminating on cycle limit wall=286 s CW 0 FATAL 0 inf 0 DT_RETRY 0 stage fallback 0 BE backward-Euler steps=0.000000e+00 stage retries=0.000000e+00 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 

### Identity check (every SDOMOM sum line)
- SDOMOM total lines 71, sum lines 421
- cause ekin: 350 lines, cells 9669613; max |dE| 0.000e+00 (nonzero in 0 lines); max |dLz| 0.000e+00 (nonzero in 0 lines); max |dKE - deint| 9.938e+19, max |dKE - deint| / |deint| 3.207e-03; max |deint| 1.845e+38
- cause other: 27 lines, cells 71; max |dE| 5.161e+19 (nonzero in 27 lines); max |dLz| 0.000e+00 (nonzero in 0 lines); max |dKE - deint| 4.313e+12, max |dKE - deint| / |deint| 3.192e-02; max |deint| 9.044e+35
- cause vceil: 44 lines, cells 67; max |dE| 1.032e+20 (nonzero in 16 lines); max |dLz| 0.000e+00 (nonzero in 0 lines); max |dKE - deint| 0.000e+00, max |dKE - deint| / |deint| 0.000e+00; max |deint| 6.476e+39

### History against Con (883d79d97)
- user.hst vs Con: rows 71 / 71, identical 71; first difference: none; largest rel difference 0.00e+00 (row 70, column 32)
- hydro.hst vs Con: rows 71 / 71, identical 71; first difference: none; largest rel difference 0.00e+00 (row 70, column 10)

### mscore.py (both loss measures as printed; no threshold, no interpretation)
```
cycles 5278 5348 71
== vceil tag: cells/cycle | sum|dp_r| sum|dp_th| sum dp_r | dKE | max|dLz| max|dE| (raw; mean per cycle)
  update1         0.7 | 2.066e+29 2.309e+27 +2.998e+28 | 2.028e+38 | 0.0e+00 1.0e+20  x0.5
  update2         0.2 | 2.452e+28 2.156e+26 -2.128e+28 | 1.056e+37 | 0.0e+00 1.3e+19  x1.0
  flux3r          0.0 | 2.255e+23 1.838e+23 -2.255e+23 | 7.959e+30 | 0.0e+00 0.0e+00  NOT PERSISTENT
  post2           0.0 | 1.093e+23 9.822e+22 -1.093e+23 | 7.059e+30 | 0.0e+00 0.0e+00  x1.0
  which2          0.0 | 1.125e+22 1.120e+22 -1.125e+22 | 4.517e+29 | 0.0e+00 0.0e+00  x0.5
== ekin tag: cells/cycle | sum|dp_r| sum|dp_th| sum dp_r | dKE | max|dLz| max|dE| (raw; mean per cycle)
  update1        36.2 | 1.520e+28 1.222e+27 +1.205e+28 | 1.097e+36 | 0.0e+00 0.0e+00  x0.5
  update2        66.0 | 3.817e+28 3.370e+27 +1.690e+27 | 3.787e+36 | 0.0e+00 0.0e+00  x1.0
  which2      18253.1 | 3.069e+26 4.182e+25 -2.199e+26 | 1.074e+34 | 0.0e+00 0.0e+00  x0.5
  post2       96927.6 | 4.943e+18 2.537e+16 +4.788e+18 | 6.902e+19 | 0.0e+00 0.0e+00  x1.0
  flux3r      20908.8 | 5.370e+23 1.327e+23 -2.472e+23 | 1.583e+31 | 0.0e+00 0.0e+00  NOT PERSISTENT
== other tag: cells/cycle | sum|dp_r| sum|dp_th| sum dp_r | dKE | max|dLz| max|dE| (raw; mean per cycle)
  update1         0.2 | 5.397e+26 3.037e+25 -4.064e+26 | 3.038e+34 | 0.0e+00 5.2e+19  x0.5
  post2           0.2 | 4.073e+05 5.623e+03 +4.073e+05 | 1.842e+12 | 0.0e+00 6.1e+12  x1.0
  which2          0.1 | 3.997e+24 4.514e+23 -3.309e+24 | 7.805e+31 | 0.0e+00 5.1e+19  x0.5
  flux3r          0.5 | 2.040e+21 1.935e+21 +1.568e+21 | 2.684e+29 | 0.0e+00 5.0e+19  NOT PERSISTENT
  update2         0.0 | 3.974e+19 1.552e+20 +3.974e+19 | 1.205e+28 | 0.0e+00 9.6e+16  x1.0
weighted |dp_r|+|dp_th| per cycle vceil mean 1.292e+29 max 2.768e+30 at cycle 5294; share of |p_pol| mean 3.83e-09 max 8.21e-08
weighted |dp_r|+|dp_th| per cycle ekin mean 4.992e+28 max 1.497e+30 at cycle 5319; share of |p_pol| mean 1.47e-09 max 4.41e-08
weighted |dp_r|+|dp_th| per cycle other mean 2.872e+26 max 1.031e+28 at cycle 5281; share of |p_pol| mean 8.54e-12 max 3.07e-10
floorcap |dp| printed at (5300, 6.714883e+26) (5325, 1.355899e+27) = 2.738e+25 per cycle
dfloor_flux_limit x1 momentum removed, this job: -2.751440e+29 (/71 cycles = 3.875e+27 per cycle)
floor_pay max 9.917e+32 CW 0
largest single cells [('5308', '1', '0', '2.516320e+12', '1.1836e+30'), ('5319', '1', '0', '2.448269e+12', '1.1774e+30'), ('5291', '0', '0', '2.448269e+12', '3.6505e+29')]
```
- Files: caltech_results/SDO_MOM/ (log.0.gz, hst, run.log, mscore.txt, id.txt, cmp.txt).
