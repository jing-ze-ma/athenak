# NOTE 2026-10-17 caltech -> raven: picstop-1017 experiment 2 on the real agcE1r restart (answer to TASK 3269c9d1)

GPU test only: agcE1r was not touched (the arms ran from a read-only copy of its agcar3d.00302.rst, md5 11340270).
- Binary: athena_gpu_he_star_m1_4775edf5_nofma, **md5 b375ccf1f8dba5304579bcf240678609**. Built with
  fin51_1015/buildq.sh (the same incremental tree and options as b16b66db_nofma; build job 4846066, 13.6 min).
  Path: builds_fin51/bin/. Not for production.
- Arms job 4846190 on **1 node x 2 H200** (hpc-sm-02-10; 2 ranks, 1 per GPU, the same launch as the agcE1r links); total 1:52.
- t_rst = 148042909.0909, **DT = 7e5**, tlim = 148742909.0909 (output2/3 dt = DT with last_time = t_rst; output5 off; ndiag 1).
  OFF took 623 s. The other arms ran with -t 00:20:46 (2x OFF).
- **b, b_nl and ab hit the -t limit before tlim** (t_end 1.48423e8 / 1.48425e8 / 1.48373e8, i.e. 46-55 % of DT). Their last
  bins are at those times, so cmpbin did not compare them with ref (time mismatch, shown below). All other arms reached tlim.
- cmpbin band: --ib 380:480 (the TASK default; photosphere index not checked separately).
- FATAL 0, NON-CONVERGED=0 (summary), DIVERGED/NaN 0, rc 0 in all 7 arms.

## stat.py
```
arm              rc    cyc  solves  pass   NC retry tail>1e-4 tail>1e-2  max_acc  checks   fails  exh(NC) ce_refus  wall_s t_end
off               0 207727    1970  5.29    0     0   12.99%   12.44%  4.1e+03       -       -        -      -/-     623 1.487429e+08
aN                0 207769    2055  6.64    0     6    6.28%    4.82%  2.2e+04 0.000000e+00 0.000000e+00 0.000000e+00 1.019000e+03/0.000000e+00     748 1.487429e+08
b                 0 208074    2661  5.59    0   101    3.47%    3.33%  1.4e+03 3.163000e+03 5.020000e+02 0.000000e+00 0.000000e+00/0.000000e+00    1258 1.484231e+08
ab                0 208017    2549  5.53    0   102    3.74%    3.29%  3.4e+03 2.551000e+03 2.000000e+00 0.000000e+00 7.130000e+02/0.000000e+00    1260 1.483732e+08
ref               0 207839    2194  9.09    0    17    5.88%    4.33%  8.1e+03       -       -        -      -/-     975 1.487429e+08
off_nl            0 207727    1970  5.29    0     0        -        -  0.0e+00       -       -        -      -/-     584 1.487429e+08
b_nl              0 208080    2673  5.59    0   102        -        -  0.0e+00 3.184000e+03 5.110000e+02 0.000000e+00 0.000000e+00/0.000000e+00    1260 1.484251e+08
```

## cmpbin.py (ref vs off aN b ab, --ib 380:480)
```
ref /resnick/groups/carnegie_poc/jingze/picstop_1017/ref t 148743000.0 148743000.0 | band i (380, 480)
off          dens max 7.1e+01 mean 2.3e-03 | band max 7.1e+01 mean 2.3e-03 ; eint max 2.0e+01 mean 1.1e-02 | band max 2.0e+01 mean 1.1e-02 ; m1_e max 4.5e+00 mean 6.8e-03 | band max 4.5e+00 mean 6.8e-03 ; m1_f1 max 8.3e-01 mean 4.8e-03 | band max 8.3e-01 mean 4.8e-03
aN           dens max 2.9e+02 mean 2.6e-03 | band max 2.9e+02 mean 2.6e-03 ; eint max 9.1e+01 mean 1.3e-02 | band max 9.1e+01 mean 1.3e-02 ; m1_e max 4.2e+00 mean 8.4e-03 | band max 4.2e+00 mean 8.4e-03 ; m1_f1 max 7.6e-01 mean 6.0e-03 | band max 7.6e-01 mean 6.0e-03
b            hydro_w: t 1.484230e+08 != ref 1.487430e+08 ; m1: t 1.484230e+08 != ref 1.487430e+08
ab           hydro_w: t 1.483730e+08 != ref 1.487430e+08 ; m1: t 1.483730e+08 != ref 1.487430e+08
```

## Summary lines per arm (`implicit transport:`, `<rad_m1> picstop:`)
```
## off
<rad_m1> implicit transport: solves=1.970000e+03 Picard iterations mean=5.287310e+00 max=1.400000e+01 NON-CONVERGED=0.000000e+00
## aN
<rad_m1> implicit transport: solves=2.055000e+03 Picard iterations mean=6.636010e+00 max=4.600000e+01 NON-CONVERGED=0.000000e+00
<rad_m1> picstop: implicit_final_resid_stop=0.000000e+00 (rel 0.000000e+00 x tol_eff) checks=0.000000e+00 failed+continued=0.000000e+00 failed-no-pass-left(NC)=0.000000e+00 solves accepted after a failed check=0.000000e+00 passes after a failed check=0.000000e+00 | implicit_conv_est_guard=0.000000e+00 refused=0.000000e+00 | implicit_conv_est_nofloor=1 refused=1.019000e+03 extra stops=0.000000e+00 
## b
<rad_m1> implicit transport: solves=2.661000e+03 Picard iterations mean=5.586246e+00 max=3.500000e+01 NON-CONVERGED=0.000000e+00
<rad_m1> picstop: implicit_final_resid_stop=1.000000e-04 (rel 0.000000e+00 x tol_eff) checks=3.163000e+03 failed+continued=5.020000e+02 failed-no-pass-left(NC)=0.000000e+00 solves accepted after a failed check=1.040000e+02 passes after a failed check=1.760000e+03 | implicit_conv_est_guard=0.000000e+00 refused=0.000000e+00 | implicit_conv_est_nofloor=0 refused=0.000000e+00 extra stops=0.000000e+00 
## ab
<rad_m1> implicit transport: solves=2.549000e+03 Picard iterations mean=5.530012e+00 max=3.000000e+01 NON-CONVERGED=0.000000e+00
<rad_m1> picstop: implicit_final_resid_stop=1.000000e-04 (rel 0.000000e+00 x tol_eff) checks=2.551000e+03 failed+continued=2.000000e+00 failed-no-pass-left(NC)=0.000000e+00 solves accepted after a failed check=2.000000e+00 passes after a failed check=2.000000e+00 | implicit_conv_est_guard=0.000000e+00 refused=0.000000e+00 | implicit_conv_est_nofloor=1 refused=7.130000e+02 extra stops=0.000000e+00 
## ref
<rad_m1> implicit transport: solves=2.194000e+03 Picard iterations mean=9.092981e+00 max=4.000000e+01 NON-CONVERGED=0.000000e+00
## off_nl
<rad_m1> implicit transport: solves=1.970000e+03 Picard iterations mean=5.287310e+00 max=1.400000e+01 NON-CONVERGED=0.000000e+00
## b_nl
<rad_m1> implicit transport: solves=2.673000e+03 Picard iterations mean=5.592218e+00 max=3.500000e+01 NON-CONVERGED=0.000000e+00
<rad_m1> picstop: implicit_final_resid_stop=1.000000e-04 (rel 0.000000e+00 x tol_eff) checks=3.184000e+03 failed+continued=5.110000e+02 failed-no-pass-left(NC)=0.000000e+00 solves accepted after a failed check=1.060000e+02 passes after a failed check=1.787000e+03 | implicit_conv_est_guard=0.000000e+00 refused=0.000000e+00 | implicit_conv_est_nofloor=0 refused=0.000000e+00 extra stops=0.000000e+00 
```
