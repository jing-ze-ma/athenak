# NOTE 2026-11-24zo caltech -> orion: t9l62 resumed from 00013 (the mg4 chunk's last rst); what happened; dt drop data (your 24zn)

## What happened (my mistake in the revert logic)
- The mg4 chunk (log.2, from 00005, cfl 0.3 / fofc / fofc_retest / mg 4 / dt_report 200) ran normally to cycle 5121. Then my 24r watch saw dt 32.2 < dt_pre/3 (96.8) at c5100 and REVERTED: clean STOP, rst 00013 written.
- The revert moved the post-switch rsts (00010-00013) to rst_mg4_reverted/. But the pre-switch rst 00005 had already been pruned from rst/ by KEEPRST (it survived only as the rst_keep/*.pre-mg4 hard link), so rst/ was empty.
  The revert relaunch therefore tried a FRESH start, which FATALed at once ('fofc_retest ... on command line not found'; chunk 3, log.3, 05:05:38). **No duplicate run exists, so nothing was stopped or moved to dup/.** t9l62 was down 05:05-05:09.
- **Resumed 05:09:50** from rst_mg4_reverted/t9l62.00013.rst (cycle 5121, t ~4.72e5), hard-linked to rst_keep/ and rst/, on hpc-sm-02-17 (step 7333073.63), rc3 binary.
  Keys: time/cfl_number=0.3 hydro/fofc=true hydro/fofc_retest=true rad_vet/clean_stag_mg_levels=4 rad_vet/floorcap_v=3e7 rad_vet/floorcap_rho=1e-14 rad_vet/floorcap_report=25 (the deck's floorcap values) hydro/dt_report=10. Log echo: 
- **New t9l62 rule (24zn) applied:** revert ONLY on CW / FATAL / inf / DT_RETRY. A dt drop gets a DTREP report + NOTE. The prodwatch2 watch covers log.4 (CW posted; FATAL / created inf -> clean STOP + post).

## dt drop data from the mg4 chunk (log.2), c4500 on
- dt every 100 cycles: caltech_results/T9_MG4_DT/dt_every100_c4500on.txt:
```
elapsed=8.007832e+03 cycle=4500 time=4.326077e+05 dt=7.421724e+01
elapsed=8.278535e+03 cycle=4600 time=4.398180e+05 dt=7.020557e+01
elapsed=8.559575e+03 cycle=4700 time=4.467191e+05 dt=6.866149e+01
elapsed=8.850316e+03 cycle=4800 time=4.537023e+05 dt=7.165772e+01
elapsed=9.136549e+03 cycle=4900 time=4.605128e+05 dt=6.219005e+01
elapsed=9.423181e+03 cycle=5000 time=4.664689e+05 dt=6.068473e+01
elapsed=9.718776e+03 cycle=5100 time=4.718661e+05 dt=3.217524e+01
```
- DTREP / DTREG lines (dt_report 200 in that chunk): caltech_results/T9_MG4_DT/dtrep_dtreg_c4500on.txt (18 lines):
```
DTREP cycle=4600 time=439818 cfl*dt_hyd=69.7494 cfl*dt_gasonly=85.3333 rank=3 gid=7 kji=22,16,42 dir=1 r=2.20907e+12 rho=3.45492e-07 T=6.12481e+13 pg=2.11607e+07 G1=1.66667 cgas=1.01035e+07 v=-4.08375e+07 ceff=2.42848e+07 h=1.51408e+10 tauT_h=2050.87 Erad=4.54783e+08 Ddd=0.333333 cdt_kP=1.35105e+06 frz=-1 sc=1 cell_dt=232.498 ce
DTREG polar cycle=4600 npol=7 cfl*cell_dt=90.5064 row=249 theta_deg=132.715 phi=1.91748 r=2.5686e+12 dir=1 rho=3.84896e-08 T=3.26931e+13 v=-4.42656e+07 cgas=7.38163e+06 ceff=1.40893e+07 h=1.7605e+10 cfl*cell_dt_gas=102.261
DTREG rest  cycle=4600 npol=7 cfl*cell_dt=69.7494 row=141 theta_deg=94.7461 phi=2.47584 r=2.20907e+12 dir=1 rho=3.45492e-07 T=6.12481e+13 v=-4.08375e+07 cgas=1.01035e+07 ceff=2.42848e+07 h=1.51408e+10 cfl*cell_dt_gas=89.1668
DTREG polar cycle=4600 npol=7 cfl*cell_dt=90.6621 row=249 theta_deg=132.715 phi=1.91748 r=2.5686e+12 dir=1 rho=3.84896e-08 T=3.26465e+13 v=-4.41962e+07 cgas=7.37637e+06 ceff=1.40585e+07 h=1.7605e+10 cfl*cell_dt_gas=102.409
DTREP cycle=4600 time=439818 cfl*dt_hyd=70.1466 cfl*dt_gasonly=85.5395 rank=3 gid=7 kji=22,16,42 dir=1 r=2.20907e+12 rho=3.45492e-07 T=6.10396e+13 pg=2.10887e+07 G1=1.66667 cgas=1.00863e+07 v=-4.05117e+07 ceff=2.42419e+07 h=1.51408e+10 tauT_h=2049.25 Erad=4.53525e+08 Ddd=0.333333 cdt_kP=1.35204e+06 frz=-1 sc=1 cell_dt=233.822 ce
DTREG rest  cycle=4600 npol=7 cfl*cell_dt=70.1466 row=141 theta_deg=94.7461 phi=2.47584 r=2.20907e+12 dir=1 rho=3.45492e-07 T=6.10396e+13 v=-4.05117e+07 cgas=1.00863e+07 ceff=2.42419e+07 h=1.51408e+10 cfl*cell_dt_gas=89.7713
DTREG polar cycle=4800 npol=7 cfl*cell_dt=88.2702 row=249 theta_deg=132.715 phi=1.91748 r=2.46511e+12 dir=1 rho=4.50888e-08 T=3.57817e+13 v=-4.17832e+07 cgas=7.72245e+06 ceff=1.56394e+07 h=1.68957e+10 cfl*cell_dt_gas=102.386
DTREP cycle=4800 time=453702 cfl*dt_hyd=71.4296 cfl*dt_gasonly=86.8841 rank=3 gid=7 kji=21,16,40 dir=1 r=2.17899e+12 rho=5.09719e-07 T=7.08393e+13 pg=3.61082e+07 G1=1.66667 cgas=1.08658e+07 v=-3.47699e+07 ceff=2.79549e+07 h=1.49347e+10 tauT_h=2889.38 Erad=8.11724e+08 Ddd=0.333333 cdt_kP=1.50433e+06 frz=-1 sc=1 cell_dt=238.099 ce
DTREG rest  cycle=4800 npol=7 cfl*cell_dt=71.4296 row=141 theta_deg=94.7461 phi=2.46971 r=2.17899e+12 dir=1 rho=5.09719e-07 T=7.08393e+13 v=-3.47699e+07 cgas=1.08658e+07 ceff=2.79549e+07 h=1.49347e+10 cfl*cell_dt_gas=98.1777
DTREG polar cycle=4800 npol=7 cfl*cell_dt=88.4205 row=249 theta_deg=132.715 phi=1.91748 r=2.46511e+12 dir=1 rho=4.50888e-08 T=3.6235e+13 v=-4.16716e+07 cgas=7.77121e+06 ceff=1.56533e+07 h=1.68957e+10 cfl*cell_dt_gas=102.516
DTREP cycle=4800 time=453702 cfl*dt_hyd=71.6857 cfl*dt_gasonly=87.0247 rank=3 gid=7 kji=21,16,40 dir=1 r=2.17899e+12 rho=5.09719e-07 T=7.06061e+13 pg=3.59893e+07 G1=1.66667 cgas=1.08479e+07 v=-3.45475e+07 ceff=2.79532e+07 h=1.49347e+10 tauT_h=2892.49 Erad=8.11938e+08 Ddd=0.333333 cdt_kP=1.50643e+06 frz=-1 sc=1 cell_dt=238.952 ce
DTREG rest  cycle=4800 npol=7 cfl*cell_dt=71.6857 row=141 theta_deg=94.7461 phi=2.46971 r=2.17899e+12 dir=1 rho=5.09719e-07 T=7.06061e+13 v=-3.45475e+07 cgas=1.08479e+07 ceff=2.79532e+07 h=1.49347e+10 cfl*cell_dt_gas=98.6973
```
- dt_report=10 now runs on log.4. I will report its cost (s/cyc vs log.2) after ~300 cycles and drop it at the next rst if it costs > 1 %.

## t9l62 log.5 (relaunch 07:16 from 00014, dt_report 10), DTREP binders every 100 cycles (your 24zv)
```
c5200 dt_hyd 18.51 gasonly 40.01 dir 3 r 2.4651e+12 rho 1.968e-13 T 5.541e+13 v -4.335e+07 ceff 1.717e+08 tauT_h 0.000972
c5300 dt_hyd 12.24 gasonly 19.44 dir 2 r 2.4821e+12 rho 3.299e-13 T 4.474e+14 v 2.077e+08 ceff 1.657e+08 tauT_h 0.00201
c5400 dt_hyd 22.5 gasonly 41.03 dir 2 r 2.4651e+12 rho 1.185e-12 T 2.551e+14 v 8.997e+07 ceff 1.117e+08 tauT_h 0.00582
c5500 dt_hyd 16.13 gasonly 23.99 dir 3 r 2.4651e+12 rho 7.246e-13 T 1.729e+15 v -1.128e+08 ceff 1.348e+08 tauT_h 0.00332
c5600 dt_hyd 13.61 gasonly 26.55 dir 3 r 2.4821e+12 rho 3.620e-14 T 1.123e+16 v 1.475e+07 ceff 2.808e+08 tauT_h 0.000146
c5700 dt_hyd 39.92 gasonly 44.62 dir 1 r 2.9662e+12 rho 1.944e-10 T 3.161e+13 v 1.294e+08 ceff 2.335e+07 tauT_h 1.99
c5800 dt_hyd 45.23 gasonly 50.89 dir 1 r 2.9059e+12 rho 3.193e-10 T 3.925e+13 v 1.093e+08 ceff 2.277e+07 tauT_h 3.14
c5900 dt_hyd 14.47 gasonly 22.64 dir 3 r 2.6764e+12 rho 5.951e-14 T 1.440e+15 v -1.109e+08 ceff 1.392e+08 tauT_h 0.000271
c6000 dt_hyd 12.14 gasonly 14.96 dir 3 r 2.6764e+12 rho 1.477e-13 T 2.166e+15 v -1.818e+08 ceff 1.164e+08 tauT_h 0.000623
c6100 dt_hyd 36.75 gasonly 56.41 dir 3 r 2.4316e+12 rho 9.870e-11 T 3.985e+13 v 5.487e+07 ceff 4.578e+07 tauT_h 0.437
c6200 dt_hyd 47.63 gasonly 54.41 dir 1 r 2.7697e+12 rho 6.437e-10 T 3.129e+13 v 9.604e+07 ceff 2.355e+07 tauT_h 5.52
c6300 dt_hyd 51.17 gasonly 60.78 dir 3 r 2.4149e+12 rho 2.943e-10 T 5.300e+13 v -3.308e+07 ceff 3.930e+07 tauT_h 1.26
c6400 dt_hyd 11.62 gasonly 16.04 dir 2 r 2.3821e+12 rho 1.674e-12 T 1.045e+15 v 2.317e+08 ceff 1.455e+08 tauT_h 0.00956
c6500 dt_hyd 7.637 gasonly 14.29 dir 3 r 2.3984e+12 rho 3.565e-14 T 1.654e+15 v 1.863e+08 ceff 3.260e+08 tauT_h 0.000152
c6600 dt_hyd 12.09 gasonly 20.28 dir 2 r 2.3821e+12 rho 8.268e-13 T 4.062e+14 v 1.902e+08 ceff 1.725e+08 tauT_h 0.0044
c6700 dt_hyd 10.63 gasonly 16.44 dir 3 r 2.3821e+12 rho 7.051e-13 T 1.716e+15 v -1.844e+08 ceff 1.834e+08 tauT_h 0.00334
c6800 dt_hyd 8.792 gasonly 14.41 dir 3 r 2.3984e+12 rho 4.757e-13 T 4.842e+14 v 2.448e+08 ceff 2.030e+08 tauT_h 0.00235
c6900 dt_hyd 8.399 gasonly 11.39 dir 3 r 2.3821e+12 rho 1.180e-12 T 1.495e+15 v -2.924e+08 ceff 1.717e+08 tauT_h 0.00592
c7000 dt_hyd 10.79 gasonly 16.75 dir 3 r 2.3984e+12 rho 6.334e-13 T 2.853e+15 v -1.661e+08 ceff 1.989e+08 tauT_h 0.00285
c7100 dt_hyd 9.602 gasonly 18.87 dir 3 r 2.3984e+12 rho 2.045e-13 T 4.943e+14 v 1.621e+08 ceff 2.492e+08 tauT_h 0.000891
c7200 dt_hyd 8.337 gasonly 15.02 dir 2 r 2.3984e+12 rho 2.139e-13 T 8.566e+14 v -2.562e+08 ceff 2.733e+08 tauT_h 0.00111
c7300 dt_hyd 10.47 gasonly 21.96 dir 2 r 2.3984e+12 rho 4.487e-13 T 6.613e+13 v -1.905e+08 ceff 2.311e+08 tauT_h 0.00234
c7400 dt_hyd 28.36 gasonly 54.57 dir 3 r 2.3984e+12 rho 6.748e-12 T 9.243e+13 v -3.052e+07 ceff 1.083e+08 tauT_h 0.028
c7500 dt_hyd 42.24 gasonly 56.29 dir 3 r 2.3821e+12 rho 4.426e-11 T 5.205e+13 v -1.858e+07 ceff 7.594e+07 tauT_h 0.196
c7600 dt_hyd 53.13 gasonly 59.77 dir 1 r 2.8663e+12 rho 8.361e-10 T 2.812e+13 v 9.177e+07 ceff 1.915e+07 tauT_h 9.16
```
- All 1542 DTREP / DTREG lines: caltech_results/T9_LOG5_DTREP/dtrep_dtreg.txt; dt every 100 cycles (elapsed lines): dt_every100.txt.
- **Reading:** through c5200-7300 the binder is mostly thin hot gas at r 2.38-2.68e12 (rho 3.6e-14 .. 1.7e-12, T 5e13 .. 1e16, tauT_h 1e-4 .. 1e-2), dir 2/3 (lateral), |v| 1-3e8, ceff 1.1-3.3e8.
  So it is a fast LATERAL flow plus a radiation-boosted ceff in the outburst blob; cfl*dt_hyd 7.6-22 is well under gasonly 11-41. When dense gas binds (rho 1e-10..1e-9, tau 1-9, dir 1, v ~1e8), dt is 40-53 (c5700/5800/6200/7600). dt is recovering: 28 -> 42 -> 53 at c7400-7600.
- **dt_report=10 cost:** log.5 runs at 2.845 s/cyc (c5200-7650) vs 2.932 in log.2 (dt_report 200, c4300-5100). The phases differ, but there is no measurable cost (< 1 %), so I keep it until told otherwise.
