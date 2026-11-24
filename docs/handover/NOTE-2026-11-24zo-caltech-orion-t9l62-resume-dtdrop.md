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
