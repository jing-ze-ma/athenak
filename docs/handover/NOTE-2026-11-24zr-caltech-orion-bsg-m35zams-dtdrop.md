# NOTE 2026-11-24zr caltech -> orion: bsg / M35ZAMS dt drops (your 24zq): dt histories, floor columns, DTREP binders

## bsg: dt every 1000 cycles (last ~15k)
```
cycle=23000 time=1.442974e+06 dt=6.024360e+01
cycle=24000 time=1.502980e+06 dt=6.022177e+01
cycle=25000 time=1.563132e+06 dt=6.008297e+01
cycle=26000 time=1.623077e+06 dt=6.068824e+01
cycle=27000 time=1.683110e+06 dt=5.922333e+01
cycle=28000 time=1.742596e+06 dt=5.234714e+01
cycle=29000 time=1.798967e+06 dt=6.003524e+01
cycle=30000 time=1.858808e+06 dt=6.001345e+01
cycle=31000 time=1.918267e+06 dt=5.999071e+01
cycle=32000 time=1.977281e+06 dt=6.039881e+01
cycle=33000 time=2.035206e+06 dt=5.176114e+01
cycle=34000 time=2.077319e+06 dt=3.218241e+01
cycle=35000 time=2.111660e+06 dt=3.278058e+01
cycle=36000 time=2.146136e+06 dt=3.632307e+01
cycle=37000 time=2.181710e+06 dt=3.210184e+01
cycle=38000 time=2.221404e+06 dt=4.616088e+01
```
## m35zams_z3: dt every 1000 cycles (last ~15k)
```
cycle=54000 time=1.041553e+05 dt=1.534326e+00
cycle=55000 time=1.056979e+05 dt=1.856368e+00
cycle=56000 time=1.070214e+05 dt=1.589882e+00
cycle=57000 time=1.085416e+05 dt=1.640182e+00
cycle=58000 time=1.101009e+05 dt=1.695892e+00
cycle=59000 time=1.114595e+05 dt=9.591549e-01
cycle=60000 time=1.122316e+05 dt=6.918412e-01
cycle=61000 time=1.129352e+05 dt=8.614701e-01
cycle=62000 time=1.139244e+05 dt=1.173690e+00
cycle=63000 time=1.152088e+05 dt=1.130959e+00
cycle=64000 time=1.165886e+05 dt=1.418600e+00
cycle=65000 time=1.179651e+05 dt=8.827438e-01
cycle=66000 time=1.191647e+05 dt=1.259213e+00
cycle=67000 time=1.201687e+05 dt=8.641185e-01
cycle=68000 time=1.217000e+05 dt=1.761222e+00
cycle=69000 time=1.233686e+05 dt=1.769237e+00
```
- bsg's drop is SUDDEN: dt 60 at c32000, 51.8 at c33000, 32.2 at c34000, then 32-36 to c37000 and 46 at c38000. M35ZAMS (Z3) went 1.5-1.7 from c53k to c58k, then 0.96 at c59k / 0.69 at c60k, and has wandered 0.86-1.76 since. (Z3 was 2.48-2.50 to c2000.)
- bsg3d.log event counters (cumulative), header + c0 + latest:
```
#  cycle eos_dfloor eos_efloor eos_tfloor eos_vceil eos_fail c2p_it fofc  efloor_de eos_tclamp   eos_tset   vceil_de  dfloor_ke m1_bcg_break m1_bcg_fall m1_bcg_keep m1_newt_fb m1_gn_switch m1_od_fall 
       0   430944        0        0        0        0      0        0  0.0000e+00          0          0  0.0000e+00  0.0000e+00  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  0.0000
   38437 52905878 113882576  3208212 489494314        0      0 1426268692  1.3628e+08          0          0  0.0000e+00  3.3700e-06  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  0.
```
- m35zams.log (cumulative), latest:
```
   66812   146276  6675204 12876483        0        0      0 9774012760  2.4581e+02          0          0  0.0000e+00  1.0352e-14  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  0.00
```
- t9l62 (your status question): log.4 ended at c5122, one cycle after the resume, on a stale STOP. The prodwatch2 FATAL branch had written that STOP at 05:05, on the dead fresh-start FATAL log.3, and it was never removed. COORD relaunched at 07:16 (log.5, step 7333073.64). I have fixed it: every launcher now does rm -f <star>/STOP before prodg, and prodwatch2 writes STOP only while the step runs. It also now ALERTS when a production step is gone before tlim.
- Next: each star is stopped cleanly and resumed with hydro/dt_report=10 (every key named); the DTREP binders follow below once ~100 cycles are in.

## m35zams_z3: resumed with hydro/dt_report=10 (11-24 07:26)
- Clean STOP: yes; resumed from m35zams.00006.rst (hard link rst_keep/m35zams.00006.rst.pre-dtr) on hpc-sm-01-15: `2026-11-24 07:24:44 chunk 1 STAR=m35zams_z3 BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc3_af836f5a3 md5=09f50abee4ffa6cd399ce15879626` + EX: time/cfl_number=0.3 hydro/fofc=true hydro/fofc_retest=true rad_vet/clean_stag_mg_levels=4 rad_vet/floorcap_lz=true hydro/dt_report=10
- DTREP binders, every other report over ~100 cycles (all 138 DTREP/DTREG lines in caltech_results/DTDROP_1124/m35zams_z3/dtrep_dtreg.txt):
```
c69250 dt_hyd 1.513 gasonly 1.514 dir 3 r 5.5971e+11 rho 1.000e-17 T 1.184e+14 v -4.142e+07 cgas 1.405e+07 ceff 1.405e+07 tauT_h 9.8e-10
c69260 dt_hyd 1.522 gasonly 1.524 dir 3 r 5.6162e+11 rho 1.000e-17 T 9.207e+11 v 5.403e+07 cgas 1.239e+06 ceff 1.323e+06 tauT_h 9.84e-10
c69270 dt_hyd 1.485 gasonly 1.487 dir 3 r 5.6162e+11 rho 1.000e-17 T 1.656e+12 v 5.498e+07 cgas 1.662e+06 ceff 1.726e+06 tauT_h 9.84e-10
c69280 dt_hyd 1.492 gasonly 1.494 dir 3 r 5.6162e+11 rho 1.000e-17 T 1.916e+12 v 5.459e+07 cgas 1.787e+06 ceff 1.847e+06 tauT_h 9.84e-10
c69290 dt_hyd 1.558 gasonly 1.564 dir 3 r 5.6162e+11 rho 1.000e-17 T 1.337e+11 v 5.338e+07 cgas 4.720e+05 ceff 6.632e+05 tauT_h 9.84e-10
c69300 dt_hyd 1.558 gasonly 1.559 dir 3 r 5.6162e+11 rho 1.000e-17 T 1.743e+13 v 4.864e+07 cgas 5.390e+06 ceff 5.410e+06 tauT_h 9.84e-10
c69310 dt_hyd 1.708 gasonly 1.709 dir 3 r 5.6162e+11 rho 1.000e-17 T 9.789e+12 v 4.525e+07 cgas 4.039e+06 ceff 4.066e+06 tauT_h 9.84e-10
c69320 dt_hyd 1.787 gasonly 1.793 dir 3 r 5.5971e+11 rho 1.000e-17 T 1.337e+11 v -4.632e+07 cgas 4.720e+05 ceff 6.539e+05 tauT_h 9.8e-10
c69330 dt_hyd 1.765 gasonly 1.772 dir 3 r 5.5971e+11 rho 1.000e-17 T 1.337e+11 v -4.689e+07 cgas 4.720e+05 ceff 6.538e+05 tauT_h 9.8e-10
c69340 dt_hyd 1.824 gasonly 1.831 dir 3 r 5.5971e+11 rho 1.000e-17 T 1.337e+11 v -4.537e+07 cgas 4.720e+05 ceff 6.538e+05 tauT_h 9.8e-10
c69350 dt_hyd 1.792 gasonly 1.799 dir 3 r 5.5971e+11 rho 1.000e-17 T 1.337e+11 v -4.618e+07 cgas 4.720e+05 ceff 6.538e+05 tauT_h 9.8e-10
c69360 dt_hyd 1.66 gasonly 1.661 dir 3 r 5.5943e+11 rho 1.013e-17 T 4.383e+12 v 4.780e+07 cgas 2.703e+06 ceff 2.741e+06 tauT_h 9.91e-10
```
- first DTREG lines: `DTREG polar cycle=69250 npol=7 cfl*cell_dt=1.931743e+00 row=4 theta_deg=8.829309e+01 phi=2.425000e-02 r=5.588872e+11 dir=3 rho=1.000000e-17 T=4.333303e+12 v=4.065206e+07 cgas=2.687410e+06 ceff=2.726319e+06 h=2.793196e+08 cfl*cell_dt_gas=1.933477e+00|DTREG rest  cycle=69250 npol=7 cfl*cell_dt=1.513318e+00 row=56 theta_deg=8.978484e+01 phi=2.375000e-02 r=5.597052e+11 dir=3 rho=1.000000e-17 T=1.183834e+14 v=-4.142355e+07 cgas=1.404655e+07 ceff=1.405400e+07 h=2.798506e+08 cfl*cell_dt_gas=1.513522e+0|`

## bsg: resumed with hydro/dt_report=10 (11-24 07:27)
- Clean STOP: yes; resumed from bsg3d.00028.rst (hard link rst_keep/bsg3d.00028.rst.pre-dtr) on hpc-sm-02-16: `2026-11-24 07:25:26 chunk 3 STAR=bsg BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc3_af836f5a3 md5=09f50abee4ffa6cd399ce15879626353 DEC` + EX: time/cfl_number=0.3 hydro/fofc=true hydro/fofc_retest=true rad_vet/clean_stag_mg_levels=4 hydro/dt_report=10
- DTREP binders, every other report over ~100 cycles (all 84 DTREP/DTREG lines in caltech_results/DTDROP_1124/bsg/dtrep_dtreg.txt):
```
c38670 dt_hyd 46.02 gasonly 38.98 dir 1 r 3.9242e+12 rho 6.152e-14 T 1.045e+15 v -9.999e+07 cgas 4.174e+07 ceff 4.174e+07 tauT_h 0.000438
c38680 dt_hyd 49.47 gasonly 40.53 dir 1 r 3.8209e+12 rho 6.235e-14 T 1.984e+13 v -6.017e+07 cgas 5.751e+06 ceff 5.751e+06 tauT_h 0.000228
c38690 dt_hyd 44.33 gasonly 40.66 dir 1 r 4.2286e+12 rho 6.152e-14 T 1.336e+15 v -9.996e+07 cgas 4.718e+07 ceff 4.718e+07 tauT_h 0.000428
c38700 dt_hyd 46.83 gasonly 41 dir 1 r 3.8209e+12 rho 6.232e-14 T 2.777e+13 v -6.283e+07 cgas 6.803e+06 ceff 6.803e+06 tauT_h 0.000228
c38710 dt_hyd 45.86 gasonly 37.66 dir 1 r 4.0981e+12 rho 6.152e-14 T 1.071e+15 v -9.999e+07 cgas 4.224e+07 ceff 4.224e+07 tauT_h 0.000428
c38720 dt_hyd 45.89 gasonly 41.35 dir 1 r 4.2286e+12 rho 6.152e-14 T 1.066e+15 v -9.999e+07 cgas 4.214e+07 ceff 4.214e+07 tauT_h 0.000421
c38730 dt_hyd 46.35 gasonly 40.37 dir 1 r 4.0112e+12 rho 6.152e-14 T 1.695e+15 v -8.758e+07 cgas 5.315e+07 ceff 5.315e+07 tauT_h 0.000434
c38740 dt_hyd 46.63 gasonly 36.02 dir 1 r 3.8209e+12 rho 6.238e-14 T 2.764e+13 v -6.315e+07 cgas 6.788e+06 ceff 6.788e+06 tauT_h 0.000228
c38750 dt_hyd 46.74 gasonly 37.33 dir 1 r 3.8209e+12 rho 6.235e-14 T 2.751e+13 v -6.301e+07 cgas 6.772e+06 ceff 6.772e+06 tauT_h 0.000228
c38760 dt_hyd 47.07 gasonly 42.66 dir 1 r 3.8209e+12 rho 6.231e-14 T 2.628e+13 v -6.267e+07 cgas 6.618e+06 ceff 6.618e+06 tauT_h 0.000228
c38770 dt_hyd 51.06 gasonly 41.16 dir 1 r 4.2503e+12 rho 6.152e-14 T 6.959e+14 v -9.368e+07 cgas 3.406e+07 ceff 3.406e+07 tauT_h 0.000421
c38780 dt_hyd 45 gasonly 35.04 dir 1 r 3.9242e+12 rho 6.152e-14 T 1.216e+15 v -9.993e+07 cgas 4.501e+07 ceff 4.501e+07 tauT_h 0.000423
```
- first DTREG lines: `DTREG polar cycle=38670 npol=7 cfl*cell_dt=84.7994 row=0 theta_deg=60.1172 phi=0.00204531 r=1.40233e+12 dir=3 rho=5.77419e-05 T=1.41913e+14 v=-1775.14 cgas=1.53792e+07 ceff=1.7594e+07 h=4.97371e+09 cfl*cell_dt_gas=97.01|DTREG rest  cycle=38670 npol=7 cfl*cell_dt=86.1643 row=7 theta_deg=61.7579 phi=0.00204531 r=1.40233e+12 dir=3 rho=5.77393e-05 T=1.41911e+14 v=-923.111 cgas=1.53791e+07 ceff=1.75939e+07 h=5.05349e+09 cfl*cell_dt_gas=98.5723|`

## Reading (24zq items 3 / 4)
- **bsg: the binder IS floor-density gas at large r, the c777 mechanism.** rho = 6.152e-14 = dfloor (deck dfloor 6.15245e-14), r 3.8-4.25e12, v_r -6e7 to **-1.0e8 = vceil (1.0e8)**, dir 1, tauT_h ~4e-4; cfl*dt_hyd 44-51.
  The other half of the binders sits at rho 6.23e-14, r 3.82e12 with |v| 6.3e7. **DFLFIX H3 started beside bsg at once** (07:28, hpc-sm-01-15 next to M35ZAMS; memory 34 + 40 GB/GPU), athena_dflfix_5b6159ad0 from bsg3d.00028 (hard link .dfltest), production keys + vet_gd_ls=0, 1000 cycles. Results go to NOTE 24zn (H3).
- **M35ZAMS (Z3): the binder is floor gas, but ABOVE the cap, not at it.** rho = 1.000e-17 (dfloor), r 5.59-5.62e11 (the thin layer just above the photosphere), **dir 3 (phi), |v| 4.1-5.5e7 > floorcap_v 3e7**, cfl*dt_hyd 1.49-1.82.
  This is the same mechanism as m35mms before floorcap_lz_spin: floorcap_lz keeps v_phi of capped floor cells.
  Your FS2 arm (floorcap_lz_spin on the Z3 state at c17560) held dt at 2.477 with the floor gas AT the cap (|v| = 3e7). So for M35ZAMS: floorcap_lz_spin would bring the binder back to "floor gas at the cap" (dt ~2.5, the dt_floor_courant item = USER decision).
  Without it, the spun-up floor layer binds dt at ~1.5-1.8. Not switched (no switch on dt alone).
- Both runs keep hydro/dt_report=10 for now; I drop it at their next rst if it costs > 1 %.

## m35zams_z3 (11-24 07:34): switched (spin)
- pre-switch (log.1): 1.0096 s/cyc (c70240-70290), dt 1.612501e+00. Clean STOP: yes. Pre-switch rst kept as rst_keep/m35zams.00007.rst.pre-spin.
- relaunch: `2026-11-24 07:32:59 chunk 2 STAR=m35zams_z3 BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc3fcs_6f2973dc8 md5=7cd2856ec083818cc35742d6cd668550 DECK=/resnick/groups/carnegie_poc/jingze/prod_rc3_1124/decks/ti` + EX: time/cfl_number=0.3 hydro/fofc=true hydro/fofc_retest=true rad_vet/clean_stag_mg_levels=4 rad_vet/floorcap_lz=true rad_vet/floorcap_lz_spin=true rad_vet/floorcap_mesh_vol=true hydro/dt_report=10
- echo: 'floorcap_lz_spin' in block 'rad_vet' is absent from the restart's inp;
- m35zams_z3 300 cycles after the switch (11-24 07:38): CW 0, FATAL 0, created inf 0, DT_RETRY 0. **before** 1.0096 s/cyc (c70240-70290), dt 1.612501e+00 -> **after** 0.2755 s/cyc (c70630-70680), dt 2.488757e+00. Last DTREP: binder c70680 dir 1 r 5.5780e+11 rho 1.000e-17 v -2.999e+07 ceff 2.795e+06 (cfl*dt_hyd 2.48875)
