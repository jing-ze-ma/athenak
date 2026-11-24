# NOTE 2026-11-24zzx caltech -> orion: Compton t9l62 arms C0 / C1 / C2 / C0b / C1r (your 24zzs) on bsg's node 02-11 (bsg paused, USER approved)

- Binary caltech-t9hot-compton-f8e55fee7 (CUDA 12.9 + g++ 13, HOPPER90, incremental on the T2 tree rel38f/build_rel): BUILD_OK md5=a2e018822f8dbc70fd0b8f27ef0e9b40; GATE PASS (gate deck copy WITHOUT vet_gd_ls): athena_t9hotcompton_f8e55fee7 md5 a2e01882 node hpc-sm-02-11 pass1=3.686e-11
- **rst: NOT rc3 00020.final_1124** (that one is at cycle 9572, t ~6.69e5, dt 97: past the small-dt phase; your note says t ~5.7e5). Used prod_rc3_1124/t9l62/rst_keep/t9l62.00015.rst.t9hot (cycle 5279, t 4.75210e5, dt 13.9 (rst header): the eruption), read only, hard-linked as t9l62.00015.rst. Args: ia123rc4_decks/IA123_t9l62_restart.args (all keys, vet_gd_ls=0) + hydro/dt_report=10; outputs 2/3/5/6/7 off. 4 GPUs per arm, one after the other (t9l62 ~104 GB/GPU on 4).
- Stored rc3 keys that the args do not name keep their rc3 values in ALL arms (same for C0 and C1): wb_cache_every=0, cs_seam_blockend/rho_guard/cs_vertex_3cell, he_ic_*, user_bcs, clean_face_hclat=0, clean_ray_source_relaxed=0, clean_stag_kpat=false, clean_stag_thinf_held=false, clean_thin_donor=0, time2_stage_local_be=false, time2_start=be, vet_col_lat_taucut=30, vet_gd_rotate_align=true, vet_rotate_avg_ema_change_form=linear, dt_retry_max=8 (the T2 decks drop these keys).

## C0 (2026-11-24 12:54:33 rc=0 wall=962 s)
- CW 0, FATAL 0, created inf 0, DT_RETRY 0; s/cyc 3.0007 (c5350-5550)
- dt every 10 cycles: 13.94 24.97 14.78 17.84 13.96 13.87 13.67 12.99 29.44 32.66 30.03 28.89 38.15 34.43 34.99 37.97 30.26 42.77 31.04 43.01 42.31 41.65 41.88 44.23 46.16 43.32 42.11 42.24 43.48 44.46 
- DTREP / DTREG every 50 cycles (all in COMPTON/C0/dtrep.txt):
```
DTREP cycle=5300 time=4.756428e+05 cfl*dt_hyd=1.782543e+01 cfl*dt_gasonly=1.782935e+01 rank=0 gid=0 kji=79,50,58 dir=3 r=2.465107e+12 rho=3.714662e-12 T=1.840681e+15 pg=6.837507e+03 G1=1.666667e+00 cgas=5.538773e+07 v=1.686993e+08 ceff=5.543703e+07 h=1.331776e+10 tauT_h=1.423448e-02 Erad=3.961774e+0
DTREG polar cycle=5300 npol=7 cfl*cell_dt=7.734394e+01 row=3 theta_deg=4.623064e+01 phi=1.810097e-01 r=2.658148e+12 dir=1 rho=2.946180e-09 T=3.283257e+13 v=6.184364e+07 cgas=7.397361e+06 ceff=8.822896e+06 h=1.821876e+10 cfl*cell_dt_gas=7.893630e+01
DTREG rest  cycle=5300 npol=7 cfl*cell_dt=1.782543e+01 row=47 theta_deg=6.169932e+01 phi=4.693981e-01 r=2.465107e+12 dir=3 rho=3.714662e-12 T=1.840681e+15 v=1.686993e+08 cgas=5.538773e+07 ceff=5.543703e+07 h=1.331776e+10 cfl*cell_dt_gas=1.782935e+01
DTREP cycle=5300 time=4.756428e+05 cfl*dt_hyd=1.804085e+01 cfl*dt_gasonly=1.804638e+01 rank=0 gid=0 kji=79,50,58 dir=3 r=2.465107e+12 rho=3.714662e-12 T=1.624963e+15 pg=6.036188e+03 G1=1.666667e+00 cgas=5.204106e+07 v=1.693511e+08 ceff=5.210899e+07 h=1.331776e+10 tauT_h=1.624155e-02 Erad=3.951490e+0
DTREG polar cycle=5300 npol=7 cfl*cell_dt=7.726963e+01 row=3 theta_deg=4.623064e+01 phi=1.810097e-01 r=2.658148e+12 dir=1 rho=2.946180e-09 T=3.258799e+13 v=6.193373e+07 cgas=7.369757e+06 ceff=8.800753e+06 h=1.821876e+10 cfl*cell_dt_gas=7.886512e+01
DTREG rest  cycle=5300 npol=7 cfl*cell_dt=1.804085e+01 row=47 theta_deg=6.169932e+01 phi=4.693981e-01 r=2.465107e+12 dir=3 rho=3.714662e-12 T=1.624963e+15 v=1.693511e+08 cgas=5.204106e+07 ceff=5.210899e+07 h=1.331776e+10 cfl*cell_dt_gas=1.804638e+01
DTREP cycle=5350 time=4.763895e+05 cfl*dt_hyd=1.449179e+01 cfl*dt_gasonly=1.449736e+01 rank=0 gid=0 kji=78,50,58 dir=3 r=2.465107e+12 rho=1.213093e-13 T=1.622858e+15 pg=1.968678e+02 G1=1.666667e+00 cgas=5.200734e+07 v=2.235827e+08 ceff=5.211314e+07 h=1.331776e+10 tauT_h=6.431792e-04 Erad=4.099502e+0
DTREG polar cycle=5350 npol=7 cfl*cell_dt=7.561188e+01 row=4 theta_deg=4.658220e+01 phi=1.810097e-01 r=2.658148e+12 dir=1 rho=2.751749e-09 T=3.285643e+13 v=6.347389e+07 cgas=7.400049e+06 ceff=8.811418e+06 h=1.821876e+10 cfl*cell_dt_gas=7.711760e+01
DTREG rest  cycle=5350 npol=7 cfl*cell_dt=1.449179e+01 row=47 theta_deg=6.169932e+01 phi=4.632622e-01 r=2.465107e+12 dir=3 rho=1.213093e-13 T=1.622858e+15 v=2.235827e+08 cgas=5.200734e+07 ceff=5.211314e+07 h=1.331776e+10 cfl*cell_dt_gas=1.449736e+01
DTREP cycle=5350 time=4.763895e+05 cfl*dt_hyd=1.502074e+01 cfl*dt_gasonly=1.502376e+01 rank=0 gid=0 kji=77,51,59 dir=3 r=2.482061e+12 rho=2.516662e-13 T=2.648729e+15 pg=6.665955e+02 G1=1.666667e+00 cgas=6.644207e+07 v=-2.022005e+08 ceff=6.649596e+07 h=1.345340e+10 tauT_h=1.091526e-03 Erad=4.066722e+
DTREG polar cycle=5350 npol=7 cfl*cell_dt=7.555207e+01 row=4 theta_deg=4.658220e+01 phi=1.810097e-01 r=2.658148e+12 dir=1 rho=2.751749e-09 T=3.262832e+13 v=6.355666e+07 cgas=7.374316e+06 ceff=8.785871e+06 h=1.821876e+10 cfl*cell_dt_gas=7.705559e+01
DTREG rest  cycle=5350 npol=7 cfl*cell_dt=1.502074e+01 row=48 theta_deg=6.205088e+01 phi=4.571263e-01 r=2.482061e+12 dir=3 rho=2.516662e-13 T=2.648729e+15 v=-2.022005e+08 cgas=6.644207e+07 ceff=6.649596e+07 h=1.345340e+10 cfl*cell_dt_gas=1.502376e+01
DTREP cycle=5400 time=4.779036e+05 cfl*dt_hyd=3.197210e+01 cfl*dt_gasonly=3.200816e+01 rank=0 gid=0 kji=77,50,58 dir=1 r=2.465107e+12 rho=9.664888e-13 T=1.265805e+15 pg=1.223387e+03 G1=1.666667e+00 cgas=4.593120e+07 v=-1.124254e+08 ceff=4.610980e+07 h=1.689568e+10 tauT_h=6.602093e-03 Erad=4.178428e+
DTREG polar cycle=5400 npol=7 cfl*cell_dt=7.256597e+01 row=5 theta_deg=4.693376e+01 phi=1.810097e-01 r=2.639991e+12 dir=1 rho=2.381367e-09 T=3.321570e+13 v=6.598849e+07 cgas=7.440396e+06 ceff=8.816476e+06 h=1.809432e+10 cfl*cell_dt_gas=7.392588e+01
DTREG rest  cycle=5400 npol=7 cfl*cell_dt=3.197210e+01 row=47 theta_deg=6.169932e+01 phi=4.571263e-01 r=2.465107e+12 dir=1 rho=9.664888e-13 T=1.265805e+15 v=-1.124254e+08 cgas=4.593120e+07 ceff=4.610980e+07 h=1.689568e+10 cfl*cell_dt_gas=3.200816e+01
DTREP cycle=5400 time=4.779036e+05 cfl*dt_hyd=3.605146e+01 cfl*dt_gasonly=3.624830e+01 rank=0 gid=0 kji=77,50,58 dir=1 r=2.465107e+12 rho=9.664888e-13 T=3.362099e+13 pg=3.249431e+01 G1=1.666667e+00 cgas=7.485652e+06 v=-1.323472e+08 ceff=8.249159e+06 h=1.689568e+10 tauT_h=5.641760e-03 Erad=4.180331e+
```
- compton report lines:
```
```
- m1-positivity line: ``
- hst tail (t dt L_bot L_top L_in Etot M_tot Lz_tot):
```
4.850275e+05 45.389 8.627974e+38 2.139591e+39 8.763705e+38 2.1540578376e+47 2.3310565725e+32 -6.9644886744e+38
4.850730e+05 45.492 8.628668e+38 2.139073e+39 8.763705e+38 2.1540572632e+47 2.3310565725e+32 -6.9644897744e+38
```

## C1 (2026-11-24 13:10:48 rc=0 wall=975 s)
- CW 0, FATAL 0, created inf 0, DT_RETRY 0; s/cyc 3.0074 (c5350-5550)
- dt every 10 cycles: 13.94 26.45 15.74 16.75 17.13 22.61 28.92 42.66 23.91 31.23 25.37 42.34 33.09 19.26 30.57 43.1 42.8 34.96 42.63 41.99 41.52 42.66 45.89 44.68 42.8 41.86 42.85 43.93 44.87 46.26 
- DTREP / DTREG every 50 cycles (all in COMPTON/C1/dtrep.txt):
```
DTREP cycle=5300 time=4.757182e+05 cfl*dt_hyd=1.485293e+01 cfl*dt_gasonly=1.485487e+01 rank=0 gid=0 kji=78,51,59 dir=2 r=2.482061e+12 rho=5.460057e-13 T=6.372609e+15 pg=3.479481e+03 G1=1.666667e+00 cgas=1.030583e+08 v=2.045123e+08 ceff=1.030985e+08 h=1.522974e+10 tauT_h=2.623230e-03 Erad=4.038165e+0
DTREG polar cycle=5300 npol=7 cfl*cell_dt=7.739349e+01 row=3 theta_deg=4.623064e+01 phi=1.810097e-01 r=2.658148e+12 dir=1 rho=2.954438e-09 T=3.284690e+13 v=6.179664e+07 cgas=7.398976e+06 ceff=8.824652e+06 h=1.821876e+10 cfl*cell_dt_gas=7.898807e+01
DTREG rest  cycle=5300 npol=7 cfl*cell_dt=1.485293e+01 row=48 theta_deg=6.205088e+01 phi=4.632622e-01 r=2.482061e+12 dir=2 rho=5.460057e-13 T=6.372609e+15 v=2.045123e+08 cgas=1.030583e+08 ceff=1.030985e+08 h=1.522974e+10 cfl*cell_dt_gas=1.485487e+01
DTREP cycle=5300 time=4.757182e+05 cfl*dt_hyd=1.736667e+01 cfl*dt_gasonly=1.740613e+01 rank=0 gid=0 kji=78,50,58 dir=3 r=2.465107e+12 rho=1.244866e-12 T=3.334409e+13 pg=4.150893e+01 G1=1.666667e+00 cgas=7.454762e+06 v=2.220809e+08 ceff=7.976381e+06 h=1.331776e+10 tauT_h=5.741939e-03 Erad=4.039057e+0
DTREG polar cycle=5300 npol=7 cfl*cell_dt=7.731479e+01 row=3 theta_deg=4.623064e+01 phi=1.810097e-01 r=2.658148e+12 dir=1 rho=2.954438e-09 T=3.258579e+13 v=6.189193e+07 cgas=7.369508e+06 ceff=8.801243e+06 h=1.821876e+10 cfl*cell_dt_gas=7.891300e+01
DTREG rest  cycle=5300 npol=7 cfl*cell_dt=1.736667e+01 row=47 theta_deg=6.169932e+01 phi=4.632622e-01 r=2.465107e+12 dir=3 rho=1.244866e-12 T=3.334409e+13 v=2.220809e+08 cgas=7.454762e+06 ceff=7.976381e+06 h=1.331776e+10 cfl*cell_dt_gas=1.740613e+01
DTREP cycle=5350 time=4.769649e+05 cfl*dt_hyd=3.465143e+01 cfl*dt_gasonly=3.469290e+01 rank=0 gid=0 kji=77,52,60 dir=1 r=2.499132e+12 rho=1.289636e-12 T=6.717361e+14 pg=8.662948e+02 G1=1.666667e+00 cgas=3.345983e+07 v=1.146587e+08 ceff=3.363710e+07 h=1.712888e+10 tauT_h=7.509563e-03 Erad=4.067619e+0
DTREG polar cycle=5350 npol=7 cfl*cell_dt=7.420428e+01 row=4 theta_deg=4.658220e+01 phi=1.810097e-01 r=2.658148e+12 dir=1 rho=2.607935e-09 T=3.312596e+13 v=6.484619e+07 cgas=7.430339e+06 ceff=8.810315e+06 h=1.821876e+10 cfl*cell_dt_gas=7.562106e+01
DTREG rest  cycle=5350 npol=7 cfl*cell_dt=3.465143e+01 row=49 theta_deg=6.240244e+01 phi=4.571263e-01 r=2.499132e+12 dir=1 rho=1.289636e-12 T=6.717361e+14 v=1.146587e+08 cgas=3.345983e+07 ceff=3.363710e+07 h=1.712888e+10 cfl*cell_dt_gas=3.469290e+01
DTREP cycle=5350 time=4.769649e+05 cfl*dt_hyd=3.670292e+01 cfl*dt_gasonly=3.691322e+01 rank=0 gid=0 kji=77,52,59 dir=1 r=2.482061e+12 rho=4.334562e-12 T=3.350628e+13 pg=1.452350e+02 G1=1.666667e+00 cgas=7.472871e+06 v=1.307855e+08 ceff=8.265090e+06 h=1.701188e+10 tauT_h=2.573757e-02 Erad=4.117801e+0
DTREG polar cycle=5350 npol=7 cfl*cell_dt=7.399974e+01 row=4 theta_deg=4.658220e+01 phi=1.810097e-01 r=2.658148e+12 dir=1 rho=2.607935e-09 T=3.258120e+13 v=6.510499e+07 cgas=7.368989e+06 ceff=8.755105e+06 h=1.821876e+10 cfl*cell_dt_gas=7.541504e+01
DTREG rest  cycle=5350 npol=7 cfl*cell_dt=3.670292e+01 row=49 theta_deg=6.240244e+01 phi=4.571263e-01 r=2.482061e+12 dir=1 rho=4.334562e-12 T=3.350628e+13 v=1.307855e+08 cgas=7.472871e+06 ceff=8.265090e+06 h=1.701188e+10 cfl*cell_dt_gas=3.691322e+01
DTREP cycle=5400 time=4.784668e+05 cfl*dt_hyd=2.997701e+01 cfl*dt_gasonly=2.999759e+01 rank=0 gid=0 kji=76,50,58 dir=3 r=2.465107e+12 rho=2.484972e-12 T=2.140003e+15 pg=5.317849e+03 G1=1.666667e+00 cgas=5.972162e+07 v=-7.346667e+07 ceff=5.981306e+07 h=1.331776e+10 tauT_h=1.266798e-02 Erad=4.274293e+
DTREG polar cycle=5400 npol=7 cfl*cell_dt=7.097220e+01 row=5 theta_deg=4.693376e+01 phi=1.810097e-01 r=2.658148e+12 dir=1 rho=2.346263e-09 T=3.293530e+13 v=6.827000e+07 cgas=7.408925e+06 ceff=8.740831e+06 h=1.821876e+10 cfl*cell_dt_gas=7.222127e+01
DTREG rest  cycle=5400 npol=7 cfl*cell_dt=2.997701e+01 row=47 theta_deg=6.169932e+01 phi=4.509904e-01 r=2.465107e+12 dir=3 rho=2.484972e-12 T=2.140003e+15 v=-7.346667e+07 cgas=5.972162e+07 ceff=5.981306e+07 h=1.331776e+10 cfl*cell_dt_gas=2.999759e+01
DTREP cycle=5400 time=4.784668e+05 cfl*dt_hyd=2.724786e+01 cfl*dt_gasonly=2.735490e+01 rank=0 gid=0 kji=77,50,58 dir=3 r=2.465107e+12 rho=1.557265e-12 T=3.367822e+13 pg=5.244590e+01 G1=1.666667e+00 cgas=7.492020e+06 v=-1.385633e+08 ceff=8.065811e+06 h=1.331776e+10 tauT_h=7.180764e-03 Erad=4.203013e+
```
- compton report lines:
```
<rad_vet> compton cycle 5300 (cumulative): cell-steps 1.473292e+09 stopped at e_floor 1.917990e+05
<rad_vet> compton cycle 5350 (cumulative): cell-steps 4.822288e+09 stopped at e_floor 6.047580e+05
<rad_vet> compton cycle 5400 (cumulative): cell-steps 8.171665e+09 stopped at e_floor 9.744290e+05
<rad_vet> compton cycle 5450 (cumulative): cell-steps 1.152140e+10 stopped at e_floor 1.296353e+06
<rad_vet> compton cycle 5500 (cumulative): cell-steps 1.487152e+10 stopped at e_floor 1.572749e+06
<rad_vet> compton cycle 5550 (cumulative): cell-steps 1.822205e+10 stopped at e_floor 1.811430e+06
```
- m1-positivity line: ``
- hst tail (t dt L_bot L_top L_in Etot M_tot Lz_tot):
```
4.856980e+05 47.45 8.639533e+38 2.132709e+39 8.763705e+38 2.1540493934e+47 2.3310565725e+32 -6.9644871813e+38
4.857457e+05 47.672 8.640334e+38 2.132231e+39 8.763705e+38 2.1540487947e+47 2.3310565725e+32 -6.9644851468e+38
```

## C2 (2026-11-24 13:31:33 rc=0 wall=1241 s)
- CW 0, FATAL 0, created inf 0, DT_RETRY 0; s/cyc 2.9971 (c5350-5600)
- dt every 10 cycles: 13.94 13.94 7.8 7.431 7.143 9.851 8.324 9.129 6.802 7.499 8.101 6.482 7.422 6.921 6.846 9.695 7.927 7.563 7.188 6.99 7.27 9.787 9.079 8.814 12.75 9.205 6.86 11.93 9.589 12.23 11.88 8.645 9.894 9.071 12.22 9.414 
- DTREP / DTREG every 50 cycles (all in COMPTON/C2/dtrep.txt):
```
DTREP cycle=5300 time=4.754673e+05 cfl*dt_hyd=9.085302e+00 cfl*dt_gasonly=9.087862e+00 rank=0 gid=0 kji=78,51,59 dir=3 r=2.482061e+12 rho=3.729475e-13 T=2.561605e+15 pg=9.553441e+02 G1=1.666667e+00 cgas=6.534020e+07 v=-1.567154e+08 ceff=6.540277e+07 h=1.345340e+10 tauT_h=1.736638e-03 Erad=4.027447e+
DTREG polar cycle=5300 npol=7 cfl*cell_dt=3.863921e+01 row=3 theta_deg=4.623064e+01 phi=1.810097e-01 r=2.658148e+12 dir=1 rho=2.930684e-09 T=3.270285e+13 v=6.191644e+07 cgas=7.382733e+06 ceff=8.810006e+06 h=1.821876e+10 cfl*cell_dt_gas=3.943502e+01
DTREG rest  cycle=5300 npol=7 cfl*cell_dt=9.085302e+00 row=48 theta_deg=6.205088e+01 phi=4.632622e-01 r=2.482061e+12 dir=3 rho=3.729475e-13 T=2.561605e+15 v=-1.567154e+08 cgas=6.534020e+07 ceff=6.540277e+07 h=1.345340e+10 cfl*cell_dt_gas=9.087862e+00
DTREP cycle=5300 time=4.754673e+05 cfl*dt_hyd=1.063465e+01 cfl*dt_gasonly=1.064576e+01 rank=0 gid=0 kji=78,51,59 dir=3 r=2.482061e+12 rho=3.729475e-13 T=2.000680e+14 pg=7.461488e+01 G1=1.666667e+00 cgas=1.826052e+07 v=-1.712994e+08 ceff=1.845870e+07 h=1.345340e+10 tauT_h=1.637939e-03 Erad=4.027413e+
DTREG polar cycle=5300 npol=7 cfl*cell_dt=3.862201e+01 row=3 theta_deg=4.623064e+01 phi=1.810097e-01 r=2.658148e+12 dir=1 rho=2.930684e-09 T=3.259241e+13 v=6.195806e+07 cgas=7.370257e+06 ceff=8.799878e+06 h=1.821876e+10 cfl*cell_dt_gas=3.941844e+01
DTREG rest  cycle=5300 npol=7 cfl*cell_dt=1.063465e+01 row=48 theta_deg=6.205088e+01 phi=4.632622e-01 r=2.482061e+12 dir=3 rho=3.729475e-13 T=2.000680e+14 v=-1.712994e+08 cgas=1.826052e+07 ceff=1.845870e+07 h=1.345340e+10 cfl*cell_dt_gas=1.064576e+01
DTREP cycle=5350 time=4.758818e+05 cfl*dt_hyd=7.541028e+00 cfl*dt_gasonly=7.541923e+00 rank=0 gid=0 kji=78,51,59 dir=3 r=2.482061e+12 rho=1.979834e-14 T=7.330013e+15 pg=1.451221e+02 G1=1.666667e+00 cgas=1.105291e+08 v=-1.570433e+08 ceff=1.105609e+08 h=1.345340e+10 tauT_h=8.497707e-05 Erad=4.042782e+
DTREG polar cycle=5350 npol=7 cfl*cell_dt=3.867742e+01 row=4 theta_deg=4.658220e+01 phi=1.810097e-01 r=2.639991e+12 dir=1 rho=2.781135e-09 T=3.305596e+13 v=6.130860e+07 cgas=7.422484e+06 ceff=8.865341e+06 h=1.809432e+10 cfl*cell_dt_gas=3.948937e+01
DTREG rest  cycle=5350 npol=7 cfl*cell_dt=7.541028e+00 row=48 theta_deg=6.205088e+01 phi=4.632622e-01 r=2.482061e+12 dir=3 rho=1.979834e-14 T=7.330013e+15 v=-1.570433e+08 cgas=1.105291e+08 ceff=1.105609e+08 h=1.345340e+10 cfl*cell_dt_gas=7.541923e+00
DTREP cycle=5350 time=4.758818e+05 cfl*dt_hyd=9.078892e+00 cfl*dt_gasonly=9.083051e+00 rank=0 gid=0 kji=78,51,59 dir=3 r=2.482061e+12 rho=1.979834e-14 T=5.522491e+14 pg=1.093361e+01 G1=1.666667e+00 cgas=3.033835e+07 v=-1.918348e+08 ceff=3.044012e+07 h=1.345340e+10 tauT_h=7.979868e-05 Erad=4.042160e+
DTREG polar cycle=5350 npol=7 cfl*cell_dt=3.865217e+01 row=4 theta_deg=4.658220e+01 phi=1.810097e-01 r=2.639991e+12 dir=1 rho=2.781135e-09 T=3.292820e+13 v=6.136917e+07 cgas=7.408126e+06 ceff=8.850619e+06 h=1.809432e+10 cfl*cell_dt_gas=3.946284e+01
DTREG rest  cycle=5350 npol=7 cfl*cell_dt=9.078892e+00 row=48 theta_deg=6.205088e+01 phi=4.632622e-01 r=2.482061e+12 dir=3 rho=1.979834e-14 T=5.522491e+14 v=-1.918348e+08 cgas=3.033835e+07 ceff=3.044012e+07 h=1.345340e+10 cfl*cell_dt_gas=9.083051e+00
DTREP cycle=5400 time=4.762384e+05 cfl*dt_hyd=7.292980e+00 cfl*dt_gasonly=7.295841e+00 rank=0 gid=0 kji=77,50,59 dir=3 r=2.482061e+12 rho=7.959702e-13 T=9.011332e+14 pg=7.172752e+02 G1=1.666667e+00 cgas=3.875421e+07 v=-2.369375e+08 ceff=3.886234e+07 h=1.340935e+10 tauT_h=3.739227e-03 Erad=4.074579e+
DTREG polar cycle=5400 npol=7 cfl*cell_dt=3.808002e+01 row=4 theta_deg=4.658220e+01 phi=1.810097e-01 r=2.639991e+12 dir=1 rho=2.648876e-09 T=3.300128e+13 v=6.243496e+07 cgas=7.416342e+06 ceff=8.839879e+06 h=1.809432e+10 cfl*cell_dt_gas=3.885607e+01
DTREG rest  cycle=5400 npol=7 cfl*cell_dt=7.292980e+00 row=47 theta_deg=6.169932e+01 phi=4.571263e-01 r=2.482061e+12 dir=3 rho=7.959702e-13 T=9.011332e+14 v=-2.369375e+08 cgas=3.875421e+07 ceff=3.886234e+07 h=1.340935e+10 cfl*cell_dt_gas=7.295841e+00
DTREP cycle=5400 time=4.762384e+05 cfl*dt_hyd=7.552899e+00 cfl*dt_gasonly=7.567978e+00 rank=0 gid=0 kji=77,51,59 dir=3 r=2.482061e+12 rho=1.424192e-12 T=3.337851e+13 pg=4.753739e+01 G1=1.666667e+00 cgas=7.458609e+06 v=-2.591926e+08 ceff=7.990978e+06 h=1.345340e+10 tauT_h=6.639782e-03 Erad=4.056008e+
```
- compton report lines:
```
<rad_vet> compton cycle 5300 (cumulative): cell-steps 1.473229e+09 stopped at e_floor 1.943600e+05
<rad_vet> compton cycle 5350 (cumulative): cell-steps 4.822063e+09 stopped at e_floor 6.250100e+05
<rad_vet> compton cycle 5400 (cumulative): cell-steps 8.171010e+09 stopped at e_floor 1.044171e+06
<rad_vet> compton cycle 5450 (cumulative): cell-steps 1.152005e+10 stopped at e_floor 1.452650e+06
<rad_vet> compton cycle 5500 (cumulative): cell-steps 1.486919e+10 stopped at e_floor 1.849640e+06
<rad_vet> compton cycle 5550 (cumulative): cell-steps 1.821846e+10 stopped at e_floor 2.232656e+06
<rad_vet> compton cycle 5600 (cumulative): cell-steps 2.156786e+10 stopped at e_floor 2.599679e+06
```
- m1-positivity line: ``
- hst tail (t dt L_bot L_top L_in Etot M_tot Lz_tot):
```
4.785096e+05 8.9354 8.894211e+38 2.207246e+39 8.763705e+38 2.1541424007e+47 2.3310565725e+32 -6.9645095107e+38
4.785186e+05 9.075 8.894181e+38 2.207173e+39 8.763705e+38 2.1541422800e+47 2.3310565725e+32 -6.9645085505e+38
```

## C0b (2026-11-24 13:23:49 rc=0 wall=218 s)
- CW 0, FATAL 0, created inf 0, DT_RETRY 0; s/cyc 
- dt every 10 cycles: 13.94 24.97 14.78 17.84 13.96 
- DTREP / DTREG every 50 cycles (all in COMPTON/C0b/dtrep.txt):
```
DTREP cycle=5300 time=4.756428e+05 cfl*dt_hyd=1.782543e+01 cfl*dt_gasonly=1.782935e+01 rank=0 gid=0 kji=79,50,58 dir=3 r=2.465107e+12 rho=3.714662e-12 T=1.840681e+15 pg=6.837507e+03 G1=1.666667e+00 cgas=5.538773e+07 v=1.686993e+08 ceff=5.543703e+07 h=1.331776e+10 tauT_h=1.423448e-02 Erad=3.961774e+0
DTREG polar cycle=5300 npol=7 cfl*cell_dt=7.734394e+01 row=3 theta_deg=4.623064e+01 phi=1.810097e-01 r=2.658148e+12 dir=1 rho=2.946180e-09 T=3.283257e+13 v=6.184364e+07 cgas=7.397361e+06 ceff=8.822896e+06 h=1.821876e+10 cfl*cell_dt_gas=7.893630e+01
DTREG rest  cycle=5300 npol=7 cfl*cell_dt=1.782543e+01 row=47 theta_deg=6.169932e+01 phi=4.693981e-01 r=2.465107e+12 dir=3 rho=3.714662e-12 T=1.840681e+15 v=1.686993e+08 cgas=5.538773e+07 ceff=5.543703e+07 h=1.331776e+10 cfl*cell_dt_gas=1.782935e+01
DTREP cycle=5300 time=4.756428e+05 cfl*dt_hyd=1.804085e+01 cfl*dt_gasonly=1.804638e+01 rank=0 gid=0 kji=79,50,58 dir=3 r=2.465107e+12 rho=3.714662e-12 T=1.624963e+15 pg=6.036188e+03 G1=1.666667e+00 cgas=5.204106e+07 v=1.693511e+08 ceff=5.210899e+07 h=1.331776e+10 tauT_h=1.624155e-02 Erad=3.951490e+0
DTREG polar cycle=5300 npol=7 cfl*cell_dt=7.726963e+01 row=3 theta_deg=4.623064e+01 phi=1.810097e-01 r=2.658148e+12 dir=1 rho=2.946180e-09 T=3.258799e+13 v=6.193373e+07 cgas=7.369757e+06 ceff=8.800753e+06 h=1.821876e+10 cfl*cell_dt_gas=7.886512e+01
DTREG rest  cycle=5300 npol=7 cfl*cell_dt=1.804085e+01 row=47 theta_deg=6.169932e+01 phi=4.693981e-01 r=2.465107e+12 dir=3 rho=3.714662e-12 T=1.624963e+15 v=1.693511e+08 cgas=5.204106e+07 ceff=5.210899e+07 h=1.331776e+10 cfl*cell_dt_gas=1.804638e+01
```
- compton report lines:
```
```
- m1-positivity line: ``
- hst tail (t dt L_bot L_top L_in Etot M_tot Lz_tot):
```
4.760913e+05 14.092 8.852796e+38 2.231549e+39 8.763705e+38 2.1541748764e+47 2.3310565725e+32 -6.9645182258e+38
4.761053e+05 14.002 8.853186e+38 2.231405e+39 8.763705e+38 2.1541746867e+47 2.3310565725e+32 -6.9645187281e+38
```

## C1ra (2026-11-24 13:31:38 rc=0 wall=469 s)
- CW 0, FATAL 0, created inf 0, DT_RETRY 0; s/cyc 
- dt every 10 cycles: 13.94 26.45 15.74 16.75 17.13 22.61 28.92 42.66 23.91 31.23 
- DTREP / DTREG every 50 cycles (all in COMPTON/C1ra/dtrep.txt):
```
DTREP cycle=5300 time=4.757182e+05 cfl*dt_hyd=1.485293e+01 cfl*dt_gasonly=1.485487e+01 rank=0 gid=0 kji=78,51,59 dir=2 r=2.482061e+12 rho=5.460057e-13 T=6.372609e+15 pg=3.479481e+03 G1=1.666667e+00 cgas=1.030583e+08 v=2.045123e+08 ceff=1.030985e+08 h=1.522974e+10 tauT_h=2.623230e-03 Erad=4.038165e+0
DTREG polar cycle=5300 npol=7 cfl*cell_dt=7.739349e+01 row=3 theta_deg=4.623064e+01 phi=1.810097e-01 r=2.658148e+12 dir=1 rho=2.954438e-09 T=3.284690e+13 v=6.179664e+07 cgas=7.398976e+06 ceff=8.824652e+06 h=1.821876e+10 cfl*cell_dt_gas=7.898807e+01
DTREG rest  cycle=5300 npol=7 cfl*cell_dt=1.485293e+01 row=48 theta_deg=6.205088e+01 phi=4.632622e-01 r=2.482061e+12 dir=2 rho=5.460057e-13 T=6.372609e+15 v=2.045123e+08 cgas=1.030583e+08 ceff=1.030985e+08 h=1.522974e+10 cfl*cell_dt_gas=1.485487e+01
DTREP cycle=5300 time=4.757182e+05 cfl*dt_hyd=1.736667e+01 cfl*dt_gasonly=1.740613e+01 rank=0 gid=0 kji=78,50,58 dir=3 r=2.465107e+12 rho=1.244866e-12 T=3.334409e+13 pg=4.150893e+01 G1=1.666667e+00 cgas=7.454762e+06 v=2.220809e+08 ceff=7.976381e+06 h=1.331776e+10 tauT_h=5.741939e-03 Erad=4.039057e+0
DTREG polar cycle=5300 npol=7 cfl*cell_dt=7.731479e+01 row=3 theta_deg=4.623064e+01 phi=1.810097e-01 r=2.658148e+12 dir=1 rho=2.954438e-09 T=3.258579e+13 v=6.189193e+07 cgas=7.369508e+06 ceff=8.801243e+06 h=1.821876e+10 cfl*cell_dt_gas=7.891300e+01
DTREG rest  cycle=5300 npol=7 cfl*cell_dt=1.736667e+01 row=47 theta_deg=6.169932e+01 phi=4.632622e-01 r=2.465107e+12 dir=3 rho=1.244866e-12 T=3.334409e+13 v=2.220809e+08 cgas=7.454762e+06 ceff=7.976381e+06 h=1.331776e+10 cfl*cell_dt_gas=1.740613e+01
DTREP cycle=5350 time=4.769649e+05 cfl*dt_hyd=3.465143e+01 cfl*dt_gasonly=3.469290e+01 rank=0 gid=0 kji=77,52,60 dir=1 r=2.499132e+12 rho=1.289636e-12 T=6.717361e+14 pg=8.662948e+02 G1=1.666667e+00 cgas=3.345983e+07 v=1.146587e+08 ceff=3.363710e+07 h=1.712888e+10 tauT_h=7.509563e-03 Erad=4.067619e+0
DTREG polar cycle=5350 npol=7 cfl*cell_dt=7.420428e+01 row=4 theta_deg=4.658220e+01 phi=1.810097e-01 r=2.658148e+12 dir=1 rho=2.607935e-09 T=3.312596e+13 v=6.484619e+07 cgas=7.430339e+06 ceff=8.810315e+06 h=1.821876e+10 cfl*cell_dt_gas=7.562106e+01
DTREG rest  cycle=5350 npol=7 cfl*cell_dt=3.465143e+01 row=49 theta_deg=6.240244e+01 phi=4.571263e-01 r=2.499132e+12 dir=1 rho=1.289636e-12 T=6.717361e+14 v=1.146587e+08 cgas=3.345983e+07 ceff=3.363710e+07 h=1.712888e+10 cfl*cell_dt_gas=3.469290e+01
DTREP cycle=5350 time=4.769649e+05 cfl*dt_hyd=3.670292e+01 cfl*dt_gasonly=3.691322e+01 rank=0 gid=0 kji=77,52,59 dir=1 r=2.482061e+12 rho=4.334562e-12 T=3.350628e+13 pg=1.452350e+02 G1=1.666667e+00 cgas=7.472871e+06 v=1.307855e+08 ceff=8.265090e+06 h=1.701188e+10 tauT_h=2.573757e-02 Erad=4.117801e+0
DTREG polar cycle=5350 npol=7 cfl*cell_dt=7.399974e+01 row=4 theta_deg=4.658220e+01 phi=1.810097e-01 r=2.658148e+12 dir=1 rho=2.607935e-09 T=3.258120e+13 v=6.510499e+07 cgas=7.368989e+06 ceff=8.755105e+06 h=1.821876e+10 cfl*cell_dt_gas=7.541504e+01
DTREG rest  cycle=5350 npol=7 cfl*cell_dt=3.670292e+01 row=49 theta_deg=6.240244e+01 phi=4.571263e-01 r=2.482061e+12 dir=1 rho=4.334562e-12 T=3.350628e+13 v=1.307855e+08 cgas=7.472871e+06 ceff=8.265090e+06 h=1.701188e+10 cfl*cell_dt_gas=3.691322e+01
```
- compton report lines:
```
<rad_vet> compton cycle 5300 (cumulative): cell-steps 1.473292e+09 stopped at e_floor 1.917990e+05
<rad_vet> compton cycle 5350 (cumulative): cell-steps 4.822288e+09 stopped at e_floor 6.047580e+05
```
- m1-positivity line: ``
- hst tail (t dt L_bot L_top L_in Etot M_tot Lz_tot):
```
4.777518e+05 25.448 8.888902e+38 2.215284e+39 8.763705e+38 2.1541525148e+47 2.3310565725e+32 -6.9645119926e+38
4.777771e+05 25.315 8.889273e+38 2.215026e+39 8.763705e+38 2.1541521759e+47 2.3310565725e+32 -6.9645124437e+38
```

## C1rb (2026-11-24 13:38:34 rc=0 wall=416 s)
- CW 0, FATAL 0, created inf 0, DT_RETRY 0; s/cyc 
- dt every 10 cycles: 25.37 42.34 33.09 19.26 30.57 43.1 42.8 34.96 42.63 
- DTREP / DTREG every 50 cycles (all in COMPTON/C1rb/dtrep.txt):
```
DTREP cycle=5400 time=4.784668e+05 cfl*dt_hyd=2.997701e+01 cfl*dt_gasonly=2.999759e+01 rank=0 gid=0 kji=76,50,58 dir=3 r=2.465107e+12 rho=2.484972e-12 T=2.140003e+15 pg=5.317849e+03 G1=1.666667e+00 cgas=5.972162e+07 v=-7.346667e+07 ceff=5.981306e+07 h=1.331776e+10 tauT_h=1.266798e-02 Erad=4.274293e+
DTREG polar cycle=5400 npol=7 cfl*cell_dt=7.097220e+01 row=5 theta_deg=4.693376e+01 phi=1.810097e-01 r=2.658148e+12 dir=1 rho=2.346263e-09 T=3.293530e+13 v=6.827000e+07 cgas=7.408925e+06 ceff=8.740831e+06 h=1.821876e+10 cfl*cell_dt_gas=7.222127e+01
DTREG rest  cycle=5400 npol=7 cfl*cell_dt=2.997701e+01 row=47 theta_deg=6.169932e+01 phi=4.509904e-01 r=2.465107e+12 dir=3 rho=2.484972e-12 T=2.140003e+15 v=-7.346667e+07 cgas=5.972162e+07 ceff=5.981306e+07 h=1.331776e+10 cfl*cell_dt_gas=2.999759e+01
DTREP cycle=5400 time=4.784668e+05 cfl*dt_hyd=2.724786e+01 cfl*dt_gasonly=2.735490e+01 rank=0 gid=0 kji=77,50,58 dir=3 r=2.465107e+12 rho=1.557265e-12 T=3.367822e+13 pg=5.244590e+01 G1=1.666667e+00 cgas=7.492020e+06 v=-1.385633e+08 ceff=8.065811e+06 h=1.331776e+10 tauT_h=7.180764e-03 Erad=4.203013e+
DTREG polar cycle=5400 npol=7 cfl*cell_dt=7.080633e+01 row=5 theta_deg=4.693376e+01 phi=1.810097e-01 r=2.658148e+12 dir=1 rho=2.346263e-09 T=3.247121e+13 v=6.849944e+07 cgas=7.356540e+06 ceff=8.691799e+06 h=1.821876e+10 cfl*cell_dt_gas=7.205270e+01
DTREG rest  cycle=5400 npol=7 cfl*cell_dt=2.724786e+01 row=47 theta_deg=6.169932e+01 phi=4.571263e-01 r=2.465107e+12 dir=3 rho=1.557265e-12 T=3.367822e+13 v=-1.385633e+08 cgas=7.492020e+06 ceff=8.065811e+06 h=1.331776e+10 cfl*cell_dt_gas=2.735490e+01
DTREG polar cycle=5450 npol=7 cfl*cell_dt=7.030252e+01 row=6 theta_deg=4.728532e+01 phi=1.810097e-01 r=2.658148e+12 dir=1 rho=2.320980e-09 T=3.306138e+13 v=6.901663e+07 cgas=7.423092e+06 ceff=8.727776e+06 h=1.821876e+10 cfl*cell_dt_gas=7.150246e+01
DTREP cycle=5450 time=480121 cfl*dt_hyd=40.4955 cfl*dt_gasonly=40.6272 rank=1 gid=2 kji=75,122,54 dir=2 r=2.39844e+12 rho=2.82409e-12 T=3.32001e+14 pg=937.602 G1=1.66667 cgas=2.35231e+07 v=8.51478e+07 ceff=2.38767e+07 h=1.47167e+10 tauT_h=0.0124741 Erad=6.2039e+07 Ddd=0.335925 cdt_kP=1.55956 frz=-1 
DTREG rest  cycle=5450 npol=7 cfl*cell_dt=40.4955 row=119 theta_deg=87.0117 phi=1.23025 r=2.39844e+12 dir=2 rho=2.82409e-12 T=3.32001e+14 v=8.51478e+07 cgas=2.35231e+07 ceff=2.38767e+07 h=1.47167e+10 cfl*cell_dt_gas=40.6272
DTREP cycle=5450 time=4.801212e+05 cfl*dt_hyd=4.330079e+01 cfl*dt_gasonly=4.344853e+01 rank=0 gid=1 kji=14,5,84 dir=1 r=2.945974e+12 rho=1.859152e-10 T=2.712838e+13 pg=5.043577e+03 G1=1.666667e+00 cgas=6.724133e+06 v=1.326926e+08 ceff=7.199836e+06 h=2.019151e+10 tauT_h=1.823591e+00 Erad=1.769499e+07
DTREG polar cycle=5450 npol=7 cfl*cell_dt=7.009950e+01 row=6 theta_deg=4.728532e+01 phi=1.810097e-01 r=2.658148e+12 dir=1 rho=2.320980e-09 T=3.229177e+13 v=6.932660e+07 cgas=7.336185e+06 ceff=8.642974e+06 h=1.821876e+10 cfl*cell_dt_gas=7.129441e+01
DTREG rest  cycle=5450 npol=7 cfl*cell_dt=4.330079e+01 row=130 theta_deg=9.087890e+01 phi=7.056312e-02 r=2.945974e+12 dir=1 rho=1.859152e-10 T=2.712838e+13 v=1.326926e+08 cgas=6.724133e+06 ceff=7.199836e+06 h=2.019151e+10 cfl*cell_dt_gas=4.344853e+01
```
- compton report lines:
```
<rad_vet> compton cycle 5400 (cumulative): cell-steps 1.473772e+09 stopped at e_floor 1.576800e+05
<rad_vet> compton cycle 5450 (cumulative): cell-steps 4.823504e+09 stopped at e_floor 4.796040e+05
```
- m1-positivity line: ``
- hst tail (t dt L_bot L_top L_in Etot M_tot Lz_tot):
```
4.808010e+05 42.261 8.793523e+38 2.184684e+39 8.763705e+38 2.1541121638e+47 2.3310565725e+32 -6.9645027982e+38
4.808433e+05 42.211 8.790464e+38 2.184231e+39 8.763705e+38 2.1541116117e+47 2.3310565725e+32 -6.9645033903e+38
```

## Comparisons
- C0b (compton binary, key off, 50 cycles) vs C0 (T2 binary), user.hst first 50 rows: BITWISE
- C1r (C1ra 0-100 + restart + C1rb, 88 rows: stopped early on zzzw) vs C1: first 100 BITWISE; C1rb rows 1-88 vs C1 rows 101-188 BITWISE
- C1 vs C2 vs C0 at equal model time:
```
common t 4.752244e+05..4.785186e+05
L_top/L_in C1-C2 rms 2.047e-04 max 5.209e-04 (rel 8.067e-05 / 2.053e-04)
L_top/L_in C0-C2 rms 1.922e-04 max 5.787e-04 (rel 7.575e-05 / 2.281e-04)
L_top/L_in C1-C0 rms 8.386e-05 max 4.411e-04 (rel 3.305e-05 / 1.738e-04)
KE_int     C1-C2 rms 1.051e+38 max 2.298e+38 (rel 6.017e-07 / 1.316e-06)
KE_int     C0-C2 rms 9.190e+38 max 1.493e+40 (rel 5.263e-06 / 8.549e-05)
KE_int     C1-C0 rms 9.112e+38 max 1.485e+40 (rel 5.218e-06 / 8.505e-05)
```
- Drift (Etot col 19, M_tot col 16, Lz_tot col 28): C0 dE/E -6.007e-05 dM/M 3.091e-16 dLz/Lz -4.710e-06; C1 dE/E -6.400e-05 dM/M -2.628e-15 dLz/Lz -5.374e-06; C2 dE/E -2.060e-05 dM/M 1.546e-15 dLz/Lz -2.014e-06; 
