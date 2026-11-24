# NOTE 2026-11-24zzf caltech -> orion: T9HOT2 bisect B / F / K / W (24zzc) and T9HOT P1 / P2 (24zy, trimmed by 24zze)

- Binary athena_t9hot_e8eaa95af (md5 fbc14236..., gate PASS 3.686e-11). B/F/K/W from rst_keep/t9l62.00013.rst (cycle 5121), 300 cycles, hydro_w + vet bin every 50 cycles (single precision), stats + B sub-boxes per dump (dumps deleted after).
  P1 / P2 from t9l62's final rst (final_1124). All keys named: time/cfl_number=0.3 hydro/fofc=true hydro/fofc_retest=true rad_vet/clean_stag_mg_levels=4 rad_vet/floorcap_v=3e7 rad_vet/floorcap_rho=1e-14 rad_vet/floorcap_report=25 hydro/dt_report=10 hydro/rad_signal_speed_phase=true + the arm key.
- Stop rule: first CW / FATAL -> scancel (B/F/K/W per 24zzc). The hydro_fofc count variable is not an output variable of this deck/binary, so it is not posted.
- Files: caltech_results/T9HOT2/<arm>/ (stats.txt, dtrep.txt, hst, log.0.gz, run.log; B: box_c<cycle>.npz).

## W (11-24 10:49)
- wellbalance_dynamic=false nlim=5421 ; 2026-11-24 10:49:41 rc=0 wall=51 s 
- CW 0, FATAL 4, created inf 0, DT_RETRY 0; s/cyc n/a
- dt every 10 cycles: 
- DTREP every 50 cycles:
```
```
- hst tail (t dt L_bot L_top L_in Etot M_tot Lz_tot):
```
```
- **Arm W FATAL at startup (10:49), as designed by the pgen, not a run failure:** `he_star_vet: needs <hydro>/etotgrav = true, wellbalance_dynamic = true and wb_x1 = true (the gravity source reads the x1 well-balanced cache)` (he_star_vet.cpp:656). So hydro/wellbalance_dynamic=false cannot be bisected on this binary. **W is not rerun**; tell me if you want another WB-off form. Its node 02-17 goes to T9HOT P1 (from t9l62.00020).

## F (11-24 10:58)
- fofc=false nlim=5421 ; 2026-11-24 10:58:09 STOP at first CW/FATAL: CLEAN WARNING cycle 5278 t2s 2: stage solve failed, the step is redone with backward Euler (stage fallback 1) 2026-11-24 10:58:11 rc=0 wall=573 s 
- CW 1, FATAL 0, created inf 0, DT_RETRY 0; s/cyc 2.9761 (c5200-5250)
- dt every 10 cycles: 33.12 32.72 35.98 43.79 30.87 18.36 9.828 14.25 22.29 39.4 46.51 46.07 33.23 17.31 12.95 10.95 6.696 
- per-dump stats:
```
cycle 5150: N_hot 0 | thin |v| max p99 1.118e+08 2.934e+07 | dense |v| max p99 1.126e+08 3.308e+07 | thin sum eint 1.569e+06 KE 1.707e+08 (cell sums, no volume)
cycle 5200: N_hot 5 | thin |v| max p99 1.552e+08 2.925e+07 | dense |v| max p99 1.227e+08 3.313e+07 | thin sum eint 1.638e+06 KE 1.736e+08 (cell sums, no volume)
cycle 5250: N_hot 0 | thin |v| max p99 2.601e+08 2.909e+07 | dense |v| max p99 1.231e+08 3.320e+07 | thin sum eint 1.592e+06 KE 1.782e+08 (cell sums, no volume)
```
- DTREP every 50 cycles:
```
c5150 dt_hyd 40.8 gasonly 40.89 dir 3 gid 0 kji 79,49,58 r 2.4651e+12 rho 1.991e-12 T 3.099e+14 v -7.466e+07 ceff 2.293e+07 tauT_h 0.00998
c5200 dt_hyd 18.95 gasonly 18.97 dir 3 gid 0 kji 78,49,59 r 2.4821e+12 rho 1.688e-11 T 1.077e+15 v -1.690e+08 ceff 4.250e+07 tauT_h 0.0761
c5250 dt_hyd 18.18 gasonly 18.19 dir 3 gid 0 kji 77,50,59 r 2.4821e+12 rho 1.755e-13 T 1.052e+15 v -1.792e+08 ceff 4.201e+07 tauT_h 0.000973
```
- hst tail (t dt L_bot L_top L_in Etot M_tot Lz_tot):
```
4.770823e+05 6.6959 2.557970e+38 2.221271e+39 8.763705e+38 2.1541614981e+47 2.3310565725e+32 2.2929342399e+38
4.770894e+05 7.1478 2.551334e+38 2.221193e+39 8.763705e+38 2.1541614020e+47 2.3310565725e+32 2.3439229339e+38
```

## K (11-24 11:05)
- floor_from_ekin=false nlim=5421 ; 2026-11-24 11:04:51 rc=0 wall=973 s 
- CW 0, FATAL 0, created inf 0, DT_RETRY 0; s/cyc 2.9972 (c5200-5400)
- dt every 10 cycles: 33.12 32.72 35.98 43.77 30.86 44.94 21.55 37.6 15.04 17.49 32.17 38.2 28.6 27.91 26.64 24.42 19.85 28.58 29.71 41.62 35.82 43.07 41.57 31.51 23.01 15.35 24.54 22.22 32.23 16.07 
- per-dump stats:
```
cycle 5150: N_hot 0 | thin |v| max p99 1.118e+08 2.929e+07 | dense |v| max p99 1.126e+08 3.308e+07 | thin sum eint 1.569e+06 KE 1.707e+08 (cell sums, no volume)
cycle 5200: N_hot 2 | thin |v| max p99 3.000e+08 2.913e+07 | dense |v| max p99 1.219e+08 3.315e+07 | thin sum eint 1.597e+06 KE 1.761e+08 (cell sums, no volume)
cycle 5250: N_hot 2 | thin |v| max p99 1.608e+08 2.898e+07 | dense |v| max p99 1.232e+08 3.320e+07 | thin sum eint 1.594e+06 KE 1.784e+08 (cell sums, no volume)
cycle 5300: N_hot 0 | thin |v| max p99 1.252e+08 2.884e+07 | dense |v| max p99 1.363e+08 3.324e+07 | thin sum eint 1.612e+06 KE 1.829e+08 (cell sums, no volume)
cycle 5350: N_hot 0 | thin |v| max p99 1.348e+08 2.863e+07 | dense |v| max p99 1.450e+08 3.331e+07 | thin sum eint 1.625e+06 KE 1.878e+08 (cell sums, no volume)
cycle 5400: N_hot 1 | thin |v| max p99 1.115e+08 2.852e+07 | dense |v| max p99 1.436e+08 3.334e+07 | thin sum eint 1.622e+06 KE 1.892e+08 (cell sums, no volume)
cycle 5421: N_hot 2 | thin |v| max p99 2.024e+08 2.848e+07 | dense |v| max p99 1.431e+08 3.336e+07 | thin sum eint 1.639e+06 KE 1.899e+08 (cell sums, no volume)
```
- DTREP every 50 cycles:
```
c5150 dt_hyd 40.78 gasonly 40.86 dir 3 gid 0 kji 79,49,58 r 2.4651e+12 rho 1.987e-12 T 3.107e+14 v -7.469e+07 ceff 2.296e+07 tauT_h 0.00996
c5200 dt_hyd 14.34 gasonly 14.34 dir 1 gid 0 kji 78,51,60 r 2.4991e+12 rho 3.385e-12 T 1.160e+16 v 2.193e+08 ceff 1.391e+08 tauT_h 0.0241
c5250 dt_hyd 27.51 gasonly 27.56 dir 3 gid 0 kji 78,50,59 r 2.4821e+12 rho 1.950e-12 T 2.000e+14 v -1.277e+08 ceff 1.852e+07 tauT_h 0.00984
c5300 dt_hyd 29.63 gasonly 29.74 dir 3 gid 0 kji 78,50,58 r 2.4651e+12 rho 6.776e-12 T 1.567e+14 v 1.182e+08 ceff 1.668e+07 tauT_h 0.0435
c5350 dt_hyd 29.54 gasonly 29.64 dir 2 gid 2 kji 75,122,54 r 2.3984e+12 rho 2.076e-12 T 4.422e+14 v 1.218e+08 ceff 2.763e+07 tauT_h 0.0116
c5400 dt_hyd 32.72 gasonly 32.77 dir 3 gid 0 kji 77,50,58 r 2.4651e+12 rho 2.994e-11 T 3.001e+14 v 9.955e+07 ceff 2.256e+07 tauT_h 0.135
```
- hst tail (t dt L_bot L_top L_in Etot M_tot Lz_tot):
```
4.814911e+05 17.294 8.644341e+38 2.176744e+39 8.763705e+38 2.1541031667e+47 2.3310565725e+32 -6.4520524272e+41
4.815082e+05 17.072 8.694694e+38 2.176539e+39 8.763705e+38 2.1541029447e+47 2.3310565725e+32 -6.4517163939e+41
```

## B (11-24 11:05)
- nlim=5421 ; 2026-11-24 11:04:48 rc=0 wall=970 s 
- CW 0, FATAL 0, created inf 0, DT_RETRY 0; s/cyc 2.9981 (c5200-5400)
- dt every 10 cycles: 33.12 32.72 35.98 43.77 30.86 44.84 23.6 32 30.15 21.63 18.76 21.61 16.43 16.8 33.53 36.09 18.31 19.75 31.14 36.46 23.98 40.07 36 43.08 42.72 31.53 23.61 18.89 37.27 25.84 
- per-dump stats:
```
cycle 5150: N_hot 0 | thin |v| max p99 1.118e+08 2.929e+07 | dense |v| max p99 1.126e+08 3.308e+07 | thin sum eint 1.569e+06 KE 1.707e+08 (cell sums, no volume)
cycle 5200: N_hot 3 | thin |v| max p99 1.309e+08 2.913e+07 | dense |v| max p99 1.218e+08 3.315e+07 | thin sum eint 1.605e+06 KE 1.758e+08 (cell sums, no volume)
cycle 5250: N_hot 2 | thin |v| max p99 1.775e+08 2.904e+07 | dense |v| max p99 1.233e+08 3.318e+07 | thin sum eint 1.598e+06 KE 1.778e+08 (cell sums, no volume)
cycle 5300: N_hot 0 | thin |v| max p99 1.666e+08 2.890e+07 | dense |v| max p99 1.294e+08 3.322e+07 | thin sum eint 1.603e+06 KE 1.813e+08 (cell sums, no volume)
cycle 5350: N_hot 0 | thin |v| max p99 1.074e+08 2.872e+07 | dense |v| max p99 1.442e+08 3.328e+07 | thin sum eint 1.604e+06 KE 1.839e+08 (cell sums, no volume)
cycle 5400: N_hot 0 | thin |v| max p99 1.506e+08 2.856e+07 | dense |v| max p99 1.419e+08 3.333e+07 | thin sum eint 1.625e+06 KE 1.891e+08 (cell sums, no volume)
cycle 5421: N_hot 0 | thin |v| max p99 1.351e+08 2.849e+07 | dense |v| max p99 1.436e+08 3.335e+07 | thin sum eint 1.620e+06 KE 1.894e+08 (cell sums, no volume)
```
- DTREP every 50 cycles:
```
c5150 dt_hyd 40.78 gasonly 40.86 dir 3 gid 0 kji 79,49,58 r 2.4651e+12 rho 1.987e-12 T 3.107e+14 v -7.469e+07 ceff 2.296e+07 tauT_h 0.00996
c5200 dt_hyd 27.28 gasonly 27.28 dir 1 gid 0 kji 78,51,59 r 2.4821e+12 rho 3.006e-13 T 6.782e+15 v 8.076e+07 ceff 1.063e+08 tauT_h 0.00124
c5250 dt_hyd 16.83 gasonly 16.83 dir 3 gid 0 kji 77,51,59 r 2.4821e+12 rho 3.462e-14 T 5.202e+15 v -1.467e+08 ceff 9.316e+07 tauT_h 0.000158
c5300 dt_hyd 29.6 gasonly 29.62 dir 1 gid 0 kji 76,51,58 r 2.4651e+12 rho 2.684e-12 T 1.842e+15 v -1.157e+08 ceff 5.549e+07 tauT_h 0.0134
c5350 dt_hyd 42.11 gasonly 42.27 dir 2 gid 2 kji 75,120,53 r 2.3821e+12 rho 1.908e-11 T 2.209e+14 v 8.454e+07 ceff 1.959e+07 tauT_h 0.0827
c5400 dt_hyd 30.82 gasonly 30.86 dir 1 gid 0 kji 78,50,59 r 2.4821e+12 rho 1.067e-12 T 9.496e+14 v 1.256e+08 ceff 3.999e+07 tauT_h 0.00717
```
- hst tail (t dt L_bot L_top L_in Etot M_tot Lz_tot):
```
4.813213e+05 18.106 8.135956e+38 2.178756e+39 8.763705e+38 2.1541053764e+47 2.3310565725e+32 5.7439052673e+39
4.813373e+05 16.06 8.184581e+38 2.178559e+39 8.763705e+38 2.1541051673e+47 2.3310565725e+32 5.7743372314e+39
```
