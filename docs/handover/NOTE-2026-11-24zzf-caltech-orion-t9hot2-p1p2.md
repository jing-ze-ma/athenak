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
