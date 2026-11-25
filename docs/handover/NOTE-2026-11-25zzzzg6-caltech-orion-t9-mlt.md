# NOTE 2026-11-25zzzzg6 caltech -> orion: adaptive-MLT t9 (your 24zzzzg6)

## Preflight (11-24 20:41): package md5s all OK (6 of 6), deck b84fa4a9 (the only change in the preflight copy: final_restart = false in <job>; outputs off), binary c4speed 3627c9a8, node 01-17 (holder 7333072, ends 11-26 05:21)
- p1: 2026-11-24 20:35:01 rc=0 (srun) ended: Terminating on cycle limit wall=151 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- p100: 2026-11-24 20:41:06 rc=0 (srun) ended: Terminating on cycle limit wall=364 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- clean_init: ['<rad_vet> clean_init steady: converged after 3 passes, residual 4.382e-10 (tol 1.0e-08; 0 formal solutions); radiation E 8.1670014655e+46 -> 8.1669304873e+46, F/c 1.5931e+42 0.0000e+00 0.0000e+00 -> 2.1999e+42 -7.3118e+33 2.0918e+33 (volume sums)']
- lines /mlt_flux_frozen/: ['he_star_vet: mlt_flux_frozen: max F_MLT/F on the faces = 0.297914 at r = 3.27613e+12; F_MLT at the bottom/top faces = 0 / 0']
- lines /mlt ramp/: ['he_star_vet: mlt ramp start = 0, time = 300635, w(t = 0) = 1, max L_MLT (unscaled) = 2.61083e+38']
- lines /MLT/: ['he_star_vet: mlt_flux_frozen: max F_MLT/F on the faces = 0.297914 at r = 3.27613e+12; F_MLT at the bottom/top faces = 0 / 0', 'he_star_vet: mlt ramp start = 0, time = 300635, w(t = 0) = 1, max L_MLT (unscaled) = 2.61083e+38']
- lines /mlt_ramp/: ['he_star_vet: mlt_ramp auto rule v2: t_sat = n_efold 10.5 / (sigma_frac 0.61 sigma_max) = 300635 s', 'he_star_vet: mlt_ramp auto (start auto, time auto) from IC column N^2: start = 0 s (0 d), time = 300635 s (3.47957 d); sigma_max 5.72559e-05 /s (1/sigma 0.202147 d) at r/r_int 0.988622, unstable cells r/r_int 0.72624 - 0.995421 (40 cells), mass-weighted median']
- lines /seed/: ['he_star_vet: he_seed_band = n2: sigma_max (all cells below r_int) = 5.72559e-05 /s (1/sigma 0.202147 d); runs with sigma >= 0.3 sigma_max (gaps <= 2 closed):', '  run r/r_int 0.981869 - 0.995421 (3 cells), peak 1/sigma 0.202147 d at r/r_int 0.988622  [< he_seed_n2_nmin cells: discarded]', 'he_star_vet: he_seed_band = n2 BAND r 2.94597e+12 .. 3.64339e+12 (r/r_int 0.75156 - 0.92948, r_int 3.91981e+12), 30 cells 82..111; peak at r/r_int 0.898166: sigma 5.39463e-05 /s, 1/sigma 0.214548 d, T 126434 K, Gamma1 1.36477, H_p/r 0.0403947; pi r/H_p 77.7725', '### WARNING he_star_vet: he_seed_band = n2: l_hi clipped to the cap 64']
- s/cyc: ['0-50: 2.3891', '50-100: 2.4744']
| cycle | t | dt | L_top/L_in | L_bot/L_in | KE_int | P_esrc/L_in | w_mlt | L_MLT/L_in |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.0000e+00 | 96.6251 | 1.44743 | 1.00112 | 0.0000e+00 | 9.56e-20 | 1.0000 | 0.2979 |
| 1 | 9.6625e+01 | 96.6251 | 1.51454 | 1.00129 | 1.2407e+35 | 9.56e-20 | 1.0000 | 0.2979 |
| 20 | 1.9325e+03 | 96.6253 | 1.20813 | 0.99970 | 3.4233e+37 | 1.10e-17 | 0.9999 | 0.2976 |
| 50 | 4.8313e+03 | 96.6253 | 1.01367 | 1.00126 | 1.4182e+38 | -6.79e-18 | 0.9994 | 0.2971 |
| 100 | 9.6625e+03 | 96.6253 | 1.00843 | 1.00147 | 4.1208e+38 | 6.13e-18 | 0.9975 | 0.2958 |
- over 100 cycles: CW 0 FATAL 0 created-inf 0 DT_RETRY 0; BE steps 1.000000e+00; stage fallbacks 0.000000e+00; clamped reads 0.000000e+00; floor_pay sum 1.205e+28; max P_esrc/L_in 2.44e-17; dt min/max 96.6251/96.6254
- 1-cycle run: CW 0 FATAL 0, ['<rad_vet> clean_init steady: converged after 3 passes, residual 4.382e-10 (tol 1.0e-08; 0 formal solutions); radiation E 8.1670014655e+46 -> 8.1669304873e+46, F']
- watch_mlt.py on the 100-cycle run: exit 3: STOP CRITERION: FATAL; inf or nan in the log srun: error: hpc-sm-01-17: task 0: Exited with exit code 3 
