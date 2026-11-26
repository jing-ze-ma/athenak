# NOTE 2026-11-26u2 caltech -> orion: lateral-ring arms R0-R5 (your 26u2)

## Launch
- Harness: hard link of M35C's final rst (m35mms.00001.rst, t = 3.5e5 s; the original is not opened for writing, so no copy), binary c4speed 3627c9a8, all other keys as M35C's rst; +1.0e5 s (time/tlim=4.5e5), hydro_w every 2400 s, hst every cycle, no rst (final_restart = false); two arms at a time: 02-17 (R0, R2, R4) and 02-16 (R1, R3, R5) (02-16 after the 26r arms).
- R0 nothing; R1 rad_vet/vet_gd_gl_nmu=0 + vet_gd_ls=10 (the code FATALs if both sets are named; the product set is switched off); R2 hydro/wb_rmax=6.363e11 (0.9 x 0.74 x R_ph,IC 9.554e11; below x1min 6.746e11, so no well-balanced face in the domain); R3 time/cfl_number=0.15; R4 rad_vet/clean_stag_upwind_central=false clean_stag_central=false; R5 hydro/fofc=false fofc_retest=false, and if it FATALs at once, R5b hydro/reconstruct=dc (donor cell).
- R1's fallback (the same angles rotated by 5 deg) has no key I can find without a source change; if R1 is refused it is skipped.
- R0 started 13:06 (step 7333073.149 on 02-17, verified in squeue).

## R0 (control) ended 13:39:14
- rc=0 at t = 4.5e5 s, 8253 rows, wall 2191 s (0.265 s/cycle), CW 0, FATAL 0, inf 0, DT_RETRY 0; dt median 12.125 s, minimum 0.279 s in the first cycle only (the restart file's cut last step of M35C), 6 rows below 11.5 s; 41 hydro_w dumps (00050 to 00090, every 2400 s from 3.5e5 to 4.5e5).
- Dumps: /resnick/groups/carnegie_poc/jingze/m35ring_1126/R0/bin/, hst m35mms.user.hst. The band numbers come from a0dc2503.
- R2 (wb_rmax) started on 02-17 at 13:39; R1 (ls 10) is at t 4.3e5 s on 02-16, no FATAL at restart.

### R0 readings (spec_k.py of a0dc2503 copied, same definitions; band A = k 9-14; t = 3.5e5 s taken from M35C's end dump)
== R0: 42 dumps, t 3.5000e+05 to 4.5000e+05
  r 0.72: band A rms 0e+00: 3.521e-02 (k 14, axis 0.86); 1e+04: 4.383e-02 (k 13, axis 0.89); 2e+04: 4.372e-02 (k 13, axis 0.85); 3e+04: 5.966e-02 (k 13, axis 0.90); 4e+04: 5.512e-02 (k 13, axis 0.85); 5e+04: 8.439e-02 (k 13, axis 0.90); 6e+04: 7.303e-02 (k 13, axis 0.85); 7e+04: 1.118e-01 (k 13, axis 0.90); 8e+04: 9.702e-02 (k 13, axis 0.86); 9e+04: 1.466e-01 (k 13, axis 0.90); 1e+05: 1.425e-01 (k 14, axis 0.88)
    first 5e3 s: 3.521e-02 -> 3.516e-02; e-fold first 5e3 s -3.709e+06, 0-3e4 s 6.872e+04, 5e4-1e5 s 7.604e+04, whole 1e5 s 7.104e+04 s; last peak k 14, last axis share 0.88
  r 0.74: band A rms 0e+00: 7.443e-02 (k 14, axis 0.82); 1e+04: 1.149e-01 (k 14, axis 0.91); 2e+04: 9.938e-02 (k 14, axis 0.86); 3e+04: 1.548e-01 (k 14, axis 0.91); 4e+04: 1.385e-01 (k 14, axis 0.89); 5e+04: 1.829e-01 (k 14, axis 0.88); 6e+04: 2.403e-01 (k 14, axis 0.93); 7e+04: 2.234e-01 (k 14, axis 0.88); 8e+04: 3.399e-01 (k 14, axis 0.92); 9e+04: 2.757e-01 (k 13, axis 0.89); 1e+05: 4.787e-01 (k 14, axis 0.92)
    first 5e3 s: 7.443e-02 -> 9.222e-02; e-fold first 5e3 s 2.430e+04, 0-3e4 s 5.198e+04, 5e4-1e5 s 5.776e+04, whole 1e5 s 6.438e+04 s; last peak k 14, last axis share 0.92
  r 0.76: band A rms 0e+00: 8.397e-02 (k 14, axis 0.92); 1e+04: 8.138e-02 (k 14, axis 0.82); 2e+04: 1.273e-01 (k 14, axis 0.93); 3e+04: 1.016e-01 (k 14, axis 0.84); 4e+04: 1.775e-01 (k 14, axis 0.92); 5e+04: 1.678e-01 (k 14, axis 0.93); 6e+04: 2.137e-01 (k 14, axis 0.89); 7e+04: 2.538e-01 (k 14, axis 0.94); 8e+04: 2.543e-01 (k 14, axis 0.88); 9e+04: 3.697e-01 (k 14, axis 0.94); 1e+05: 2.887e-01 (k 14, axis 0.90)
    first 5e3 s: 8.397e-02 -> 9.977e-02; e-fold first 5e3 s 3.020e+04, 0-3e4 s 7.522e+04, 5e4-1e5 s 7.634e+04, whole 1e5 s 6.581e+04 s; last peak k 14, last axis share 0.90
- Reading: at 0.74 R_ph,IC the band grows from 0.074 to 0.479 km/s in 1e5 s: e-fold 5.2e4 s (0-3e4 s), 5.8e4 s (last 5e4 s), 6.4e4 s (whole); M35C itself had 5.1e4 s (2.6e5-3.5e5 s). At 0.72: 7.1e4 s, at 0.76: 6.6e4 s (M35C 5.7e4, 5.0e4). Peak k 13-14 at all three radii (M35C 14), axis share 0.85-0.94. The harness reproduces M35C's growth and band within the scatter. The band rms is not monotonic (the 1e4 s samples jump by +/- 20-30 %), so the e-folds have about that scatter; an arm within 25 % of R0 is not distinguishable.
- The first 5e3 s: at 0.74 0.074 -> 0.092, 0.76 0.084 -> 0.100, 0.72 unchanged; one dump pair only (the restart's first dump is at 3.552e5 s).

## R1
- 2026-11-25 13:42:26 rc=0 ended: Terminating on time limit wall=2199 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- keys: keys: rad_vet/vet_gd_gl_nmu=0 rad_vet/vet_gd_ls=10 common: time/tlim=4.5e5 output2/dt=2400 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=; 8253 rows, dt median 12.1253, min 0.5578, rows below 11.5: 5
== R1: 42 dumps, t 3.5000e+05 to 4.5000e+05
  r 0.72: band A rms 0e+00: 3.521e-02 (k 14, axis 0.86); 1e+04: 4.383e-02 (k 13, axis 0.89); 2e+04: 4.372e-02 (k 13, axis 0.85); 3e+04: 5.966e-02 (k 13, axis 0.90); 4e+04: 5.512e-02 (k 13, axis 0.85); 5e+04: 8.439e-02 (k 13, axis 0.90); 6e+04: 7.304e-02 (k 13, axis 0.85); 7e+04: 1.118e-01 (k 13, axis 0.90); 8e+04: 9.703e-02 (k 13, axis 0.86); 9e+04: 1.466e-01 (k 13, axis 0.90); 1e+05: 1.425e-01 (k 14, axis 0.88)
    first 5e3 s: 3.521e-02 -> 3.516e-02; e-fold first 5e3 s -3.709e+06, 0-3e4 s 6.872e+04, 5e4-1e5 s 7.606e+04, whole 1e5 s 7.104e+04 s; last peak k 14, last axis share 0.88
  r 0.74: band A rms 0e+00: 7.443e-02 (k 14, axis 0.82); 1e+04: 1.149e-01 (k 14, axis 0.91); 2e+04: 9.938e-02 (k 14, axis 0.86); 3e+04: 1.548e-01 (k 14, axis 0.91); 4e+04: 1.385e-01 (k 14, axis 0.89); 5e+04: 1.829e-01 (k 14, axis 0.88); 6e+04: 2.404e-01 (k 14, axis 0.93); 7e+04: 2.234e-01 (k 14, axis 0.88); 8e+04: 3.399e-01 (k 14, axis 0.92); 9e+04: 2.758e-01 (k 13, axis 0.89); 1e+05: 4.788e-01 (k 14, axis 0.92)
    first 5e3 s: 7.443e-02 -> 9.222e-02; e-fold first 5e3 s 2.430e+04, 0-3e4 s 5.198e+04, 5e4-1e5 s 5.775e+04, whole 1e5 s 6.437e+04 s; last peak k 14, last axis share 0.92
  r 0.76: band A rms 0e+00: 8.397e-02 (k 14, axis 0.92); 1e+04: 8.138e-02 (k 14, axis 0.82); 2e+04: 1.273e-01 (k 14, axis 0.93); 3e+04: 1.016e-01 (k 14, axis 0.84); 4e+04: 1.775e-01 (k 14, axis 0.92); 5e+04: 1.678e-01 (k 14, axis 0.93); 6e+04: 2.137e-01 (k 14, axis 0.89); 7e+04: 2.538e-01 (k 14, axis 0.94); 8e+04: 2.544e-01 (k 14, axis 0.88); 9e+04: 3.697e-01 (k 14, axis 0.94); 1e+05: 2.887e-01 (k 14, axis 0.90)
    first 5e3 s: 8.397e-02 -> 9.977e-02; e-fold first 5e3 s 3.020e+04, 0-3e4 s 7.522e+04, 5e4-1e5 s 7.631e+04, whole 1e5 s 6.580e+04 s; last peak k 14, last axis share 0.90

### R1 checked against R0 (hst and v_r maps; added to the R1 entry)
- The key was in use: R1's log prints "vet_gd direction set: 120 directions ... level-symmetric LQ10" and rebuilds the closure arrays (per-shell halo, ragged band intensity 1.406 GB against 1.477 GB of R0's GL7x18, 126 directions); the restart does not ignore it.
hst rows R0 8253 R1 8253; |R1/R0 - 1| of L_top, KE_int, KEr_int, L_mid, Etot at rows: 
  row 0 t 3.50000e+05: t diff 0.0e+00; L_top 4.74e-02, KE_int 1.67e-10, KEr_int 3.45e-10, L_mid 3.00e-13, Etot 4.46e-10
  row 1 t 3.50001e+05: t diff 0.0e+00; L_top 4.66e-02, KE_int 1.67e-09, KEr_int 1.24e-09, L_mid 5.24e-14, Etot 1.32e-09
  row 10 t 3.50078e+05: t diff -2.8e-09; L_top 5.99e-05, KE_int 3.33e-07, KEr_int 1.24e-07, L_mid 4.43e-11, Etot 7.18e-09
  row 100 t 3.51170e+05: t diff 5.2e-08; L_top 5.99e-06, KE_int 2.26e-06, KEr_int 1.50e-06, L_mid 1.19e-08, Etot 7.35e-09
  row 1000 t 3.62083e+05: t diff 8.5e-07; L_top 1.29e-05, KE_int 6.66e-05, KEr_int 2.41e-05, L_mid 2.71e-06, Etot 6.17e-09
  row 3000 t 3.86334e+05: t diff 2.0e-04; L_top 1.41e-03, KE_int 2.63e-04, KEr_int 1.13e-03, L_mid 3.16e-05, Etot 4.63e-08
  row 8000 t 4.46959e+05: t diff 1.4e-01; L_top 1.27e-02, KE_int 1.40e-01, KEr_int 1.84e-01, L_mid 1.59e-02, Etot 1.83e-06
  row 8252 t 4.50000e+05: t diff 0.0e+00; L_top 6.20e-03, KE_int 1.01e-01, KEr_int 1.16e-01, L_mid 8.45e-03, Etot 9.37e-07
  rms over all rows: L_top 6.94e-03, KE_int 4.38e-02, KEr_int 6.76e-02
map rms of v_r difference / rms of v_r (R0), after the shell mean; R = r/R_ph,IC: 0.74 0.80 0.87 0.90 0.93 0.95 0.98 1.00 1.03
  t 3.5521e+05 / 3.5521e+05: 3.00e-07 (rms 9.87e-02); 9.80e-08 (rms 8.35e-02); 2.69e-05 (rms 8.90e-01); 1.95e-05 (rms 1.21e+01); 5.73e-05 (rms 3.22e+01); 2.44e-04 (rms 3.68e+01); 5.54e-03 (rms 2.49e+01); 7.37e-02 (rms 3.88e+01); 1.88e-01 (rms 3.87e+01)
  t 3.5761e+05 / 3.5761e+05: 3.47e-07 (rms 1.14e-01); 3.54e-06 (rms 7.75e-02); 3.85e-05 (rms 8.31e-01); 1.55e-04 (rms 1.26e+01); 3.53e-04 (rms 3.23e+01); 7.34e-04 (rms 3.74e+01); 1.93e-02 (rms 2.55e+01); 1.03e-01 (rms 3.58e+01); 3.46e-01 (rms 4.56e+01)
  t 3.6241e+05 / 3.6241e+05: 2.20e-06 (rms 1.23e-01); 5.54e-05 (rms 8.73e-02); 5.97e-04 (rms 1.23e+00); 2.10e-03 (rms 1.03e+01); 2.55e-03 (rms 3.18e+01); 3.76e-03 (rms 3.59e+01); 2.83e-02 (rms 2.36e+01); 1.65e-01 (rms 4.92e+01); 5.38e-01 (rms 4.94e+01)
  t 3.6960e+05 / 3.6960e+05: 2.86e-05 (rms 1.06e-01); 5.16e-04 (rms 8.06e-02); 1.34e-03 (rms 9.44e-01); 1.10e-02 (rms 1.07e+01); 1.45e-02 (rms 3.24e+01); 2.24e-02 (rms 3.50e+01); 6.02e-02 (rms 2.36e+01); 1.58e-01 (rms 4.25e+01); 3.92e-01 (rms 3.15e+01)
  t 3.8160e+05 / 3.8160e+05: 3.58e-04 (rms 1.55e-01); 6.28e-03 (rms 8.73e-02); 4.58e-02 (rms 1.19e+00); 1.83e-01 (rms 8.04e+00); 1.37e-01 (rms 2.96e+01); 1.26e-01 (rms 3.76e+01); 2.07e-01 (rms 2.43e+01); 2.40e-01 (rms 4.70e+01); 2.66e-01 (rms 1.51e+01)
  t 4.0560e+05 / 4.0560e+05: 1.06e-02 (rms 1.82e-01); 1.02e-01 (rms 1.14e-01); 6.74e-01 (rms 7.82e-01); 1.16e+00 (rms 1.01e+01); 9.88e-01 (rms 3.05e+01); 9.53e-01 (rms 3.84e+01); 1.13e+00 (rms 2.57e+01); 9.93e-01 (rms 4.90e+01); 1.19e+00 (rms 3.81e+01)
  t 4.2961e+05 / 4.2961e+05: 1.92e-02 (rms 3.57e-01); 7.63e-01 (rms 7.36e-02); 1.10e+00 (rms 8.56e-01); 1.41e+00 (rms 9.18e+00); 1.31e+00 (rms 3.28e+01); 1.27e+00 (rms 3.77e+01); 1.40e+00 (rms 2.64e+01); 1.75e+00 (rms 3.36e+01); 2.23e+00 (rms 2.84e+01)
  t 4.5000e+05 / 4.5000e+05: 2.59e-02 (rms 4.91e-01); 7.49e-01 (rms 8.68e-02); 1.17e+00 (rms 1.00e+00); 1.12e+00 (rms 1.56e+01); 1.34e+00 (rms 3.15e+01); 1.33e+00 (rms 3.82e+01); 1.72e+00 (rms 2.22e+01); 1.13e+00 (rms 7.34e+01); 1.28e+00 (rms 4.08e+01)
- Reading: the convection zone and photosphere of R1 and R0 separate (chaos) and are uncorrelated by 4.05e5 s: at 0.95 R_ph,IC the map difference over rms v_r goes 2.4e-4 (3.55e5), 2.2e-2 (3.7e5), 0.13 (3.8e5), 0.95 (4.06e5), 1.3 (4.5e5); KE_int differs by 14 % at 4.47e5 s. The ring at 0.74 stays the same: 3e-7, 3e-5, 4e-4, 1.1e-2, 2.6e-2 at the same times, up to 1e5 s later, after the convection zone has decorrelated. At 0.80 it is 0.75 by the end (that shell lies in the transition). So the ring's phase and amplitude at 0.74 do not follow the zone above it: it is not driven from above in this restart, which agrees with an ongoing deep mode.
- L_top differs by 4.7e-2 in the first row (the new ray set needs one closure rebuild) and 6e-5 ten rows later.

## R2
- 2026-11-25 14:16:33 rc=0 ended: Terminating on time limit wall=2238 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- keys: keys: hydro/wb_rmax=6.3630e11 common: time/tlim=4.5e5 output2/dt=2400 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30; 8254 rows, dt median 12.1249, min 0.5578, rows below 11.5: 5
== R2: 42 dumps, t 3.5000e+05 to 4.5000e+05
  r 0.72: band A rms 0e+00: 3.521e-02 (k 14, axis 0.86); 1e+04: 4.392e-02 (k 13, axis 0.89); 2e+04: 4.395e-02 (k 13, axis 0.85); 3e+04: 6.003e-02 (k 13, axis 0.90); 4e+04: 5.566e-02 (k 13, axis 0.85); 5e+04: 8.526e-02 (k 13, axis 0.90); 6e+04: 7.393e-02 (k 13, axis 0.85); 7e+04: 1.133e-01 (k 13, axis 0.90); 8e+04: 9.829e-02 (k 13, axis 0.86); 9e+04: 1.490e-01 (k 13, axis 0.90); 1e+05: 1.441e-01 (k 14, axis 0.88)
    first 5e3 s: 3.521e-02 -> 3.519e-02; e-fold first 5e3 s -9.882e+06, 0-3e4 s 6.760e+04, 5e4-1e5 s 7.553e+04, whole 1e5 s 7.026e+04 s; last peak k 14, last axis share 0.88
  r 0.74: band A rms 0e+00: 7.443e-02 (k 14, axis 0.82); 1e+04: 1.149e-01 (k 14, axis 0.91); 2e+04: 9.930e-02 (k 14, axis 0.86); 3e+04: 1.547e-01 (k 14, axis 0.91); 4e+04: 1.383e-01 (k 14, axis 0.89); 5e+04: 1.830e-01 (k 14, axis 0.88); 6e+04: 2.397e-01 (k 14, axis 0.93); 7e+04: 2.240e-01 (k 14, axis 0.88); 8e+04: 3.390e-01 (k 14, axis 0.92); 9e+04: 2.769e-01 (k 13, axis 0.89); 1e+05: 4.785e-01 (k 14, axis 0.92)
    first 5e3 s: 7.443e-02 -> 9.220e-02; e-fold first 5e3 s 2.433e+04, 0-3e4 s 5.210e+04, 5e4-1e5 s 5.777e+04, whole 1e5 s 6.437e+04 s; last peak k 14, last axis share 0.92
  r 0.76: band A rms 0e+00: 8.397e-02 (k 14, axis 0.92); 1e+04: 8.128e-02 (k 14, axis 0.82); 2e+04: 1.269e-01 (k 14, axis 0.93); 3e+04: 1.012e-01 (k 14, axis 0.84); 4e+04: 1.765e-01 (k 14, axis 0.92); 5e+04: 1.665e-01 (k 14, axis 0.93); 6e+04: 2.125e-01 (k 14, axis 0.89); 7e+04: 2.512e-01 (k 14, axis 0.94); 8e+04: 2.529e-01 (k 14, axis 0.88); 9e+04: 3.655e-01 (k 14, axis 0.94); 1e+05: 2.867e-01 (k 14, axis 0.90)
    first 5e3 s: 8.397e-02 -> 9.973e-02; e-fold first 5e3 s 3.027e+04, 0-3e4 s 7.601e+04, 5e4-1e5 s 7.651e+04, whole 1e5 s 6.619e+04 s; last peak k 14, last axis share 0.90

## R3
- 2026-11-25 14:55:24 rc=0 ended: Terminating on time limit wall=4378 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- keys: keys: time/cfl_number=0.15 common: time/tlim=4.5e5 output2/dt=2400 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30; 16500 rows, dt median 6.0626, min 0.5578, rows below 11.5: 16497
== R3: 42 dumps, t 3.5000e+05 to 4.5000e+05
  r 0.72: band A rms 0e+00: 3.521e-02 (k 14, axis 0.86); 1e+04: 4.383e-02 (k 13, axis 0.89); 2e+04: 4.375e-02 (k 13, axis 0.85); 3e+04: 5.972e-02 (k 13, axis 0.90); 4e+04: 5.521e-02 (k 13, axis 0.85); 5e+04: 8.456e-02 (k 13, axis 0.90); 6e+04: 7.319e-02 (k 13, axis 0.85); 7e+04: 1.121e-01 (k 13, axis 0.90); 8e+04: 9.734e-02 (k 13, axis 0.86); 9e+04: 1.472e-01 (k 13, axis 0.90); 1e+05: 1.431e-01 (k 14, axis 0.88)
    first 5e3 s: 3.521e-02 -> 3.517e-02; e-fold first 5e3 s -4.267e+06, 0-3e4 s 6.856e+04, 5e4-1e5 s 7.578e+04, whole 1e5 s 7.083e+04 s; last peak k 14, last axis share 0.88
  r 0.74: band A rms 0e+00: 7.443e-02 (k 14, axis 0.82); 1e+04: 1.149e-01 (k 14, axis 0.91); 2e+04: 9.949e-02 (k 14, axis 0.86); 3e+04: 1.550e-01 (k 14, axis 0.91); 4e+04: 1.388e-01 (k 14, axis 0.89); 5e+04: 1.834e-01 (k 14, axis 0.88); 6e+04: 2.411e-01 (k 14, axis 0.93); 7e+04: 2.241e-01 (k 14, axis 0.88); 8e+04: 3.412e-01 (k 14, axis 0.92); 9e+04: 2.770e-01 (k 13, axis 0.89); 1e+05: 4.810e-01 (k 14, axis 0.92)
    first 5e3 s: 7.443e-02 -> 9.224e-02; e-fold first 5e3 s 2.427e+04, 0-3e4 s 5.184e+04, 5e4-1e5 s 5.759e+04, whole 1e5 s 6.417e+04 s; last peak k 14, last axis share 0.92
  r 0.76: band A rms 0e+00: 8.397e-02 (k 14, axis 0.92); 1e+04: 8.147e-02 (k 14, axis 0.82); 2e+04: 1.275e-01 (k 14, axis 0.93); 3e+04: 1.018e-01 (k 14, axis 0.84); 4e+04: 1.779e-01 (k 14, axis 0.92); 5e+04: 1.682e-01 (k 14, axis 0.93); 6e+04: 2.145e-01 (k 14, axis 0.89); 7e+04: 2.549e-01 (k 14, axis 0.94); 8e+04: 2.555e-01 (k 14, axis 0.88); 9e+04: 3.716e-01 (k 14, axis 0.94); 1e+05: 2.903e-01 (k 14, axis 0.90)
    first 5e3 s: 8.397e-02 -> 9.983e-02; e-fold first 5e3 s 3.008e+04, 0-3e4 s 7.487e+04, 5e4-1e5 s 7.608e+04, whole 1e5 s 6.558e+04 s; last peak k 14, last axis share 0.90

## R4
- 2026-11-25 14:53:28 rc=0 ended: Terminating on time limit wall=2212 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- keys: keys: rad_vet/clean_stag_upwind_central=false rad_vet/clean_stag_central=false common: time/tlim=4.5e5 output2/dt=2400 output3/dt=1e30 output5/dt=1e30; 8253 rows, dt median 12.1253, min 0.5578, rows below 11.5: 5
== R4: 42 dumps, t 3.5000e+05 to 4.5000e+05
  r 0.72: band A rms 0e+00: 3.521e-02 (k 14, axis 0.86); 1e+04: 4.383e-02 (k 13, axis 0.89); 2e+04: 4.372e-02 (k 13, axis 0.85); 3e+04: 5.966e-02 (k 13, axis 0.90); 4e+04: 5.512e-02 (k 13, axis 0.85); 5e+04: 8.439e-02 (k 13, axis 0.90); 6e+04: 7.303e-02 (k 13, axis 0.85); 7e+04: 1.118e-01 (k 13, axis 0.90); 8e+04: 9.702e-02 (k 13, axis 0.86); 9e+04: 1.466e-01 (k 13, axis 0.90); 1e+05: 1.425e-01 (k 14, axis 0.88)
    first 5e3 s: 3.521e-02 -> 3.516e-02; e-fold first 5e3 s -3.709e+06, 0-3e4 s 6.872e+04, 5e4-1e5 s 7.604e+04, whole 1e5 s 7.104e+04 s; last peak k 14, last axis share 0.88
  r 0.74: band A rms 0e+00: 7.443e-02 (k 14, axis 0.82); 1e+04: 1.149e-01 (k 14, axis 0.91); 2e+04: 9.938e-02 (k 14, axis 0.86); 3e+04: 1.548e-01 (k 14, axis 0.91); 4e+04: 1.385e-01 (k 14, axis 0.89); 5e+04: 1.829e-01 (k 14, axis 0.88); 6e+04: 2.403e-01 (k 14, axis 0.93); 7e+04: 2.234e-01 (k 14, axis 0.88); 8e+04: 3.399e-01 (k 14, axis 0.92); 9e+04: 2.757e-01 (k 13, axis 0.89); 1e+05: 4.787e-01 (k 14, axis 0.92)
    first 5e3 s: 7.443e-02 -> 9.222e-02; e-fold first 5e3 s 2.430e+04, 0-3e4 s 5.198e+04, 5e4-1e5 s 5.777e+04, whole 1e5 s 6.438e+04 s; last peak k 14, last axis share 0.92
  r 0.76: band A rms 0e+00: 8.397e-02 (k 14, axis 0.92); 1e+04: 8.138e-02 (k 14, axis 0.82); 2e+04: 1.273e-01 (k 14, axis 0.93); 3e+04: 1.016e-01 (k 14, axis 0.84); 4e+04: 1.775e-01 (k 14, axis 0.92); 5e+04: 1.678e-01 (k 14, axis 0.93); 6e+04: 2.137e-01 (k 14, axis 0.89); 7e+04: 2.538e-01 (k 14, axis 0.94); 8e+04: 2.543e-01 (k 14, axis 0.88); 9e+04: 3.697e-01 (k 14, axis 0.94); 1e+05: 2.887e-01 (k 14, axis 0.90)
    first 5e3 s: 8.397e-02 -> 9.977e-02; e-fold first 5e3 s 3.020e+04, 0-3e4 s 7.522e+04, 5e4-1e5 s 7.634e+04, whole 1e5 s 6.581e+04 s; last peak k 14, last axis share 0.90

## 27d: do the keys of R3 and R4 act?
- R3 (time/cfl_number=0.15): 16500 hst rows against R0's 8253 for the same 1e5 s, mean dt half of R0's (it acts). R4 (clean_stag_upwind_central=false, clean_stag_central=false): both keys are FALSE by default and true in the stored deck, so the command-line value equals the default and the log's list of non-default keys simply omits them (no line); the run differs from R0 at once (below), so they act.

### R4 against R0 (hst at equal t; v_r maps, rms of the difference over the rms of v_r, after the shell mean)
hst rows R0 8253 R4 8253; |ARM/R0 - 1| at equal t (R0 interpolated linearly in t):
  t 3.50009e+05: L_top 4.63e-09, KE_int 2.20e-13, KEr_int 1.76e-13, L_mid 4.15e-14, Etot 6.66e-16
  t 3.50102e+05: L_top 1.50e-08, KE_int 2.74e-12, KEr_int 4.02e-12, L_mid 2.80e-13, Etot 2.41e-14
  t 3.51000e+05: L_top 1.23e-07, KE_int 2.66e-09, KEr_int 1.01e-08, L_mid 2.55e-09, Etot 4.99e-13
  t 3.59997e+05: L_top 1.89e-05, KE_int 5.30e-08, KEr_int 2.05e-08, L_mid 4.92e-09, Etot 7.96e-12
  t 3.80005e+05: L_top 3.68e-07, KE_int 3.56e-06, KEr_int 7.80e-06, L_mid 4.37e-07, Etot 5.28e-11
  t 4.00000e+05: L_top 3.54e-05, KE_int 1.67e-04, KEr_int 3.73e-04, L_mid 1.12e-05, Etot 2.34e-11
  t 4.19994e+05: L_top 6.02e-05, KE_int 4.99e-04, KEr_int 3.03e-03, L_mid 7.40e-05, Etot 8.99e-09
  t 4.46995e+05: L_top 7.72e-03, KE_int 1.36e-01, KEr_int 1.75e-01, L_mid 1.88e-03, Etot 3.35e-06
  t 4.49990e+05: L_top 3.88e-03, KE_int 1.41e-01, KEr_int 1.66e-01, L_mid 9.11e-04, Etot 3.62e-06
map rms of v_r difference / rms of v_r (R0), after the shell mean; R = r/R_ph,IC: 0.74 0.80 0.87 0.90 0.93 0.95 0.98 1.00 1.03
  t 3.5521e+05 / 3.5521e+05: 2.69e-07 (rms 9.87e-02); 9.73e-08 (rms 8.35e-02); 2.68e-07 (rms 8.90e-01); 3.44e-07 (rms 1.21e+01); 2.21e-07 (rms 3.22e+01); 2.98e-07 (rms 3.68e+01); 1.01e-05 (rms 2.49e+01); 7.46e-05 (rms 3.88e+01); 5.20e-03 (rms 3.87e+01)
  t 3.5761e+05 / 3.5761e+05: 2.45e-07 (rms 1.14e-01); 1.77e-07 (rms 7.75e-02); 6.41e-07 (rms 8.31e-01); 5.01e-07 (rms 1.26e+01); 5.87e-07 (rms 3.23e+01); 8.93e-07 (rms 3.74e+01); 3.79e-05 (rms 2.55e+01); 1.21e-04 (rms 3.58e+01); 4.79e-03 (rms 4.56e+01)
  t 3.6241e+05 / 3.6241e+05: 2.31e-07 (rms 1.23e-01); 2.28e-07 (rms 8.73e-02); 7.46e-07 (rms 1.23e+00); 2.87e-06 (rms 1.03e+01); 3.91e-06 (rms 3.18e+01); 4.25e-06 (rms 3.59e+01); 2.54e-05 (rms 2.36e+01); 3.01e-04 (rms 4.92e+01); 3.80e-03 (rms 4.94e+01)
  t 3.6960e+05 / 3.6960e+05: 2.80e-07 (rms 1.06e-01); 5.59e-07 (rms 8.06e-02); 2.33e-06 (rms 9.44e-01); 2.41e-05 (rms 1.07e+01); 2.33e-05 (rms 3.24e+01); 2.68e-05 (rms 3.50e+01); 8.70e-05 (rms 2.36e+01); 5.13e-04 (rms 4.25e+01); 1.04e-02 (rms 3.15e+01)
  t 3.8160e+05 / 3.8160e+05: 7.99e-07 (rms 1.55e-01); 7.91e-06 (rms 8.73e-02); 5.47e-05 (rms 1.19e+00); 1.95e-04 (rms 8.04e+00); 1.56e-04 (rms 2.96e+01); 1.34e-04 (rms 3.76e+01); 4.43e-04 (rms 2.43e+01); 3.91e-04 (rms 4.70e+01); 1.93e-03 (rms 1.51e+01)
  t 4.0560e+05 / 4.0560e+05: 1.89e-04 (rms 1.82e-01); 1.63e-03 (rms 1.14e-01); 1.93e-02 (rms 7.82e-01); 5.76e-02 (rms 1.01e+01); 3.09e-02 (rms 3.05e+01); 1.58e-02 (rms 3.84e+01); 5.59e-02 (rms 2.57e+01); 3.90e-02 (rms 4.90e+01); 8.31e-02 (rms 3.81e+01)
  t 4.2961e+05 / 4.2961e+05: 7.09e-04 (rms 3.57e-01); 2.22e-02 (rms 7.36e-02); 8.14e-02 (rms 8.56e-01); 4.77e-01 (rms 9.18e+00); 3.81e-01 (rms 3.28e+01); 3.81e-01 (rms 3.77e+01); 5.10e-01 (rms 2.64e+01); 5.35e-01 (rms 3.36e+01); 6.02e-01 (rms 2.84e+01)
  t 4.5000e+05 / 4.5000e+05: 6.01e-03 (rms 4.91e-01); 1.52e-01 (rms 8.68e-02); 6.52e-01 (rms 1.00e+00); 1.01e+00 (rms 1.56e+01); 1.00e+00 (rms 3.15e+01); 1.09e+00 (rms 3.82e+01); 1.13e+00 (rms 2.22e+01); 9.84e-01 (rms 7.34e+01); 1.21e+00 (rms 4.08e+01)

### R3 against R0 (hst at equal t; v_r maps, rms of the difference over the rms of v_r, after the shell mean)
hst rows R0 8253 R3 16500; |ARM/R0 - 1| at equal t (R0 interpolated linearly in t):
  t 3.50009e+05: L_top 0.00e+00, KE_int 0.00e+00, KEr_int 0.00e+00, L_mid 0.00e+00, Etot 0.00e+00
  t 3.50100e+05: L_top 1.47e-05, KE_int 6.61e-06, KEr_int 1.05e-05, L_mid 2.77e-07, Etot 9.27e-11
  t 3.51003e+05: L_top 3.07e-05, KE_int 3.75e-05, KEr_int 6.53e-05, L_mid 1.03e-06, Etot 6.57e-10
  t 3.60001e+05: L_top 4.91e-05, KE_int 3.64e-05, KEr_int 5.08e-05, L_mid 2.73e-05, Etot 1.83e-09
  t 3.80002e+05: L_top 1.92e-04, KE_int 1.80e-04, KEr_int 1.29e-04, L_mid 2.31e-05, Etot 2.93e-10
  t 4.00003e+05: L_top 1.48e-03, KE_int 1.92e-03, KEr_int 1.07e-02, L_mid 8.89e-05, Etot 5.90e-08
  t 4.19997e+05: L_top 2.70e-03, KE_int 1.64e-02, KEr_int 2.05e-02, L_mid 4.87e-03, Etot 1.17e-07
  t 4.46998e+05: L_top 2.14e-03, KE_int 1.00e-01, KEr_int 1.55e-01, L_mid 1.71e-02, Etot 7.08e-07
  t 4.49987e+05: L_top 6.97e-04, KE_int 1.08e-01, KEr_int 1.50e-01, L_mid 3.95e-03, Etot 2.95e-07
map rms of v_r difference / rms of v_r (R0), after the shell mean; R = r/R_ph,IC: 0.74 0.80 0.87 0.90 0.93 0.95 0.98 1.00 1.03
  t 3.5521e+05 / 3.5520e+05: 6.64e-04 (rms 9.87e-02); 1.53e-03 (rms 8.35e-02); 1.61e-03 (rms 8.90e-01); 1.72e-03 (rms 1.21e+01); 1.41e-03 (rms 3.22e+01); 1.25e-03 (rms 3.68e+01); 2.37e-03 (rms 2.49e+01); 7.66e-03 (rms 3.88e+01); 3.20e-02 (rms 3.87e+01)
  t 3.5761e+05 / 3.5761e+05: 6.44e-04 (rms 1.14e-01); 1.09e-03 (rms 7.75e-02); 2.04e-03 (rms 8.31e-01); 1.55e-03 (rms 1.26e+01); 1.40e-03 (rms 3.23e+01); 1.34e-03 (rms 3.74e+01); 3.83e-03 (rms 2.55e+01); 1.01e-02 (rms 3.58e+01); 5.03e-02 (rms 4.56e+01)
  t 3.6241e+05 / 3.6240e+05: 1.56e-03 (rms 1.23e-01); 4.16e-03 (rms 8.73e-02); 2.99e-03 (rms 1.23e+00); 5.29e-03 (rms 1.03e+01); 4.02e-03 (rms 3.18e+01); 3.67e-03 (rms 3.59e+01); 9.00e-03 (rms 2.36e+01); 2.93e-02 (rms 4.92e+01); 1.20e-01 (rms 4.94e+01)
  t 3.6960e+05 / 3.6960e+05: 1.46e-03 (rms 1.06e-01); 1.14e-03 (rms 8.06e-02); 1.78e-03 (rms 9.44e-01); 6.06e-03 (rms 1.07e+01); 3.80e-03 (rms 3.24e+01); 4.52e-03 (rms 3.50e+01); 6.53e-03 (rms 2.36e+01); 1.34e-02 (rms 4.25e+01); 8.93e-02 (rms 3.15e+01)
  t 3.8160e+05 / 3.8160e+05: 1.70e-03 (rms 1.55e-01); 2.31e-03 (rms 8.73e-02); 3.90e-03 (rms 1.19e+00); 1.94e-02 (rms 8.04e+00); 2.24e-02 (rms 2.96e+01); 1.97e-02 (rms 3.76e+01); 2.85e-02 (rms 2.43e+01); 4.48e-02 (rms 4.70e+01); 5.64e-02 (rms 1.51e+01)
  t 4.0560e+05 / 4.0560e+05: 3.97e-03 (rms 1.82e-01); 2.83e-02 (rms 1.14e-01); 1.97e-01 (rms 7.82e-01); 5.83e-01 (rms 1.01e+01); 3.77e-01 (rms 3.05e+01); 3.55e-01 (rms 3.84e+01); 5.45e-01 (rms 2.57e+01); 4.23e-01 (rms 4.90e+01); 5.72e-01 (rms 3.81e+01)
  t 4.2961e+05 / 4.2961e+05: 1.21e-02 (rms 3.57e-01); 3.04e-01 (rms 7.36e-02); 6.01e-01 (rms 8.56e-01); 1.33e+00 (rms 9.18e+00); 1.13e+00 (rms 3.28e+01); 1.18e+00 (rms 3.77e+01); 1.28e+00 (rms 2.64e+01); 1.28e+00 (rms 3.36e+01); 1.46e+00 (rms 2.84e+01)
  t 4.5000e+05 / 4.5000e+05: 1.86e-02 (rms 4.91e-01); 5.44e-01 (rms 8.68e-02); 9.53e-01 (rms 1.00e+00); 1.17e+00 (rms 1.56e+01); 1.35e+00 (rms 3.15e+01); 1.32e+00 (rms 3.82e+01); 1.24e+00 (rms 2.22e+01); 1.11e+00 (rms 7.34e+01); 1.20e+00 (rms 4.08e+01)

- Reading: R4 starts at 1e-7 (hst) and 5e-3 at 1.03 R_ph,IC in the first dump (the advective face flux matters where the light streams, at the top), grows to order one in the convection zone by 4.3e5 to 4.5e5 s (KE_int 14 %); the ring at 0.74 R_ph,IC differs by 6e-3 at the end. R3 differs at once at 1e-3 everywhere (the time step), the convection zone decorrelates (order one from 4.06e5 s), KE_int 11 %; the ring at 0.74 differs by 1.9 % at the end. Both keys act, and as in R1 the ring at 0.74 stays the same while the zone above decorrelates.

## R5
- 2026-11-25 15:32:15 rc=0 ended: Terminating on time limit wall=2196 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- keys: keys: hydro/fofc=false hydro/fofc_retest=false common: time/tlim=4.5e5 output2/dt=2400 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30; 8253 rows, dt median 12.1253, min 0.5578, rows below 11.5: 5
== R5: 42 dumps, t 3.5000e+05 to 4.5000e+05
  r 0.72: band A rms 0e+00: 3.521e-02 (k 14, axis 0.86); 1e+04: 4.383e-02 (k 13, axis 0.89); 2e+04: 4.372e-02 (k 13, axis 0.85); 3e+04: 5.966e-02 (k 13, axis 0.90); 4e+04: 5.512e-02 (k 13, axis 0.85); 5e+04: 8.439e-02 (k 13, axis 0.90); 6e+04: 7.303e-02 (k 13, axis 0.85); 7e+04: 1.118e-01 (k 13, axis 0.90); 8e+04: 9.702e-02 (k 13, axis 0.86); 9e+04: 1.466e-01 (k 13, axis 0.90); 1e+05: 1.425e-01 (k 14, axis 0.88)
    first 5e3 s: 3.521e-02 -> 3.516e-02; e-fold first 5e3 s -3.709e+06, 0-3e4 s 6.872e+04, 5e4-1e5 s 7.604e+04, whole 1e5 s 7.104e+04 s; last peak k 14, last axis share 0.88
  r 0.74: band A rms 0e+00: 7.443e-02 (k 14, axis 0.82); 1e+04: 1.149e-01 (k 14, axis 0.91); 2e+04: 9.938e-02 (k 14, axis 0.86); 3e+04: 1.548e-01 (k 14, axis 0.91); 4e+04: 1.385e-01 (k 14, axis 0.89); 5e+04: 1.829e-01 (k 14, axis 0.88); 6e+04: 2.403e-01 (k 14, axis 0.93); 7e+04: 2.234e-01 (k 14, axis 0.88); 8e+04: 3.399e-01 (k 14, axis 0.92); 9e+04: 2.757e-01 (k 13, axis 0.89); 1e+05: 4.787e-01 (k 14, axis 0.92)
    first 5e3 s: 7.443e-02 -> 9.222e-02; e-fold first 5e3 s 2.430e+04, 0-3e4 s 5.198e+04, 5e4-1e5 s 5.776e+04, whole 1e5 s 6.438e+04 s; last peak k 14, last axis share 0.92
  r 0.76: band A rms 0e+00: 8.397e-02 (k 14, axis 0.92); 1e+04: 8.138e-02 (k 14, axis 0.82); 2e+04: 1.273e-01 (k 14, axis 0.93); 3e+04: 1.016e-01 (k 14, axis 0.84); 4e+04: 1.775e-01 (k 14, axis 0.92); 5e+04: 1.678e-01 (k 14, axis 0.93); 6e+04: 2.137e-01 (k 14, axis 0.89); 7e+04: 2.538e-01 (k 14, axis 0.94); 8e+04: 2.543e-01 (k 14, axis 0.88); 9e+04: 3.697e-01 (k 14, axis 0.94); 1e+05: 2.887e-01 (k 14, axis 0.90)
    first 5e3 s: 8.397e-02 -> 9.977e-02; e-fold first 5e3 s 3.020e+04, 0-3e4 s 7.522e+04, 5e4-1e5 s 7.634e+04, whole 1e5 s 6.581e+04 s; last peak k 14, last axis share 0.90

## R5 (hydro/fofc=false, fofc_retest=false) and the summary
- R5 was accepted at restart (no FATAL), so R5b (donor cell) was not run. The log shows fofc_retest 0 and "retested 0" each cycle: FOFC is off. rc=0, wall 2196 s, CW / FATAL / DT_RETRY 0.
- It acts: against R0 (hst at equal t; v_r maps as for R3 / R4):
hst rows R0 8253 R5 8253; |ARM/R0 - 1| at equal t (R0 interpolated linearly in t):
  t 3.50009e+05: L_top 8.01e-12, KE_int 4.22e-15, KEr_int 1.17e-13, L_mid 1.55e-13, Etot 1.33e-15
  t 3.50102e+05: L_top 2.65e-10, KE_int 1.64e-14, KEr_int 1.53e-12, L_mid 2.50e-12, Etot 4.44e-16
  t 3.51000e+05: L_top 3.94e-08, KE_int 1.29e-11, KEr_int 6.31e-12, L_mid 2.83e-12, Etot 2.74e-14
  t 3.59997e+05: L_top 5.47e-06, KE_int 4.35e-08, KEr_int 4.24e-08, L_mid 9.26e-09, Etot 7.35e-13
  t 3.80005e+05: L_top 4.74e-07, KE_int 5.22e-07, KEr_int 5.57e-08, L_mid 7.29e-08, Etot 5.16e-12
  t 4.00000e+05: L_top 3.01e-06, KE_int 1.56e-05, KEr_int 3.64e-05, L_mid 8.02e-07, Etot 6.64e-10
  t 4.19994e+05: L_top 1.58e-04, KE_int 2.70e-04, KEr_int 8.72e-04, L_mid 3.64e-05, Etot 2.76e-09
  t 4.46995e+05: L_top 1.95e-04, KE_int 2.17e-03, KEr_int 1.39e-02, L_mid 2.94e-03, Etot 1.12e-07
  t 4.49990e+05: L_top 3.52e-03, KE_int 1.67e-02, KEr_int 3.94e-02, L_mid 3.09e-03, Etot 3.49e-07
map rms of v_r difference / rms of v_r (R0), after the shell mean; R = r/R_ph,IC: 0.74 0.80 0.87 0.90 0.93 0.95 0.98 1.00 1.03
  t 3.5521e+05 / 3.5521e+05: 2.70e-07 (rms 9.87e-02); 9.59e-08 (rms 8.35e-02); 5.50e-08 (rms 8.90e-01); 1.39e-07 (rms 1.21e+01); 6.66e-08 (rms 3.22e+01); 1.81e-07 (rms 3.68e+01); 2.04e-06 (rms 2.49e+01); 6.98e-05 (rms 3.88e+01); 4.66e-02 (rms 3.87e+01)
  t 3.5761e+05 / 3.5761e+05: 2.40e-07 (rms 1.14e-01); 1.43e-07 (rms 7.75e-02); 3.62e-07 (rms 8.31e-01); 2.65e-07 (rms 1.26e+01); 2.45e-07 (rms 3.23e+01); 3.89e-07 (rms 3.74e+01); 3.76e-05 (rms 2.55e+01); 1.17e-04 (rms 3.58e+01); 6.99e-02 (rms 4.56e+01)
  t 3.6241e+05 / 3.6241e+05: 2.24e-07 (rms 1.23e-01); 2.06e-07 (rms 8.73e-02); 3.74e-07 (rms 1.23e+00); 1.05e-06 (rms 1.03e+01); 1.09e-06 (rms 3.18e+01); 1.20e-06 (rms 3.59e+01); 1.86e-05 (rms 2.36e+01); 1.30e-02 (rms 4.92e+01); 7.83e-02 (rms 4.94e+01)
  t 3.6960e+05 / 3.6960e+05: 2.78e-07 (rms 1.06e-01); 4.40e-07 (rms 8.06e-02); 9.01e-07 (rms 9.44e-01); 8.66e-06 (rms 1.07e+01); 7.11e-06 (rms 3.24e+01); 1.69e-05 (rms 3.50e+01); 7.36e-05 (rms 2.36e+01); 1.32e-03 (rms 4.25e+01); 1.13e-01 (rms 3.15e+01)
  t 3.8160e+05 / 3.8160e+05: 3.64e-07 (rms 1.55e-01); 5.17e-06 (rms 8.73e-02); 3.59e-05 (rms 1.19e+00); 2.19e-04 (rms 8.04e+00); 1.42e-04 (rms 2.96e+01); 1.23e-04 (rms 3.76e+01); 3.44e-04 (rms 2.43e+01); 3.26e-04 (rms 4.70e+01); 2.07e-02 (rms 1.51e+01)
  t 4.0560e+05 / 4.0560e+05: 1.20e-05 (rms 1.82e-01); 1.81e-04 (rms 1.14e-01); 1.74e-03 (rms 7.82e-01); 5.52e-03 (rms 1.01e+01); 3.63e-03 (rms 3.05e+01); 3.01e-03 (rms 3.84e+01); 6.51e-03 (rms 2.57e+01); 8.14e-03 (rms 4.90e+01); 1.51e-01 (rms 3.81e+01)
  t 4.2961e+05 / 4.2961e+05: 1.50e-04 (rms 3.57e-01); 3.44e-03 (rms 7.36e-02); 1.29e-02 (rms 8.56e-01); 8.20e-02 (rms 9.18e+00); 5.92e-02 (rms 3.28e+01); 6.93e-02 (rms 3.77e+01); 1.18e-01 (rms 2.64e+01); 2.20e-01 (rms 3.36e+01); 2.61e-01 (rms 2.84e+01)
  t 4.5000e+05 / 4.5000e+05: 1.47e-03 (rms 4.91e-01); 4.03e-02 (rms 8.68e-02); 2.23e-01 (rms 1.00e+00); 5.54e-01 (rms 1.56e+01); 6.50e-01 (rms 3.15e+01); 6.26e-01 (rms 3.82e+01); 9.02e-01 (rms 2.22e+01); 8.29e-01 (rms 7.34e+01); 5.34e-01 (rms 4.08e+01)
- Band readings (a0dc's definitions): equal to R0 to the 4th digit at all three radii (e-fold at 0.74 R_ph,IC 5.198e4 / 5.776e4 / 6.438e4 s; rms at 1e5 s as R0; peak k 14, axis 0.92); full table below.
== R5: 42 dumps, t 3.5000e+05 to 4.5000e+05
  r 0.72: band A rms 0e+00: 3.521e-02 (k 14, axis 0.86); 1e+04: 4.383e-02 (k 13, axis 0.89); 2e+04: 4.372e-02 (k 13, axis 0.85); 3e+04: 5.966e-02 (k 13, axis 0.90); 4e+04: 5.512e-02 (k 13, axis 0.85); 5e+04: 8.439e-02 (k 13, axis 0.90); 6e+04: 7.303e-02 (k 13, axis 0.85); 7e+04: 1.118e-01 (k 13, axis 0.90); 8e+04: 9.702e-02 (k 13, axis 0.86); 9e+04: 1.466e-01 (k 13, axis 0.90); 1e+05: 1.425e-01 (k 14, axis 0.88)
    first 5e3 s: 3.521e-02 -> 3.516e-02; e-fold first 5e3 s -3.709e+06, 0-3e4 s 6.872e+04, 5e4-1e5 s 7.604e+04, whole 1e5 s 7.104e+04 s; last peak k 14, last axis share 0.88
  r 0.74: band A rms 0e+00: 7.443e-02 (k 14, axis 0.82); 1e+04: 1.149e-01 (k 14, axis 0.91); 2e+04: 9.938e-02 (k 14, axis 0.86); 3e+04: 1.548e-01 (k 14, axis 0.91); 4e+04: 1.385e-01 (k 14, axis 0.89); 5e+04: 1.829e-01 (k 14, axis 0.88); 6e+04: 2.403e-01 (k 14, axis 0.93); 7e+04: 2.234e-01 (k 14, axis 0.88); 8e+04: 3.399e-01 (k 14, axis 0.92); 9e+04: 2.757e-01 (k 13, axis 0.89); 1e+05: 4.787e-01 (k 14, axis 0.92)
    first 5e3 s: 7.443e-02 -> 9.222e-02; e-fold first 5e3 s 2.430e+04, 0-3e4 s 5.198e+04, 5e4-1e5 s 5.776e+04, whole 1e5 s 6.438e+04 s; last peak k 14, last axis share 0.92
  r 0.76: band A rms 0e+00: 8.397e-02 (k 14, axis 0.92); 1e+04: 8.138e-02 (k 14, axis 0.82); 2e+04: 1.273e-01 (k 14, axis 0.93); 3e+04: 1.016e-01 (k 14, axis 0.84); 4e+04: 1.775e-01 (k 14, axis 0.92); 5e+04: 1.678e-01 (k 14, axis 0.93); 6e+04: 2.137e-01 (k 14, axis 0.89); 7e+04: 2.538e-01 (k 14, axis 0.94); 8e+04: 2.543e-01 (k 14, axis 0.88); 9e+04: 3.697e-01 (k 14, axis 0.94); 1e+05: 2.887e-01 (k 14, axis 0.90)
    first 5e3 s: 8.397e-02 -> 9.977e-02; e-fold first 5e3 s 3.020e+04, 0-3e4 s 7.522e+04, 5e4-1e5 s 7.634e+04, whole 1e5 s 6.581e+04 s; last peak k 14, last axis share 0.90

### Summary of 26u2
| arm | change | e-fold at 0.74 (first 3e4 / last 5e4 / whole) | band A rms at +1e5 s, 0.74 | peak k | v_r map difference at the end: 0.74 / 0.95 R_ph,IC |
|---|---|---|---|---|---|
| R0 | none | 5.20e4 / 5.78e4 / 6.44e4 s | 0.4787 km/s | 14 | - |
| R1 | LQ10, 120 rays | 5.20e4 / 5.78e4 / 6.44e4 | 0.4788 | 14 | 2.6e-2 / 1.33 |
| R2 | wb off | 5.21e4 / 5.78e4 / 6.44e4 | 0.4785 | 14 | (band 1 %) |
| R3 | cfl 0.15 | in the band table, within 0.5 % (orion) | | 14 | 1.9e-2 / 1.32 |
| R4 | central fluxes off | within 0.5 % | | 14 | 6.0e-3 / 1.09 |
| R5 | FOFC off | 5.20e4 / 5.78e4 / 6.44e4 | as R0 | 14 | 1.5e-3 / 0.63 |
- None of the five changes the ring by more than the 25 % (e-fold) or 20 % (peak k) of the fixed reading; every key acts (R1, R3, R4, R5 shown; R2's band moves by 1 %), and in each arm the convection zone decorrelates from R0 while the ring at 0.74 R_ph,IC stays within 0.2-2.6 % of R0. So none of these terms makes the ring, and it is not driven from the zone above.

## R6
- 2026-11-25 16:15:05 rc=0 ended: Terminating on time limit wall=2234 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- keys: keys: hydro/rsolver=hlle common: time/tlim=4.5e5 output2/dt=2400 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30; 8253 rows, dt median 12.1262, min 0.5578, rows below 11.5: 5
== R6: 42 dumps, t 3.5000e+05 to 4.5000e+05
  r 0.72: band A rms 0e+00: 3.521e-02 (k 14, axis 0.86); 1e+04: 5.147e-03 (k 12, axis 0.83); 2e+04: 1.572e-03 (k 12, axis 0.71); 3e+04: 7.690e-04 (k 6, axis 0.55); 4e+04: 7.043e-04 (k 6, axis 0.61); 5e+04: 4.852e-04 (k 6, axis 0.43); 6e+04: 4.577e-04 (k 6, axis 0.58); 7e+04: 4.450e-04 (k 7, axis 0.52); 8e+04: 3.315e-04 (k 6, axis 0.53); 9e+04: 3.773e-04 (k 6, axis 0.60); 1e+05: 2.824e-04 (k 7, axis 0.54)
    first 5e3 s: 3.521e-02 -> 1.103e-02; e-fold first 5e3 s -4.487e+03, 0-3e4 s -8.074e+03, 5e4-1e5 s -1.055e+05, whole 1e5 s -3.279e+04 s; last peak k 7, last axis share 0.54
  r 0.74: band A rms 0e+00: 7.443e-02 (k 14, axis 0.82); 1e+04: 9.507e-03 (k 14, axis 0.87); 2e+04: 2.196e-03 (k 12, axis 0.76); 3e+04: 9.933e-04 (k 6, axis 0.68); 4e+04: 7.752e-04 (k 6, axis 0.58); 5e+04: 4.914e-04 (k 6, axis 0.57); 6e+04: 4.817e-04 (k 6, axis 0.54); 7e+04: 3.283e-04 (k 6, axis 0.52); 8e+04: 3.635e-04 (k 6, axis 0.57); 9e+04: 2.911e-04 (k 8, axis 0.52); 1e+05: 2.688e-04 (k 6, axis 0.59)
    first 5e3 s: 7.443e-02 -> 2.332e-02; e-fold first 5e3 s -4.487e+03, 0-3e4 s -6.875e+03, 5e4-1e5 s -6.752e+04, whole 1e5 s -2.560e+04 s; last peak k 6, last axis share 0.59
  r 0.76: band A rms 0e+00: 8.397e-02 (k 14, axis 0.92); 1e+04: 7.477e-03 (k 14, axis 0.82); 2e+04: 1.998e-03 (k 14, axis 0.88); 3e+04: 8.506e-04 (k 12, axis 0.76); 4e+04: 4.576e-04 (k 6, axis 0.71); 5e+04: 3.050e-04 (k 6, axis 0.60); 6e+04: 2.144e-04 (k 6, axis 0.52); 7e+04: 1.773e-04 (k 6, axis 0.52); 8e+04: 1.533e-04 (k 6, axis 0.47); 9e+04: 1.170e-04 (k 6, axis 0.51); 1e+05: 1.117e-04 (k 6, axis 0.44)
    first 5e3 s: 8.397e-02 -> 2.202e-02; e-fold first 5e3 s -3.891e+03, 0-3e4 s -6.674e+03, 5e4-1e5 s -5.020e+04, whole 1e5 s -2.048e+04 s; last peak k 6, last axis share 0.44
