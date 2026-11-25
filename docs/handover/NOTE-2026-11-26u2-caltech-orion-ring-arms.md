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
