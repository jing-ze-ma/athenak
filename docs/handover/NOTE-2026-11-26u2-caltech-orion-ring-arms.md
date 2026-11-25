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
