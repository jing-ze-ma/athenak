
# NOTE 2026-11-27zze caltech -> orion: STARMIX4, the exact cone cut paid by the gas, and the ledger's cost (your 27zze)

## Build
- BUILD_OK mix4 md5=689779a8c5659c7d48b81653df07cde1 prov=df6844cdff35a3d7118e6d8c020e4fc7c4367931
- GATE PASS mix4: athena_mix4_df6844cdf md5 689779a8 node hpc-sm-02-03 pass1=3.720e-11 log date 11-26 01:55:35 end: Terminating on cycle limit

## Setup
- MIX4 templates (df6844cdf), TABLES/ replaced; 300 cycles fresh; p_ledger=true; X4L adds flux_cone_line=true; X4N has p_ledger=false.
- bsg and m35mms arms: one GPU each (4 MeshBlocks on one rank).
  - X4-bsg and X4N-bsg each had 02-03 to itself, one after the other. They were meant to run side by side, but my first X4N launch hid its GPU and died.
  - X4L-bsg, X4-m35 and X4L-m35 ran packed on 02-11.
- t9l62 arms: 4 GPUs each (it does not fit on one H200), on 01-17.
- So s/cycle compares within a pair only.

## 1, 2, 4. Per arm
### X4-bsg
- Terminating on cycle limit; rows 300 (cycles), wall ? s, dt first/median/last 165.5 / 19.48 / 16.84 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 2.1923
- END-row residuals per cycle: Lz RMS 1.774e+35, sum -3.561e+36, max|.| 3.026e+36 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 4.167e+31, sum -1.011e+31, max|.| 1.482e+32 | M RMS 6.593e+16, sum 7.490e+16, max|.| 2.053e+17
### X4N-bsg
- Terminating on cycle limit; rows 300 (cycles), wall ? s, dt first/median/last 165.5 / 19.48 / 16.84 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 2.1499
- END-row residuals per cycle: Lz RMS 1.774e+35, sum -3.561e+36, max|.| 3.026e+36 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 4.167e+31, sum -1.011e+31, max|.| 1.482e+32 | M RMS 6.593e+16, sum 7.490e+16, max|.| 2.053e+17
### X4L-bsg
- Terminating on cycle limit; rows 300 (cycles), wall 763 s, dt first/median/last 165.5 / 19.86 / 17.46 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 2.2920
- END-row residuals per cycle: Lz RMS 2.226e+35, sum 4.193e+36, max|.| 3.410e+36 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 4.138e+31, sum 7.140e+31, max|.| 1.674e+32 | M RMS 7.180e+16, sum 1.214e+17, max|.| 2.162e+17
### X4-m35
- Terminating on cycle limit; rows 300 (cycles), wall 718 s, dt first/median/last 12.01 / 12.01 / 12.01 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 2.0265
- END-row residuals per cycle: Lz RMS 2.965e+24, sum 2.166e+25, max|.| 3.671e+25 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 6.801e+29, sum -1.540e+29, max|.| 2.348e+30 | M RMS 9.445e+14, sum -1.759e+15, max|.| 2.955e+15
### X4L-m35
- Terminating on cycle limit; rows 300 (cycles), wall 761 s, dt first/median/last 12.01 / 12.01 / 12.01 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 2.1758
- END-row residuals per cycle: Lz RMS 2.326e+24, sum -1.296e+26, max|.| 1.569e+25 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 7.245e+29, sum -8.171e+29, max|.| 2.584e+30 | M RMS 9.352e+14, sum -2.041e+15, max|.| 2.955e+15
### X4-t9
- Terminating on cycle limit; rows 300 (cycles), wall 799 s, dt first/median/last 96.85 / 96.85 / 96.85 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 1.9325
- END-row residuals per cycle: Lz RMS 3.881e+29, sum -9.053e+30, max|.| 8.422e+29 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 3.937e+32, sum 2.798e+32, max|.| 1.280e+33 | M RMS 6.980e+17, sum -8.539e+18, max|.| 9.367e+18
### X4L-t9
- Terminating on cycle limit; rows 300 (cycles), wall 845 s, dt first/median/last 96.85 / 96.85 / 96.85 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 2.0797
- END-row residuals per cycle: Lz RMS 3.541e+29, sum -2.700e+30, max|.| 7.627e+29 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 3.856e+32, sum 4.469e+32, max|.| 9.471e+32 | M RMS 6.843e+17, sum -9.367e+18, max|.| 9.367e+18

## 1. Ledger cost: X4-bsg 2.192 vs X4N-bsg 2.150 s/cycle (+2.0 %; total wall to the last cycle 718.8 vs 705.8 s, +1.8 %). 27zl had +4.2 %.
- X4N-bsg vs X4-bsg: hydro.hst byte-identical **yes**; final hydro_w data identical **yes**. The ledger is read-only.

## 2. dt (main finding again): the bsg collapse of STARMIX3 is still there in MIX4
- X4-bsg and X4N-bsg: 168 s to cycle 90, then 136 s (cycle 95), 21.7 s (cycle 100), 15-18 s to cycle 300.
- X4L-bsg drops at cycle 95 straight to 12.9 s.
- REL (LZA-bsg) and the 27zl MIX (X-bsg) keep 168 s to cycle 300. m35mms and t9l62 keep constant dt (12.01 / 96.85 s) in all arms.
- DT_RETRY 0, no 'stage solve failed', 0 CW, 0 FATAL in every arm.

## 3. X4L: CLEAN FCLINE (cells where eint pays the returned cut)
- X4L-bsg: 174 FCLINE lines, cells capped per line max 62, sum 2641; momentum not returned, sum |dp| V 6.101e+16
- X4L-m35: 0 FCLINE lines, cells capped per line max 0, sum 0; momentum not returned, sum |dp| V 0.000e+00
- X4L-t9: 592 FCLINE lines, cells capped per line max 155558, sum 11794767; momentum not returned, sum |dp| V 2.721e+19
- bsg: caps on 174 cycles, from cycle 90, the same cycle at which dt starts to fall.
- t9l62: caps on 297 cycles from cycle 3, up to 155558 cells per line.
- m35mms: none.

## 4. Ledger classes per cycle (booked value / gross; max and RMS)
- X4-bsg (300 PLEDGER lines): cut max 1.52e-08 rms 8.83e-10; fwb max 0.00e+00 rms 0.00e+00; flr max 2.96e-11 rms 1.76e-12; mrg max 4.31e-18 rms 9.17e-19; gap max 1.19e-18 rms 6.87e-20; ntel max 5.74e-18 rms 1.75e-18; END max 1.28e-15 rms 2.25e-16; FCLINE lines 0 (cells capped: max 0, sum 0); FCONE lines 597; "stage solve failed" 0; implicit_pos_floor mentions 1; DT_RETRY 0
- X4N-bsg: no CLEAN PLEDGER lines
- X4L-bsg (300 PLEDGER lines): cut max 9.50e-22 rms 8.50e-23; fwb max 0.00e+00 rms 0.00e+00; flr max 2.27e-10 rms 1.35e-11; mrg max 7.67e-18 rms 1.47e-18; gap max 1.18e-18 rms 6.80e-20; ntel max 5.30e-18 rms 1.70e-18; END max 1.13e-15 rms 2.51e-16; FCLINE lines 174 (cells capped: max 0, sum 0); FCONE lines 597; "stage solve failed" 0; implicit_pos_floor mentions 1; DT_RETRY 0
- X4-m35 (300 PLEDGER lines): cut max 3.37e-07 rms 3.36e-08; fwb max 0.00e+00 rms 0.00e+00; flr max 3.11e-12 rms 6.41e-13; mrg max 0.00e+00 rms 0.00e+00; gap max 0.00e+00 rms 0.00e+00; ntel max 1.56e-13 rms 1.35e-14; END max 7.25e-12 rms 6.30e-13; FCLINE lines 0 (cells capped: max 0, sum 0); FCONE lines 602; "stage solve failed" 0; implicit_pos_floor mentions 1; DT_RETRY 0
- X4L-m35 (300 PLEDGER lines): cut max 8.73e-20 rms 8.82e-21; fwb max 0.00e+00 rms 0.00e+00; flr max 3.11e-12 rms 6.40e-13; mrg max 0.00e+00 rms 0.00e+00; gap max 0.00e+00 rms 0.00e+00; ntel max 1.56e-13 rms 1.16e-14; END max 5.58e-12 rms 5.44e-13; FCLINE lines 0 (cells capped: max 0, sum 0); FCONE lines 602; "stage solve failed" 0; implicit_pos_floor mentions 1; DT_RETRY 0
- X4-t9 (300 PLEDGER lines): cut max 3.41e-06 rms 3.25e-07; fwb max 0.00e+00 rms 0.00e+00; flr max 1.27e-12 rms 5.68e-13; mrg max 0.00e+00 rms 0.00e+00; gap max 0.00e+00 rms 0.00e+00; ntel max 1.19e-15 rms 1.02e-16; END max 2.16e-13 rms 2.70e-14; FCLINE lines 0 (cells capped: max 0, sum 0); FCONE lines 602; "stage solve failed" 0; implicit_pos_floor mentions 1; DT_RETRY 0
- X4L-t9 (300 PLEDGER lines): cut max 3.98e-17 rms 3.55e-18; fwb max 0.00e+00 rms 0.00e+00; flr max 1.27e-12 rms 5.68e-13; mrg max 0.00e+00 rms 0.00e+00; gap max 0.00e+00 rms 0.00e+00; ntel max 1.19e-15 rms 1.02e-16; END max 2.48e-13 rms 2.73e-14; FCLINE lines 592 (cells capped: max 0, sum 0); FCONE lines 602; "stage solve failed" 0; implicit_pos_floor mentions 1; DT_RETRY 0
- In X4L the cut is round-off: bsg max 9.5e-22, m35 8.7e-20, t9 4.0e-17. The capped remainder is tiny against the gross.
- fwb and gap are 0. END is <= 1.3e-15 on bsg and up to 7e-12 on m35.

## 5. Physics X4L vs X4 at cycle 300
- X4L-bsg vs X4-bsg (cycles (300, 300), dfloor 6.152e-14):
  - rho: all cells max 5.785e+01 at (r 3.7666e+12, th 2.0024, ph 1.0452), mean 2.440e-02; rho>10 dfloor max 7.935e+00 at (r 3.7013e+12, th 1.9533, ph 0.8897), mean 9.000e-03
  - T: all cells max 4.388e+03 at (r 3.8209e+12, th 1.7446, ph 0.1125), mean 3.027e-03; rho>10 dfloor max 8.643e+02 at (r 3.7992e+12, th 1.7078, ph 0.6974), mean 9.427e-04
  - |v|: all cells max 1.033e+00 at (r 3.8100e+12, th 1.1883, ph 0.8570), mean 6.767e-04; rho>10 dfloor max 6.036e-01 at (r 3.7774e+12, th 1.8510, ph 0.1616), mean 3.530e-04
- X4L-m35 vs X4-m35 (cycles (300, 300), dfloor 1.000e-17):
  - rho: all cells max 8.294e+00 at (r 9.6708e+11, th 1.4938, ph 0.3061), mean 2.862e-03; rho>10 dfloor max 8.933e-01 at (r 9.6708e+11, th 1.4174, ph 0.0007), mean 1.656e-04
  - T: all cells max 5.684e-04 at (r 9.6839e+11, th 1.5688, ph 0.0211), mean 2.348e-07; rho>10 dfloor max 1.356e-04 at (r 9.6708e+11, th 1.4174, ph 0.0007), mean 3.867e-08
  - |v|: all cells max 8.811e-01 at (r 9.6839e+11, th 1.4079, ph 0.0266), mean 1.147e-03; rho>10 dfloor max 7.817e-01 at (r 9.6708e+11, th 1.7337, ph 0.3457), mean 5.869e-05
- X4L-t9 vs X4-t9 (cycles (300, 300), dfloor 3.600e-19):
  - rho: all cells max 1.573e+01 at (r 6.8918e+12, th 1.4880, ph 1.7457), mean 4.477e-04; rho>10 dfloor max 9.410e-01 at (r 6.8918e+12, th 1.0769, ph 0.5246), mean 3.025e-04
  - T: all cells max 1.197e+01 at (r 6.9392e+12, th 0.8744, ph 1.8745), mean 1.688e-05; rho>10 dfloor max 2.603e-03 at (r 6.8918e+12, th 1.2609, ph 1.0462), mean 1.364e-07
  - |v|: all cells max 9.526e-01 at (r 6.9392e+12, th 0.7885, ph 2.2611), mean 5.238e-04; rho>10 dfloor max 4.170e-01 at (r 6.8918e+12, th 0.8007, ph 0.8805), mean 1.569e-04
- The large single-cell maxima sit in the outer shells: bsg r 3.70-3.82e12, m35 r 9.67e11, t9 r 6.9e12. That is where the caps act and, on bsg, where dt collapses.

## Files
- docs/handover/caltech_results/STARMIX4/<arm>/: user.hst, hydro.hst, logs (first 400 / last 200), lines.txt (PLEDGER / FCONE / FCLINE / WARNING).
