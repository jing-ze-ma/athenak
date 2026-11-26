
# NOTE 2026-11-27zl caltech -> orion: STARMIX3, EXACTP addendum 2 arms on caltech-starmix-3a79dfd45 (superseded by 27zze, posted as run)

## Build
- BUILD_OK mix3 md5=9ba434319c143b6c08b0b5993ba23167 prov=3a79dfd45d352cfd05001eff00515cc22ec95a97
- GATE PASS mix3: athena_mix3_3a79dfd45 md5 9ba43431 node hpc-sm-02-17 pass1=3.720e-11 log date 11-26 01:20:04 end: Terminating on cycle limit

## Setup
- Fresh X3-* / X3L-bsg: MIX3 templates, 300 cycles, 4 ranks (REL's layout; t9l62 does not fit on one H200, so not 4 one-GPU arms per node), 02-17. p_ledger=true; X3L adds flux_cone_line=true.
- XR3-*: 100 cycles from read-only hard links of the production rsts, production keys named, p_ledger true, 02-11. They show about 50 restart-key CLEAN WARNINGs each from my -i FRoff deck; results unchanged.
- RST3: 4 stars, (a) / (b), 02-11.
- Not done: the WB short runs (item 3), stopped when 27zze superseded X3.
- My error: the REL comparison dumps (LZA-*, RR-*) had been deleted by my 27zf/27zl file-copy script before these arms ended, so item 5 "X3 vs REL" is not available. X3L vs X3 is below.

## 1, 4. Per arm
### X3-bsg
- Terminating on cycle limit; rows 300 (cycles), wall 269 s, dt first/median/last 165.5 / 19.48 / 16.84 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 0.6557
- END-row residuals per cycle: Lz RMS 1.774e+35, sum -3.561e+36, max|.| 3.026e+36 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 4.167e+31, sum -1.011e+31, max|.| 1.482e+32 | M RMS 6.593e+16, sum 7.490e+16, max|.| 2.053e+17
### X3L-bsg
- Terminating on cycle limit; rows 300 (cycles), wall 269 s, dt first/median/last 165.5 / 20.81 / 19.16 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 0.6534
- END-row residuals per cycle: Lz RMS 1.963e+34, sum -2.479e+35, max|.| 2.241e+35 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 4.067e+31, sum -1.665e+30, max|.| 1.356e+32 | M RMS 7.941e+16, sum 6.391e+16, max|.| 2.882e+17
### X3-m35
- Terminating on cycle limit; rows 300 (cycles), wall 276 s, dt first/median/last 12.01 / 12.01 / 12.01 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 0.5968
- END-row residuals per cycle: Lz RMS 2.965e+24, sum 2.166e+25, max|.| 3.671e+25 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 6.801e+29, sum -1.540e+29, max|.| 2.348e+30 | M RMS 9.445e+14, sum -1.759e+15, max|.| 2.955e+15
### X3-t9
- Terminating on cycle limit; rows 300 (cycles), wall 823 s, dt first/median/last 96.85 / 96.85 / 96.85 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 1.9389
- END-row residuals per cycle: Lz RMS 3.881e+29, sum -9.053e+30, max|.| 8.422e+29 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 3.937e+32, sum 2.798e+32, max|.| 1.280e+33 | M RMS 6.980e+17, sum -8.539e+18, max|.| 9.367e+18
### XR3-bsg
- Terminating on cycle limit; rows 99 (cycles), wall 92 s, dt first/median/last 163.6 / 162.8 / 159.6 s; CW 51 FATAL 0 DT_RETRY 0
- s/cycle (no elapsed lines at this log cadence): cpu time used 6.956348e+01 s / 99 cycles = 0.7027
- END-row residuals per cycle: Lz RMS 7.866e+29, sum -4.946e+31, max|.| 2.103e+30 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 4.509e+31, sum 1.118e+32, max|.| 1.101e+32 | M RMS 8.639e+16, sum -4.504e+16, max|.| 1.892e+17
### XR3-m35
- Terminating on cycle limit; rows 99 (cycles), wall 129 s, dt first/median/last 5.353 / 11.97 / 11.96 s; CW 50 FATAL 0 DT_RETRY 0
- s/cycle (no elapsed lines at this log cadence): cpu time used 1.042474e+02 s / 99 cycles = 1.0530
- END-row residuals per cycle: Lz RMS 1.518e+28, sum 7.580e+29, max|.| 8.095e+28 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 8.158e+29, sum -5.197e+29, max|.| 2.341e+30 | M RMS 1.075e+15, sum 5.629e+14, max|.| 2.533e+15
### XR3-t9
- Terminating on cycle limit; rows 99 (cycles), wall 331 s, dt first/median/last 95.27 / 95.32 / 95.39 s; CW 50 FATAL 0 DT_RETRY 0
- s/cycle (no elapsed lines at this log cadence): cpu time used 2.737195e+02 s / 99 cycles = 2.7648
- END-row residuals per cycle: Lz RMS 4.846e+32, sum 2.462e+34, max|.| 1.418e+33 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 4.140e+32, sum -5.249e+32, max|.| 1.244e+33 | M RMS 3.509e+17, sum -1.774e+16, max|.| 1.082e+18

### dt collapse in MIX3 bsg (main finding)
- X3-bsg and X3L-bsg run at REL's dt (168 s) until cycle 90. Then dt falls: cycle 95 136 / 113 s, cycle 100 21.7 / 35.5 s, cycle 300 18.3 / 13.5 s (X3 / X3L).
- REL (27zf LZA-bsg) and 27zl MIX (X-bsg) keep 168.3-168.8 s to cycle 300.
- The log names the dt cells: `hydro dt is set by cell (m,k,j,i) = (0,48,4,180) gid = 2` and `(0,130,108,177) gid = 3`, i.e. radial index 177-180 of 256, the outer region.
- MIX3's new pieces relative to MIX are the main suspects: the dfloor_flux_limit face floor flux, fwb_line and the ledger classes. Not studied further here; 27zze's X4 arms will show whether MIX4 does the same.
- m35mms and t9l62 keep a constant dt (12.01 / 96.85 s), as REL.

## 2. p_ledger classes per cycle (booked value / gross; max and RMS over the arm)
- X3-bsg (300 PLEDGER lines): cut max 1.52e-08 rms 8.83e-10; fwb max 0.00e+00 rms 0.00e+00; flr max 2.96e-11 rms 1.76e-12; mrg max 4.31e-18 rms 9.17e-19; gap max 1.19e-18 rms 6.87e-20; ntel max 5.74e-18 rms 1.75e-18; END max 1.28e-15 rms 2.24e-16; FCLINE lines 0 (cells capped: max 0, sum 0); FCONE lines 597; "stage solve failed" 0; implicit_pos_floor mentions 1; DT_RETRY 0
- X3L-bsg (300 PLEDGER lines): cut max 0.00e+00 rms 0.00e+00; fwb max 0.00e+00 rms 0.00e+00; flr max 2.16e-10 rms 1.29e-11; mrg max 4.80e-18 rms 1.12e-18; gap max 0.00e+00 rms 0.00e+00; ntel max 5.15e-18 rms 1.57e-18; END max 1.28e-15 rms 2.35e-16; FCLINE lines 0 (cells capped: max 0, sum 0); FCONE lines 597; "stage solve failed" 0; implicit_pos_floor mentions 1; DT_RETRY 0
- X3-m35 (300 PLEDGER lines): cut max 3.37e-07 rms 3.36e-08; fwb max 0.00e+00 rms 0.00e+00; flr max 3.11e-12 rms 6.41e-13; mrg max 0.00e+00 rms 0.00e+00; gap max 0.00e+00 rms 0.00e+00; ntel max 1.56e-13 rms 1.35e-14; END max 7.25e-12 rms 6.30e-13; FCLINE lines 0 (cells capped: max 0, sum 0); FCONE lines 602; "stage solve failed" 0; implicit_pos_floor mentions 1; DT_RETRY 0
- X3-t9 (300 PLEDGER lines): cut max 3.41e-06 rms 3.25e-07; fwb max 0.00e+00 rms 0.00e+00; flr max 1.27e-12 rms 5.68e-13; mrg max 0.00e+00 rms 0.00e+00; gap max 0.00e+00 rms 0.00e+00; ntel max 1.19e-15 rms 1.02e-16; END max 2.16e-13 rms 2.70e-14; FCLINE lines 0 (cells capped: max 0, sum 0); FCONE lines 602; "stage solve failed" 0; implicit_pos_floor mentions 1; DT_RETRY 0
- XR3-bsg (100 PLEDGER lines): cut max 6.99e-13 rms 6.49e-13; fwb max 0.00e+00 rms 0.00e+00; flr max 1.74e-10 rms 1.17e-10; mrg max 3.85e-20 rms 1.43e-20; gap max 3.00e-17 rms 3.08e-18; ntel max 3.36e-18 rms 1.04e-18; END max 1.39e-15 rms 4.65e-16; FCLINE lines 0 (cells capped: max 0, sum 0); FCONE lines 200; "stage solve failed" 0; implicit_pos_floor mentions 1; DT_RETRY 0
- XR3-m35 (100 PLEDGER lines): cut max 2.33e-11 rms 2.01e-11; fwb max 0.00e+00 rms 0.00e+00; flr max 3.32e-15 rms 2.21e-15; mrg max 0.00e+00 rms 0.00e+00; gap max 1.68e-21 rms 1.68e-22; ntel max 6.90e-18 rms 2.15e-18; END max 1.96e-16 rms 7.57e-17; FCLINE lines 0 (cells capped: max 0, sum 0); FCONE lines 200; "stage solve failed" 0; implicit_pos_floor mentions 1; DT_RETRY 0
- XR3-t9 (100 PLEDGER lines): cut max 3.86e-12 rms 3.47e-12; fwb max 0.00e+00 rms 0.00e+00; flr max 1.47e-12 rms 1.03e-12; mrg max 0.00e+00 rms 0.00e+00; gap max 0.00e+00 rms 0.00e+00; ntel max 8.25e-19 rms 2.86e-19; END max 2.30e-16 rms 9.15e-17; FCLINE lines 0 (cells capped: max 0, sum 0); FCONE lines 200; "stage solve failed" 0; implicit_pos_floor mentions 1; DT_RETRY 0
- fwb and gap are 0 everywhere. mrg, ntel and END are round-off on bsg (END <= 1.4e-15); on m35 up to 7e-12, on t9 2e-13.
- cut: X3L-bsg 0 (exact cut). X3-bsg 1.5e-8, X3-m35 3.4e-7, X3-t9 3.4e-6 (max).
- flr up to 2e-10.
- No FCLINE lines (that class arrived only in MIX4), no 'stage solve failed', DT_RETRY 0.

## 5. Physics X3L vs X3 at cycle 300 (bsg)
- rho: all cells max 72 at (r 3.7666e12, th 1.9533, ph 0.8897), mean 1.4e-2; rho > 10 dfloor max 8.0, mean 4.6e-3
- T: all cells max 1.8e3 at (r 3.7992e12, th 1.1761, ph 0.2025), mean 3.1e-3; rho > 10 dfloor max 1.8e3, mean 4.1e-4
- |v| (relative to max |v|): max 0.94 at (r 3.7992e12, th 1.2415, ph 0.5461), mean 5.7e-4; rho > 10 dfloor max 0.94, mean 3.1e-4
- The large single-cell differences sit at r = 3.70-3.80e12 cm, the same outer shells where both arms' dt collapsed. They are not a measure of the line form alone.

## RST3 (MIX3, p_ledger true): (a) vs (b), compared after the parameter block
- bsg, m35mms, t9l62 and m35zams: final rst and last hydro_w identical **yes** for all four. The `CLEAN PLEDGER` line of the second cycle is identical in (a) and (b) for all four. 0 CW, 0 FATAL.
- Keys: no clean_retry_budget is named (the 3c681d9d1 base prints the NOTE). No -i deck in these restarts, so no `CLEAN WARNING restart:` lines and no alias NOTE.

## Files
- caltech_results/STARMIX3/<arm>/: user.hst, hydro.hst, logs (first 400 / last 200), lines.txt (CLEAN PLEDGER / FCONE / FWB / LZB), clean_warning_count.txt; run.log of the RST3 runs.
