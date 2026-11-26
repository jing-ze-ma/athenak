
# NOTE 2026-11-27zl caltech -> orion: STARMIX, the merged fix binary against the release on three stars (your 27zl)

## Build
- BUILD_OK mix md5=bbc5dcbd21ed225881e3dbea1a9df39f prov=bd1503befeead2aa536b3574961eb0836d2392f9
- GATE PASS mix: athena_mix_bd1503bef md5 bbc5dcbd node hpc-sm-02-09 pass1=3.720e-11 log date 11-25 22:59:02 end: Terminating on cycle limit
- worktree wt_orion1121/mix at bd1503bef; CUDA 12.9 + g++13, star_vet.

## Setup
- MIX fresh arms (X-*, MIXNL-bsg): templates of the MIX branch with TABLES/ replaced; hydro/p_ledger=true named (false in MIXNL).
- REL fresh arms (R-*): the 27zf LZA-* arms. Same REL binary and aa2646b92 template; the two LZ keys are named at their REL default (false), which changes nothing.
- Restart pairs RR (REL) / XR (MIX, p_ledger true): 100 cycles from READ-ONLY hard links of bsg 00009, m35mms 00011, t9l62 00012.
  - Both arms name: the production run's command-line keys, time/tlim=1e9 (the rsts are at tlim), rad_vet/clean_verbose=1, and -i a 2-line deck with job/final_restart=false.
  - **The XR arms also name rad_vet/clean_retry_budget=true.** All production rsts store 0, which MIX FATALs on (my 27zr). true is MIX's only behaviour, so it moves no result.
  - The XR arms print about 50 `CLEAN WARNING restart: <key> = ... from the restart file (the deck does not name it)` lines each (MIX's restart warning, 27zj). That is because of the small -i deck; results unchanged. The RR arms print none (REL has no such warning).
- s/cycle notes: the fresh arms use the elapsed lines at the template cadence (ndiag 50); the restart arms use cpu time / cycles. XR-bsg overlapped nothing; RR-bsg ran its first ~8 minutes next to non-timing restart-key runs on 02-09.

## Per-arm readings (items 1-5)
```
### X-bsg
- Terminating on cycle limit; rows 300 (cycles), wall 284 s, dt first/median/last 165.5 / 168.6 / 168.8 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 0.6922
- END-row residuals per cycle: Lz RMS 9.840e+29, sum -8.588e+31, max|.| 3.906e+30 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 4.482e+31, sum 4.628e+31, max|.| 1.219e+32 | M RMS 7.124e+16, sum 8.679e+16, max|.| 1.892e+17
- CLEAN LZB: 300 lines, ledger first 1.3667e+30 last -1.3887e+34 (growth -1.3889e+34); sum of Lz_hflx*dt over the arm 4.2268e+36
- CLEAN FCONE_LZ: 20 lines, max |X| n/a; first `CLEAN FCONE_LZ cycle 0 t2s 1 momentum removed (signed, V/c^2) x1 5.129269e+30 x2 -2.396935e+14 x3 0.000000e+00`; last `CLEAN FCONE_LZ cycle 10 t2s 2 momentum removed (signed, V/c^2) x1 6.703830e+22 x2 -1.545125e+14 x3 0.000000e+00`
- CLEAN LZSPIN: 3 lines, max |X| n/a; first `CLEAN LZSPIN cycle 1 t2s 3: L_z of the cut and capped momentum X 2.107034e+30, the gas's as dOmega 2.973114e-26 (I 7.086962e+55)`; last `CLEAN LZSPIN cycle 2 t2s 3: L_z of the cut and capped momentum X 1.420739e+30, the gas's as dOmega 2.004722e-26 (I 7.086962e+55)`
- p_ledger columns in user.hst: ['Pr_tot', 'Pr_hflx', 'Pr_ntel', 'Pr_src', 'Pr_rad', 'Pr_flr', 'Pr_mrg', 'Pr_gap'] (per cycle in the committed user.hst); means: Pr_tot 4.3149e+35, Pr_hflx -1.7854e+34, Pr_ntel 1.5725e+17, Pr_src -1.7908e+34, Pr_rad 6.1851e+31, Pr_flr 1.6125e+26, Pr_mrg 7.6240e+16, Pr_gap 0.0000e+00
- RMS of Pr_mrg 1.7630e+17, Pr_rad 6.2410e+31
- CLEAN PLEDGER END lines: 300; last `CLEAN PLEDGER cycle 299 dt 1.688312e+02 hflx -1.789440e+34 ntel 1.748184e+18 src -1.796201e+34 rad 5.095578e+31 flr 2.604726e+26 mrg 0.000000e+00 gap 0.000000e+00 P 3.9640519344829167e+35 END -3.735466e+20 gross 1.726527e+36`
### MIXNL-bsg
- Terminating on cycle limit; rows 300 (cycles), wall 277 s, dt first/median/last 165.5 / 168.6 / 168.8 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 0.6644
- END-row residuals per cycle: Lz RMS 9.840e+29, sum -8.588e+31, max|.| 3.906e+30 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 4.482e+31, sum 4.628e+31, max|.| 1.219e+32 | M RMS 7.124e+16, sum 8.679e+16, max|.| 1.892e+17
- CLEAN LZB: 300 lines, ledger first 1.3667e+30 last -1.3887e+34 (growth -1.3889e+34); sum of Lz_hflx*dt over the arm 4.2268e+36
- CLEAN FCONE_LZ: 20 lines, max |X| n/a; first `CLEAN FCONE_LZ cycle 0 t2s 1 momentum removed (signed, V/c^2) x1 5.129269e+30 x2 -2.396935e+14 x3 0.000000e+00`; last `CLEAN FCONE_LZ cycle 10 t2s 2 momentum removed (signed, V/c^2) x1 6.703830e+22 x2 -1.545125e+14 x3 0.000000e+00`
- CLEAN LZSPIN: 3 lines, max |X| n/a; first `CLEAN LZSPIN cycle 1 t2s 3: L_z of the cut and capped momentum X 2.107034e+30, the gas's as dOmega 2.973114e-26 (I 7.086962e+55)`; last `CLEAN LZSPIN cycle 2 t2s 3: L_z of the cut and capped momentum X 1.420739e+30, the gas's as dOmega 2.004722e-26 (I 7.086962e+55)`
### X-m35
- Terminating on cycle limit; rows 300 (cycles), wall 267 s, dt first/median/last 12.01 / 12.01 / 12.01 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 0.5992
- END-row residuals per cycle: Lz RMS 9.832e+23, sum 9.195e+24, max|.| 9.572e+24 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 7.372e+29, sum 7.763e+29, max|.| 2.107e+30 | M RMS 1.047e+15, sum -2.744e+15, max|.| 3.448e+15
- CLEAN LZB: 300 lines, ledger first 6.9817e+32 last -2.6308e+34 (growth -2.7006e+34); sum of Lz_hflx*dt over the arm 1.5465e+24
- CLEAN FCONE_LZ: none
- CLEAN LZSPIN: none
- p_ledger columns in user.hst: ['Pr_tot', 'Pr_hflx', 'Pr_ntel', 'Pr_src', 'Pr_rad', 'Pr_flr', 'Pr_mrg', 'Pr_gap'] (per cycle in the committed user.hst); means: Pr_tot -4.0888e+31, Pr_hflx -2.1500e+33, Pr_ntel -1.2760e+15, Pr_src -2.1500e+33, Pr_rad -1.9035e+28, Pr_flr 4.3523e+20, Pr_mrg 0.0000e+00, Pr_gap 0.0000e+00
- RMS of Pr_mrg 0.0000e+00, Pr_rad 1.9393e+28
- CLEAN PLEDGER END lines: 300; last `CLEAN PLEDGER cycle 299 dt 1.200762e+01 hflx -2.157169e+33 ntel 1.920317e+17 src -2.157157e+33 rad -1.277075e+28 flr 4.520572e+20 mrg 0.000000e+00 gap 0.000000e+00 P -6.1861181676119516e+31 END 1.661406e+17 gross 6.271005e+31`
### X-t9
- Terminating on cycle limit; rows 300 (cycles), wall 817 s, dt first/median/last 96.85 / 96.85 / 96.85 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 2.0291
- END-row residuals per cycle: Lz RMS 3.128e+27, sum -5.922e+28, max|.| 3.604e+28 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 3.862e+32, sum 2.519e+32, max|.| 1.067e+33 | M RMS 6.813e+17, sum -8.287e+18, max|.| 9.367e+18
- CLEAN LZB: 300 lines, ledger first 7.8092e+30 last 3.4998e+38 (growth 3.4998e+38); sum of Lz_hflx*dt over the arm 2.5599e+25
- CLEAN FCONE_LZ: 602 lines, max |X| n/a; first `CLEAN FCONE_LZ cycle 0 t2s 0 momentum removed (signed, V/c^2) x1 3.631451e+29 x2 9.085045e+17 x3 0.000000e+00`; last `CLEAN FCONE_LZ cycle 299 t2s 3 momentum removed (signed, V/c^2) x1 5.468444e+28 x2 4.905105e+22 x3 0.000000e+00`
- CLEAN LZSPIN: 597 lines, max |X| n/a; first `CLEAN LZSPIN cycle 1 t2s 3: L_z of the cut and capped momentum X 1.782568e+30, the gas's as dOmega 2.610936e-27 (I 6.827314e+56)`; last `CLEAN LZSPIN cycle 299 t2s 3: L_z of the cut and capped momentum X -7.615526e+32, the gas's as dOmega -1.115512e-24 (I 6.826935e+56)`
- p_ledger columns in user.hst: ['Pr_tot', 'Pr_hflx', 'Pr_ntel', 'Pr_src', 'Pr_rad', 'Pr_flr', 'Pr_mrg', 'Pr_gap'] (per cycle in the committed user.hst); means: Pr_tot -3.9180e+35, Pr_hflx -2.5201e+35, Pr_ntel -3.2398e+17, Pr_src -2.5128e+35, Pr_rad -7.3329e+32, Pr_flr 4.1175e+23, Pr_mrg 0.0000e+00, Pr_gap 0.0000e+00
- RMS of Pr_mrg 0.0000e+00, Pr_rad 7.3465e+32
- CLEAN PLEDGER END lines: 300; last `CLEAN PLEDGER cycle 299 dt 9.684942e+01 hflx -2.528718e+35 ntel 0.000000e+00 src -2.521268e+35 rad -7.147626e+32 flr 4.395273e+23 mrg 0.000000e+00 gap 0.000000e+00 P -2.7407469332180895e+35 END 1.180592e+21 gross 2.897544e+35`
### RR-bsg
- Terminating on cycle limit; rows 99 (cycles), wall 89 s, dt first/median/last 163.6 / 162.8 / 159.6 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle (no elapsed lines at this log cadence): cpu time used 6.580058e+01 s / 99 cycles = 0.6647
- END-row residuals per cycle: Lz RMS 7.603e+29, sum -4.767e+31, max|.| 2.087e+30 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 4.124e+31, sum -5.309e+31, max|.| 1.075e+32 | M RMS 7.865e+16, sum 1.441e+17, max|.| 1.937e+17
- CLEAN LZB: 100 lines, ledger first -5.6851e+42 last -5.5826e+42 (growth 1.0248e+41); sum of Lz_hflx*dt over the arm -2.8105e+37
- CLEAN FCONE_LZ: none
- CLEAN LZSPIN: none
### XR-bsg
- Terminating on cycle limit; rows 99 (cycles), wall 95 s, dt first/median/last 163.6 / 162.8 / 159.6 s; CW 51 FATAL 0 DT_RETRY 0
- s/cycle (no elapsed lines at this log cadence): cpu time used 7.218900e+01 s / 99 cycles = 0.7292
- END-row residuals per cycle: Lz RMS 8.415e+29, sum -4.505e+31, max|.| 2.501e+30 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 4.288e+31, sum -4.914e+31, max|.| 1.031e+32 | M RMS 6.962e+16, sum 7.206e+16, max|.| 2.477e+17
- CLEAN LZB: 100 lines, ledger first -5.6851e+42 last -5.5807e+42 (growth 1.0433e+41); sum of Lz_hflx*dt over the arm -5.4104e+36
- CLEAN FCONE_LZ: 200 lines, max |X| n/a; first `CLEAN FCONE_LZ cycle 29657 t2s 2 momentum removed (signed, V/c^2) x1 3.650496e+22 x2 -8.221955e+20 x3 0.000000e+00`; last `CLEAN FCONE_LZ cycle 29756 t2s 3 momentum removed (signed, V/c^2) x1 6.215786e+25 x2 -2.425332e+23 x3 0.000000e+00`
- CLEAN LZSPIN: 200 lines, max |X| n/a; first `CLEAN LZSPIN cycle 29657 t2s 2: L_z of the cut and capped momentum X -1.625910e+32, the gas's as dOmega -2.295185e-24 (I 7.084004e+55)`; last `CLEAN LZSPIN cycle 29756 t2s 3: L_z of the cut and capped momentum X 1.343363e+36, the gas's as dOmega 1.896363e-20 (I 7.083890e+55)`
- p_ledger columns in user.hst: ['Pr_tot', 'Pr_hflx', 'Pr_ntel', 'Pr_src', 'Pr_rad', 'Pr_flr', 'Pr_mrg', 'Pr_gap'] (per cycle in the committed user.hst); means: Pr_tot -1.8727e+34, Pr_hflx -1.7923e+34, Pr_ntel 1.2653e+17, Pr_src -1.7966e+34, Pr_rad 4.2492e+31, Pr_flr 2.6996e+26, Pr_mrg 8.5693e+15, Pr_gap 0.0000e+00
- RMS of Pr_mrg 1.4668e+16, Pr_rad 4.2492e+31
- CLEAN PLEDGER END lines: 100; last `CLEAN PLEDGER cycle 29756 dt 1.595842e+02 hflx -1.792319e+34 ntel 0.000000e+00 src -1.796656e+34 rad 4.269058e+31 flr 2.692198e+26 mrg 0.000000e+00 gap 0.000000e+00 P -2.6891584114901729e+34 END 3.804641e+19 gross 1.098050e+36`
### RR-m35
- Terminating on cycle limit; rows 99 (cycles), wall ? s, dt first/median/last 5.353 / 11.97 / 11.96 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle (no elapsed lines at this log cadence): cpu time used 1.019723e+02 s / 99 cycles = 1.0300
- END-row residuals per cycle: Lz RMS 1.905e+27, sum 1.918e+28, max|.| 6.442e+27 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 8.519e+29, sum -4.239e+28, max|.| 2.219e+30 | M RMS 8.877e+14, sum -9.148e+14, max|.| 2.393e+15
- CLEAN LZB: 100 lines, ledger first -1.4022e+41 last -1.4129e+41 (growth -1.0670e+39); sum of Lz_hflx*dt over the arm -8.0680e+25
- CLEAN FCONE_LZ: none
- CLEAN LZSPIN: none
### XR-m35
- Terminating on cycle limit; rows 99 (cycles), wall 130 s, dt first/median/last 5.353 / 11.97 / 11.96 s; CW 50 FATAL 0 DT_RETRY 0
- s/cycle (no elapsed lines at this log cadence): cpu time used 1.034062e+02 s / 99 cycles = 1.0445
- END-row residuals per cycle: Lz RMS 1.617e+27, sum -9.872e+27, max|.| 5.706e+27 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 6.343e+29, sum -1.938e+29, max|.| 1.890e+30 | M RMS 8.309e+14, sum -8.444e+14, max|.| 2.252e+15
- CLEAN LZB: 100 lines, ledger first -1.4022e+41 last -1.4079e+41 (growth -5.7110e+38); sum of Lz_hflx*dt over the arm 1.1843e+25
- CLEAN FCONE_LZ: 200 lines, max |X| n/a; first `CLEAN FCONE_LZ cycle 73280 t2s 2 momentum removed (signed, V/c^2) x1 8.503267e+22 x2 -4.532910e+20 x3 0.000000e+00`; last `CLEAN FCONE_LZ cycle 73379 t2s 3 momentum removed (signed, V/c^2) x1 5.869126e+25 x2 8.750848e+22 x3 0.000000e+00`
- CLEAN LZSPIN: 200 lines, max |X| n/a; first `CLEAN LZSPIN cycle 73280 t2s 2: L_z of the cut and capped momentum X -1.137424e+34, the gas's as dOmega -6.180356e-20 (I 1.840386e+53)`; last `CLEAN LZSPIN cycle 73379 t2s 3: L_z of the cut and capped momentum X -1.216570e+36, the gas's as dOmega -6.610403e-18 (I 1.840386e+53)`
- p_ledger columns in user.hst: ['Pr_tot', 'Pr_hflx', 'Pr_ntel', 'Pr_src', 'Pr_rad', 'Pr_flr', 'Pr_mrg', 'Pr_gap'] (per cycle in the committed user.hst); means: Pr_tot -1.1861e+31, Pr_hflx -2.1578e+33, Pr_ntel -2.1535e+15, Pr_src -2.1571e+33, Pr_rad -3.2808e+29, Pr_flr 3.2218e+20, Pr_mrg 0.0000e+00, Pr_gap 0.0000e+00
- RMS of Pr_mrg 0.0000e+00, Pr_rad 3.2824e+29
- CLEAN PLEDGER END lines: 100; last `CLEAN PLEDGER cycle 73379 dt 1.196494e+01 hflx -2.157805e+33 ntel 0.000000e+00 src -2.157134e+33 rad -3.468324e+29 flr 3.194264e+20 mrg 0.000000e+00 gap 0.000000e+00 P 1.9460222872429911e+32 END -3.473401e+17 gross 5.599111e+34`
### RR-t9
- Terminating on cycle limit; rows 99 (cycles), wall 327 s, dt first/median/last 95.27 / 95.32 / 95.39 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle (no elapsed lines at this log cadence): cpu time used 2.711122e+02 s / 99 cycles = 2.7385
- END-row residuals per cycle: Lz RMS 6.291e+32, sum 3.489e+34, max|.| 1.485e+33 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 4.252e+32, sum 4.385e+32, max|.| 9.263e+32 | M RMS 4.122e+17, sum -1.073e+17, max|.| 1.125e+18
- CLEAN LZB: 100 lines, ledger first -1.5307e+45 last -1.5524e+45 (growth -2.1672e+43); sum of Lz_hflx*dt over the arm -2.9860e+45
- CLEAN FCONE_LZ: none
- CLEAN LZSPIN: none
### XR-t9
- Terminating on cycle limit; rows 99 (cycles), wall 342 s, dt first/median/last 95.27 / 95.32 / 95.39 s; CW 50 FATAL 0 DT_RETRY 0
- s/cycle (no elapsed lines at this log cadence): cpu time used 2.755277e+02 s / 99 cycles = 2.7831
- END-row residuals per cycle: Lz RMS 6.440e+32, sum 3.680e+34, max|.| 1.490e+33 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 3.975e+32, sum 7.759e+31, max|.| 9.688e+32 | M RMS 3.963e+17, sum -3.070e+17, max|.| 9.860e+17
- CLEAN LZB: 100 lines, ledger first -1.5308e+45 last -1.5577e+45 (growth -2.6962e+43); sum of Lz_hflx*dt over the arm -2.9860e+45
- CLEAN FCONE_LZ: 200 lines, max |X| n/a; first `CLEAN FCONE_LZ cycle 29228 t2s 2 momentum removed (signed, V/c^2) x1 9.603482e+26 x2 -4.166644e+26 x3 0.000000e+00`; last `CLEAN FCONE_LZ cycle 29327 t2s 3 momentum removed (signed, V/c^2) x1 1.070579e+28 x2 -3.635515e+26 x3 0.000000e+00`
- CLEAN LZSPIN: 200 lines, max |X| n/a; first `CLEAN LZSPIN cycle 29228 t2s 2: L_z of the cut and capped momentum X 2.242720e+38, the gas's as dOmega 2.645629e-19 (I 8.477076e+56)`; last `CLEAN LZSPIN cycle 29327 t2s 3: L_z of the cut and capped momentum X 2.543025e+40, the gas's as dOmega 3.003972e-17 (I 8.465542e+56)`
- p_ledger columns in user.hst: ['Pr_tot', 'Pr_hflx', 'Pr_ntel', 'Pr_src', 'Pr_rad', 'Pr_flr', 'Pr_mrg', 'Pr_gap'] (per cycle in the committed user.hst); means: Pr_tot 8.6014e+35, Pr_hflx -2.5309e+35, Pr_ntel -2.2281e+18, Pr_src -2.5244e+35, Pr_rad -6.3460e+32, Pr_flr 6.0138e+25, Pr_mrg 0.0000e+00, Pr_gap 0.0000e+00
- RMS of Pr_mrg 0.0000e+00, Pr_rad 6.3461e+32
- CLEAN PLEDGER END lines: 100; last `CLEAN PLEDGER cycle 29327 dt 9.538710e+01 hflx -2.530857e+35 ntel -2.475369e+19 src -2.524378e+35 rad -6.413587e+32 flr 4.904220e+25 mrg 0.000000e+00 gap 0.000000e+00 P 9.2602131610873333e+35 END 4.270421e+21 gross 6.003334e+37`
```
- The R-* rows are the LZA-* rows in the 27zf NOTE.

## Item 6: physics MIX vs REL at the last cycle
- X-bsg vs LZA-bsg (cycles (300, 300), dfloor 6.152e-14):
  - rho: all cells max 1.287e-02 at (r 3.7883e+12, th 1.2783, ph 0.8693), mean 1.476e-05; rho>10 dfloor max 4.551e-03 at (r 3.7774e+12, th 1.2129, ph 0.8693), mean 1.117e-05
  - T: all cells max 7.475e-05 at (r 3.7448e+12, th 1.2088, ph 0.8529), mean 1.482e-06; rho>10 dfloor max 7.475e-05 at (r 3.7448e+12, th 1.2088, ph 0.8529), mean 1.382e-06
  - |v|: all cells max 4.812e-02 at (r 5.4026e+12, th 1.2415, ph 0.9715), mean 4.572e-03; rho>10 dfloor max 1.813e-03 at (r 3.7774e+12, th 1.2783, ph 0.8733), mean 6.165e-06
- X-m35 vs LZA-m35 (cycles (300, 300), dfloor 1.000e-17):
  - rho: all cells max 1.111e-02 at (r 9.7102e+11, th 1.4024, ph 0.0402), mean 1.709e-05; rho>10 dfloor max 1.292e-03 at (r 9.6839e+11, th 1.7392, ph 0.3088), mean 1.377e-05
  - T: all cells max 2.456e-05 at (r 9.3230e+11, th 1.3969, ph 0.1548), mean 9.879e-07; rho>10 dfloor max 2.456e-05 at (r 9.3230e+11, th 1.3969, ph 0.1548), mean 8.220e-07
  - |v|: all cells max 9.440e-02 at (r 9.7895e+11, th 1.7446, ph 0.1698), mean 1.890e-04; rho>10 dfloor max 7.236e-04 at (r 9.6839e+11, th 1.6792, ph 0.1739), mean 5.423e-06
- X-t9 vs LZA-t9 (cycles (300, 300), dfloor 3.600e-19):
  - rho: all cells max 4.602e-03 at (r 6.9870e+12, th 1.5739, ph 2.1875), mean 1.130e-06; rho>10 dfloor max 2.463e-04 at (r 6.8918e+12, th 1.3591, ph 0.6044), mean 1.440e-06
  - T: all cells max 1.018e-04 at (r 3.3569e+13, th 2.3531, ph 1.5677), mean 8.808e-06; rho>10 dfloor max 3.727e-05 at (r 6.1760e+12, th 2.2181, ph 1.0462), mean 1.027e-06
  - |v|: all cells max 8.365e-02 at (r 5.6130e+13, th 2.3531, ph 0.1687), mean 5.559e-04; rho>10 dfloor max 3.615e-04 at (r 6.8918e+12, th 1.3284, ph 0.5001), mean 2.369e-06
- XR-bsg vs RR-bsg (cycles (29757, 29757), dfloor 6.152e-14):
  - rho: all cells max 2.207e-01 at (r 3.6035e+12, th 1.2947, ph 0.6238), mean 8.788e-06; rho>10 dfloor max 7.803e-02 at (r 3.6361e+12, th 2.0923, ph 0.6197), mean 3.280e-06
  - T: all cells max 6.507e-03 at (r 3.5926e+12, th 1.3438, ph 0.5911), mean 6.174e-06; rho>10 dfloor max 2.395e-03 at (r 3.6470e+12, th 1.5728, ph 0.4929), mean 1.330e-06
  - |v|: all cells max 2.367e-01 at (r 3.6579e+12, th 1.5197, ph 0.5338), mean 1.392e-04; rho>10 dfloor max 5.551e-02 at (r 3.6361e+12, th 2.0923, ph 0.6156), mean 1.336e-06
- XR-m35 vs RR-m35 (cycles (73380, 73380), dfloor 1.000e-17):
  - rho: all cells max 5.442e+00 at (r 9.7763e+11, th 1.4692, ph 0.2379), mean 1.208e-03; rho>10 dfloor max 5.442e+00 at (r 9.7763e+11, th 1.4692, ph 0.2379), mean 1.125e-03
  - T: all cells max 5.779e-01 at (r 9.5017e+11, th 1.6956, ph 0.3402), mean 9.455e-05; rho>10 dfloor max 5.779e-01 at (r 9.5017e+11, th 1.6956, ph 0.3402), mean 6.939e-05
  - |v|: all cells max 3.681e-01 at (r 9.8829e+11, th 1.5701, ph 0.3361), mean 1.192e-03; rho>10 dfloor max 3.567e-01 at (r 9.8695e+11, th 1.5660, ph 0.3347), mean 2.944e-04
- XR-t9 vs RR-t9 (cycles (29328, 29328), dfloor 3.600e-19):
  - rho: all cells max 2.093e+00 at (r 4.5681e+12, th 2.1261, ph 0.1749), mean 7.142e-05; rho>10 dfloor max 2.093e+00 at (r 4.5681e+12, th 2.1261, ph 0.1749), mean 7.142e-05
  - T: all cells max 1.738e-01 at (r 2.1501e+13, th 1.5309, ph 1.0278), mean 4.189e-05; rho>10 dfloor max 1.738e-01 at (r 2.1501e+13, th 1.5309, ph 1.0278), mean 4.189e-05
  - |v|: all cells max 2.446e-01 at (r 4.5369e+12, th 1.8807, ph 1.1137), mean 2.943e-05; rho>10 dfloor max 2.446e-01 at (r 4.5369e+12, th 1.8807, ph 1.1137), mean 2.943e-05

## Byte checks

## MIXNL
- MIXNL-bsg vs X-bsg: hydro.hst byte-identical **yes**; final hydro_w data (after the parameter block) **yes**. p_ledger is read-only.

## RST: (a) 2 cycles vs (b) 1 + rst + 1, MIX with p_ledger true and clean_retry_budget=true named
- RST-MIX-bsg: CW 0 FATAL 0; final rst bsg3d.00010.rst vs bsg3d.00011.rst identical: yes; hydro_w bsg3d.hydro_w.00117.bin vs bsg3d.hydro_w.00117.bin identical: yes
- RST-MIX-m35mms: CW 0 FATAL 0; final rst m35mms.00012.rst vs m35mms.00013.rst identical: yes; hydro_w m35mms.hydro_w.00483.bin vs m35mms.hydro_w.00483.bin identical: yes
- RST-MIX-m35zams: CW 0 FATAL 0; final rst m35zams.00012.rst vs m35zams.00013.rst identical: yes; hydro_w m35zams.hydro_w.00483.bin vs m35zams.hydro_w.00483.bin identical: yes
- RST-MIX-t9l62: CW 0 FATAL 0; final rst t9l62.00013.rst vs t9l62.00014.rst identical: yes; hydro_w t9l62.hydro_w.00252.bin vs t9l62.hydro_w.00252.bin identical: yes
- All 4 stars: final rst and last hydro_w identical. 0 CW (no -i deck in these restarts, so no restart-key warning), 0 FATAL.
- A first pass without clean_retry_budget FATALed at start-up on all 4 stars (27zr). It is kept in lzab_1127/failed_v1/retrybudget.

## Short reading
- Fresh, MIX vs REL at 300 cycles, where rho > 10 dfloor:
  - bsg: rho 4.6e-3, T 7.5e-5, |v| 1.8e-3
  - m35mms: rho 1.3e-3, T 2.5e-5, |v| 7.2e-4
  - t9l62: rho 2.5e-4, T 3.7e-5, |v| 3.6e-4
  - Residual RMS for L_z, E and M are within 10-40 % of REL.
- DFL in MIX: on the bsg restart, the sum of Lz_hflx*dt falls from -2.81e37 (RR) to -5.41e36 (XR), as in 27zf DFL1.
- p_ledger: max |END| / gross is 4.9e-16 (bsg restart), 2.2e-16 (bsg fresh), 2.6e-15 (m35 fresh), 4.1e-15 (t9 fresh), 6.2e-18 (m35 restart), 7.1e-17 (t9 restart), from each arm's last PLEDGER line. RMS of mrg / rad are in the rows above.
- Cost: X-bsg 0.692 vs LZA-bsg 0.658 s/cycle (+5 %); X-t9 2.029 vs 1.980 (+2.5 %); X-m35 0.599 vs 0.569 (+5 %). The restart pairs are not clean timings.

## Files
- docs/handover/caltech_results/STARMIX/<arm>/: user.hst, hydro.hst, logs, lines.txt, clean_warning_count.txt. Readings: docs/handover/caltech_results/STARMIX/readings_rst.txt (RST rows) and the 27zf readings file for the arm rows.
