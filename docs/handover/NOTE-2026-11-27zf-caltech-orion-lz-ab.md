
# NOTE 2026-11-27zf caltech -> orion: L_z-exact keys A/B, dfloor_flux_limit fix, p_ledger, LZD, restart sanity (your 27zf)

## Builds (CUDA 12.9 + g++13, star_vet, HOPPER90, -ffp-contract=off), tiny gate each
- BUILD_OK pled md5=3ae1965592fedd42240a71d726680328 prov=174233149540fe31190d9c79c910a4f88653dce0
- GATE PASS pled: athena_pled_174233149 md5 3ae19655 node hpc-sm-02-11 pass1=3.720e-11 log date 11-25 22:51:00 end: Terminating on cycle limit
- BUILD_OK lzd md5=fe538ddba2e8753b458b328381766caa prov=ca9048dd42fd9b7969e82cc0191412c6f352ee41
- GATE PASS lzd: athena_lzd_ca9048dd4 md5 fe538ddb node hpc-sm-02-11 pass1=3.720e-11 log date 11-25 22:51:01 end: Terminating on cycle limit
- BUILD_OK dfl md5=35fa6d4bfd9310a27042332654654825 prov=2640c8d33eb0d4a7981aa2f38e155a3831820ec8
- GATE PASS dfl: athena_dfl_2640c8d33 md5 35fa6d4b node hpc-sm-02-11 pass1=3.720e-11 log date 11-25 22:50:46 end: Terminating on cycle limit

## Setup
- REL: rel6 01deec542. 4 ranks (production layout) on 02-17 (fresh arms) and 02-09 (restart arms); same node type in each pair.
- Fresh decks: the aa2646b92 templates, TABLES/ replaced by our copies; bsg's profile is the production linhao_profile_full.txt. Restarts: READ-ONLY hard links of bsg 00009, m35mms 00011, t9l62 00012, M35ZAMS 00011.
- Named in every arm: time/nlim, output1/dcycle=1, output2/dcycle=<last cycle>, output3/5/6/7 dt=1e30, and both LZ keys (A false, B true).
  - Fresh arms: output2/dt=1e30 dropped. The release FATALs on a command-line key it never reads, and dt is not read next to dcycle.
- Restarts additionally name:
  - the production run's command-line keys (its chunk EX list);
  - time/tlim=1e9 (every rst is at its run's tlim);
  - rad_vet/clean_verbose=1 (the bsg rst stores 0, and then no CLEAN LZB lines are printed);
  - -i a 2-line deck with job/final_restart=false.
- Dropped per 27zo: LZRA-t9/LZRB-t9 (27zl RR-t9/XR-t9 repeat them). The other arms ran.

## Per-arm readings (items 1-4, 6)
- Per arm: cycles, wall, s/cycle, dt; residual RMS/sum/max per cycle for L_z, E and M; ledger growth from CLEAN LZB; FCONE_LZ / LZSPIN counts and first/last lines.
- In full: docs/handover/caltech_results/LZAB/readings_all_arms.txt.
- s/cycle notes:
  - The median of the last elapsed intervals only exists at the log cadence (ndiag 50 in the templates).
  - The restart arms have cpu time / cycles instead.
  - PLED1's 7.8 s/cycle (wall 802 s vs DFL0 89 s) overlapped two CPU builds on 02-09 and is not a clean timing. The MIX arm with p_ledger on (27zl X-bsg) costs +4 % against MIXNL.
```
### LZA-bsg
- Terminating on cycle limit; rows 300 (cycles), wall 268 s, dt first/median/last 165.5 / 168.6 / 168.8 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 0.6582
- END-row residuals per cycle: Lz RMS 9.743e+29, sum -1.033e+32, max|.| 3.674e+30 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 4.118e+31, sum -2.889e+30, max|.| 1.149e+32 | M RMS 6.982e+16, sum 6.368e+16, max|.| 1.937e+17
- CLEAN LZB: 300 lines, ledger first 4.8032e+30 last -1.1818e+34 (growth -1.1823e+34); sum of Lz_hflx*dt over the arm 1.0216e+30
- CLEAN FCONE_LZ: none
- CLEAN LZSPIN: none
### LZB-bsg
- Terminating on cycle limit; rows 300 (cycles), wall 277 s, dt first/median/last 165.5 / 168.6 / 168.8 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 0.6579
- END-row residuals per cycle: Lz RMS 9.503e+29, sum -8.084e+31, max|.| 3.614e+30 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 4.129e+31, sum -3.766e+31, max|.| 1.114e+32 | M RMS 7.073e+16, sum 4.629e+16, max|.| 1.892e+17
- CLEAN LZB: 300 lines, ledger first 9.4083e+29 last -1.1823e+34 (growth -1.1824e+34); sum of Lz_hflx*dt over the arm -9.1382e+29
- CLEAN FCONE_LZ: 34 lines, max |X| n/a; first `CLEAN FCONE_LZ cycle 0 t2s 1 momentum removed (signed, V/c^2) x1 7.426480e+30 x2 1.231453e+14 x3 0.000000e+00`; last `CLEAN FCONE_LZ cycle 17 t2s 2 momentum removed (signed, V/c^2) x1 2.311362e+23 x2 -2.304503e+14 x3 0.000000e+00`
- CLEAN LZSPIN: 3 lines, max |X| n/a; first `CLEAN LZSPIN cycle 1 t2s 3: L_z of the cut and capped momentum X 2.117066e+30, the gas's as dOmega 2.987269e-26 (I 7.086962e+55)`; last `CLEAN LZSPIN cycle 2 t2s 3: L_z of the cut and capped momentum X 1.424001e+30, the gas's as dOmega 2.009325e-26 (I 7.086962e+55)`
### LZD-bsg
- Terminating on cycle limit; rows 300 (cycles), wall 279 s, dt first/median/last 165.5 / 168.6 / 168.8 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 0.6606
- END-row residuals per cycle: Lz RMS 9.503e+29, sum -8.084e+31, max|.| 3.614e+30 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 4.129e+31, sum -3.766e+31, max|.| 1.114e+32 | M RMS 7.073e+16, sum 4.629e+16, max|.| 1.892e+17
- CLEAN LZB: 300 lines, ledger first 9.4083e+29 last -1.1823e+34 (growth -1.1824e+34); sum of Lz_hflx*dt over the arm -9.1382e+29
- CLEAN FCONE_LZ: 34 lines, max |X| n/a; first `CLEAN FCONE_LZ cycle 0 t2s 1 momentum removed (signed, V/c^2) x1 7.426480e+30 x2 1.231453e+14 x3 0.000000e+00`; last `CLEAN FCONE_LZ cycle 17 t2s 2 momentum removed (signed, V/c^2) x1 2.311362e+23 x2 -2.304503e+14 x3 0.000000e+00`
- CLEAN LZSPIN: 3 lines, max |X| n/a; first `CLEAN LZSPIN cycle 1 t2s 3: L_z of the cut and capped momentum X 2.117066e+30, the gas's as dOmega 2.987269e-26 (I 7.086962e+55)`; last `CLEAN LZSPIN cycle 2 t2s 3: L_z of the cut and capped momentum X 1.424001e+30, the gas's as dOmega 2.009325e-26 (I 7.086962e+55)`
### LZA-m35
- Terminating on cycle limit; rows 300 (cycles), wall 257 s, dt first/median/last 12.01 / 12.01 / 12.01 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 0.5685
- END-row residuals per cycle: Lz RMS 1.345e+24, sum -6.666e+24, max|.| 1.784e+25 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 7.448e+29, sum -1.893e+29, max|.| 2.199e+30 | M RMS 9.527e+14, sum -1.689e+15, max|.| 2.533e+15
- CLEAN LZB: 300 lines, ledger first 6.9797e+32 last -2.6276e+34 (growth -2.6974e+34); sum of Lz_hflx*dt over the arm 1.7798e+24
- CLEAN FCONE_LZ: none
- CLEAN LZSPIN: none
### LZB-m35
- Terminating on cycle limit; rows 300 (cycles), wall 261 s, dt first/median/last 12.01 / 12.01 / 12.01 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 0.5701
- END-row residuals per cycle: Lz RMS 1.345e+24, sum -6.666e+24, max|.| 1.784e+25 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 7.448e+29, sum -1.893e+29, max|.| 2.199e+30 | M RMS 9.527e+14, sum -1.689e+15, max|.| 2.533e+15
- CLEAN LZB: 300 lines, ledger first 6.9797e+32 last -2.6276e+34 (growth -2.6974e+34); sum of Lz_hflx*dt over the arm 1.7798e+24
- CLEAN FCONE_LZ: none
- CLEAN LZSPIN: none
### LZA-t9
- Terminating on cycle limit; rows 300 (cycles), wall 821 s, dt first/median/last 96.85 / 96.85 / 96.85 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 1.9803
- END-row residuals per cycle: Lz RMS 2.271e+27, sum -6.143e+28, max|.| 6.762e+27 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 3.914e+32, sum -7.275e+31, max|.| 1.103e+33 | M RMS 6.754e+17, sum -9.512e+18, max|.| 9.367e+18
- CLEAN LZB: 300 lines, ledger first 8.4153e+35 last 1.5131e+38 (growth 1.5047e+38); sum of Lz_hflx*dt over the arm -2.5567e+24
- CLEAN FCONE_LZ: none
- CLEAN LZSPIN: none
### LZB-t9
- Terminating on cycle limit; rows 300 (cycles), wall 822 s, dt first/median/last 96.85 / 96.85 / 96.85 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 6 elapsed intervals 1.9801
- END-row residuals per cycle: Lz RMS 2.291e+27, sum -8.899e+28, max|.| 9.690e+27 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 4.109e+32, sum 1.427e+32, max|.| 1.446e+33 | M RMS 6.866e+17, sum -8.971e+18, max|.| 9.367e+18
- CLEAN LZB: 300 lines, ledger first 7.8108e+30 last 3.4956e+38 (growth 3.4956e+38); sum of Lz_hflx*dt over the arm -5.3434e+24
- CLEAN FCONE_LZ: 602 lines, max |X| n/a; first `CLEAN FCONE_LZ cycle 0 t2s 0 momentum removed (signed, V/c^2) x1 4.332853e+31 x2 -1.447376e+18 x3 0.000000e+00`; last `CLEAN FCONE_LZ cycle 299 t2s 3 momentum removed (signed, V/c^2) x1 1.904117e+29 x2 -3.954313e+21 x3 0.000000e+00`
- CLEAN LZSPIN: 597 lines, max |X| n/a; first `CLEAN LZSPIN cycle 1 t2s 3: L_z of the cut and capped momentum X 1.782309e+30, the gas's as dOmega 2.610557e-27 (I 6.827314e+56)`; last `CLEAN LZSPIN cycle 299 t2s 3: L_z of the cut and capped momentum X -1.119682e+33, the gas's as dOmega -1.640095e-24 (I 6.826935e+56)`
### LZRA-bsg
- Terminating on cycle limit; rows 199 (cycles), wall 147 s, dt first/median/last 163.6 / 163.2 / 169 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 3 elapsed intervals 0.5830
- END-row residuals per cycle: Lz RMS 2.198e+30, sum -3.375e+32, max|.| 4.647e+30 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 3.774e+31, sum -3.545e+31, max|.| 1.075e+32 | M RMS 7.831e+16, sum 3.603e+16, max|.| 2.297e+17
- CLEAN LZB: 200 lines, ledger first -5.6851e+42 last -5.5688e+42 (growth 1.1626e+41); sum of Lz_hflx*dt over the arm -5.6144e+37
- CLEAN FCONE_LZ: none
- CLEAN LZSPIN: none
### LZRB-bsg
- Terminating on cycle limit; rows 199 (cycles), wall 148 s, dt first/median/last 163.6 / 163.2 / 169 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 3 elapsed intervals 0.5806
- END-row residuals per cycle: Lz RMS 2.213e+30, sum -3.343e+32, max|.| 4.632e+30 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 3.875e+31, sum 4.720e+31, max|.| 9.708e+31 | M RMS 8.318e+16, sum 7.656e+16, max|.| 1.801e+17
- CLEAN LZB: 200 lines, ledger first -5.6851e+42 last -5.5725e+42 (growth 1.1260e+41); sum of Lz_hflx*dt over the arm -5.6511e+37
- CLEAN FCONE_LZ: 400 lines, max |X| n/a; first `CLEAN FCONE_LZ cycle 29657 t2s 2 momentum removed (signed, V/c^2) x1 4.610215e+25 x2 -1.065701e+23 x3 0.000000e+00`; last `CLEAN FCONE_LZ cycle 29856 t2s 3 momentum removed (signed, V/c^2) x1 2.108813e+26 x2 -8.598559e+23 x3 0.000000e+00`
- CLEAN LZSPIN: 400 lines, max |X| n/a; first `CLEAN LZSPIN cycle 29657 t2s 2: L_z of the cut and capped momentum X -2.037949e+32, the gas's as dOmega -2.876832e-24 (I 7.084004e+55)`; last `CLEAN LZSPIN cycle 29856 t2s 3: L_z of the cut and capped momentum X 9.410635e+36, the gas's as dOmega 1.328497e-19 (I 7.083673e+55)`
### LZRA-m35
- Terminating on cycle limit; rows 199 (cycles), wall 224 s, dt first/median/last 5.353 / 11.96 / 11.96 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 3 elapsed intervals 0.9627
- END-row residuals per cycle: Lz RMS 1.999e+27, sum 6.330e+28, max|.| 7.213e+27 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 7.628e+29, sum -2.473e+29, max|.| 2.219e+30 | M RMS 9.781e+14, sum -3.518e+14, max|.| 2.955e+15
- CLEAN LZB: 200 lines, ledger first -1.4022e+41 last -1.4256e+41 (growth -2.3382e+39); sum of Lz_hflx*dt over the arm -1.2518e+26
- CLEAN FCONE_LZ: none
- CLEAN LZSPIN: none
### LZRB-m35
- Terminating on cycle limit; rows 199 (cycles), wall 230 s, dt first/median/last 5.353 / 11.96 / 11.96 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle median over the last 3 elapsed intervals 0.9703
- END-row residuals per cycle: Lz RMS 2.253e+27, sum 3.012e+28, max|.| 6.146e+27 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 7.396e+29, sum 3.093e+29, max|.| 2.827e+30 | M RMS 9.613e+14, sum 2.815e+14, max|.| 2.955e+15
- CLEAN LZB: 200 lines, ledger first -1.4022e+41 last -1.4180e+41 (growth -1.5774e+39); sum of Lz_hflx*dt over the arm -4.6341e+26
- CLEAN FCONE_LZ: 400 lines, max |X| n/a; first `CLEAN FCONE_LZ cycle 73280 t2s 2 momentum removed (signed, V/c^2) x1 5.174089e+24 x2 7.849905e+22 x3 0.000000e+00`; last `CLEAN FCONE_LZ cycle 73479 t2s 3 momentum removed (signed, V/c^2) x1 1.022280e+26 x2 2.394810e+23 x3 0.000000e+00`
- CLEAN LZSPIN: 400 lines, max |X| n/a; first `CLEAN LZSPIN cycle 73280 t2s 2: L_z of the cut and capped momentum X -1.132778e+34, the gas's as dOmega -6.155111e-20 (I 1.840386e+53)`; last `CLEAN LZSPIN cycle 73479 t2s 3: L_z of the cut and capped momentum X -5.667343e+34, the gas's as dOmega -3.079420e-19 (I 1.840393e+53)`
### DFL0
- Terminating on cycle limit; rows 99 (cycles), wall 89 s, dt first/median/last 163.6 / 162.8 / 159.6 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle (no elapsed lines at this log cadence): cpu time used 6.658633e+01 s / 99 cycles = 0.6726
- END-row residuals per cycle: Lz RMS 7.603e+29, sum -4.767e+31, max|.| 2.087e+30 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 4.124e+31, sum -5.309e+31, max|.| 1.075e+32 | M RMS 7.865e+16, sum 1.441e+17, max|.| 1.937e+17
- CLEAN LZB: 100 lines, ledger first -5.6851e+42 last -5.5826e+42 (growth 1.0248e+41); sum of Lz_hflx*dt over the arm -2.8105e+37
- CLEAN FCONE_LZ: none
- CLEAN LZSPIN: none
- Lz_hflx / E_hflx / M_hflx means -1.7508e+33 / 7.4264e+27 / 1.1398e+07; last 10 rows:
  - t 4.939577e+06: Lz_hflx -1.7512e+33 E_hflx 7.4058e+27 M_hflx 1.7597e+07
  - t 4.939737e+06: Lz_hflx -1.7507e+33 E_hflx 7.4054e+27 M_hflx 1.7959e+07
  - t 4.939896e+06: Lz_hflx -1.7500e+33 E_hflx 7.4050e+27 M_hflx 1.8321e+07
  - t 4.940056e+06: Lz_hflx -1.7493e+33 E_hflx 7.4046e+27 M_hflx 1.8676e+07
  - t 4.940216e+06: Lz_hflx -1.7489e+33 E_hflx 7.4042e+27 M_hflx 1.9040e+07
  - t 4.940376e+06: Lz_hflx -1.7481e+33 E_hflx 7.4038e+27 M_hflx 1.9394e+07
  - t 4.940536e+06: Lz_hflx -1.7474e+33 E_hflx 7.4034e+27 M_hflx 1.9740e+07
  - t 4.940695e+06: Lz_hflx -1.7469e+33 E_hflx 7.4030e+27 M_hflx 2.0077e+07
  - t 4.940855e+06: Lz_hflx -1.7463e+33 E_hflx 7.4026e+27 M_hflx 2.0425e+07
  - t 4.941014e+06: Lz_hflx -1.7455e+33 E_hflx 7.4022e+27 M_hflx 2.0753e+07
### DFL1
- Terminating on cycle limit; rows 99 (cycles), wall 95 s, dt first/median/last 163.6 / 162.8 / 159.6 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle (no elapsed lines at this log cadence): cpu time used 6.717507e+01 s / 99 cycles = 0.6785
- END-row residuals per cycle: Lz RMS 8.487e+29, sum -5.009e+31, max|.| 2.309e+30 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 4.189e+31, sum 2.662e+31, max|.| 1.225e+32 | M RMS 8.111e+16, sum 1.441e+17, max|.| 3.198e+17
- CLEAN LZB: 100 lines, ledger first -5.6851e+42 last -5.5826e+42 (growth 1.0247e+41); sum of Lz_hflx*dt over the arm -5.3920e+36
- CLEAN FCONE_LZ: none
- CLEAN LZSPIN: none
- Lz_hflx / E_hflx / M_hflx means -3.3581e+32 / 1.8313e+27 / 1.1345e+07; last 10 rows:
  - t 4.939577e+06: Lz_hflx -3.7290e+32 E_hflx 1.8383e+27 M_hflx 1.7456e+07
  - t 4.939737e+06: Lz_hflx -3.6750e+32 E_hflx 1.8384e+27 M_hflx 1.7814e+07
  - t 4.939896e+06: Lz_hflx -3.7226e+32 E_hflx 1.8374e+27 M_hflx 1.8169e+07
  - t 4.940056e+06: Lz_hflx -3.6736e+32 E_hflx 1.8369e+27 M_hflx 1.8517e+07
  - t 4.940216e+06: Lz_hflx -3.7183e+32 E_hflx 1.8372e+27 M_hflx 1.8863e+07
  - t 4.940376e+06: Lz_hflx -3.7230e+32 E_hflx 1.8373e+27 M_hflx 1.9218e+07
  - t 4.940536e+06: Lz_hflx -3.6980e+32 E_hflx 1.8369e+27 M_hflx 1.9569e+07
  - t 4.940695e+06: Lz_hflx -3.6709e+32 E_hflx 1.8372e+27 M_hflx 1.9905e+07
  - t 4.940855e+06: Lz_hflx -3.7069e+32 E_hflx 1.8368e+27 M_hflx 2.0243e+07
  - t 4.941014e+06: Lz_hflx -3.5864e+32 E_hflx 1.8383e+27 M_hflx 2.0560e+07
### PLED1
- Terminating on cycle limit; rows 99 (cycles), wall 802 s, dt first/median/last 163.6 / 162.8 / 159.6 s; CW 0 FATAL 0 DT_RETRY 0
- s/cycle (no elapsed lines at this log cadence): cpu time used 7.741968e+02 s / 99 cycles = 7.8202
- END-row residuals per cycle: Lz RMS 7.603e+29, sum -4.767e+31, max|.| 2.087e+30 | E (dEtot + E_hflx dt - (L_in - L_top) dt) RMS 4.124e+31, sum -5.309e+31, max|.| 1.075e+32 | M RMS 7.865e+16, sum 1.441e+17, max|.| 1.937e+17
- CLEAN LZB: 100 lines, ledger first -5.6851e+42 last -5.5826e+42 (growth 1.0248e+41); sum of Lz_hflx*dt over the arm -2.8105e+37
- CLEAN FCONE_LZ: none
- CLEAN LZSPIN: none
- p_ledger columns in user.hst: ['Pr_tot', 'Pr_hflx', 'Pr_ntel', 'Pr_src', 'Pr_rad', 'Pr_flr', 'Pr_mrg', 'Pr_gap'] (per cycle in the committed user.hst); means: Pr_tot -1.8727e+34, Pr_hflx -1.7923e+34, Pr_ntel -7.3969e+16, Pr_src -1.7966e+34, Pr_rad 4.2492e+31, Pr_flr 2.6992e+26, Pr_mrg 8.6423e+15, Pr_gap 0.0000e+00
- RMS of Pr_mrg 1.3607e+16, Pr_rad 4.2492e+31
- CLEAN PLEDGER END lines: 100; last `CLEAN PLEDGER cycle 29756 dt 1.595807e+02 hflx -1.792319e+34 ntel -1.849521e+18 src -1.796656e+34 rad 4.269058e+31 flr 2.691733e+26 mrg 0.000000e+00 gap 0.000000e+00 P -2.6891461970346939e+34 END -5.372614e+20 gross 1.098050e+36`
```

## Short reading (you score from the files)
- **LZ keys, B vs A:**
  - m35mms fresh: A and B are bitwise identical (dumps and every residual). No FCONE_LZ / LZSPIN line appears in either, so the cut never acts there.
  - bsg fresh: the cut acts on cycles 0-17 only (34 FCONE_LZ, 3 LZSPIN lines). L_z residual RMS is 9.5e29 vs 9.7e29. The ledger grows by -1.182e34 in both arms; the sum of Lz_hflx*dt is only about 1e30, so the ledger growth is not the boundary torque.
  - t9l62 fresh: the cut acts every stage (602 FCONE_LZ, 597 LZSPIN lines). The ledger grows 3.50e38 in B against 1.50e38 in A.
  - R arms (developed convection): the ledger grows +1.13e41 (B) vs +1.16e41 (A) on bsg, and -1.58e39 vs -2.34e39 on m35mms. Again the boundary torque is far smaller (Lz_hflx*dt sum -5.7e37 and -4.6e26).
  - So in my reading the ledger in B is not just the radial boundaries' torque. Please check the ledger definition against these numbers.
- **Physics B vs A at the last cycle:** fresh bsg/t9: rho <= 3.5e-3 / 5.3e-5 and T <= 4e-5 / 2.3e-5 where rho > 10 dfloor. The R arms differ more after 200 cycles of developed convection: bsg rho 0.25 and m35mms rho 16 at single cells, means 1.5e-5 / 4.7e-3. This is the chaotic growth of a small change, not a bias.
- **DFL:** the sum of Lz_hflx*dt over 100 cycles falls from -2.81e37 (DFL0) to -5.39e36 (DFL1), and E_hflx / M_hflx follow (rows in the readings). Physics DFL1 vs DFL0: rho max 2.1e-2 and T max 4.0e-3 where rho > 10 dfloor.

## Item 5: physics difference at the last cycle (B vs A; DFL1 vs DFL0)
- LZB-bsg vs LZA-bsg (cycles (300, 300), dfloor 6.152e-14):
  - rho: all cells max 1.128e-02 at (r 3.7992e+12, th 1.8878, ph 0.1779), mean 7.162e-06; rho>10 dfloor max 3.522e-03 at (r 3.7774e+12, th 1.9328, ph 0.1779), mean 4.068e-06
  - T: all cells max 4.146e-05 at (r 3.7448e+12, th 1.9410, ph 0.5216), mean 3.724e-07; rho>10 dfloor max 4.146e-05 at (r 3.7448e+12, th 1.9410, ph 0.5216), mean 1.656e-07
  - |v|: all cells max 4.563e-02 at (r 5.4243e+12, th 1.2333, ph 0.7629), mean 4.369e-03; rho>10 dfloor max 1.355e-03 at (r 3.7774e+12, th 1.8837, ph 0.4316), mean 1.740e-06
- LZB-m35 vs LZA-m35 (cycles (300, 300), dfloor 1.000e-17):
  - rho: all cells max 0.000e+00 at (r 6.7509e+11, th 1.3969, ph 0.0007), mean 0.000e+00; rho>10 dfloor max 0.000e+00 at (r 6.7509e+11, th 1.3969, ph 0.0007), mean 0.000e+00
  - T: all cells max 0.000e+00 at (r 6.7509e+11, th 1.3969, ph 0.0007), mean 0.000e+00; rho>10 dfloor max 0.000e+00 at (r 6.7509e+11, th 1.3969, ph 0.0007), mean 0.000e+00
  - |v|: all cells max 0.000e+00 at (r 6.7509e+11, th 1.3969, ph 0.0007), mean 0.000e+00; rho>10 dfloor max 0.000e+00 at (r 6.7509e+11, th 1.3969, ph 0.0007), mean 0.000e+00
- LZB-t9 vs LZA-t9 (cycles (300, 300), dfloor 3.600e-19):
  - rho: all cells max 1.147e-03 at (r 6.9870e+12, th 2.0954, ph 2.9729), mean 6.465e-08; rho>10 dfloor max 5.349e-05 at (r 6.8918e+12, th 1.9359, ph 0.0644), mean 3.950e-08
  - T: all cells max 6.833e-05 at (r 2.1208e+13, th 0.7885, ph 2.7949), mean 5.893e-06; rho>10 dfloor max 2.267e-05 at (r 6.5690e+12, th 1.0830, ph 2.3102), mean 2.854e-07
  - |v|: all cells max 6.816e-02 at (r 5.6130e+13, th 0.7885, ph 0.1687), mean 5.375e-04; rho>10 dfloor max 3.587e-05 at (r 6.8918e+12, th 2.2488, ph 0.2608), mean 2.502e-08
- LZRB-bsg vs LZRA-bsg (cycles (29857, 29857), dfloor 6.152e-14):
  - rho: all cells max 2.458e-01 at (r 3.6470e+12, th 2.0923, ph 0.6034), mean 2.948e-05; rho>10 dfloor max 2.458e-01 at (r 3.6470e+12, th 2.0923, ph 0.6034), mean 1.474e-05
  - T: all cells max 4.573e-02 at (r 3.5926e+12, th 1.7324, ph 0.6525), mean 6.267e-05; rho>10 dfloor max 4.573e-02 at (r 3.5926e+12, th 1.7324, ph 0.6525), mean 7.706e-06
  - |v|: all cells max 1.205e-01 at (r 3.6470e+12, th 1.7610, ph 0.5256), mean 8.612e-05; rho>10 dfloor max 5.600e-02 at (r 3.6470e+12, th 2.0923, ph 0.6034), mean 3.249e-06
- LZRB-m35 vs LZRA-m35 (cycles (73480, 73480), dfloor 1.000e-17):
  - rho: all cells max 1.614e+01 at (r 9.7895e+11, th 1.5551, ph 0.3143), mean 5.611e-03; rho>10 dfloor max 1.605e+01 at (r 9.7763e+11, th 1.5538, ph 0.3143), mean 4.735e-03
  - T: all cells max 3.558e-01 at (r 9.5405e+11, th 1.3969, ph 0.3334), mean 1.137e-04; rho>10 dfloor max 3.558e-01 at (r 9.5405e+11, th 1.3969, ph 0.3334), mean 6.701e-05
  - |v|: all cells max 5.198e-01 at (r 9.6970e+11, th 1.5101, ph 0.1507), mean 2.233e-03; rho>10 dfloor max 5.198e-01 at (r 9.6970e+11, th 1.5101, ph 0.1507), mean 6.872e-04
- DFL1 vs DFL0 (cycles (29757, 29757), dfloor 6.152e-14):
  - rho: all cells max 3.609e-02 at (r 3.7557e+12, th 1.8919, ph 0.6974), mean 1.420e-06; rho>10 dfloor max 2.134e-02 at (r 3.6361e+12, th 1.8019, ph 0.7343), mean 7.615e-07
  - T: all cells max 1.576e-02 at (r 3.5926e+12, th 1.3438, ph 0.5911), mean 4.911e-06; rho>10 dfloor max 3.973e-03 at (r 3.6253e+12, th 1.8019, ph 0.7343), mean 5.149e-07
  - |v|: all cells max 5.478e-02 at (r 3.7122e+12, th 1.1801, ph 0.7220), mean 8.789e-05; rho>10 dfloor max 4.208e-03 at (r 3.6035e+12, th 1.3519, ph 0.6156), mean 1.580e-07

## Byte checks
- LZD-bsg vs LZB-bsg: user.hst byte-identical **yes**; final hydro_w data (after the parameter block) **yes**. The whole file differs only in the header, which holds the command line with the 2 named keys.
- PLED1 vs DFL0 hydro.hst: byte-identical **yes** (md5 0f10b0ca).

## RST: (a) 2 cycles straight vs (b) 1 cycle + rst + 1 cycle, compared after the parameter block
- RST-DFL-bsg: CW 0 FATAL 0; final rst bsg3d.00010.rst vs bsg3d.00011.rst identical: yes; hydro_w bsg3d.hydro_w.00117.bin vs bsg3d.hydro_w.00117.bin identical: yes
- RST-DFL-m35mms: CW 0 FATAL 0; final rst m35mms.00012.rst vs m35mms.00013.rst identical: yes; hydro_w m35mms.hydro_w.00483.bin vs m35mms.hydro_w.00483.bin identical: yes
- RST-DFL-m35zams: CW 0 FATAL 0; final rst m35zams.00012.rst vs m35zams.00013.rst identical: yes; hydro_w m35zams.hydro_w.00483.bin vs m35zams.hydro_w.00483.bin identical: yes
- RST-DFL-t9l62: CW 0 FATAL 0; final rst t9l62.00013.rst vs t9l62.00014.rst identical: yes; hydro_w t9l62.hydro_w.00252.bin vs t9l62.hydro_w.00252.bin identical: yes
- RST-LZD-bsg: CW 0 FATAL 0; final rst bsg3d.00010.rst vs bsg3d.00011.rst identical: yes; hydro_w bsg3d.hydro_w.00117.bin vs bsg3d.hydro_w.00117.bin identical: yes
- RST-LZD-m35mms: CW 0 FATAL 0; final rst m35mms.00012.rst vs m35mms.00013.rst identical: yes; hydro_w m35mms.hydro_w.00483.bin vs m35mms.hydro_w.00483.bin identical: yes
- RST-LZD-m35zams: CW 0 FATAL 0; final rst m35zams.00012.rst vs m35zams.00013.rst identical: yes; hydro_w m35zams.hydro_w.00483.bin vs m35zams.hydro_w.00483.bin identical: yes
- RST-LZD-t9l62: CW 0 FATAL 0; final rst t9l62.00013.rst vs t9l62.00014.rst identical: yes; hydro_w t9l62.hydro_w.00252.bin vs t9l62.hydro_w.00252.bin identical: yes
- RST-PLED-bsg: CW 0 FATAL 0; final rst bsg3d.00010.rst vs bsg3d.00011.rst identical: yes; hydro_w bsg3d.hydro_w.00117.bin vs bsg3d.hydro_w.00117.bin identical: yes
- All 9 pairs: final rst and last hydro_w identical, 0 CW, 0 FATAL.

## Files
- docs/handover/caltech_results/LZAB/<arm>/: user.hst, hydro.hst, log first 400 / last 200 lines, lz_lines.txt (CLEAN LZB/FCONE/LZSPIN/PLEDGER), clean_warning_count.txt, run.log. The RST arms have run.log + CW count.
