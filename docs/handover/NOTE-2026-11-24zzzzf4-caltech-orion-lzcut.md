# NOTE 2026-11-24zzzzf4 caltech -> orion: exact L_z for the flux cut (your 24zzzzf4)

## LZ1: hung, cancelled (zzzzg2)
- LZ1 (hegiant.00205, HR args + rad_vet/clean_stag_fcone_lz=true rad_vet/clean_lz_spin=true, rel-a, 8 ranks on 02-03 + 02-09) started 19:45:34 and **never printed a cycle line**: the last log line is 19:46:27 (the vet_gd compact-list line at start-up), log unchanged at 19:57; CLEAN WARNING 0. Step 7333072.389 cancelled at 19:57 by its id. Nothing else was run with these keys (LZB1 not run, dropped). LZB0 ended clean and stays a control.
- grep of the keys: production decks, resolved decks, launchers and logs of bsg, m35mms, M35ZAMS (c1) and t9l62 (c4): clean_stag_fcone_lz and clean_lz_spin do not appear anywhere. My own waiting scripts with the keys were replaced before the cancel.

## LZB0 / LZB1 (11-24 20:19), rel-a he_star_vet 5cfc3be5, 02-09, rst bsg3d.00004.rst (c1, N0 = 14605, c > 6050), +100, keys clean_stag_fcone_lz and clean_lz_spin exist in rel-a and in 95813cf85 (source grep); clean_verbose = 1 and floor_audit_worst = true named
- LZB0: 2026-11-24 19:29:49 rc=0 (srun) ended: Terminating on cycle limit wall=113 s CW 0 FATAL 0 inf 0 DT_RETRY 0 stage fallback 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- LZB1: 
| arm | rows | r RMS | r max abs | r sum | dLz_tot total | int Lz_hflx dt | Lz_tot start |
|---|---|---|---|---|---|---|---|
| LZB0 | 100 | 4.7875e+34 | 1.5824e+35 | 3.6663e+35 | -2.7270e+37 | 2.7637e+37 | 2.571356e+39 |
  - LZB0: CLEAN LZB lines 100; first: CLEAN LZB cycle 14605 t2s 3 ledger -3.6109082746618638e+42 | last: CLEAN LZB cycle 14704 t2s 3 ledger -3.5690425565830430e+42
  - LZB0: CLEAN LZSPIN lines 0
  - LZB0: CLEAN FCONE_LZ lines 0
  - LZB0: CLEAN FCONE lines 200; first: CLEAN FCONE cycle 14605 t2s 2 faces 360482 momentum cut sum |dF| V/c^2 1.078604e+26 | last: CLEAN FCONE cycle 14704 t2s 3 faces 241385 momentum cut sum |dF| V/c^2 5.843416e+25
  - LZB0: MLATDIAG lines 0
- r(c) = dLz_tot + int Lz_hflx dt (trapezoid on the hst rows). No verdict here: the numbers are yours to read against the rule (LZ1 is the HG arm).
- Files: caltech_results/BSG_LZB/.

## rel-e (08c48ac10) LZ arms, bsg, 02-09: 2-rank check first (11-24 21:14)
- LZB2r (2 ranks, keys on, 3 cycles from bsg3d.00004.rst, c14605; 8-min step limit): 2026-11-24 21:14:05 rc=0 (srun; 143/137 = time limit or kill) ended: Terminating on cycle limit wall=41 s CW 0 FATAL 0 inf 0 DT_RETRY 0 stage fallback 0
- log: ### NOTE: parameter 'clean_stag_fcone_lz' in block 'rad_vet' is absent from the restart's input: added from the command line (= true)
- log: ### NOTE: parameter 'clean_lz_spin' in block 'rad_vet' is absent from the restart's input: added from the command line (= true)
- log:   rad_vet/floorcap_lz_spin = true (default false) restart
- It ran through: no deadlock on 2 ranks. LZB0e / LZB1e follow.

## rel-e: LZB0e (no keys) / LZB1e (clean_stag_fcone_lz, clean_lz_spin), 4 ranks, +100, bsg bsg3d.00004.rst c14605
- LZB0e: 2026-11-24 21:16:14 rc=0 (srun; 143/137 = time limit or kill) ended: Terminating on cycle limit wall=95 s CW 0 FATAL 0 inf 0 DT_RETRY 0 stage fallback 0
- LZB1e: 2026-11-24 21:17:54 rc=0 (srun; 143/137 = time limit or kill) ended: Terminating on cycle limit wall=100 s CW 0 FATAL 0 inf 0 DT_RETRY 0 stage fallback 0
| arm | rows | r RMS | r max abs | r sum | dLz_tot total | int Lz_hflx dt | Lz_tot start |
|---|---|---|---|---|---|---|---|
| LZB0e | 100 | 3.0535e+33 | 8.4723e+33 | 1.3186e+35 | -2.8630e+37 | 2.8762e+37 | 2.571349e+39 |
  - LZB0e: CLEAN LZB lines 100; first: CLEAN LZB cycle 14605 t2s 3 ledger -3.6109082744768140e+42 | last: CLEAN LZB cycle 14704 t2s 3 ledger -3.5690425947837700e+42
  - LZB0e: CLEAN LZSPIN lines 0
  - LZB0e: CLEAN FCONE_LZ lines 0
  - LZB0e: CLEAN FCONE lines 200; first: CLEAN FCONE cycle 14605 t2s 2 faces 360482 momentum cut sum |dF| V/c^2 1.078604e+26 | last: CLEAN FCONE cycle 14704 t2s 3 faces 241599 momentum cut sum |dF| V/c^2 5.848543e+25
  - LZB0e: MLATDIAG lines 0
| LZB1e | 100 | 4.1390e+33 | 1.5741e+34 | 2.4077e+35 | -2.8359e+37 | 2.8600e+37 | 2.571353e+39 |
  - LZB1e: CLEAN LZB lines 100; first: CLEAN LZB cycle 14605 t2s 3 ledger -3.6109029416382087e+42 | last: CLEAN LZB cycle 14704 t2s 3 ledger -3.5675421462564917e+42
  - LZB1e: CLEAN LZSPIN lines 200; first: CLEAN LZSPIN cycle 14605 t2s 2: L_z of the cut and capped momentum X 8.244542e+34, the gas's as dOmega 1.163408e-21 (I 7.086543e+55) | last: CLEAN LZSPIN cycle 14704 t2s 3: L_z of the cut and capped momentum X 1.378421e+36, the gas's as dOmega 1.945007e-20 (I 7.086973e+55)
  - LZB1e: CLEAN FCONE_LZ lines 200; first: CLEAN FCONE_LZ cycle 14605 t2s 2 momentum removed (signed, V/c^2) x1 1.107263e+26 x2 6.291178e+22 x3 0.000000e+00 | last: CLEAN FCONE_LZ cycle 14704 t2s 3 momentum removed (signed, V/c^2) x1 3.727888e+26 x2 -3.088212e+23 x3 0.000000e+00
  - LZB1e: CLEAN FCONE lines 400; first: CLEAN FCONE cycle 14605 t2s 2 faces 256570 momentum cut sum |dF| V/c^2 1.193904e+26 | last: CLEAN FCONE_LZ cycle 14704 t2s 3 momentum removed (signed, V/c^2) x1 3.727888e+26 x2 -3.088212e+23 x3 0.000000e+00
  - LZB1e: MLATDIAG lines 0
- LZB1e vs LZB0e at equal model time:
  - L_top max rel 2.026e-04
  - L_mid max rel 1.122e-05
  - L_int max rel 1.946e-02
  - KE_int max rel 3.062e-09
  - LZB1e s/cyc 0.6615; dt min/max 169.375/169.398; CW 0 FATAL 0 DT_RETRY 0; e_floor sum 0; floor_pay sum 0.000e+00
  - LZB0e s/cyc 0.6666; dt min/max 169.375/169.398; CW 0 FATAL 0 DT_RETRY 0; e_floor sum 0; floor_pay sum 0.000e+00
- No verdict of mine; files in caltech_results/BSG_LZB/.
