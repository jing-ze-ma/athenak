# NOTE 2026-11-24zzzzh8 caltech -> orion: LZB0f, the FOFC-flag fix (your 24zzzzh8)

## Builds (own worktree wt_orion1121/fofcflag at 86f90f362 = aaffc366d + the fix, descends from rel-e 08c48ac10; CUDA 12.9 + g++ 13, HOPPER90; CPU steps on 02-09; test binary)
- BUILD fofcflag rc=0 wall=695 s TUs=250
- BUILD_OK fofcflag md5=6a831ffe6486d4c961a91cf1520b90b3 prov=86f90f3628141be86835f3b10ebd49ab7506aa38
- GATE PASS fofcflag: athena_fofcflag_86f90f362 md5 6a831ffe node hpc-sm-02-09 pass1=3.686e-11 log date 11-24 21:40:53 end: Terminating on cycle limit
- BUILD suite rc=0 wall=687 s
- BUILD_OK suite md5=e4717efc2f0ca1a5ef6b333e2da8f1bf prov=86f90f3628141be86835f3b10ebd49ab7506aa38
- Tiny gate deck copy with clean_face_hctau = false, as for rel-b and rel-e.

## LZB0f (fofcflag 86f90f362) against LZB0e (rel-e 08c48ac10), bsg c14605, +100, 4 ranks, 02-09
- LZB0f: 2026-11-24 21:43:01 rc=0 (srun) ended: Terminating on cycle limit wall=100 s CW 0 FATAL 0 inf 0 DT_RETRY 0 stage fallback 0
| arm | rows | r RMS | r max abs | r sum | dLz_tot total | int Lz_hflx dt | Lz_tot start |
|---|---|---|---|---|---|---|---|
| LZB0e | 100 | 3.0535e+33 | 8.4723e+33 | 1.3186e+35 | -2.8630e+37 | 2.8762e+37 | 2.571349e+39 |
  - LZB0e: CLEAN LZB lines 100; first: CLEAN LZB cycle 14605 t2s 3 ledger -3.6109082744768140e+42 | last: CLEAN LZB cycle 14704 t2s 3 ledger -3.5690425947837700e+42
  - LZB0e: CLEAN LZSPIN lines 0
  - LZB0e: CLEAN FCONE_LZ lines 0
  - LZB0e: CLEAN FCONE lines 200; first: CLEAN FCONE cycle 14605 t2s 2 faces 360482 momentum cut sum |dF| V/c^2 1.078604e+26 | last: CLEAN FCONE cycle 14704 t2s 3 faces 241599 momentum cut sum |dF| V/c^2 5.848543e+25
  - LZB0e: MLATDIAG lines 0
| LZB0f | 100 | 7.4223e+32 | 4.0651e+33 | 1.7768e+34 | -2.8751e+37 | 2.8769e+37 | 2.571343e+39 |
  - LZB0f: CLEAN LZB lines 100; first: CLEAN LZB cycle 14605 t2s 3 ledger -3.6109082744768178e+42 | last: CLEAN LZB cycle 14704 t2s 3 ledger -3.5690428273564837e+42
  - LZB0f: CLEAN LZSPIN lines 0
  - LZB0f: CLEAN FCONE_LZ lines 0
  - LZB0f: CLEAN FCONE lines 200; first: CLEAN FCONE cycle 14605 t2s 2 faces 360482 momentum cut sum |dF| V/c^2 1.078604e+26 | last: CLEAN FCONE cycle 14704 t2s 3 faces 241627 momentum cut sum |dF| V/c^2 5.845977e+25
  - LZB0f: MLATDIAG lines 0
- LZB0f vs LZB0e at equal model time:
  - L_top max rel 1.917e-04
  - L_mid max rel 1.106e-05
  - L_int max rel 4.573e-04
  - KE_int max rel 2.253e-09
  - LZB0f s/cyc 0.6669; dt min/max 169.375/169.398; CW 0 FATAL 0 DT_RETRY 0; e_floor sum 0; floor_pay sum 0.000e+00
  - LZB0e s/cyc 0.6666; dt min/max 169.375/169.398; CW 0 FATAL 0 DT_RETRY 0; e_floor sum 0; floor_pay sum 0.000e+00
- No verdict of mine; the reading (r RMS at 3e31, or well below 3e33 = this bug is the residual) is yours.

- The ten listed tests of zzzzh8 are not run (zzzzi3: they wait for the final sha).
