# NOTE 2026-11-24zzzzc2 caltech -> orion: nsys profiles and the t9l62 merged-base arms (your 24zzzzc2), sdorder3 30d01765

## PRF_t9q (11-24 17:42)
TAG=PRF_t9q BIN md5=30d01765970051212176cc2eb81ea034 RST=/resnick/groups/carnegie_poc/jingze/prod_rc4v_1124/t9l62/rst_keep/t9l62.stop1124b.rst N0=2553 NC=12 HOLDER=7333072 NODE=hpc-sm-02-03 NR=4 NSYS=/central/software9/spack/opt/spack/linux-rhel9-x86_64/none-n
run rc=0 (srun) ended: Terminating on cycle limit wall=122 s CW 0 FATAL 0
main loop time (before the final outputs) = 31.3736
nsys_r0.sqlite: cycles 11 wall 2413.0 ms/cyc  busy 2244.5  kernel sum 2244.5  launches 5361/cyc  keys 130
  ms/cyc     %   n/cyc    us/l  reg  blk     grid lmem  label | key
  463.97  20.7     646   718.2  130  128    16128    0   | radvet::RadiationVet::VetGdSweep#1
nsys_r0.sqlite: cycles 11 wall 2413.0 ms/cyc  GPU idle 141.5 ms/cyc (5.9 %)  gaps 6346/cyc
idle by host activity (ms/cyc): MPI 75.2  MPI(part) 31.5  sync 9.9  memcpy-api(part) 9.1  sync(part) 4.7  launch(part) 4.7  launch 2.5  memcpy-api 2.0  cuda-other 1.2  cuda-other(part) 0.5  host-compute 0.1
idle by gap size   (ms/cyc, n/cyc): <10us 13.6 (5008)  10-100us 22.9 (1304)  0.1-1ms 6.4 (25)  >1ms 98.5 (9)
memcpy (kind: n/cyc, MB/cyc, ms/cyc; 1=HtoD 2=DtoH 8=DtoD 10=PtoP): 1: 632, 13.9, 1.10  2: 80, 13.0, 0.43  8: 21, 7743.6, 3.71  10: 1984, 5339.7, 50.06
host time in (ms/cyc, overlapping GPU work or not): sync 2085.7  MPI 312.5
largest gaps (us, host activity, previous kernel -> next):
    57863.3 MPI(part)      M -> M
    42174.4 MPI            M -> M
    40439.1 MPI            M -> M
    38022.4 MPI            M -> M
    36232.0 MPI            M -> M
    35948.0 MPI            M -> M
- Text files: caltech_results/PRF_1124/PRF_t9q/ (res.txt, kname_r0..3.txt, gaps_r0..3.txt, stats_r0.txt, CLEAN PASS/cycle lines, hst). Plain s/cyc of the same state in the CAND arm: 2.582 s/cyc (t9l62_R, b8b0d79c). Only the kernel SHARES of this run count (nsys slows the step). nsys from cuda-12.9.0.

## PRF_m35 (11-24 17:45)
TAG=PRF_m35 BIN md5=30d01765970051212176cc2eb81ea034 RST=/resnick/groups/carnegie_poc/jingze/prod_rc4v_1124/m35mms/rst_keep/m35mms.stop1124b.rst N0=6959 NC=12 HOLDER=7333072 NODE=hpc-sm-02-03 NR=4 NSYS=/central/software9/spack/opt/spack/linux-rhel9-x86_64/none
run rc=0 (srun) ended: Terminating on cycle limit wall=70 s CW 0 FATAL 0
main loop time (before the final outputs) = 16.3626
nsys_r0.sqlite: cycles 11 wall 1207.5 ms/cyc  busy 996.4  kernel sum 996.4  launches 6671/cyc  keys 133
  ms/cyc     %   n/cyc    us/l  reg  blk     grid lmem  label | key
   99.00   9.9    1035    95.7   61  352     8192    0   | radvet::RadiationVet::CleanStagLineSolve#1
nsys_r0.sqlite: cycles 11 wall 1207.5 ms/cyc  GPU idle 187.6 ms/cyc (15.5 %)  gaps 9418/cyc
idle by host activity (ms/cyc): MPI 119.8  MPI(part) 22.3  sync 18.7  memcpy-api(part) 9.3  sync(part) 5.7  memcpy-api 3.8  cuda-other 2.5  launch 2.1  launch(part) 2.1  cuda-other(part) 1.2  host-compute 0.2
idle by gap size   (ms/cyc, n/cyc): <10us 18.6 (6997)  10-100us 44.7 (2386)  0.1-1ms 5.3 (20)  >1ms 119.0 (15)
memcpy (kind: n/cyc, MB/cyc, ms/cyc; 1=HtoD 2=DtoH 8=DtoD 10=PtoP): 1: 1858, 9.9, 2.14  2: 262, 7.8, 0.71  8: 21, 2556.4, 1.27  10: 1866, 4721.0, 44.71
host time in (ms/cyc, overlapping GPU work or not): sync 851.8  MPI 314.4
largest gaps (us, host activity, previous kernel -> next):
    33573.1 MPI            M -> M
    32946.4 MPI            M -> M
    32260.0 MPI            M -> M
    32015.6 MPI            M -> M
    31798.8 MPI            M -> M
    31602.6 MPI            M -> M
- Text files: caltech_results/PRF_1124/PRF_m35/ (res.txt, kname_r0..3.txt, gaps_r0..3.txt, stats_r0.txt, CLEAN PASS/cycle lines, hst). Plain s/cyc of the same state in the CAND arm: 1.097 s/cyc (m35mms_R; m35mms state = stop1124b c6959, the file that arm started from). Only the kernel SHARES of this run count (nsys slows the step). nsys from cuda-12.9.0.

## MLQ arms: t9l62 quiet stop rst (c2553), candidate keys, merge on (MLQ1) vs off (MLQ0) (11-24 18:04)
- **USER decision to make: the merge (hydro/merge_lat_rc 2.022e12, f 2) changes the lateral resolution of the innermost 27 of 512 shells (2 x 2 lateral cells merged). A clean run is not a claim that the answer is unchanged.**
| arm | cycles (from the rst) | rows | dt min/mean/max | its t2s2 (sum, solves) | t2s3 |
|---|---|---|---|---|---|
| MLQ0 | 0-11 | 11 | 96.84/96.84/96.84 | 5.6, 2.09 | 5.8, 2.00 |
| MLQ0 | 11-100 | 89 | 96.84/96.84/96.84 | 5.4, 3.00 | 9.6, 2.57 |
| MLQ1 | 0-11 | 11 | 96.84/147.45/152.51 | 6.3, 2.00 | 6.3, 2.00 |
| MLQ1 | 11-100 | 89 | 152.51/152.51/152.51 | 6.0, 2.00 | 8.5, 2.20 |
- MLQ0: loop 255.7 s over 100 cycles = 2.5571 s/cyc mean; s per 1000 s of model time 26.672; CW 0 FATAL 0 DT_RETRY 0; BE 0 fallback 0 retry 0; floor_pay 4.648e+27; L_top/L_in end 1.015120; t end 2.569318e+05
  - DTREP: 20 lines; binders (gid, kji, dir: count) {('0', '54,3,3', '3'): 7, ('0', '52,3,3', '3'): 4, ('0', '57,3,3', '3'): 4}; last r 1.690906e+12 rho 9.981467e-05 v -3.231816e+01 ceff 2.279700e+07 h 7.358921e+09 cfl*dt_hyd 9.684052e+01
  - DTREG polar: 20 lines; binders (gid, kji, dir: count) {(None, None, '3'): 20}; last r 1.690906e+12 rho 9.981467e-05 v -3.231816e+01 ceff 2.279700e+07 h 7.358921e+09 cfl*dt_hyd 9.684052e+01
  - DTREG rest: 20 lines; binders (gid, kji, dir: count) {(None, None, '3'): 20}; last r 1.690906e+12 rho 9.981465e-05 v -3.313580e+01 ceff 2.279700e+07 h 7.666177e+09 cfl*dt_hyd 1.008839e+02
- MLQ1: loop 239.7 s over 100 cycles = 2.3971 s/cyc mean; s per 1000 s of model time 15.876; CW 0 FATAL 0 DT_RETRY 0; BE 0 fallback 0 retry 0; floor_pay 9.164e+27; L_top/L_in end 1.017342; t end 2.624432e+05
  - DTREP: 20 lines; binders (gid, kji, dir: count) {('7', '10,93,3', '3'): 4, ('7', '6,97,3', '3'): 3, ('5', '124,101,3', '3'): 2}; last r 1.69091e+12 rho 9.98146e-05 v -2.7533 ceff 2.2797e+07 h 9.13513e+09 cfl*dt_hyd 152.511
  - DTREG polar: 20 lines; binders (gid, kji, dir: count) {(None, None, '3'): 20}; last r 1.690906e+12 rho 9.981451e-05 v -2.506798e+00 ceff 2.279705e+07 h 7.358921e+09 cfl*dt_hyd 9.684044e+01
  - DTREG rest: 20 lines; binders (gid, kji, dir: count) {(None, None, '3'): 20}; last r 1.690906e+12 rho 9.981452e-05 v -1.352055e+00 ceff 2.279705e+07 h 7.666177e+09 cfl*dt_hyd 1.008838e+02
- MLATDIAG (every 10 cycles; 20 lines): relative drift of M 2.22e-16, E 2.01e-06; |dP|/|P0| 2.13e+00; |dL|/|L0| 4.76e-01; cycle step 10 max M 3.09e-16 E 2.78e-07
- MLATDIAG prints the totals only of the arm that names merge_lat_diag (MLQ1); MLQ0 has none, so "vs MLQ0" drifts come from the hst below.

| cycles of MLQ0 (from the rst) | t | L_top/L_in MLQ0 / MLQ1 | rel L_top | rel L_bot | rel KE_int | rel Etot | rel M_tot |
|---|---|---|---|---|---|---|---|
| 2 | 2.47538e+05 | 1.014099 / 1.014127 | 2.689e-05 | 3.349e-04 | 3.413e-06 | -1.172e-11 | 4.441e-16 |
| 5 | 2.47829e+05 | 1.014393 / 1.014366 | -2.640e-05 | 3.329e-04 | 1.788e-06 | -1.665e-11 | -1.443e-15 |
| 10 | 2.48313e+05 | 1.015463 / 1.015334 | -1.272e-04 | 2.510e-04 | 2.407e-06 | 1.453e-11 | -7.772e-16 |
| 20 | 2.49281e+05 | 1.012759 / 1.013177 | 4.125e-04 | 3.324e-04 | 1.669e-07 | -4.441e-11 | 2.220e-16 |
| 50 | 2.52187e+05 | 1.013499 / 1.013491 | -8.743e-06 | 4.386e-04 | -5.546e-06 | 6.994e-11 | 1.776e-15 |
| 99 | 2.56932e+05 | 1.015120 / 1.015109 | -1.054e-05 | 4.290e-04 | -1.490e-05 | 7.645e-11 | 2.887e-15 |
- max |relative difference| at equal model time: cycles 0-10 and the rest apart:
  - L_top: cycles 0-10 1.347e-04, cycles 11-end 4.125e-04
  - L_bot: cycles 0-10 3.349e-04, cycles 11-end 4.411e-04
  - L_in: cycles 0-10 2.220e-16, cycles 11-end 2.220e-16
  - KE_int: cycles 0-10 4.280e-06, cycles 11-end 1.591e-05
  - Etot: cycles 0-10 3.315e-11, cycles 11-end 1.472e-10
  - M_tot: cycles 0-10 4.663e-15, cycles 11-end 4.219e-15
  - Lz_tot: cycles 0-10 3.726e-07, cycles 11-end 7.426e-07
- MLQ1r: part 1 50 rows vs MLQ1 rows 1-50: BITWISE; part 2 50 rows vs MLQ1 rows 50-99: BITWISE (printed hst digits)
- radial E and T profiles vs MLQ0 at equal model time: **not available** (no bins are kept and the hst has no profile); the request says "hst / existing bins", and there are none.
- MLQ0: 2026-11-24 17:51:35 rc=0 (srun) ended: Terminating on cycle limit wall=320 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- MLQ1: 2026-11-24 17:56:35 rc=0 (srun) ended: Terminating on cycle limit wall=300 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- MLQ1ra: 2026-11-24 18:01:24 rc=0 (srun) ended: Terminating on cycle limit wall=289 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- MLQ1rb: 2026-11-24 18:04:35 rc=0 (srun) ended: Terminating on cycle limit wall=191 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- merge NOTE lines of MLQ1's log:
```
### NOTE: parameter 'merge_lat_rc' in block 'hydro' is absent from the restart's input: added from the command line (= 2.022e12)
### NOTE: parameter 'merge_lat_f' in block 'hydro' is absent from the restart's input: added from the command line (= 2)
### NOTE: parameter 'merge_lat_blockend' in block 'hydro' is absent from the restart's input: added from the command line (= exchange)
### NOTE: parameter 'merge_lat_diag' in block 'hydro' is absent from the restart's input: added from the command line (= 10)
MLAT merged lateral cells: 1 shell(s): x1v < 2.022e+12 f 2; mom rigid, rad fine
MLAT init: 884736 non-uniform gas groups, 0 non-uniform radiation groups -> averaged, ghosts refilled
MLATDIAG h cyc 2560 t 2.4825965697e+05 dt 1.525118e+02 M 2.3310565725424865e+32 E 1.3352846316879129e+47 P -7.4114538118e+32 2.3784055268e+34 3.0528546728e+33 L -1.3734661930871166e+45 4.2906567203172836e+44 -3.025401274
MLATDIAG h cyc 2560 t 2.4825965697e+05 dt 1.525118e+02 M 2.3310565725424865e+32 E 1.3352846407255275e+47 P -7.4677898501e+32 -2.8161327717e+34 3.0619477718e+33 L -1.3784697268373582e+45 4.2906309107513319e+44 -3.03538494
```

## PRF_zam (11-24 18:06)
TAG=PRF_zam BIN md5=30d01765970051212176cc2eb81ea034 RST=/resnick/groups/carnegie_poc/jingze/prod_rc4v_1124/m35zams/rst_keep/m35zams.stop1124b.rst N0=22155 NC=12 HOLDER=7333072 NODE=hpc-sm-02-03 NR=4 NSYS=/central/software9/spack/opt/spack/linux-rhel9-x86_64/n
run rc=0 (srun) ended: Terminating on cycle limit wall=47 s CW 0 FATAL 0
main loop time (before the final outputs) = 7.120714e+00
elapsed=9.658395e+00 cycle=22160
nsys_r0.sqlite: cycles 5 wall 640.4 ms/cyc  busy 477.5  kernel sum 477.5  launches 8139/cyc  keys 128
  ms/cyc     %   n/cyc    us/l  reg  blk     grid lmem  label | key
   42.90   9.0     408   105.1  130  128     2016    0   | radvet::RadiationVet::VetGdSweep#1
nsys_r0.sqlite: cycles 5 wall 640.4 ms/cyc  GPU idle 113.1 ms/cyc (17.7 %)  gaps 11311/cyc
idle by host activity (ms/cyc): MPI(part) 29.5  MPI 27.8  sync 18.4  memcpy-api(part) 16.9  sync(part) 8.7  memcpy-api 4.4  launch(part) 2.7  cuda-other 2.0  launch 1.7  cuda-other(part) 0.7  host-compute 0.4
idle by gap size   (ms/cyc, n/cyc): <10us 21.3 (8199)  10-100us 51.1 (3079)  0.1-1ms 17.2 (27)  >1ms 23.6 (5)
memcpy (kind: n/cyc, MB/cyc, ms/cyc; 1=HtoD 2=DtoH 8=DtoD 10=PtoP): 1: 1871, 7.3, 2.04  2: 364, 5.0, 0.86  8: 42, 1591.9, 0.89  10: 2774, 10083.6, 94.11
host time in (ms/cyc, overlapping GPU work or not): sync 435.3  MPI 336.1
largest gaps (us, host activity, previous kernel -> next):
     8311.0 sync(part)     M -> M
     7840.1 MPI(part)      M -> M
     7044.2 memcpy-api(part) M -> M
     7019.6 memcpy-api(part) M -> M
     7006.4 memcpy-api(part) M -> M
     6974.3 sync(part)     M -> M
- Text files: caltech_results/PRF_1124/PRF_zam/ (res.txt, kname_r0..3.txt, gaps_r0..3.txt, stats_r0.txt, CLEAN PASS/cycle lines, hst). Plain s/cyc of the same state in the CAND arm: 0.297 s/cyc (m35zams_R). Only the kernel SHARES of this run count (nsys slows the step). nsys from cuda-12.9.0.
