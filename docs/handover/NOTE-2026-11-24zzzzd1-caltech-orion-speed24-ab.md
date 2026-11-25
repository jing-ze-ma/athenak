# NOTE 2026-11-24zzzzd1 caltech -> orion: H200 A/B of speed24 e7b2eea39 (your 24zzzzd1)

## Extra, no run: MLATDIAG nfb in the bsg logs (read-only grep, 11-24 18:19)
- nfb is the running total since the log's start (hydro_mergelat.cpp: nfallback += per call); a group falls back when its rotational energy would exceed half the remaining thermal energy, and its omega is then set to 0.
- c1 bsg (prod_c1_1124/bsg/log.0; printed every 1000 cycles, two lines per print):
```
cycle t nfb(line 1) nfb(line 2)
0 0.0000000000e+00 0 0
0 0.0000000000e+00 0 0
1000 1.6884532047e+05 0 0
2000 3.3812388789e+05 0 0
3000 5.0752915826e+05 0 0
4000 6.7697849784e+05 0 0
5000 8.4646420869e+05 0 0
6000 1.0121137691e+06 0 0
7000 1.1622801965e+06 101260 101260
8000 1.3277475984e+06 105044 105044
9000 1.4969114017e+06 575785 575791
```
- Other bsg production logs (max nfb in the log, last printed cycle):
```
prod_c1_1124/bsg/log.0 lines 22 max nfb 575791 last cyc 9000
prod_fresh_1121b/bsg/log.0 lines 38 max nfb 556069 last cyc 35000
prod_fresh_1121/bsg/log.0 lines 13 max nfb 2136263 last cyc 10000
prod_ia11_1123/bsg/log.0 lines 68 max nfb 357879 last cyc 32000
prod_ia9_1123/bsg/log.0 lines 48 max nfb 337496 last cyc 22000
prod_rc3_1124/bsg/log.0 lines 6 max nfb 1951991 last cyc 1000
prod_rc3_1124/bsg/log.1 lines 8 max nfb 58197905 last cyc 4000
prod_rc3_1124/bsg/log.2 lines 68 max nfb 2888321104 last cyc 38000
prod_rc3_1124/bsg/log.3 lines 2 max nfb 51505833 last cyc 39000
prod_rc3_1124/bsg/log.4 lines 32 max nfb 397460 last cyc 55000
prod_rc4v_1124/bsg/log.0 lines 6 max nfb 0 last cyc 1000
prod_rc4v_1124/bsg/log.1 lines 10 max nfb 0 last cyc 6000
```
- No reading attached: what it means for L_z is yours. The watchers do not count nfb.

## Build (own worktree wt_orion1121/speed24 at e7b2eea39, CUDA 12.9 + g++ 13, HOPPER90, CPU step on 02-09)
- BUILD speed24 rc=0 wall=402 s TUs=250
- BUILD_OK speed24 md5=601f709cb98bf3632f981d79bb9fe5ac prov=e7b2eea39f39130a190fa1f2e0dc18e0e3d8f36a
- GATE PASS speed24: athena_speed24_e7b2eea39 md5 601f709c node hpc-sm-02-09 pass1=3.686e-11 log date 11-24 18:23:41 end: Terminating on cycle limit
- Diff 95813cf85..e7b2eea39 (src, tools, docs):  docs/vet_clean/KEYS.md, src/hydro/hydro.cpp, src/hydro/hydro_newdt.cpp, src/rad_vet/rad_vet.cpp, src/rad_vet/rad_vet.hpp, src/rad_vet/rad_vet_clean.cpp, src/rad_vet/rad_vet_clean.hpp, src/rad_vet/rad_vet_clean_solver.cpp, src/rad_vet/rad_vet_dtretry_visit.hpp, src/rad_vet/rad_vet_fluid.cpp, tools/dtretry_gen.py, 11 files changed, 285 insertions(+), 55 deletions(-),
- Plan: 17 arms (the note's base / tip / base / tip taken literally for F: 12 runs; then G 2, H 3; G off = F tip run 1), plus 2 nsys. The final bins of each arm stay on disk.

## F: A1+A2 base / tip / base / tip (11-24 19:26)
- F_m35_b1: 2026-11-24 18:54:45 rc=0 (srun) ended: Terminating on cycle limit wall=137 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- F_m35_b2: 2026-11-24 18:59:15 rc=0 (srun) ended: Terminating on cycle limit wall=136 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- F_m35_t1: 2026-11-24 18:56:58 rc=0 (srun) ended: Terminating on cycle limit wall=131 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- F_m35_t2: 2026-11-24 19:01:31 rc=0 (srun) ended: Terminating on cycle limit wall=135 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- F_t9q_b1: 2026-11-24 19:10:19 rc=0 (srun) ended: Terminating on cycle limit wall=320 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- F_t9q_b2: 2026-11-24 19:21:06 rc=0 (srun) ended: Terminating on cycle limit wall=333 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- F_t9q_t1: 2026-11-24 19:15:32 rc=0 (srun) ended: Terminating on cycle limit wall=297 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- F_t9q_t2: 2026-11-24 19:26:20 rc=0 (srun) ended: Terminating on cycle limit wall=313 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- F_zam_b1: 2026-11-24 19:02:21 rc=0 (srun) ended: Terminating on cycle limit wall=50 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- F_zam_b2: 2026-11-24 19:04:09 rc=0 (srun) ended: Terminating on cycle limit wall=56 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- F_zam_t1: 2026-11-24 19:03:10 rc=0 (srun) ended: Terminating on cycle limit wall=48 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- F_zam_t2: 2026-11-24 19:04:58 rc=0 (srun) ended: Terminating on cycle limit wall=48 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 

| arm | rows | s/cyc (last two elapsed prints) | loop s | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | CW/FATAL/inf/DT_RETRY | floor_pay sum | L_top/L_in end | mem max GB |
|---|---|---|---|---|---|---|---|---|---|---|
| F_m35_b1 | 100 | 1.0957 (c7000-7050) | 111.1 | 12.01/12.01/12.01 | 36.97, 5.180 | 18.24, 4.020 | 0/0/0/0 | 8.74e+25 | 0.988138 | see mem.log |
| F_m35_b2 | 100 | 1.0974 (c7000-7050) | 111.1 | 12.01/12.01/12.01 | 36.97, 5.180 | 18.24, 4.020 | 0/0/0/0 | 8.74e+25 | 0.988138 | see mem.log |
| F_m35_t1 | 100 | 1.0381 (c7000-7050) | 105.2 | 12.01/12.01/12.01 | 36.97, 5.180 | 18.24, 4.020 | 0/0/0/0 | 8.74e+25 | 0.988138 | see mem.log |
| F_m35_t2 | 100 | 1.0383 (c7000-7050) | 105.2 | 12.01/12.01/12.01 | 36.97, 5.180 | 18.24, 4.020 | 0/0/0/0 | 8.74e+25 | 0.988138 | see mem.log |
| F_t9q_b1 | 100 | 2.5821 (c2600-2650) | 255.6 | 96.84/96.84/96.84 | 5.43, 2.900 | 9.15, 2.510 | 0/0/0/0 | 4.65e+27 | 1.015120 | see mem.log |
| F_t9q_b2 | 100 | 2.5874 (c2600-2650) | 256.2 | 96.84/96.84/96.84 | 5.43, 2.900 | 9.15, 2.510 | 0/0/0/0 | 4.65e+27 | 1.015120 | see mem.log |
| F_t9q_t1 | 100 | 2.2406 (c2600-2650) | 221.5 | 96.84/96.84/96.84 | 5.43, 2.900 | 9.15, 2.510 | 0/0/0/0 | 4.65e+27 | 1.015120 | see mem.log |
| F_t9q_t2 | 100 | 2.2409 (c2600-2650) | 221.7 | 96.84/96.84/96.84 | 5.43, 2.900 | 9.15, 2.510 | 0/0/0/0 | 4.65e+27 | 1.015120 | see mem.log |
| F_zam_b1 | 100 | 0.2951 (c22240-22250) | 33.0 | 3.481/3.481/3.481 | 18.96, 2.000 | 16.63, 2.240 | 0/0/0/0 | 0.00e+00 | 1.000146 | see mem.log |
| F_zam_b2 | 100 | 0.2953 (c22240-22250) | 32.9 | 3.481/3.481/3.481 | 18.96, 2.000 | 16.63, 2.240 | 0/0/0/0 | 0.00e+00 | 1.000146 | see mem.log |
| F_zam_t1 | 100 | 0.2856 (c22240-22250) | 32.2 | 3.481/3.481/3.481 | 18.96, 2.000 | 16.63, 2.240 | 0/0/0/0 | 0.00e+00 | 1.000146 | see mem.log |
| F_zam_t2 | 100 | 0.2876 (c22240-22250) | 32.3 | 3.481/3.481/3.481 | 18.96, 2.000 | 16.63, 2.240 | 0/0/0/0 | 0.00e+00 | 1.000146 | see mem.log |

History identity (printed digits; user and hydro):
- F_m35_b1 vs F_m35_b2: user 100 rows IDENTICAL; hydro 100 rows IDENTICAL
- F_m35_b1 vs F_m35_t1: user 100 rows IDENTICAL; hydro 100 rows IDENTICAL
- F_m35_b2 vs F_m35_t2: user 100 rows IDENTICAL; hydro 100 rows IDENTICAL
- F_zam_b1 vs F_zam_b2: user 100 rows IDENTICAL; hydro 100 rows IDENTICAL
- F_zam_b1 vs F_zam_t1: user 100 rows IDENTICAL; hydro 100 rows IDENTICAL
- F_zam_b2 vs F_zam_t2: user 100 rows IDENTICAL; hydro 100 rows IDENTICAL
- F_t9q_b1 vs F_t9q_b2: user 100 rows IDENTICAL; hydro 100 rows IDENTICAL
- F_t9q_b1 vs F_t9q_t1: user 100 rows IDENTICAL; hydro 100 rows IDENTICAL
- F_t9q_b2 vs F_t9q_t2: user 100 rows IDENTICAL; hydro 100 rows IDENTICAL

## nsys of the tip (12 cycles, the same cut as the base profiles of zzzzc2); us per launch, rank 0
- t9q base: nsys_r0.sqlite: cycles 11 wall 2413.0 ms/cyc  busy 2244.5  kernel sum 2244.5  launches 5361/cyc  keys 130
```
  203.29   9.1      19 10961.9   80  128   132096    0   | radvet::RadiationVet::CleanStagFSolve#6
  103.22   4.6      19  5566.0   75  128   132096    0   | radvet::RadiationVet::CleanStagFCone#4
   63.93   2.8      19  3447.1   80  128   132096    0   | radvet::RadiationVet::CleanStagFSolve#3
   35.49   1.6       2 17747.4  126  128   131072    0   | radvet::RadiationVet::FloorGasCap#2
   13.94   0.6      25   563.8   64  128   132096    0   | radvet::RadiationVet::CleanStagFSolve#8
   13.75   0.6      37   370.6   80  128   131072    0   | radvet::RadiationVet::CleanStagFSolve#5
   12.71   0.6      19   685.2   78  128   131072    0   | radvet::RadiationVet::CleanStagFSolve#4
   11.60   0.5      19   625.4   80  128   132096    0   | radvet::RadiationVet::CleanStagFSolve#7
   11.26   0.5       6  1822.0   68  256    32768    0   | radvet::RadiationVet::CleanStagFSolve#1
   10.29   0.5       6  1664.8   80  128   131072    0   | radvet::RadiationVet::CleanStagFCone#1
    8.11   0.4      12   656.3   91   32     4096    0   | radvet::RadiationVet::CleanStagFSolve#2
    5.98   0.3      19   322.6   56  128   132096    0   | radvet::RadiationVet::CleanStagFCone#2
nsys_r0.sqlite: cycles 11 wall 2413.0 ms/cyc  GPU idle 141.5 ms/cyc (5.9 %)  gaps 6346/cyc
idle by gap size   (ms/cyc, n/cyc): <10us 13.6 (5008)  10-100us 22.9 (1304)  0.1-1ms 6.4 (25)  >1ms 98.5 (9)
1.4,750606620,2973,252474.5,6930.0,560,40004728,2248482.7,collectives,MPI_Allreduce
```
- t9q tip: nsys_r0.sqlite: cycles 11 wall 2017.7 ms/cyc  busy 1936.0  kernel sum 1936.0  launches 5373/cyc  keys 130
```
   63.85   3.3      19  3443.1   80  128   132096    0   | radvet::RadiationVet::CleanStagFSolve#3
   24.59   1.3      19  1326.1   80  128   132096    0   | radvet::RadiationVet::CleanStagFSolve#6
   13.93   0.7      25   563.4   64  128   132096    0   | radvet::RadiationVet::CleanStagFSolve#8
   13.74   0.7      37   370.4   80  128   131072    0   | radvet::RadiationVet::CleanStagFSolve#5
   12.70   0.7      19   685.0   78  128   131072    0   | radvet::RadiationVet::CleanStagFSolve#4
   11.59   0.6      19   624.8   80  128   132096    0   | radvet::RadiationVet::CleanStagFSolve#7
   11.25   0.6       6  1820.4   68  256    32768    0   | radvet::RadiationVet::CleanStagFSolve#1
    9.74   0.5       6  1574.9   80  128   131072    0   | radvet::RadiationVet::CleanStagFCone#1
    8.33   0.4      19   449.4   76  128   132096    0   | radvet::RadiationVet::CleanStagFCone#4
    8.11   0.4      12   656.0   91   32     4096    0   | radvet::RadiationVet::CleanStagFSolve#2
    5.98   0.3      19   322.3   56  128   132096    0   | radvet::RadiationVet::CleanStagFCone#2
nsys_r0.sqlite: cycles 11 wall 2017.7 ms/cyc  GPU idle 55.3 ms/cyc (2.7 %)  gaps 6364/cyc
idle by gap size   (ms/cyc, n/cyc): <10us 13.6 (5032)  10-100us 23.6 (1306)  0.1-1ms 6.1 (23)  >1ms 12.1 (2)
0.1,68016747,2899,23462.1,7060.0,890,15099133,290035.5,collectives,MPI_Allreduce
```
- m35 base: nsys_r0.sqlite: cycles 11 wall 1207.5 ms/cyc  busy 996.4  kernel sum 996.4  launches 6671/cyc  keys 133
```
   39.69   4.0      34  1173.5   80  128    43344    0   | radvet::RadiationVet::CleanStagFSolve#3
   37.20   3.7      34  1100.0   80  128    43344    0   | radvet::RadiationVet::CleanStagFSolve#6
   16.49   1.7       2  8243.5  126  128    43008    0   | radvet::RadiationVet::FloorGasCap#2
    8.39   0.8      68   124.0   80  128    43008    0   | radvet::RadiationVet::CleanStagFSolve#5
    7.74   0.8      11   686.5   68  256    16384    0   | radvet::RadiationVet::CleanStagFSolve#1
    7.71   0.8      34   227.9   78  128    43008    0   | radvet::RadiationVet::CleanStagFSolve#4
    7.17   0.7      45   159.0   64  128    43344    0   | radvet::RadiationVet::CleanStagFSolve#8
    6.92   0.7      23   307.0   91   32     1344    0   | radvet::RadiationVet::CleanStagFSolve#2
    6.66   0.7      34   196.8   75  128    43344    0   | radvet::RadiationVet::CleanStagFCone#4
    6.43   0.6      34   190.1   80  128    43344    0   | radvet::RadiationVet::CleanStagFSolve#7
    5.97   0.6      11   529.8   80  128    43008    0   | radvet::RadiationVet::CleanStagFCone#1
    3.71   0.4      34   109.7   56  128    43344    0   | radvet::RadiationVet::CleanStagFCone#2
nsys_r0.sqlite: cycles 11 wall 1207.5 ms/cyc  GPU idle 187.6 ms/cyc (15.5 %)  gaps 9418/cyc
idle by gap size   (ms/cyc, n/cyc): <10us 18.6 (6997)  10-100us 44.7 (2386)  0.1-1ms 5.3 (20)  >1ms 119.0 (15)
5.1,1225115264,2958,414170.1,3190.0,710,32392138,2914376.4,collectives,MPI_Allreduce
```
- m35 tip: nsys_r0.sqlite: cycles 11 wall 1060.4 ms/cyc  busy 955.6  kernel sum 955.6  launches 6688/cyc  keys 133
```
   39.65   4.1      34  1172.5   80  128    43344    0   | radvet::RadiationVet::CleanStagFSolve#3
   12.70   1.3      34   375.6   80  128    43344    0   | radvet::RadiationVet::CleanStagFSolve#6
    8.39   0.9      68   124.0   80  128    43008    0   | radvet::RadiationVet::CleanStagFSolve#5
    7.74   0.8      11   686.2   68  256    16384    0   | radvet::RadiationVet::CleanStagFSolve#1
    7.71   0.8      34   227.9   78  128    43008    0   | radvet::RadiationVet::CleanStagFSolve#4
    7.18   0.8      45   159.1   64  128    43344    0   | radvet::RadiationVet::CleanStagFSolve#8
    6.92   0.7      23   307.0   91   32     1344    0   | radvet::RadiationVet::CleanStagFSolve#2
    6.42   0.7      34   189.8   80  128    43344    0   | radvet::RadiationVet::CleanStagFSolve#7
    5.96   0.6      11   529.0   80  128    43008    0   | radvet::RadiationVet::CleanStagFCone#1
    3.94   0.4      34   116.6   76  128    43344    0   | radvet::RadiationVet::CleanStagFCone#4
    3.71   0.4      34   109.6   56  128    43344    0   | radvet::RadiationVet::CleanStagFCone#2
nsys_r0.sqlite: cycles 11 wall 1060.4 ms/cyc  GPU idle 81.4 ms/cyc (7.7 %)  gaps 9436/cyc
idle by gap size   (ms/cyc, n/cyc): <10us 18.5 (7007)  10-100us 44.1 (2409)  0.1-1ms 4.8 (16)  >1ms 14.0 (3)
0.3,55191123,2823,19550.5,3060.0,730,29372780,553901.1,collectives,MPI_Allreduce
```
