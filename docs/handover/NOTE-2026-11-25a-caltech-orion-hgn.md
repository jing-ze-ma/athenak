# NOTE 2026-11-25a caltech -> orion: HG profile and old-vs-new time on H200 (your 25a)

- **Layout:** 4 H200 do not fit: rel-f on 02-03 with 4 ranks (6 blocks per GPU) stopped in cycle 3 with `Kokkos ERROR: Cuda memory space failed to allocate 4.254 GiB (label=cln_rk)` at 140.7-141.1 GB per GPU (posted from chk4). All arms below: 8 ranks on 02-03 + 02-09 (3 blocks per GPU). State: hegiant.00205.rst (the newest HG rst at Caltech), N0 794784, t 6.6096e7, dt 90.0. Keys: HR's args (the newest HG keys run here), thinf true, no outputs, ndiag 1, hydro/dt_report = 10 (DTREP).
- Old binary: athena_gpu_he_star_m1_prod1118_1d197280 (1d1972805, md5 26a473d9, CUDA 12.9, built 11-13) exists; HGO = it at the same rst with the rst's own stored keys + table paths only (the rst was written by the old scheme; see the run line). It has no final_restart key, so it writes its final rst into its own dir (kept, nothing is deleted).

## HGN (11-24 23:58): rel-f 80d16d39, 60 cycles
- 2026-11-24 23:58:54 rc=0 (srun) ended: Terminating on cycle limit wall=250 s CW 0 FATAL 0
### HGN: 61 elapsed prints, cycles 794784-794844
- s/cycle last 50: mean 3.4793, median 3.3874, min 3.2791, max 6.2698; whole loop 2.098805e+02 s
- cycles with wall > 1.25 x median (cycle, s): [(794785, 6.411), (794844, 6.27)]; spacing: [59]
- dt per cycle (min/mean/max): 90.0000/90.0000/90.0000; first 6 and last 3: ['90.0000', '90.0000', '90.0000', '90.0000', '90.0000', '90.0000', '90.0000', '90.0000', '90.0000']
- stage t2s 2: its (sum) mean 16.00 max 22; solves mean 4.567; lines 60
- stage t2s 3: its (sum) mean 14.95 max 20; solves mean 5.633; lines 60
- DTREP lines 12; last: DTREP cycle=794840 time=6.6101e+07 cfl*dt_hyd=92.1702 cfl*dt_gasonly=124.533 rank=5 gid=15 kji=66,18,23 dir=1 r=7.34434e+11 rho=5.81522e-06 T=5.31723e+13 pg=2.30817e+08 G1=1.66667 cgas=8.13346e+06 v=4.13655e+06 ceff=1.27621e+07 h=5.19184e+09 tauT_h=10484 Erad=1.26549e+09 Ddd=0.333333 cdt_kP=4.38354e
- lines starting "CLEAN FWB": 0; first: -
- lines starting "CLEAN FCONE": 0; first: -
- lines starting "CLEAN GNSW": 5; first: CLEAN GNSW cycle 794793 t2s 2 pass 5: stalled, the gas Newton slope off (1 so far) | last: CLEAN GNSW cycle 794822 t2s 2 pass 5: stalled, the gas Newton slope off (5 so far)
- lines starting "CLEAN XT": 0; first: -
- lines starting "CLEAN CLOS": 0; first: -
- lines starting "<rad_vet> vet_gd": 18; first: <rad_vet> vet_gd_rotate_xfade = true: the closure is cross-faded from the current to the next ray set over each rotation block (two sweeps per build) | last: <rad_vet> vet_gd_rotate_xfade_fuse (rank 0): fused builds 6.000000e+01, two-sweep fallbacks 0.000000e+00, merged exchanges 0.000000e+00, out-of-slab accesses (other set) 0.000000e+00
- lines starting "<rad_vet> vet_rotate": 3; first: <rad_vet> vet_rotate_advance = true: with vet_rotate_avg the rotated sets do not repeat; the averaged last-K window advances (no effect without vet_rotate_avg) | last: <rad_vet> vet_rotate_avg_ema = 128: the closure is an exponential moving average over the rotated sets (a = 1/N per block; change test eps 0.05); vet_rotate_avg is not used
- lines starting "<rad_vet> vet_gd_rotate": 3; first: <rad_vet> vet_gd_rotate_xfade = true: the closure is cross-faded from the current to the next ray set over each rotation block (two sweeps per build) | last: <rad_vet> vet_gd_rotate_xfade_fuse (rank 0): fused builds 6.000000e+01, two-sweep fallbacks 0.000000e+00, merged exchanges 0.000000e+00, out-of-slab accesses (other set) 0.000000e+00
- CW 0 FATAL 0 inf 0 DT_RETRY 0; fallback 0.000000e+00; BE 0.000000e+00
- end lines of the rays / split timers: ['<rad_vet> vet_gd cs band (rank 0): compact, intensity values sent 1.036078e+09 of dense 3.119790e+11 (3.320985e-03); list builds 7.000000e+00 (3.428292e+00 s)', '<rad_vet> vet_gd ragged band (rank 0): out-of-slab accesses 0.000000e+00', '<rad_vet> vet_gd_band_ring (rank 0): sweeps 6.300000e+01, dense band saves 0.000000e+00, save store 1.297930e-01 GB', '<rad_vet> vet_gd_rotate_xfade_fuse (rank 0): fused builds 6.000000e+01, two-sweep fallbacks 0.000000e+00, merged exchanges 0.000000e+00, out-of-slab accesses (other set) 0.000000e+00', '<rad_vet> (rank 0) vet_gd halo lists 0.000000e+00 made, 0.000000e+00 MB of 0.000000e+00 (key 4096, GPU: <= free/2); exchanges by list 0.000000e+00, dense 0.000000e+00', '<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 2.087758e+02 | failed stage attempts 0.000000e+00 | their BE redo 0.000000e+00 | auto_be steps 0.000000e+00 | other BE steps 0.000000e+00 | discarded dt_retry a']
- peak memory per GPU (MiB, nvidia-smi samples): {'\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00hpc-sm-02-09:0': 5, 'hpc-sm-02-03:0': 75591, 'hpc-sm-02-03:1': 75931, 'hpc-sm-02-03:2': 75917, 'hpc-sm-02-03:3': 75651, 'hpc-sm-02-09:0': 75595, 'hpc-sm-02-09:1': 76011, 'hpc-sm-02-09:2': 75997, 'hpc-sm-02-09:3': 75587}

## HGO (11-25 00:05): old HG production binary 1d1972805, its own (rst-stored) keys, 60 cycles
- 2026-11-25 00:05:07 rc=0 (srun) ended: Terminating on cycle limit wall=361 s CW 0 FATAL 0
### HGO: 61 elapsed prints, cycles 794784-794844
- s/cycle last 50: mean 4.9378, median 3.2600, min 1.9292, max 46.4885; whole loop 2.521928e+02 s
- cycles with wall > 1.25 x median (cycle, s): [(794785, 7.858), (794786, 5.821), (794787, 10.198), (794788, 4.583), (794794, 7.122), (794795, 5.913), (794797, 6.998), (794798, 4.283), (794806, 8.15), (794807, 8.796), (794809, 5.0), (794810, 8.548), (794811, 4.478), (794813, 8.276), (794814, 4.501), (794815, 7.996), (794816, 4.572), (794817, 4.532), (794818, 8.133), (794819, 4.705)]; spacing: [1, 1, 1, 6, 1, 2, 1, 8, 1, 2, 1, 1]
- dt per cycle (min/mean/max): 22.5000/67.8689/90.0000; first 6 and last 3: ['90.0000', '90.0000', '90.0000', '90.0000', '90.0000', '90.0000', '22.5000', '22.5000', '22.5000']
- DTREP lines 0; last: -
- lines starting "CLEAN FWB": 0; first: -
- lines starting "CLEAN FCONE": 0; first: -
- lines starting "CLEAN GNSW": 0; first: -
- lines starting "CLEAN XT": 0; first: -
- lines starting "CLEAN CLOS": 0; first: -
- lines starting "<rad_vet> vet_gd": 0; first: -
- lines starting "<rad_vet> vet_rotate": 0; first: -
- lines starting "<rad_vet> vet_gd_rotate": 0; first: -
- CW 0 FATAL 0 inf 0 DT_RETRY 2; fallback 8.000000e+00; BE 2.000000e+01
- end lines of the rays / split timers: []
- peak memory per GPU (MiB, nvidia-smi samples): {'hpc-sm-02-03:0': 66979, 'hpc-sm-02-03:1': 67271, 'hpc-sm-02-03:2': 67257, 'hpc-sm-02-03:3': 66977, 'hpc-sm-02-09:0': 66653, 'hpc-sm-02-09:1': 67103, 'hpc-sm-02-09:2': 67093, 'hpc-sm-02-09:3': 66655, 'hpc-sm-02-c-sm-02-03:1': 66843}

## HGN and HGO: the per-cycle statistics above include the final print (which contains the final bins / rst) and the start-up cycle; corrected:
- HGN: s/cycle without the first and the last print: mean 3.4484 median 3.3946 min 3.2791 max 4.1656 (58 cycles); last 50 before the final print: mean 3.4372 median 3.3874; model time advanced 5.4000e+03 s in 60 cycles (mean dt 90.00); wall of the main loop 2.098805e+02 s = 38.87 s per 1000 s model
  dt sequence: 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90
  cycles above 1.15 x median: [(794794, 4.17), (794807, 3.93), (794820, 3.92), (794833, 3.95)]
- HGO: s/cycle without the first and the last print: mean 4.1682 median 3.0284 min 1.9292 max 10.1981 (58 cycles); last 50 before the final print: mean 4.1505 median 3.2600; model time advanced 4.0500e+03 s in 60 cycles (mean dt 67.50); wall of the main loop 2.521928e+02 s = 62.27 s per 1000 s model
  dt sequence: 90 90 90 90 90 90 90 90 90 90 90 90 90 90 22.5 22.5 22.5 22.5 22.5 22.5 22.5
  cycles above 1.15 x median: [(794786, 5.82), (794787, 10.2), (794788, 4.58), (794794, 7.12), (794795, 5.91), (794797, 7.0), (794798, 4.28), (794806, 8.15), (794807, 8.8), (794809, 5.0), (794810, 8.55), (794811, 4.48), (794812, 3.52), (794813, 8.28), (794814, 4.5), (794815, 8.0), (794816, 4.57), (794817, 4.53), (794818, 8.13), (794819, 4.7), (794820, 6.44), (794821, 4.66), (794822, 4.49), (794823, 3.55), (794824, 9.88), (794825, 4.74), (794828, 3.78)]
- HGO retry/fallback lines: ['DT_RETRY cycle=794824 t=6.609960e+07 dt 9.000000e+01 -> 4.500000e+01 reason=eint<=0 resid=5.000000e+00 cut 1\n', 'DT_RETRY cycle=794824 t=6.609960e+07 dt 4.500000e+01 -> 2.250000e+01 reason=eint<=0 resid=2.500000e+01 cut 2\n']
- HGO's last print includes writing its final rst (46 s); the old binary has no final_restart key. HGO is NOT at the same dt as HGN: it runs 90 s for the first cycles, then drops to 22.5 (DT_RETRY 2, stage fallback 8, BE 20) with the old scheme's own rule, so its cycles are partly retried steps and its s/cycle is not at the same state as HGN's (compare per model time, and read dt from the sequence).

## PRF_hgn: hung, cancelled (no profile)
- The nsys run (rel-f, 8 ranks on 02-03 + 02-09, HR's args + clean_verbose 2, nsys from CUDA 12.9, -t cuda,nvtx,mpi) started 00:05:31 and stopped advancing: its run.log's last line is the end of cycle 794785 (written 00:06:46; the elapsed print of 794785 at 26.3 s), nothing after it for 26 min; memory log still running. Step 7333072.544 cancelled at 00:33 by its id. No nsys report file was written, so there is no partial profile. HGN without nsys passed the same cycle (250 s for 60 cycles on the same 8 ranks).
- Not retried: the cause is not known (a hang under nsys with the MPI trace on two nodes is the first suspect, not proved). If you want it, options: nsys without the mpi trace, or one node if a smaller cut fits (4 H200 do not: OOM), or 12 cycles without the 13-cycle rotation block. Tell me which.
- HGN, HGO and the rotation cadence are posted above; the kernel table (res.txt, kname, gaps) does not exist for HG.

## HGV (25d item 1): HGN's command + clean_verbose = 2, 30 cycles, no nsys, 8 ranks
- 2026-11-25 01:07:17 rc=0 (srun) ended: Terminating on cycle limit wall=150 s CW 0 FATAL 0
### HGV: 31 elapsed prints, cycles 794784-794814
- s/cycle last 30: mean 3.5504, median 3.2893, min 3.1669, max 6.3270; whole loop 1.035006e+02 s
- cycles with wall > 1.25 x median (cycle, s): [(794785, 6.327), (794814, 6.204)]; spacing: [29]
- dt per cycle (min/mean/max): 90.0000/90.0000/90.0000; first 6 and last 3: ['90.0000', '90.0000', '90.0000', '90.0000', '90.0000', '90.0000', '90.0000', '90.0000', '90.0000']
- stage t2s 2: its (sum) mean 16.53 max 22; solves mean 4.900; lines 30
- stage t2s 3: its (sum) mean 15.07 max 20; solves mean 5.633; lines 30
- CLEAN PASS lines per cycle: mean 10.53 max 14 (30 cycles)
- DTREP lines 6; last: DTREP cycle=794810 time=6.60983e+07 cfl*dt_hyd=90.2453 cfl*dt_gasonly=123.272 rank=5 gid=15 kji=66,19,23 dir=1 r=7.34434e+11 rho=5.62572e-06 T=5.33763e+13 pg=2.24152e+08 G1=1.66667 cgas=8.14904e+06 v=4.30047e+06 ceff=1.29586e+07 h=5.19184e+09 tauT_h=9979.73 Erad=1.28501e+09 Ddd=0.333333 cdt_kP=4.042
- lines starting "CLEAN FWB": 60; first: CLEAN FWB cycle 794784 t2s 2 faces 151577909 radiation momentum sum |dp| 5.092018e+29 x1 faces sum dp 8.701600e+27 | last: CLEAN FWB cycle 794813 t2s 3 faces 151547295 radiation momentum sum |dp| 6.242986e+29 x1 faces sum dp 1.516717e+28
- lines starting "CLEAN FCONE": 60; first: CLEAN FCONE cycle 794784 t2s 2 faces 36237663 momentum cut sum |dF| V/c^2 1.557552e+29 | last: CLEAN FCONE cycle 794813 t2s 3 faces 31459105 momentum cut sum |dF| V/c^2 3.815734e+28
- lines starting "CLEAN GNSW": 3; first: CLEAN GNSW cycle 794793 t2s 2 pass 5: stalled, the gas Newton slope off (1 so far) | last: CLEAN GNSW cycle 794807 t2s 2 pass 5: stalled, the gas Newton slope off (3 so far)
- lines starting "CLEAN XT": 60; first: CLEAN XT cycle 794784 t2s 2 write-back closure 2.479e-03 | last: CLEAN XT cycle 794813 t2s 3 write-back closure 4.764e-03
- lines starting "CLEAN CLOS": 0; first: -
- lines starting "<rad_vet> vet_gd": 18; first: <rad_vet> vet_gd_rotate_xfade = true: the closure is cross-faded from the current to the next ray set over each rotation block (two sweeps per build) | last: <rad_vet> vet_gd_rotate_xfade_fuse (rank 0): fused builds 3.000000e+01, two-sweep fallbacks 0.000000e+00, merged exchanges 0.000000e+00, out-of-slab accesses (other set) 0.000000e+00
- lines starting "<rad_vet> vet_rotate": 3; first: <rad_vet> vet_rotate_advance = true: with vet_rotate_avg the rotated sets do not repeat; the averaged last-K window advances (no effect without vet_rotate_avg) | last: <rad_vet> vet_rotate_avg_ema = 128: the closure is an exponential moving average over the rotated sets (a = 1/N per block; change test eps 0.05); vet_rotate_avg is not used
- lines starting "<rad_vet> vet_gd_rotate": 3; first: <rad_vet> vet_gd_rotate_xfade = true: the closure is cross-faded from the current to the next ray set over each rotation block (two sweeps per build) | last: <rad_vet> vet_gd_rotate_xfade_fuse (rank 0): fused builds 3.000000e+01, two-sweep fallbacks 0.000000e+00, merged exchanges 0.000000e+00, out-of-slab accesses (other set) 0.000000e+00
- CW 0 FATAL 0 inf 0 DT_RETRY 0; fallback 0.000000e+00; BE 0.000000e+00
- end lines of the rays / split timers: ['<rad_vet> vet_gd cs band (rank 0): compact, intensity values sent 5.203564e+08 of dense 1.589212e+11 (3.274305e-03); list builds 5.000000e+00 (2.515761e+00 s)', '<rad_vet> vet_gd ragged band (rank 0): out-of-slab accesses 0.000000e+00', '<rad_vet> vet_gd_band_ring (rank 0): sweeps 3.300000e+01, dense band saves 0.000000e+00, save store 1.206439e-01 GB', '<rad_vet> vet_gd_rota
Traceback (most recent call last):
  File "/resnick/groups/carnegie_poc/jingze/b1_1124/hgn.py", line 30, in <module>
    for l in open(R + a + '/mem.log', errors='replace'):
FileNotFoundError: [Errno 2] No such file or directory: '/resnick/groups/carnegie_poc/jingze/b1_1124/hgn/HGV/mem.log'
- CLEAN PASS lines: 316; first: CLEAN PASS cycle 794784 t2s 2 pass 0 it 9 de 1.343e-02 dv 5.851e-04 tol 1.000e-05 fe 1.000e-10 fv 3.118e-10
- passes per cycle and stage (number of CLEAN PASS lines): c794784 {2: 7, 3: 6}; c794785 {2: 5, 3: 6}; c794786 {2: 5, 3: 9}; c794787 {2: 6, 3: 6}; c794788 {2: 5, 3: 9}; c794789 {2: 6, 3: 5}; c794790 {2: 5, 3: 5}; c794791 {2: 5, 3: 5}; c794792 {2: 6, 3: 6}; c794793 {2: 7, 3: 6}; c794794 {2: 7, 3: 5}; c794795 {2: 4, 3: 5}; c794796 {2: 4, 3: 5}; c794797 {2: 4, 3: 6}; c794798 {2: 4, 3: 5}; c794799 {2: 7, 3: 5}; c794800 {
- user.hst vs HGN: 30 rows in common, printed-bitwise
- All CLEAN PASS lines and per-cycle prints: caltech_results/HGN/HGV/lines.txt.gz.
- The hung nsys run, one line: no utilisation was logged; its memory samples (nvidia-smi every 2 s) stay at 73.6-74.1 GB per GPU on all 8 for the whole 26 min (no growth, no release), the log ended inside cycle 794786, no stack was taken.
