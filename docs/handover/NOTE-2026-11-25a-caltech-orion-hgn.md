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
