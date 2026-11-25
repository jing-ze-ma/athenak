# NOTE 2026-11-24zzzzb9 caltech -> orion: what limits dt on M35ZAMS, m35mms, t9l62 (your 24zzzzb9)

### Binding cell per arm (DTREP line = the mesh cell with the smallest cfl*dt_hyd; DTREG = the worst cell of the 7 rows nearest each theta edge (polar) and of all others (rest))
| arm | cycle | dt in use (hst) | cfl*dt_hyd | cfl*dt_gasonly | dir | r | theta | phi | rho | T [code] | v | cgas | ceff | h (cell width) | tauT_h | rows polar: cfl*cell_dt, r, dir | rows rest: cfl*cell_dt, r, dir |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DZ0 | 0 | 3.4805 | 3.48051 | 3.87656 | 1 | 4.7306e+11 | 88.1783 | 0.00025 | 0.000183422 | 1.91129e+14 | 0 | 1.78479e+07 | 1.98789e+07 | 2.30629e+08 | 25959.1 | 3.480508e+00, 4.730553e+11, 1 | 3.48051, 4.73055e+11, 1 |
| DZ0 | 20 | 3.4805 | 3.48051 | 3.87656 | 1 | 4.7306e+11 | 88.7808 | 0.05475 | 0.000183422 | 1.91129e+14 | -0.112009 | 1.78479e+07 | 1.98789e+07 | 2.30629e+08 | 25959.1 | 3.48051, 4.73055e+11, 1 | 3.48051, 4.73055e+11, 1 |
| DM0 | 0 | 12.0076 | 12.0076 | 14.1642 | 3 | 6.7509e+11 | 80.0391 | 0.0374973 | 0.000174261 | 2.21248e+14 | 0 | 1.92028e+07 | 2.26516e+07 | 9.06639e+08 | 87985.3 | 12.0076, 6.75093e+11, 3 | 12.0272, 6.75093e+11, 3 |
| DM0 | 20 | 12.0076 | 12.0076 | 14.1642 | 3 | 6.7509e+11 | 99.9609 | 0.0784035 | 0.000174261 | 2.21248e+14 | 0.00429395 | 1.92028e+07 | 2.26516e+07 | 9.06639e+08 | 87985.3 | 12.0076, 6.75093e+11, 3 | 12.0272, 6.75093e+11, 3 |
| DZs | 22155 | 3.4805 | 3.48051 | 3.80097 | 1 | 4.7306e+11 | 89.6127 | 0.03825 | 0.000183423 | 1.91129e+14 | 0.152465 | 1.78479e+07 | 1.98789e+07 | 2.30629e+08 | 25959.3 | 3.48051, 4.73055e+11, 1 | 3.48051, 4.73055e+11, 1 |
| DZs | 22175 | 3.4805 | 3.48051 | 3.7741 | 1 | 4.7306e+11 | 89.6701 | 0.03775 | 0.000183423 | 1.91129e+14 | 0.137165 | 1.78479e+07 | 1.98789e+07 | 2.30629e+08 | 25959.3 | 3.48051, 4.73055e+11, 1 | 3.48051, 4.73055e+11, 1 |
| DMs | 6280 | 10.5826 | 10.5981 | 10.5982 | 1 | 9.1353e+11 | 92.9297 | 0.218848 | 3.47201e-09 | 1.94481e+13 | 2.93801e+07 | 5.69328e+06 | 5.69359e+06 | 1.23904e+09 | 4.60392 | 12.0072, 6.75093e+11, 3 | 10.592, 9.13526e+11, 1 |
| DMs | 6300 | 10.7108 | 10.7167 | 10.7169 | 1 | 9.1229e+11 | 92.9297 | 0.218848 | 3.27703e-09 | 1.90454e+13 | 2.90038e+07 | 5.63403e+06 | 5.6344e+06 | 1.23736e+09 | 4.43771 | 12.0072, 6.75093e+11, 3 | 10.7167, 9.12288e+11, 1 |
| DT15 | 5280 | 13.9435 | 22.0042 | 22.0095 | 1 | 2.4651e+12 | 61.6993 | 0.463262 | 8.43787e-13 | 5.6033e+15 | -1.33658e+08 | 9.66377e+07 | 9.66932e+07 | 1.68957e+10 | 0.00478027 | 77.4412, 2.65815e+12, 1 | 21.7512, 2.46511e+12, 1 |
| DT15 | 5300 | 14.7812 | 1.804085e+01 | 1.804638e+01 | 3 | 2.4651e+12 | 6.169932e+01 | 4.693981e-01 | 3.714662e-12 | 1.624963e+15 | 1.693511e+08 | 5.204106e+07 | 5.210899e+07 | 1.331776e+10 | 1.624155e-02 | 7.726963e+01, 2.658148e+12, 1 | 1.804085e+01, 2.465107e+12, 3 |
| DTq | 2555 | 96.8409 | 96.8409 | 123.301 | 3 | 1.6909e+12 | 45.176 | 0.303728 | 9.9815e-05 | 1.92347e+14 | -32.6728 | 1.79047e+07 | 2.27969e+07 | 7.35892e+09 | 373621 | 96.8409, 1.69091e+12, 3 | 100.884, 1.69091e+12, 3 |
| DTq | 2575 | 96.8408 | 96.8408 | 123.301 | 3 | 1.6909e+12 | 45.176 | 0.316 | 9.98149e-05 | 1.92347e+14 | -33.8057 | 1.79047e+07 | 2.27969e+07 | 7.35892e+09 | 373621 | 96.8408, 1.69091e+12, 3 | 100.884, 1.69091e+12, 3 |

### What sets it (over all printed cycles of each arm)
- DZ0: 11 lines; dir counts {'1': 11} (1 r, 2 theta, 3 phi); binder cells (gid,kji: count) {('0', '3,3,3'): 2, ('1', '18,32,3'): 1, ('3', '35,39,3'): 1}; r 4.7306e+11 - 4.7306e+11; rho 1.83e-04 - 1.83e-04; floor flag {'0': 11}; |v| max 1.31e+00, cgas 1.78e+07 - 1.78e+07, ceff 1.99e+07 - 1.99e+07 (ceff > cgas: radiation signal in the speed); cell width h 2.306e+08 - 2.306e+08; cell_dt/cfl*dt_hyd ratio 1.0000
- DM0: 11 lines; dir counts {'3': 11} (1 r, 2 theta, 3 phi); binder cells (gid,kji: count) {('0', '30,3,3'): 2, ('1', '14,130,3'): 2, ('2', '100,3,3'): 2}; r 6.7509e+11 - 6.7509e+11; rho 1.74e-04 - 1.74e-04; floor flag {'0': 11}; |v| max 1.87e-02, cgas 1.92e+07 - 1.92e+07, ceff 2.27e+07 - 2.27e+07 (ceff > cgas: radiation signal in the speed); cell width h 9.066e+08 - 9.066e+08; cell_dt/cfl*dt_hyd ratio 1.0000
- DZs: 11 lines; dir counts {'1': 11} (1 r, 2 theta, 3 phi); binder cells (gid,kji: count) {('2', '15,53,3'): 3, ('2', '16,51,3'): 2, ('2', '16,53,3'): 2}; r 4.7306e+11 - 4.7306e+11; rho 1.83e-04 - 1.83e-04; floor flag {'0': 11}; |v| max 4.50e-01, cgas 1.78e+07 - 1.78e+07, ceff 1.99e+07 - 1.99e+07 (ceff > cgas: radiation signal in the speed); cell width h 2.306e+08 - 2.306e+08; cell_dt/cfl*dt_hyd ratio 1.0000
- DMs: 10 lines; dir counts {'1': 10} (1 r, 2 theta, 3 phi); binder cells (gid,kji: count) {('3', '35,40,226'): 8, ('3', '35,40,225'): 2}; r 9.1229e+11 - 9.1353e+11; rho 3.28e-09 - 3.84e-09; floor flag {'0': 10}; |v| max 2.95e+07, cgas 5.57e+06 - 5.71e+06, ceff 5.57e+06 - 5.71e+06 (ceff > cgas: radiation signal in the speed); cell width h 1.237e+09 - 1.239e+09; cell_dt/cfl*dt_hyd ratio 1.0000
- DT15: 10 lines; dir counts {'1': 2, '3': 5, '2': 3} (1 r, 2 theta, 3 phi); binder cells (gid,kji: count) {('0', '78,50,58'): 4, ('0', '78,51,59'): 3, ('0', '79,50,58'): 2}; r 2.4651e+12 - 2.4821e+12; rho 7.82e-13 - 1.91e-11; floor flag {'0': 10}; |v| max 2.79e+08, cgas 1.47e+07 - 9.66e+07, ceff 1.50e+07 - 9.67e+07 (ceff > cgas: radiation signal in the speed); cell width h 1.327e+10 - 1.690e+10; cell_dt/cfl*dt_hyd ratio 1.0000
- DTq: 10 lines; dir counts {'3': 10} (1 r, 2 theta, 3 phi); binder cells (gid,kji: count) {('0', '52,3,3'): 7, ('0', '53,3,3'): 2, ('0', '54,3,3'): 1}; r 1.6909e+12 - 1.6909e+12; rho 9.98e-05 - 9.98e-05; floor flag {'0': 10}; |v| max 3.39e+01, cgas 1.79e+07 - 1.79e+07, ceff 2.28e+07 - 2.28e+07 (ceff > cgas: radiation signal in the speed); cell width h 7.359e+09 - 7.359e+09; cell_dt/cfl*dt_hyd ratio 1.0000

### T2 production logs (read only): binder over the first chunk, where DTREP was on (every 10 cycles)
- m35mms: 412 lines, cycles 0-2040; dir counts {'3': 412}; first r 6.75093e+11 dir 3 rho 0.000174261 h 9.06639e+08 cfl*dt_hyd 12.0076; last r 6.75093e+11 dir 3 rho 0.000174263 h 9.06639e+08 cfl*dt_hyd 12.0076
- m35zams: 1352 lines, cycles 0-6740; dir counts {'1': 1352}; first r 4.73055e+11 dir 1 rho 0.000183422 h 2.30629e+08 cfl*dt_hyd 3.48051; last r 4.73055e+11 dir 1 rho 0.000183423 h 2.30629e+08 cfl*dt_hyd 3.48051
- t9l62: 98 lines, cycles 0-470; dir counts {'3': 98}; first r 1.69091e+12 dir 3 rho 9.98281e-05 h 7.35892e+09 cfl*dt_hyd 96.8523; last r 1.690906e+12 dir 3 rho 9.983636e-05 h 7.358921e+09 cfl*dt_hyd 9.684894e+01
- bsg: 268 lines, cycles 0-1320; dir counts {'3': 268}; first r 1.40233e+12 dir 3 rho 5.62878e-05 h 4.97371e+09 cfl*dt_hyd 167.506; last r 1.40233e+12 dir 3 rho 5.7607e-05 h 4.97371e+09 cfl*dt_hyd 169.288

### What these lines do NOT give
- The share of cells within 5 % of the limit, and the minimum dt per direction over the mesh: DTREP prints one cell (the binder) and DTREG the worst cell of the theta-edge rows and of the rest. So the "second limiting class" (e.g. the radial-only dt) is not in the log and I do not estimate it. cfl*dt_gasonly is the same cell without the radiation signal, not the second-worst cell.
- The same cell in all printed cycles means the binder is fixed; a changing (gid, kji) is listed above by its count.
- DZ0: 2026-11-24 17:33:49 rc=0 (srun) ended: Terminating on cycle limit wall=31 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean st
- DM0: 2026-11-24 17:34:34 rc=0 (srun) ended: Terminating on cycle limit wall=45 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean st
- DZs: 2026-11-24 17:35:06 rc=0 (srun) ended: Terminating on cycle limit wall=32 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean st
- DMs: 2026-11-24 17:36:10 rc=0 (srun) ended: Terminating on cycle limit wall=64 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean st
- DT15: 2026-11-24 17:38:33 rc=0 (srun) ended: Terminating on cycle limit wall=143 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean st
- DTq: 2026-11-24 17:40:43 rc=0 (srun) ended: Terminating on cycle limit wall=130 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean st
- All DTREP / DTREG lines: caltech_results/DT_LIMIT/<arm>/dtrep.txt (printed every 5 cycles).
