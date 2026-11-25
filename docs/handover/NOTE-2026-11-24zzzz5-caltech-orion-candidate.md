# NOTE 2026-11-24zzzj caltech -> orion: stop of RSG + 4 stars and the test batch (your 24zzzh / 24zzzi)
# NOTE 2026-11-24zzzz5 caltech -> orion: candidate timing (your 24zzzzr-zzzzu): start state, what was changed, what is missing

- **24zzzzq items 1-3 cancelled** (steps 7333072.189, .191 and 7333073.125 cancelled by id; scripts stopped by recorded PID). The PS t9l62 arm ran on. Item 4 is in NOTE zzzz4.
- **PL (t9l62.00031) not run: the file no longer exists** (the dtdiag NOTE 23zz6 already said so on 11-23: '00031 no longer exists'). prod_fresh_1121b/t9l62/rst_p1_seg still holds 00037 (cycle 30085, t 1.5776e6, dt 117), 00038 and 00039; not substituted, as you said. Tell me if you want one of them.
- Binary of all arms: caltech-sdorder2-09d1e4df4 b8b0d79c45351fecfd04cee19c996605. Candidate keys named on the command line of every restart (read from CK() in lib.sh): time2_gtime_slope=true, thinf=true, lat_sums=true, tie_gid=true, opac_drag=true, allow_clamp=false, vceil=3.0e8, clean_pc_sweep=rb, mg 4 for m35mms else 2, cfl 0.3, dt_report=0. Fresh decks: copies of prod_rc4v_1124/decks with exactly these lines changed (from -> to), nothing else: thinf false->true, lat_sums / tie_gid / opac_drag 0->true, allow_clamp true->false (bsg, m35mms, t9l62), vceil bsg 1.0e8->3.0e8, mg m35mms 2->4, dt_report 10->0, nlim -1->100, outputs off; **added (absent in the deck): rad_vet/time2_gtime_slope = true**; final_restart not named (fresh: FATAL on the command line).
- My script error on the first try: two arms started without a restart file (a variable clash); their steps are FATAL 'Either an input or restart file must be specified' and are discarded; note the binary returns exit code 0 there, so I read the end from the log line. M35ZAMS_R and ZM2 ran on one node for ~10 s; both were cancelled and rerun one after the other.

## PL: late t9l62 restarts (t9l62.00031 is gone; not substituted; orion picks)
- Present in prod_fresh_1121b/t9l62/rst_p1_seg/ (28.0 GB each, nx1 512, 8 blocks 512x128x128, closure vet_col, clean_stag_flim false):
| file | cycle | t [s] | dt [s] | written (11-21) | binary that wrote it (chunks.txt) |
|---|---|---|---|---|---|
| t9l62.00037 | 30085 | 1.57756e6 | 117.3 | 17:15 (chunk 5, final) | athena_p1_bf4f917d4 (md5 d38d411c) |
| t9l62.00038 | 30321 | 1.59842e6 | 127.5 | 18:05 (chunk 6) | athena_p1_bf4f917d4 (md5 d38d411c) |
| t9l62.00039 | 30948 | 1.64163e6 | 41.8 | 18:38 (chunk 6) | athena_p1_bf4f917d4 (md5 d38d411c) |
- The decks name the J4 / old clean scheme (stored block <rad_m1>, read as <rad_vet>); the scheme names are in the rst headers (clean_stag_flim = false). 00040 and 00044-00047 also exist (00040 in rst_p1_seg, 00044-00047 in rst/).
- **Readable by b8b0d79c?** Load-only run (0 cycles, time/nlim = the file's own cycle, T2 restart args + candidate keys) on 02-16 (timing-suspect, not timed):
```
ld37 2026-11-24 15:30:46 rc=0 (srun) ended: Terminating on cycle limit wall=32 s CW 0 FATAL 0 clamp sweeps) 0
37 cycle 30085: 2026-11-24 15:30:46 rc=0 (srun) ended: Terminating on cycle limit wall=32 s CW 0 FATAL 0 clamp sweeps) 0 | 
ld38 2026-11-24 15:31:17 rc=0 (srun) ended: Terminating on cycle limit wall=30 s CW 0 FATAL 0 clamp sweeps) 0
38 cycle 30321: 2026-11-24 15:31:17 rc=0 (srun) ended: Terminating on cycle limit wall=30 s CW 0 FATAL 0 clamp sweeps) 0 | 
ld39 2026-11-24 15:31:47 rc=0 (srun) ended: Terminating on cycle limit wall=30 s CW 0 FATAL 0 clamp sweeps) 0
39 cycle 30948: 2026-11-24 15:31:47 rc=0 (srun) ended: Terminating on cycle limit wall=30 s CW 0 FATAL 0 clamp sweeps) 0 | 
```
- Exit status recorded: srun's code and the log's last state are in each line.

## Candidate timing: results (11-24 15:43). Binary b8b0d79c45351fecfd04cee19c996605 (caltech-sdorder2-09d1e4df4); settings as in the header (keys named on the command line for restarts; fresh decks with the lines listed there)
| arm | rows | s/cyc last 50 (cycles) | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | passes/cyc | CW/FATAL/inf/DT_RETRY | BE/fallback/retry | clamped reads | floor_pay sum | L_top/L_in end | end |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| bsg_R | 100 | 0.6598 (c6950-7000) ; loop 68.7 s | 147/149.3/153 | 15.2, 3.45 | 14.7, 3.78 | - | 0/0/0/0 | 0/0/0 | 0.000000e+00 | 0.00e+00 | 1.09672 | rc=0 (srun) ended: Terminating on cycle limit wall=93 s CW 0 FATAL 0 clamp sweep |
| bsg_F | 101 | 0.7289 (c50-100, loop time) ; loop 74.0 s | 164.6/167.3/168.5 | 20.7, 3.26 | 25.8, 3.94 | - | 0/0/0/0 | 1/0/0 | 0.000000e+00 | 0.00e+00 | 1.04016 | rc=0 (srun) ended: Terminating on cycle limit wall=149 s CW 0 FATAL 0 clamp swee |
| m35zams_R | 131 | 0.2971 (c22200-22250) ; loop 32.7 s | 3.481/3.481/3.481 | 18.9, 2.00 | 15.1, 2.15 | - | 0/0/0/0 | 0/0/0 | 0.000000e+00 | 0.00e+00 | 1.00015 | rc=0 (srun) ended: Terminating on cycle limit wall=53 s CW 0 FATAL 0 clamp sweep |
| m35zams_F | 101 | 0.2825 (c50-100, loop time) ; loop 29.7 s | 3.481/3.481/3.481 | 15.0, 2.33 | 15.4, 2.17 | - | 0/0/0/0 | 1/0/0 | 0.000000e+00 | 0.00e+00 | 1.00072 | rc=0 (srun) ended: Terminating on cycle limit wall=74 s CW 0 FATAL 0 clamp sweep |
| m35mms_R | 100 | 1.0972 (c7000-7050) ; loop 111.4 s | 12.01/12.01/12.01 | 37.3, 5.19 | 17.9, 4.04 | - | 0/0/0/0 | 0/0/0 | 0.000000e+00 | 8.73e+25 | 0.98817 | rc=0 (srun) ended: Terminating on cycle limit wall=140 s CW 0 FATAL 0 clamp swee |
| m35mms_F | 101 | 0.6543 (c50-100, loop time) ; loop 66.5 s | 12.01/12.01/12.01 | 8.6, 2.00 | 10.5, 2.27 | - | 0/0/0/0 | 1/0/0 | 0.000000e+00 | 0.00e+00 | 1.00160 | rc=0 (srun) ended: Terminating on cycle limit wall=154 s CW 0 FATAL 0 clamp swee |
| t9l62_R | 100 | 2.5820 (c2600-2650) ; loop 255.6 s | 96.84/96.84/96.84 | 5.4, 2.90 | 9.2, 2.51 | 5.41 | 0/0/0/0 | 0/0/0 | 0.000000e+00 | 3.69e+28 | 1.01512 | rc=0 (srun) ended: Terminating on cycle limit wall=315 s CW 0 FATAL 0 clamp swee |
| t9l62_F | 101 | 2.4168 (c50-100, loop time) ; loop 238.3 s | 96.85/96.85/96.85 | 8.0, 2.74 | 5.9, 1.99 | - | 1/0/0/0 | 1/0/0 | 0.000000e+00 | 1.46e+30 | 1.00957 | rc=0 (srun) ended: Terminating on cycle limit wall=543 s CW 1 FATAL 0 clamp swee |
| PO | 25 | 3.2653 (mean of loop, 25 cyc) ; loop 81.6 s | 13.94/22.54/27.66 | 14.1, 5.04 | 8.2, 4.08 | 9.12 | 0/0/0/0 | 0/0/0 | 0.000000e+00 | 1.52e+27 | 2.55212 | rc=0 (srun) ended: Terminating on cycle limit wall=149 s CW 0 FATAL 0 clamp swee |
| PQ | 25 | 2.5068 (mean of loop, 25 cyc) ; loop 62.7 s | 96.84/96.84/96.84 | 5.8, 2.60 | 6.7, 2.04 | 4.64 | 0/0/0/0 | 0/0/0 | 0 | 8.48e+27 | 1.01285 | rc=0 (srun) ended: Terminating on cycle limit wall=126 s CW 0 FATAL 0 clamp swee |
| PL | 1 | 4.2732 (c30100-30150) ; loop 271.7 s | 96.37/96.37/96.37 | 24.4, 10.12 | 17.2, 5.91 | 16.03 | 0/0/0/0 | 0/0/0 | 0.000000e+00 | 2.02e+23 | 1.59571 | rc=0 (srun) ended: Terminating on cycle limit wall=309 s CW 0 FATAL 0 clamp swee |
| ZM2 | 113 | 0.2998 (c22200-22250) ; loop 33.2 s | 3.481/3.481/3.481 | 18.9, 2.00 | 15.1, 2.15 | - | 0/0/0/0 | 0/0/0 | 0.000000e+00 | 0.00e+00 | 1.00015 | rc=0 (srun) ended: Terminating on cycle limit wall=50 s CW 0 FATAL 0 clamp sweep |
| ZM3 | 100 | 0.2988 (c22200-22250) ; loop 32.9 s | 3.481/3.481/3.481 | 14.5, 2.00 | 15.9, 2.41 | - | 0/0/0/0 | 0/0/0 | 0.000000e+00 | 0.00e+00 | 1.00015 | rc=0 (srun) ended: Terminating on cycle limit wall=48 s CW 0 FATAL 0 clamp sweep |
| ZM4 | 100 | 0.3013 (c22200-22250) ; loop 33.3 s | 3.481/3.481/3.481 | 13.6, 2.00 | 16.1, 2.50 | - | 0/0/0/0 | 0/0/0 | 0.000000e+00 | 0.00e+00 | 1.00015 | rc=0 (srun) ended: Terminating on cycle limit wall=54 s CW 0 FATAL 0 clamp sweep |

- Measured loop time per cycle: restart arms = elapsed over the last 50 cycles of the prints; **fresh arms = (main-loop time minus the time to cycle 50)/50**, because the last elapsed print of a fresh run includes the final rst write (its plain last-50 reading, 1.35 / 4.37 s for bsg / t9l62, is not a step time). M35ZAMS_R, ZM2 and ZM3 hst rows > 100 include rows of my cancelled overlapping try (same cycles repeated); the table values above use the prints only.
- **t9l62 fresh, CLEAN WARNING 1:** 'clean_init steady: not converged (residual 6.167e-08, tol 1.0e-08)' after 50 passes; the production T2 fresh run (thinf=false etc.) converged in 3 passes (residual 1.3e-9). It is the only CW in the candidate arms; the run went on and ended on the cycle limit.
- bsg_F, M35ZAMS_F, m35mms_F, t9l62_F show 1 backward-Euler step (the first step of a fresh run); no fallback, retry or DT_RETRY anywhere.
- PL (00037, older-scheme state p1_bf4f917d4 written 11-21; single 65-cycle run: the first 40 leave the scheme change behind; the profile is the last 10): dt 96.4, 16 pass lines per cycle against 4.6 (PQ) and 9.1 (PO); s/cyc 4.27 by the last-50 window (includes the first 40's part), hst rows 1 (its hst cadence is a time one). Its stored clean_conv_max is 1e-4 (others 1e-5), so its tol column reads 1.0e-4: this is a state with a looser stored tolerance, not the candidate's conv 1e-5; it was not named, so it is as stored.

### Per-pass means over the last 10 cycles (CLEAN PASS: it, de, dv, tol; tags = the last field: opac or the fe/fv value)
```
PO: cycles 5294-5303 (last 10), 9.3 pass lines per cycle
  t2s 2 pass  0: n  10 (of 10 cycles) | it 6.4 | de 8.18e-03 | dv 9.30e-07 | tol 1.00e-05 | tags 6.747e-16 6.795e-16 6.827e-16
  t2s 2 pass  1: n  10 (of 10 cycles) | it 1.3 | de 1.89e-02 | dv 4.96e-08 | tol 1.00e-05 | tags opac
  t2s 2 pass  2: n  10 (of 10 cycles) | it 3.4 | de 7.79e-04 | dv 1.90e-08 | tol 1.00e-05 | tags opac
  t2s 2 pass  3: n  10 (of 10 cycles) | it 1.9 | de 3.81e-05 | dv 5.78e-11 | tol 1.00e-05 | tags 1.054e-15 1.068e-15 1.113e-15
  t2s 2 pass  4: n  10 (of 10 cycles) | it 1.3 | de 3.00e-06 | dv 8.85e-12 | tol 1.00e-05 | tags 1.054e-15 1.068e-15 1.113e-15
  t2s 2 pass  5: n   1 (of 10 cycles) | it 0.0 | de 4.17e-12 | dv 1.05e-14 | tol 1.00e-05 | tags 6.813e-16
  t2s 3 pass  0: n  10 (of 10 cycles) | it 1.1 | de 1.52e-02 | dv 9.97e-07 | tol 1.00e-05 | tags 6.359e-16 6.535e-16 7.004e-16
  t2s 3 pass  1: n  10 (of 10 cycles) | it 2.5 | de 1.64e-03 | dv 2.58e-07 | tol 1.00e-05 | tags opac
  t2s 3 pass  2: n  10 (of 10 cycles) | it 2.7 | de 6.14e-05 | dv 1.43e-08 | tol 1.00e-05 | tags opac
  t2s 3 pass  3: n  10 (of 10 cycles) | it 1.0 | de 6.98e-06 | dv 5.36e-09 | tol 1.00e-05 | tags 6.331e-16 6.535e-16 7.079e-16
  t2s 3 pass  4: n   2 (of 10 cycles) | it 3.5 | de 2.50e-06 | dv 1.15e-09 | tol 1.00e-05 | tags 6.331e-16 7.593e-16
PQ: cycles 2568-2577 (last 10), 5.1 pass lines per cycle
  t2s 2 pass  0: n  10 (of 10 cycles) | it 0.0 | de 1.24e-18 | dv 6.92e-09 | tol 1.00e-11 | tags 2.001e-13
  t2s 2 pass  1: n  10 (of 10 cycles) | it 2.0 | de 2.03e-08 | dv 2.32e-11 | tol 1.01e-08 | tags opac
  t2s 2 pass  2: n  10 (of 10 cycles) | it 4.0 | de 1.75e-09 | dv 8.53e-13 | tol 9.95e-09 | tags opac
  t2s 3 pass  0: n  10 (of 10 cycles) | it 2.0 | de 4.06e-09 | dv 6.93e-09 | tol 2.03e-09 | tags 2.001e-13 2.002e-13
  t2s 3 pass  1: n  10 (of 10 cycles) | it 5.2 | de 1.16e-09 | dv 1.67e-12 | tol 1.55e-09 | tags opac
  t2s 3 pass  2: n   1 (of 10 cycles) | it 4.0 | de 1.40e-10 | dv 7.87e-14 | tol 8.45e-10 | tags opac
PL: cycles 30140-30149 (last 10), 15.0 pass lines per cycle
  t2s 2 pass  0: n  10 (of 10 cycles) | it 1.7 | de 2.35e-02 | dv 4.58e-05 | tol 1.00e-04 | tags 1.100e-12 1.354e-12 1.435e-12
  t2s 2 pass  1: n  10 (of 10 cycles) | it 1.0 | de 7.68e-03 | dv 2.10e-05 | tol 1.00e-04 | tags opac
  t2s 2 pass  2: n  10 (of 10 cycles) | it 1.0 | de 2.18e-02 | dv 1.46e-05 | tol 1.00e-04 | tags opac
  t2s 2 pass  3: n  10 (of 10 cycles) | it 1.1 | de 8.17e-03 | dv 1.09e-05 | tol 1.00e-04 | tags 1.165e-12 1.214e-12 1.248e-12
  t2s 2 pass  4: n  10 (of 10 cycles) | it 2.1 | de 1.99e-02 | dv 5.73e-06 | tol 1.00e-04 | tags 1.239e-12 1.256e-12 1.328e-12
  t2s 2 pass  5: n  10 (of 10 cycles) | it 6.2 | de 1.97e-03 | dv 5.22e-06 | tol 1.00e-04 | tags 1.247e-12 1.284e-12 1.329e-12
  t2s 2 pass  6: n  10 (of 10 cycles) | it 4.5 | de 9.43e-04 | dv 2.71e-06 | tol 1.00e-04 | tags 1.251e-12 1.287e-12 1.329e-12
  t2s 2 pass  7: n   8 (of 10 cycles) | it 3.8 | de 1.04e-03 | dv 2.38e-06 | tol 1.00e-04 | tags 1.251e-12 1.288e-12 1.329e-12
  t2s 2 pass  8: n   6 (of 10 cycles) | it 5.0 | de 6.54e-04 | dv 2.26e-06 | tol 1.00e-04 | tags 1.251e-12 1.288e-12 1.329e-12
  t2s 2 pass  9: n   5 (of 10 cycles) | it 3.2 | de 3.42e-04 | dv 9.77e-07 | tol 1.00e-04 | tags 1.251e-12 1.288e-12 1.329e-12
  t2s 2 pass 10: n   3 (of 10 cycles) | it 1.0 | de 4.92e-04 | dv 5.91e-07 | tol 1.00e-04 | tags 1.251e-12 1.288e-12 2.396e-12
  t2s 2 pass 11: n   2 (of 10 cycles) | it 2.0 | de 1.61e-04 | dv 7.70e-07 | tol 1.00e-04 | tags 1.251e-12 1.288e-12
  t2s 2 pass 12: n   1 (of 10 cycles) | it 1.0 | de 2.53e-04 | dv 2.63e-07 | tol 1.00e-04 | tags 1.251e-12
  t2s 2 pass 13: n   1 (of 10 cycles) | it 1.0 | de 1.75e-05 | dv 2.33e-07 | tol 1.00e-04 | tags 1.251e-12
  t2s 3 pass  0: n  10 (of 10 cycles) | it 2.9 | de 1.01e-02 | dv 5.76e-05 | tol 1.00e-04 | tags 1.418e-12 1.440e-12 1.471e-12
  t2s 3 pass  1: n  10 (of 10 cycles) | it 2.1 | de 1.77e-03 | dv 8.31e-06 | tol 1.00e-04 | tags opac
  t2s 3 pass  2: n  10 (of 10 cycles) | it 3.6 | de 1.20e-03 | dv 5.11e-06 | tol 1.00e-04 | tags opac
  t2s 3 pass  3: n   8 (of 10 cycles) | it 4.4 | de 3.35e-03 | dv 4.56e-06 | tol 1.00e-04 | tags 1.482e-12 1.496e-12 1.501e-12
  t2s 3 pass  4: n   8 (of 10 cycles) | it 2.9 | de 6.68e-04 | dv 3.98e-06 | tol 1.00e-04 | tags 1.482e-12 1.496e-12 1.539e-12
  t2s 3 pass  5: n   5 (of 10 cycles) | it 2.4 | de 2.86e-04 | dv 6.95e-07 | tol 1.00e-04 | tags 1.496e-12 1.567e-12 1.570e-12
  t2s 3 pass  6: n   3 (of 10 cycles) | it 1.3 | de 1.01e-05 | dv 9.02e-08 | tol 1.00e-04 | tags 1.496e-12 1.570e-12
```
- PS t9l62 (100 cycles, 00015): rb 3.0685 vs sym 3.1325 s/cyc, its t2s2 14.7 -> 13.1, t2s3 8.7 -> 8.2, dt mean -6.3 % (different trajectories), cost per simulated second 1.089 x, L_top at matched time max 1.0e-4; not harmless by the bounds (cost and dt). Earlier PS m35mms: cost 1.114 x.

## QUESTION before any fresh t9l62 start: the candidate clean_init CW
- Candidate fresh t9l62 (b8b0d79c, thinf / lat_sums / tie_gid / opac_drag true, gtime_slope true, allow_clamp false, vceil 3e8): 'CLEAN WARNING clean_init steady: not converged (residual 6.167e-08, tol 1.0e-08)', 50 passes. The T2 production fresh t9l62 (thinf = false, the old values): converged in 3 passes, residual 1.3e-9, no CW. The other three stars (bsg, M35ZAMS, m35mms) have no such CW in the candidate fresh start.
- The run continued and finished 100 cycles (CW 1 only). I have not looked into the cause; **what explains it, and is a fresh t9l62 start on these settings acceptable? Please answer before any fresh t9l62 start.**

## t9l62 fresh: which key makes clean_init steady warn (your 24zzzzz), nlim = 1, candidate binary b8b0d79c, deck copies (one line changed back per arm)
| arm | changed back | steady init | CW | end |
|---|---|---|---|---|
| A | thinf=false | converged after 3 passes, residual 1.343e-09 (tol 1.0e-08; 0 formal so | 0 | srun rc rc=0, Terminating on cycle limit; last log line: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 0 | failed |
| B | lat_sums=false | NOT converged after 50 passes, residual 6.166e-08 (tol 1.0e-08; 0 form
not converged (residual 6.166e-08, tol 1.0e-08) | 1 | srun rc rc=0, Terminating on cycle limit; last log line: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 0 | failed |
| C | tie_gid=false | NOT converged after 50 passes, residual 6.167e-08 (tol 1.0e-08; 0 form
not converged (residual 6.167e-08, tol 1.0e-08) | 1 | srun rc rc=0, Terminating on cycle limit; last log line: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 0 | failed |
| D | opac_drag=false | NOT converged after 50 passes, residual 6.167e-08 (tol 1.0e-08; 0 form
not converged (residual 6.167e-08, tol 1.0e-08) | 1 | srun rc rc=0, Terminating on cycle limit; last log line: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 0 | failed |
- **Arm A (thinf = false) converges in 3 passes (1.343e-9), CW 0, like the T2 production fresh; B, C, D each warn after 50 passes (6.17e-8).** So thinf = true is the key that makes the steady init not converge to 1e-8 (the production T2 and the candidate differ in thinf, allow_clamp, gtime_slope, the binary and dt_report; only thinf matters for the init). E (gtime_slope = false) is running, F (sym) skipped as you said; E's result follows. Production deck vs candidate deck differences: lat_sums, tie_gid, opac_drag, thinf, allow_clamp (true -> false), time2_gtime_slope (added), vceil (3.0e8 already in both), outputs, dt_report, nlim. No other differing line. Not tested: allow_clamp and the binary on their own (A differs from the T2 production by those and still converges).
- **Arm E** (time2_gtime_slope = false): NOT converged after 50 passes, residual 6.167e-08, CW 1 (srun rc 0, cycle limit). So of A-E only A (thinf = false) converges; B, C, D, E all warn. F skipped (zzzza1).
