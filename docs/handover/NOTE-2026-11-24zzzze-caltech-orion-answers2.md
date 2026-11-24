# NOTE 2026-11-24zzzj caltech -> orion: stop of RSG + 4 stars and the test batch (your 24zzzh / 24zzzi)
# NOTE 2026-11-24zzzze caltech -> orion: answers to 24zzzzb / zzzzc / zzzzd

## S2 (zzzzd 1): it ran; the empty post was my script's error
- The scripts' star-name variable overwrote the results directory, so the poster looked in m35zams/S2off. The arms themselves are fine (Slurm: t_S2off COMPLETED 0:0 1:53, t_S2on COMPLETED 0:0 1:52; log ends 'Terminating on cycle limit'; 600 CLEAN lines, 300 SDO lines, 300 hst rows each). The S2 table is reposted below (PASS); clean lines are now all_lines.txt.gz (the empty 20-byte files are removed) and every post checks file sizes and that a verdict exists (a failed check is written into the post).
### S2off last 20 lines of log.0
```
MeshBlock-cycles = 1200
cpu time used  = 9.287979e+01
main loop time (before the final outputs) = 9.209982e+01
zone-cycles/cpu_second = 2.032129e+07
particle-updates/cpu_second = 0.000000e+00
<hydro>/split_bval_fast: energy-only refreshes 0, full refreshes 0 (rank 0)
<rad_vet> closure = vet_col: 0.000000e+00 tensor builds (0.000000e+00 skipped), 0.000000e+00 s (rank 0, fenced), 0.000000e+00 ms per build
<rad_vet> vet_col_lat: 3.000000e+02 sweeps, 0.000000e+00 s (rank 0, fenced), 0.000000e+00 ms per sweep; lateral reads clamped to the ghost band (rank 0, all sweeps): 0.00
<rad_vet> vet_gd: rank 0 fenced seconds: source 1.205400e+00, sweeps 1.825563e+01 (of which per-shell halo+walls 1.800697e+01, of which MPI wait 4.979713e+00), tables (ro
<rad_vet> vet_gd halo (rank 0): exchanges 6.330000e+04, bytes sent 9.829219e+11 (1.552799e+07 per exchange), post 0.000000e+00 s, partners 3; band-clamped lateral reads (
<rad_vet> vet_gd global band (rank 0): compact, intensity values sent 1.078837e+11 of dense 3.777829e+11 (2.855707e-01); list builds 0.000000e+00 (0.000000e+00 s)
<rad_vet> vet_gd ragged band (rank 0): out-of-slab accesses 0.000000e+00
<rad_vet> vet_gd_band_ring (rank 0): sweeps 3.000000e+02, dense band saves 0.000000e+00, save store 2.176951e+00 GB
<rad_vet> (rank 0) vet_gd halo lists 0.000000e+00 made, 0.000000e+00 MB of 0.000000e+00 (key 4096, GPU: <= free/2); exchanges by list 0.000000e+00, dense 0.000000e+00
<rad_vet> (rank 0) exact-list batches 0.000000e+00 (group 0, 0.000000e+00 s, list reallocs 0.000000e+00); table prefetch 1 swapped in 0.000000e+00 of the rotations
<rad_vet> implicit energy solve: cells=0.000000e+00 iterations mean=0.000000e+00 max=0.000000e+00 FAILURES=0.000000e+00 (no bracket: 0.000000e+00)
<rad_vet> vet-positivity: implicit_pos_gas=1 (frac 1.000000e-03) cell-solves=0.000000e+00 energy moved rad->gas=0.000000e+00 | E floor raises cell-solves=0.000000e+00 ene
<rad_vet> time_scheme=hesdirk2: stage steps=3.000000e+02 backward-Euler steps=0.000000e+00 stage fallbacks=0.000000e+00 vet clips=0.000000e+00 face_fct calls=6.000000e+02
<rad_vet> time2_auto_be=16: skipped-stage BE steps=0.000000e+00 stage retries=0.000000e+00 (succeeded 0.000000e+00)
<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 8.812223e+01 | failed stage attempts 0.000000e+00 | their BE redo 0.000000e+00 | auto_be steps 0.000000e+00 | 
```
### S2on last 20 lines of log.0
```
MeshBlock-cycles = 1200
cpu time used  = 9.343978e+01
main loop time (before the final outputs) = 9.243166e+01
zone-cycles/cpu_second = 2.019950e+07
particle-updates/cpu_second = 0.000000e+00
<hydro>/split_bval_fast: energy-only refreshes 0, full refreshes 0 (rank 0)
<rad_vet> closure = vet_col: 0.000000e+00 tensor builds (0.000000e+00 skipped), 0.000000e+00 s (rank 0, fenced), 0.000000e+00 ms per build
<rad_vet> vet_col_lat: 3.000000e+02 sweeps, 0.000000e+00 s (rank 0, fenced), 0.000000e+00 ms per sweep; lateral reads clamped to the ghost band (rank 0, all sweeps): 0.00
<rad_vet> vet_gd: rank 0 fenced seconds: source 1.183640e+00, sweeps 1.829962e+01 (of which per-shell halo+walls 1.806647e+01, of which MPI wait 4.995134e+00), tables (ro
<rad_vet> vet_gd halo (rank 0): exchanges 6.330000e+04, bytes sent 9.829219e+11 (1.552799e+07 per exchange), post 0.000000e+00 s, partners 3; band-clamped lateral reads (
<rad_vet> vet_gd global band (rank 0): compact, intensity values sent 1.078837e+11 of dense 3.777829e+11 (2.855707e-01); list builds 0.000000e+00 (0.000000e+00 s)
<rad_vet> vet_gd ragged band (rank 0): out-of-slab accesses 0.000000e+00
<rad_vet> vet_gd_band_ring (rank 0): sweeps 3.000000e+02, dense band saves 0.000000e+00, save store 2.176951e+00 GB
<rad_vet> (rank 0) vet_gd halo lists 0.000000e+00 made, 0.000000e+00 MB of 0.000000e+00 (key 4096, GPU: <= free/2); exchanges by list 0.000000e+00, dense 0.000000e+00
<rad_vet> (rank 0) exact-list batches 0.000000e+00 (group 0, 0.000000e+00 s, list reallocs 0.000000e+00); table prefetch 1 swapped in 0.000000e+00 of the rotations
<rad_vet> implicit energy solve: cells=0.000000e+00 iterations mean=0.000000e+00 max=0.000000e+00 FAILURES=0.000000e+00 (no bracket: 0.000000e+00)
<rad_vet> vet-positivity: implicit_pos_gas=1 (frac 1.000000e-03) cell-solves=0.000000e+00 energy moved rad->gas=0.000000e+00 | E floor raises cell-solves=0.000000e+00 ene
<rad_vet> time_scheme=hesdirk2: stage steps=3.000000e+02 backward-Euler steps=0.000000e+00 stage fallbacks=0.000000e+00 vet clips=0.000000e+00 face_fct calls=6.000000e+02
<rad_vet> time2_auto_be=16: skipped-stage BE steps=0.000000e+00 stage retries=0.000000e+00 (succeeded 0.000000e+00)
<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 8.863439e+01 | failed stage attempts 0.000000e+00 | their BE redo 0.000000e+00 | auto_be steps 0.000000e+00 | 
```
## Z arms (zzzzd 2): my mistakes
- ZC (13:48) FATAL 'final_restart not found' (exit 1): the key was on the command line of a fresh deck. ZA / ZB (13:52): no directory was made (same variable error) so the step died at once (CANCELLED / FAILED, 9-10 s). Fixed: no final_restart on the command line; ZA runs alone and its first 5 cycles are checked before ZB and ZC. They run on 01-15 (7333036) now; the post follows.
## HG (zzzzd 3, zzzzc 1)
- HG arms: P1 and H0 on 02-09, H1 on 02-03, HT on 01-17 (7333072, 4 ranks each), +6-cycle vet dump in H0 and HT (output3/variable=vet, dcycle 6; the rest is deleted); seam_vet.py after. final_restart=false is kept only if the +1 step check passes with it. R used for seam_vet: R_ph 61.78 Rsun = 4.298e12 cm (header of hegiant.00139).
- hegiant.00205: 33776730112 of 68152836794 B at 13:56; the coordinator confirms the md5.
## 7x36 ray arm (zzzzc 3)
- It ran: preflight 13:46-13:47 (log ends on the cycle limit, 71 s), G2b 13:47-13:50 (t_G2b COMPLETED 0:0 2:33). NOTE zzzw has its table. **G2's abort line, quoted (pre_G2/log.0):** `  what():  Kokkos ERROR: Cuda memory space failed to allocate 60.41 GiB (label="vet_arena").` (and the same in G2/log.0 line 73).
## S5 / S3 / S1 (zzzzd 4)
- S3off ended (COMPLETED 0:0, cycle limit, c6959+300, CW 0). S3on ended (COMPLETED 0:0, cycle limit, CW 0); S5off / S5on running (about 250 of 300 cycles at 13:56, CW 0). S1 not started: it waits for S3 and S2 (zzzzb). S4 waits for S5.
