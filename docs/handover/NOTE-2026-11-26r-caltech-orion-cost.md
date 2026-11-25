# NOTE 2026-11-26r caltech -> orion: module cost of the star runs (your 26r)

## bsg (11-25 12:30): rst bsg_bsg3d.00008.rst (hard link, c29127, state in rst), next3 c3f133b44 ad2f2007, own keys none; node 02-16
- bsg plain: 39 cycles timed, wall per cycle mean 0.7081 median 0.6368 min 0.5754 max 3.5351 s; dt 169.4271 (first) 164.6824 (last); its per cycle 32.9, solves per cycle 6.95; stage fallbacks 0.000000e+00, backward-Euler steps 0.000000e+00, DT_RETRY 0, CLEAN WARNING 0, FATAL 0; peak MiB per GPU 0:20389 1:20389 2:20389 3:20389; start cycle 29127 t 4.838566e+06
- per cycle wall [s]: 0.657 0.654 0.652 0.635 0.650 0.670 0.667 0.680 0.629 0.638 0.642 0.639 0.637 0.643 0.667 0.630 0.633 0.641 0.641 0.632 0.642 0.635 0.634 0.583 0.633 0.631 0.631 0.632 0.637 0.632 0.579 0.579 0.651 0.639 0.607 0.575 0.635 0.587 3.535
- <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 2.668885e+01 | failed stage attempts 0.000000e+00 | their BE redo 0.000000e+00 | auto_be steps 0.000000e+00 | other BE steps 0.000000e+00 | discarded dt_retry attempts 0.000000e+00
- <rad_vet> vet_gd_band_ring (rank 0): sweeps 4.000000e+01, dense band saves 0.000000e+00, save store 6.394461e-02 GB
- mpit arm (26 cycles, MPI per call site, no fences): bsg mpit rc=0 ended: Terminating on cycle limit wall=41 s mt files 4
- nsys arm: bsg nsys rc=0 ended: Terminating on cycle limit wall=70 s reports 4  
- files: caltech_results/COST26R/bsg/ (kname_r0..3, gaps_r0..3, stats_r0). Window: kname / gaps C0 = 2 (cycles 0 and 1 of the 25 dropped).
