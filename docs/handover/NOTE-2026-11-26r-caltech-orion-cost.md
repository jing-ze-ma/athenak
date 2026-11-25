# NOTE 2026-11-26r caltech -> orion: module cost of the star runs (your 26r)

## bsg (11-25 12:30): rst bsg_bsg3d.00008.rst (hard link, c29127, state in rst), next3 c3f133b44 ad2f2007, own keys none; node 02-16
- bsg plain: 39 cycles timed, wall per cycle mean 0.7081 median 0.6368 min 0.5754 max 3.5351 s; dt 169.4271 (first) 164.6824 (last); its per cycle 32.9, solves per cycle 6.95; stage fallbacks 0.000000e+00, backward-Euler steps 0.000000e+00, DT_RETRY 0, CLEAN WARNING 0, FATAL 0; peak MiB per GPU 0:20389 1:20389 2:20389 3:20389; start cycle 29127 t 4.838566e+06
- per cycle wall [s]: 0.657 0.654 0.652 0.635 0.650 0.670 0.667 0.680 0.629 0.638 0.642 0.639 0.637 0.643 0.667 0.630 0.633 0.641 0.641 0.632 0.642 0.635 0.634 0.583 0.633 0.631 0.631 0.632 0.637 0.632 0.579 0.579 0.651 0.639 0.607 0.575 0.635 0.587 3.535
- <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 2.668885e+01 | failed stage attempts 0.000000e+00 | their BE redo 0.000000e+00 | auto_be steps 0.000000e+00 | other BE steps 0.000000e+00 | discarded dt_retry attempts 0.000000e+00
- <rad_vet> vet_gd_band_ring (rank 0): sweeps 4.000000e+01, dense band saves 0.000000e+00, save store 6.394461e-02 GB
- mpit arm (26 cycles, MPI per call site, no fences): bsg mpit rc=0 ended: Terminating on cycle limit wall=41 s mt files 4
- nsys arm: bsg nsys rc=0 ended: Terminating on cycle limit wall=70 s reports 4  
- files: caltech_results/COST26R/bsg/ (kname_r0..3, gaps_r0..3, stats_r0). Window: kname / gaps C0 = 2 (cycles 0 and 1 of the 25 dropped).

## m35zams (11-25 12:34): rst m35zams_m35zams.00010.rst (hard link, c223416, state in rst), next3 c3f133b44 ad2f2007, own keys none; node 02-16
- m35zams plain: 39 cycles timed, wall per cycle mean 0.2998 median 0.2715 min 0.2589 max 1.0571 s; dt 3.4805 (first) 3.4805 (last); its per cycle 35.0, solves per cycle 4.40; stage fallbacks 0.000000e+00, backward-Euler steps 0.000000e+00, DT_RETRY 0, CLEAN WARNING 0, FATAL 0; peak MiB per GPU 0:15535 1:15535 2:15535 3:15535; start cycle 223416 t 7.776028e+05
- per cycle wall [s]: 0.271 0.270 0.294 0.293 0.262 0.269 0.299 0.265 0.267 0.295 0.278 0.271 0.269 0.298 0.310 0.297 0.296 0.261 0.292 0.264 0.267 0.289 0.299 0.291 0.267 0.262 0.266 0.266 0.261 0.265 0.296 0.299 0.291 0.267 0.259 0.297 0.307 0.263 1.057
- <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 1.380257e+01 | failed stage attempts 0.000000e+00 | their BE redo 0.000000e+00 | auto_be steps 0.000000e+00 | other BE steps 0.000000e+00 | discarded dt_retry attempts 0.000000e+00
- <rad_vet> vet_gd_band_ring (rank 0): sweeps 4.000000e+01, dense band saves 0.000000e+00, save store 2.176951e+00 GB
- mpit arm (26 cycles, MPI per call site, no fences): m35zams mpit rc=0 ended: Terminating on cycle limit wall=25 s mt files 4
- nsys arm: m35zams nsys rc=0 ended: Terminating on cycle limit wall=51 s reports 4  
- files: caltech_results/COST26R/m35zams/ (kname_r0..3, gaps_r0..3, stats_r0). Window: kname / gaps C0 = 2 (cycles 0 and 1 of the 25 dropped).

## m35mms (11-25 12:41): rst m35mms_m35mms.00010.rst (hard link, c65876, state in rst), next3 c3f133b44 ad2f2007, own keys none; node 02-16
- m35mms plain: 39 cycles timed, wall per cycle mean 0.9261 median 0.8689 min 0.8047 max 3.2721 s; dt 11.9928 (first) 11.9925 (last); its per cycle 40.8, solves per cycle 7.95; stage fallbacks 0.000000e+00, backward-Euler steps 0.000000e+00, DT_RETRY 0, CLEAN WARNING 0, FATAL 0; peak MiB per GPU 0:26753 1:26513 2:26513 3:26513; start cycle 65876 t 7.776063e+05
- per cycle wall [s]: 0.921 0.921 0.934 0.867 0.900 0.860 0.911 0.891 0.896 0.865 0.909 0.866 0.824 0.869 0.876 0.870 0.807 0.874 0.863 0.805 0.860 0.811 0.899 0.825 0.872 0.825 0.848 0.829 0.831 0.812 0.890 0.875 0.860 0.905 0.875 0.822 0.809 0.870 3.272
- <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 3.578093e+01 | failed stage attempts 0.000000e+00 | their BE redo 0.000000e+00 | auto_be steps 0.000000e+00 | other BE steps 0.000000e+00 | discarded dt_retry attempts 0.000000e+00
- <rad_vet> vet_gd_band_ring (rank 0): sweeps 4.000000e+01, dense band saves 0.000000e+00, save store 7.167220e-02 GB
- mpit arm (26 cycles, MPI per call site, no fences): m35mms mpit rc=0 ended: Terminating on cycle limit wall=55 s mt files 4
- nsys arm: m35mms nsys rc=0 ended: Terminating on cycle limit wall=75 s reports 4  
- files: caltech_results/COST26R/m35mms/ (kname_r0..3, gaps_r0..3, stats_r0). Window: kname / gaps C0 = 2 (cycles 0 and 1 of the 25 dropped).

## t9l62 (11-25 12:54): rst t9l62_t9l62.00009.rst (hard link, c22866, state in rst), next3 c3f133b44 ad2f2007, own keys rad_vet/time2_gtime_slope=true rad_vet/clean_stag_thinf=true rad_vet/clean_stag_lat_sums=true rad_vet/clean_stag_tie_gid=true rad_vet/clean_stag_opac_drag=true rad_vet/vet_gd_allow_clamp=false hydro/vceil=3.0e8 rad_vet/clean_pc_sweep=rb rad_vet/clean_stag_mg_levels=2 time/cfl_number=0.3 rad_vet/clean_conv_max=1e-5; node 02-16
- t9l62 plain: 39 cycles timed, wall per cycle mean 2.7281 median 2.5148 min 2.4077 max 9.1700 s; dt 96.2045 (first) 96.1725 (last); its per cycle 43.7, solves per cycle 5.35; stage fallbacks 0.000000e+00, backward-Euler steps 0.000000e+00, DT_RETRY 0, CLEAN WARNING 0, FATAL 0; peak MiB per GPU 0:73381 1:73179 2:73179 3:73383; start cycle 22866 t 2.073641e+06
- per cycle wall [s]: 2.555 2.433 2.490 2.843 2.669 2.515 2.412 2.483 2.445 2.564 2.741 2.557 2.478 2.685 2.434 2.480 2.437 2.505 2.434 2.479 2.539 2.673 2.490 2.611 2.413 2.458 2.926 2.670 2.474 2.554 2.754 2.735 2.623 2.408 2.706 2.598 2.484 2.471 9.170
- <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 1.038470e+02 | failed stage attempts 0.000000e+00 | their BE redo 0.000000e+00 | auto_be steps 0.000000e+00 | other BE steps 0.000000e+00 | discarded dt_retry attempts 0.000000e+00
- <rad_vet> vet_gd_band_ring (rank 0): sweeps 4.000000e+01, dense band saves 0.000000e+00, save store 2.024421e-01 GB
- mpit arm (26 cycles, MPI per call site, no fences): t9l62 mpit rc=0 ended: Terminating on cycle limit wall=128 s mt files 4
- nsys arm: t9l62 nsys rc=137 ended:  wall=381 s reports 0 12:52:40 HANG: nsys step 7333035.262 cancelled (no log progress for 4 min) t9l62 kernel timer (nsys did not finish): Terminating on cycle limit 
- files: caltech_results/COST26R/t9l62/ (kname_r0..3, gaps_r0..3, stats_r0, kernel timer kt files: nsys did not finish; the timer fences after every kernel). Window: kname / gaps C0 = 2 (cycles 0 and 1 of the 25 dropped).
