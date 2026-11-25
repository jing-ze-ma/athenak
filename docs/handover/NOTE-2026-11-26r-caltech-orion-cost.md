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

## Cost per module, all stars (26r + 26s)
Setting: binary next3 c3f133b44 (ad2f2007), 4 ranks on 02-16, one H200 each (as production); rst = hard link of the newest complete production rst: bsg3d.00008 (c29127, t 4.84e6, bsg's last-but-one), m35zams.00010 (c223416, t 7.78e5), m35mms.00010 (c65876, t 7.78e5), t9l62.00009 (c22866, t 2.07e6): late, developed states. No vet rotation in these stars, so the rays are rebuilt every cycle (there is no cycle without a ray rebuild).

Plain 40 cycles (wall per cycle mean / median; the first timed cycle of bsg, m35mms, t9l62 is 3-9 s because of the set-up after the restart):
- bsg plain: 39 cycles timed, wall per cycle mean 0.7081 median 0.6368 min 0.5754 max 3.5351 s; dt 169.4271 (first) 164.6824 (last); its per cycle 32.9, solves per cycle 6.95; stage fallbacks 0.000000e+00, backward-Euler steps 0.000000e+00, DT_RETRY 0, CLEAN WARNING 0, FATAL 0; peak MiB per GPU 0:20389 1:20389 2:20389 3:20389; start cycle 29127 t 4.838566e+06
- m35zams plain: 39 cycles timed, wall per cycle mean 0.2998 median 0.2715 min 0.2589 max 1.0571 s; dt 3.4805 (first) 3.4805 (last); its per cycle 35.0, solves per cycle 4.40; stage fallbacks 0.000000e+00, backward-Euler steps 0.000000e+00, DT_RETRY 0, CLEAN WARNING 0, FATAL 0; peak MiB per GPU 0:15535 1:15535 2:15535 3:15535; start cycle 223416 t 7.776028e+05
- m35mms plain: 39 cycles timed, wall per cycle mean 0.9261 median 0.8689 min 0.8047 max 3.2721 s; dt 11.9928 (first) 11.9925 (last); its per cycle 40.8, solves per cycle 7.95; stage fallbacks 0.000000e+00, backward-Euler steps 0.000000e+00, DT_RETRY 0, CLEAN WARNING 0, FATAL 0; peak MiB per GPU 0:26753 1:26513 2:26513 3:26513; start cycle 65876 t 7.776063e+05
- t9l62 plain: 39 cycles timed, wall per cycle mean 2.7281 median 2.5148 min 2.4077 max 9.1700 s; dt 96.2045 (first) 96.1725 (last); its per cycle 43.7, solves per cycle 5.35; stage fallbacks 0.000000e+00, backward-Euler steps 0.000000e+00, DT_RETRY 0, CLEAN WARNING 0, FATAL 0; peak MiB per GPU 0:73381 1:73179 2:73179 3:73383; start cycle 22866 t 2.073641e+06

Arms per star: plain (40 cycles), MPI timer (libmpit, no fences), nsys (bsg, m35zams, m35mms; t9l62's hung: the nsys step was cancelled by its id after 4 min without log progress, no reports, so its numbers come from the kernel timer), named kernel timer (kt, global fencing; 25 cycles minus a 6-cycle run so that start-up is out; it lists the kernels by name, which nsys lumps as "void" / "bool").

How the modules are cut (names of the kernel timer): rays = vet_vgd*, vet_vlat*, vet_lay*, cln_gd*, cln_raymom, cln_rayj, cln_raycut, cln_clositer, vet_t2_avg*; hydro = hflux*, hyd*, HydroNudt, h_update, hydro_am, FOFC*, dfl*, duale*, wbcache, spsrc, star_*, etotgrav*, HistSums*; rad-vet = all other cln_*, cls_*, vet_*; others = the rest (view copies and fills, SendBuff / RecvBuff / Seam* of the hydro boundary exchange). Pack / copy / unpack = cln_hd_*, cln_ghost*, vet_vgd_hl_*, vet_vgd_hxb_*, vet_vgd_halo, vet_vgd_unpack, vet_vgd_shl, SendBuff, RecvBuff. Kernels that could be calculation or communication (counted as calculation here): vet_vgd_ring (the ring sweep and shell copies are one kernel), cln_cslat and cln_cslat_sh, cln_cstie*, cln_tiet, cln_fwb*, cln_slot, cln_apply_st_sh, cln_fs3_sh; their sizes are printed per star. Remarks: gravity and the hydro boundary values appear partly in others.


### bsg
kernel timer (fenced, 4 ranks, 25-cycle minus 6-cycle run, per 19 cycles (start-up and end removed)): ms per cycle, mean of the ranks
| module | calculation | pack / copy / unpack / ghost / request | Kokkos view copies, init | sum | share of kernel time |
|---|---|---|---|---|---|
| hydro | 51.94 | 0.00 | 0.00 | 51.94 | 8.3 % |
| rays | 126.37 | 4.31 | 0.00 | 130.68 | 20.8 % |
| rad-vet | 424.71 | 19.18 | 0.00 | 443.89 | 70.6 % |
| others | 0.62 | 1.01 | 0.40 | 2.02 | 0.3 % |
- kernel sum 628.5 ms per cycle (fenced run wall per cycle includes the fence: see below); ambiguous names (cannot be split with certainty; counted as calculation above): cln_cslat 15.20; vet_vgd_ring 13.00; cln_fwb 7.73; cln_fwbs 5.74; cln_apply_st_sh 3.73; cln_cslat_sh 3.23; cln_cstie 2.04; cln_tiet 2.00; cln_fwbq 1.92; cln_cstie2 1.22; cln_slot 1.03; cln_fwbp 0.81
- largest "others" names: Kokkos::ViewFill-4D 0.64; SendBuff 0.51; RecvBuff 0.50; SeamFluxSendMPI 0.28; SeamStateRecvRaw 0.25; Kokkos::ViewFill-1D 0.10; Kokkos::View::initialization [cln_rmc] via memset 0.07; Kokkos::View::initialization [vetf2in] via memset 0.04
MPI timer (26 minus 6 cycles, per 20 cycles, 4 ranks, no fences): host time in MPI, ms per cycle (mean of the ranks), by module of the last kernel before the call: hydro 0.15, rays 2.04, rad-vet 25.17, others 3.19
  - sites (module, call): rad-vet Waitall cln_hd_pack 10.53; rad-vet Allreduce cln_lztot 5.45; rad-vet Waitall cln_apply_st 2.82; rays Waitall vet_vgd_halo 1.48; rad-vet Waitall cln_cslat 1.45; rad-vet Allreduce cln_extrap 1.39; others Test SendBuff 1.12; others Bcast Kokkos::View::initialization [_mirror] via memset 0.84; rad-vet Isend cln_hd_pack 0.66; rad-vet Waitall cln_fs3 0.61; others Wait SendBuff 0.61; rad-vet Waitall cln_hd_unpack 0.59; others Test SeamFluxSendMPI 0.54; rad-vet Allreduce cln_smax 0.50
  - non-MPI host waits in the same timer (fences 313.73, deep copies 126.02 ms per cycle, these overlap GPU work)
nsys (kernels unfenced; per real cycle: values divided by 1.14 = launches of CleanGdFaceSums#6 per printed cycle; printed cycles 21): mean of 4 ranks: wall 652.1, GPU busy 584.7, GPU idle 52.4 ms per cycle
  - idle by host activity (ms per cycle, mean of 4 ranks): MPI 12.34, MPI(part) 7.33, cuda-other 1.68, cuda-other(part) 0.81, host-compute 0.18, launch 1.31, launch(part) 1.51, memcpy-api 2.52, memcpy-api(part) 8.27, sync 11.86, sync(part) 4.57
- bsg plain: 39 cycles timed, wall per cycle mean 0.7081 median 0.6368 min 0.5754 max 3.5351 s; dt 169.4271 (first) 164.6824 (last); its per cycle 32.9, solves per cycle 6.95; stage fallbacks 0.000000e+00, backward-Euler steps 0.000000e+00, DT_RETRY 0, CLEAN WARNING 0, FATAL 0; peak MiB per GPU 0:20389 1:20389 2:20389 3:20389; start cycle 29127 t 4.838566e+06

### m35zams
kernel timer (fenced, 4 ranks, 25-cycle minus 6-cycle run, per 19 cycles (start-up and end removed)): ms per cycle, mean of the ranks
| module | calculation | pack / copy / unpack / ghost / request | Kokkos view copies, init | sum | share of kernel time |
|---|---|---|---|---|---|
| hydro | 18.52 | 0.00 | 0.00 | 18.52 | 7.3 % |
| rays | 66.52 | 0.00 | 0.00 | 66.52 | 26.2 % |
| rad-vet | 154.99 | 12.76 | 0.00 | 167.76 | 66.0 % |
| others | 0.02 | 0.81 | 0.49 | 1.32 | 0.5 % |
- kernel sum 254.1 ms per cycle (fenced run wall per cycle includes the fence: see below); ambiguous names (cannot be split with certainty; counted as calculation above): vet_vgd_ring 11.28; cln_cslat 4.36; cln_apply_st_sh 2.54; cln_fwb 1.81; cln_fwbs 1.77; cln_cslat_sh 1.71; cln_slot 0.97; cln_tiet 0.84; cln_fwbq 0.78; cln_cstie 0.64; cln_cstie2 0.60; cln_fwbp 0.36
- largest "others" names: RecvBuff 0.42; SendBuff 0.40; Kokkos::ViewCopy-4D 0.38; Kokkos::ViewFill-4D 0.29; Kokkos::View::initialization [cln_rmc] via memset 0.04; Kokkos::ViewFill-1D 0.04; radvetbc_x1 0.03; Kokkos::View::initialization [vet_vgd_rawrec_mirror] via memset 0.02
MPI timer (26 minus 6 cycles, per 20 cycles, 4 ranks, no fences): host time in MPI, ms per cycle (mean of the ranks), by module of the last kernel before the call: hydro 0.09, rays 19.27, rad-vet 17.02, others 2.47
  - sites (module, call): rays Waitall vet_vgdgb_kpack 14.84; rad-vet Waitall cln_hd_pack 6.02; rad-vet Allreduce cln_extrap 3.56; rad-vet Allreduce cln_lztot 2.66; rays Waitall vet_vgdgb_wflag 2.41; rad-vet Waitall cln_apply_st 2.21; others Isend SendBuff 1.57; rays Waitall vet_vgdgb_pack 1.33; others Test SendBuff 1.32; rad-vet Waitall cln_cslat 0.89; rad-vet Isend cln_hd_pack 0.53; rays Isend vet_vgdgb_kpack 0.42; others Bcast Kokkos::View::initialization [_mirror] via memset 0.37; others Wait SendBuff 0.36
  - non-MPI host waits in the same timer (fences 85.56, deep copies 35.52 ms per cycle, these overlap GPU work)
nsys (kernels unfenced; per real cycle: values divided by 8.00 = launches of CleanGdFaceSums#6 per printed cycle; printed cycles 1.5): mean of 4 ranks: wall 291.4, GPU busy 217.0, GPU idle 50.8 ms per cycle
  - idle by host activity (ms per cycle, mean of 4 ranks): MPI 18.03, MPI(part) 6.18, cuda-other 1.04, cuda-other(part) 0.40, host-compute 0.17, launch 0.84, launch(part) 1.40, memcpy-api 2.19, memcpy-api(part) 7.88, sync 8.83, sync(part) 3.86
- m35zams plain: 39 cycles timed, wall per cycle mean 0.2998 median 0.2715 min 0.2589 max 1.0571 s; dt 3.4805 (first) 3.4805 (last); its per cycle 35.0, solves per cycle 4.40; stage fallbacks 0.000000e+00, backward-Euler steps 0.000000e+00, DT_RETRY 0, CLEAN WARNING 0, FATAL 0; peak MiB per GPU 0:15535 1:15535 2:15535 3:15535; start cycle 223416 t 7.776028e+05

### m35mms
kernel timer (fenced, 4 ranks, 25-cycle minus 6-cycle run, per 19 cycles (start-up and end removed)): ms per cycle, mean of the ranks
| module | calculation | pack / copy / unpack / ghost / request | Kokkos view copies, init | sum | share of kernel time |
|---|---|---|---|---|---|
| hydro | 52.42 | 0.00 | 0.00 | 52.42 | 6.2 % |
| rays | 134.93 | 7.75 | 0.00 | 142.68 | 16.8 % |
| rad-vet | 627.54 | 23.57 | 0.00 | 651.12 | 76.4 % |
| others | 0.09 | 1.20 | 4.21 | 5.50 | 0.6 % |
- kernel sum 851.7 ms per cycle (fenced run wall per cycle includes the fence: see below); ambiguous names (cannot be split with certainty; counted as calculation above): cln_cslat 20.43; vet_vgd_ring 12.86; cln_fwb 9.77; cln_fwbs 7.99; cln_apply_st_sh 5.09; cln_cslat_sh 3.86; cln_cstie 2.77; cln_tiet 2.60; cln_fwbq 2.21; cln_cstie2 1.54; cln_slot 1.11; cln_fwbp 0.89
- largest "others" names: Kokkos::ViewCopy-4D 2.38; Kokkos::ViewFill-4D 0.77; RecvBuff 0.61; SendBuff 0.60; Kokkos::View::initialization [rst-tail-ahead] via memset 0.42; Kokkos::View::initialization [rst-cc-in] via memset 0.21; Kokkos::View::initialization [cc_outvar] via memset 0.14; Kokkos::ViewFill-1D 0.08
MPI timer (26 minus 6 cycles, per 20 cycles, 4 ranks, no fences): host time in MPI, ms per cycle (mean of the ranks), by module of the last kernel before the call: hydro 0.08, rays 11.29, rad-vet 28.64, others 0.02
  - sites (module, call): rad-vet Waitall cln_hd_pack 12.70; rays Waitall vet_vgd_ring 9.24; rad-vet Allreduce cln_lztot 6.25; rad-vet Waitall cln_apply_st 3.55; rad-vet Waitall cln_cslat 1.70; rad-vet Allreduce cln_extrap 1.40; rays Waitall vet_vgd_halo 1.24; others Test SendBuff 1.20; rad-vet Isend cln_hd_pack 0.76; rad-vet Waitall cln_fs3 0.66; rays Isend vet_vgd_ring 0.57; others Wait SendBuff 0.52; rad-vet Allreduce cln_smax 0.32; rad-vet Waitall cln_hd_unpack 0.31
  - non-MPI host waits in the same timer (fences 350.99, deep copies 118.99 ms per cycle, these overlap GPU work)
nsys (kernels unfenced; per real cycle: values divided by 1.00 = launches of CleanGdFaceSums#6 per printed cycle; printed cycles 24): mean of 4 ranks: wall 890.9, GPU busy 810.1, GPU idle 61.1 ms per cycle
  - idle by host activity (ms per cycle, mean of 4 ranks): MPI 16.40, MPI(part) 7.55, cuda-other 2.22, cuda-other(part) 1.05, host-compute 0.17, launch 1.80, launch(part) 1.65, memcpy-api 2.85, memcpy-api(part) 8.05, sync 15.57, sync(part) 3.80
- m35mms plain: 39 cycles timed, wall per cycle mean 0.9261 median 0.8689 min 0.8047 max 3.2721 s; dt 11.9928 (first) 11.9925 (last); its per cycle 40.8, solves per cycle 7.95; stage fallbacks 0.000000e+00, backward-Euler steps 0.000000e+00, DT_RETRY 0, CLEAN WARNING 0, FATAL 0; peak MiB per GPU 0:26753 1:26513 2:26513 3:26513; start cycle 65876 t 7.776063e+05

### t9l62
kernel timer (fenced, 4 ranks, 25-cycle minus 6-cycle run, per 19 cycles (start-up and end removed)): ms per cycle, mean of the ranks
| module | calculation | pack / copy / unpack / ghost / request | Kokkos view copies, init | sum | share of kernel time |
|---|---|---|---|---|---|
| hydro | 131.05 | 0.00 | 0.00 | 131.05 | 5.2 % |
| rays | 882.40 | 22.26 | 0.00 | 904.67 | 36.2 % |
| rad-vet | 1420.64 | 36.63 | 0.00 | 1457.27 | 58.3 % |
| others | 0.05 | 1.88 | 5.34 | 7.28 | 0.3 % |
- kernel sum 2500.3 ms per cycle (fenced run wall per cycle includes the fence: see below); ambiguous names (cannot be split with certainty; counted as calculation above): vet_vgd_ring 72.14; cln_cslat 51.65; cln_fwbs 16.41; cln_fwb 15.70; cln_apply_st_sh 11.66; cln_cslat_sh 7.83; cln_tiet 7.54; cln_cstie 5.74; cln_fwbq 2.82; cln_cstie2 2.68; cln_slot 1.13; cln_fwbp 0.87
- largest "others" names: Kokkos::ViewCopy-4D 2.85; Kokkos::ViewFill-4D 1.14; RecvBuff 0.96; SendBuff 0.92; Kokkos::View::initialization [d_out_var_mirror] via memset 0.38; Kokkos::View::initialization [cc_outvar] via memset 0.21; Kokkos::ViewFill-1D 0.20; Kokkos::View::initialization [vet_vgd_tab_w3] via memset 0.13
MPI timer (26 minus 6 cycles, per 20 cycles, 4 ranks, no fences): host time in MPI, ms per cycle (mean of the ranks), by module of the last kernel before the call: hydro 0.09, rays 68.66, rad-vet 36.07, others 1.41
  - sites (module, call): rays Waitall vet_vgd_ring 65.38; rad-vet Waitall cln_hd_copy 17.18; rad-vet Allreduce cln_lztot 5.51; rad-vet Waitall cln_apply_st 4.40; rays Waitall vet_vgd_halo 2.01; rad-vet Waitall cln_cslat 1.96; rad-vet Allreduce cln_extrap 1.84; rad-vet Waitall cln_fwbp 1.15; others Wait RecvBuff 0.99; rays Isend vet_vgd_ring 0.87; others Wait SendBuff 0.83; rad-vet Allreduce vet_t2_ffct_r 0.80; others Isend SendBuff 0.79; rad-vet Waitall cln_fs3 0.71
  - non-MPI host waits in the same timer (fences 832.92, deep copies 250.62 ms per cycle, these overlap GPU work)
nsys: the run hung, no kname / gaps (kernel timer and MPI timer only)
- t9l62 plain: 39 cycles timed, wall per cycle mean 2.7281 median 2.5148 min 2.4077 max 9.1700 s; dt 96.2045 (first) 96.1725 (last); its per cycle 43.7, solves per cycle 5.35; stage fallbacks 0.000000e+00, backward-Euler steps 0.000000e+00, DT_RETRY 0, CLEAN WARNING 0, FATAL 0; peak MiB per GPU 0:73381 1:73179 2:73179 3:73383; start cycle 22866 t 2.073641e+06
- Sum check (nsys stars, mean of 4 ranks, ms per cycle; wall = GPU kernels + GPU copies + idle): bsg 652.1 = 584.7 + 52.4 + 15.0 residual (2.3 %); m35zams 291.4 = 217.0 + 50.8 + 23.6 (8 %); m35mms 890.9 = 810.1 + 61.1 + 19.7 (2.2 %). The residual is GPU memcpy / memset time, which the kernel sum does not count. The nsys wall is 2-7 % above the plain median (bsg 0.637 s, m35zams 0.272 s, m35mms 0.869 s): nsys overhead. For t9l62 (no nsys): plain median 2.515 s; kernel sum 2500 ms (fenced), so GPU-idle and MPI are not separated there; the MPI timer gives 106 ms per cycle of host time in MPI (mean of the ranks).
- m35zams's printed nsys cycle count is 1.5 (the cycle detection merged 8 real cycles); values divided by 8.00, the launches of CleanGdFaceSums#6 (1 per cycle).
- Files: caltech_results/COST26R/cost.txt, kname_r*, gaps_r*, stats_r0 (nsys stars), kt*.gz and mt*.gz (per rank; kt6 / mt6 = the 6-cycle arms), cost26.py (the map; edit the regexes to regroup).
