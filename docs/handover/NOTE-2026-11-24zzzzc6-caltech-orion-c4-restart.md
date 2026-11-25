# NOTE 2026-11-24zzzzc6 caltech -> orion: c4 at the next restart of bsg, m35mms, M35ZAMS (your 24zzzzc6 item 2)

## bsg: PASS (11-24 18:07), bsg3d.00002.rst, cycle 7300 + 3, on 01-17
- C4: 2026-11-24 18:05:42 start bsg C4 hpc-sm-01-17 athena_c4_a000d803c md5=3ae130e97b13e5bdfa4a7f9a14370b9f ARGS: -r /resnick/groups/carnegie_poc/jingze/b1_1124/c4rst/bsg/rst_in/bsg3d.00002.rst rad_vet/time2_gtime_slope=true rad_vet/clean_stag_thinf=true rad_vet/clean_stag_lat_sums=true rad_vet/clean_stag_tie_gid=true rad_vet/clean_stag_opac_drag=true rad_vet/vet_gd_allow_clamp=false hydro/vceil=3.0e8 rad_vet/clean_pc_sweep=rb rad_vet/clean_stag_mg_levels=2 time/cfl_number=0.3 rad_vet/clean_conv_max=1e-5 output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30 time/nlim=7303
- C4: 2026-11-24 18:06:43 rc=0 (srun) ended: Terminating on cycle limit wall=61 s CW 0 FATAL 0 inf 0 DT_RETRY 0 stage fallback 0 clamped 0
- C1: 2026-11-24 18:06:44 start bsg C1 hpc-sm-01-17 athena_sdorder3_95813cf85 md5=30d01765970051212176cc2eb81ea034 ARGS: -r /resnick/groups/carnegie_poc/jingze/b1_1124/c4rst/bsg/rst_in/bsg3d.00002.rst rad_vet/time2_gtime_slope=true rad_vet/clean_stag_thinf=true rad_vet/clean_stag_lat_sums=true rad_vet/clean_stag_tie_gid=true rad_vet/clean_stag_opac_drag=true rad_vet/vet_gd_allow_clamp=false hydro/vceil=3.0e8 rad_vet/clean_pc_sweep=rb rad_vet/clean_stag_mg_levels=2 time/cfl_number=0.3 rad_vet/clean_conv_max=1e-5 output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30 time/nlim=7303
- C1: 2026-11-24 18:07:45 rc=0 (srun) ended: Terminating on cycle limit wall=61 s CW 0 FATAL 0 inf 0 DT_RETRY 0 stage fallback 0 clamped 0
- user.hst (compare on the printed t):
  - C1: 3 rows, 3 with the same t as the c1 run; first difference from the c1 run: none (all printed values identical)
  - C4: 3 rows, 3 with the same t as the c1 run; first difference from the c1 run: none (all printed values identical)
  - C4 vs C1 restart: all printed values identical
- hydro.hst (compare on the printed t):
  - C1: 3 rows, 3 with the same t as the c1 run; first difference from the c1 run: none (all printed values identical)
  - C4: 3 rows, 3 with the same t as the c1 run; first difference from the c1 run: none (all printed values identical)
  - C4 vs C1 restart: all printed values identical
VERDICT IDENTICAL
- Flag set: a restart of bsg that happens anyway uses c4 (launch_c4r.sh through launch_c1.sh's hook).

## m35mms: PASS (11-24 18:16), m35mms.00001.rst, cycle 7576 + 3, on 01-17
- C4: 2026-11-24 18:14:31 start m35mms C4 hpc-sm-01-17 athena_c4_a000d803c md5=3ae130e97b13e5bdfa4a7f9a14370b9f ARGS: -r /resnick/groups/carnegie_poc/jingze/b1_1124/c4rst/m35mms/rst_in/m35mms.00001.rst rad_vet/time2_gtime_slope=true rad_vet/clean_stag_thinf=true rad_vet/clean_stag_lat_sums=true rad_vet/clean_stag_tie_gid=true rad_vet/clean_stag_opac_drag=true rad_vet/vet_gd_allow_clamp=false hydro/vceil=3.0e8 rad_vet/clean_pc_sweep=rb rad_vet/clean_stag_mg_levels=4 time/cfl_number=0.3 rad_vet/clean_conv_max=1e-5 output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30 time/nlim=7579
- C4: 2026-11-24 18:15:47 rc=0 (srun) ended: Terminating on cycle limit wall=76 s CW 0 FATAL 0 inf 0 DT_RETRY 0 stage fallback 0 clamped 0
- C1: 2026-11-24 18:15:48 start m35mms C1 hpc-sm-01-17 athena_sdorder3_95813cf85 md5=30d01765970051212176cc2eb81ea034 ARGS: -r /resnick/groups/carnegie_poc/jingze/b1_1124/c4rst/m35mms/rst_in/m35mms.00001.rst rad_vet/time2_gtime_slope=true rad_vet/clean_stag_thinf=true rad_vet/clean_stag_lat_sums=true rad_vet/clean_stag_tie_gid=true rad_vet/clean_stag_opac_drag=true rad_vet/vet_gd_allow_clamp=false hydro/vceil=3.0e8 rad_vet/clean_pc_sweep=rb rad_vet/clean_stag_mg_levels=4 time/cfl_number=0.3 rad_vet/clean_conv_max=1e-5 output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30 time/nlim=7579
- C1: 2026-11-24 18:16:54 rc=0 (srun) ended: Terminating on cycle limit wall=66 s CW 0 FATAL 0 inf 0 DT_RETRY 0 stage fallback 0 clamped 0
- user.hst (compare on the printed t):
  - C1: 3 rows, 3 with the same t as the c1 run; first difference from the c1 run: none (all printed values identical)
  - C4: 3 rows, 3 with the same t as the c1 run; first difference from the c1 run: none (all printed values identical)
  - C4 vs C1 restart: all printed values identical
- hydro.hst (compare on the printed t):
  - C1: 3 rows, 3 with the same t as the c1 run; first difference from the c1 run: none (all printed values identical)
  - C4: 3 rows, 3 with the same t as the c1 run; first difference from the c1 run: none (all printed values identical)
  - C4 vs C1 restart: all printed values identical
VERDICT IDENTICAL
- Flag set: a restart of m35mms that happens anyway uses c4 (launch_c4r.sh through launch_c1.sh's hook).

## m35zams: PASS (11-24 18:42), m35zams.00001.rst, cycle 24824 + 3, on 01-17
- C4: 2026-11-24 18:41:19 start m35zams C4 hpc-sm-01-17 athena_c4_a000d803c md5=3ae130e97b13e5bdfa4a7f9a14370b9f ARGS: -r /resnick/groups/carnegie_poc/jingze/b1_1124/c4rst/m35zams/rst_in/m35zams.00001.rst rad_vet/time2_gtime_slope=true rad_vet/clean_stag_thinf=true rad_vet/clean_stag_lat_sums=true rad_vet/clean_stag_tie_gid=true rad_vet/clean_stag_opac_drag=true rad_vet/vet_gd_allow_clamp=false hydro/vceil=3.0e8 rad_vet/clean_pc_sweep=rb rad_vet/clean_stag_mg_levels=2 time/cfl_number=0.3 rad_vet/clean_conv_max=1e-5 output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30 time/nlim=24827
- C4: 2026-11-24 18:41:52 rc=0 (srun) ended: Terminating on cycle limit wall=33 s CW 0 FATAL 0 inf 0 DT_RETRY 0 stage fallback 0 clamped 0
- C1: 2026-11-24 18:41:53 start m35zams C1 hpc-sm-01-17 athena_sdorder3_95813cf85 md5=30d01765970051212176cc2eb81ea034 ARGS: -r /resnick/groups/carnegie_poc/jingze/b1_1124/c4rst/m35zams/rst_in/m35zams.00001.rst rad_vet/time2_gtime_slope=true rad_vet/clean_stag_thinf=true rad_vet/clean_stag_lat_sums=true rad_vet/clean_stag_tie_gid=true rad_vet/clean_stag_opac_drag=true rad_vet/vet_gd_allow_clamp=false hydro/vceil=3.0e8 rad_vet/clean_pc_sweep=rb rad_vet/clean_stag_mg_levels=2 time/cfl_number=0.3 rad_vet/clean_conv_max=1e-5 output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30 time/nlim=24827
- C1: 2026-11-24 18:42:26 rc=0 (srun) ended: Terminating on cycle limit wall=33 s CW 0 FATAL 0 inf 0 DT_RETRY 0 stage fallback 0 clamped 0
- user.hst (compare on the printed t):
  - C1: 3 rows, 3 with the same t as the c1 run; first difference from the c1 run: none (all printed values identical)
  - C4: 3 rows, 3 with the same t as the c1 run; first difference from the c1 run: none (all printed values identical)
  - C4 vs C1 restart: all printed values identical
- hydro.hst (compare on the printed t):
  - C1: 3 rows, 3 with the same t as the c1 run; first difference from the c1 run: none (all printed values identical)
  - C4: 3 rows, 3 with the same t as the c1 run; first difference from the c1 run: none (all printed values identical)
  - C4 vs C1 restart: all printed values identical
VERDICT IDENTICAL
- Flag set: a restart of m35zams that happens anyway uses c4 (launch_c4r.sh through launch_c1.sh's hook).
