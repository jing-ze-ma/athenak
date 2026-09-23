#!/bin/bash
# Step 1 CPU gates: implicit_halo_overlap (binary ovl) vs new (m1-int merge) and vs off
W=/viper/ptmp2/jinma/m1int_0924; R=$W/scripts/run.sh
V="rad_m1/closure=vet_sc rad_m1/vet_tensor=full"
S="time/tlim=100"
X2="meshblock/nx2=16"; X4="meshblock/nx2=8"
B3="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=30"
HM="rad_m1/implicit_halo_mpi=true"; PP="rad_m1/implicit_krylov_pipe=true"; OV="rad_m1/implicit_halo_overlap=true"
# off: ovl binary, switch absent, vs the merged binary (runs of the merge gate: e_new2/4, ...)
echo "2	$R o_e2hm   ovl 2 slab2d_def $S $X2 $HM"
echo "2	$R n_e2hm   new 2 slab2d_def $S $X2 $HM"
echo "4	$R o_e4     ovl 4 slab2d_def $S $X4"
echo "4	$R o_e34hm  ovl 4 box3d_def $B3 $HM"
echo "4	$R n_e34hm  new 4 box3d_def $B3 $HM"
echo "4	$R o_pvh4hm ovl 4 slab2d_plm_vimp_h2 $S $X4 $HM"
# on vs off (HM)
echo "2	$R o_e2ov   ovl 2 slab2d_def $S $X2 $HM $OV"
echo "4	$R o_e4hm   ovl 4 slab2d_def $S $X4 $HM"
echo "4	$R o_e4ov   ovl 4 slab2d_def $S $X4 $HM $OV"
echo "2	$R o_v2hm   ovl 2 slab2d_def $S $X2 $V $HM"
echo "2	$R o_v2ov   ovl 2 slab2d_def $S $X2 $V $HM $OV"
echo "4	$R o_e34ov  ovl 4 box3d_def $B3 $HM $OV"
echo "4	$R o_v34hm  ovl 4 box3d_def $B3 $V $HM"
echo "4	$R o_v34ov  ovl 4 box3d_def $B3 $V $HM $OV"
echo "4	$R o_pvh4ov ovl 4 slab2d_plm_vimp_h2 $S $X4 $HM $OV"
echo "4	$R o_pv4hm  ovl 4 slab2d_plm_vimp $S $X4 $HM"
echo "4	$R o_pv4ov  ovl 4 slab2d_plm_vimp $S $X4 $HM $OV"
echo "4	$R o_b4hm   ovl 4 box3d_plm_vimp_h2 $B3 $HM"
echo "4	$R o_b4ov   ovl 4 box3d_plm_vimp_h2 $B3 $HM $OV"
# with the pipe (P1/P2 red = 0: bitwise except the start-up op)
echo "4	$R o_e4ph   ovl 4 slab2d_def $S $X4 $HM $PP"
echo "4	$R o_e4pho  ovl 4 slab2d_def $S $X4 $HM $PP $OV"
echo "4	$R o_b4ph   ovl 4 box3d_plm_vimp_h2 $B3 $HM $PP"
echo "4	$R o_b4pho  ovl 4 box3d_plm_vimp_h2 $B3 $HM $PP $OV"
# restart (overlap on, 2 ranks, plm+vimp+hesdirk2): to t=40 with rst at t=20
echo "2	$R o_rsfull ovl 2 slab2d_plm_vimp_h2 time/tlim=40 output3/dt=20 output2/dt=40 output5/dt=40 $X2 $HM $OV && $W/scripts/rst.sh o_rsfull o_rsrst 2"
