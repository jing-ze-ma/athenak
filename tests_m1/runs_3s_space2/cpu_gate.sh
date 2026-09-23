#!/bin/bash
# CPU gates of implicit_enthalpy (runs_3s_space2): neutrality (base vs new, switch off),
# plm/central on the 2-D He slab (200 s) and the 3-D box (60 cycles), 1 block / 2 blocks /
# 2 ranks.  Same setup as tests_m1/runs_3p_fastdefault.
W=/viper/ptmp2/jinma/space2_3s/gate; R=$W/cpu_run.sh
V="rad_m1/closure=vet_sc rad_m1/vet_tensor=full"
S="time/tlim=200"
P="rad_m1/implicit_enthalpy=plm"
B3="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=60"
$R base e_base 1 slab2d_def $S &
$R new  e_off  1 slab2d_def $S &
$R new  e_plm  1 slab2d_def $S $P &
$R new  e_cen  1 slab2d_def $S rad_m1/implicit_enthalpy=central &
$R new  e_plm_mb 1 slab2d_def $S $P meshblock/nx2=16 &
$R new  e_plm_m2 2 slab2d_def $S $P meshblock/nx2=16 &
$R new  e_plm_nost 1 slab2d_def $S $P rad_m1/implicit_op_stencil=false &
$R base v_base 1 slab2d_def $S $V &
$R new  v_off  1 slab2d_def $S $V &
$R new  v_plm  1 slab2d_def $S $V $P &
$R new  v_plm_m2 2 slab2d_def $S $V $P meshblock/nx2=16 &
$R base e3_base 1 box3d_def $B3 &
$R new  e3_off 1 box3d_def $B3 &
$R new  e3_plm 1 box3d_def $B3 $P &
$R new  e3_plm_m2 2 box3d_def $B3 $P &
$R new  v3_off 1 box3d_def $B3 $V &
$R new  v3_plm 1 box3d_def $B3 $V $P &
$R new  v3_plm_m2 2 box3d_def $B3 $V $P &
wait
echo CPU GATE DONE
