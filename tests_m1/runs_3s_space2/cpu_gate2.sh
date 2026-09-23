#!/bin/bash
# the implicit_enthalpy arms of cpu_gate.sh, with the switch written in the input file
# (a command-line override cannot add a key the file does not have)
W=/viper/ptmp2/jinma/space2_3s/gate; R=$W/cpu_run.sh
V="rad_m1/closure=vet_sc rad_m1/vet_tensor=full"
S="time/tlim=200"
B3="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=60"
$R new  e_plm  1 slab2d_plm $S &
$R new  e_cen  1 slab2d_central $S &
$R new  e_plm_mb 1 slab2d_plm $S meshblock/nx2=16 &
$R new  e_plm_m2 2 slab2d_plm $S meshblock/nx2=16 &
$R new  e_plm_nost 1 slab2d_plmnost $S &
$R new  v_plm  1 slab2d_plm $S $V &
$R new  v_plm_m2 2 slab2d_plm $S $V meshblock/nx2=16 &
$R new  e3_plm 1 box3d_plm $B3 &
$R new  e3_plm_m2 2 box3d_plm $B3 &
$R new  v3_plm 1 box3d_plm $B3 $V &
$R new  v3_plm_m2 2 box3d_plm $B3 $V &
wait
echo CPU GATE2 DONE
