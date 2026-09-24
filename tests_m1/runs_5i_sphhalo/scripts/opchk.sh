#!/bin/bash
# the implicit_op_check matrix on the sp wedge: closures x split geometries x theta bc
W=/viper/ptmp2/jinma/sphhalo_0924; X=$W/bin/athena_new_none_cpu; R=$W/scripts/run.sh
declare -A G
G[t2]="2 meshblock/nx2=6 meshblock/nx3=12"
G[p2]="2 meshblock/nx2=12 meshblock/nx3=6"
G[q2]="2 meshblock/nx2=6 meshblock/nx3=6"
G[q4]="4 meshblock/nx2=6 meshblock/nx3=6"
G[r4]="4 meshblock/nx2=6 meshblock/nx3=6 mesh/ix2_bc=reflect mesh/ox2_bc=reflect"
G[t4]="4 meshblock/nx2=3 meshblock/nx3=12"
for c in edd m1 vc; do
  for g in t2 p2 q2 q4 r4 t4; do
    set -- ${G[$g]}; np=$1; shift
    $R chk_${c}_$g $X $np $W/inp/chk_$c.athinput "$@" rad_m1/implicit_op_check=-2 &
  done
  wait
done
