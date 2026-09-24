#!/bin/bash
# CPU share of the vet_col build in the M1 step (correctness-machine number only; the cost
# comparison proper is the apudev GPU job): He wedge grid at 96 x 32 x 32, one block, 1 rank
cd /viper/ptmp2/jinma/vetcol_0924
H=$PWD/inp/hewedge.athinput
G="mesh/nx2=32 mesh/nx3=32 time/nlim=20"
for r in 1 2; do
  ./run.sh cc_edd_$r new none 1 $H $G rad_m1/closure=eddington
  ./run.sh cc_vc_$r new none 1 $H $G rad_m1/closure=vet_col
done
for t in cc_edd_1 cc_vc_1 cc_edd_2 cc_vc_2; do
  echo "== $t: $(grep -h 'cpu time used\|vet_col: .*builds' cpu/$t/log.txt | tr '\n' ' ')"
done
