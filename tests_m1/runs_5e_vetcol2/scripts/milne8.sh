#!/bin/bash
cd /viper/ptmp2/jinma/vetcol2_0924
M=$PWD/inp2/milne_vc.athinput
for n in 32 64 128 256; do
  L=800; [ $n = 256 ] && L=1600
  N="mesh/nx1=$n meshblock/nx1=$n time/nlim=$L rad_m1/vet_col_nmu=8"
  ./run.sh ml8o$n new none 1 $M $N &
  ./run.sh ml8q$n new none 1 $M $N rad_m1/vet_col_surface_q=true &
done
wait
