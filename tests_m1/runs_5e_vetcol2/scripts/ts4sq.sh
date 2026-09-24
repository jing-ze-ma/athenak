#!/bin/bash
cd /viper/ptmp2/jinma/vetcol2_0924
I=$PWD/inp2/sph_atm_vc.athinput; M=$PWD/inp2/milne_vc.athinput
SQ="rad_m1/vet_col_surface_q=true"
for n in 32 64 128 256; do
  L=800; [ $n = 256 ] && L=1600
  N="mesh/nx1=$n meshblock/nx1=$n time/nlim=$L"
  ./run.sh s4q$n new none 1 $I $N $SQ &
  ./run.sh s4qt$n new none 1 $I $N $SQ rad_m1/vet_col_team=true &
  ./run.sh mlo$n new none 1 $M $N &
  ./run.sh mlq$n new none 1 $M $N $SQ &
done
./run.sh s4t64 new none 1 $I mesh/nx1=64 meshblock/nx1=64 rad_m1/vet_col_team=true &
./run.sh s4n64 new none 1 $I mesh/nx1=64 meshblock/nx1=64 &
./run.sh s4o64 ref none 1 $I mesh/nx1=64 meshblock/nx1=64 &
wait
