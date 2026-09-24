#!/bin/bash
cd /viper/ptmp2/jinma/vetcol_0924
Q="rad_m1/marshak_q=0.59740391"
for n in 32 64 128 256; do
  L=800; [ $n = 256 ] && L=1600
  ./run.sh q4v$n new none 1 $PWD/inp/sph_atm_vc.athinput mesh/nx1=$n meshblock/nx1=$n $Q time/nlim=$L &
  ./run.sh q4m$n ref none 1 $PWD/inp/sph_atm.athinput mesh/nx1=$n meshblock/nx1=$n $Q time/nlim=$L &
  ./run.sh q4e$n ref none 1 $PWD/inp/sph_atm.athinput mesh/nx1=$n meshblock/nx1=$n rad_m1/closure=eddington $Q time/nlim=$L &
done
./run.sh q4v2_128 new none 1 $PWD/inp/sph_atm_vc.athinput mesh/nx1=128 meshblock/nx1=128 rad_m1/vet_col_nsub=2 $Q &
wait
