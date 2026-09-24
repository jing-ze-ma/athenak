#!/bin/bash
cd /viper/ptmp2/jinma/vetcol_0924
for n in 16 32 64 128; do
  ./run.sh a4v$n new none 1 $PWD/inp/sph_atm_vc.athinput mesh/nx1=$n meshblock/nx1=$n &
  ./run.sh a4m$n ref none 1 $PWD/inp/sph_atm.athinput mesh/nx1=$n meshblock/nx1=$n &
  ./run.sh a4e$n ref none 1 $PWD/inp/sph_atm.athinput mesh/nx1=$n meshblock/nx1=$n rad_m1/closure=eddington &
done
for n in 32 64; do
  ./run.sh a4vs$n new none 1 $PWD/inp/sph_atm_str_vc.athinput mesh/nx1=$n meshblock/nx1=$n &
  ./run.sh a4v2_$n new none 1 $PWD/inp/sph_atm_vc.athinput mesh/nx1=$n meshblock/nx1=$n rad_m1/vet_col_nsub=2 &
done
wait
