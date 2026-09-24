#!/bin/bash
cd /viper/ptmp2/jinma/s2_0924
for n in 16 32 64; do ./run.sh sa$n new none 1 $PWD/inp/sph_atm.athinput mesh/nx1=$n meshblock/nx1=$n & done
for n in 32 64; do ./run.sh sas$n new none 1 $PWD/inp/sph_atm_str.athinput mesh/nx1=$n meshblock/nx1=$n & done
./run.sh sae32 new none 1 $PWD/inp/sph_atm.athinput mesh/nx1=32 meshblock/nx1=32 rad_m1/closure=eddington &
./run.sh r6A new none 1 $PWD/inp/sph_atm_r6.athinput output4/dcycle=1000 &
wait
