#!/bin/bash
# the three pRT emission runs, sequential (login node, nice); no core dumps
ulimit -c 0
source /viper/ptmp2/jinma/prt_venv/bin/activate
export OMP_NUM_THREADS=1
cd /viper/ptmp2/jinma/w121prod_0929/synth2_emis
for a in "w1x 1x_eq" "w10x 10x_eq" "w1x 1x_nodiss"; do
  set -- $a
  nice -n 10 python -u emis_prt.py $1 $2 12 > logs/prt_$1_$2.log 2>&1
done
echo ALLDONE >> logs/prt_w1x_1x_nodiss.log
