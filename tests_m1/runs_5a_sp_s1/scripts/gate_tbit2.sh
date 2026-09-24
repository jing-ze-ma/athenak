#!/bin/bash -l
W=/viper/ptmp2/jinma/s1_0924; R=$W/run.sh
for c in ref new; do
  $R t_vimp2_$c $c box_convection 2 /viper/ptmp2/jinma/s0_0923/inp/slab2d_vimp.athinput time/tlim=50 problem/vpert=1.0e-2 meshblock/nx2=16 &
  $R t_spblast_$c $c sp_test 2 $W/inp/sp_blast.athinput time/nlim=100 &
done
wait
for t in vimp2 spblast; do
  echo "== $t: $(python3 $W/cmp.py $W/cpu/t_${t}_ref $W/cpu/t_${t}_new | tail -n1) rc $(tail -n1 $W/cpu/t_${t}_ref/log.txt) $(tail -n1 $W/cpu/t_${t}_new/log.txt)"
done
