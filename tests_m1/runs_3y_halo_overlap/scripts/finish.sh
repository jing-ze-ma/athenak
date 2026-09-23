#!/bin/bash
# tests_m1/runs_3y_halo_overlap: every table of the README from the run trees
W=/viper/ptmp2/jinma/m1int_0924
python3 $W/scripts/tsum.py $W/runs/s12 $W/runs/w2 $W/runs/f4 $W/runs/f8
python3 $W/scripts/profsum.py $W/runs/prof
for a in g2_hm g2_hmo g2_ph g2_pho; do
  echo "== $a"; python3 $W/scripts/halogap.py $W/runs/prof/${a}_n20/prof/rank_0_kernel_trace.csv 30 | grep -E "total idle|Halo"
done
