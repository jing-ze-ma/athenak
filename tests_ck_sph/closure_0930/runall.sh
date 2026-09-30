#!/bin/bash
# runall.sh <bin> <suffix> [extra keys]  -- the full matrix for one binary/key set
B=$1; SUF=$2; shift 2; EX="$@"
cd /viper/ptmp2/jinma/cksph_test_0930; . keys.sh
U="mesh/use_grid_stretch_r_poly=false"
for T in 2000 3000; do
 ./run.sh $B w_T${T}_$SUF base_w121.athinput problem/rt_test_T=$T $SIMP $EX &
 ./run.sh $B o_T${T}_$SUF base_w121.athinput problem/rt_test_T=$T $OLD $SIMP $EX &
done; wait
./run.sh $B tr_x1.2_$SUF base_w121.athinput problem/rt_test_T=2000 problem/rt_test_rho=1e-14 problem/rt_test_rstep=1 problem/rt_test_rho_top=1e-14 $U mesh/x1min=1.78127e10 mesh/x1max=1.95129e10 $SIMP $EX &
for xs in "1.7 1.46886e10" "3 1.95129e10"; do set -- $xs
 ./run.sh $B tr_x$1_$SUF base_w121.athinput problem/rt_test_T=2000 problem/rt_test_rho=1e-14 problem/rt_test_rstep=1 problem/rt_test_rho_top=1e-14 $U mesh/x1max=$2 $SIMP $EX &
done
for xs in "1.2 1.78127e10" "1.7 1.49658e10"; do set -- $xs
 ./run.sh $B st_x$1_$SUF base_w121.athinput problem/rt_test_T=2000 problem/rt_test_rho=1e-1 problem/rt_test_rstep=$2 problem/rt_test_rho_top=1e-14 $U mesh/x1max=1.95129e10 $SIMP $EX &
done; wait
for r in runs/*_$SUF; do nice python3 ana/exact.py $r/dump.txt $(basename $r) 60 > $r/exact.txt 2>&1; done
