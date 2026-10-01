#!/bin/bash
# runv.sh <bin> <suffix> [extra keys] -- the VARIANTS.md matrix for one binary/key set:
# runall.sh's static cases + grad (T gradient, deep diffusion) + nw10 (Newton, dt 10 s)
B=$1; SUF=$2; shift 2; EX="$@"
cd /viper/ptmp2/jinma/cksph_test_0930; . keys.sh
U="mesh/use_grid_stretch_r_poly=false"
for T in 2000 3000; do
 ./run.sh $B w_T${T}_$SUF base_w121.athinput problem/rt_test_T=$T $SIMP $EX &
 ./run.sh $B o_T${T}_$SUF base_w121.athinput problem/rt_test_T=$T $OLD $SIMP $EX &
done
./run.sh $B grad_$SUF base_w121.athinput problem/rt_test_T=2500 problem/rt_test_dT=-0.3 problem/rt_test_nwave=0.5 $SIMP $EX &
./run.sh $B nw10_$SUF base_w121.athinput problem/rt_test_T=2000 problem/rt_test_dt=10 time/nlim=10 problem/rt_test_out=rt.txt $SIMP $EX &
./run.sh $B tr_x1.2_$SUF base_w121.athinput problem/rt_test_T=2000 problem/rt_test_rho=1e-14 problem/rt_test_rstep=1 problem/rt_test_rho_top=1e-14 $U mesh/x1min=1.78127e10 mesh/x1max=1.95129e10 $SIMP $EX &
for xs in "1.7 1.46886e10" "3 1.95129e10"; do set -- $xs
 ./run.sh $B tr_x$1_$SUF base_w121.athinput problem/rt_test_T=2000 problem/rt_test_rho=1e-14 problem/rt_test_rstep=1 problem/rt_test_rho_top=1e-14 $U mesh/x1max=$2 $SIMP $EX &
done
for xs in "1.2 1.78127e10" "1.7 1.49658e10"; do set -- $xs
 ./run.sh $B st_x$1_$SUF base_w121.athinput problem/rt_test_T=2000 problem/rt_test_rho=1e-1 problem/rt_test_rstep=$2 problem/rt_test_rho_top=1e-14 $U mesh/x1max=1.95129e10 $SIMP $EX &
done; wait
for r in w_T2000 w_T3000 o_T2000 o_T3000 tr_x1.2 tr_x1.7 tr_x3 st_x1.2 st_x1.7; do
 R=runs/${r}_$SUF; I0=60
 nice python3 ana/exact.py $R/dump.txt $(basename $R) $I0 > $R/exact.txt 2>&1 &
done; wait
