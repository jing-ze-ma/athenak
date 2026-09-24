#!/bin/bash
cd /viper/ptmp2/jinma/s2_0924
I=$PWD/inp/sph_atm_bin.athinput
N="mesh/nx1=32 meshblock/nx1=32 time/nlim=50 rad_m1/implicit_lin_tol=1.0e-10 rad_m1/implicit_tol=1.0e-10 rad_m1/implicit_maxit=100 output3/dcycle=10 problem/atm_seed=1.0e-6"
S="rad_m1/kappa_p=0.0 rad_m1/kappa_e=0.0 rad_m1/kappa_f=0.0 rad_m1/kappa_s=25.0"
C="mesh/use_spherical_polar=false mesh/x2min=0 mesh/x2max=0.2 output1/slice_x2=0.1 output2/slice_x2=0.1"
./run.sh zc_mem0 new3 none 1 $I $N $S $C rad_m1/dbg_trans_memory=0.0 &
./run.sh zc_pass new3 none 1 $I $N $S $C rad_m1/implicit_closure_lag=pass &
./run.sh zc_rel new3 none 1 $I $N $S $C rad_m1/implicit_closure_lag=pass rad_m1/implicit_closure_relax=0.3 &
./run.sh zc_ker new3 none 1 $I $N $S $C rad_m1/closure=kershaw &
./run.sh zc_base new3 none 1 $I $N $S $C &
wait
