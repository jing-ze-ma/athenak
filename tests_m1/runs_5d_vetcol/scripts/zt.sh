#!/bin/bash
cd /viper/ptmp2/jinma/vetcol_0924
I=$PWD/inp/sph_atm_bin_vc.athinput
N="mesh/nx1=32 meshblock/nx1=32 time/nlim=50 rad_m1/implicit_lin_tol=1.0e-10 rad_m1/implicit_tol=1.0e-10 rad_m1/implicit_maxit=100 output3/dcycle=10 problem/atm_seed=1.0e-6"
N2="mesh/nx1=32 meshblock/nx1=32 time/nlim=400 rad_m1/implicit_lin_tol=1.0e-10 rad_m1/implicit_tol=1.0e-10 rad_m1/implicit_maxit=100 output3/dcycle=50 problem/atm_seed=1.0e-6"
S="rad_m1/kappa_p=0.0 rad_m1/kappa_e=0.0 rad_m1/kappa_f=0.0 rad_m1/kappa_s=25.0"
C="mesh/use_spherical_polar=false mesh/x2min=0 mesh/x2max=0.2 output1/slice_x2=0.1 output2/slice_x2=0.1"
./run.sh zc_m1 new none 1 $I $N $S $C rad_m1/closure=m1 &
./run.sh zc_edd new none 1 $I $N $S $C rad_m1/closure=eddington &
./run.sh zc_vc new none 1 $I $N $S $C rad_m1/closure=vet_col &
./run.sh zc_vcL new none 1 $I $N2 $S $C rad_m1/closure=vet_col &
./run.sh zc_vcF new none 1 $I $N2 $S $C rad_m1/closure=vet_col rad_m1/vet_col_axis=flux &
./run.sh zs_m1 new none 1 $I $N $S rad_m1/closure=m1 &
./run.sh zs_vcL new none 1 $I $N2 $S rad_m1/closure=vet_col &
./run.sh zs_vcF new none 1 $I $N2 $S rad_m1/closure=vet_col rad_m1/vet_col_axis=flux &
wait
