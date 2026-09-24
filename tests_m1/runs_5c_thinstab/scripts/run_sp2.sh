#!/bin/bash -l
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
W=/viper/ptmp2/jinma/thinstab_0924; I=$W/sp2/sph_atm_bin.athinput; X=$W/bin/athena_sp2ctr_cpu
N="mesh/nx1=32 meshblock/nx1=32 time/nlim=300 rad_m1/implicit_lin_tol=1.0e-10 rad_m1/implicit_tol=1.0e-10 rad_m1/implicit_maxit=100 output3/dcycle=25 problem/atm_seed=1.0e-6"
S="rad_m1/kappa_p=0.0 rad_m1/kappa_e=0.0 rad_m1/kappa_f=0.0 rad_m1/kappa_s=25.0"
C="mesh/use_spherical_polar=false mesh/x2min=0 mesh/x2max=0.2 output1/slice_x2=0.1 output2/slice_x2=0.1"
r() { d=$W/sp2/$1; shift; rm -rf $d; mkdir -p $d; cd $d
  OMP_NUM_THREADS=1 nice mpirun -np 1 --bind-to none $X -i $I -d $d "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt; }
r zc_base $N $S $C &
r zc_ctr $N $S $C rad_m1/implicit_closure_thin_relax=1.5 &
r sph_base $N $S &
r sph_ctr $N $S rad_m1/implicit_closure_thin_relax=1.5 &
wait
