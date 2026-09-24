#!/bin/bash -l
# runs_5c: the 09-22 seeded 2-D He slab (runs_3b7 A0 s_base input), 50 s, serial CPU
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
W=/viper/ptmp2/jinma/thinstab_0924; I=$W/he/he_slab_m1_2d.athinput; X=$W/bin/athena_s1bc_cpu
A="problem/vpert=1.0e-2 time/tlim=50.0 output2/dt=1.0e30 output5/dt=1.0e30 output3/dt=1.0e30"
r() { d=$W/he/$1; shift; rm -rf $d; mkdir -p $d; cd $d
  OMP_NUM_THREADS=1 nice $X -i $I -d $d $A "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt; }
r H0_pass &
r H1_step rad_m1/implicit_closure_lag=step &
r H2_step_ctr rad_m1/implicit_closure_lag=step rad_m1/implicit_closure_thin_relax=1.5 &
r H3_pass_and rad_m1/implicit_accel=anderson &
wait
