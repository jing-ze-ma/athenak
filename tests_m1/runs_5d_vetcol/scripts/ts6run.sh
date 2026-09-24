#!/bin/bash
# T-S6: one step, dump of the first vet_col build (the analytic initial state)
cd /viper/ptmp2/jinma/vetcol_0924
I=$PWD/inp/sph_atm_vc.athinput
B="time/nlim=1 problem/atm_seed=0.0 rad_m1/vet_col_dump_every=1 output3/dcycle=100000"
FS="rad_m1/kappa_p=1.0e-12 rad_m1/kappa_e=1.0e-12 rad_m1/kappa_f=1.0e-12 rad_m1/kappa_s=0.0"
EX="rad_m1/kappa_p=0.0 rad_m1/kappa_e=0.0 rad_m1/kappa_f=10.0 rad_m1/kappa_s=0.0"
TK="rad_m1/kappa_p=0.0 rad_m1/kappa_e=0.0 rad_m1/kappa_f=1.0e4 rad_m1/kappa_s=0.0"
run() { # name case nr nc np
  n=$1; shift; c=$1; shift; nr=$1; nc=$2; np=$3
  eval K=\$$c
  ./run.sh $n new none 1 $I $B $K mesh/nx1=$nr meshblock/nx1=$nr rad_m1/vet_col_ncore=$nc \
    rad_m1/vet_col_nsub=$np rad_m1/vet_col_dump=$PWD/cpu/$n/dump
}
for c in ${CASES:-FS EX TK}; do
  for nr in 32 64 128 256; do run t6_${c}_r${nr} $c $nr 8 1 & done
  for nc in 2 4 16 32; do run t6_${c}_c${nc} $c 128 $nc 1 & done
  for np in 2 4; do run t6_${c}_p${np} $c 128 8 $np & done
  wait
done
