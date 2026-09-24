#!/bin/bash
cd /viper/ptmp2/jinma/thinstab_0924
I=$PWD/inp/uni.athinput
R="./run.sh"
P="rad_m1/implicit_closure_lag=pass"
$R c_pass_and10 s0 $I $P rad_m1/implicit_accel=anderson rad_m1/implicit_anderson_m=10 rad_m1/implicit_maxit=200 &
$R c_tlim05 s0 $I rad_m1/implicit_trans_limit=lp rad_m1/implicit_trans_fmax=0.5 &
$R c_tlim05_odn s0 $I rad_m1/implicit_trans_limit=lp rad_m1/implicit_trans_fmax=0.5 rad_m1/implicit_offdiag=none &
for ks in 0.25 0.5 2; do
  $R and_k$ks s0 $I $P rad_m1/implicit_accel=anderson rad_m1/kappa_s=$ks &
  $R rel03_k$ks s0 $I $P rad_m1/implicit_closure_relax=0.3 rad_m1/kappa_s=$ks &
  $R base_pass_k$ks s0 $I rad_m1/kappa_s=$ks $P &
done
$R and_nyq s0 $I $P rad_m1/implicit_accel=anderson problem/atm_seed_k=0 &
$R and_ker s0 $I $P rad_m1/implicit_accel=anderson rad_m1/closure=kershaw &
$R step_k2 s0 $I rad_m1/kappa_s=2 &
wait
