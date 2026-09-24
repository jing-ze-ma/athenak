#!/bin/bash
# runs_5c step 1: growth vs tau_cell, k, nu, closure, offdiag (binary s0 = HEAD + pgen seed)
cd /viper/ptmp2/jinma/thinstab_0924
I=$PWD/inp/uni.athinput
R="./run.sh"
for ks in 0.25 0.5 1 2 4 8; do $R tau_k$ks s0 $I rad_m1/kappa_s=$ks & done
for k in 2 4 0; do $R kmode_$k s0 $I problem/atm_seed_k=$k & done
wait
for c in 1 3 30 3000 30000; do $R nu_$c s0 $I rad_m1/implicit_cfl=$c & done
$R cl_ker s0 $I rad_m1/closure=kershaw &
$R cl_edd s0 $I rad_m1/closure=eddington &
$R od_none s0 $I rad_m1/implicit_offdiag=none &
$R od_lag s0 $I rad_m1/implicit_offdiag=lagged &
$R mem0 s0 $I rad_m1/dbg_trans_memory=0.0 &
$R ny4 s0 $I mesh/nx2=4 meshblock/nx2=4 mesh/x2max=0.5 &
wait
for c in pass_and pass_rel03 pass_rel01 tlim_lp tlim_lp05 pass; do :; done
$R c_pass s0 $I rad_m1/implicit_closure_lag=pass &
$R c_pass_rel01 s0 $I rad_m1/implicit_closure_lag=pass rad_m1/implicit_closure_relax=0.1 &
$R c_pass_rel03 s0 $I rad_m1/implicit_closure_lag=pass rad_m1/implicit_closure_relax=0.3 &
$R c_pass_and s0 $I rad_m1/implicit_closure_lag=pass rad_m1/implicit_accel=anderson &
$R c_pass_and10 s0 $I rad_m1/implicit_closure_lag=pass rad_m1/implicit_accel=anderson rad_m1/implicit_anderson_m=10 rad_m1/implicit_maxit=200 &
$R c_tlim s0 $I rad_m1/implicit_trans_limit=lp &
$R c_tlim05 s0 $I rad_m1/implicit_trans_limit=lp rad_m1/implicit_trans_fmax=0.5 &
$R c_tlim05_odn s0 $I rad_m1/implicit_trans_limit=lp rad_m1/implicit_trans_fmax=0.5 rad_m1/implicit_offdiag=none &
wait
