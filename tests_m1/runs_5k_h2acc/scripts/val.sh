#!/bin/bash -l
# runs_5h relaxation tests (T-S4 atmosphere from the Eddington start, t = 6.4).  usage: val.sh TV|TV2|T1
W=/viper/ptmp2/jinma/h2acc_0924; I=/viper/ptmp2/jinma/sph2_0924/inp; A=$W/scripts/ana.py
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
run() { local n=$1 x=$2 np=$3; shift 3; local d=$W/cpu/val/$n; rm -rf $d; mkdir -p $d
  (cd $d && OMP_NUM_THREADS=1 nice mpirun -np $np --oversubscribe --bind-to none $W/bin/athena_${x}_none_cpu -d $d "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt) & }
G="rad_m1/time2_dbg_g3=0.15"
for g in "$@"; do case $g in
TV)  # vet_col with absorption (gas held): base (vet_gas old) vs new (start); m1 closure exp
  O="time/nlim=-1 time/tlim=6.4"
  for c in 0.3 0.15 0.075 0.0375; do
    run tv_base_$c base 1 -i $I/atmvc_h2.athinput $O time/cfl_number=$c
    run tv_new_$c new 1 -i $I/atmvc_h2.athinput $O time/cfl_number=$c
    run tv_exp_$c exp 1 -i $I/atmvc_h2.athinput $O time/cfl_number=$c $G
    run tv_be_$c base 1 -i $I/atmvc_be.athinput $O time/cfl_number=$c
    run tvm_new_$c new 1 -i $I/atm_h2.athinput $O time/cfl_number=$c mesh/nx1=32 meshblock/nx1=32
    run tvm_exp_$c exp 1 -i $I/atm_h2.athinput $O time/cfl_number=$c mesh/nx1=32 meshblock/nx1=32 $G
  done
  run tv_base_ref base 1 -i $I/atmvc_h2.athinput $O time/cfl_number=0.009375
  run tv_new_ref new 1 -i $I/atmvc_h2.athinput $O time/cfl_number=0.009375
  run tvm_new_ref new 1 -i $I/atm_h2.athinput $O time/cfl_number=0.009375 mesh/nx1=32 meshblock/nx1=32
  wait
  for a in base new; do echo "TV vet_col h2 $a (own ref): $(python3 $A order 'tab/sa.m1.*.tab' m1_e $W/cpu/val/tv_${a}_ref $W/cpu/val/tv_${a}_{0.3,0.15,0.075,0.0375})"; done
  echo "TV vet_col h2 exp (new ref): $(python3 $A order 'tab/sa.m1.*.tab' m1_e $W/cpu/val/tv_new_ref $W/cpu/val/tv_exp_{0.3,0.15,0.075,0.0375})"
  echo "TV vet_col be (new ref): $(python3 $A order 'tab/sa.m1.*.tab' m1_e $W/cpu/val/tv_new_ref $W/cpu/val/tv_be_{0.3,0.15,0.075,0.0375})"
  echo "TV vet_col base vs new refs: $(python3 $A order 'tab/sa.m1.*.tab' m1_e $W/cpu/val/tv_new_ref $W/cpu/val/tv_base_ref)"
  for a in new exp; do echo "TV m1 h2 $a: $(python3 $A order 'tab/sa.m1.*.tab' m1_e $W/cpu/val/tvm_new_ref $W/cpu/val/tvm_${a}_{0.3,0.15,0.075,0.0375})"; done;;
TV2)  # pure scattering relaxation: the stiff thin-top order, new vs exp tableau
  O="time/nlim=-1 time/tlim=6.4 rad_m1/kappa_p=0.0 rad_m1/kappa_e=0.0 rad_m1/kappa_f=0.0 rad_m1/kappa_s=100.0 rad_m1/implicit_tol=1.0e-12 rad_m1/implicit_lin_tol=1.0e-13"
  for c in 0.3 0.15 0.075 0.0375; do
    run tw_new_$c new 1 -i $I/atmvc_h2.athinput $O time/cfl_number=$c
    run tw_base_$c base 1 -i $I/atmvc_h2.athinput $O time/cfl_number=$c
    run tw_exp_$c exp 1 -i $I/atmvc_h2.athinput $O time/cfl_number=$c $G
  done
  run tw_ref new 1 -i $I/atmvc_h2.athinput $O time/cfl_number=0.009375
  wait
  for a in base new exp; do echo "TV2 vet_col $a: $(python3 $A order 'tab/sa.m1.*.tab' m1_e $W/cpu/val/tw_ref $W/cpu/val/tw_${a}_{0.3,0.15,0.075,0.0375})"; done;;
T1)  # transient spherical diffusion (T-S1 shell, stretched r, tau 2000), t = 100, vet_col
  O="time/nlim=-1 time/tlim=100.0 rad_m1/closure=vet_col"
  for c in 100.0 50.0 25.0 12.5; do
    run t1_new_$c new 4 -i $I/sym_h2.athinput $O rad_m1/implicit_cfl=$c
    run t1_exp_$c exp 4 -i $I/sym_h2.athinput $O rad_m1/implicit_cfl=$c $G
  done
  run t1_ref new 4 -i $I/sym_h2.athinput $O rad_m1/implicit_cfl=3.125
  wait
  for a in new exp; do echo "T1 vet_col $a: $(python3 $A order 'tab/sd.j0k0.*.tab' m1_e $W/cpu/val/t1_ref $W/cpu/val/t1_${a}_{100.0,50.0,25.0,12.5} --pert 0.0)"; done;;
esac; done
echo "VAL $* DONE"
