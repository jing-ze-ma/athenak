#!/bin/bash -l
# T3b (design sect. 11 and 13): the opacity-jump steady state, default f_source = cell
# against the interface form f_source = wb, for the two jumps straddling tau_cell = 1.
#   case A: tau_cell 1    -> 1e3   (kappa_f = 64,    flux 1e-6, tlim 1e5)
#   case B: tau_cell 1e-3 -> 1     (kappa_f = 0.064, flux 1e-3, tlim 2000)
R=/viper/u2/jinma/ATHENAK/bench/wt_m1gates
X=$R/build/src/athena
IN=$R/inputs/tests/rad_m1_jump.athinput
O=$R/tests_m1/runs_open
A="rad_m1/kappa_f=64.0 problem/jump_flux=1.0e-6 time/tlim=1.0e5 output1/dt=5.0e4 output2/dt=5.0e4"
B="rad_m1/kappa_f=0.064 problem/jump_flux=1.0e-3 time/tlim=2000.0 output1/dt=1000.0 output2/dt=1000.0"
run () { local n=$1; shift; local d=$O/t3b_$n; rm -rf $d; mkdir -p $d
  (cd $d && $X -i $IN "$@" > run.log 2>&1; echo "exit=$?" >> run.log); }
run A_cell $A rad_m1/f_source=cell &
run A_wb   $A rad_m1/f_source=wb &
run B_cell $B rad_m1/f_source=cell &
run B_wb   $B rad_m1/f_source=wb &
wait
echo T3B_DONE
