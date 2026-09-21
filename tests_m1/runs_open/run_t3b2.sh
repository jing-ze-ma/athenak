#!/bin/bash -l
# T3b with the half-cell ghost offset of the Dirichlet BC removed
# (<problem>/exact_ghost = true): what is left is the JUMP error alone.
R=/viper/u2/jinma/ATHENAK/bench/wt_m1gates
X=$R/build/src/athena
IN=$R/inputs/tests/rad_m1_jump.athinput
O=$R/tests_m1/runs_open
A="rad_m1/kappa_f=64.0 problem/jump_flux=1.0e-6 time/tlim=1.0e5 output1/dt=5.0e4 output2/dt=5.0e4"
B="rad_m1/kappa_f=0.064 problem/jump_flux=1.0e-3 time/tlim=2000.0 output1/dt=1000.0 output2/dt=1000.0"
G="problem/exact_ghost=true"
run () { local n=$1; shift; local d=$O/t3bx_$n; rm -rf $d; mkdir -p $d
  (cd $d && $X -i $IN $G "$@" > run.log 2>&1; echo "exit=$?" >> run.log); }
run A_cell $A rad_m1/f_source=cell &
run A_wb   $A rad_m1/f_source=wb &
run B_cell $B rad_m1/f_source=cell &
run B_wb   $B rad_m1/f_source=wb &
run A_cell_n128 $A rad_m1/f_source=cell mesh/nx1=128 meshblock/nx1=128 &
run A_wb_n128   $A rad_m1/f_source=wb   mesh/nx1=128 meshblock/nx1=128 &
run A_none $A rad_m1/thick_flux=none &
run A_scaled $A rad_m1/thick_flux=scaled &
wait
echo T3BX_DONE
