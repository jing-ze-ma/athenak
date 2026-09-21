#!/bin/bash -l
R=/viper/u2/jinma/ATHENAK/athenak
X=$R/build_cpu_m1/src/athena
IN=$R/inputs/tests/rad_m1_jump.athinput
O=$R/tests_m1/runs_1c
# case A (tau_cell 1 -> 1e3): the thick side needs L^2/D ~ 5e4 to reach the DISCRETE
# steady state, so tlim = 1e5 with dumps at 5e4 and 1e5 to show it has stopped moving.
A="rad_m1/kappa_f=64.0 problem/jump_flux=1.0e-6 time/tlim=1.0e5 output1/dt=5.0e4 output2/dt=5.0e4"
run () { local n=$1; shift; local d=$O/t3bA_$n; rm -rf $d; mkdir -p $d
  (cd $d && $X -i $IN $A "$@" > run.log 2>&1; echo "exit=$?" >> run.log); }
run ap_unified rad_m1/ap_form=unified &
run ap_alpha2  rad_m1/ap_form=alpha2 &
run none       rad_m1/thick_flux=none &
run scaled     rad_m1/thick_flux=scaled &
wait
echo T3BA_DONE
