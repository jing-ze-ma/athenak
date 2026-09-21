#!/bin/bash -l
# T3b redefined (design sect. 10): steady state across a 1e3 opacity jump, for two
# jump positions straddling tau_cell = 1.
#   case A: tau_cell 1    -> 1e3   (kappa_f = 64,    flux 1e-6)
#   case B: tau_cell 1e-3 -> 1     (kappa_f = 0.064, flux 1e-3)
R=/viper/u2/jinma/ATHENAK/athenak
X=$R/build_cpu_m1/src/athena
IN=$R/inputs/tests/rad_m1_jump.athinput
O=$R/tests_m1/runs_1c

run () {
  local n=$1; shift
  local d=$O/t3b_$n
  rm -rf $d; mkdir -p $d
  (cd $d && $X -i $IN "$@" > run.log 2>&1; echo "exit=$?" >> run.log)
}

A="rad_m1/kappa_f=64.0 problem/jump_flux=1.0e-6 time/tlim=2000.0 output1/dt=1000.0 output2/dt=1000.0"
B="rad_m1/kappa_f=0.064 problem/jump_flux=1.0e-3 time/tlim=2000.0 output1/dt=1000.0 output2/dt=1000.0"

for F in unified alpha2; do
  run A_ap_$F   $A rad_m1/ap_form=$F
  run B_ap_$F   $B rad_m1/ap_form=$F
done
run A_none   $A rad_m1/thick_flux=none
run A_scaled $A rad_m1/thick_flux=scaled
run B_none   $B rad_m1/thick_flux=none
run B_scaled $B rad_m1/thick_flux=scaled
echo T3B_DONE
