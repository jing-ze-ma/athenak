#!/bin/bash -l
# Milestone 1c-B, STEP 0: reproduce the 1c-A default-scheme numbers on the post-merge
# binary.  T3 (tau_cell 1e3, plm), T4 dynamic 512 with and without the split, T4b.
R=/viper/u2/jinma/ATHENAK/athenak
X=$R/build_cpu_m1/src/athena
O=$R/tests_m1/runs_1cB/step0
rm -rf $O; mkdir -p $O
run () { local n=$1; local IN=$2; shift 2; local d=$O/$n; mkdir -p $d
  (cd $d && $X -i $IN "$@" > run.log 2>&1; echo "exit=$?" >> run.log); }

IT=$R/inputs/tests/rad_m1_thick_pulse.athinput
IP=$R/inputs/tests/rad_m1_advect_pulse.athinput
IU=$R/inputs/tests/rad_m1_advect_uniform.athinput
TL="time/tlim=4.8e-6 output1/dt=1.2e-6 output2/dt=1.2e-6"
N512="mesh/nx1=512 meshblock/nx1=512"
DYN="rad_m1/kappa_p=500 problem/pulse_v=3.0e7"

run t3_tau1e3_plm $IT rad_m1/kappa_s=128000 rad_m1/reconstruct=plm \
    time/tlim=12.0 output1/dt=1.2 &
run t4_s512       $IP $TL $N512 rad_m1/kappa_p=500 &
run t4_d512_split $IP $TL $N512 $DYN &
run t4_d512_nosplit $IP $TL $N512 $DYN rad_m1/advect_split=false &
run t4b_full      $IU &
run t4b_ovc       $IU rad_m1/source_form=ovc &
wait
echo STEP0_RUNS_DONE
