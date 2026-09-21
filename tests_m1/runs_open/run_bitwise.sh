#!/bin/bash -l
# BITWISE gate: the default path (f_source = cell) against the pristine rt-integration
# HEAD binary, on the same input files.  $1 = binary, $2 = output subdirectory.
R=/viper/u2/jinma/ATHENAK/bench/wt_m1gates
X=$1
O=$R/tests_m1/runs_open/bw/$2
IT=$R/inputs/tests/rad_m1_thick_pulse.athinput
IP=$R/inputs/tests/rad_m1_advect_pulse.athinput
IU=$R/inputs/tests/rad_m1_advect_uniform.athinput
IM=$R/inputs/tests/rad_m1_marshak.athinput
IJ=$R/inputs/tests/rad_m1_jump.athinput
rm -rf $O; mkdir -p $O
run () { local n=$1; local IN=$2; shift 2; local d=$O/$n; mkdir -p $d
  (cd $d && $X -i $IN "$@" > run.log 2>&1; echo "exit=$?" >> run.log); }
run t3_tau1e3 $IT mesh/nghost=3 rad_m1/kappa_s=128000 time/tlim=12.0 \
    output1/dt=1.2 output2/dt=1.2 &
run tophat $IT rad_m1/kappa_s=1280000 problem/m1_test=tophat time/tlim=100.0 \
    output1/dt=10.0 output2/dt=10.0 &
run t4_d512 $IP time/tlim=4.8e-6 output1/dt=1.2e-6 output2/dt=1.2e-6 \
    mesh/nx1=512 meshblock/nx1=512 rad_m1/kappa_p=500 problem/pulse_v=3.0e7 &
run t4b_full $IU &
run t6_n128 $IM &
run t3b_A $IJ rad_m1/kappa_f=64.0 problem/jump_flux=1.0e-6 time/tlim=2000.0 \
    output1/dt=1000.0 output2/dt=1000.0 &
wait
echo BW_DONE_$2
