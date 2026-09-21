#!/bin/bash -l
R=/viper/u2/jinma/ATHENAK/bench/wt_m1gates
X=$R/build/src/athena
O=$R/tests_m1/runs_open/reg
IP=$R/inputs/tests/rad_m1_advect_pulse.athinput
TL="time/tlim=4.8e-6 output1/dt=1.2e-6 output2/dt=1.2e-6 mesh/nx1=512 meshblock/nx1=512"
DYN="rad_m1/kappa_p=500 problem/pulse_v=3.0e7"
run () { local n=$1; shift; local d=$O/$n; rm -rf $d; mkdir -p $d
  (cd $d && $X -i $IP $TL "$@" > run.log 2>&1; echo "exit=$?" >> run.log); }
for F in cell wb; do
  run ${F}_t4_d512   rad_m1/f_source=$F $DYN &
  run ${F}_t4_d512ns rad_m1/f_source=$F $DYN rad_m1/advect_split=false &
done
wait
echo REG3_DONE
