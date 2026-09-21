#!/bin/bash -l
# The runs of run_reg.sh that were lost to the inode quota, with output2 (the tab
# series) throttled to the same cadence as output1.
R=/viper/u2/jinma/ATHENAK/bench/wt_m1gates
X=$R/build/src/athena
O=$R/tests_m1/runs_open/reg
IT=$R/inputs/tests/rad_m1_thick_pulse.athinput
IP=$R/inputs/tests/rad_m1_advect_pulse.athinput
run () { local n=$1; local IN=$2; shift 2; local d=$O/$n; rm -rf $d; mkdir -p $d
  (cd $d && $X -i $IN "$@" > run.log 2>&1; echo "exit=$?" >> run.log); }
for F in cell wb; do
  S="rad_m1/f_source=$F"
  run ${F}_t3_tau1e6 $IT $S mesh/nghost=3 rad_m1/kappa_s=128000000 time/tlim=3000.0 \
      output1/dt=300.0 output2/dt=300.0 &
  run ${F}_tophat $IT $S rad_m1/kappa_s=1280000 problem/m1_test=tophat \
      time/tlim=100.0 output1/dt=10.0 output2/dt=10.0 &
  run ${F}_t4_s512 $IP $S time/tlim=4.8e-6 output1/dt=1.2e-6 output2/dt=1.2e-6 \
      mesh/nx1=512 meshblock/nx1=512 rad_m1/kappa_p=500 &
done
wait
echo REG2_DONE
