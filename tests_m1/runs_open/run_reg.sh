#!/bin/bash -l
# Milestone 1d regression: the DEFAULT path (f_source = cell) must reproduce the
# milestone-1c/1c-B dumps BITWISE, and the same set is re-run with f_source = wb so
# that the two can be tabulated side by side.
R=/viper/u2/jinma/ATHENAK/bench/wt_m1gates
X=$R/build/src/athena
O=$R/tests_m1/runs_open/reg
IT=$R/inputs/tests/rad_m1_thick_pulse.athinput
IP=$R/inputs/tests/rad_m1_advect_pulse.athinput
IU=$R/inputs/tests/rad_m1_advect_uniform.athinput
IM=$R/inputs/tests/rad_m1_marshak.athinput
rm -rf $O; mkdir -p $O
run () { local n=$1; local IN=$2; shift 2; local d=$O/$n; mkdir -p $d
  (cd $d && $X -i $IN "$@" > run.log 2>&1; echo "exit=$?" >> run.log); }
TL="time/tlim=4.8e-6 output1/dt=1.2e-6 output2/dt=1.2e-6"
N512="mesh/nx1=512 meshblock/nx1=512"
DYN="rad_m1/kappa_p=500 problem/pulse_v=3.0e7"
for F in cell wb; do
  S="rad_m1/f_source=$F"
  run ${F}_t3_tau1e1  $IT $S mesh/nghost=3 rad_m1/kappa_s=1280 time/tlim=0.5 \
      output1/dt=0.05 &
  run ${F}_t3_tau1e3  $IT $S mesh/nghost=3 rad_m1/kappa_s=128000 time/tlim=12.0 \
      output1/dt=1.2 &
  run ${F}_t3_tau1e6  $IT $S mesh/nghost=3 rad_m1/kappa_s=128000000 time/tlim=3000.0 \
      output1/dt=300.0 &
  run ${F}_t3_nyq     $IT $S mesh/nghost=3 rad_m1/kappa_s=128000 \
      problem/nyquist_amp=1.0e-3 time/tlim=20.0 output1/dt=2.0 &
  run ${F}_tophat     $IT $S rad_m1/kappa_s=1280000 problem/m1_test=tophat \
      time/tlim=100.0 output1/dt=10.0 &
  run ${F}_t4_s512    $IP $S $TL $N512 rad_m1/kappa_p=500 &
  run ${F}_t4_d512    $IP $S $TL $N512 $DYN &
  run ${F}_t4_d512ns  $IP $S $TL $N512 $DYN rad_m1/advect_split=false &
  run ${F}_t4b_full   $IU $S &
  run ${F}_t4b_ovc    $IU $S rad_m1/source_form=ovc &
  run ${F}_t6_n64     $IM $S mesh/nx1=64 meshblock/nx1=64 &
  run ${F}_t6_n128    $IM $S &
done
wait
echo REG_DONE
