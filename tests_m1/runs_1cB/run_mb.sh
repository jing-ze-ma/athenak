#!/bin/bash -l
# The second 1c-A gap: hydro ghost zones at MeshBlock boundaries when the module does
# not write hydro back.  T4 dynamic, 512 cells, 1 MeshBlock vs 2, with and without
# gas_feedback.  All four must agree to truncation level.
R=/viper/u2/jinma/ATHENAK/athenak
X=$R/build_cpu_m1/src/athena
IN=$R/inputs/tests/rad_m1_advect_pulse.athinput
O=$R/tests_m1/runs_1cB
TL="time/tlim=4.8e-6 output1/dt=4.8e-6 output2/dt=4.8e-6"
DYN="rad_m1/kappa_p=500 problem/pulse_v=3.0e7 mesh/nx1=512"
run () { local n=$1; shift; local d=$O/mb_$n; rm -rf $d; mkdir -p $d
  (cd $d && $X -i $IN $TL $DYN "$@" > run.log 2>&1; echo "exit=$?" >> run.log); }
run fb1_1blk meshblock/nx1=512 &
run fb1_2blk meshblock/nx1=256 &
run fb1_4blk meshblock/nx1=128 &
run fb0_1blk meshblock/nx1=512 rad_m1/gas_feedback=false &
run fb0_2blk meshblock/nx1=256 rad_m1/gas_feedback=false &
run fb0_4blk meshblock/nx1=128 rad_m1/gas_feedback=false &
run cp0_1blk meshblock/nx1=512 rad_m1/coupling=false &
run cp0_2blk meshblock/nx1=256 rad_m1/coupling=false &
wait
echo MB_DONE
