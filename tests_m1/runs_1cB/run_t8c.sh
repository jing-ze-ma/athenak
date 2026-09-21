#!/bin/bash -l
# T8c: conservation over 200 cycles of the T4 dynamic case.
R=/viper/u2/jinma/ATHENAK/athenak
X=$R/build_cpu_m1/src/athena
IN=$R/inputs/tests/rad_m1_advect_pulse.athinput
O=$R/tests_m1/runs_1cB
d=$O/t8c; rm -rf $d; mkdir -p $d
(cd $d && $X -i $IN time/nlim=200 time/tlim=1.0 output1/dt=3.1e-7 output2/dt=3.1e-7 \
   mesh/nx1=512 meshblock/nx1=512 rad_m1/kappa_p=500 problem/pulse_v=3.0e7 \
   > run.log 2>&1; echo "exit=$?" >> run.log)
echo T8C_DONE
