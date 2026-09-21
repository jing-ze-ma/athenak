#!/bin/bash -l
# T2 (design sect. 8): Hayes & Norman / HERACLES shadow test, 280 x 80, 10 light
# crossing times, M1 against the closure = eddington control (which must fill the
# shadow in).  The ARRIVAL time needs a fine output cadence but only just past one
# crossing, so it comes from a separate short run; the depth comes from the long one.
R=/viper/u2/jinma/ATHENAK/bench/wt_m1gates
X=$R/build/src/athena
IN=$R/inputs/tests/rad_m1_shadow.athinput
O=$R/tests_m1/runs_open
TC=3.3356e-11
run () { local n=$1; shift; local d=$O/t2_$n; rm -rf $d; mkdir -p $d
  (cd $d && $X -i $IN "$@" > run.log 2>&1; echo "exit=$?" >> run.log); }
run m1        rad_m1/closure=m1        output1/dt=1.66780e-11 output2/dt=3.3356e-10 &
run eddington rad_m1/closure=eddington output1/dt=1.66780e-11 output2/dt=3.3356e-10 &
run arr_m1    rad_m1/closure=m1        time/tlim=4.0027e-11 \
              output1/dt=1.33424e-12 output2/dt=4.0027e-11 &
wait
echo T2_DONE
