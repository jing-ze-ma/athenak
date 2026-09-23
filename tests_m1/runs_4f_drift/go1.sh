#!/bin/bash
# build (incremental) then run jobs1.txt, 12 at a time
cd /viper/ptmp2/jinma/wt_drift/build_cpu && nice -n 10 make -j 12 > ../build_cpu.make2.log 2>&1 || exit 1
cd /viper/ptmp2/jinma/wt_drift/tests_m1/runs_4f_drift
xargs -P 12 -L 1 ./run_t7.sh < jobs1.txt
echo ALLDONE
