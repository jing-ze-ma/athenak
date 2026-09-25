#!/bin/bash
# CPU gates of m1-coarse2: gates_gc.py (default bitwise base vs new; mg_gc and
# rho_direct fixed point), rst_gc.py (restart), tst/test_suite/rad_m1 (from the snapshot)
W=/viper/ptmp2/jinma/coarse2_0925; T=$W/src_${1:-c2}/tests_m1; C=${1:-c2}
until [ -f $W/bin/athena_${C}_box_convection_cpu ]; do sleep 60; done
source /etc/profile.d/modules.sh; module purge; module load gcc/14 openmpi/5.0 cmake/4.0
cp /viper/ptmp2/jinma/wt_coarse2/tests_m1/runs_5p_coarse2/*.py $T/runs_5p_coarse2/ 2>/dev/null || { mkdir -p $T/runs_5p_coarse2; cp /viper/ptmp2/jinma/wt_coarse2/tests_m1/runs_5p_coarse2/*.py $T/runs_5p_coarse2/; }
python3 $T/runs_5p_coarse2/gates_gc.py run $W/bin/athena_base_box_convection_cpu $W/bin/athena_${C}_box_convection_cpu $W/gates > $W/gates_run.log 2>&1
python3 $T/runs_5p_coarse2/gates_gc.py eval $W/gates > $W/gates_eval.txt 2>&1
ATHENAK_M1_DATA=/viper/ptmp2/jinma/faces_0924/m1data python3 $T/runs_5p_coarse2/rst_gc.py $W/bin/athena_${C}_box_convection_cpu $W/rstgc > $W/rstgc.log 2>&1
cd $W/src_${C}/tst && ATHENAK_M1_DATA=/viper/ptmp2/jinma/faces_0924/m1data python3 -m pytest -q test_suite/rad_m1 > $W/tst_rad_m1.log 2>&1
python3 $W/scripts/refgates.py run > $W/refgates_run.log 2>&1; python3 $W/scripts/refgates.py eval > $W/refgates_eval.txt 2>&1
echo CPU_GATES_DONE
