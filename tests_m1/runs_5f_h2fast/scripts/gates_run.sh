#!/bin/bash -l
source /etc/profile.d/modules.sh; module purge; module load gcc/14 openmpi/5.0
cd /viper/ptmp2/jinma/wt_h2fast/tests_m1/gates
export OMP_NUM_THREADS=1
nice -n 10 python3 gates.py run /viper/ptmp2/jinma/h2fast_0924/bin/athena_v4_box /viper/ptmp2/jinma/h2fast_0924/cpu/gates > /viper/ptmp2/jinma/h2fast_0924/cpu/gates_run.log 2>&1
python3 gates.py eval /viper/ptmp2/jinma/h2fast_0924/cpu/gates
echo GATES_DONE
