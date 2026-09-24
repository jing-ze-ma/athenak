#!/bin/bash -l
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
W=/viper/ptmp2/jinma/sph2_0924; G=/viper/ptmp2/jinma/wt_sph2/tests_m1/gates
export OMP_NUM_THREADS=1
nice python3 $G/gates.py run $W/bin/athena_new_box_cpu $W/gates > $W/gates_run.log 2>&1
python3 $G/gates.py eval $W/gates > $W/gates_eval.txt 2>&1
echo GATESPY DONE
