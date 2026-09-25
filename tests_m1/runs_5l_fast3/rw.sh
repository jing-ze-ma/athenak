#!/bin/bash
source /etc/profile.d/modules.sh; module purge; module load gcc/14 openmpi/5.0
export OMP_NUM_THREADS=1
export RW_EXE=/viper/ptmp2/jinma/fast3/bin/athena_n2_none_cpu
cd /viper/ptmp2/jinma/wt_fast3/tests_m1/runs_4a_accel
for lt in 1.0e-10 1.0e-9 1.0e-8; do
  export RW_RUNS=/viper/ptmp2/jinma/fast3/rw/r$lt; mkdir -p $RW_RUNS
  nice -n 10 python3 run_radwave.py run /viper/ptmp2/jinma/fast3/rw/lists/lt$lt.txt --np 16
done
echo RW DONE
