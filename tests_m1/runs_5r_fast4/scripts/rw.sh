#!/bin/bash
# radwave order set for MR k = 1, 2, 4 (binary $1, tag $2)
source /etc/profile.d/modules.sh; module purge; module load gcc/14 openmpi/5.0
export OMP_NUM_THREADS=1
export RW_EXE=$1
cd /viper/ptmp2/jinma/wt_fast4/tests_m1/runs_4a_accel
for k in ${KS:-1 2 4}; do
  export RW_RUNS=/viper/ptmp2/jinma/fast4_0925/rw/$2/k$k; mkdir -p $RW_RUNS
  nice -n 10 python3 run_radwave.py run /viper/ptmp2/jinma/fast4_0925/rw/lists/k$k.txt --np 16 > $RW_RUNS.log 2>&1
done
echo RW DONE
