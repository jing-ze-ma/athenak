#!/bin/bash -l
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 cmake/4.0 >/dev/null 2>&1
export ATHENAK_M1_DATA=/viper/ptmp2/jinma/thinstab_0924/m1data
cd /viper/ptmp2/jinma/wt_thinstab/tst
for t in test_suite/rad_m1/test_rad_m1_slab_cpu.py test_suite/rad_m1/test_rad_m1_opcheck_mpicpu.py test_suite/rad_m1/test_rad_m1_restart_mpicpu.py; do
  echo "=== $t"; python3 run_test_suite.py --test $t 2>&1 | tail -15
done
