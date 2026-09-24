#!/bin/bash -l
source /etc/profile.d/modules.sh; module purge; module load gcc/14 openmpi/5.0 cmake/4.0
export ATHENAK_M1_DATA=/viper/ptmp2/jinma/faces_0924/m1data
cd /viper/ptmp2/jinma/m1def2_0924/src_new/tst
nice -n 10 python3 -m pytest -q test_suite/rad_m1 2>&1 | tail -15
echo PYTEST_DONE
