#!/bin/bash
# usage: rw.sh LIST NP BUILD RUNSDIR  (radwave cases on the login node, nice 10)
#   BUILD = base | new (the build_mpi of /viper/ptmp2/jinma/accel_0923/<BUILD>)
source /etc/profile.d/modules.sh; module purge; module load gcc/14 openmpi/5.0
export OMP_NUM_THREADS=1
export RW_EXE=/viper/ptmp2/jinma/accel_0923/$3/build_mpi/src/athena
export RW_RUNS=/viper/ptmp2/jinma/accel_0923/rw/$4
mkdir -p $RW_RUNS
cd /viper/ptmp2/jinma/wt_accel/tests_m1/runs_4a_accel
nice -n 10 python3 run_radwave.py run $1 --np $2
echo "RW DONE $1 $3 $4"
