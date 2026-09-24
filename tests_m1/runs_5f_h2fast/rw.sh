#!/bin/bash
# usage: rw.sh LIST NP BIN RUNSDIR  (radwave cases on the login node, nice 10)
#   BIN = the tag of /viper/ptmp2/jinma/h2fast_0924/bin/athena_<BIN>_mpi
source /etc/profile.d/modules.sh; module purge; module load gcc/14 openmpi/5.0
export OMP_NUM_THREADS=1
export RW_EXE=/viper/ptmp2/jinma/h2fast_0924/bin/athena_$3_mpi
export RW_RUNS=/viper/ptmp2/jinma/h2fast_0924/rw/$4
mkdir -p $RW_RUNS
cd /viper/ptmp2/jinma/wt_h2fast/tests_m1/runs_5f_h2fast
nice -n 10 python3 run_radwave.py run $1 --np $2
echo "RW DONE $1 $3 $4"
