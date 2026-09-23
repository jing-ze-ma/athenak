#!/bin/bash
# usage: rw.sh LIST NP  (radwave cases with the new MPI build, login node nice 10)
source /etc/profile.d/modules.sh; module purge; module load gcc/14 openmpi/5.0
export OMP_NUM_THREADS=1
export RW_EXE=${RW_EXE:-/viper/ptmp2/jinma/h2_3x/new/build_mpi/src/athena}
cd /viper/ptmp2/jinma/wt_time2b/tests_m1/runs_3x_hesdirk2
nice -n 10 python3 run_radwave.py run $1 --np $2
echo "RW DONE $1"
