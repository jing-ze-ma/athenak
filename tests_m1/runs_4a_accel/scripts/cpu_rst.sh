#!/bin/bash -l
# usage: cpu_rst.sh <src rundir> <new rundir> [overrides]  restart from the t=100 file (rst 00001)
module purge; module load gcc/14 openmpi/5.0
W=/viper/ptmp2/jinma/accel_0923
s=$1; n=$2; shift 2
d=$W/cpu/$n; rm -rf $d; mkdir -p $d; cd $d
OMP_NUM_THREADS=1 nice -n 10 mpirun -np 1 $W/new/build_boxcpu/src/athena \
  -r $W/cpu/$s/rst/m1slab.00001.rst -d $d "$@" > log.txt 2>&1
echo "rc=$?" >> log.txt
