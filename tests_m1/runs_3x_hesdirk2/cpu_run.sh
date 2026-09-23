#!/bin/bash -l
# usage: cpu_run.sh <base|new> <rundir> <np> <inp> [overrides...]
module purge; module load gcc/14 openmpi/5.0
W=/viper/ptmp2/jinma/h2_3x
b=$1; n=$2; np=$3; inp=$4; shift 4
X=$W/$b/build_boxcpu/src/athena
d=$W/cpu/$n; rm -rf $d; mkdir -p $d; cd $d
OMP_NUM_THREADS=1 nice -n 10 mpirun -np $np --oversubscribe $X \
  -i $W/inp/$inp.athinput -d $d problem/vpert=1.0e-2 "$@" > log.txt 2>&1
echo "rc=$?" >> log.txt
