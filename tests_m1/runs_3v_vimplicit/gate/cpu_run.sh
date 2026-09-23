#!/bin/bash -l
# usage: cpu_run.sh <base|new> <rundir> <np> <inp> [overrides...]   (base = 89dc28d7 build)
module purge; module load gcc/14 openmpi/5.0
W=/viper/ptmp2/jinma/vimp_3v
b=$1; n=$2; np=$3; inp=$4; shift 4
if [ "$b" == "base" ]; then X=/viper/ptmp2/jinma/space2_3s/new/build_boxcpu/src/athena; else X=$W/new/build_boxcpu/src/athena; fi
d=$W/gate/cpu/$n; rm -rf $d; mkdir -p $d; cd $d
OMP_NUM_THREADS=1 nice -n 10 mpirun -np $np --oversubscribe $X \
  -i $W/gate/inp/$inp.athinput -d $d problem/vpert=1.0e-2 "$@" > log.txt 2>&1
echo "rc=$?" >> log.txt
