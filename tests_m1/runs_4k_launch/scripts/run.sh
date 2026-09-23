#!/bin/bash -l
# usage: run.sh <name> <base|new> <nranks> <input: cpu/<x>.athinput name> [overrides...]
name=$1; which=$2; np=$3; inp=$4; shift 4
source /etc/profile.d/modules.sh 2>/dev/null; module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
W=/viper/ptmp2/jinma/launch_0923
exe=$W/bin/athena_${which}_boxcpu
d=$W/cpu/$name; rm -rf $d; mkdir -p $d; cd $d
OMP_NUM_THREADS=1 nice mpirun -np $np --oversubscribe --bind-to none $exe \
  -i $W/cpu/$inp.athinput -d $d problem/vpert=1.0e-2 "$@" > log.txt 2>&1
echo "rc=$?" >> log.txt
