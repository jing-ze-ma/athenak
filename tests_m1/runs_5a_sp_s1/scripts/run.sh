#!/bin/bash -l
# usage: run.sh <name> <ref|new> <problem> <nranks> <input-path> [overrides...]   (CPU)
name=$1; which=$2; prob=$3; np=$4; inp=$5; shift 5
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
W=/viper/ptmp2/jinma/s1_0924
exe=$W/bin/athena_${which}_${prob}_cpu
d=$W/cpu/$name; rm -rf $d; mkdir -p $d; cd $d
OMP_NUM_THREADS=1 nice mpirun -np $np --oversubscribe --bind-to none $exe \
  -i $inp -d $d "$@" > log.txt 2>&1
echo "rc=$?" >> log.txt
