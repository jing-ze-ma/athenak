#!/bin/bash -l
# usage: run.sh <name> <binary tag> <input> [overrides...]   (CPU, serial)
name=$1; tag=$2; inp=$3; shift 3
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
W=/viper/ptmp2/jinma/thinstab_0924
d=$W/runs/$name; rm -rf $d; mkdir -p $d; cd $d
OMP_NUM_THREADS=1 nice mpirun -np 1 --bind-to none $W/bin/athena_${tag}_cpu \
  -i $inp -d $d "$@" > log.txt 2>&1
echo "rc=$?" >> log.txt
