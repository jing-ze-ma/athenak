#!/bin/bash -l
# usage: run.sh <name> <ref|new> <nranks> <input> [overrides...]
name=$1; which=$2; np=$3; inp=$4; shift 4
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
R=/viper/ptmp2/jinma/sctb_0923
exe=$R/bin/athena_${which}_cpu
rm -rf $R/cpu/$name; mkdir -p $R/cpu/$name && cd $R/cpu/$name
if [ "$np" = 0 ]; then
  nice -n 10 $exe -i $R/cpu/$inp "$@" > log.txt 2>&1
else
  nice -n 10 mpirun -np $np --oversubscribe --bind-to none $exe -i $R/cpu/$inp "$@" > log.txt 2>&1
fi
echo "exit $?" >> log.txt
