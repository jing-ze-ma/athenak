#!/bin/bash -l
# run.sh <name> <exe> <np> <input> [overrides...]   (CPU, runs/<name>)
source /etc/profile.d/modules.sh 2>/dev/null; module purge >/dev/null 2>&1
module load gcc/14 openmpi/5.0 >/dev/null 2>&1
W=/viper/ptmp2/jinma/sphhalo_0924; n=$1 exe=$2 np=$3 inp=$4; shift 4
d=$W/runs/$n; rm -rf $d; mkdir -p $d; cd $d
OMP_NUM_THREADS=1 nice mpirun -np $np --oversubscribe --bind-to none $exe $([[ $inp == *.rst ]] && echo -r || echo -i) $inp -d $d "$@" > log.txt 2>&1
echo "rc=$?" >> log.txt
