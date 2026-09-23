#!/bin/bash
# usage: run.sh <dir> <rstfile> <nlim> <input>
module load gcc/14 openmpi/5.0
cd /viper/ptmp2/jinma/rstjolt_0923
mkdir -p $1
nice -n 10 mpirun -np 6 ./athena_cpu -r $2 -i $4 -d $1 time/nlim=$3 \
  output3/dt=1e30 output4/dt=1e30 > $1/run.log 2>&1
echo done $? >> $1/run.log
