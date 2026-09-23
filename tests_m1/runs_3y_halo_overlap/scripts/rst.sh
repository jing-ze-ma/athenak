#!/bin/bash -l
# usage: rst.sh <fullrun> <rstrun> <np>: restart the ovl binary from <fullrun>'s first rst file
source /etc/profile.d/modules.sh 2>/dev/null; module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
W=/viper/ptmp2/jinma/m1int_0924
d=$W/cpu/$2; rm -rf $d; mkdir -p $d; cd $d
rf=$(ls $W/cpu/$1/rst/*.00001.rst)
OMP_NUM_THREADS=1 nice mpirun -np $3 --oversubscribe --bind-to none $W/bin/athena_ovl_boxcpu -r $rf -d $d > log.txt 2>&1
echo "rc=$?" >> log.txt
