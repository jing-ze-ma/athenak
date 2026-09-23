#!/bin/bash -l
# usage: rundhj.sh <name> <ref|new> <hyd|mhd> <nranks> [overrides...]
# 3 cycles of the production input from the newest-at-setup production restart (read in
# place): hyd = bench/cs_hyd4_prod rst 00113 (t=1.723252e7, rot 56.50, cycle 873580),
# mhd = bench/cs_mhd_prod4 rst 00081 (t=1.235251e7, rot 40.50, cycle 821798).
name=$1; which=$2; sys=$3; np=$4; shift 4
source /etc/profile.d/modules.sh 2>/dev/null; module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
W=/viper/ptmp2/jinma/m1int_0924; B=/viper/u2/jinma/ATHENAK/bench
if [ $sys = hyd ]; then P=$B/cs_hyd4_prod; R=$P/rst/dhj.00113.rst; N=873583
else P=$B/cs_mhd_prod4; R=$P/rst/dhj.00081.rst; N=821801; fi
exe=$W/bin/athena_${which}_dhjcpu
d=$W/cpu/$name; rm -rf $d; mkdir -p $d; cd $d
OMP_NUM_THREADS=1 nice mpirun -np $np --oversubscribe --bind-to none $exe \
  -r $R -i $P/deep_hot_jupiter.athinput -d $d time/nlim=$N output1/dt=1.0 \
  output2/dt=1.0 output3/dt=1.0e20 output4/dt=1.0e20 "$@" > log.txt 2>&1
echo "rc=$?" >> log.txt
