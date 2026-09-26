#!/bin/bash
# gate (A) resolution study: nx1 = 48, 96, 192 on a 4x4 lateral wedge, 1 rank, t = 3
source /etc/profile.d/modules.sh; module purge; module load gcc/14 openmpi/5.0
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/sprhd_0926; B=$W/bin/athena_$1_none_cpu; O=$W/cpu/Ares_$1
for n in 48 96 192; do
  mkdir -p $O/n$n; (cd $O/n$n && timeout 5000 mpirun -np 1 $B -i $W/inp/gateA.athinput -d . time/nlim=-1 time/tlim=3.0 mesh/nx1=$n meshblock/nx1=$n mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 > log.txt 2>&1; echo rc=$? >> log.txt) &
done
wait; echo ARES DONE
