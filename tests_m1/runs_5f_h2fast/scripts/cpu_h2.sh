#!/bin/bash -l
# hesdirk2 CPU gates with the v3 defaults (no lever keys named): He slab 2-D plm+vimp
source /etc/profile.d/modules.sh; module purge; module load gcc/14 openmpi/5.0
W=/viper/ptmp2/jinma/h2fast_0924; E=$W/bin/athena_v4_box; I=$W/inp/slab2d_plm_vimp_h2.athinput
V="rad_m1/closure=vet_sc rad_m1/vet_tensor=full"
P="problem/vpert=1.0e-2 output3/dt=100.0"
run() { # name np args...
  local n=$1 np=$2; shift 2; local d=$W/cpu/$n; rm -rf $d; mkdir -p $d
  (cd $d && OMP_NUM_THREADS=1 nice -n 10 mpirun -np $np --oversubscribe --bind-to none $E "$@" -d $d > log.txt 2>&1; echo "rc=$?" >> log.txt) }
run Zh2E 1 -i $I $P time/tlim=200 &
run Zh2V 2 -i $I $P time/tlim=200 meshblock/nx2=16 $V &
run Zh2E1000 1 -i $I $P time/tlim=1000 &
run Zh2V1000 2 -i $I $P time/tlim=1000 meshblock/nx2=16 $V &
run Zh2Efail 1 -i $I $P time/tlim=200 rad_m1/time2_dbg_fail=50 &
wait
run RZh2E 1 -r $W/cpu/Zh2E/rst/m1slab.00001.rst time/tlim=200 &
run RZh2V 2 -r $W/cpu/Zh2V/rst/m1slab.00001.rst time/tlim=200 &
wait
echo CPU_H2_DONE
