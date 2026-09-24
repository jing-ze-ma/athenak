#!/bin/bash -l
# the M1-closure radwave on the CARTESIAN slab (twin of T5 m1): is its hesdirk2 order sp-specific?
W=/viper/ptmp2/jinma/sph2_0924; I=$W/inp; X=$W/bin/athena_new_none_cpu
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
run() { local n=$1; shift; local d=$W/cpu/$n; rm -rf $d; mkdir -p $d
  (cd $d && OMP_NUM_THREADS=1 nice mpirun -np 1 $X -d $d "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt) & }
O="time/nlim=-1 output1/dt=0.592156525463792"
for c in 0.4 0.2 0.1 0.05; do for a in be h2; do run t5c_${a}_$c -i $I/rwc_$a.athinput $O time/cfl_number=$c; done; done
run t5c_ref -i $I/rwc_h2.athinput $O time/cfl_number=0.0125
wait
for a in be h2; do echo "T5C m1 cartesian $a: $(python3 $W/scripts/ana.py order 'tab/*.tab' dens $W/cpu/t5c_ref $W/cpu/t5c_${a}_{0.4,0.2,0.1,0.05} --pert 1.0)"; done
echo T5C DONE
