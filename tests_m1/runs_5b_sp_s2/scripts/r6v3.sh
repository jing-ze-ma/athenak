#!/bin/bash -l
# risk 6 v3: stage A2 (hot gas, held, cfl 100, 1000 steps) -> A3 (held, cfl 0.3, 10 steps:
# brings dt down, since a restart keeps the last dt) -> B3 (gas free, 100 steps)
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
W=/viper/ptmp2/jinma/s2_0924; X=$W/bin/athena_new_none_cpu
d=$W/cpu/r6A3; rm -rf $d; mkdir -p $d; cd $d
OMP_NUM_THREADS=1 mpirun -np 1 $X -r $W/cpu/r6A2/rst/r6.00001.rst -d $d time/nlim=1010 \
  time/cfl_number=0.3 output3/dcycle=10 > log.txt 2>&1; echo "rc=$?" >> log.txt
R=$(ls $d/rst/r6.*.rst | tail -1)
run() { local d=$W/cpu/$1; shift; rm -rf $d; mkdir -p $d; cd $d
  OMP_NUM_THREADS=1 mpirun -np 1 $X -r $R -d $d time/nlim=1110 time/cfl_number=0.3 \
    problem/atm_hold=false output1/dcycle=50 output2/dcycle=50 output4/dcycle=10 \
    "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt; }
run r6B3_wb rad_m1/dbg_gas_force=true rad_m1/force_reference=wb_arad &
run r6B3_none rad_m1/dbg_gas_force=true rad_m1/force_reference=none &
run r6B3_ctl rad_m1/dbg_gas_force=false &
wait
