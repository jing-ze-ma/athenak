#!/bin/bash -l
# risk 6, stage B: restart the static radiative-equilibrium state (r6A, cycle 400) with
# the gas free to move, 100 steps
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
W=/viper/ptmp2/jinma/s2_0924; X=$W/bin/athena_new_none_cpu; R=$W/cpu/r6A/rst/r6.00001.rst
run() { local d=$W/cpu/$1; shift; rm -rf $d; mkdir -p $d; cd $d
  OMP_NUM_THREADS=1 mpirun -np 1 $X -r $R -d $d time/nlim=1100 problem/atm_hold=false output1/dcycle=100 output2/dcycle=100 \
    "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt; }
run r6B_wb rad_m1/dbg_gas_force=true rad_m1/force_reference=wb_arad &
run r6B_wb0 rad_m1/dbg_gas_force=true rad_m1/force_reference=wb_arad hydro/sp_wellbalanced_src=false &
run r6B_none rad_m1/dbg_gas_force=true rad_m1/force_reference=none &
run r6B_ctl rad_m1/dbg_gas_force=false &
wait
