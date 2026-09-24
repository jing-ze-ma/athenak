#!/bin/bash -l
# risk 6 v2: hot uniform gas (T = 1, c_s ~ 1.3 << c = 100).  Stage A: held gas, the
# hydro step made irrelevant by the hold (cfl_number = 100, the static uniform gas has
# balanced fluxes), 1000 steps to radiative equilibrium.  Stage B: restart with the gas
# free, cfl_number = 0.3, 100 steps: wb_arad (face form), none, and the control.
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
W=/viper/ptmp2/jinma/s2_0924; X=$W/bin/athena_new_none_cpu
cd $W && ./run.sh r6A2 new none 1 $W/inp/sph_atm_r6.athinput problem/atm_temp=1.0 time/cfl_number=100.0 output4/dcycle=1000
R=$W/cpu/r6A2/rst/r6.00001.rst
run() { local d=$W/cpu/$1; shift; rm -rf $d; mkdir -p $d; cd $d
  OMP_NUM_THREADS=1 mpirun -np 1 $X -r $R -d $d time/nlim=1100 time/cfl_number=0.3 \
    problem/atm_hold=false output1/dcycle=100 output2/dcycle=100 output4/dcycle=10 \
    "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt; }
run r6B2_wb rad_m1/dbg_gas_force=true rad_m1/force_reference=wb_arad &
run r6B2_none rad_m1/dbg_gas_force=true rad_m1/force_reference=none &
run r6B2_ctl rad_m1/dbg_gas_force=false &
wait
