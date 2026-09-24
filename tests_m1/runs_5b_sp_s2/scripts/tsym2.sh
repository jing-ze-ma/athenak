#!/bin/bash
# T-sym with gas and restart, M1 closure, c = 100 (with c = 1 the gas reaches v ~ c and
# the O(v/c) moment system leaves its domain: f -> 1 in a tau = 2000 medium)
cd /viper/ptmp2/jinma/s2_0924
X=$PWD/bin/athena_new_none_cpu
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
M="rad_m1/closure=m1 rad_m1/implicit_precond=line rad_m1/c_light=100.0"
./run.sh symg_m1 new none 4 $PWD/inp/sph_sym_gas.athinput $M &
./run.sh rstA_m1 new none 4 $PWD/inp/sph_sym_gas_rst.athinput $M &
wait
d=$PWD/cpu/rstB_m1; rm -rf $d; mkdir -p $d; cd $d
OMP_NUM_THREADS=1 mpirun -np 4 --oversubscribe --bind-to none $X -r $d/../rstA_m1/rst/sd.00001.rst -d $d > log.txt 2>&1
echo "rc=$?" >> log.txt
