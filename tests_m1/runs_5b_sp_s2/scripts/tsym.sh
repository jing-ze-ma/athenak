#!/bin/bash
# T-sym and restart with the M1 closure on the wedge (S1 inputs, closure = m1)
cd /viper/ptmp2/jinma/s2_0924
X=$PWD/bin/athena_new2_none_cpu
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
M="rad_m1/closure=m1 rad_m1/implicit_precond=line"
./run.sh sym_m1d new2 none 4 $PWD/inp/sph_sym.athinput $M rad_m1/implicit_lin_tol=1.0e-14 rad_m1/implicit_tol=1.0e-12 &
./run.sh sym_m1a new2 none 4 $PWD/inp/sph_sym.athinput $M &
./run.sh symg_m1 new2 none 4 $PWD/inp/sph_sym_gas.athinput $M &
./run.sh rstA_m1 new2 none 4 $PWD/inp/sph_sym_gas_rst.athinput $M &
wait
d=$PWD/cpu/rstB_m1; rm -rf $d; mkdir -p $d; cd $d
OMP_NUM_THREADS=1 mpirun -np 4 --oversubscribe --bind-to none $X -r $PWD/../rstA_m1/rst/sd.00001.rst -d $d > log.txt 2>&1
echo "rc=$?" >> log.txt
