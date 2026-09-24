#!/bin/bash
cd /viper/ptmp2/jinma/vetcol_0924
X=$PWD/bin/athena_new_none_cpu
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
M="rad_m1/closure=vet_col"
./run.sh sym_vcd new none 4 $PWD/inp/sph_sym.athinput $M rad_m1/implicit_lin_tol=1.0e-14 rad_m1/implicit_tol=1.0e-12 &
./run.sh sym_vca new none 4 $PWD/inp/sph_sym.athinput $M &
./run.sh symg_vc new none 4 $PWD/inp/sph_sym_gas.athinput $M rad_m1/c_light=100.0 &
./run.sh rstA_vc new none 4 $PWD/inp/sph_sym_gas_rst.athinput $M rad_m1/c_light=100.0 &
./run.sh fast64 new none 1 $PWD/inp/sph_atm_vc.athinput mesh/nx1=64 meshblock/nx1=64 time/nlim=200 &
./run.sh slow64 new none 1 $PWD/inp/sph_atm_vc_slow.athinput mesh/nx1=64 meshblock/nx1=64 time/nlim=200 &
wait
d=$PWD/cpu/rstB_vc; rm -rf $d; mkdir -p $d; cd $d
OMP_NUM_THREADS=1 mpirun -np 4 --oversubscribe --bind-to none $X -r $PWD/../rstA_vc/rst/sd.00001.rst -d $d > log.txt 2>&1
echo "rc=$?" >> log.txt
