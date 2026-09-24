#!/bin/bash
# restart gate with vet_col_surface_q = true (every-step rebuild), 4 ranks, restart at 250
cd /viper/ptmp2/jinma/vetcol2_0924
X=$PWD/bin/athena_new_none_cpu
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
./run.sh rstA_sq new none 4 $PWD/inp2/sph_sym_gas_rst.athinput rad_m1/closure=vet_col rad_m1/c_light=100.0
d=$PWD/cpu/rstB_sq; rm -rf $d; mkdir -p $d; cd $d
OMP_NUM_THREADS=1 mpirun -np 4 --oversubscribe --bind-to none $X -r $PWD/../rstA_sq/rst/sd.00001.rst -d $d > log.txt 2>&1
echo "rc=$?" >> log.txt
cd ..; python3 ../cmp.py rstB_sq rstA_sq | tail -n1
grep -h "Marshak q" rstA_sq/log.txt | head -1
