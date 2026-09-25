#!/bin/bash
source /etc/profile.d/modules.sh; module purge; module load gcc/14 openmpi/5.0
I=/viper/ptmp2/jinma/vetcol2_0924/inp/sph_atm_vc.athinput
E=/viper/ptmp2/jinma/fast3/bin/athena_n2_none_cpu
for lt in 1.0e-10 1.0e-9 1.0e-8; do
  d=/viper/ptmp2/jinma/fast3/ts4/l$lt; rm -rf $d; mkdir -p $d
  (cd $d && OMP_NUM_THREADS=1 nice mpirun -np 1 --oversubscribe --bind-to none $E -i $I -d $d \
    mesh/nx1=64 meshblock/nx1=64 time/nlim=800 \
    rad_m1/implicit_lin_tol=$lt > log.txt 2>&1) &
done
wait
cd /viper/ptmp2/jinma/wt_fast3/tests_m1/runs_5d_vetcol
python3 ts4v.py /viper/ptmp2/jinma/fast3/ts4 l1.0e-10 l1.0e-9 l1.0e-8
for lt in 1.0e-10 1.0e-9 1.0e-8; do grep -h "inner iterations mean\|NON-CONV" /viper/ptmp2/jinma/fast3/ts4/l$lt/log.txt | cut -c1-160; done
