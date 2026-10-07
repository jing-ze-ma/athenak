#!/bin/bash
# stage-1 input, quarter resolution, 60 cycles, old (1e8d0ff8) vs new CPU binary; diff hst
source /etc/profile.d/modules.sh
module purge >/dev/null 2>&1
module load gcc/14 openmpi/5.0
INP=/viper/ptmp2/jinma/accretor_1006/bin/ry_per_accretor_1e8d0ff8.athinput
LOW="mesh/nx1=112 mesh/nx3=512 mesh/x2min=1.5462527 mesh/x2max=1.5953400 meshblock/nx1=112 meshblock/nx3=64 time/nlim=60 output1/dt=0.0001 output2/dt=100"
for arm in old new; do
  if [ $arm = old ]; then B=/viper/ptmp2/jinma/accretor_1006/bin/athena_ryper_cpu_1e8d0ff8; else B=/viper/ptmp2/jinma/accretor_1006/build_cpu/src/athena; fi
  D=/viper/ptmp2/jinma/accretor_1006/tests2/bw_$arm
  mkdir -p $D && cd $D || exit 1
  nice -n 10 mpirun -np 8 --oversubscribe $B -i $INP -d $D $LOW > out.log 2>&1
done
cmp /viper/ptmp2/jinma/accretor_1006/tests2/bw_old/ryper.user.hst /viper/ptmp2/jinma/accretor_1006/tests2/bw_new/ryper.user.hst && echo USER_HST_IDENTICAL
cmp /viper/ptmp2/jinma/accretor_1006/tests2/bw_old/ryper.hydro.hst /viper/ptmp2/jinma/accretor_1006/tests2/bw_new/ryper.hydro.hst && echo HYDRO_HST_IDENTICAL
