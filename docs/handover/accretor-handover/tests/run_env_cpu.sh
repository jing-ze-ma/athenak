#!/bin/bash
# usage: [RES=q|h] run_env_cpu.sh <name> <nranks> <override string...>   (login node, nice)
# stage-2 envelope input at quarter (141 x 4 x 512, R_acc = face 29) or half resolution
# (282 x 4 x 1024, R_acc = face 58), 8 blocks
NAME=$1; NP=$2; shift 2
BIN=/viper/ptmp2/jinma/accretor_1006/build_cpu/src/athena
INP=/viper/ptmp2/jinma/wt_accretor/inputs/hydro/ry_per_accretor_env.athinput
D=/viper/ptmp2/jinma/accretor_1006/tests2/$NAME
mkdir -p $D && cd $D || exit 1
source /etc/profile.d/modules.sh
module purge >/dev/null 2>&1
module load gcc/14 openmpi/5.0
md5sum $BIN > md5.txt
cp $INP input.athinput
if [ "${RES:-q}" = h ]; then
LOW="mesh/nx1=282 mesh/nx3=1024 mesh/x2min=1.5585244804918964 mesh/x2max=1.5830681731030002 meshblock/nx1=282 meshblock/nx3=128"
else
LOW="mesh/nx1=141 mesh/nx3=512 mesh/x2min=1.5462526341887264 mesh/x2max=1.5953400194010667 meshblock/nx1=141 meshblock/nx3=64"
fi
echo "$@" > keys.txt
echo "RES=${RES:-q} $LOW" >> keys.txt
/usr/bin/time -v nice -n 10 mpirun -np $NP --oversubscribe $BIN -i $INP -d $D $LOW "$@" > out.log 2> err.log
echo "exit $?" >> out.log
