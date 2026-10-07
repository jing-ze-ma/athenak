#!/bin/bash
# usage: run_cpu.sh <name> <nranks> <override string...>   (login node, nice)
NAME=$1; NP=$2; shift 2
BIN=/viper/ptmp2/jinma/accretor_1006/build_cpu/src/athena
INP=/viper/ptmp2/jinma/wt_accretor/inputs/hydro/ry_per_accretor.athinput
D=/viper/ptmp2/jinma/accretor_1006/tests/$NAME
mkdir -p $D && cd $D || exit 1
source /etc/profile.d/modules.sh
module purge >/dev/null 2>&1
module load gcc/14 openmpi/5.0
md5sum $BIN > md5.txt
echo "$@" > keys.txt
# quarter resolution: 112 x 4 x 512 (dphi 0.01227), 8 blocks
LOW="mesh/nx1=112 mesh/nx3=512 mesh/x2min=1.5462527 mesh/x2max=1.5953400 meshblock/nx1=112 meshblock/nx3=64"
/usr/bin/time -v nice -n 10 mpirun -np $NP --oversubscribe $BIN -i $INP -d $D $LOW "$@" > out.log 2> err.log
echo "exit $?" >> out.log
