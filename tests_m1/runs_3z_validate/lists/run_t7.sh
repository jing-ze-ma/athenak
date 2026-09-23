#!/bin/bash
# usage: run_t7.sh <name> <input> [overrides]   (one run, in t7/<name>)
V=/viper/ptmp2/jinma/validate_0923; X=$V/wt/build_mpi/src/athena
n=$1; inp=$2; shift 2
d=$V/t7/$n; rm -rf $d; mkdir -p $d; cp $V/t7/ref_*.txt $d/
cd $d && OMP_NUM_THREADS=1 nice -n 10 $X -i $inp "$@" > run.log 2>&1; echo "exit=$?" >> run.log
