#!/bin/bash
# run.sh <name> [param=value ...]
# One serial 1-D M1 column run in tests_m1/runs_2a_diag/<name>.
cd /viper/u2/jinma/ATHENAK/athenak/tests_m1/runs_2a_diag
N=$1; shift
rm -rf $N; mkdir -p $N
../../build_cpu_box/src/athena -i ../../inputs/hydro/he_box_m1_1d.athinput -d $N \
  problem/m1_top_bc=dark output2/dt=1.0e30 output3/dt=1.0e30 "$@" > $N/run.log 2>&1
echo "$N done rc=$?"
