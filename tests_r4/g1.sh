#!/bin/bash -l
# G1: the plane-parallel BOX regression.  box_w8 in the production mode-3 configuration,
# shrunk to nx2 = nx3 = 16 / meshblock 134x8x8, 50 cycles, serial CPU.
set -e
ROOT=/viper/u2/jinma/ATHENAK/bench/wt_he4
IN=/viper/u2/jinma/ATHENAK/bench/hestar_fecz/box_w8/he_box_w8.athinput
EXE=$1
DIR=$2
MODE=$3          # 3 (production) or 0 (explicit sweep)
cd $ROOT/tests_r4/$DIR
rm -rf bin cbin_hydro_w_2 rst rt_surface.bin rt_profile.bin *.hst *.log column_used.txt
OPT=""
if [ "$MODE" = "0" ]; then OPT="problem/rt_implicit_column=0 problem/rt_col3_skip_sweep=false"; fi
$EXE -i $IN mesh/nx2=16 mesh/nx3=16 meshblock/nx1=134 meshblock/nx2=8 meshblock/nx3=8 \
  time/nlim=50 $OPT > run.log 2>&1
echo "exit $?"
