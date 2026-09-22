#!/bin/bash -l
G=/viper/u2/jinma/ATHENAK/athenak/tests_gate_merge
IN=/viper/u2/jinma/ATHENAK/bench/hestar_fecz/box_w8/he_box_w8.athinput
EX=/viper/u2/jinma/ATHENAK/athenak/build_cpu_box/src/athena
for M in 3 0; do
  D=$G/g1_m1_2adiag_m${M}
  rm -rf $D; mkdir -p $D; cd $D
  OPT=""
  if [ "$M" = "0" ]; then OPT="problem/rt_implicit_column=0 problem/rt_col3_skip_sweep=false"; fi
  $EX -i $IN mesh/nx2=16 mesh/nx3=16 meshblock/nx1=134 meshblock/nx2=8 meshblock/nx3=8 time/nlim=50 $OPT > run.log 2>&1
  echo "RUN m$M exit=$?"
done
echo G0_DONE
