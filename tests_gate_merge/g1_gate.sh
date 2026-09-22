#!/bin/bash -l
G=/viper/u2/jinma/ATHENAK/athenak/tests_gate_merge
IN=/viper/u2/jinma/ATHENAK/bench/hestar_fecz/box_w8/he_box_w8.athinput
# wait for builds
until grep -q "BUILD_DONE" $G/build_box_new.log && grep -q "BUILD_DONE" $G/build_box_ref.log; do sleep 10; done
grep -h "MAKE_EXIT\|BUILD_DONE" $G/build_box_new.log $G/build_box_ref.log
cp /viper/u2/jinma/ATHENAK/athenak/build_cpu_box/src/athena $G/athena_box_new || exit 1
cp /viper/u2/jinma/ATHENAK/bench/wt_ref/build_cpu_box/src/athena $G/athena_box_ref || exit 1
for v in ref new; do
  for M in 3 0; do
    D=$G/g1_${v}_m${M}
    rm -rf $D; mkdir -p $D; cd $D
    OPT=""
    if [ "$M" = "0" ]; then OPT="problem/rt_implicit_column=0 problem/rt_col3_skip_sweep=false"; fi
    $G/athena_box_$v -i $IN mesh/nx2=16 mesh/nx3=16 meshblock/nx1=134 meshblock/nx2=8 meshblock/nx3=8 time/nlim=50 $OPT > run.log 2>&1
    echo "RUN $v m$M exit=$?"
  done
done
echo GATE_RUNS_DONE
