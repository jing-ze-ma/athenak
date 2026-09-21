#!/bin/bash
# Gate (d): relax the same column with the semi-implicit apply and with ck_implicit, then
# dump T(p) for three columns out of each final state.  Same pattern as
# tests_ck_sph/rerun_dumps.sh: the correlated-k column dump is one-shot and one column per
# run, so the comparison is three short restarts per arm.
set -u
R=/viper/u2/jinma/ATHENAK/athenak
N=$R/build_ckimp_new/src/athena
IN=$R/inputs/tests/dhj_ck_implicit.athinput
B="output3/dt=1e30 problem/rt_use_cons=true problem/ck_spherical=true"
NCYC=${1:-600}

for arm in off on; do
  d=d_relax_$arm; rm -rf $d; mkdir -p $d
  EX=""
  [ $arm = on ] && EX="problem/ck_implicit=true problem/ck_impl_verbose=true"
  ( cd $d && $N -i $IN $B time/nlim=$NCYC output4/dt=1e30 $EX > run.log 2>&1 )
  echo "$arm exit $? $(tail -2 $d/run.log | head -1)"
done

for arm in off on; do
  RST=$(ls -1 d_relax_$arm/rst/*.rst | tail -1)
  for col in "night 0 2" "day92 1 2" "day38 0 4"; do
    set -- $col
    D=d_col_${1}_$arm; rm -rf $D; mkdir -p $D
    ( cd $D && $N -r ../$RST -t 00:00:20 time/nlim=999999 \
        problem/ck_dump_file=col.txt problem/ck_dump_m=$2 problem/ck_dump_k=$3 \
        output1/dt=1e30 output3/dt=1e30 output4/dt=1e30 > run.log 2>&1 )
    rm -rf $D/bin $D/rst
    echo "$D $(wc -l < $D/col.txt 2>/dev/null) lines"
  done
done
echo "relax done"
