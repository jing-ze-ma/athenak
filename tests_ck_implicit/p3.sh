#!/bin/bash
# tests_ck_implicit phase 3: the FROZEN EXCHANGE OPERATOR and the WARM START.
set -u
R=/viper/u2/jinma/ATHENAK/athenak
REF=$R/build_ckfast_ref/src/athena
NEW=$R/build_ckfast/src/athena
IN=$R/inputs/tests/dhj_ck_implicit.athinput
BASE="output1/dt=1e30 output2/dt=1e30 output4/dt=1e30 problem/rt_use_cons=true"
SPH="problem/ck_spherical=true problem/ck_beam_sph=true"

run() {  # run <dir> <binary> <extra args...>
  local d=$1; shift
  local b=$1; shift
  rm -rf "$d"; mkdir -p "$d"
  ( cd "$d" && /usr/bin/time -f "%e" -o t.txt $b -i $IN $BASE "$@" > run.log 2>&1 )
  echo "$? $d $(cat $d/t.txt 2>/dev/null)"
}

echo "=== (a) ck_implicit off: new vs reference, 20 cycles, one dump per cycle"
for sp in false true; do
  run p3a_ref_sp$sp $REF time/nlim=20 output3/dt=1e-6 problem/ck_spherical=$sp \
      problem/ck_beam_sph=true
  run p3a_new_sp$sp $NEW time/nlim=20 output3/dt=1e-6 problem/ck_spherical=$sp \
      problem/ck_beam_sph=true
done

echo "=== (b,d) implicit: phase-2 path vs frozen operator vs frozen+warm, 20 cycles"
run p3b_p2   $NEW time/nlim=20 output3/dt=1e-6 $SPH problem/ck_implicit=true \
    problem/ck_impl_verbose=true
run p3b_frz  $NEW time/nlim=20 output3/dt=1e-6 $SPH problem/ck_implicit=true \
    problem/ck_impl_verbose=true problem/ck_impl_frozen_op=true
run p3b_frzn $NEW time/nlim=20 output3/dt=1e-6 $SPH problem/ck_implicit=true \
    problem/ck_impl_verbose=true problem/ck_impl_frozen_op=true \
    problem/ck_impl_frozen_cof=false
run p3b_warm $NEW time/nlim=20 output3/dt=1e-6 $SPH problem/ck_implicit=true \
    problem/ck_impl_verbose=true problem/ck_impl_frozen_op=true \
    problem/ck_impl_warm=true
run p3b_warmonly $NEW time/nlim=20 output3/dt=1e-6 $SPH problem/ck_implicit=true \
    problem/ck_impl_verbose=true problem/ck_impl_warm=true

echo "=== (c) cost, 100 cycles, no dumps"
run p3c_semi  $NEW time/nlim=100 output3/dt=1e30 $SPH
run p3c_p2    $NEW time/nlim=100 output3/dt=1e30 $SPH problem/ck_implicit=true \
    problem/ck_impl_verbose=true
run p3c_frz   $NEW time/nlim=100 output3/dt=1e30 $SPH problem/ck_implicit=true \
    problem/ck_impl_verbose=true problem/ck_impl_frozen_op=true
run p3c_frzn  $NEW time/nlim=100 output3/dt=1e30 $SPH problem/ck_implicit=true \
    problem/ck_impl_verbose=true problem/ck_impl_frozen_op=true \
    problem/ck_impl_frozen_cof=false
run p3c_warm  $NEW time/nlim=100 output3/dt=1e30 $SPH problem/ck_implicit=true \
    problem/ck_impl_verbose=true problem/ck_impl_frozen_op=true \
    problem/ck_impl_warm=true
echo "=== done"
