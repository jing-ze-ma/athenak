#!/bin/bash
# tests_ck_implicit phase 4: the JACOBIAN REUSE (problem/ck_impl_reuse_jac) and the
# BETTER INITIAL GUESS (problem/ck_impl_seed).  See README_phase4.md.
set -u
R=/viper/u2/jinma/ATHENAK/athenak
NEW=$R/build_ckq/src/athena
REF=/viper/u2/jinma/ATHENAK/bench/ckq_0922/ref_src/build_ref/src/athena
IN=$R/inputs/tests/dhj_ck_implicit.athinput
INS=$R/inputs/tests/dhj_ck_spherical.athinput
BASE="output1/dt=1e30 output2/dt=1e30 output4/dt=1e30 problem/rt_use_cons=true"
SPH="problem/ck_spherical=true problem/ck_beam_sph=true"
IMP="problem/ck_implicit=true problem/ck_impl_verbose=true problem/ck_impl_frozen_op=true"

run() {  # run <dir> <binary> <input> <extra args...>
  local d=$1; shift
  local b=$1; shift
  local i=$1; shift
  rm -rf "$d"; mkdir -p "$d"
  ( cd "$d" && /usr/bin/time -f "%e" -o t.txt $b -i $i "$@" > run.log 2>&1 )
  echo "$? $d $(cat $d/t.txt 2>/dev/null)"
}

if [ "${1:-all}" = a ] || [ "${1:-all}" = all ]; then
echo "=== (a) default off: new vs HEAD reference, dhj_ck_spherical, 20 cycles"
for sp in false true; do
  run p4a_ref_sp$sp $REF $INS time/nlim=20 output3/dt=1e-6 problem/ck_spherical=$sp
  run p4a_new_sp$sp $NEW $INS time/nlim=20 output3/dt=1e-6 problem/ck_spherical=$sp
done
fi

if [ "${1:-all}" = b ] || [ "${1:-all}" = all ]; then
echo "=== (b) 20 cycles with dumps: state and passes, each lever against the baseline"
run p4b_base $NEW $IN time/nlim=20 output3/dt=1e-6 $BASE $SPH $IMP
for v in "problem/ck_impl_reuse_jac=1:r1" "problem/ck_impl_reuse_jac=2:r2" \
         "problem/ck_impl_seed=1:s1" "problem/ck_impl_seed=2:s2" \
         "problem/ck_impl_reuse_jac=1 problem/ck_impl_seed=2:b1" \
         "problem/ck_impl_reuse_jac=2 problem/ck_impl_seed=2:b2"; do
  a=${v%:*}; n=${v##*:}
  run p4b_$n $NEW $IN time/nlim=20 output3/dt=1e-6 $BASE $SPH $IMP $a
done
fi

if [ "${1:-all}" = c ] || [ "${1:-all}" = all ]; then
echo "=== (c) cost, 100 cycles, no dumps"
run p4c_semi $NEW $IN time/nlim=100 output3/dt=1e30 $BASE $SPH
run p4c_base $NEW $IN time/nlim=100 output3/dt=1e30 $BASE $SPH $IMP
for v in "problem/ck_impl_reuse_jac=1:r1" "problem/ck_impl_reuse_jac=2:r2" \
         "problem/ck_impl_seed=1:s1" "problem/ck_impl_seed=2:s2" \
         "problem/ck_impl_reuse_jac=1 problem/ck_impl_seed=2:b1" \
         "problem/ck_impl_reuse_jac=2 problem/ck_impl_seed=2:b2"; do
  a=${v%:*}; n=${v##*:}
  run p4c_$n $NEW $IN time/nlim=100 output3/dt=1e30 $BASE $SPH $IMP $a
done
fi
echo "=== done"
