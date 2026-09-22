#!/bin/bash
# Bitwise gate: HEAD (ref) vs the FOP-tag build (new), 20 cycles, hst + bin dumps.
# Usage: ./gate.sh
set -u
REF=/viper/u2/jinma/ATHENAK/bench/cktag_0922/athena.cpu.head
NEW=/viper/u2/jinma/ATHENAK/bench/cktag_0922/athena.cpu.tag
IN1=/viper/u2/jinma/ATHENAK/athenak/inputs/tests/dhj_ck_spherical.athinput
IN2=/viper/u2/jinma/ATHENAK/athenak/inputs/tests/dhj_ck_implicit.athinput
BASE="time/nlim=20 output1/dt=1e-30 output3/dt=1e-30 output4/dt=1e30"
run() { local d=$1; shift; local b=$1; shift; local in=$1; shift
  rm -rf "$d"; mkdir -p "$d"
  ( cd "$d" && $b -i $in $BASE "$@" > run.log 2>&1 ); echo "rc=$? $d"; }
gate() { local name=$1; shift; local in=$1; shift
  run ${name}_ref $REF $in "$@"
  run ${name}_new $NEW $in "$@"
  local ok=1
  for f in $(cd ${name}_ref && ls *.hst bin/*.bin 2>/dev/null); do
    cmp -s ${name}_ref/$f ${name}_new/$f || { echo "DIFF $name $f"; ok=0; }
  done
  local n=$(cd ${name}_ref && ls *.hst bin/*.bin 2>/dev/null | wc -l)
  echo "GATE $name: $([ $ok = 1 ] && echo BITWISE || echo FAIL) ($n files)"; }

gate a_sph_off $IN1
gate b_sph_on  $IN1 problem/ck_spherical=true
gate c_imp     $IN2 problem/ck_implicit=true problem/rt_use_cons=true
gate d_imp_fop $IN2 problem/ck_implicit=true problem/rt_use_cons=true problem/ck_impl_frozen_op=true
gate e_imp_sph $IN2 problem/ck_implicit=true problem/rt_use_cons=true problem/ck_spherical=true \
                    problem/ck_impl_frozen_op=true
echo GATEDONE
