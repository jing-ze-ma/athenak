#!/bin/bash
# T10 (the radiation-modified acoustic wave) with and without the TRANSVERSE
# realizability limiter.  Serial CPU, build_cpu_m1.
set -u
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=$ROOT/build_cpu_m1/src/athena
INP=$ROOT/tests_m1/runs_3b8/rad_m1_radwave.athinput
OUT=$ROOT/tests_m1/runs_3b8
TLIM2=1.184313050927584
TLIMD=0.8374357893586236
ODT2=0.0370097828
ODTD=0.0261698684
IMPL="rad_m1/transport=implicit rad_m1/implicit_solver=bicgstab \
rad_m1/implicit_offdiag=operator rad_m1/implicit_closure_lag=step"
run() { local name=$1; shift; local dir=$OUT/$name
  rm -rf "$dir"; mkdir -p "$dir"
  (cd "$dir" && timeout 1200 "$EXE" -i "$INP" -d . "$@" > log.txt 2>&1)
  echo "=== $name rc=$?"; }
for L in off on; do
  [ $L = on ] && LIM="rad_m1/implicit_trans_limit=lp" || LIM=""
  run a_x1_$L mesh/nx2=1 meshblock/nx2=1 problem/radwave_dir=x1 \
    rad_m1/transport=implicit_x1 time/tlim=$TLIM2 output1/dt=$ODT2
  run c_x2_$L mesh/nx1=4 mesh/nx2=64 meshblock/nx1=4 meshblock/nx2=64 \
    problem/radwave_dir=x2 $IMPL $LIM time/tlim=$TLIM2 output1/dt=$ODT2
  run e_xy_$L mesh/nx1=64 mesh/nx2=64 meshblock/nx1=64 meshblock/nx2=64 \
    problem/radwave_dir=xy $IMPL $LIM time/tlim=$TLIMD output1/dt=$ODTD
done
python3 "$ROOT/tests_m1/t10_radwave.py" \
  --arm $OUT/a_x1_off,x1,a_x1_off, \
  --arm $OUT/c_x2_off,x2,c_x2_off,a_x1_off \
  --arm $OUT/e_xy_off,xy,e_xy_off, \
  --arm $OUT/c_x2_on,x2,c_x2_on,a_x1_off \
  --arm $OUT/e_xy_on,xy,e_xy_on,
