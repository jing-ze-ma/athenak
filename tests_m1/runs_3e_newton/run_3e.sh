#!/bin/bash
# MILESTONE 3e -- ANDERSON ACCELERATION of the lagged-closure Picard iteration
# (<rad_m1>/implicit_accel = anderson).  Serial CPU, build_cpu_box3.
#   usage: bash tests_m1/runs_3e_newton/run_3e.sh <arm> [overrides...]
set -u
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${EXE:-$ROOT/build_cpu_box3/src/athena}
INP=${INP:-$ROOT/tests_m1/runs_3e_newton/he_slab_m1_2d.athinput}
OUT=$ROOT/tests_m1/runs_3e_newton
BASE="rad_m1/implicit_offdiag=operator"
name=$1; shift
dir=$OUT/$name
rm -rf "$dir"; mkdir -p "$dir"
(cd "$dir" && timeout 1500 "$EXE" -i "$INP" -d . $BASE "$@" > log.txt 2>&1)
echo "=== $name rc=$?"; tail -2 "$dir/log.txt"
