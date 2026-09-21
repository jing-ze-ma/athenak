#!/bin/bash
# MILESTONE 3b PHASE G -- the TRANSVERSE realizability limiter
# (<rad_m1>/implicit_trans_limit = lp).  Serial CPU, build_cpu_box.
#   usage: bash tests_m1/runs_3b8/run_3b8.sh <arm> [overrides...]
set -u
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${EXE:-$ROOT/build_cpu_box/src/athena}
INP=$ROOT/tests_m1/runs_3b8/he_slab_m1_2d.athinput
OUT=$ROOT/tests_m1/runs_3b8
BASE="rad_m1/implicit_offdiag=operator rad_m1/implicit_closure_lag=step"
name=$1; shift
dir=$OUT/$name
rm -rf "$dir"; mkdir -p "$dir"
(cd "$dir" && timeout 1200 "$EXE" -i "$INP" -d . $BASE "$@" > log.txt 2>&1)
echo "=== $name rc=$?"; tail -2 "$dir/log.txt"
