#!/bin/bash
# tests_m1/runs_3d_edd -- the EDDINGTON-consistent He slab column (V3edd).
#   usage: bash run_3d.sh <arm> <athinput> [overrides...]
set -u
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
EXE=$HERE/athena_56c97566
name=$1; shift
inp=$1; shift
dir=$HERE/$name
rm -rf "$dir"; mkdir -p "$dir"
(cd "$dir" && timeout 1500 "$EXE" -i "$inp" -d . "$@" > log.txt 2>&1)
echo "=== $name rc=$?"; tail -2 "$dir/log.txt"
