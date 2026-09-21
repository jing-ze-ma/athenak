#!/bin/bash
# tests_m1/runs_3g_newton_T -- MILESTONE 3g gates.
#   usage: bash run.sh <exe> <name> <athinput> [overrides...]
set -u
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
EXE=$1; shift
name=$1; shift
inp=$1; shift
dir=$HERE/$name
rm -rf "$dir"; mkdir -p "$dir"
(cd "$dir" && timeout 1500 "$EXE" -i "$inp" -d . "$@" > log.txt 2>&1)
rc=$?
echo "=== $name rc=$rc"
grep -E "zone-cycles/cpu_second|cpu time used" "$dir/log.txt" | tail -2
grep -E "implicit transport:|gas coupling:|eos_cache:" "$dir/log.txt" | tail -4
