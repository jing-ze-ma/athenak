#!/bin/bash
# one order run: one.sh <rundir> <input>   (CPU binary $ATHENA_CPU, OMP 1, nice)
R=$1; IN=$2; shift 2
mkdir -p "$R"
cd "$R" || exit 1
t0=$(date +%s.%N)
OMP_NUM_THREADS=1 nice -n 10 "$ATHENA_CPU" -i "$IN" -d "$R" "$@" > run.log 2>&1
rc=$?
t1=$(date +%s.%N)
echo "rc=$rc wall=$(echo "$t1 - $t0" | bc)" > wall.txt
echo "$ATHENA_CPU $*" > cmd.txt
