#!/bin/bash
# run_one.sh <dir> <binary> <tmax_s> : run, analyse, prune dumps
d=$1; bin=$2; tmax=$3
cd $d || exit 1
[ -f metrics.json ] && exit 0
mkdir .lock 2>/dev/null || exit 0
rm -rf bin run.*hst run.log run-errs.dat out.log
t0=$(date +%s)
timeout $tmax $bin -i run.athinput > out.log 2>&1
rc=$?
t1=$(date +%s)
echo "{\"rc\": $rc, \"wall_s\": $((t1-t0)), \"timeout\": $([ $rc -eq 124 ] && echo true || echo false)}" > status.json
python3 /resnick/groups/carnegie_poc/jingze/lowbeta_1102/analyze.py $d > analyze.err 2>&1
rm -rf bin core.* 2>/dev/null
