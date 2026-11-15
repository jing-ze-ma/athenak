#!/bin/bash
# run_all.sh <athena> [npar=24]: every run of runs.txt (1 rank each), arm CC1, then score.
H=$(cd "$(dirname "$0")" && pwd); B=$(readlink -f $1); NP=${2:-24}
KEYS="scheme = clean;clean_layout = cell;clean_closure = extrap;clean_maxpic = 4;clean_pc_sweeps = 2"
while read -r dk f t; do
  grep -q '^rc=0' $H/suite/$dk/CC1_x$f/run.log 2>/dev/null && continue
  while [ $(jobs -rp | wc -l) -ge $NP ]; do sleep 20; done
  $H/run_ladder.sh CC1 $B "$KEYS" $dk $f $t &
done < $H/runs.txt
wait
cd $H && python3 score_all.py > score.txt
for l in suite/*/CC1_x*/run.log; do
  echo "$l | $(grep '^CLEAN' $l | tail -1) | $(grep '^wall' $l) | $(tail -1 $l)"
done > clean_lines.txt
