#!/bin/bash
# Resolution scan of the cs field loop (cs_test iprob 11) to t = 2.0, hst every 0.0625.
set -u
H="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$H/out" || exit 1
for n in 16 32 64; do
  nb=$(( n / 8 ))
  "$H/athena_diag" -i "$H/loop_hist.athinput" -d . job/basename="n$n" \
    mesh/nx1=$n mesh/nx2=$n mesh/nx3=$n \
    meshblock/nx1=$n meshblock/nx2=$n meshblock/nx3=$n \
    problem/conv_nband=$nb > "$H/logs/n$n.log" 2>&1 &
done
wait
