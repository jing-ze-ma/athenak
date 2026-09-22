#!/bin/bash
# Field-amplitude arm at n = 32: does the L1(B) growth rate scale with the SPURIOUS
# VELOCITY (which is the Maxwell-stress truncation, so v ~ b0c^2) or not?
set -u
H="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$H/out" || exit 1
for b in 0.15 0.45; do
  tag="b$b"
  "$H/athena_diag" -i "$H/loop_hist.athinput" -d . job/basename="$tag" \
    mesh/nx1=32 mesh/nx2=32 mesh/nx3=32 \
    meshblock/nx1=32 meshblock/nx2=32 meshblock/nx3=32 \
    problem/conv_nband=4 problem/b0c=$b > "$H/logs/$tag.log" 2>&1 &
done
wait
