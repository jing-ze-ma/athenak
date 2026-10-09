#!/bin/bash
# Delta test runner poller (TASK-2026-10-08-delta-test-runner): every 15 min fetch the radiation branches and
# exit (printing them) as soon as a batch file not listed in seen.txt appears.
B=/work/nvme/bivj/jma20/delta_1008; R=/u/jma20/ATHENAK/athenak
while true; do
  new=""
  for br in sp-blend2-1008 rad-beam-1008; do
    git -C $R fetch -q origin $br 2>/dev/null || continue
    for f in $(git -C $R ls-tree --name-only origin/$br docs/handover/ | grep -E 'TASK-2026-10-08-delta-.*-batch.*\.md$'); do
      grep -qxF "$br $f" $B/runner/seen.txt || new="$new$br $f\n"
    done
  done
  [ -n "$new" ] && { printf "NEW BATCH(ES) $(date '+%F %T'):\n$new"; exit 0; }
  sleep 900
done
