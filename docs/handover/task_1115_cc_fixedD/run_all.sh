#!/bin/bash
# run_all.sh <athena> [npar=32]: every line of runs.txt (one 1-rank run each), then score.
cd "$(dirname "$0")" || exit 1
export ATHENA=$(readlink -f $1) OMP_NUM_THREADS=1
grep -v -F -f <(grep -l '^rc=0' runs/*/log.txt 2>/dev/null | sed 's,/log.txt$, \&\&,') runs.txt \
  | xargs -P ${2:-32} -d '\n' -I{} bash -c {}
python3 rwfarm.py score cases.txt runs > score.txt
python3 rwfarm.py table cases.txt runs > table.md
for d in runs/*/; do echo "$(basename $d) | $(grep '^CLEAN' $d/log.txt | tail -1) |" \
  "$(grep -o 'E-floor created energy dropped[^,]*' $d/log.txt | tail -1) | $(tail -1 $d/log.txt)"
done > clean_lines.txt
