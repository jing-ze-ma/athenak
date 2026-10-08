#!/bin/bash
# usage: bash metrics.sh <rundir> ...
# s/cycle from the elapsed= lines (time/ndiag=1), skipping the first C0SKIP (default 4) cycles of the
# run; then the rad_m1 solver summary lines and the last hydro.hst row.
for d in "$@"; do
  L=$d/run.log
  echo "== $d"
  grep -a "^elapsed=" "$L" | sed 's/elapsed=\([^ ]*\) cycle=\([0-9]*\).*/\1 \2/' | \
    awk -v s="${C0SKIP:-4}" 'NR==1{c00=$2} $2>=c00+s{if(!n){e0=$1;k0=$2} n++; e1=$1;k1=$2}
      END{if(k1>k0) printf "s/cycle %.4f over %d cycles (%d-%d)\n",(e1-e0)/(k1-k0),k1-k0,k0,k1;
          else print "no cycles"}'
  grep -a "Picard iterations mean\|inner iterations mean\|breakdowns=\|positivity fallbacks" "$L"
  grep -a "final 7-point\|FATAL" "$L"
  grep -a -o "|F| > c E scaled back.*" "$L"
  tail -1 "$d"/*.hydro.hst
done
