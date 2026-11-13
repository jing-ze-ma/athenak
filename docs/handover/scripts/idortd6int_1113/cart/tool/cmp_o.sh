#!/bin/bash
# idortd6int-1113: D6 (p, d6int build) vs FIN (f, wt_fin build) per row, as Raven jobs/cmp_d6.sh.  usage: cmp_o.sh <W>
# env BINCONV = <worktree>/vis/python (cmpout.py reads .bin files with bin_convert).
# (1) row lines: FLAG_RC (D6 rc!=0 or FATAL), FLAG_NC (D6 NON-CONVERGED > 0), FLAG_PICARD (D6 Picard mean > 1.25 x FIN)
# (2) cmpout.py FIN -> D6 outputs (BITWISE, or DIFF with maxrel of the differing data)
W=$1; C=$(cd "$(dirname "$0")/.." && pwd)
pic() { grep -h -o 'Picard iterations mean=[^ ]* max=[^ ]* NON-CONVERGED=[^ ]*' "$1/all.log" 2>/dev/null | tail -1 | sed 's/Picard iterations //'; }
row() { [ -d "$1" ] || { echo MISSING; return; }; echo "rc=$(cat $1/rc.txt 2>/dev/null) fatal=$(grep -c FATAL $1/all.log 2>/dev/null) dtr=$(grep -c '^DT_RETRY cycle' $1/all.log 2>/dev/null) $(pic $1) t=$(cut -d' ' -f1 $1/time.txt 2>/dev/null)"; }
dirs=$(cd "$W" && ls -d runs/*/*/*/p runs/*/*/p runs_omp/*/*/*/p 2>/dev/null | sed 's#/p$##' | sort -u)
echo "== (1) rows: D6 (p) | FIN (f)"
for r in $dirs; do a=$(row $W/$r/p); b=$(row $W/$r/f); fl=""
  echo "$a" | grep -q 'rc=0 fatal=0' || fl="$fl FLAG_RC"
  na=$(echo "$a" | grep -o 'NON-CONVERGED=[0-9.e+]*' | cut -d= -f2); awk -v x="${na:-0}" 'BEGIN{exit !(x+0 > 0)}' && fl="$fl FLAG_NC"
  ma=$(echo "$a" | grep -o 'mean=[0-9.e+-]*' | cut -d= -f2); mb=$(echo "$b" | grep -o 'mean=[0-9.e+-]*' | cut -d= -f2)
  [ -n "$ma" ] && [ -n "$mb" ] && awk -v x=$ma -v y=$mb 'BEGIN{exit !(x > 1.25*y)}' && fl="$fl FLAG_PICARD"
  echo "$r | D6 $a | FIN $b |$fl"; done
echo "== (2) outputs FIN -> D6 (cmpout.py)"
for r in $dirs; do [ -d "$W/$r/f" ] && [ -d "$W/$r/p" ] && echo "$r $(python3 $C/tool/cmpout.py $W/$r/f $W/$r/p 2>&1 | tail -1)"; done
