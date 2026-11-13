#!/bin/bash
# idortd6int-1113 (Orion CPU): one IDORT row, arms concurrently on pre-resolved decks (no key processing here).
#   q = prod-1118 = D6 deck + $EXTRA keys on the prod-1118 build (fork/prod-1118 1d1972805)  binary: $EXE_Q
# usage: one_q.sh <row> <nt>   env: C (= .../idortd6int_1113/cart of this branch) W (run root) EXE_Q EXTRA [NLIM]
# Decks were resolved on viper with tool/resolve (mkp.py + timekeys.sh, CAND=d6 / fin); identical to the Raven d6fast / FIN decks
# up to the shock-reference path (now /orion/ptmp/jinma/idortd6int_1113/refs/...). nt is the core reservation; the binary is serial
# (OMP_NUM_THREADS=$nt only matters for an OpenMP build).
r=$1; nt=$2; d="${W:?}/runs/${r:?}"; mkdir -p "$d"
arm() { x=$1; case $x in q) dk=$C/decks_d6/$r; B=$EXE_Q;; *) return;; esac
  { [ -f "$d/$x/rc.txt" ] && [ "$(cat "$d/$x/rc.txt")" = 0 ]; } && return
  mkdir "$d/.lock_$x" 2>/dev/null || return
  mkdir -p "$d/$x"; cd "$d/$x" || return; rm -f rc.txt all.log run.log; cp "$dk/in.athinput" "$dk/timekeys.txt" .
  t0=$(date +%s); md5sum "$B" > bin.txt
  OMP_NUM_THREADS=$nt OMP_PROC_BIND=false "$B" -i in.athinput -d "$d/$x" ${EXTRA:?} ${NLIM:+time/nlim=$NLIM} > run.log 2>&1; echo $? > rc.txt
  cp run.log all.log; echo "$(( $(date +%s) - t0 )) s nt=$nt" > time.txt
  rm -rf "${d:?}/${x:?}/rst"; rmdir "$d/.lock_$x"; }
for x in ${ARMS:-q}; do arm "$x" & done; wait
cd "$d" || exit 1
pic() { grep -h -o 'Picard iterations mean=[^ ]* max=[^ ]* NON-CONVERGED=[^ ]*' "$1/all.log" | tail -1 | sed 's/Picard iterations //'; }
line="$r |"; for x in ${ARMS:-q}; do line="$line $x rc=$(cat $x/rc.txt) fatal=$(grep -c FATAL $x/all.log) dtr=$(grep -c '^DT_RETRY cycle' $x/all.log) unk=$(grep -c -i 'unread\|not read' $x/all.log) $(pic $x) t=$(cut -d' ' -f1 $x/time.txt) |"; done
echo "$line" > sum_q.txt; cat sum_q.txt
echo "$EXTRA" > "$d/q/extra.txt"
