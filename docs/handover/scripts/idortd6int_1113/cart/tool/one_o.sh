#!/bin/bash
# idortd6int-1113 (Orion CPU): one IDORT row, arms concurrently on pre-resolved decks (no key processing here).
#   p = D6 keys on the d6int-1113 build (tip, >= 62abe7ec1) deck: $C/decks_d6/<row>/in.athinput   binary: $EXE_P
#   f = FIN keys on the FIN build (421a3d024 wt_fin)  deck: $C/decks_fin/<row>/in.athinput  binary: $EXE_F
# usage: one_o.sh <row> <nt>   env: C (= .../idortd6int_1113/cart of this branch) W (run root) EXE_P EXE_F [ARMS="p f"] [NLIM]
# Decks were resolved on viper with tool/resolve (mkp.py + timekeys.sh, CAND=d6 / fin); identical to the Raven d6fast / FIN decks
# up to the shock-reference path (now /orion/ptmp/jinma/idortd6int_1113/refs/...). nt is the core reservation; the binary is serial
# (OMP_NUM_THREADS=$nt only matters for an OpenMP build).
r=$1; nt=$2; d="${W:?}/runs/${r:?}"; mkdir -p "$d"
arm() { x=$1; case $x in p) dk=$C/decks_d6/$r; B=$EXE_P;; f) dk=$C/decks_fin/$r; B=$EXE_F;; *) return;; esac
  { [ -f "$d/$x/rc.txt" ] && [ "$(cat "$d/$x/rc.txt")" = 0 ]; } && return
  mkdir "$d/.lock_$x" 2>/dev/null || return
  mkdir -p "$d/$x"; cd "$d/$x" || return; rm -f rc.txt all.log run.log; cp "$dk/in.athinput" "$dk/timekeys.txt" .
  t0=$(date +%s); md5sum "$B" > bin.txt
  OMP_NUM_THREADS=$nt OMP_PROC_BIND=false "$B" -i in.athinput -d "$d/$x" ${NLIM:+time/nlim=$NLIM} > run.log 2>&1; echo $? > rc.txt
  cp run.log all.log; echo "$(( $(date +%s) - t0 )) s nt=$nt" > time.txt
  rm -rf "${d:?}/${x:?}/rst"; rmdir "$d/.lock_$x"; }
for x in ${ARMS:-p f}; do arm "$x" & done; wait
cd "$d" || exit 1
pic() { grep -h -o 'Picard iterations mean=[^ ]* max=[^ ]* NON-CONVERGED=[^ ]*' "$1/all.log" | tail -1 | sed 's/Picard iterations //'; }
line="$r |"; for x in ${ARMS:-p f}; do line="$line $x rc=$(cat $x/rc.txt) fatal=$(grep -c FATAL $x/all.log) dtr=$(grep -c '^DT_RETRY cycle' $x/all.log) unk=$(grep -c -i 'unread\|not read' $x/all.log) $(pic $x) t=$(cut -d' ' -f1 $x/time.txt) |"; done
echo "$line" > sum.txt; cat sum.txt
