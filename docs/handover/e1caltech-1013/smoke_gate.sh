#!/bin/bash
# e1caltech-1013 smoke_gate.sh SMOKE_DIR: the production gate. Compares DATA only (never whole bin files: the headers carry
# nlim / file counters / offsets that legitimately differ between a continuous run and a chain -- your BSG false positive).
#   chaincmp.py C L L : every hst DATA row of C (by time) present in the chain, identical text, none missing / extra;
#                       every bin dump of C (stem, time) present in the chain, every variable's data array byte-identical.
#   cmplast.py C L    : last dump per output, every variable's data array byte-identical.
#   hst files         : rows without '#' lines, cmp.
#   logs              : rc 0 (link stdout), FATAL 0, NaN 0, DIVERGED 0, NON-CONVERGED 0 in C, L (both links), he_ic_equil ran.
# Writes SMOKE_DIR/VERIFIED (and prints GATE PASS) only if ALL pass; else prints GATE FAIL and the reasons. Never edits run dirs.
S=$1; W=$(dirname "$(readlink -f "$0")"); PY=${PY:-python3}; bad=""
cd "$S" || exit 1
$PY $W/chaincmp.py C L L | tee gate_chaincmp.txt; grep -q "CHAINCMP BITWISE" gate_chaincmp.txt || bad="$bad chaincmp"
[ "$(grep -c 'missing 0 \[\], extra 0 \[\], differing 0' gate_chaincmp.txt)" -ge 2 ] || bad="$bad hst-rows"
nd=$(grep -o "^bin: C dumps [0-9]*" gate_chaincmp.txt | awk '{print $4}'); [ "${nd:-0}" -ge 4 ] || bad="$bad too-few-dumps(${nd:-0})"
$PY $W/cmplast.py C L | tee gate_cmplast.txt | tail -2; grep -q "CMPLAST BITWISE" gate_cmplast.txt || bad="$bad cmplast"
for h in $(cd C; ls *.hst); do cmp -s <(grep -v '^#' C/$h) <(grep -v '^#' L/$h) && echo "HST $h BITWISE" || bad="$bad hst:$h"; done
for f in C/run.log.* L/run.log.*; do
  [ -f "$f" ] || continue
  grep -q FATAL $f && bad="$bad FATAL:$f"; grep -q DIVERGED $f && bad="$bad DIVERGED:$f"
  tail -n 60 $f | grep -qiw nan && bad="$bad NaN:$f"
  n=$(grep -o "NON-CONVERGED=[0-9.e+-]*" $f | tail -1 | cut -d= -f2); [ -n "$n" ] && awk -v x="$n" 'BEGIN{exit !(x+0>0)}' && bad="$bad NC=$n:$f"
done
grep -q he_ic_equil C/run.log.* || bad="$bad no-he_ic_equil-line"
[ -f C/STOP ] && bad="$bad STOP:C"; [ -f L/STOP ] && bad="$bad STOP:L"
if [ -z "$bad" ]; then echo "GATE PASS $(date +%F_%T)" | tee VERIFIED; else echo "GATE FAIL:$bad"; rm -f VERIFIED; exit 1; fi
