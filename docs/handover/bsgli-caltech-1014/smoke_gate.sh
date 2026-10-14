#!/bin/bash
# bsgli-caltech-1014 smoke_gate.sh SMOKE_DIR: the production gate (= e1caltech-1013/smoke_gate.sh with the Coff check added and the
# he_ic_equil check replaced by the Linhao IC/BC lines). It compares DATA only, never whole bin files: the headers carry nlim / file
# counters / offsets that legitimately differ between a continuous run and a chain (your 10-13 BSG false positive).
#   chaincmp.py C L L : every hst DATA row of C (by time) present in the chain, identical text, none missing / extra;
#                       every bin dump of C (stem, time) present in the chain, every variable's data array byte-identical.
#   cmplast.py C L    : last dump per output, every variable's data array byte-identical.
#   cmplast.py C Coff : the same for halo_ipc off (vet_gd_halo_ipc = true is kept only if this is BITWISE).
#   hst files         : rows without '#' lines, cmp (C vs L and C vs Coff).
#   logs              : FATAL 0, NaN 0, DIVERGED 0, NON-CONVERGED 0 in C, L (both links), Coff; he_ic_linhao + he_bc_linhao lines in C.
# Writes SMOKE_DIR/VERIFIED (and prints GATE PASS) only if ALL pass; else prints GATE FAIL and the reasons, exit 1. Never edits run dirs.
S=$1; W=$(dirname "$(readlink -f "$0")"); PY=${PY:-python3}; bad=""
export BINCONV_DIR=${BINCONV_DIR:-$W}
cd "$S" || exit 1
$PY $W/chaincmp.py C L L | tee gate_chaincmp.txt; grep -q "CHAINCMP BITWISE" gate_chaincmp.txt || bad="$bad chaincmp"
[ "$(grep -c 'missing 0 \[\], extra 0 \[\], differing 0' gate_chaincmp.txt)" -ge 2 ] || bad="$bad hst-rows"
nd=$(grep -o "^bin: C dumps [0-9]*" gate_chaincmp.txt | awk '{print $4}'); [ "${nd:-0}" -ge 4 ] || bad="$bad too-few-dumps(${nd:-0})"
$PY $W/cmplast.py C L | tee gate_cmplast.txt | tail -2; grep -q "CMPLAST BITWISE" gate_cmplast.txt || bad="$bad cmplast"
$PY $W/cmplast.py C Coff | tee gate_coff.txt | tail -2; grep -q "CMPLAST BITWISE" gate_coff.txt || bad="$bad coff(halo_ipc)"
for h in $(cd C; ls *.hst); do
  cmp -s <(grep -v '^#' C/$h) <(grep -v '^#' L/$h) && echo "HST $h BITWISE" || bad="$bad hst:$h"
  cmp -s <(grep -v '^#' C/$h) <(grep -v '^#' Coff/$h) && echo "COFF HST $h BITWISE" || bad="$bad coffhst:$h"
done
for f in C/run.log.* L/run.log.* Coff/run.log.*; do
  [ -f "$f" ] || continue
  grep -q FATAL $f && bad="$bad FATAL:$f"; grep -q DIVERGED $f && bad="$bad DIVERGED:$f"
  tail -n 60 $f | grep -qiw nan && bad="$bad NaN:$f"
  n=$(grep -o "NON-CONVERGED=[0-9.e+-]*" $f | tail -1 | cut -d= -f2); [ -n "$n" ] && awk -v x="$n" 'BEGIN{exit !(x+0>0)}' && bad="$bad NC=$n:$f"
done
grep -q "he_ic_linhao:" C/run.log.* || bad="$bad no-he_ic_linhao-line"
grep -q "he_bc_linhao" C/run.log.* || bad="$bad no-he_bc_linhao-line"
for d in C L Coff; do [ -f $d/STOP ] && bad="$bad STOP:$d"; [ -f $d/NOGO ] && bad="$bad NOGO:$d"; done
if [ -z "$bad" ]; then echo "GATE PASS $(date +%F_%T)" | tee VERIFIED; else echo "GATE FAIL:$bad"; rm -f VERIFIED; exit 1; fi
