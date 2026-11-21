#!/bin/bash
# apply_v207.sh (INTEG11): keys of v20.7 for HG / RSG / R136 = the v20.5 keys unchanged (apply_v205.sh: v20.4 option B +
# clean_stag_flim=false named), as v20.6. v20.7 adds the LEDGER fix (no key) and clean_pass_damp (default 0 = off, not
# named here; the USER chooses N later), so the production restarts run bitwise to v20.6 / v20.5.
# Nothing is appended to common_scheme.txt: this script regenerates the three keys files from the
# current uniform keys (which carry the v20.5 block) and reports any difference; it changes nothing when v20.5 was applied.
# DRY=<dir>: write the three keys files to <dir> only. Run apply_v205.sh first if the v20.5 block is missing.
S=/orion/ptmp/jinma/clean_1115/stag2; U=$S/uniform/keys
grep -q '^# --- v20.5 ' $U/common_scheme.txt || { echo "apply_v207: v20.5 block missing in common_scheme.txt; run apply_v205.sh first"; exit 1; }
gen() {  # gen OUTDIR: common + <star>_problem -> OUTDIR/<star>_latest.txt (first value of a key wins)
  for s in hg rsg r136; do
    cat $U/common_scheme.txt $U/${s}_problem.txt | grep -v '^#' | grep -v '^[[:space:]]*$' | awk -F= '!seen[$1]++' > $1/${s}_latest.txt; done; }
T=${DRY:-$(mktemp -d)}; mkdir -p $T; gen $T
for s in hg rsg r136; do
  if cmp -s $T/${s}_latest.txt $S/${s}_latest.txt; then echo "$s: keys unchanged (= v20.5)"
  else echo "$s: differs from $S/${s}_latest.txt:"; diff $S/${s}_latest.txt $T/${s}_latest.txt; fi; done
[ -n "$DRY" ] && exit 0
for s in hg rsg r136; do f=$S/${s}_latest.txt; cmp -s $T/${s}_latest.txt $f && continue
  cp $f $f.bak_v205; cp $T/${s}_latest.txt $f; a=$S/integ/keys/k1_${s}.args; cp $a $a.bak_v205; tr '\n' ' ' < $f > $a; done
rm -rf $T; echo "apply_v207 done"
