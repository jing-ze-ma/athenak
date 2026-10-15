#!/bin/bash
# seed_cr2.sh CASE=hg|e1 (login node or CPU job): COPY the newest VALID production rst into $W/$CASE/seed/rst (read-only use of the
# production dir; nothing there is written, linked or moved). All CR2 arms (gate, window, ROL2U chain) start from this one copy.
set -e
CASE=$1; source "$(dirname "$0")/cr2_common.sh"
mkdir -p $SEED/rst
[ -n "$(seedrst)" ] && { echo "seed exists: $(seedrst)"; exit 0; }
P=$(python3 $SCR/rst_pick.py $PROD/rst | head -1); read -r F T NC TL <<< "$P"
[ -n "$F" ] || { echo "no valid production rst in $PROD/rst"; exit 1; }
age=$(( $(date +%s) - $(stat -c %Y "$F") )); [ $age -lt 300 ] && { echo "$F only $age s old: retry later"; exit 1; }
cp -p "$F" $SEED/rst/                                     # a COPY (production prune keeps 3; never hard-link production files)
cmp "$F" $SEED/rst/$(basename "$F") && echo "copied $F (t $T, cycle $NC) -> $SEED/rst ($(du -h $SEED/rst/*.rst | cut -f1))"
md5sum $SEED/rst/*.rst > $SEED/rst.md5
echo "$F $T $NC $TL" > $SEED/seed.txt
echo "par_end offset of the seed rst: $(parend $SEED/rst/$(basename "$F")) (limit 45056)"
