#!/bin/bash
# run_all.sh <bindir> [npar=16]: every suite line (two 1-rank runs each) -> table.txt
PKG=$(cd "$(dirname "$0")" && pwd); BD=$(readlink -f $1)
grep . $PKG/suite.txt | xargs -d '\n' -P ${2:-16} -I{} bash -c '
  IFS="|" read -r n p d k l <<< "{}"
  '$PKG'/tools/cmp_ab.sh '$PKG' '$BD' "$n" "$p" "$d" "$k" "$l"' | sort > $PKG/table.txt
awk '{print $2}' $PKG/table.txt | sort | uniq -c
