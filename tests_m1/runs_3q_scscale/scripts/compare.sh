#!/bin/bash
# usage: compare.sh <dirA> <dirB> : byte comparison of every output file; the history
# files separately; for the others the number of differing bytes and the last offset
cd /viper/ptmp2/jinma/scscale_0923/cpu
a=$1; b=$2; n=0; same=0; out=""; hst=""
for f in $(cd $a && find . -type f \( -name "*.bin" -o -name "*.rst" -o -name "*.hst" \) | sort); do
  n=$((n+1))
  if cmp -s $a/$f $b/$f; then same=$((same+1)); r=same; else
    nd=$(cmp -l $a/$f $b/$f 2>/dev/null | wc -l); mo=$(cmp -l $a/$f $b/$f 2>/dev/null | tail -1 | awk '{print $1}')
    out="$out $(basename $f):${nd}B@<=$mo"; r=DIFF; fi
  case $f in *.hst) hst="$hst $(basename $f .hst)=$r";; esac
done
echo "$a vs $b: $same/$n identical; hst:$hst;$out | $(tail -1 $a/log.txt) $(tail -1 $b/log.txt)"
