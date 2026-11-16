#!/bin/bash
# cmp.sh ARM BASE: files (bin, tab, hst) of each run of ARM that differ from arm BASE
cd "$(dirname "$0")"
for d in cc/suite/*/${1}_x* runs/*/${1}_* runs/thin/${1}_*_32 cs/*_$1 fb/runsf/flz/?${1}_32 rhse/$1 rw/*_$1 rst/runs/${1}_*; do
  [ -d $d ] || continue
  o=$(echo $d | sed "s#/${1}_x#/${2}_x#; s#/${1}_\([a-z0-9]*\)_32\$#/${2}_\1_32#; s#/${1}_\([0-9]*\)\$#/${2}_\1#; s#_${1}\$#_${2}#; s#/q${1}_#/q${2}_#; s#/Q${1}_#/Q${2}_#; s#rhse/${1}\$#rhse/${2}#; s#runs/${1}_#runs/${2}_#")
  n=0; b=0
  for f in $(cd $d && find . -name '*.bin' -o -name '*.tab' -o -name '*.hst' | grep -v '/rst/'); do
    n=$((n+1)); cmp -s $d/$f $o/$f || b=$((b+1)); done
  echo "$d files $n differ $b"
done | sort -u
echo "CLEAN WARNING: $(cat $(find cc/suite runs cs fb/runsf rhse rw rst/runs -path "*$1*" \( -name run.log -o -name log.txt -o -name '?.log' \)) 2>/dev/null | grep -c 'CLEAN WARNING')"
