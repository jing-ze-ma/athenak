#!/bin/bash
# addm2.sh: prepend W189 m2 dumps (complete: older than 2 min) not yet in queue3.txt
O=/resnick/groups/carnegie_poc/jingze/fe_mhd_1103; B=/resnick/groups/carnegie_poc/jingze/r1planets_1101/wasp189b/runs/m2/bin
new=$(for f in $(find $B -name 'dhj.mhd_w_bcc.*.bin' -mmin +2 | sort -r); do n=$(basename $f .bin); n=${n##*.}; [ $n = 00000 ] && continue
  [ $((10#$n)) -le 3 ] && continue; grep -qw "w189m2_$n" $O/queue3.txt || echo "w189 $f w189m2_$n"; done)
[ -n "$new" ] && { { echo "$new"; cat $O/queue3.txt; } > $O/queue3.new && mv $O/queue3.new $O/queue3.txt; echo "$new"; }
