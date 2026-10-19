#!/bin/bash
# usage: collect.sh <athenak checkout of bsg-files-1009 (has docs/handover/csedge-1019/csedge_diag.py + vis/python)>
# Runs the CSEDGE diag on each arm's dumps so far and copies the TASK files into <wt>/docs/handover/csedge-1019/caltech/<arm>/,
# plus the restart-diff (TASK B) results into <wt>/docs/handover/csedge-1019/caltech/rdiff/. Then commit + push from <wt>.
W=${1:?checkout}; P=/resnick/groups/carnegie_poc/jingze/rsg_g25_prod; C=$P/csedge_1019
for a in B F H N V B24; do D=$C/$a; O=$W/docs/handover/csedge-1019/caltech/$a; [ -d $D ] || continue; mkdir -p $O
  cp -f $D/rsg.user.hst $D/rsg.hydro.hst $D/cmdline.txt $O/ 2>/dev/null; tail -30 $D/run.log > $O/run.log.tail30 2>/dev/null
  python3 $W/docs/handover/csedge-1019/csedge_diag.py $W/vis/python $D/bin rsg > $O/diag.txt 2>&1
  # s/cycle from the elapsed lines (ndiag=100): last-minus-first elapsed / cycles
  grep -E '^elapsed' $D/run.log | awk -F'[= ]' 'NR==1{e0=$2;c0=$4} {e=$2;c=$4} END{if(c>c0) printf "s/cycle %.3f over cycles %d-%d (t=%s)\n",(e-e0)/(c-c0),c0,c,$6; else print "no cycles yet"}' > $O/s_per_cycle.txt
  J=$(ls $C/cse_$a.*.out 2>/dev/null | tail -1); [ -n "$J" ] && cp -f $J $O/slurm.out
  echo "$a: $(cat $O/s_per_cycle.txt); dumps $(ls $D/bin/rsg.hydro_w.*.bin 2>/dev/null | wc -l)"; done
R=$P/det/rdiff_1019; O=$W/docs/handover/csedge-1019/caltech/rdiff; mkdir -p $O
cp -f $R/rdiff.sh $R/compare.sh $O/ ; for d in R2 CV RV; do [ -f $R/$d/cmdline.txt ] && cp -f $R/$d/cmdline.txt $O/cmdline_$d.txt; done
J=$(ls $R/rdiff.*.out 2>/dev/null | tail -1); [ -n "$J" ] && grep -vE '^ *(Loading|mpfr|zlib|xz|libxml)' $J > $O/rdiff.out
bash $R/compare.sh > $O/compare.txt 2>&1; cat $O/compare.txt
