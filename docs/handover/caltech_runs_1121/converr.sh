#!/bin/bash
# converr.sh OUT [OUTDIR]  (orion CONVERR 2026-11-20): run converr.py on the CONVERR arms under OUT (the armc.sh run root).
# Needs ~25 GB RAM and ~10-20 min per star on a compute node: run it inside a holder, e.g.
#   srun --overlap --jobid=$HOLDER -w $NODE -N1 -n1 -c 8 -t 90 bash converr.sh $OUT
# Writes OUTDIR/converr_{rsg,hg}.txt (+ _3c) (default OUTDIR = $OUT/CONVERR); commit those text files only.
# The bins can be deleted after the text files are committed.
O=$1; R=${2:-$O/CONVERR}; mkdir -p $R; H=$(dirname "$(readlink -f "$0")")
for s in rsg hg; do
  ref=$O/CET3_$s; [ -f $ref/a/eos_table.txt ] && ls $ref/a/bin/*.m1.*.bin >/dev/null 2>&1 || ref=$O/CET2_$s
  [ -f $ref/a/eos_table.txt ] && ls $ref/a/bin/*.m1.*.bin >/dev/null 2>&1 || ref=$O/CET1_$s
  runs=; for t in CEP CEP0 CET1 CET2; do d=$O/${t}_$s; [ $d != $ref ] && ls $d/a/bin/*.m1.*.bin >/dev/null 2>&1 && runs="$runs $d"; done
  python3 $H/converr.py $ref $runs > $R/converr_$s.txt 2>&1
  [ -d $O/CET3c_$s/a/bin ] && [ -d $O/CEPc_$s/a/bin ] && python3 $H/converr.py $O/CET3c_$s $O/CEPc_$s > $R/converr_${s}_3c.txt 2>&1
done
ls -l $R
