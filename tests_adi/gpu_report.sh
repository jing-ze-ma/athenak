#!/bin/bash
cd "$(dirname "$0")"
for d in G_exp_16 G_sts_16 G_adi_16 G_exp_32 G_sts_32 G_adi_32; do
  [ -d $d/bin ] || continue
  echo -n "$d "; python3 measure_cs.py $d
  grep -o "|sum V de|/sum V|de| = [0-9.e+-]*" $d/log.txt 2>/dev/null | tail -1
done
for d in GS_sts_3.0e-3 GS_adi_3.0e-3 GS_sts_1.0e-5 GS_adi_1.0e-5; do
  [ -d $d/bin ] || continue
  echo "=== $d"; grep -o "max z_i = [0-9.e+-]*" $d/log.txt | head -1
  grep -o "substages = [0-9]*" $d/log.txt | head -1
  python3 amp_hist.py $d | tail -2
done
